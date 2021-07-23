#' helper_functions.R 
#' 
#' functions used to help process and analyze ouput from BEER simulations

#' integrate_vector()
#' 
#' Function to perform trapezoidal approximation given two vectors. 
#' 
#' @param x numeric vector of x values
#' @param y numeric vector with the same length as x
#' 
#' @return numeric value of the trapezoidal approximation 
integrate_vector <- function(x, y){
    
    if (length(y) != length(x)) {
        stop("The length of the vectors must be equal.")
    }
    
    n_points <- length(y)
    sorted_ind <- sort(x, index.return = TRUE)$ix
    sorted_x <- x[sorted_ind]
    sorted_y <- y[sorted_ind]
    
    # trapezoidal rule
    sum(0.5 * (sorted_y[-1] + sorted_y[-n_points]) * 
            (sorted_x[-1] - sorted_x[-n_points]))
}

#' get_roc()
#' 
#' Function to calculate the ppv, sens, and spec for various
#' cutoffs between 0 and 1. 
#' 
#' @param data data frame  with columns `prop_enriched` for the threshold and 
#' `Z` for the true enrichment status, and optional column `extra_info`
#' @param min_cutoff minimum cutoff, default to 0
#' @param max_cutoff maximum cutoff, default ot 1 - 1e-6
#' @param extra_info boolean indicating whether extra information should be used
#' in the classification at each cutoff. If `TRUE`, then the column `extra_info`
#' must be present in the data frame. 
#' 
#' @return data fram with columns `cutoff`, `ppv`, `sens`, and `spec`. 
get_roc <- function(data, min_cutoff = 0, max_cutoff = 1 - 1e-6, 
                    extra_info = FALSE){
    
    cutoffs <- seq(min_cutoff, max_cutoff, length.out = 1000)
    
    # Calculate ppv, sens, spec for each cutoff
    ppv <- sapply(cutoffs, function(x) {
        if(extra_info){
            predict <- (data$prop_enriched >= x & data$extra_info)
        } else {
            predict <- (data$prop_enriched >= x)
        }
        sum((data$Z == 1) & predict)/sum(predict)
    })
    
    spec <- sapply(cutoffs, function(x){
        if(extra_info){
            predict <- (data$prop_enriched >= x & data$extra_info)
        } else {
            predict <- (data$prop_enriched >= x)
        }
        sum((data$Z == 0) & !predict)/sum(data$Z == 0)
    })
    
    sens <- sapply(cutoffs, function(x){
        if(extra_info){
            predict <- (data$prop_enriched >= x & data$extra_info)
        } else {
            predict <- (data$prop_enriched >= x)
        }
        sum((data$Z == 1) & predict)/sum(data$Z == 1)
    })
    
    # Get AUC for ROC
    npoints <- length(sens)
    area_roc <- integrate_vector(1-spec, sens)
    
    return(cutoffs = data.frame(cutoff = cutoffs,
                                ppv = ppv,
                                sens = sens,
                                spec = spec, 
                                area_roc = area_roc))
}

#' get_legend()
#' 
#' Function given [here](https://stackoverflow.com/questions/12539348/ggplot-separate-legend-and-plot)
#' to extract the legend for plotting purposes. 
#' 
#' @param plot ggplot
#' @return ggplot legend
get_legend <- function(myggplot){
    tmp <- ggplot_gtable(ggplot_build(myggplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
}

#' mysqrt_trans()
#' 
#' Function given [here](https://stackoverflow.com/questions/47944992/ggplot2-removes-zero-when-using-scale-x-sqrt)
#' to add zero to the plot after sqrt transforming the x-axis
#' 
#' @param plot ggplot
#' @return ggplot legend
mysqrt_trans <- function() {
    scales::trans_new("mysqrt", 
                      transform = base::sqrt,
                      inverse = function(x) ifelse(x<0, 0, x^2),
                      domain = c(0, Inf))
}

#' penriched_fit()
#' 
#' Function that returns a data frame with the point estimate, and 95\% 
#' confidence intervals. 
#' 
#' @param model logistic regression model
#' @param covariates data frame of covariates
#' @return data frame with columns for the covariate, point estimate, lower CI
#' and upper CI.
penriched_fit <- function(model, covariates){
    
    # Predict based on the model
    prediction <- predict(model, covariates, type = "link", se.fit = TRUE)
    pred_lower <- prediction$fit - 1.96*prediction$se.fit
    pred_upper <- prediction$fit + 1.96*prediction$se.fit
    
    # Transform logit to probabilities
    ppred <- 1/(1 + exp(-prediction$fit))
    plower <- 1/(1 + exp(-pred_lower))
    pupper <- 1/(1 + exp(-pred_upper))
    
    # Return covariates with prediction + 95 CI added
    bind_cols(covariates, 
              data.frame(predict_p = ppred, 
                         lower_ci = plower, 
                         upper_ci = pupper))
}

#' as_df()
#' 
#' Function to convert a PhIPData object into a tidy data frame 
#' @param phip_obj PhIPData object
#' @return a tibble 
as_df <- function(phip_obj, metadata = FALSE){
    
    n_samples <- ncol(phip_obj)
    n_peps <- nrow(phip_obj)
    
    assay_df <- vapply(assays(phip_obj), as.vector, numeric(prod(dim(phip_obj))))
    metadata_df <- if(metadata){
        data.frame(est_c = rep(phip_obj$c, each = n_peps), 
                   est_pi = rep(phip_obj$pi, each = n_peps))
    } else { NULL }
    
    bind_cols(data.frame(sample = rep(colnames(phip_obj), each = n_peps), 
                         peptide = rep(rownames(phip_obj), times = n_samples),
                         beads = rep(phip_obj$group, each = n_peps), 
                         n = rep(librarySize(phip_obj), each = n_peps)),
              metadata_df, 
              as_tibble(assay_df))
}


