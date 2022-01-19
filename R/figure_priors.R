#' figure_priors.R
#' 
#' Code to generate figure:
#' - Figure S16: prior.pdf

# Set-up -----------
if(!"here" %in% installed.packages()){
    install.packages(here)
}

# Figure S16: prior.pdf -----------
pdf(here::here("figures", "priors.pdf"), width = 10, height = 4)
par(mfrow = c(1,3),las=1,yaxs="i",mar=c(5,1,3,1),cex.axis=1.3)
curve(dbeta(x,2,300),from=0,to=0.05,n=501,
      xlab="",ylab="",yaxt="n",ylim=c(0,115),lwd=3,col="red")
curve(dgamma(x,1.25,0.1),from=0,to=80,n=501,
      xlab="",ylab="",yaxt="n",ylim=c(0,0.065),lwd=3,col="blue")
curve(dbeta(x,80,20),from=0.5,to=1,n=501,
      xlab="",ylab="",yaxt="n",ylim=c(0,10.5),lwd=3,col="green3")
dev.off()