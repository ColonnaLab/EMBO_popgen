#---------------------------------------------------
# Function to simulate data
#---------------------------------------------------
simData <- function(f, N){
  #simulate genotypes
  g <- rbinom(N, 2, f);
  
  #determin counts n_0, n_1, n_2
  n <- tabulate(g+1, nbin=3); 
  return(n);
};

#---------------------------------------------------
# Function to calculate log-likelihood and MLE
# Note: only proportional to full likelihood
#---------------------------------------------------
calcLL <- function(n, f){
  #Note: indexes in R are 1,2,3 for genotypes 0,1,2
  return((n[2] + 2*n[3]) * log(f) + (n[2] + 2*n[1]) * log(1-f));
};

calcMLE <- function(n){
  #return maximum likelihood estimate
  return((n[2]+2*n[3])/(2*sum(n)));
}

#---------------------------------------------------
# Function to calculate prior and posterior density
#---------------------------------------------------
priorDensity <- function(f, alpha, beta){
  return(dbeta(f, alpha, beta));
};

posteriorDensity <- function(f, n, alpha, beta){
  #Note: indexes in R are 1,2,3 for genotypes 0,1,2
  return(dbeta(f, alpha+n[2]+2*n[3], beta+n[2]+2*n[1]));
};

calcPriorPosterior <- function(n, alpha, beta, length.out=100000){
  #Note: beta distribution not defined at f=0 and f=1 if alpha,beta < 1
  f <- seq(0, 1, length.out=length.out)[2:(length.out-1)];
  prior <- priorDensity(f, alpha, beta);
  
  #get posterior
  posterior <- posteriorDensity(f, n, alpha, beta);

  #return
  return(data.frame(f=f, prior=prior, posterior=posterior));
}

#---------------------------------------------------
# Function to plot
#---------------------------------------------------
plotPriorPosterior <- function(densities, col='red', ylim=c(0, max(max(densities$prior), max(densities$posterior))), true_f=NA, MLE=NA, plotPrior=TRUE, plotPosterior=TRUE, plotMLE=TRUE, plotLegend=TRUE, add=FALSE, lwd=1.5){
  #open plot
  if(!add){
    plot(0, type='n', xlab="Allele frequency f", ylab="Density", xlim=c(0,1), ylim=ylim);
  }  
  
  #add MLE
  if(!is.na(MLE)){
    lines(rep(MLE, 2), par("usr")[3:4], lty=2, col=col);
    axis(side=3, at=MLE, labels=quote(hat("f")), col = NA, tick=FALSE, col.axis=col);
  }
  
  #plot true_f
  if(!is.na(true_f)){
    lines(rep(true_f, 2), par("usr")[3:4], lty=1, col='black');
    axis(side=3, at=true_f, labels=expression(f^'*'), col = NA, tick=FALSE);
  }
  
  #plot
  if(plotPrior){
    lines(densities$f, densities$prior, col=col, lty=2, lwd=lwd);
  }
  if(plotPosterior){
    lines(densities$f, densities$posterior, col=col, lty=1, lwd=lwd);
  }
  
  #add legend
  if(plotLegend){
    if(which(densities$posterior==max(densities$posterior))/length(densities$posterior) < 0.5){
      legend('topright', bty='n', lty=c(2,1), legend=c("Prior", "Posterior"), col=col, lwd=lwd);
    } else {
      legend('topleft', bty='n', lty=c(2,1), legend=c("Prior", "Posterior"), col=col, lwd=lwd);
    }
  }
};
  
#---------------------------------------------------
# Simulate and plot
#---------------------------------------------------
true_f <- 0.3;
n <- simData(true_f, N=5);

#or use data from example in script
#n <- c(12,11,2); 

#set prior
alpha <- 0.8;
beta <- 10;

#calculate and plot
densities <- calcPriorPosterior(n, alpha, beta);
MLE <- calcMLE(n);
plotPriorPosterior(densities, true_f=true_f, MLE=MLE);

