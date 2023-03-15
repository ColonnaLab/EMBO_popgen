#---------------------------------------------------
# Function to simulate data
#---------------------------------------------------
simData <- function(f, N){
  #simulate genotypes
  g <- rbinom(N, 2, f);
  
  #determin counts n_0, n_1, n_2
  n <- tabulate(g+1, nbin=3); 
  return(n);
}

#---------------------------------------------------
# Function to calculate log-likelihood and MLE
# Note: only proportional to full likelihood
#---------------------------------------------------
calcLL <- function(n, f){
  #Note: indexes in R are 1,2,3 for genotypes 0,1,2
  return((n[2] + 2*n[3]) * log(f) + (n[2] + 2*n[1]) * log(1-f));
}

calcMLE <- function(n){
  #return maximum likelihood estimate
  return((n[2]+2*n[3])/(2*sum(n)));
}

#---------------------------------------------------
# Function to calculate prior and posterior density
#---------------------------------------------------
priorDensity <- function(f, alpha, beta, log = FALSE){
  return(dbeta(f, alpha, beta, log = log));
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
# Functions to run MCMC
#---------------------------------------------------
runMCMC <- function(n, alpha, beta, f_initial, f_sd = 0.1, nIterations = 10000){
  cat(paste0("Running MCMC with ", nIterations, " iterations ..."));
  # prepare storage
  f <- numeric(nIterations);
  f[1] <- f_initial;
  
  #run MCMC chain
  for(i in 2:nIterations){
    # propose new f, mirror at 0 and 1
    new_f <- abs(rnorm(1, f[i-1], f_sd));
    if(new_f > 1){ new_f <- 2 - new_f; }
    
    # calculate hastings ratio
    h_log <- calcLL(n, new_f) - calcLL(n, f[i-1]) + priorDensity(new_f, alpha, beta, log = TRUE) - priorDensity(f[i-1], alpha, beta, log = TRUE);
    
    # accept or reject
    if(log(runif(1)) < h_log){
      f[i] <- new_f;
    } else {
      f[i] <- f[i-1];
    }
  }
  cat(" done!\n");
  cat(paste0("Acceptance rate was ", sum(f[1:(length(f)-1)] != f[2:length(f)]) / length(f)), ".\n");
  return(f);
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
n <- simData(true_f, N=25);

#set prior
alpha <- 1/2;
beta <- 1/2;

#calculate and plot
densities <- calcPriorPosterior(n, alpha, beta);
MLE <- calcMLE(n);
plotPriorPosterior(densities, true_f=true_f, MLE=MLE);

# run MCMC and add to plot
mcmc <- runMCMC(n, alpha, beta, true_f)
lines(density(mcmc), col='blue')
