#-----------------------------
#function to simulate data
#-----------------------------
simulate <- function(N, f){
  g <- 0:2;
  p <- choose(2,g) * f^(2-g) * (1-f)^g;
  if(sum(p<0) > 0){
    print(f);
    print(p);
  }
  G <- c(sample(0:2, N, replace=TRUE, prob=p), g); #add a 0, 1 and 2 to make sure each category exists
  return(tabulate(G+1)-1);
};

simulateMany <- function(N, fPrime){
  nPrime <- sapply(fPrime, simulate, N=N);
  return(data.frame(f=fPrime, n_0=nPrime[1,], n_1=nPrime[2,], n_2=nPrime[3,]));
};

simFromPrior <- function(len, N, alpha, beta){
  fPrime <- rbeta(len, alpha, beta);
  return(simulateMany(N, fPrime));
};


#-----------------------------
#function to get prior and analytical posterior
#-----------------------------
posterior <- function(n, f, alpha, beta){
  return(dbeta(f, alpha + 2*n[1] + n[2] + 1, beta + n[2] + 2*n[3] + 1));
}

prior <- function(f, alpha, beta){
  return(dbeta(f, alpha, beta));
}

#-----------------------------
#function to do ABC
#-----------------------------
calcDistances <- function(n, simulations){
  #we use Euclidean norm
  d <- sqrt(colSums((t(simulations[,2:4]) - n)^2));
  return(d);
}

retainSimulations <- function(n, simulations, delta=5, maxToRetain=5000){
  #calculate distances
  d <- calcDistances(n, simulations);
  
  #retain sims below threshold
  numToRetain <- min(sum(d <= delta), maxToRetain);
  retained <- simulations[d <= delta,];
  return(retained[1:numToRetain,]);
}

retainBestSimulations <- function(n, simulations, numToRetain=5000){
  #calculate distances
  d <- calcDistances(n, simulations);
  
  #sort to get best sims
  dd <- sort(d, index.return=TRUE);
  
  #retain best sims
  retained <- simulations[dd$ix[1:numToRetain]];
  return(retained);
}

ABC_GLM <- function(n, retained, sd=0.0001, nPoints=200){
  f <- seq(0, 1, length.out=nPoints);
  
  #prepare matrices
  numSim <- length(retained[,1]);
  P <- matrix(retained$f, nrow=1);
  X <- matrix(c(rep(1, numSim), retained$f), byrow=FALSE, ncol=2);
  S <- as.matrix(t(retained[,2:dim(retained)[2]]));
  
  #fit glm
  C_c0 <- S %*% X %*% solve(t(X) %*% X);
  c0 <- C_c0[,1];
  C <- C_c0[,2];
  R <- t(S) - X %*% t(C_c0);
  Sigma_s <- 1/(numSim - 1) * (t(R) %*% R)
  
  #calculate posterior
  Sigma_s_inv <- solve(Sigma_s + diag(0.01,length(n)));
  Sigma_theta_inv <- matrix(1/sd);
  T <- solve(t(C) %*% Sigma_s_inv %*% C + Sigma_theta_inv);
  T_inv <- solve(T);
  
  c_j <- numeric(numSim);
  t_j <- numeric(numSim);
  for(j in 1:numSim){
    v_j <- t(C) %*% Sigma_s_inv %*% (n - c0) + Sigma_theta_inv * retained$f[j];
    t_j[j] <- T %*% v_j;
    c_j[j] <- -0.5 * (retained$f[j] * Sigma_theta_inv * retained$f[j] - t(v_j) %*% T %*% v_j);
  }
  c_j <- exp(c_j - max(c_j));
  
  posterior <- numeric(nPoints);
  for(i in 1:nPoints){
    aa <- f[i]-t_j;
    posterior[i] <- sum(c_j * exp(-0.5 * aa * as.numeric(T_inv) * aa));
  }
  
  #normalize
  d <- 1/(nPoints-1);
  A <- sum(posterior[2:(nPoints-1)]) * d + (posterior[1] + posterior[nPoints]) * d/2;
  posterior <- posterior / A;
  
  #return
  return(data.frame(f=f, posterior=posterior));
}

#---------------------------
# Functions to plot ABC results
#---------------------------
L1 <- function(x, trueP, abcP){
  return(sum(abs(trueP - abcP)) * (x[2]-x[1]) / 2);
}

openPlot <- function(n, alpha, beta, main="", delta=NA, colors=NA){
  #get true prior and analytic posterior
  f <- seq(0, 1, length.out=1000);
  posteriorDens <- posterior(n, f, alpha, beta);
  priorDens <- prior(f, alpha, beta);
  
  #open plot
  par(xaxs='i', yaxs='i', las=1);
  plot(0, type='n', xlim=c(0,1), ylim=c(0, max(posteriorDens)*1.05), xlab="Alelle frequency f", ylab="Density", main=main);
  
  #plot prior and analytic posterior
  polygon(c(f, f[1]), c(posteriorDens, posteriorDens[1]), border=NA, col='gray92');
  lines(f, priorDens, lwd=1.25, col='black', lty=2);
  
  #legend
  if(sum(is.na(delta)) == 0){
   legend('topright', bty='n', col=colors, lwd=1.5, legend=as.expression(lapply(1:length(delta), function(x) bquote(delta==.(delta[x])))))
  }
}

plotPosteriorDensityFromSimulations <- function(simulations, col='red', lwd=1.5, lty=1){
  dens <- density(simulations$f, adjust=1.25);
  lines(dens$x, dens$y, col=col, lwd=lwd, lty=lty);
}

plot_ABC_GLM_Posterior<-function(reg, col='red'){
  lines(reg$f, reg$posterior, col=col, lwd=1.5, lty=1);
}

#-----------------------------
# Calculating coverage property
#-----------------------------
coverageProperty <- function(pseudo_obs, simulations, delta=0){
  #prepare storage
  coverage <- numeric(dim(pseudo_obs)[1]);
  
  #conduct ABC estimations
  for(i in 1:dim(pseudo_obs)[1]){
    n_tmp <- as.numeric(pseudo_obs[i,2:4]);
    retained <- retainSimulations(n_tmp, simulations, delta=delta, maxToRetain=5000);
    coverage[i] <- sum(retained$f < pseudo_obs$f[i]) / length(retained$f);
  }
  
  #return
  return(coverage);
}

coverageProperty_ABC_GLM <- function(pseudo_obs, simulations, delta=0){
  #prepare storage
  coverage <- numeric(dim(pseudo_obs)[1]);
  
  #conduct ABC estimations
  for(i in 1:dim(pseudo_obs)[1]){
    n_tmp <- as.numeric(pseudo_obs[i,2:4]);
    retained <- retainSimulations(n_tmp, simulations, delta=delta, maxToRetain=5000);
    reg <- ABC_GLM(n_tmp, retained, sd=0.0001);
    coverage[i] <- sum(reg$posterior[reg$f<pseudo_obs$f[i]]) / sum(reg$posterior);
  }
  
  #return
  return(coverage);
}

#-----------------------------
#Run ABC inference
#-----------------------------
true_f <- 1/3; #true allele frequency
alpha <- 1/2; beta <- alpha; #prior

# observed data: either simulate or use values from script.
N <- 25; #number of individuals
n <- simulate(N, true_f);
#n <- c(2, 11, 12);
N <- sum(n);

#general graphics settings
par(mfrow=c(1,2), xaxs='i', yaxs='i', las=1, mgp=c(2.5,0.66, 0), mar=c(3.5,4,1.5,1));

# 0) Generate simulations
#-----------------------------
#generate simulations
simulations <- simFromPrior(len=500000, N, alpha, beta);

#visualize simulations
par(mfrow=c(1,3), xaxs='i', yaxs='i', las=1, mgp=c(2.5,0.66, 0), mar=c(3.5,4,1.5,1));
these <- 1:10000
plot(simulations$f[these], simulations$n_0[these], xlab="f", ylab="n_0");
plot(simulations$f[these], simulations$n_1[these], xlab="f", ylab="n_1");
plot(simulations$f[these], simulations$n_2[these], xlab="f", ylab="n_2");

# 1) ABC-REJ
#------------------
delta <- c(1, 5, 10, 15);
colors <- c("dodgerblue", "purple", "orange2", "red");

#open plot
par(mfrow=c(1,1), xaxs='i', yaxs='i', las=1, mgp=c(2.5,0.66, 0), mar=c(3.5,4,1.5,1));
openPlot(n, alpha, beta, main="ABC-REJ", delta=delta, colors=colors);

#get posterior for different thresholds delta
for(i in 1:length(delta)){
  retained <- retainSimulations(n, simulations, delta=delta[i], maxToRetain=5000);
  plotPosteriorDensityFromSimulations(retained, col=colors[i]);
}

# 2) ABC-GLM
#------------------
#use simulations and delta values from 1) ABC-REJ
#open plot
openPlot(n, alpha, beta, main="ABC-GLM", delta=delta, colors);

#get posterior for different thresholds delta
for(i in 1:length(delta)){
  retained <- retainSimulations(n, simulations, delta=delta[i], maxToRetain=5000);
  reg <- ABC_GLM(n, retained);
  plot_ABC_GLM_Posterior(reg, col=colors[i]);
}

#-----------------------------
# Asses coverage property
#-----------------------------
#function to plot coverage
plotcoverage <- function(coverage, main){
  ks <- ks.test(coverage, "punif")
  hist(coverage, xlab="Coverage", main=paste0("ABC-REJ (p = ", round(ks$p.value, 4), ")"));
}

#make pseuso-observed simulations
pseudo_obs <- simFromPrior(1000, N, alpha, beta);

#using simulations conducted above under 1)
coverage <- coverageProperty(pseudo_obs, simulations, delta=1);
hist(coverage, xlab="Coverage", main="ABC-REJ");
ks.test(coverage, "punif")

#coverage property of ABC-GLM using simulations conducted above under 1)
coverage <- coverageProperty_ABC_GLM(pseudo_obs, simulations, delta=15);
hist(coverage, xlab="Coverage", main="ABC-GLM");
ks.test(coverage, "punif")