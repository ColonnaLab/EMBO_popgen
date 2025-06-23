
sims <- read.csv("mosquito-task2.csv", head=T)

# check prior distributions
# x11() # only needed if on command line
par(mfrow=c(2,2))
hist(sims$N1)
hist(sims$N2)
hist(sims$T_split)
hist(sims$MigRate)

# remove simulations with NaN for some summary stats!

# find useful summary stats which correlate with T_split
cor(sims$Fst, sims$T_split)
cor(sims$dxy, sims$T_split)
cor(sims$segsites1, sims$T_split)
cor(sims$segsites2, sims$T_split)
cor(sims$pi1, sims$T_split)
cor(sims$pi2, sims$T_split)
cor(sims$tajima1, sims$T_split)
cor(sims$tajima2, sims$T_split)

# load observed summary stats
obs <- read.csv("mosquito-observed.csv", head=T)

# check if simulated retained summary stats contain the observed one
quantile(sims$Fst); cat(obs$Fst)
quantile(sims$segsites1); cat(obs$segsites1)
quantile(sims$segsites2); cat(obs$segsites2)

# merge obs with retained sims to scale them
sumstats <- scale(rbind(obs[c(1,3,4)],sims[,c(5,7,8)]))
# testing without real data
# sumstats <- sims[,c(5,7,8)]

library(abc)

est <- abc(target=sumstats[1,], param=sims$T_split, sumstat=sumstats[-1,], tol=0.05, method="rejection")

# check distances in the acceptance region
hist(est$dist)
abline(v=max(est$dist[which(est$region)]), lty=2)

# posterior
# x11() # only needed if on command line
par(mfrow=c(2,1))
hist(est$unadj.values, freq=FALSE, xlim=range(sims$T_split), col=rgb(0,0,1,1/4), main="Posterior probability", xlab="Split time")
# MAP
map <- mean(est$unadj.values)
abline(v=map, lty=2)

# confidence intervals
hpd <- quantile(x=est$unadj.values, probs=c(0.025,0.975))
abline(v=hpd, lty=3)

# prior
hist(sims$T_split, freq=FALSE, xlim=range(sims$T_split), col=rgb(1,0,0,1/4))







