

# observations: MMFM

thetas <- seq(0,1,0.01)

# likelihood
like <- c()
for (theta in thetas) like <- c(like, (theta^3)*(1-theta)^1)

plot(thetas,like, type="l")

# MLE
thetas[which.max(like)]

# LRT: H0: theta<=0.05
tD <- 2*log(max(like)/like[6])
tD
# p-value
1-pchisq(tD,1)

# priors

# skeptical
belief1 <- 1/thetas
belief1[1] <- NA
belief1 <- belief1/sum(belief1, na.rm=T)

par(mfrow=c(2,3))
plot(thetas,like,type="l")
plot(thetas,belief1,type="l")
plot(thetas,belief1*like,type="l")
thetas[which.max(belief1*like)]

# very skeptical
belief2 <- 1/thetas^3
belief2[1] <- NA
belief2 <- belief2/sum(belief2, na.rm=T)

#par(mfrow=c(1,3))
plot(thetas,like,type="l")
plot(thetas,belief2,type="l")
plot(thetas,belief2*like,type="l")
thetas[which.max(belief2*like)]




