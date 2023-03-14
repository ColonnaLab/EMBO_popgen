# WFsim.r  simple simulations of the Wright-Fisher model

# Initialize popsize (N), number of samples (samp), number of 
# generations (ngen) and starting frequency (startfreq)

# Note this is for a haploid population of size N.

N<-200
nsamp<-100
ngen<-200
startfreq<-.5

#Wright-Fisher is simply recurrent binomial sampling over generations

x<-rbinom(nsamp,N,startfreq)
p<-x/N
pcur<-p

for (i in 1:ngen){
x<-rbinom(nsamp,N,p)
p<-x/N
pcur<-rbind(pcur,p)
}

#Now plot these trajectories
gen<-seq(1,ngen+1)
for (k in 1:nsamp){
plot(gen,pcur[,k],type="l",ylim=c(0,1))
par(new=TRUE)
}


