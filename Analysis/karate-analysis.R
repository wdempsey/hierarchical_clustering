## August 16, 2015
## Analysis of Zachary

### Data from Karate Network

setwd("C:/users/hcrane/dropbox/network sampling/hierarchical clustering/data")
data = read.table("zachary-weighted.txt")

## Projected Karate Network
nr <- nrow(data)
proj <- data
for(i in 1:nr){
	for(j in 1:nr){
		if(data[i,j]>0){proj[i,j] <- 1}else{proj[i,j] <- 0}
	}
}



## "True" Clustering
B0  = c(1,1,1,1,1,1,1,1,2,2,1,1,1,1,2,2,1,1,2,1,2,1,2,2,2,2,2,2,2,2,2,2,2,2)
B0 <- rep(1,34); B0[19] <- 2
B0 <- rep(1,34)
pp <- p.max.sbm(proj,B0)
ll <- lambda.max(data,B0,same=TRUE)
log.lik.sbm(proj,B0,pp$p.in,pp$p.out)
log.lik(data,B0,ll$lambda.in,ll$lambda.out)

log.lik.sbm(proj,rep(1,34),p.in=c(.11),p.out=c(1/2))

B0 <- c(1,2,1,2,1,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1)
## Dempsey max
B0 <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,1,1,2,1,2,1,2,2,2,2,2,2,2,2,2,2,2,2)
pp <- p.max.sbm(proj,B0)
ll <- lambda.max(data,B0)
log.lik.sbm(proj,B0,pp$p.in,pp$p.out)
log.lik(data,B0,ll$lambda.in,ll$lambda.out)

B0 <- opt$B.s
ll <- lambda.max(data,B0)
log.posterior(data,B0,1,1,ll$lambda.in,ll$lambda.out)


B0 <- opt$B.s
begin <- proc.time()[3]
opt <- search.local.global(data,B0,iter.loc=50,iter.gl=10,k=2,same=TRUE)
opt
end <- proc.time()[3]
(end-begin)/60

begin <- proc.time()[3]
opt.sbm <- search.local.global.sbm(proj,rep(1,34),iter.loc=20,iter.gl=40)
opt.sbm
end <- proc.time()[3]
(end-begin)/60

B00 <- rep(1,34)
B00[[8] <- 2
log.lik(data,B0,alpha.out=1/2,alpha.in=1/2,p=1)
log.lik(data,B00,alpha.out=1/2,alpha.in=1/2,p=1)