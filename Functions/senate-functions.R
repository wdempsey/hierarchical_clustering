## September 7, 2015
## Senate voting network clustering

#hc directory
setwd("C:/users/hcrane/dropbox/network sampling/hierarchical clustering/data/")
#wd directory (just a guess)
setwd("C:/users/wdempsey/dropbox/network sampling/hierarchical clustering/data/")

senate <- read.csv("senate.txt",header = FALSE)
votes <- senate[,4:ncol(senate)]

##Functions for getting a weighted and unweighted network from the above 
##'votes' data

##remove: returns a list of indices for missing votes, i.e., 0 entries
##in a vector 'x'.  These are fed into the 'network' function below to
##extract the relevant network data.

remove <- function(x){
	l <- length(x)
	rem <- c()
	for(i in 1:l){
		if(x[i] == 0){
			rem <- c(rem,i)
		}
	}
	rem
}

##network: takes a matrix of 'votes' and outputs a weighted network 'weight'
##along with a network 'counts' for which the (i,j) entry is the number of
##bills on which senators i and j both voted.

network  <- function(M){
	r <- nrow(M)
	G <- N <- matrix(0,r,r)
	
	for(i in 1:(r-1)){
		for(j in (i+1):r){
			re <- unique(c(remove(M[i,]),remove(M[j,])))
			if(is.null(re)){mat <- abs(M[i,]-M[j,])}else{
			mat <- abs(M[i,-re]-M[j,-re])}
			N[i,j] <- length(mat)
			G[i,j] <- N[i,j] - sum(mat)
		}
	}
list(weight = G, counts = N)
}

##project: projects a weighted network 'M' with counts 'N' to an unweighted network 
##by thresholding edge to 0 if the fraction of votes M[i,j]/N[i,j] <= p and 
##1 otherwise.

project <- function(M,N,p){
	P <- matrix(0,nrow=nrow(M),ncol= ncol(M))
	for(i in 1:(nrow(M)-1)){
		for(j in (i+1):nrow(M)){
			if(M[i,j]/N[i,j] > p){P[i,j] <- 1}
		}
	}
	P
}

## Garbage from before that should probably be ignored

senate <- t(matrix(c(R1,R2,R3,R4,R5,R6,R7,R8,R9,R10,R11,R12,R13,R14,R15,R16,R17,R18,R19,R20,
R21,R22,R23,R24,R25,R26,R27,R28,R29,R30,R31,R32,R33,R34,R35,R36,R37,R38,R39,R40,R41,R42,R43,R44,R45,R46,R47,R48,R49,R50,
D1,D2,D3,D4,D5,D6,D7,D8,D9,D10,D11,D12,D13,D14,D15,D16,D17,D18,D19,D20,
D21,D22,D23,D24,D25,D26,D27,D28,D29,D30,D31,D32,D33,D34,D35,D36,D37,D38,D39,D40,D41,D42,D43,D44,D45,D46,D47,D48,D49,D50),ncol=100))

sen.net <- network(senate)

