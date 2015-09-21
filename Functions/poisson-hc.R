## August 16, 2015
## Clustering in interaction networks
####################################################

## Code for fitting the Poisson model to the Zachary data set
## Markov chain search functions from JASA paper are below.

## computes log of x(x+1)...(x+n-1)
log.factorial <- function(x,n){
	if(n==0){y <- 0}else{
	if(x<0){y <- log(-x)}else{if(x>0){y <- log(x)}else{y <- 0}}
	if(n>1){
	for(i in 1:(n-1)){
		y <- y + log(x+i)
	}}}
y
}


## Given a partition of labeled elements, will return an unlabeled partition
## in the form of a symmetric boolean matrix.

boolean.part <- function(partition){
M = max(partition)
N = length(partition)
B = matrix(0,N,N)
for(i in 1:M){
d = rep(0,N)
for(j in 1:N){
if(partition[j]==i){d[j]=1}
}
B = B + outer(d,d)
}
B
}


log.lik <- function(G,B,lambda.in=1,lambda.out=1){
	nr <- nrow(G)
	lik <- 0
	for(i in 1:(nr-1)){
		for(j in 2:nr){
			if(B[i]==B[j]){if(lambda.in[B[i]]>0){lik <- lik + G[i,j]*log(lambda.in[B[i]])-lambda.in[B[i]]}else{lik <- log(0)}}else{
				if(lambda.out>0){lik <- lik + G[i,j]*log(lambda.out) - lambda.out}else{lik <- log(0)}
			}
		}
	}
lik
}

## stochastic blockmodel
log.lik.sbm <- function(G,B,p.in=1,p.out=1){
	nr <- nrow(G)
	lik <- 0
	for(i in 1:(nr-1)){
		for(j in 2:nr){
			if(B[i]==B[j]){if(p.in[B[i]]>0){lik <- lik + G[i,j]*log(p.in[B[i]]) + (1-G[i,j])*log(1-p.in[B[i]])}else{lik <- log(0)}}else{
				if(p.out>0){lik <- lik + G[i,j]*log(p.out) + (1-G[i,j])*log(p.out)}else{lik <- log(0)}
			}
		}
	}
lik
}

lambda.max <- function(G,B){
	nr <- nrow(G)
	counts.in <- rep(0,max(B))
	counts.out <- 0
	lambda.in <- rep(0,max(B))
	lambda.out <- 0
	part <- boolean.part(B)
	part.comp <- 1-part
	G.part <- G*part
	G.comp <- G*part.comp

	for(i in 1:(nr-1)){
		for(j in 2:nr){
			counts.in[B[i]] <- counts.in[B[i]] + part[i,j]
			counts.out <- counts.out + part.comp[i,j]
			lambda.in[B[i]] <- lambda.in[B[i]] + G.part[i,j]
			lambda.out <- lambda.out + G.comp[i,j]
		}
	}
	for(i in 1:max(B)){
		if(counts.in[i] > 0){lambda.in[i] <- lambda.in[i]/counts.in[i]}else{lambda.in[i] <- 0}
	}
	if(counts.out > 0){lambda.out <- lambda.out/counts.out}else{lambda.out <- 0}
list(lambda.in = lambda.in, lambda.out = lambda.out)
}


p.max.sbm <- function(G,B){
	nr <- nrow(G)
	counts.in <- rep(0,max(B))
	counts.out <- 0
	p.in <- lambda.in <- rep(0,max(B))
	p.out <- lambda.out <- 0
	part <- boolean.part(B)
	part.comp <- 1-part
	G.part <- G*part
	G.comp <- G*part.comp

	for(i in 1:(nr-1)){
		for(j in 2:nr){
			counts.in[B[i]] <- counts.in[B[i]] + part[i,j]
			counts.out <- counts.out + part.comp[i,j]
			lambda.in[B[i]] <- lambda.in[B[i]] + G.part[i,j]
			lambda.out <- lambda.out + G.comp[i,j]
		}
	}
	for(i in 1:max(B)){
		if(counts.in[i] > 0){p.in[i] <- lambda.in[i]/counts.in[i]}else{p.in[i] <- 0}
	}
	if(counts.out > 0){p.out <- lambda.out/counts.out}else{p.out <- 0}
list(p.in = p.in, p.out = p.out)
}
####################################################################################################
## Search for optimal

## optimal Zachary
##B0  = c(1,1,1,1,1,1,1,1,2,2,1,1,1,1,2,2,1,1,2,1,2,1,2,2,2,2,2,2,2,2,2,2,2,2)

search.local <- function(G,B0,iter = 100,alpha.search=1,k=2){
		n <- length(B0)
		B.star <-  B0
		lam <- lambda.max(G,B0)
		p <- log.lik(G,B0,lam$lambda.in,lam$lambda.out)		
		C <- B0

		for(i in 1:iter){
			C0 <- Step.local(C,alpha.search,k)
			lam <- lambda.max(G,C0)
			q <- log.lik(G,C0,lam$lambda.in,lam$lambda.out)
			if(q>p){p <- q; B.star <- C0; C <- C0}else{if(q>log(0)){if(runif(1)>exp(q-p)){C <- C0}}
			}
		}
	list(B.s = B.star, p.s = p)
}



search.global <- function(G,B0,iter = 100,alpha.search=1,k=2){
		n <- length(B0)
		B.star <-  B0
		lam <- lambda.max(G,B0)
		p <- log.lik(G,B0,lam$lambda.in,lam$lambda.out)		
		C <- B0



		for(i in 1:iter){
			C0 <- Step(C,alpha.search,k)
			lam <- lambda.max(G,C0)
			q <- log.lik(G,C0,lam$lambda.in,lam$lambda.out)
			if(q>p){p <- q; B.star <- C0}
			C <- C0
			
		}
	list(B.s = B.star, p.s = p)
}



search.local.global <- function(G,B0,iter.loc = 100, iter.gl = 100,alpha.search=1,k=2){
		n <- length(B0)
		B.star <-  B0
		lam <- lambda.max(G,B0)
		p <- log.lik(G,B0,lam$lambda.in,lam$lambda.out)		
		C <- B0


		for(i in 1:iter.gl){
			C0 <- Step(C,alpha.search,k)
			out = search.local(G,C0,iter.loc,alpha.search,k)
			q <- out$p.s
			if(q>p){p = q; B.star <- out$B.s; C <- B.star}
		}
	list(B.s = B.star, p.s = p)
}

## SBM search

search.local <- function(G,B0,iter = 100,alpha.search=1,k=2){
		n <- length(B0)
		B.star <-  B0
		lam <- lambda.max(G,B0)
		p <- log.lik(G,B0,lam$lambda.in,lam$lambda.out)		
		C <- B0

		for(i in 1:iter){
			C0 <- Step.local(C,alpha.search,k)
			lam <- lambda.max(G,C0)
			q <- log.lik(G,C0,lam$lambda.in,lam$lambda.out)
			if(q>p){p <- q; B.star <- C0; C <- C0}else{if(q>log(0)){if(runif(1)>exp(q-p)){C <- C0}}
			}
		}
	list(B.s = B.star, p.s = p)
}



search.global <- function(G,B0,iter = 100,alpha.search=1,k=2){
		n <- length(B0)
		B.star <-  B0
		lam <- lambda.max(G,B0)
		p <- log.lik(G,B0,lam$lambda.in,lam$lambda.out)		
		C <- B0



		for(i in 1:iter){
			C0 <- Step(C,alpha.search,k)
			lam <- lambda.max(G,C0)
			q <- log.lik(G,C0,lam$lambda.in,lam$lambda.out)
			if(q>p){p <- q; B.star <- C0}
			C <- C0
			
		}
	list(B.s = B.star, p.s = p)
}



search.local.global <- function(G,B0,iter.loc = 100, iter.gl = 100,alpha.search=1,k=2){
		n <- length(B0)
		B.star <-  B0
		lam <- lambda.max(G,B0)
		p <- log.lik(G,B0,lam$lambda.in,lam$lambda.out)		
		C <- B0


		for(i in 1:iter.gl){
			C0 <- Step(C,alpha.search,k)
			out = search.local(G,C0,iter.loc,alpha.search,k)
			q <- out$p.s
			if(q>p){p = q; B.star <- out$B.s; C <- B.star}
		}
	list(B.s = B.star, p.s = p)
}



	
## Ewens partition

pitman.ewens <- function(n,alpha=0,theta=1){
	N <- 1;
	b <- c(1);
	pi <- c(1);
	
	if(n > 1){
		for(j in 2:n){
			new <- sample(N+1,1,prob=c(b-alpha,alpha*N+theta))
			pi <- c(pi,new);
			if(new == N+1){N <- N+1; b <- c(b,1)}else{b[new] <- b[new]+1}
		}
	}
pi
}

## ## adds element to B according to ewens probability
## k=0 corresponds to no upper bound on number of blocks

cocktail <- function(B,alpha=1,k=0){
	M <- max(B);
	b <- block.sizes(B);
	if(k==0){
		pr <- c(b,alpha)}else{
			pr <- c(b+alpha,(k-M)*alpha);
		}
	
	new <- sample(M+1,1,prob=pr)
	c(B,new);
}


ListPart <- function(L){
	n <- length(L);
	entry <- c(L[1]);
	E <- 1;
	part <- c(1);
	if(n > 1){
	for(i in 2:n){
		new <- 1;
		for(j in 1:E){
			if(L[i]==entry[j]){new <- 0; part[i] <- j}
		}
		if(new==1){entry <- c(entry,L[i]); E <- E+1; part[i] <- E};
	}}
	part;
}

ListPart.seq <- function(L){
	C <- ncol(L)
	S.part <- L;
	for(i in 1:C){
		r <- recused(L[,i])
		if(length(r)>0){	
			S.part[-r,i] <- ListPart(L[-r,i])
		}else{S.part[,i] <- ListPart(L[,i])}
	}
S.part
}

## Makes a local move in the state space

Step.local <- function(B,alpha=1,k=0){
	n <- length(B);
	U <- sample(n,1);
	B.n <- B;

	B.n[U] <- B[n]; B.n[n] <- B[U];
	B.n <- ListPart(B.n[-c(n)]);
	B.n <- cocktail(B.n,alpha,k);
	
	B.new <- B.n;
	B.new[U] <- B.n[n]; B.new[n] <- B.n[U];
	ListPart(B.new);
}

Step <- function(B, alpha=1, k=2, P=-1){
	n <- length(B);
	Bp <- c(1);
	if(sum(P)==-1){P <- rep(1,n)}
	B.n <- boolean.part(B)*boolean.part(P);

	for(i in 2:n){
		b <- B.n[i,1:(i-1)];
		N <- max(Bp);
		M <- rep(0,N+1);
		for(j in 1:N){
			C <- boolean.part(c(Bp,j));
			M[j] <- sum(b*C[i,1:(i-1)]);
		}
		M <- M + c(rep(alpha/k,N),alpha-alpha*N/k)
		Bp[i] <- sample(N+1,1,prob=M);
	}
	Bp
}


block.sizes <- function(part){
	M <- max(part)
	L <- length(part)
	s <- rep(0,M)
	for(i in 1:L){
		s[part[i]] <- s[part[i]]+1
	}
s
}

	
