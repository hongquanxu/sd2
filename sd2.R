## SD2.R : Supplemental R codes for
## Tian, Y. and Xu, H. (2025+). A Stratified L2-Discrepancy with Application to Space-Filling Designs. Journal of the Royal Statistical Society, Series B.
##
##
## key functions: sd2=sd2.w, sd2.y,  sd2.mix
##    	sd2=sd2.w : stratified L2-discrepancy using general weights as in Theorem 1
##		sd2.y: stratified L2-discrepancy using exponential weights (7) or special weights in Theorem 4
##    	sd2.mix: general stratified L2-discrepancy in Section 7
##
## The following R files reproduce the numerical results in Tian and Xu (2025+)
## Examples.R : Examples 2, 3, 5
## Table1.R : Table 1 and Figure 4
## Table2.R : Table 2
## Figure5.R : Figure 5 and Figure S1
##  
###
## Date: 7/12/25

## function to generate a full s^k factorial design
full.fd=function(s,k)
{ # generate s^k full factorial design, 4/26/20
	if(k==1) return(matrix(0:(s-1),ncol=1))
	x0=full.fd(s,k-1)
	xk=cbind(0,x0)
	for(i in 1:(s-1)){
		xk=rbind(xk, cbind(i,x0))
	}
	dimnames(xk)=list(1:s^k, 1:k)
	xk
}

### functions for the stratified L2-discrepancy 
sd.kernel=function(q, k, w=rep(1,k), w0=1)
{ # w is a vector of k weights, 3/7/21
# return a q^k * q^k matrix K(x,y)
	Fd = full.fd(q, k) # q^k full factorial
	Ker = matrix(0, q^k, q^k)
	for(x in 1:q^k) for(y in 1:x){
		k1 = w0
		for(j in 1:k){
			if(Fd[x,j] != Fd[y,j]) break  # no need to continue
			k1 = k1 + w[j]/q^j
		}  
		Ker[x,y] = Ker[y,x] = k1
	}
	Ker   # include the first term 1
}

sd2 = sd2.w=function(x, q=2, k=floor(log(nrow(x)+1e-8, q)), w=rep(1,k), w0=1)
{ # stratified L2-discrepancy with weights, Theorem 1
# x: N*n design, within [0,1) or scaled to [0,1)
# q, k: stratification parameters. Each dimension is divided into q^k intervals.
# w: a vector of weights with length k

	N=nrow(x); n=ncol(x); 
	if(max(x)-min(x)>=1){ # 
		s = max(x) - min(x) + 1	 # number of levels
		x = x - min(x)  # coded as 0:(s-1)
	# change levels (0:(s-1)) to (0,1)
		x = (x+0.5)/s
	}
	
	if(max(x)-min(x)<1){ # 0<=x<1
		x = floor(q^k * x) + 1 # convert to [0,q^k-1] +1
	}	
	Kw = sd.kernel(q, k, w, w0)
	
	res = 0
	for(a in 1:N) for(b in a:N){
		pk = 1
		for(j in 1:n) pk = pk * Kw[x[a,j], x[b,j]] 
		if(b>a) res = res + 2* pk  # (b,a)
		else res = res + pk  # a==b
	}
	res1= w0
	for(j in 1:k) res1 = res1 + w[j]/q^(2*j)
	res/N^2 -(res1)^n
}

sd2.y=function(x, q=2, k=floor(log(nrow(x)+1e-8, q)), y=NULL, adjust=F)
{ # stratified L2-discrepancy with exponenital weighting scheme (7) or special weights in Theorem 4 
# x: N*n design, y can be a complex number
# q, k: stratification parameters. Each dimension is divided into q^k intervals.
# 0 <y <1

	if(is.null(y))	y = 2/(ncol(x)+1)	# set default y=2/(m+1), 7/8/25
	if(adjust == T){	# adjust weights as in Theorem 4
		w=(q^2*y)^(1:k); 
		w[k]=(q^2*y)^k/(1-y)
	}
	else 	w=(y)^(1:k); 	# exponenital weighting scheme (7)
	sd2.w(x,q,k, w=w, w0=1)
}

sd2.lb=function(n, m, s, p, w=rep(1,p))
{ # lower bound in Theorem 5. 
 # w is the vector of weights. Assume w0=1.
 
	g=function(k, w){
		res = 1
		if(k < p) for(i in 1:(p-k)) res = res + w[i]/s^i
		res
	}
	
	lambda=n/s^p
	n0 = (lambda-1)*m/(n-1)
	nk = lambda * s^((1:p)-1) * (s-1)*m/(n-1)	# a vector of length p
	lb1 = - ( 1 + sum(w[1:p] * s^(-2*(1:p))) )^m
	g0w = g(0, w)
	lb2 = 1
	for(k in 1:p) lb2 = lb2 * g(k,w)^nk[k] 
	lb2 = ( g0w^m + (n-1) * g0w^n0 * lb2 )/n
	lb1 + lb2
}

### general stratified L2-discrepancy with mixed levels; see Section 7

sd.kernel.mix=function(ss,  w=rep(1,k), w0=1)
{ # w is a vector of k weights
## ss=c(s1,s2,...) is a vector of levels with length k, 7/17/24
## stratification by s1, s1*s2, s1*s2*s3, etc for each dimension
# return a qk * qk matrix K(x,y), where qk=prod(ss)
# Note: sd.kernel.mix(rep(q,k)) is equivalent to sd.kernel(q,k)
# q=3; k=2; system.time(x2<-sd.kernel.mix(rep(q,k))); system.time(x1<-sd.kernel(q,k)); range(x2-x1) 
	k = length(ss)
	qk=prod(ss)
	Ker = matrix(0, qk, qk)
	for(x in 1:qk) for(y in 1:x){
		k1 = w0
		for(j in 1:k){
#			if(Fd[x,j] != Fd[y,j]) break  # no need to continue
#			k1 = k1 + w[j]/q^j
		# collapse qk levels to prod(ss[1:j]) levels
			ssj = prod(ss[1:j])
			sdiv = qk / ssj
			xj = floor((x-1)/sdiv)
			yj = floor((y-1)/sdiv)
			if(xj != yj) break  # no need to continue
			k1 = k1 + w[j]/ssj
		}  
		Ker[x,y] = Ker[y,x] = k1
	}
	Ker   # include the first term 1
}

sd2.mix=function(x, ss, w=rep(1,k), w0=1)
{ # general stratified L2-discrepancy with mixed levels; see Section 7
# x: N*n design, ss: levels
# w is a vector of k weights
# ss=c(s1,s2,...) is a vector of levels with length k, 7/17/24
# stratification by s1, s1*s2, s1*s2*s3, etc for each dimension

# sd2.mix is an extend of sd2.w, 7/17/24
# x=soa32x8; sd2(x,2,3); sd2.mix(x, rep(2,3))
# x=soa18x4; sd2(x,2,3); sd2.mix(x, c(3,2));  sd2.mix(x, c(2,3))
	k = length(ss)
	N=nrow(x); n=ncol(x); 
	if(max(x)-min(x)>=1){ # 3/8/21
		s = max(x) - min(x) + 1	 # number of levels
		x = x - min(x)  # coded as 0:(s-1)
	# change levels (0:(s-1)) to (0,1)
		x = (x+0.5)/s
	}
	
	if(max(x)-min(x)<1){ # 0<=x<1
		x = floor(prod(ss) * x) + 1 # convert to [0,q^k-1] +1
	}	
	Kw = sd.kernel.mix(ss, w, w0)
	
	res = 0
	for(a in 1:N) for(b in a:N){
		pk = 1
		for(j in 1:n) pk = pk * Kw[x[a,j], x[b,j]] 
		if(b>a) res = res + 2* pk  # (b,a)
		else res = res + pk  # a==b
	}
	res1= w0
	for(j in 1:k) res1 = res1 + w[j]/prod(ss[1:j])^2
	res/N^2 -(res1)^n
}


