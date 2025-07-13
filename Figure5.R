## Figure5.R : R code for generating Figure 5 and Figure S1
## Tian, Y. and Xu, H. (2025+). A Stratified L2-Discrepancy with Application to Space-Filling Designs. Journal of the Royal Statistical Society, Series B.
##
## require SD2.R and R packages UniPro, MaxPro, SLHD, parallel

## It may take a few minutes to run.

## Date: 7/12/25

## install UniPro package from GitHub
if (!require(devtools)) install.packages("devtools")
if (!require(UniPro)) devtools::install_github("oonyambu/UniPro")

library(UniPro)
library(MaxPro)
library(SLHD)
require(parallel) # mclapply, works for MacBook Pro

source("sd2.R")		# sd2()

get_design <- function(type, n, k)
{	# levels: 0, ..., n-1
	x = switch(type,
      	MaxPro = MaxPro::MaxProLHD(n, k)$Design * n  - 0.5,
		LHD = replicate(k, sample(n)-1),
    		UniPro = UniPro::UniPro(n, k)$Design,		# slow for large designs
  		Maximin = SLHD::maximinSLHD(1, n, k)$Design -1,		
  	)
  	x
}

eval_design=function(x, y=1)
{ # 3 criteria only
	maxpro=MaxPro::MaxProMeasure( x/nrow(x) )  # scale x to [0,1]
	upd2 = UniPro::unipromeasure( round(x) )/1000 # restore unipro value, x needs to be integers 
	w2=sd2.y(x, q=2, y=y)
	c(maxpro=maxpro, unipro=upd2, sd2=w2^0.5) 
}



plot_crt3=function(n=23, k=5, nreps=10, designs=c('LHD', "Maximin", 'MaxPro',  'UniPro'), y=1)
{ 
	## generate design and compute the crteria values
	df = NULL
	for(i in 1:nreps){
		a = lapply(designs, function(type) get_design(type, n, k))
		res=t(sapply(a, eval_design, y=y))
		res = as.data.frame(res)
		dimnames(res)[[1]]=designs
		res$design=designs
		df = rbind(df, res)
	} 
	
	## make scatter plots
	cols = as.numeric(factor(df$design))
	des = sort(unique(df$design))
	vars = dimnames(df)[[2]]

	for(i in 1:2) for(j in (i+1):3){
		cor1= cor(df[,i], df[,j])
		main=paste0("LHD(", n, ",", k, ") (", round(cor1,2), ")"	) 
		plot(df[,i], df[,j], col=cols, pch=cols, xlab=vars[i], ylab=vars[j], main=main)
		legend("bottomright", legend=des, col=1:7, pch=1:7)
	}	
	invisible(df)
}


Figure5=function(n=23, nreps=50)
{ 
	# comparing 3 criteria and 4 types of designs 
#	file=paste0("cmp-lhd", n,".pdf")
#	pdf(file, w=8.4, h=5.6)
	par(mfrow=c(2,3))
	par(mar=c(5,4, 2, 1)+0.1)		# default c(5,4,4,2)+0.1
	a=plot_crt3(n=n, k=5,  nreps=nreps, designs=c('LHD', "Maximin", 'MaxPro',  'UniPro'))
	a=plot_crt3(n=n, k=n-1, nreps=nreps, designs=c('LHD', "Maximin", 'MaxPro',  'UniPro'))
#	dev.off()	
}


cmp_sd2y=function(n=20, k=10, nreps=10, designs=c('LHD', "Maximin", 'MaxPro'))
{ 
	require(parallel) # mclapply, works for MacBook Pro
 	mc.cores = if (.Platform$OS.type == "windows")   1   else  parallel::detectCores()
#  	mc.cores = if (.Platform$OS.type == "unix")   parallel::detectCores() else 1

	## generate design and compute the crteria values
	fn = function(i)  simplify2array( lapply(designs, function(type) sd2.y(get_design(type, n, k))^0.5 ))
	res = mclapply(1:nreps, fn, mc.cores = mc.cores)
	df = t(simplify2array(res))
	df = as.data.frame(df)
	dimnames(df)[[2]]=designs
	main=paste0("LHD(", n, ",", k, ")"	) 
	boxplot(df, main=main, ylab="sd2")
	df
} 

Figure_S1=function()
{ # examine sd2.y() with y=2/(m+1) for m=c(10,20,50,100)
# UniPro is slow and not effective in constructing 101*100 UPDs

#	system.time( get_design("MaxPro", 101, 100)) # 6s
#	system.time( get_design("Maximin", 101, 100)) # 10s
#	system.time( get_design("UniPro", 101, 100)) # 64s

	par(mfrow=c(2,4))
	par(mar=c(5,4, 2, 1)+0.1)		# default c(5,4,4,2)+0.1
	nreps=10; mm=c(10,20,50,100); designs=c('LHD', "Maximin", 'MaxPro')
	for(m in mm) 	a=cmp_sd2y(n=101, k=m,  nreps=nreps, designs=designs); 
	for(m in mm)		a=cmp_sd2y(n=200, k=m,  nreps=nreps, designs=designs); 
}

## Draw Figure 5 
Figure5()

## 
## Draw Figure S1 
## It may take a few minutes to run
date() 
Figure_S1()
date()

