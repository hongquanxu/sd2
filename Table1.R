## Table1.R : R codes for generating Table 1 and Figure 4
## Tian, Y. and Xu, H. (2025+). A Stratified L2-Discrepancy with Application to Space-Filling Designs. Journal of the Royal Statistical Society, Series B.
##
## require: sd2.R 
##
## Date: 7/12/25

library(DiceDesign) # discrepancyCriteria requires N>n
library(lhs)

source("sd2.R")		# sd2.y() 

# Table 1 and Figure 4: 64x63 random designs
N=64; n=63; 
x1=randomLHS(N,n); 
x2=1-x1*0.25; 
x3=(x1-1/2)*0.25+1/2; 
x4 = (x2+0.125)%%1; 
res=matrix(0,4,6)
designs=list(x1,x2,x3,x4)
par(mfrow=c(2,2)) # 
for(i in 1:4){
	x=designs[[i]]; # print(range(x))
	plot(x[,1], x[,2], xlim=c(0,1), ylim=c(0,1));
	Disc = discrepancyCriteria(x)
	sd_s2 = sd2.y(x,2,y=NULL)^(1/2) 		# default y=2/(m+1)
	sd_s3 = sd2.y(x,3,y=NULL)^(1/2)		# default y=2/(m+1)
	
	res[i,]=c(Disc$DisL2star, Disc$DisC2, Disc$DisW2, Disc$DisMix2, sd_s2, sd_s3) 
}
dimnames(res)[[2]]=c("D*","CD", "WD", "MD", "SD(s=2)", "SD(s=3)")
dimnames(res)[[1]]=c("P1", "P2", "P3", "P4")
round(log(res), 4)		# Table 1, the numbers differ slightly as the designs are random
