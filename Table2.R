## Table2.R : R code for generating Table 2
## Tian, Y. and Xu, H. (2025+). A Stratified L2-Discrepancy with Application to Space-Filling Designs. Journal of the Royal Statistical Society, Series B.
##
## require: sd2.R and four 19x18 LHDs under the data folder
##
## Date: 7/12/25

library(DiceDesign) # discrepancyCriteria requires N>n

source("sd2.R")		# sd2()


# Table 2: 19x18 LHDs from Sun, Wang and Xu (2019)
Maximin =read.table("data/maximin19x18.txt", header = F) -1
MaxPro =read.table("data/maxpro19x18.txt", header = F) -1 
Uniform =read.table("data/ud19x18.txt", header = F) -1 
UPD=read.table("data/upd19x18.txt", header = F) -1 

designs=list(Maximin, MaxPro,  Uniform, UPD)
res=matrix(0,4,5)
for(i in 1:4){
	x=designs[[i]]; # print(range(x))
	x = (x+0.5)/nrow(x)		# transform to [0,1]
	Disc = discrepancyCriteria(x)
	sd_s2 = sd2(x,2)^(1/2) 
	sd_s3 = sd2(x,3)^(1/2)
	
	res[i,]=c(Disc$DisC2, Disc$DisW2, Disc$DisMix2, sd_s2, sd_s3) 
}
dimnames(res)[[2]]=c("CD", "WD", "MD", "SD(s=2)", "SD(s=3)")
dimnames(res)[[1]]=c("Maximin",   "MaxPro",  "Uniform", "UPD")
round((res), 4)		# Table 2
