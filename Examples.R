## Examples.R : R codes for reproducing Examples 2, 3, 5  
## Tian, Y. and Xu, H. (2025+). A Stratified L2-Discrepancy with Application to Space-Filling Designs. Journal of the Royal Statistical Society, Series B.
## require: sd2.R 
##
## Date: 7/12/25

library(DiceDesign) # discrepancyCriteria requires N>n
library(lhs)

source("sd2.R")		# sd2.y() 

### Codes to generate Tables and Figures 
##
###
Example2=function()
{  ## Example 2
	par(mfrow=c(1,2)) # 
	P1=rbind( c(0,0), c(1,1), c(2,4), c(3,5), c(4,2), c(5,3), c(6,6), c(7,7) )/8
	plot(P1+1/16, xlim=c(0,1), ylim=c(0,1))
	print( sqrt(sd2(P1)) )	# 0.2415566
	sqrt(sd2(P1+1/16)) # same SD2 with a small shift
	
	P2=rbind( c(0,0), c(2,3), c(3,6), c(1,5), c(6,2), c(4,1), c(5,4), c(7,7) )/8 
	plot(P2+1/16, xlim=c(0,1), ylim=c(0,1))
	print( sqrt(sd2(P2))	 ) # 0.138878
	sqrt(sd2(P2+1/16))	# same SD2 with a small shift
}


Table3=function(q=3, k=2)
{	## Table 3
	library(lhs)
	gf=lhs::create_galois_field(q^k); 
	x8=gf$times[,-1];	# remove col 1
	x8	# Table 3: GSOA(9,8,3^2,1)
}

Example5=function(q=3, k=2)
{ ## Example 5 using Table3()
	x8 = Table3(q=q, k=k)	# Table 3: GSOA(9,8,3^2,1)
	
	sd2.lb(q^k, (q^k-1), q, k)		# lower bound 1.148028
	print(sd2(x8, q, k))	#  sd2=1.148028, achieving the lower bound

	# GSOA(9,4,3^2,2) 
	sd2.lb(9, 4, q, k)		# lower bound 0.07583258
	print( sd2(x8[,c(1,3,4,5)], q, k) ) # with sd2= 0.07583258
}

Table4=function(no=1)
{	## Table 4: mixed-level designs
	## SOA(18,4,6,2+) from He, Cheng and Tang (2018)
#	library(DoE.base)	# GWLP
x1=matrix(c(0,0,0,2,2,2,4,4,4,1,1,1,3,3,3,5,5,5,0,2,4,0,2,4,0,2,4,1,3,5,1,3,5,1,3,5,0,2,4,2,4,0,4,0,2,1,3,5,3,5,1,5,1,3,0,2,4,4,0,2,2,4,0,1,3,5,5,1,3,3,5,1), 18, 4); #  

	# a GSOA(18,4,6,1) with the same GWLP as x1		
x2=matrix(c(0,0,1,2,2,2,5,5,5,1,1,0,3,3,3,4,4,4,1,2,4,0,2,4,0,3,5,0,2,4,1,3,5,1,3,5,0,3,5,3,5,1,4,0,3,1,2,4,2,4,0,5,1,2,1,3,5,4,1,3,2,4,1,0,2,4,5,0,2,3,5,0), 18, 4)
	
	if(no==1) x=x1 else	x = x2

#	write.table(t(x), row=F, col=F)
	x
}


	## Example 2
	Example2() # SD values for 2 designs
	
	# Example 3
	sd.kernel(2,3)	
	round( sd.kernel(3,2), 2)	


	## Example 5 and Table 3
	Table3()
	Example5()  # SD values for 2 designs

	## Table 4: mixed-level designs
	( x1 = Table4(no=1) )			# SOA(18,4,6,2+) and SD value
	sd2.mix(x1, c(3,2))^(1/2)		# 0.183595 with s1=3, s2=2
	
	( x2 = Table4(no=2) )	# GSOA(18,4,6,1) and SD value
	sd2.mix(x2, c(3,2))^(1/2)	#  0.1983595  with s1=3, s2=2

