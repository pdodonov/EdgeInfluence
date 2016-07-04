######################################
#Function for estimating Edge Contrast
#By Pavel Dodonov - pdodonov@gmail.com (feel free to email!)
#When using this function, please cite the following manuscript: Dodonov P., Harper K. A. & Silva-Matos D. M. (2013) The role of edge contrast and forest structure in edge influence: vegetation and microclimate at edges in the Brazilian cerrado. Plant Ecol. 214, 1345-59..
#And also please use citation(msm) to cite the msm library, which is used to generate the weighting function.
#This function, as is, considers four different contrasts at different distances from the edge, but can be easily modified to consider more or less contrasts.
#Parameters of the function:
#x: a sequence of numbers used to calculate the area below the contrast curve. Has to follow a random uniform distribution (or another random distribution) in order for the montecarlo approach to work.
#contrasts: a vector with four numbers: contrast_up_to_distance_1 ... contrast_up_to_distance_4, distance_1 ... distance_4 - e.g.: c(1,0.8,0.1,1,5,23,35,46)
#SD: the standard deviation used for the weighting function
#dist: up to what distance from edge are contrsts considered?
#The code below may be used to calculate the contrast for more than one edge (by means of a loop)

##############Specifying the objects
####Here you specify the objects used by the function
data <- read.table("clipboard",header=T,row.names=1) #A data.frame containing the set of contrasts for different edges; one edge per row.

dist.max <- 40 #will be the default of the dist argument in the function. Must be included here if you want to analyze more than one edge.
standdev <- 15 #will be the default of the SD argument in the function. Must be included here if you want to analyze more than one edge.
Nreplicate <- 1000 #number of random numbers used in the monte-carlo calculation of edge contrast.
######################
library(msm) #Needed to generate the truncated normal.
sqrtpi <- sqrt(pi)
#Here is the function to create the weighting function.
dweighting <-function (x,dist=dist.max, SD=standdev) {   #This is the function that created the PDF of the weighting factor
	#If you're wondering: it's a truncated normal with mean=0 and SD = dist.max/2.
	toone <- 1 / max(dtnorm(x,mean=0,sd=SD, lower=0, upper=dist))
	weighting.factor <- dtnorm(x,mean=0,sd=SD, lower=0, upper=dist)*toone
	return(weighting.factor)
}
#And here is the function that multiplies it by the contrasts at different distance
dcontrast <- function (x, contrasts=contrastes, dist=dist.max, SD=standdev) {
	SDcontr <- SD
	weighting.factor <- dweighting(x, SD <- SDcontr)
	contr1 <- contrasts[1]
	contr2 <- contrasts[2]
	contr3 <- contrasts[3]
	contr4 <- contrasts[4]
	dist1 <- contrasts[5]
	dist2 <- contrasts[6]
	dist3 <- contrasts[7]
	dist4 <- contrasts[8]
	contrast.unweighted <- numeric(1)
	contrast.weighted <- numeric(1)
	contrast.unweighted <- as.numeric(ifelse (x >= 0 & x < dist1, contr1, 
		ifelse (x >= dist1 & x < dist2, contr2,
		ifelse (x >= dist2 & x < dist3, contr3,
		ifelse (x >= dist3 & x < dist4, contr4, "error")))))
	contrast.weighted <- contrast.unweighted * weighting.factor
	return(contrast.weighted)
}
#And here is a code to calculate an area below the second function using the Monte Carlo approach
answer <- numeric(nrow(data))
for (i in 1:nrow(data)) {
	contrastes <- as.numeric(data[i,])
	replicates <- runif(Nreplicate,min=0,max=dist.max)
	weighted.contrast <- mean(dcontrast(replicates))
	weighting.factor.contrast <- mean(dweighting(replicates))
	final.contrast <- weighted.contrast / weighting.factor.contrast
	a <- seq(0,dist.max,by=dist.max/100)
	plot(a,dcontrast(a),type="l", main=i)
	answer[i] <- final.contrast 
	}
names(answer) <- row.names(data)
answer


		

