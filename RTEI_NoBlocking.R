# This is a function for depth of edge influence - DEI analysis, based on permutation of critical values obtained using edge and interior values.
### Parameters included in the function: 
#x (dataset with both edge and interior values); 
#Ntrans (number of edge transects); 
#Nedge (number of edge measurements in each transect); 
#Nint (number of interior measurements in each transect); 
#Nperm (how many permutations for each distance), 
#alpha (the desired significance level), 
#probs (whether you want the significance values at each distance (TRUE) or just whether it is significant at the desired alpha level (FALSE; 1 for significant).
###Required input: one file with edge values, with the first column having the distances (or a rank of distances), and the second column, which may have any name, with the variables.
###Output: five columns, and one row for each distance from edge (the first row beign the one closest to the edge). 
#Column 1: mean value at the distance; 
#Column 2: MEI, calculated as (edge - int)/(edge + int); 
#Column 3: one-tailed significance for that distance being significantly lower than the interior; 
#Column 4: one-tailed significance for it being higher than the interior; 
#Column 5: two-tailed significance. In probs=FALSE, only the one-tailed results are given. 
###Please cite the following material when using this script:
#Harper, K.A. and S.E. Macdonald. 2011. Quantifying distance of edge influence: a comparison of methods and a new randomization method. Ecosphere. 2: article 94. 
#Dodonov P., Harper K. A. and Silva-Matos D. M. (2013) The role of edge contrast and forest structure in edge influence: vegetation and microclimate at edges in the Brazilian cerrado. Plant Ecol. 214, 1345-59.
###Instructions:
#1) Copy the entire function into the console in R, or open it as a script
#2) Open your dataset as an object in R, using, e.g., dataset = read.csv(file.choose(),sep="",dec=".")
#3) Analyze your data by typing RTEI(dataset). 
#4) If your number of transects and of edge and interior plots differs from the default values, you have to input them into the function with e.g. RTEI(dataset, Ntrans=7, Nedge=10, Nint=4)

RTEI <- function (x, Ntrans=5, Nedge=12, Nint=3, Nperm=10000, alpha=0.05, probs=TRUE) {
	# Step 1: sort according to distances from the edge
	coisa <- x[order(x[,1]),]
	# Step 2: Separate the dataset into two, one for the edge values and one for the interior values.
	edge <- Ntrans*Nedge  #number of edge measurement
	int <- Ntrans*Nint   #number of interior measurements
	coisa.edge <- coisa[1:edge,] 
	#Edge values are transformed into a matrix
	coisa.edge <- t(matrix(data=coisa.edge[,2], nrow=Ntrans, ncol=Nedge))
	coisa.int <- coisa[edge+1:int,2]
	#Create an object with each row containg edge values and interir values
	joined <- matrix (nrow=Nedge, ncol=Ntrans+Ntrans*Nint)
	for (i in 1:Nedge) joined[i,]=c(coisa.edge[i,],coisa.int) 
	# Create a function that calculates the mean difference of a vector, permutes it and recalculates it.
	permute.MEI <- function (x) {
		permutation <- numeric(0)
		perm.MEI <- numeric(0)
		mean.edge <- mean(x[1:Ntrans],na.rm=TRUE)
		mean.int <- mean(x[(Ntrans+1):length(x)],na.rm=TRUE)
		real.MEI <- (mean.edge - mean.int) / (mean.edge + mean.int) 
		Nperm <- Nperm-1
		for (i in 1:Nperm) {
			permutation <- sample (x)
			perm.mean.edge <- mean(permutation[1:Ntrans],na.rm=TRUE)
			perm.mean.int <- mean(permutation[(Ntrans+1):length(permutation)],na.rm=TRUE)
			perm.MEI[i] <- (perm.mean.edge - perm.mean.int)/(perm.mean.edge + perm.mean.int)
		}
		perm.MEI <- c(perm.MEI,real.MEI)
		if (probs==TRUE) {
			prob.below <- length(perm.MEI[perm.MEI<=real.MEI])/length(perm.MEI)
			prob.above <- length(perm.MEI[perm.MEI>=real.MEI])/length(perm.MEI)
			output <- c(mean(x[1:Ntrans],na.rm=TRUE),real.MEI,prob.below,prob.above,min(prob.below,prob.above)*2)
			output
		}	else {
			lim.inf <- quantile(perm.MEI,(alpha/2),na.rm=TRUE)
			lim.sup <- quantile(perm.MEI,(1-alpha/2),na.rm=TRUE)
			sign.MEI.inf <- ifelse(real.MEI<=lim.inf,1,0)
			sign.MEI.sup <- ifelse (real.MEI>=lim.sup,1,0)
			output <- c(mean(x[1:Ntrans],na.rm=TRUE),real.MEI,sign.MEI.inf,sign.MEI.sup)
			return(output)
		}
	}
	answer <- apply(joined,1,permute.MEI)
	answer <- t(answer)
	colnames(answer) <- c("Mean (edge)", "MEI", "p(edge>=int)", "p(edge<=int)", "p(edge=int)")
	rownames(answer) <- unique(coisa[,1])[1:Nedge]
	interior <- aggregate(coisa.int,list(dist=sort(rep(1:Nint, Ntrans))), FUN=mean, na.rm=T)
	interior <- data.frame(interior[,2],rep(NA,Nint),rep(NA,Nint),rep(NA,Nint),rep(NA,Nint))
	colnames(interior) <- c("Mean (edge)", "MEI", "p(edge>=int)", "p(edge<=int)", "p(edge=int)")
	SDs <- aggregate(coisa[,2], list(dist=coisa[,1]), sd, na.rm=T)
	SDs <- SDs[,2]
	distances <- unique(coisa[,1])
	answer2 <- rbind(answer,interior)
	row.names(answer2) <- unique(coisa[,1])
	answer2$error <- SDs
	answer2 <- answer2[,c(1,6,2,3,4,5)]
	return(answer2)
	}
