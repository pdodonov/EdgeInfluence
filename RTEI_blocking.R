##############################################################
# Randomization Test for Assessing Edge Influence (RTEI)     #
#   including a blocking factor                              #
##############################################################
# This text is analogous to a mixed model or a paired t-test #
#   in that it assumes that sampling points are not          #
#   independent, but grouped according to some factor.       #
##############################################################
# The grouping factor is usually the transect where the      #
#   samples were collected.                                  #
##############################################################

# Data for the analysis: should have three columns:
# 1) Transect - representing the transect or blocking factor;
# 2) Distance - the distance corresponding to each sample;
# 3) Value - The value of the response variable for which 
#      distance of edge influence (DEI) is to be assessed.

# First, let's simulate some data with eight transects and several distances:

dists <- c(0,5,10,15,20,30,40,50,60,80,100,120,150,180)
Ndists <- length(dists)
Ntrans <- 8
dists2 <- rep(dists, times=Ntrans)
transects <- rep(c(1:8), each=Ndists)

set.seed(17)

response <- 15+5*exp(-dists2/15) + rnorm(length(dists2),0,0.5)

plot(response~dists2, xaxt="n")
axis(at=dists, side=1)

# In this simulated data, the response variable is higher at the edge
#  and stabilizes around 30 m. This is a realistic scenario for 
#  Neotropical vegetation.

# Now let us introduce random variation, with some transects having
#  much higher values than other transects

response2 <- response + transects*5

plot(response2 ~ dists2, pch=16) # The pattern is much less clear

plot(response2 ~ dists2, pch=21, bg=transects) # But it can be seen within each transect

# Combine into one data frame

data.full <- data.frame(Transect = transects, Distance = dists2, Value = response2)

# Randomizations will be performed within each transect
# Test statistics: mean (edge - interior) differences (so NOT difference in edge and interior means)

distance <- dists2
block <- transects
y <- response2
dist.int.min <- 120
Nperm <- 5000

RTEI.blocking <- function(block, distance, y, dist.int.min, Nperm) {
    dists <- unique(distance)
    Ndists <- length(dists)
    Ntrans <- length(unique(block))
    result <- matrix(ncol=8, nrow=length(dists))
    colnames(result) <- c("Distance", "Mean", "SD", "MEI", "Signif_below", "Signif_above", "Signif_abs", "Signif_2tailed")
    which.int <- which(distance >= dist.int.min)
    y.int <- y[which.int]
    trans.int <- block[which.int]
    mean.int <- aggregate(y.int ~ trans.int, FUN=mean)[,2]
    for(i in 1:Ndists) {
        print(dists[i])
        y.edge <- y[distance == dists[i]]
        trans.edge <- block[distance == dists[i]]
        trans.edge <- block[distance == dists[i]]
        mean.edge <- aggregate(y.edge ~ trans.edge, FUN=mean)[,2]
        dif.real <- mean(mean.edge - mean.int)
        dif.perm <- numeric(Nperm)
        dif.perm[1] <- dif.real
        if (dists[i] < dist.int.min) {
            for (j in 2:Nperm) {
                dif.temp <- numeric(Ntrans)
                for(k in 1:Ntrans) {
                    trans.now <- unique(block)[k]
                    y.edge.now <- y.edge[trans.edge == trans.now]
                    y.int.now <- y.int[trans.int == trans.now]
                    Nedge.now <- length(y.edge.now)
                    data.now <- c(y.edge.now, y.int.now)
                    data.now.rand <- sample(data.now)
                    dif.temp[k] <- mean(data.now.rand[1:Nedge.now]) - mean(data.now.rand[(1+Nedge.now):length(data.now.rand)])  
                }
                dif.perm[j] <- mean(dif.temp)
            }
        }
        result[i,"Distance"] <- dists[i]
        result[i, "Mean"] <- mean(y.edge)
        result[i, "SD"] <- sd(y.edge)
        if (dists[i] < dist.int.min) {
            result[i, "MEI"] <- (mean(y.edge) - mean(y.int)) / (mean(y.edge) + mean(y.int))
            result[i, "Signif_below"] <- sum(dif.perm <= dif.real) / Nperm
            result[i, "Signif_above"] <- sum(dif.perm >= dif.real) / Nperm
            result[i, "Signif_abs"] <- sum(abs(dif.perm) >= abs(dif.real)) / Nperm
            result[i, "Signif_2tailed"] <- min(c(result[i, "Signif_below"], result[i, "Signif_above"]))*2
        }
    }
    return(result)
}


# Testing
RTEI.blocking(block = data.full$Transect, distance=data.full$Distance, y=data.full$Value, dist.int.min=120,Nperm=5000)

# Seems to work :-)                                                                                         




































