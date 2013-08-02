## Capwire simulations

## Matthew Pennell and Craig Miller
## July 19, 2013

## Purpose: 
## Validate the Partitioning Method under the assumptions of two versus three capture classes


library(capwire)
## make sure that the latest version from github is installed
library(multicore)

## Simulating conditions:
## 1. Two discrete capture classes
## 2. Three discrete capture classes
## 3. Bimodal (two means with some variance)
## 4. Trimodal (three means with some variance)


## Simulating discrete classes
## Takes 3 arguments
## 1. alpha -- vector of differences in capturability
## 2. n.class -- number of individuals in each class (must be longer than alpha by 1)
## should be in order of capture class (i.e. first number, the easiest to capture)
## 3. n.samples -- total sample size

simCapwireD <- function(alpha, n.class, n.samples){
	
	if (length(alpha) != (length(n.class) -1)){
		return(print("number of alphas must correspond to number of classes"))
	}
	
	all.ind <- c(1:(sum(n.class)))
	
	## get sampling probabilities
	prob.samp <- rep(1, n.class[1])
	
	for (i in 2:length(n.class)){
		tmp <- rep(alpha[i-1], n.class[i])
		prob.samp <- c(prob.samp, tmp)
	}
	
	## draw from the population
	samp <- sample(all.ind, size=n.samples, replace=TRUE, prob=prob.samp)
	
	## get the number of individuals sampled
	ind.samp <- unique(samp)
	
	data <- data.frame()
	for (j in 1:length(ind.samp)){
		x <- length(samp[samp == ind.samp[j]])
		y <- c(ind.samp[j], x)
		data <- rbind(data, y)
	}
	
	data <- indToClass(data)
	return(data)

}




## Simulate from continuous distributions


## Create distribution function first

# Takes 3 arguments
## 1. alpha -- vector of differences in capturability
## 2. sd -- standard distribution of the normal dists (assumed to be the same across class)
## 3. n.class -- number of individuals in each class (must be longer than alpha by 1)
## should be in order of capture class (i.e. first number, the easiest to capture)



capRatesMM <- function(alpha, sd, n.class){
	
	## create a distribution for the first class
	
	cap.dist <- rnorm(10000 * n.class[1], mean=10, sd=sd)
	
	## add each rate class with values proportional to their relative numbers
	
	for (i in 2:length(n.class)){
		
		class.i <- rnorm(10000 * n.class[i], mean=10*alpha[i-1], sd=sd)
		cap.dist <- c(cap.dist, class.i)
	}
	
	## truncate the distribution at 0
	
	trunc.dist <- cap.dist[which(cap.dist >= 0)]
	
	## divide by 10 so that we remove the effect of the truncation on the relative values
	
	trunc.dist <- trunc.dist/10
	
	## output a function which will be used to simulate data
	
	rate.fxn <- function(n){
		
		sample(trunc.dist, size=n, replace=TRUE)
		
	}
	
	return(rate.fxn)
}





## A wrapper function which also samples from the distribution

## Takes 4 arguments
## The first three will just be fed into capRatesMM (used internally)

## 1. alpha -- vector of differences in capturability
## 2. sd -- standard distribution of the normal dists (assumed to be the same across class)
## 3. n.class -- number of individuals in each class (must be longer than alpha by 1)
## should be in order of capture class (i.e. first number, the easiest to capture)
## 4. n.samples -- total sample size


simCapwireC <- function(alpha, sd, n.class, n.samples){
	
	## make sampling fxn according to parameters
	
	rate.dist <- capRatesMM(alpha, sd, n.class)
	
	## simulate data set
	
	data <- simCapture(n=sum(n.class), s=n.samples, dist.func=rate.dist)
	
	return(data)
}








## Simulate then fit data to both tirm and part
## This is a general function for both discrete and continuous
## If sd == 0 (default): discrete
## If sd != 0: continuous

partTest <- function(alpha, sd=0, n.class, n.samples, max.pop){
	
	## Simulate discrete data
	
	if (sd == 0){
	
		sim.data <- simCapwireD(alpha, n.class, n.samples)
	}
	
	## Simulate continuous data
	
	if (sd != 0){
		
		sim.data <- simCapwireC(alpha, sd, n.class, n.samples)
	}
	
	## fit Tirm
	
	fit.tirm <- fitTirm(sim.data, max.pop)
	
	tirm.pop <- fit.tirm$ml.pop.size
	
	## fit Part
	
	fit.part <- fitTirmPartition(sim.data, max.pop)
	
	part.pop <- fit.part$ml.pop.size
	
	
	## get number of classes
	
	nom.class <- length(n.class)
	
	## add extra alpha/pop class for bookkeeping purposes
	
	if (length(alpha) == 1){
		alpha <- c(alpha, 0)
	}
	
	if (length(n.class) == 2){
		n.class <- c(n.class, 0)
	}
	
	output <- c(sum(n.class), tirm.pop, part.pop, n.samples, nom.class, alpha[1], alpha[2], n.class[1], n.class[2], n.class[3], sd)
	
	names(output) <- c("Pop.size", "Tirm.est", "Part.est", "Samples", "No.classes", "Alpha.1", "Alpha.2", "N.class.1", "N.class2", "N.class3", "Capture.sd")
	
	return(output)
	
}






## Simulations start here


## Two discrete populations

alpha <- c(2, 4, 6, 8, 10)

pop.size <- c(50, 100, 200)

## specify population subdivision. Proportion of individuals in first class
#sub.div <- c(0.25, 0.50, 0.75)

sample.size <- c(100, 200, 300, 400, 500)

two.dis <- data.frame()

for (i in 1:length(alpha)){
	for (j in 1:length(pop.size)){
			for (m in 1:length(sample.size)){
				
				class.one <- round(0.5 * pop.size[j])
				class.two <- pop.size[j] - class.one
				n.class <- c(class.one, class.two)
				
				tmp <- mclapply(c(1:100), function(x) partTest(alpha[i], sd=0, n.class=n.class, n.samples=sample.size[m], max.pop=300), mc.cores=24)
				
				tmp2 <- do.call(rbind, tmp)
				
				two.dis <- rbind(two.dis, tmp2)
				
				
		}
	}
}

## write the file
write.csv(two.dis, file="two_discrete_v2.csv")






# Three discrete classes


alpha <- c(2, 2, 4, 4, 6)
alpha.two <- c(3, 4, 6, 8, 10)

pop.size <- c(50, 100, 200)

sample.size <- c(100, 200, 300, 400, 500)

three.dis <- data.frame()

for (i in 1:length(alpha)){
	for (j in 1:length(pop.size)){
			for (m in 1:length(sample.size)){
				
				class.one <- round(1/3 * pop.size[j])
				class.two <- round(1/3 * pop.size[j])
				class.three <- pop.size[j] - class.one - class.two
				n.class <- c(class.one, class.two, class.three)
				
				alpha.one <- alpha[i]
				alpha.sec <- alpha.two[i]
				alpha.all <- c(alpha.one, alpha.sec)
				
				tmp <- mclapply(c(1:100), function(x) partTest(alpha.all, sd=0, n.class=n.class, n.samples=sample.size[m], max.pop=300), mc.cores=24)
				
				tmp2 <- do.call(rbind, tmp)
				
				three.dis <- rbind(three.dis, tmp2)
				
				
			}
	}
}

## write the file
write.csv(three.dis, file="three_discrete.csv")







# Two cont populations

alpha <- c(2, 4, 6, 8, 10)

pop.size <- c(50, 100, 200)

sample.size <- c(100, 200, 300, 400, 500)

two.cont <- data.frame()

for (i in 1:length(alpha)){
	for (j in 1:length(pop.size)){
			for (m in 1:length(sample.size)){
				
				class.one <- round(0.5 * pop.size[j])
				class.two <- pop.size[j] - class.one
				n.class <- c(class.one, class.two)
								
				tmp <- mclapply(c(1:100), function(x) partTest(alpha[i], sd=0.5, n.class=n.class, n.samples=sample.size[m], max.pop=300), mc.cores=24)
				
				tmp2 <- do.call(rbind, tmp)
				
				two.cont <- rbind(two.cont, tmp2)	
				
			}
	}
}

## write the file
write.csv(two.cont, file="two_continuous.csv")






## Three cont classes

alpha <- c(2, 2, 4, 4, 6)
alpha.two <- c(3, 4, 6, 8, 10)

pop.size <- c(50, 100, 200)

sample.size <- c(100, 200, 300, 400, 500)

three.cont <- data.frame()

for (i in 1:length(alpha)){
	for (j in 1:length(pop.size)){
			for (m in 1:length(sample.size)){
				
				class.one <- round(1/3 * pop.size[j])
				class.two <- round(1/3 * pop.size[j])
				class.three <- pop.size[j] - class.one - class.two
				n.class <- c(class.one, class.two, class.three)
				
				alpha.one <- alpha[i]
				alpha.sec <- alpha.two[i]
				alpha.all <- c(alpha.one, alpha.sec)
								
				tmp <- mclapply(c(1:100), function(x) partTest(alpha.all, sd=0.5, n.class=n.class, n.samples=sample.size[m], max.pop=300), mc.cores=24)
				
				tmp2 <- do.call(rbind, tmp)
				
				three.cont <- rbind(three.cont, tmp2)
				
				
		}
	}
}

## write the file
write.csv(three.cont, file="three_continuous.csv")







