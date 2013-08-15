## Processing and plotting sims validating the use of the part method in capwire

## Matthew Pennell Aug 5/2013
library(ggplot2)
library(colorbrewer)

twod <- read.csv(file="two_discrete.csv", as.is=TRUE) 
thrd <- read.csv(file="three_discrete.csv", as.is=TRUE)
#twoc <- read.csv(file="two_continuous.csv", as.is=TRUE)
#thrc <- read.csv(file="three_continuous.csv", as.is=TRUE)


## fxn to calculate mean squared error 
## x is a vector of values
## y is the known value
calcMse <- function(x, y){
	
	mse <- 1/length(x) * sum(sapply(x, function(z) return((y - z)^2)))
	
	mse <- mse^(1/2) / y
	
	mse
	
}

## fxn to calculate biase
## x is a vector of values
## y is the known value

calcBias <- function(x, y){
	
	bias <- mean((x - y) / y)
	
	bias
}


## Combine all the data into a single data frame

sim.res <- rbind(twod, thrd)

## Get rid of first column
sim.res <- sim.res[,-1]

## Get make index for each simulation
## i is the rownumber
simIndex <- function(i, data){
	
	d <- data[i,]
	
	index <- d[c("Pop.size", "Samples", "No.classes", "Alpha.1", "Alpha.2", "Capture.sd")]
	
	rownames(index) <- NULL
	
	index
	
}


sim.index <- lapply(c(1:nrow(sim.res)), function(x) simIndex(x, sim.res))

sim.index <- do.call(rbind, sim.index)

sim.index <- as.matrix(sim.index)

index <- unique(sim.index)




## get all sims by index
## i is the index 
## indices is the indices for all data sets
sim.by.index <- function(i, sim.index, sim.res){
	
	
	sbi <- apply(sim.index, MARGIN=1, function(x) all(x == i))
	
	sbi <- sim.res[sbi,]
	
	sbi
	
}


## gather all sims identified by the index

sbi <- lapply(c(1:nrow(index)), function(x) sim.by.index(index[x,], sim.index, sim.res))



## for each matrix that is output, calculate bias and MSE

calcError <- function(x){
	
	y <- unique(x[,"Pop.size"])
	
	s <- unique(x[,"Samples"])
	
	a1 <- unique(x[,"Alpha.1"])
	
	a2 <- unique(x[,"Alpha.2"])
	
	n <- unique(x[,"No.classes"])
	
	sd <- unique(x[,"Capture.sd"])
		
	tirm.bias <- calcBias(x[,"Tirm.est"], y)
	
	part.bias <- calcBias(x[,"Part.est"], y)
	
	cond.bias <- calcBias(x[,"Cond.est"], y)
	
	tirm.mse <- calcMse(x[,"Tirm.est"], y)
	
	part.mse <- calcMse(x[,"Part.est"], y)
	
	cond.mse <- calcMse(x[,"Cond.est"], y)
	
	if (n == 2){
		alpha <- paste("[", a1, "]", sep=" ")
	}
	if (n == 3){
		alpha <- paste("[", paste(a1, a2, sep=","), "]", sep=" ")
	}
	
	if (sd == 0){
		type <- "Disc"
	} else {
		type <- "Cont"
	}
	
	res <- cbind.data.frame(y, "TIRM", tirm.bias, tirm.mse, s, alpha, n, type)
	res2 <- cbind.data.frame(y, "PART", part.bias, part.mse, s, alpha, n, type)
	res3 <- cbind.data.frame(y, "Cond. PART", cond.bias, cond.mse, s, alpha, n, type)
	names(res) <- names(res2) <- names(res3) <- c("Pop.size", "Estimator", "Bias", "MSE", "Samples", "Alpha", "N.class", "Type")
	
	res <- rbind(res, res2, res3)
	
	res
}



## Calculate across all simulated data sets

res <- lapply(sbi, function(z) calcError(z))

data <- do.call(rbind, res)

write.csv(data, file="part_valid.csv")







## Note: all plotting functions can be run without running the above as long as a file called part_valid.csv is in working directory

## N.B.: Change the facet labels so that they are alpha = [1,2,3]




## fxn for converting the names of the facet labels


## plotting

plot.bias.two <- function(){
	
	data <- read.csv(file="part_valid.csv")
	
	
	dd <- data[which(data[,"N.class"] == 2),]
	
	## change the levels of the facets
	dd <- within(dd, Alpha <- factor(Alpha, levels=c("[ 2 ]", "[ 4 ]", "[ 6 ]", "[ 8 ]", "[ 10 ]"), labels=c("alpha==1:2", "alpha==1:4", "alpha==1:6", "alpha==1:8", "alpha==1:10")))
	
	## change the levels of the estimators
	dd <- within(dd, Estimator <- factor(Estimator, levels=c("TIRM", "PART", "Cond. PART")))
	
			
	dd <- within(dd, Pop.size <- factor(Pop.size, levels=c(50,100,200), labels=c("N==50", "N==100", "N==200")))

	
	
	p <- ggplot(dd, aes(Samples, Bias, color=Estimator), environment=environment()) + geom_smooth(method="loess") + geom_point(size=3) + geom_hline(aes(yintercept=0, alpha=0.85)) + facet_grid(Alpha ~ Pop.size, labeller = label_parsed) + ggtitle("Comparison of bias - 2 rate classes")+ scale_colour_brewer(palette="Dark2")
	
	print(p)
	
}


plot.bias.three <- function(){
	
	
	data <- read.csv(file="part_valid.csv")
	
	
	dd <- data[which(data[,"N.class"] == 3),]
	
	## change the levels of the facets
	dd <- within(dd, Alpha <- factor(Alpha, levels=c("[ 2,4 ]", "[ 4,6 ]", "[ 4,8 ]", "[ 6,10 ]", "[ 6,12 ]"), labels=c("alpha==1:2:4", "alpha==1:4:6", "alpha==1:4:8", "alpha==1:6:10", "alpha==1:6:12")))
	
	
	## convert the levels of popsize to N = x
	dd <- within(dd, Pop.size <- factor(Pop.size, levels=c(50,100,200), labels=c("N==50", "N==100", "N==200")))

	## change the levels of the estimators
	dd <- within(dd, Estimator <- factor(Estimator, levels=c("TIRM", "PART", "Cond. PART")))
	

	
	p <- ggplot(dd, aes(Samples, Bias, color=Estimator), environment=environment()) + geom_smooth(method="loess") + geom_point(size=3) + geom_hline(aes(yintercept=0, alpha=0.85)) + facet_grid(Alpha ~ Pop.size, labeller = label_parsed) + ggtitle("Comparison of bias - 3 rate classes") + scale_colour_brewer(palette="Dark2")
	
	print(p)

}



plot.mse.two <- function(){
	
	data <- read.csv(file="part_valid.csv")
	
	
	dd <- data[which(data[,"N.class"] == 2),]
	
	
	## change the levels of the facets
	dd <- within(dd, Alpha <- factor(Alpha, levels=c("[ 2 ]", "[ 4 ]", "[ 6 ]", "[ 8 ]", "[ 10 ]"), labels=c("alpha==1:2", "alpha==1:4", "alpha==1:6", "alpha==1:8", "alpha==1:10")))
		
			
	dd <- within(dd, Pop.size <- factor(Pop.size, levels=c(50,100,200), labels=c("N==50", "N==100", "N==200")))

	## change the levels of the estimators
	dd <- within(dd, Estimator <- factor(Estimator, levels=c("TIRM", "PART", "Cond. PART")))
	

	p <- ggplot(dd, aes(Samples, MSE, color=Estimator), environment=environment()) + geom_smooth(method="loess") + geom_point(size=3) +  facet_grid(Alpha ~ Pop.size, labeller = label_parsed) + ggtitle("Comparison of errors - 2 rate classes") + ylab(expression(sqrt("MSE"))) + scale_colour_brewer(palette="Dark2")
	
	print(p)
	
}



plot.mse.three <- function(){
	
	
	data <- read.csv(file="part_valid.csv")
	
	
	dd <- data[which(data[,"N.class"] == 3),]
	
	## change the levels of the facets
	dd <- within(dd, Alpha <- factor(Alpha, levels=c("[ 2,4 ]", "[ 4,6 ]", "[ 4,8 ]", "[ 6,10 ]", "[ 6,12 ]"), labels=c("alpha==1:2:4", "alpha==1:4:6", "alpha==1:4:8", "alpha==1:6:10", "alpha==1:6:12")))
	
	
	## convert the levels of popsize to N = x
	dd <- within(dd, Pop.size <- factor(Pop.size, levels=c(50,100,200), labels=c("N==50", "N==100", "N==200")))

	## change the levels of the estimators
	dd <- within(dd, Estimator <- factor(Estimator, levels=c("TIRM", "PART", "Cond. PART")))
	

	
	p <- ggplot(dd, aes(Samples, MSE, color=Estimator), environment=environment()) + geom_smooth(method="loess") + geom_point(size=3) +  facet_grid(Alpha ~ Pop.size, labeller = label_parsed) + ggtitle("Comparison of errors - 3 rate classes") + ylab(expression(sqrt("MSE"))) + scale_colour_brewer(palette="Dark2")
	
	print(p)

}






dev.off()

pdf(file="bias_2_class.pdf", height=7, width=9)
plot.bias.two()
dev.off()

pdf(file="bias_3_class.pdf", height=7, width=9)
plot.bias.three()
dev.off()

pdf(file="mse_2_class.pdf", height=7, width=9)
plot.mse.two()
dev.off()

pdf(file="mse_3_class.pdf", height=7, width=9)
plot.mse.three()
dev.off()




