install.packages("foreach")
install.packages("doMC")
install.packages("ggtheme")

### libraries
library(ggplot2)
library(data.table)
library(foreach)  ### you will probably need to install this package
library(doMC)     ### you will probably need to install this package

registerDoMC(2)

## A first attempt.
gen1 <- c(1,0, 1,0, 1,0, 1,0, 1,0) ### initialize our population
gen2 <- sample(gen1, replace=T) ### The sample function draws samples from the input vector.

table(gen1) ## this function tabulates the allele counts for gen1..
table(gen2) ### and gen2.

table(gen1)/length(gen1) ## calculate frequencies
table(gen2)/length(gen2)

## Streamlining things a bit with `rbinom`
rbinom(n=1, size=10, prob=.5) ### this returns the number of successes.
rbinom(n=1, size=10, prob=.5)/10 ### this returns the frequency of the 'A' allele.

## One generation of drift using `rbinom`
gen1 <- .5
gen2 <- rbinom(n=1, size=10, prob=gen1)/10
gen2

## Multiple generations
gen1 <- .5
gen2 <- rbinom(n=1, size=10, prob=gen1)/10
gen3 <- rbinom(n=1, size=10, prob=gen2)/10
gen4 <- rbinom(n=1, size=10, prob=gen3)/10

### First we define our parameters
nGens <- 10
popSize=10
startingAlleleFreq=.5


wf_fun <- function(nGens, popSize, startingAlleleFreq, locus) {
  ### next, we generate a data.table that will hold our output.
tmp <- data.table(gen=c(1:nGens), af=-1, popSize=popSize, locus=locus) ### we set `af` equal to -1 because it is just a placeholder. We include popSize for bookkeeping
tmp[gen==1]$af <- startingAlleleFreq ### we need to initialize the very first generation at the specified allele frequency

### we run the for loop. We start at generation two because we need to use generation 1 as the starting generation
### this style of "for" loop can look in the "past" (i-1) and, if we wanted, the "future" (i+1)
for(i in 2:nGens) {
  tmp[gen==i]$af <- rbinom(1, popSize, tmp[gen==(i-1)]$af)/popSize
}
tmp
} 
  
wf_fun(nGens=100, popSize=10000, startingAlleleFreq = 0.5)

tmp2 <- foreach(popSize.i = c(100, 500, 1000, 10000, 100000, 1000000), .combine="rbind")%do%{
  foreach(locus.i=c(1:10), .combine="rbind")%do%{
    wf_fun(popSize = popSize.i, nGens=10, startingAlleleFreq = 0.5, locus = locus.i)
  }
}



ggplot(data=tmp2, aes(x=gen, y=af, group=locus)) + geom_line() + facet_grid(~popSize)
  ylab("Allele Frequency") + xlab("Generation")
  


## loops that do not care about the past and cannot see it. But, it is useful for independent simulations
test_fun <- function(popSize) {
  return(data.table(initialPopSize=popSize, newPopSize=popSize * 2))
}

tmp2 <- foreach(popSize.i=c(100, 500), .combine="rbind")%do%{
  test_fun(popSize=popSize.i)
}

ggplot(data=tmp2, aes(x=popSize, y=locus.i)) + geom_line() + ylim(0,1)


foreach(locus.i=c(1:10),combine="rbind")%do%{
  
### Your turn: combine the "for" and "foreach" loops, and likely a function that you write, 
### to simulate different population sizes and a few hundred loci per population size for several hundred generations


