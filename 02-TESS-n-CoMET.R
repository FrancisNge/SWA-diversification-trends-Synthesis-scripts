## template script for the TESS and CoMET analyses


#This method treats the number of specation-rate shifts, extinction-rate shifts,
#and mass-extinction events |as well as the parameters associated with these events|
#as random variables, and estimates their joint posterior distribution.

library(phyloch)
library(TESS)

#Selecting Piptostigmateae
tree<-read.tree("Xylopioa.168t.prune.beast.tree")
#remove_tip<-as.character(liste_especes[c(1:23,24:29,53),1])
#Piptostigmateae_coded<-drop.tip(tree,remove_tip)
plot.phylo(tree, cex=0.5)

#Preparing CoMET
{numExpectedMassExtinctions <- 1
  numExpectedRateChanges <- 1
  #Mean and standard deviation of the lognormal prior on the speciation rate in real space
  speciationPriorMu <- 0.2
  speciationPriorSigma <- 0.5
  #Same for extinction rate
  extinctionPriorMu <- 0.05
  extinctionPriorSigma <- 0.5
  
  # Transform the priors on the speciation rate into log space
  speciationRatePriorMean <- log((speciationPriorMu^2)
                                 /sqrt(speciationPriorSigma^2+
                                         speciationPriorMu^2))
  speciationRatePriorStDev <- sqrt( log(1+speciationPriorSigma^2
                                        /(speciationPriorMu^2)))
  
  # Transform the priors on the extinction rate into log space
  extinctionRatePriorMean <- log((extinctionPriorMu^2)
                                 /sqrt(extinctionPriorSigma^2+
                                         extinctionPriorMu^2))
  extinctionRatePriorStDev <- sqrt( log(1+extinctionPriorSigma^2
                                        /(extinctionPriorMu^2)))
  
  #beta distribution
  expectedSurvivalProbability <- 0.5
  pMassExtinctionPriorShape2 <- 100
  pMassExtinctionPriorShape1 <- - pMassExtinctionPriorShape2 *expectedSurvivalProbability /(expectedSurvivalProbability - 1)
  
  # Plot the density function of our beta distribution.
  curve(dbeta(x,shape1=pMassExtinctionPriorShape1,shape2=pMassExtinctionPriorShape2),n=1001,
        xlab="survival probability",ylab="density",las=1)
  # Plot the 95% prior interval on the survival probability.
  abline(v = qbeta(c(0.025,0.975),shape1=pMassExtinctionPriorShape1,
                   shape2=pMassExtinctionPriorShape2),lty=2)
  
  samplingFraction <- 0.9
}

#Running COMET
{
  tess.analysis(tree,
                empiricalHyperPriors = FALSE,
                initialSpeciationRate = speciationPriorMu,
                speciationRatePriorMean = speciationRatePriorMean,
                speciationRatePriorStDev = speciationRatePriorStDev,
                initialExtinctionRate = extinctionPriorMu,
                extinctionRatePriorMean = extinctionRatePriorMean,
                extinctionRatePriorStDev = extinctionRatePriorStDev,
                samplingProbability = samplingFraction,
                numExpectedRateChanges = numExpectedRateChanges,
                numExpectedMassExtinctions = numExpectedMassExtinctions,
                pMassExtinctionPriorShape1 = pMassExtinctionPriorShape1,
                pMassExtinctionPriorShape2 = pMassExtinctionPriorShape2,
                dir = "tess_analysis",
                MIN_ESS = 1000)}


#Results
output<- tess.process.output("tess_analysis",numExpectedRateChanges = numExpectedRateChanges,
                               numExpectedMassExtinctions = numExpectedMassExtinctions)
  
par(mfrow=c(3,2)) 
tess.plot.output(output,
                   fig.types = c("speciation rates",
                                 "speciation shift times",
                                 "extinction rates",
                                 "extinction shift times",
                                 "mass extinction Bayes factors",
                                 "mass extinction times"),
                   xlab = "millions years ago",lwd=1,las=1)

layout.mat <- matrix(1:3,nrow=2,ncol=1,byrow=TRUE)
layout(layout.mat)
tess.plot.output(output,fig.types = c("extinction shift times", "mass extinction times"))



effectiveSize(output$numSpeciationCategories)
## var value, what is the cut-off?

##### PLOT

#PLOT convergence

layout.mat <- matrix(1:6,nrow=3,ncol=2,byrow=TRUE)
layout(layout.mat)
tess.plot.singlechain.diagnostics(output,
                                  parameters = c("speciation rates",
                                                 "extinction rates",
                                                 "mass extinction times"),
                                  las=2)
