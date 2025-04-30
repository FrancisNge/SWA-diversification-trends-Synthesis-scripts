## Template script to plot BAMM results and run STRAPP
## input tree file and csv traits (SWA vs SEA) files are specific to each plant group

library(BAMMtools)
library(geiger)
library(ape)

tree <- read.tree("xxx")
Ntip(tree)

edata <- getEventData(tree, eventdata = "Corymbia.nu.event_data.txt", burnin=0.2)
shift_probs <- summary(edata)

##Error in Casuarinaceae subsp. (edata doesn't like full stops), go to text and replace ' with nothing.

#Rate through time
plotRateThroughTime(edata, ratetype="speciation")
plotRateThroughTime(edata, ratetype="extinction")

#export as 4 x5 frmae

plot(tree)

##Refer to Bamm tutorial 2015 Phyloseminar
plot.bammdata(edata, legend=T, spex='netdiv')

##did the MCMC runs converge?
mcmcout <- read.csv("Casuarinaceae.nu.mcmc_out.txt")

plot(mcmcout$logLik ~ mcmcout$generation)

#discard the first 10% of samples as burnin:
burnstart <- floor(0.1 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]

#it is good to check the effective sample sizes of the log-likelihood and the number of shift events present in each sample. 
#We’ll do this using the coda library for R:
library(coda)
effectiveSize(postburn$N_shifts)
effectiveSize(postburn$logLik)
#In general, we want these to be at least 200 (and 200 is on the low side, but might be reasonable for very large datasets).



plot.bammdata(edata, lwd=1, legend=T)
axis(1, at=max(branching.times(tree))-0:80, labels=0:80)
css <- credibleShiftSet(edata, expectedNumberOfShifts=1, threshold=5, set.limit = 0.95)
css$number.distinct
summary(css)
plot.credibleshiftset(css)

#plot best shift configuration (i.e. most frequent shift configuration)
best <- getBestShiftConfiguration(edata, expectedNumberOfShifts=1, threshold=5)
plot.bammdata(best, lwd=1.25)
addBAMMshifts(best, cex=2)




##STRAPP
# CSV input file
tempfile = read.csv("Allo.Casuarina_only.Geohisse_traits.SW.match.csv",  sep = "\t", header = FALSE)
trait.vector = tempfile$V2
lapply(trait.vector, as.numeric)
names(trait.vector) <- row.names(tempfile)
names(trait.vector) <- tree$tip.label


edata_diversification_ladderized <- getEventData(tree, eventdata = "Loganiaceae.nu.event_data.txt", burnin=0.1)            


# If you get a warning about negative log values, you can put logrates = FALSE
library(parallel)
environment_vs_diversification = traitDependentBAMM(edata, pm.sorted, 1000, rate = "net diversification", return.full = FALSE, method = "spearman", logrates = FALSE, two.tailed = TRUE, traitorder = NA, nthreads = 4) 
#environment_vs_speciation = traitDependentBAMM(edata, trait.vector, 1000, rate = "speciation", return.full = FALSE, method = "spearman", logrates = TRUE, two.tailed = TRUE, traitorder = NA, nthreads = 4) 
#environment_vs_extinction = traitDependentBAMM(edata, trait.vector, 1000, rate = "extinction", return.full = FALSE, method = "spearman", logrates = TRUE, two.tailed = TRUE, traitorder = NA, nthreads = 4) 

#can change method, spearman, pearson, mann-whitney, kruskal

environment_vs_diversification.speciation = traitDependentBAMM(edata, pm.sorted, 1000, rate = "speciation", return.full = FALSE, method = "spearman", logrates = FALSE, two.tailed = TRUE, traitorder = NA, nthreads = 4) 

#output txt files
write.table(environment_vs_diversification, file="Goodeniaceae.geo.STRAPP.txt", quote=FALSE, sep="\t")
write.table(environment_vs_diversification.speciation, file="Goodeniaceae.geo.speciation.STRAPP.txt", quote=FALSE, sep="\t")



