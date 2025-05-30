##Template script for running GeoHiSSE

#### GeoHiSSE code ####

# For more information, see:
# Caetano, D. S., O'Meara, B. C., & Beaulieu, J. M. (2018). 
# Hidden state models improve state‐dependent diversification approaches, 
# including biogeographical models. Evolution, 72(11), 2308-2324.

library(ape)
library(hisse)
library(parallel)

group = "xxx" # name the group
tree <- read.tree(file.choose()) # load tree file 
dist <- read.csv(file.choose()) # load distribution file 

plot(tree)
Ntip(tree)

# Preparing data - areas have to be as 0 (11 - widespread), 
# 1 (10, endemic of first area) 
# and 2 (01, endemic of second area)

areas <- as.data.frame(rep(1, nrow(dist)))
dist <- cbind(dist, areas)
colnames(dist)[4] <- "area"


#####

#warnings()

#####
for (i in 1:length(dist$area)){
  if (dist[i, "southwest"] == 1 && dist[i, "non_southwest"]  == 1){
    dist[i, "area"] = 0 
  }
  if (dist[i, "southwest"] == 0 && dist[i, "non_southwest"]  == 1){
    dist[i, "area"] = 1
  }
  if (dist[i, "southwest"] == 1 && dist[i, "non_southwest"]  == 0){
    dist[i, "area"] = 2
  }
}
states<-dist[,c("Species", "area")]
##change X to "species' in own data

table(states$area) # check if species-richness in each range make sense

# 2 - cr "endemic"
# 1 - non-cr "endemic"
# 0 - widespread

# Load sampling fraction for the group
#sf<-c(1,1,1) # e.g. if it's fully sampled
sf<-c(0.9,0.9,0.9) # e.g. if it's fully sampled

#
phy=tree
dat=states

# We used the same 18 models of Caetano et al. (2018) (plus a second set of models including jump dispersal) - see their original publication for more information

###############################################################################
## Block of GeoSSE-like models.
## Here extirpation is linked to range reduction.
###############################################################################

## Model 1 - Dispersal parameters vary only, no range-dependent diversification. 
speciation <- c(1,1,1)
extirpation <- c(1,1)
trans.rate <- TransMatMakerGeoHiSSE(hidden.areas=0, include.jumps=FALSE
                                    , separate.extirpation=FALSE) 
mod1 <- GeoHiSSE(phy, dat, f=sf, speciation=speciation, extirpation=extirpation,
                 hidden.areas=FALSE, trans.rate=trans.rate, assume.cladogenetic=TRUE) 


## Model 2. Canonical GeoSSE model, range effect on diversification 
speciation <- c(1,2,3)
extirpation <- c(1,2)
trans.rate <- TransMatMakerGeoHiSSE(hidden.areas=0, include.jumps=FALSE
                                    , separate.extirpation=FALSE) 
mod2 <- GeoHiSSE(phy, dat, f=sf, speciation=speciation , extirpation=extirpation,
                 hidden.areas=FALSE , trans.rate=trans.rate, assume.cladogenetic=TRUE)



## Model 3. Heterogeneous diversification, not tied to range evolution.
## Assumes three distinct diversification rates.
## Dispersion parameters across hidden areas are the same.
speciation <- c(1,1,1,2,2,2,3,3,3)
extirpation <- c(1,1,2,2,3,3)
trans.rate <- TransMatMakerGeoHiSSE(hidden.areas=2, make.null=TRUE
                                    , include.jumps=FALSE, separate.extirpation=FALSE) 
mod3 <- GeoHiSSE(phy, dat, f=sf, speciation=speciation, extirpation=extirpation,
                 hidden.areas=TRUE, trans.rate=trans.rate, assume.cladogenetic=TRUE) 



## Model 4. Heterogeneous diversification, tied to range evolution. 
## Assumes 6 distinct diversification rates.
speciation <- c(1,2,3,4,5,6)
extirpation <- c(1,2,3,4)
trans.rate <- TransMatMakerGeoHiSSE(hidden.areas=1, include.jumps=FALSE , 
                                    separate.extirpation=FALSE)
mod4 <- GeoHiSSE(phy, dat, f=sf, speciation=speciation
                 , extirpation=extirpation, hidden.areas=TRUE
                 , trans.rate=trans.rate, assume.cladogenetic=TRUE)


## Model 5. Heterogeneous diversification, not tied to range evolution. 
## Assumes 5 distinct diversification rates.
speciation <- c(rep(1,3), rep(2,3), rep(3,3), rep(4,3), rep(5,3)) 
extirpation <- c(rep(1,2), rep(2,2), rep(3,2), rep(4,2), rep(5,2)) 
trans.rate <- TransMatMakerGeoHiSSE(hidden.areas=4, make.null=TRUE, 
                                    include.jumps=FALSE, separate.extirpation=FALSE) 
mod5 <- GeoHiSSE(phy, dat, f=sf, speciation=speciation, extirpation=extirpation,
                 hidden.areas=TRUE, trans.rate=trans.rate, assume.cladogenetic=TRUE)



## Model 6. Heterogeneous diversification, not tied to range evolution. 
## Assumes two distinct diversification rates.
speciation <- c(1,1,1,2,2,2)
extirpation <- c(1,1,2,2)
trans.rate <- TransMatMakerGeoHiSSE(hidden.areas=1, include.jumps=FALSE , 
                                    separate.extirpation=FALSE)
mod6 <- GeoHiSSE(phy, dat, f=sf, speciation=speciation
                 , extirpation=extirpation, hidden.areas=TRUE
                 , trans.rate=trans.rate, assume.cladogenetic=TRUE)



###############################################################################
## Block of GeoSSE+extinction models.
## Here extirpation is NOT linked to range reduction.
## Range reduction is different from the extinction of an endemic lineage.
###############################################################################
## Model 7 - Dispersal parameters vary only, no range-dependent diversification. 
speciation <- c(1,1,1)
extirpation <- c(1,1)
trans.rate <- TransMatMakerGeoHiSSE(hidden.areas=0, include.jumps=FALSE
                                    , separate.extirpation=TRUE) 
mod7 <- GeoHiSSE(phy, dat, f=sf, speciation=speciation, extirpation=extirpation,
                 hidden.areas=FALSE , trans.rate=trans.rate, assume.cladogenetic=TRUE)



## Model 8. GeoSSE model, with range effect on diversification speciation <- c(1,2,3)
extirpation <- c(1,2)
trans.rate <- TransMatMakerGeoHiSSE(hidden.areas=0, include.jumps=FALSE
                                    , separate.extirpation=TRUE) 
mod8 <- GeoHiSSE(phy, dat, f=sf, speciation=speciation, extirpation=extirpation,
                 hidden.areas=FALSE, trans.rate=trans.rate, assume.cladogenetic=TRUE)


## Model 9. Heterogeneous diversification, not tied to range evolution.
## Assumes three distinct diversification rates.
speciation <- c(1,1,1,2,2,2,3,3,3)
extirpation <- c(1,1,2,2,3,3)
trans.rate <- TransMatMakerGeoHiSSE(hidden.areas=2, make.null=TRUE,include.jumps=FALSE,
                                    separate.extirpation=TRUE)
mod9 <- GeoHiSSE(phy, dat, f=sf, speciation=speciation, extirpation=extirpation,
                 hidden.areas=TRUE, trans.rate=trans.rate, assume.cladogenetic=TRUE)



## Model 10. Heterogeneous diversification, tied to range evolution.
## Assumes 6 distinct diversification rates.
speciation <- c(1,2,3,4,5,6)
extirpation <- c(1,2,3,4)
trans.rate <- TransMatMakerGeoHiSSE(hidden.areas=1, include.jumps=FALSE, 
                                    separate.extirpation=TRUE) 
mod10 <- GeoHiSSE(phy, dat, f=sf, speciation=speciation, extirpation=extirpation,
                  hidden.areas=TRUE, trans.rate=trans.rate, assume.cladogenetic=TRUE)



## Model 11. Heterogeneous diversification, not tied to range evolution. 
## Assumes 5 distinct diversification rates.
speciation <- c(rep(1,3), rep(2,3), rep(3,3), rep(4,3), rep(5,3)) 
extirpation <- c(rep(1,2), rep(2,2), rep(3,2), rep(4,2), rep(5,2)) 
trans.rate <- TransMatMakerGeoHiSSE(hidden.areas=4, make.null=TRUE, include.jumps=FALSE,
                                    separate.extirpation=TRUE) 
mod11 <- GeoHiSSE(phy, dat, f=sf, speciation=speciation, extirpation=extirpation,
                  hidden.areas=TRUE, trans.rate=trans.rate, assume.cladogenetic=TRUE)

#start 8:09, 8:17 (8 mins 15 params)
#1:23 start, 1:28

## Model 12. Heterogeneous diversification, not tied to range evolution. 
## Assumes two distinct diversification rates.
speciation <- c(1,1,1,2,2,2)
extirpation <- c(1,1,2,2)
trans.rate <- TransMatMakerGeoHiSSE(hidden.areas=1, include.jumps=FALSE , 
                                    separate.extirpation=TRUE)
mod12 <- GeoHiSSE(phy, dat, f=sf, speciation=speciation , extirpation=extirpation,
                  hidden.areas=TRUE, trans.rate=trans.rate, assume.cladogenetic=TRUE)


###############################################################################
## Block of anagenetic geographic models (MuSSE).
## Here models emulate GeoSSE (or GeoHiSSE) but changes only happen along branches.
###############################################################################
## Model 13. Transitions only. No character effect on diversification
speciation <- c(1,1,1)
extirpation <- c(1,1,1)
trans.rate <- TransMatMakerGeoHiSSE(hidden.areas=0, include.jumps=FALSE
                                    , separate.extirpation=TRUE) 
mod13 <- GeoHiSSE(phy, dat, f=sf, speciation=speciation, 
                  extirpation=extirpation, hidden.areas=FALSE, trans.rate=trans.rate, 
                  assume.cladogenetic=FALSE) 



## Model 14. Character effect on diversification.
speciation <- c(1,2,3)
extirpation <- c(1,2,3)
trans.rate <- TransMatMakerGeoHiSSE(hidden.areas=0, include.jumps=FALSE
                                    , separate.extirpation=TRUE) 
mod14 <- GeoHiSSE(phy, dat, f=sf, speciation=speciation , extirpation=extirpation,
                  hidden.areas=FALSE, trans.rate=trans.rate, assume.cladogenetic=FALSE)


## Model 15. No character effect on diversification.
speciation <- c(1,1,1,2,2,2,3,3,3)
extirpation <- c(1,1,1,2,2,2,3,3,3)
trans.rate <- TransMatMakerGeoHiSSE(hidden.areas=2, include.jumps=FALSE
                                    , separate.extirpation=TRUE, make.null=TRUE) 
mod15 <- GeoHiSSE(phy, dat, f=sf, speciation=speciation, extirpation=extirpation,
                  hidden.areas=TRUE, trans.rate=trans.rate, assume.cladogenetic=FALSE)


## Model 16. Character effect on diversification, with a hidden state
speciation <- c(1,2,3,4,5,6)
extirpation <- c(1,2,3,4,5,6)
trans.rate <- trans.rate <- TransMatMakerGeoHiSSE(hidden.areas=1, include.jumps=FALSE
                                                  , separate.extirpation=TRUE)
mod16 <- GeoHiSSE(phy, dat, f=sf, speciation=speciation, extirpation=extirpation,
                  hidden.areas=TRUE, trans.rate=trans.rate, assume.cladogenetic=FALSE)


## Model 17. No character effect on diversification, multiple shifts
speciation <- c(rep(1,3), rep(2,3), rep(3,3), rep(4,3), rep(5,3))
extirpation <- c(rep(1,3), rep(2,3), rep(3,3), rep(4,3), rep(5,3))
trans.rate <- TransMatMakerGeoHiSSE(hidden.areas=4, include.jumps=FALSE,
                                    separate.extirpation=TRUE, make.null=TRUE) 
mod17 <- GeoHiSSE(phy, dat, f=sf, speciation=speciation, extirpation=extirpation,
                  hidden.areas=TRUE, trans.rate=trans.rate, assume.cladogenetic=FALSE) 


## Model 18. No character effect on diversification, multiple shifts.
speciation <- c(rep(1,3), rep(2,3))
extirpation <- c(rep(1,3), rep(2,3))
trans.rate <- TransMatMakerGeoHiSSE(hidden.areas=1, include.jumps=FALSE, 
                                    separate.extirpation=TRUE, make.null=TRUE)
mod18 <- GeoHiSSE(phy, dat, f=sf, speciation=speciation, extirpation=extirpation
                  , hidden.areas=TRUE, trans.rate=trans.rate, assume.cladogenetic=FALSE)



#################### JUMP MODELS ##########################
#################### JUMP MODELS ##########################
#################### JUMP MODELS ##########################

## argument "include.jumps" set to TRUE

###############################################################################
## Block of GeoSSE-like models.
## Here extirpation is linked to range reduction.
###############################################################################
## Model 19 - Dispersal parameters vary only, no range-dependent diversification. 
speciation <- c(1,1,1)
extirpation <- c(1,1)
trans.rate <- TransMatMakerGeoHiSSE(hidden.areas=0, include.jumps=TRUE 
                                    , separate.extirpation=FALSE) 
mod19 <- GeoHiSSE(phy, dat, f=sf, speciation=speciation, extirpation=extirpation,
                  hidden.areas=FALSE, trans.rate=trans.rate, assume.cladogenetic=TRUE) 


## Model 20. Canonical GeoSSE model, range effect on diversification 
speciation <- c(1,2,3)
extirpation <- c(1,2)
trans.rate <- TransMatMakerGeoHiSSE(hidden.areas=0, include.jumps=TRUE
                                    , separate.extirpation=FALSE) 
mod20 <- GeoHiSSE(phy, dat, f=sf, speciation=speciation , extirpation=extirpation,
                  hidden.areas=FALSE , trans.rate=trans.rate, assume.cladogenetic=TRUE)


## Model 21. Heterogeneous diversification, not tied to range evolution.
## Assumes three distinct diversification rates.
## Dispersion parameters across hidden areas are the same.
speciation <- c(1,1,1,2,2,2,3,3,3)
extirpation <- c(1,1,2,2,3,3)
trans.rate <- TransMatMakerGeoHiSSE(hidden.areas=2, make.null=TRUE
                                    , include.jumps=TRUE, separate.extirpation=FALSE) 
mod21 <- GeoHiSSE(phy, dat, f=sf, speciation=speciation, extirpation=extirpation,
                  hidden.areas=TRUE, trans.rate=trans.rate, assume.cladogenetic=TRUE) 


## Model 22. Heterogeneous diversification, tied to range evolution. 
## Assumes 6 distinct diversification rates.
speciation <- c(1,2,3,4,5,6)
extirpation <- c(1,2,3,4)
trans.rate <- TransMatMakerGeoHiSSE(hidden.areas=1, include.jumps=TRUE , 
                                    separate.extirpation=FALSE)
mod22 <- GeoHiSSE(phy, dat, f=sf, speciation=speciation
                  , extirpation=extirpation, hidden.areas=TRUE
                  , trans.rate=trans.rate, assume.cladogenetic=TRUE)


## Model 23. Heterogeneous diversification, not tied to range evolution. 
## Assumes 5 distinct diversification rates.
speciation <- c(rep(1,3), rep(2,3), rep(3,3), rep(4,3), rep(5,3)) 
extirpation <- c(rep(1,2), rep(2,2), rep(3,2), rep(4,2), rep(5,2)) 
trans.rate <- TransMatMakerGeoHiSSE(hidden.areas=4, make.null=TRUE, 
                                    include.jumps=TRUE, separate.extirpation=FALSE) 
mod23 <- GeoHiSSE(phy, dat, f=sf, speciation=speciation, extirpation=extirpation,
                  hidden.areas=TRUE, trans.rate=trans.rate, assume.cladogenetic=TRUE)


## Model 24. Heterogeneous diversification, not tied to range evolution. 
## Assumes two distinct diversification rates.
speciation <- c(1,1,1,2,2,2)
extirpation <- c(1,1,2,2)
trans.rate <- TransMatMakerGeoHiSSE(hidden.areas=1, include.jumps=TRUE , 
                                    separate.extirpation=FALSE)
mod24 <- GeoHiSSE(phy, dat, f=sf, speciation=speciation
                  , extirpation=extirpation, hidden.areas=TRUE
                  , trans.rate=trans.rate, assume.cladogenetic=TRUE)

###############################################################################
## Block of GeoSSE+extinction models.
## Here extirpation is NOT linked to range reduction.
## Range reduction is different from the extinction of an endemic lineage.
###############################################################################
## Model 25 - Dispersal parameters vary only, no range-dependent diversification. 
speciation <- c(1,1,1)
extirpation <- c(1,1)
trans.rate <- TransMatMakerGeoHiSSE(hidden.areas=0, include.jumps=TRUE
                                    , separate.extirpation=TRUE) 
mod25 <- GeoHiSSE(phy, dat, f=sf, speciation=speciation, extirpation=extirpation,
                  hidden.areas=FALSE , trans.rate=trans.rate, assume.cladogenetic=TRUE)


## Model 26. GeoSSE model, with range effect on diversification speciation <- c(1,2,3)
extirpation <- c(1,2)
trans.rate <- TransMatMakerGeoHiSSE(hidden.areas=0, include.jumps=TRUE
                                    , separate.extirpation=TRUE) 
mod26 <- GeoHiSSE(phy, dat, f=sf, speciation=speciation, extirpation=extirpation,
                  hidden.areas=FALSE, trans.rate=trans.rate, assume.cladogenetic=TRUE)


## Model 27. Heterogeneous diversification, not tied to range evolution.
## Assumes three distinct diversification rates.
speciation <- c(1,1,1,2,2,2,3,3,3)
extirpation <- c(1,1,2,2,3,3)
trans.rate <- TransMatMakerGeoHiSSE(hidden.areas=2, make.null=TRUE,include.jumps=TRUE,
                                    separate.extirpation=TRUE)
mod27 <- GeoHiSSE(phy, dat, f=sf, speciation=speciation, extirpation=extirpation,
                  hidden.areas=TRUE, trans.rate=trans.rate, assume.cladogenetic=TRUE)


## Model 28. Heterogeneous diversification, tied to range evolution.
## Assumes 6 distinct diversification rates.
speciation <- c(1,2,3,4,5,6)
extirpation <- c(1,2,3,4)
trans.rate <- TransMatMakerGeoHiSSE(hidden.areas=1, include.jumps=TRUE, 
                                    separate.extirpation=TRUE) 
mod28 <- GeoHiSSE(phy, dat, f=sf, speciation=speciation, extirpation=extirpation,
                  hidden.areas=TRUE, trans.rate=trans.rate, assume.cladogenetic=TRUE)


## Model 29. Heterogeneous diversification, not tied to range evolution. 
## Assumes 5 distinct diversification rates.
speciation <- c(rep(1,3), rep(2,3), rep(3,3), rep(4,3), rep(5,3)) 
extirpation <- c(rep(1,2), rep(2,2), rep(3,2), rep(4,2), rep(5,2)) 
trans.rate <- TransMatMakerGeoHiSSE(hidden.areas=4, make.null=TRUE, include.jumps=TRUE,
                                    separate.extirpation=TRUE) 
mod29 <- GeoHiSSE(phy, dat, f=sf, speciation=speciation, extirpation=extirpation,
                  hidden.areas=TRUE, trans.rate=trans.rate, assume.cladogenetic=TRUE)


## Model 30. Heterogeneous diversification, not tied to range evolution. 
## Assumes two distinct diversification rates.
speciation <- c(1,1,1,2,2,2)
extirpation <- c(1,1,2,2)
trans.rate <- TransMatMakerGeoHiSSE(hidden.areas=1, include.jumps=TRUE , 
                                    separate.extirpation=TRUE)
mod30 <- GeoHiSSE(phy, dat, f=sf, speciation=speciation , extirpation=extirpation,
                  hidden.areas=TRUE, trans.rate=trans.rate, assume.cladogenetic=TRUE)



###############################################################################
## Block of anagenetic geographic models (MuSSE).
## Here models emulate GeoSSE (or GeoHiSSE) but changes only happen along branches.
###############################################################################
## Model 31. Transitions only. No character effect on diversification
speciation <- c(1,1,1)
extirpation <- c(1,1,1)
trans.rate <- TransMatMakerGeoHiSSE(hidden.areas=0, include.jumps=TRUE
                                    , separate.extirpation=TRUE) 
mod31 <- GeoHiSSE(phy, dat, f=sf, speciation=speciation, 
                  extirpation=extirpation, hidden.areas=FALSE, trans.rate=trans.rate, 
                  assume.cladogenetic=FALSE) 


## Model 32. Character effect on diversification.
speciation <- c(1,2,3)
extirpation <- c(1,2,3)
trans.rate <- TransMatMakerGeoHiSSE(hidden.areas=0, include.jumps=TRUE
                                    , separate.extirpation=TRUE) 
mod32 <- GeoHiSSE(phy, dat, f=sf, speciation=speciation , extirpation=extirpation,
                  hidden.areas=FALSE, trans.rate=trans.rate, assume.cladogenetic=FALSE)



## Model 33. No character effect on diversification.
speciation <- c(1,1,1,2,2,2,3,3,3)
extirpation <- c(1,1,1,2,2,2,3,3,3)
trans.rate <- TransMatMakerGeoHiSSE(hidden.areas=2, include.jumps=TRUE
                                    , separate.extirpation=TRUE, make.null=TRUE) 
mod33 <- GeoHiSSE(phy, dat, f=sf, speciation=speciation, extirpation=extirpation,
                  hidden.areas=TRUE, trans.rate=trans.rate, assume.cladogenetic=FALSE)


## Model 34. Character effect on diversification, with a hidden state
speciation <- c(1,2,3,4,5,6)
extirpation <- c(1,2,3,4,5,6)
trans.rate <- trans.rate <- TransMatMakerGeoHiSSE(hidden.areas=1, include.jumps=TRUE
                                                  , separate.extirpation=TRUE)
mod34 <- GeoHiSSE(phy, dat, f=sf, speciation=speciation, extirpation=extirpation,
                  hidden.areas=TRUE, trans.rate=trans.rate, assume.cladogenetic=FALSE)


## Model 35. No character effect on diversification, multiple shifts
speciation <- c(rep(1,3), rep(2,3), rep(3,3), rep(4,3), rep(5,3))
extirpation <- c(rep(1,3), rep(2,3), rep(3,3), rep(4,3), rep(5,3))
trans.rate <- TransMatMakerGeoHiSSE(hidden.areas=4, include.jumps=TRUE,
                                    separate.extirpation=TRUE, make.null=TRUE) 
mod35 <- GeoHiSSE(phy, dat, f=sf, speciation=speciation, extirpation=extirpation,
                  hidden.areas=TRUE, trans.rate=trans.rate, assume.cladogenetic=FALSE) 


## Model 36. No character effect on diversification, multiple shifts.
speciation <- c(rep(1,3), rep(2,3))
extirpation <- c(rep(1,3), rep(2,3))
trans.rate <- TransMatMakerGeoHiSSE(hidden.areas=1, include.jumps=TRUE, 
                                    separate.extirpation=TRUE, make.null=TRUE)
mod36 <- GeoHiSSE(phy, dat, f=sf, speciation=speciation, extirpation=extirpation
                  , hidden.areas=TRUE, trans.rate=trans.rate, assume.cladogenetic=FALSE)


save.image("CHAM.SW.36.Geohisse.Rdata")

##########################################################################################
############################## aics ######################################################
##########################################################################################

##### extracting aics #######
list.geohisse <- list(model1 = mod1, 
                     model7 = mod7, model13 = mod13, model19 = mod19) ###
     


### original
list.geohisse <- list(model1 = mod1, 
                      model2 = mod2, 
                      model3 = mod3, 
                      model4 = mod4, 
                      model5 = mod5,
                      model6 = mod6, 
                      model7 = mod7, 
                      model8 = mod8, 
                      model9 = mod9, 
                      model10 = mod10, 
                      model11 = mod11, 
                      model12 = mod12, 
                      model13 = mod13,
                      model14 = mod14, 
                      model15 = mod15, 
                      model16 = mod16, 
                      model17 = mod17,
                      model18 = mod18, 
                      model19 = mod19, 
                      model20 = mod20, 
                      model21 = mod21, 
                      model22 = mod22, 
                      model23 = mod23,
                      model24 = mod24,
                      model25 = mod25, 
                      model26 = mod26, 
                      model27 = mod27, 
                      model28 = mod28, 
                      model29 = mod29, 
                      model30 = mod30, 
                      model31 = mod31,
                      model32 = mod32, 
                      model33 = mod33, 
                      model34 = mod34, 
                      model35 = mod35,
                      model36 = mod36) 


saveRDS(list.geohisse, file=paste(group, "GeoHiSSE.models", ".RData", sep="_"))

##########################################################################################
############################## recons ####################################################
##########################################################################################

recon.mod1 <- MarginReconGeoSSE(phy = mod1$phy, data = mod1$data, f = mod1$f,
                                pars = mod1$solution, hidden.areas = mod1$hidden.areas,
                                root.type = mod1$root.type, root.p = mod1$root.p,
                                aic = mod1$AIC, n.cores = 4)
recon.mod2 <- MarginReconGeoSSE(phy = mod2$phy, data = mod2$data, f = mod2$f,
                                pars = mod2$solution, hidden.areas = mod2$hidden.areas, 
                                root.type = mod2$root.type, root.p = mod2$root.p,
                                aic = mod2$AIC, n.cores = 4)

#######################################
##original
recon.mod1 <- MarginReconGeoSSE(phy = mod1$phy, data = mod1$data, f = mod1$f,
                                pars = mod1$solution, hidden.areas = mod1$hidden.areas,
                                root.type = mod1$root.type, root.p = mod1$root.p,
                                aic = mod1$AIC, n.cores = 4)
recon.mod2 <- MarginReconGeoSSE(phy = mod2$phy, data = mod2$data, f = mod2$f,
                                pars = mod2$solution, hidden.areas = mod2$hidden.areas, 
                                root.type = mod2$root.type, root.p = mod2$root.p,
                                aic = mod2$AIC, n.cores = 4)
recon.mod3 <- MarginReconGeoSSE(phy = mod3$phy, data = mod3$data, f = mod3$f,
                                pars = mod3$solution, hidden.areas = mod3$hidden.areas,
                                root.type = mod3$root.type, root.p = mod3$root.p,
                                aic = mod3$AIC, n.cores = 4)
recon.mod4 <- MarginReconGeoSSE(phy = mod4$phy, data = mod4$data, f = mod4$f,
                                pars = mod4$solution, hidden.areas = mod4$hidden.areas, 
                                root.type = mod4$root.type, root.p = mod4$root.p,
                                aic = mod4$AIC, n.cores = 4)
recon.mod5 <- MarginReconGeoSSE(phy = mod5$phy, data = mod5$data, f = mod5$f,
                                pars = mod5$solution, hidden.areas = mod5$hidden.areas, 
                                root.type = mod5$root.type, root.p = mod5$root.p,
                                aic = mod5$AIC, n.cores = 4)
recon.mod6 <- MarginReconGeoSSE(phy = mod6$phy, data = mod6$data, f = mod6$f,
                                pars = mod6$solution, hidden.areas = mod6$hidden.areas, 
                                root.type = mod6$root.type, root.p = mod6$root.p,
                                aic = mod6$AIC, n.cores = 4)
recon.mod7 <- MarginReconGeoSSE(phy = mod7$phy, data = mod7$data, f = mod7$f,
                                pars = mod7$solution, hidden.areas = mod7$hidden.areas, 
                                root.type = mod7$root.type, root.p = mod7$root.p,
                                aic = mod7$AIC, n.cores = 4)
recon.mod8 <- MarginReconGeoSSE(phy = mod8$phy, data = mod8$data, f = mod8$f,
                                pars = mod8$solution, hidden.areas = mod8$hidden.areas, 
                                root.type = mod8$root.type, root.p = mod8$root.p,
                                aic = mod8$AIC, n.cores = 4)
recon.mod9 <- MarginReconGeoSSE(phy = mod9$phy, data = mod9$data, f = mod9$f,
                                pars = mod9$solution, hidden.areas = mod9$hidden.areas, 
                                root.type = mod9$root.type, root.p = mod9$root.p,
                                aic = mod9$AIC, n.cores = 4)
recon.mod10 <- MarginReconGeoSSE(phy = mod10$phy, data = mod10$data, f = mod10$f,
                                 pars = mod10$solution, hidden.areas = mod10$hidden.areas, 
                                 root.type = mod10$root.type, root.p = mod10$root.p,
                                 aic = mod10$AIC, n.cores = 4)
recon.mod11 <- MarginReconGeoSSE(phy = mod11$phy, data = mod11$data, f = mod11$f,
                                 pars = mod11$solution, hidden.areas = mod11$hidden.areas, 
                                 root.type = mod11$root.type, root.p = mod11$root.p,
                                 aic = mod11$AIC, n.cores = 4)
recon.mod12 <- MarginReconGeoSSE(phy = mod12$phy, data = mod12$data, f = mod12$f,
                                 pars = mod12$solution, hidden.areas = mod12$hidden.areas, 
                                 root.type = mod12$root.type, root.p = mod12$root.p,
                                 aic = mod12$AIC, n.cores = 4)
recon.mod13 <- MarginReconGeoSSE(phy = mod13$phy, data = mod13$data, f = mod13$f,
                                 pars = mod13$solution, hidden.areas = mod13$hidden.areas, 
                                 root.type = mod13$root.type, root.p = mod13$root.p,
                                 aic = mod13$AIC, n.cores = 4)
recon.mod14 <- MarginReconGeoSSE(phy = mod14$phy, data = mod14$data, f = mod14$f,
                                 pars = mod14$solution, hidden.areas = mod14$hidden.areas, 
                                 root.type = mod14$root.type, root.p = mod14$root.p,
                                 aic = mod14$AIC, n.cores = 4)
recon.mod15 <- MarginReconGeoSSE(phy = mod15$phy, data = mod15$data, f = mod15$f,
                                 pars = mod15$solution, hidden.areas = mod15$hidden.areas, 
                                 root.type = mod15$root.type, root.p = mod15$root.p,
                                 aic = mod15$AIC, n.cores = 4)
recon.mod16 <- MarginReconGeoSSE(phy = mod16$phy, data = mod16$data, f = mod16$f,
                                 pars = mod16$solution, hidden.areas = mod16$hidden.areas, 
                                 root.type = mod16$root.type, root.p = mod16$root.p,
                                 aic = mod16$AIC, n.cores = 4)
recon.mod17 <- MarginReconGeoSSE(phy = mod17$phy, data = mod17$data, f = mod17$f,
                                 pars = mod17$solution, hidden.areas = mod17$hidden.areas, 
                                 root.type = mod17$root.type, root.p = mod17$root.p,
                                 aic = mod17$AIC, n.cores = 4)
recon.mod18 <- MarginReconGeoSSE(phy = mod18$phy, data = mod18$data, f = mod18$f,
                                 pars = mod18$solution, hidden.areas = mod18$hidden.areas, 
                                 root.type = mod18$root.type, root.p = mod18$root.p,
                                 aic = mod18$AIC, n.cores = 4)
recon.mod19 <- MarginReconGeoSSE(phy = mod19$phy, data = mod19$data, f = mod19$f,
                                 pars = mod19$solution, hidden.areas = mod19$hidden.areas,
                                 root.type = mod19$root.type, root.p = mod19$root.p,
                                 aic = mod19$AIC, n.cores = 4)
recon.mod20 <- MarginReconGeoSSE(phy = mod20$phy, data = mod20$data, f = mod20$f,
                                 pars = mod20$solution, hidden.areas = mod20$hidden.areas, 
                                 root.type = mod20$root.type, root.p = mod20$root.p,
                                 aic = mod2$AIC, n.cores = 4)
recon.mod21 <- MarginReconGeoSSE(phy = mod21$phy, data = mod21$data, f = mod21$f,
                                 pars = mod21$solution, hidden.areas = mod21$hidden.areas,
                                 root.type = mod21$root.type, root.p = mod21$root.p,
                                 aic = mod21$AIC, n.cores = 4)
recon.mod22 <- MarginReconGeoSSE(phy = mod22$phy, data = mod22$data, f = mod22$f,
                                 pars = mod22$solution, hidden.areas = mod22$hidden.areas, 
                                 root.type = mod22$root.type, root.p = mod22$root.p,
                                 aic = mod22$AIC, n.cores = 4)
recon.mod23 <- MarginReconGeoSSE(phy = mod23$phy, data = mod23$data, f = mod23$f,
                                 pars = mod23$solution, hidden.areas = mod23$hidden.areas, 
                                 root.type = mod23$root.type, root.p = mod23$root.p,
                                 aic = mod23$AIC, n.cores = 4)
recon.mod24 <- MarginReconGeoSSE(phy = mod24$phy, data = mod24$data, f = mod24$f,
                                 pars = mod24$solution, hidden.areas = mod24$hidden.areas, 
                                 root.type = mod24$root.type, root.p = mod24$root.p,
                                 aic = mod24$AIC, n.cores = 4)
recon.mod25 <- MarginReconGeoSSE(phy = mod25$phy, data = mod25$data, f = mod25$f,
                                 pars = mod25$solution, hidden.areas = mod25$hidden.areas, 
                                 root.type = mod25$root.type, root.p = mod25$root.p,
                                 aic = mod25$AIC, n.cores = 4)
recon.mod26 <- MarginReconGeoSSE(phy = mod26$phy, data = mod26$data, f = mod26$f,
                                 pars = mod26$solution, hidden.areas = mod26$hidden.areas, 
                                 root.type = mod26$root.type, root.p = mod26$root.p,
                                 aic = mod26$AIC, n.cores = 4)
recon.mod27 <- MarginReconGeoSSE(phy = mod27$phy, data = mod27$data, f = mod27$f,
                                 pars = mod27$solution, hidden.areas = mod27$hidden.areas, 
                                 root.type = mod27$root.type, root.p = mod27$root.p,
                                 aic = mod27$AIC, n.cores = 4)
recon.mod28 <- MarginReconGeoSSE(phy = mod28$phy, data = mod28$data, f = mod28$f,
                                 pars = mod28$solution, hidden.areas = mod28$hidden.areas, 
                                 root.type = mod28$root.type, root.p = mod28$root.p,
                                 aic = mod28$AIC, n.cores = 4)
recon.mod29 <- MarginReconGeoSSE(phy = mod29$phy, data = mod29$data, f = mod29$f,
                                 pars = mod29$solution, hidden.areas = mod29$hidden.areas, 
                                 root.type = mod29$root.type, root.p = mod29$root.p,
                                 aic = mod29$AIC, n.cores = 4)
recon.mod30 <- MarginReconGeoSSE(phy = mod30$phy, data = mod30$data, f = mod30$f,
                                 pars = mod30$solution, hidden.areas = mod30$hidden.areas, 
                                 root.type = mod30$root.type, root.p = mod30$root.p,
                                 aic = mod30$AIC, n.cores = 4)
recon.mod31 <- MarginReconGeoSSE(phy = mod31$phy, data = mod31$data, f = mod31$f,
                                 pars = mod31$solution, hidden.areas = mod31$hidden.areas, 
                                 root.type = mod31$root.type, root.p = mod31$root.p,
                                 aic = mod31$AIC, n.cores = 4)
recon.mod32 <- MarginReconGeoSSE(phy = mod32$phy, data = mod32$data, f = mod32$f,
                                 pars = mod32$solution, hidden.areas = mod32$hidden.areas, 
                                 root.type = mod32$root.type, root.p = mod32$root.p,
                                 aic = mod32$AIC, n.cores = 4)
recon.mod33 <- MarginReconGeoSSE(phy = mod33$phy, data = mod33$data, f = mod33$f,
                                 pars = mod33$solution, hidden.areas = mod33$hidden.areas, 
                                 root.type = mod33$root.type, root.p = mod33$root.p,
                                 aic = mod33$AIC, n.cores = 4)
recon.mod34 <- MarginReconGeoSSE(phy = mod34$phy, data = mod34$data, f = mod34$f,
                                 pars = mod34$solution, hidden.areas = mod34$hidden.areas, 
                                 root.type = mod34$root.type, root.p = mod34$root.p,
                                 aic = mod34$AIC, n.cores = 4)
recon.mod35 <- MarginReconGeoSSE(phy = mod35$phy, data = mod35$data, f = mod35$f,
                                 pars = mod35$solution, hidden.areas = mod35$hidden.areas, 
                                 root.type = mod35$root.type, root.p = mod35$root.p,
                                 aic = mod35$AIC, n.cores = 4)
recon.mod36 <- MarginReconGeoSSE(phy = mod36$phy, data = mod36$data, f = mod36$f,
                                 pars = mod36$solution, hidden.areas = mod36$hidden.areas, 
                                 root.type = mod36$root.type, root.p = mod36$root.p,
                                 aic = mod36$AIC, n.cores = 4)


#####################################
#############
recon.models <- list(recon.mod1, recon.mod2, recon.mod3, recon.mod4, 
                     recon.mod5, recon.mod6, recon.mod7, recon.mod8, 
                     recon.mod9, recon.mod10, recon.mod11, recon.mod12, 
                     recon.mod13, recon.mod14, recon.mod15, recon.mod16, 
                     recon.mod17, recon.mod18, recon.mod19, recon.mod20,
                     recon.mod21, recon.mod22, recon.mod23, recon.mod24,
                     recon.mod25, recon.mod26, recon.mod27, recon.mod28, 
                     recon.mod29, recon.mod30, recon.mod31, recon.mod32, 
                     recon.mod33, recon.mod34, recon.mod35, recon.mod36)

######

saveRDS(recon.models, file=paste("SW.GeoHiSSE_recon.models", group ,".RData", sep="_"))

recon.models <- load("Xylopia_GeoHiSSE.models_.RData")

###own additions

head(recon.mod4)

pdf("Xylopia.GeoH.mod4.netdiv.pdf", width=12, height=20)
plot.geohisse.states(x = recon.models, rate.param = "net.div", type = "phylogram", 
                     show.tip.label = TRUE, legend =TRUE)
dev.off()  # Turn off PDF


############################################################################################################
##########################################################################################

##GetAICWeights(recon.models, criterion="aic")

GetAICWeights(list(model1z = recon.mod1, model2z = recon.mod2, model3z = recon.mod3, model4z = recon.mod4, model5z = recon.mod5, model6z = recon.mod6, model7z = recon.mod7, model8z = recon.mod8, model9z = recon.mod9, model10z=recon.mod10, model11z=recon.mod11, model12z=recon.mod12, model13z=recon.mod13, model14z=recon.mod14, model15z=recon.mod15, model16z=recon.mod16, model17z=recon.mod17, model18z=recon.mod18, model19z=recon.mod19, model20z=recon.mod20, model21z=recon.mod21, model22z=recon.mod22, model23z=recon.mod23, model24z=recon.mod24, model25z=recon.mod25, model26z=recon.mod26, model27z=recon.mod27, model28z=recon.mod28, model29z=recon.mod29, model30z=recon.mod30, model31z=recon.mod31, model32z=recon.mod32, model33z=recon.mod33, model34z=recon.mod34, model35z=recon.mod35, model36z=recon.mod36), criterion="aic")

write.table((GetAICWeights(list(model1z = recon.mod1, model2z = recon.mod2, model3z = recon.mod3, model4z = recon.mod4, model5z = recon.mod5, model6z = recon.mod6, model7z = recon.mod7, model8z = recon.mod8, model9z = recon.mod9, model10z=recon.mod10, model11z=recon.mod11, model12z=recon.mod12, model13z=recon.mod13, model14z=recon.mod14, model15z=recon.mod15, model16z=recon.mod16, model17z=recon.mod17, model18z=recon.mod18, model19z=recon.mod19, model20z=recon.mod20, model21z=recon.mod21, model22z=recon.mod22, model23z=recon.mod23, model24z=recon.mod24, model25z=recon.mod25, model26z=recon.mod26, model27z=recon.mod27, model28z=recon.mod28, model29z=recon.mod29, model30z=recon.mod30, model31z=recon.mod31, model32z=recon.mod32, model33z=recon.mod33, model34z=recon.mod34, model35z=recon.mod35, model36z=recon.mod36), criterion="aic")), 
            file="AICweights.txt")


#### get rates for best model
model.ave.rates <- GetModelAveRates(x = recon.mod14, type = "tips")

head(model.ave.rates)

write.table(model.ave.rates, file='rates.txt')



