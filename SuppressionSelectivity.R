
# About this code ---- 

# This code was written primarily by Rae Fadlovich with sections (noted in text) adapted from code written by Timothy E. Walsworth 

----
  # This annotated code accompanies the article "Selectivity of invasive species suppression efforts influences control efficacy" in the Journal of Applied Ecology accepted for publication 11/2024

  # This code can be used to evaluate the impact of selectivity and effort on invasive species control efficacy 
  # As written, it provides projected population estimates for the invasive common carp population in Utah Lake, but we have heavily annotated it so that it can be adapted to other contexts  

  ----
  
# What data is needed to run this code ---- 

  # Population estimates
    # Our abundance-at-age estimates come from Walsworth, Wallace, et al., 2023 and the full model is described in Walsworth et al., 2020
      # To adapt this model for your own system, you will need abundance estimates for different age-, stage-, or size-classes of your target species 
    # Population estimates file "inputs/UTLCarp_Natage_Medians.csv"

  # Selectivity estimates 
    # We provide a JAGS script for assessing selectivity if you have population estimates and age-based survey data 
    # If you already have selectivity and population estimates for your population, you can skip ahead to the simulation section
    # We calculate selectivity from catch density data collected during standardized annual monitoring in 2020 and 2021 
    # Catch density files "inputs/UTL_CarpSeine_DensityAge_byHaul_2021.csv" "inputs/UTL_CarpSeine_DensityAge_byHaul_2020.csv"

  # Data for this project has been archived on Dryad and can alternatively be found in the rfadlovi/SuppressionSelectivity GitHub repository.


# Requires Packages and Programs 
  
  # Selectivity estimates
library(R2jags)
library(coda)
library(lattice)
library(MCMCvis)
library(DescTools)
library(bayestestR)

  # Simulation model 
library(MQMF)

  # Plotting the manuscript figures 
library(rcartocolor)
library(paletteer) 


# Selectivity Estimates -------------------------------------------------------- 


# Set working directory to Github folder 
setwd("~/Github/SuppressionSelectivity")

# Commented out as they show up above. Legacy from GitHub where these codes are separate files. 
# library(R2jags)
# library(coda)
# library(lattice)
# library(MCMCvis)
# library(DescTools)
# library(bayestestR)

### Calculate observed selectivity (Or skip ahead if already calc.)---------------------------------------------

# Read in the density-by-age by haul data for 2020 
dat20<-read.csv("inputs/UTL_CarpSeine_DensityAge_byHaul_2020.csv",header=T)

# Read in the density-by-age by haul data for 2021
dat21<-read.csv("inputs/UTL_CarpSeine_DensityAge_byHaul_2021.csv",header=T)

# Read in the abundance estimates 
Nage<-read.csv("inputs/UTLCarp_Natage_Medians.csv") # From 2022 model run - medians 

# Catchability coefficient 
q<-0.000334 # q from 2022 model run

# Function to calculate selectivity at age by gear  
sag<-function(catchat, Nat, q){
  # give some storage 
  s<-rep(NA, length(catchat))
  # loop through every row in csv 
  for (i in 1:length(catchat)){
    # main function for selectivity
    s[i]<-(log(1-(catchat[i]/Nat[i])))/(-q) 
  }
  # return the selectivity and name it
  return("selectivity"=s)
  
}

#### Calculate observed selectivity for 2020 -----------------------------------

# Carp density  
c.at20<-dat20$density

# Pull out just 2020 median abundance 
databund20<-Nage$X2020

# Adding Nage to main df - median at age
dat20$N<-Nage[dat20$age+1,13]

# Add in a seine haul number column rep
# Note that the 49 is hard coded. 
haul20<-rep(1:49,each=8)
dat20$Haul<-haul20


# Redundant but I want to have data org and code for function seperate so its double the fun 
N.at20<-dat20$N

# Run selectivity function with 2020 data 
sag20<-sag(catchat=c.at20,Nat=dat20$N,q=q)

# Add this selectivity into the main data 
dat20$Selectivity<-sag20

# Pull out the max selectivity for each haul 
maxselecs20<-aggregate(dat20$Selectivity,by=list(dat20$Haul),max)

# Make this vector long so it can be added onto dat 
ssvec20<-rep(maxselecs20$x,each=8)
# Add it on to dat
dat20$MaxSelectivity<-ssvec20

# Scale the selectivity values by max 
scale.selec20<-dat20$Selectivity/dat20$MaxSelectivity
# Replace NAs that resulted from 0/
scale.selec20[is.na(scale.selec20)]<- 0

# Add to data frame 
dat20$ScaleSelec<-scale.selec20

# Extract the large mesh stuff 
lmesh20<-dat20[dat20$mesh == "largemesh",]

# Extract the small mesh stuff 
smesh20<-dat20[dat20$mesh == "smallmesh",]

# Change large mesh values that will cause error in the JAGS model
lmesh20$ScaleSelec[lmesh20$ScaleSelec==0]<- .001
lmesh20$ScaleSelec[lmesh20$ScaleSelec==1]<- .999
lmesh_obs_20<-lmesh20$ScaleSelec
# Write out csv file to save data for later 
# write.csv(lmesh_obs_20, file="outputs/lmesh_scaled_selectivity_2020.csv")

# Change small mesh values that will cause error in the JAGS model
smesh20$ScaleSelec[smesh20$ScaleSelec==0]<- .001
smesh20$ScaleSelec[smesh20$ScaleSelec==1]<- .999
smesh_obs_20<-smesh20$ScaleSelec
# Write out csv file to save data for later 
# write.csv(smesh_obs_20, file="outputs/smesh_scaled_selectivity_2020.csv")



#### Calculate observed selectivity for 2021 -----------------------------------

# Carp density  
c.at21<-dat21$density

# Pull out just 2020 median abundance 
databund21<-Nage$X2021

# Adding Nage to main df - median number at age
dat21$N<-Nage[dat21$age+1,14]

# Add in a seine haul number column rep
# Note that the 45 is hard coded and will need to be changed if you have a different number of hauls.
haul21<-rep(1:45,each=8)

dat21$Haul<-haul21


# redundant but I want to have data org and code for function seperate so its double the fun 
N.at21<-dat21$N

# Run the selectivity function with 2021 data
sag21<-sag(catchat=c.at21,Nat=dat21$N,q=q)

# Add this selectivity into the main data 
dat21$Selectivity<-sag21

# Pull out the max selectivty
maxselecs21<-aggregate(dat21$Selectivity,by=list(dat21$Haul),max)

# Make this vector long so it can be added onto dat 
ssvec21<-rep(maxselecs21$x,each=8)

dat21$MaxSelectivity<-ssvec21

# Scale the selectivity values by max 
scale.selec21<-dat21$Selectivity/dat21$MaxSelectivity

scale.selec21[is.na(scale.selec21)]<- 0

# Add to data frame 
dat21$ScaleSelec<-scale.selec21

# Extract the large mesh stuff 
lmesh21<-dat21[dat21$mesh == "largemesh",]

# Extract the small mesh stuff 
smesh21<-dat21[dat21$mesh == "smallmesh",]

# Change large mesh values that will cause error in the JAGS model
lmesh21$ScaleSelec[lmesh21$ScaleSelec==0]<- .001
lmesh21$ScaleSelec[lmesh21$ScaleSelec==1]<- .999
lmesh_obs_21<-lmesh21$ScaleSelec
# Write out csv file to save data for later 
# write.csv(lmesh_obs_21, file="outputs/lmesh_scaled_selectivity_2021.csv")

# Change small mesh values that will cause error in the JAGS model
smesh21$ScaleSelec[smesh21$ScaleSelec==0]<- .001
smesh21$ScaleSelec[smesh21$ScaleSelec==1]<- .999
smesh_obs_21<-smesh21$ScaleSelec
# Write out csv file to save data for later 
# write.csv(smesh_obs_21, file="outputs/smesh_scaled_selectivity_2021.csv")


#### Organize observed selectivity data  ----------------------------------------------------

# If you have already calculated observed (scaled) selectivity and saved your outputs, you can start here 

# Read in the large mesh observations from 2020 - not necessary if running straight through 
# lmesh_obs_20<-read.csv("outputs/lmesh_scaled_selectivity_2020.csv")

# If you read in the file - pull out the vector 
# lmesh_obs_20<-lmesh_obs_20$x

# Read in the small mesh observations from 2020 - not necessary if running straight through 
# smesh_obs_20<-read.csv("outputs/smesh_scaled_selectivity_2020.csv") 

# If you read in the file - pull out the vector 
# smesh_obs_20<-smesh_obs_20$x

# Read in the large mesh observations from 2021 - not necessary if running straight through 
# lmesh_obs_21<-read.csv("lmesh_scaled_selectivity_2021.csv")
# If you read in the file - pull out the vector 
# lmesh_obs_21<-lmesh_obs_21$x


# Read in the small mesh observations from 2021
# smesh_obs_21<-read.csv("smesh_scaled_selectivity_2021.csv")
# We just want the observations 
# smesh_obs_21<-smesh_obs_21$x


# Combine 2020 and 2021 data into one file 

# For the large mesh 
lmesh_obs_20_21<-c(lmesh_obs_20, lmesh_obs_21) 

# For the small mesh 
smesh_obs_20_21<-c(smesh_obs_20, smesh_obs_21)


# Assign ages for each value in vector 

# Ages for the large 
age.l.20.21<-rep(0:7, length.out=length(lmesh_obs_20_21))
# Ages for the small
age.s.20.21<-rep(0:7, length.out=length(smesh_obs_20_21))


### Set up the JAGS model -------------------------------------------------------

# Set up the JAGS model to run with the large mesh data 

mod.data.l.20.21<- list("obs"=lmesh_obs_20_21,"age"=age.l.20.21+1, "nobs"=length(age.l.20.21),
                        "nages"=8)

# Set up the JAGS model to run with the small mesh data 

mod.data.s.20.21<- list("obs"=smesh_obs_20_21,"age"=age.s.20.21+1, "nobs"=length(age.s.20.21),
                        "nages"=8)

# Model parameters - are the same for small and large 
mod.params<- c("s50","s95","predS","vari" ,"alph","beta")

# JAGS model is in the bugs file 
mod.loc<-"selectivity.bug"

# Number of chains
mc.chains<-3

# Samples to pull 
mc.pull<-500000 # You can reduce this number to speed things up for trial runs where you test other variables

# Thinning rate
mc.thin<-500

# Run the models 
#! Please note that running these models can take a few hours on a 'normal' work machine. 
# Run the Large
mod.out.l.20.21<-jags(data=mod.data.l.20.21,inits=NULL,parameters.to.save=mod.params,
                      model.file=mod.loc,n.chains=mc.chains, n.iter=mc.pull,
                      n.thin=mc.thin)

# Save the run if you want to use the outputs later
# save(mod.out.l.20.21, file="outputs/mo.l.20.21.rdata")

# Run the small 
mod.out.s.20.21<-jags(data=mod.data.s.20.21,inits=NULL,parameters.to.save=mod.params,
                      model.file=mod.loc,n.chains=mc.chains, n.iter=mc.pull,
                      n.thin=mc.thin)

# Save the small run if you want to use the outputs later
# save(mod.out.s.20.21, file="outputs/mo.s.20.21.rdata")

### Summarize Outputs -----------------------------------------------------------

# Read in the model outputs if you have already run and saved them 
# Load the large
# load("outputs/mo.l.20.21.rdata")
# # Load the small
# load("outputs/mo.s.20.21.rdata")

# Use the coda package to create a MCMC object for large mesh
mod.out.fit.l<-as.mcmc(mod.out.l.20.21)

# Obtain model summary 
lsum<-summary(mod.out.fit.l)

# Extract uantiles 
lsumq<-lsum$quantiles

# Use the coda package to create a MCMC object for small mesh
mod.out.fit.s<-as.mcmc(mod.out.s.20.21)

# Obtain model summary 
ssum<-summary(mod.out.fit.s)

# Extract uantiles 
ssumq<-ssum$quantiles


# Function for predicted selectivity
pred.selec<- function(s50,s95,age){
  s<-rep(1,length(age))

  for (i in 1:length(age)){
    s[i]<- 1 / (1+exp(-log(19)*((age[i]-(s50-1))/((s95-1)-(s50-1)))))
  }
  return(s)
}


# Set up age vector 
age200<-seq(0,7,length.out = 200)

# This is for the median large mesh s50 s95 values 
l_med_pred<-pred.selec(age = age200, s50 = lsumq[26,3], 
                       s95 = lsumq[27,3])

# This is for the lower ci small mesh 
l_lci_pred<-pred.selec(age = age200, s50 = lsumq[26,1], 
                       s95 = lsumq[27,1])

# This is for the upper ci large mesh 
l_uci_pred<-pred.selec(age = age200, s50 = lsumq[26,5], 
                       s95 = lsumq[27,5])


# Calculate the betas # Previous attempts to turn this into a loop or function have failed - may update later to make more digestible 


# Get the simulation list from the large model output
mo.list.l<-mod.out.l.20.21$BUGSoutput$sims.list


# age0
age0l<-rbeta(100,mo.list.l$alph[1,1],mo.list.l$beta[1,1])  
#grab values from each sample and c 
for(i in 2:nrow(mo.list.l$alph)){
  age0l<-c(age0l,rbeta(100,mo.list.l$alph[i,1],mo.list.l$beta[i,1]))
}
# 1
age1l<-rbeta(100,mo.list.l$alph[1,2],mo.list.l$beta[1,2])  
#grab values from each sample and c 
for(i in 2:nrow(mo.list.l$alph)){
  age1l<-c(age1l,rbeta(100,mo.list.l$alph[i,2],mo.list.l$beta[i,2]))
}
# 2
age2l<-rbeta(100,mo.list.l$alph[1,3],mo.list.l$beta[1,3])  
#grab values from each sample and c 
for(i in 2:nrow(mo.list.l$alph)){
  age2l<-c(age2l,rbeta(100,mo.list.l$alph[i,3],mo.list.l$beta[i,3]))
}
# 3
age3l<-rbeta(100,mo.list.l$alph[1,4],mo.list.l$beta[1,4])  
#grab values from each sample and c 
for(i in 2:nrow(mo.list.l$alph)){
  age3l<-c(age3l,rbeta(100,mo.list.l$alph[i,4],mo.list.l$beta[i,4]))
}
# 4
age4l<-rbeta(100,mo.list.l$alph[1,5],mo.list.l$beta[1,5])  
#grab values from each sample and c 
for(i in 2:nrow(mo.list.l$alph)){
  age4l<-c(age4l,rbeta(100,mo.list.l$alph[i,5],mo.list.l$beta[i,5]))
}
# 5
age5l<-rbeta(100,mo.list.l$alph[1,6],mo.list.l$beta[1,6])  
#grab values from each sample and c 
for(i in 2:nrow(mo.list.l$alph)){
  age5l<-c(age5l,rbeta(100,mo.list.l$alph[i,6],mo.list.l$beta[i,6]))
}
# 6
age6l<-rbeta(100,mo.list.l$alph[1,7],mo.list.l$beta[1,7])  
#grab values from each sample and c 
for(i in 2:nrow(mo.list.l$alph)){
  age6l<-c(age6l,rbeta(100,mo.list.l$alph[i,7],mo.list.l$beta[i,7]))
}
#7
age7l<-rbeta(100,mo.list.l$alph[1,8],mo.list.l$beta[1,8])  
#grab values from each sample and c 
for(i in 2:nrow(mo.list.l$alph)){
  age7l<-c(age7l,rbeta(100,mo.list.l$alph[i,8],mo.list.l$beta[i,8]))
}


# Now get all the betas for all the smalls 

# Get the simulation list from the small model output
mo.list.s<-mod.out.s.20.21$BUGSoutput$sims.list

# Age 0 
age0s<-rbeta(100,mo.list.s$alph[1,1],mo.list.s$beta[1,1])  
#grab values from each sample and c 
for(i in 2:nrow(mo.list.s$alph)){
  age0s<-c(age0s,rbeta(100,mo.list.s$alph[i,1],mo.list.s$beta[i,1]))
}
# 1
age1s<-rbeta(100,mo.list.s$alph[1,2],mo.list.s$beta[1,2])  
#grab values from each sample and c 
for(i in 2:nrow(mo.list.s$alph)){
  age1s<-c(age1s,rbeta(100,mo.list.s$alph[i,2],mo.list.s$beta[i,2]))
}
# 2
age2s<-rbeta(100,mo.list.s$alph[1,3],mo.list.s$beta[1,3])  
#grab values from each sample and c 
for(i in 2:nrow(mo.list.s$alph)){
  age2s<-c(age2s,rbeta(100,mo.list.s$alph[i,3],mo.list.s$beta[i,3]))
}
# 3
age3s<-rbeta(100,mo.list.s$alph[1,4],mo.list.s$beta[1,4])  
#grab values from each sample and c 
for(i in 2:nrow(mo.list.s$alph)){
  age3s<-c(age3s,rbeta(100,mo.list.s$alph[i,4],mo.list.s$beta[i,4]))
}
# 4
age4s<-rbeta(100,mo.list.s$alph[1,5],mo.list.s$beta[1,5])  
#grab values from each sample and c 
for(i in 2:nrow(mo.list.s$alph)){
  age4s<-c(age4s,rbeta(100,mo.list.s$alph[i,5],mo.list.s$beta[i,5]))
}
# 5
age5s<-rbeta(100,mo.list.s$alph[1,6],mo.list.s$beta[1,6])  
#grab values from each sample and c 
for(i in 2:nrow(mo.list.s$alph)){
  age5s<-c(age5s,rbeta(100,mo.list.s$alph[i,6],mo.list.s$beta[i,6]))
}
# 6
age6s<-rbeta(100,mo.list.s$alph[1,7],mo.list.s$beta[1,7])  
#grab values from each sample and c 
for(i in 2:nrow(mo.list.s$alph)){
  age6s<-c(age6s,rbeta(100,mo.list.s$alph[i,7],mo.list.s$beta[i,7]))
}
# 7
age7s<-rbeta(100,mo.list.s$alph[1,8],mo.list.s$beta[1,8])  
#grab values from each sample and c 
for(i in 2:nrow(mo.list.s$alph)){
  age7s<-c(age7s,rbeta(100,mo.list.s$alph[i,8],mo.list.s$beta[i,8]))
}



# Calculate the CI - can change method if desired 
lci0h=ci(age0l,method = 'HDI') 
lci1h=ci(age1l,method = 'HDI') 
lci2h=ci(age2l,method = 'HDI') 
lci3h=ci(age3l,method = 'HDI')
lci4h=ci(age4l,method = 'HDI')
lci5h=ci(age5l,method = 'HDI')
lci6h=ci(age6l,method = 'HDI')
lci7h=ci(age7l,method = 'HDI')

# lower CI 
l_lci_posth<-c(lci0h$CI_low,lci1h$CI_low,lci2h$CI_low,lci3h$CI_low,lci4h$CI_low,
               lci5h$CI_low,lci6h$CI_low,lci7h$CI_low)
# upper CI 
l_uci_posth<-c(lci0h$CI_high,lci1h$CI_high,lci2h$CI_high,lci3h$CI_high,lci4h$CI_high,
               lci5h$CI_high,lci6h$CI_high,lci7h$CI_high)


# CI for the small mesh - can change method if desired

sci0h=ci(age0s,method = 'HDI') 
sci1h=ci(age1s,method = 'HDI') 
sci2h=ci(age2s,method = 'HDI') 
sci3h=ci(age3s,method = 'HDI')
sci4h=ci(age4s,method = 'HDI')
sci5h=ci(age5s,method = 'HDI')
sci6h=ci(age6s,method = 'HDI')
sci7h=ci(age7s,method = 'HDI')

# Combine the lower CI 
s_lci_posth<-c(sci0h$CI_low,sci1h$CI_low,sci2h$CI_low,sci3h$CI_low,sci4h$CI_low,
               sci5h$CI_low,sci6h$CI_low,sci7h$CI_low)
# Upper CI 
s_uci_posth<-c(sci0h$CI_high,sci1h$CI_high,sci2h$CI_high,sci3h$CI_high,sci4h$CI_high,
               sci5h$CI_high,sci6h$CI_high,sci7h$CI_high)

# Get the pred posteriors 
s_lci_pred<-pred.selec(age = age200, s50 = ssumq[26,1], 
                       s95 = ssumq[27,1])

# Upper CI for the Small mesh 
s_uci_pred<-pred.selec(age = age200, s50 = ssumq[26,5], 
                       s95 = ssumq[27,5])

# Median predicted selectivity
s_med_pred<-pred.selec(age = age200, s50 = ssumq[26,3], 
                       s95 = ssumq[27,3])



# These are the outputs that I read in for the plotting section 
# If you are running this straight through, you do not necessarily need to save. But you may want to do so as it allows you to clear the rest of your working directory

save(l_lci_posth, l_uci_posth, s_lci_posth, s_uci_posth,
     age200, s_lci_pred, s_uci_pred, l_lci_pred, l_uci_pred,
     l_med_pred, s_med_pred, file = "outputs/largesmall_jags_plotinfo2.RData")



# Simulations ------------------------------------------------------------------

  # Running simulations of projected carp population with different selectivity and effort 

  # This code also runs simulations for knife-edge and domed selectivity scenarios. 
  # While we do not include them in our manuscript as they were distracting to our overall goals, we have included them here as they may be helpful for practicioners with these selectivity shapes in their removal programs


# Rae Fadlovich modification of Timothy E. Walsworth code 

###load packages#### 
library(MQMF)

#..........................................
### Load in outputs #################
#..........................................

# Load in the JAGS output - only necessary if you are not running straight through. But it is not a bad idea to clear and refresh environment here 


# # large
#  load("outputs/mo.l.20.21.rdata")

# Load population model estimates

## Our model parameters come from the 2023 carp model  (see Walsworth, Wallace, et al., 2023 for a model description)
 load("inputs/UTL_VaryQ_boundSR22023-01-06.RData")

# Remove unneeded things to clean up space 
rm(abund.table, abund.table.dn, abund.table.up,abundmc, abundquants,
   adjust.abundtable, biomass, biomassH, biomassL, biomassmc, 
   biomassquants, biomassspn, cis, dens.out, dens.tab,
   estdensmc, harvmc, mod.out, moddat, Nat, Nat.cv, Nat.se, 
   Natmcmd, Natuse, subdat)


##................................................
### Set up selectivity and variance arrays ####
##................................................

# From the JAGS model 

# Pull out the mean selectivity 
predS_mean<-mod.out.l.20.21$BUGSoutput$mean$predS 
# Pull the mean variance
predS_vari<-mod.out.l.20.21$BUGSoutput$mean$vari 


# Calculate selectivity by shifting the mean S50 and s95 by 1 year

# mean s50 from the JAGS model
s501<-mod.out.l.20.21$BUGSoutput$mean$s50 
# mean s95 from the JAGS model
s951<-mod.out.l.20.21$BUGSoutput$mean$s95 

# Find the s50 and s95 CI  
# mod.out.l.20.21$BUGSoutput$summary


# Calculate simulated log selectivities 

# Shift the s50 by 1 again and again to get a shifted vector
s50_shift<-(rep(s501,11)-c(0:10)) 
# Shifted vector for s95
s95_shift<-(rep(s951,11)-c(0:10)) 


# Function to calc selectivity off s50 and s95
pred.selec.s5095c.loop<- function(s50,s95,age){ 
  s<-matrix(NA, nrow = length(age), ncol = length(s50))
  for (j in seq_along(s50)){
    for (i in seq_along(age)){
      s[i,j]<- 1 / (1+exp(-log(19)*((age[i]-(s50[j]))/((s95[j])-(s50[j])))))
    }
  }
  
  (s)
}

# Calculate shifted selectivity values for all ages
s_log<-pred.selec.s5095c.loop(s50_shift,s95_shift, c(1:nage)) 

# If you decide to add more shapes/ scenarios onto these, you can just loop starting with new selectivity scenario index :)  



#...........................................................#
#calculate domed selectivity - needs MQMF - HADDON 2011 
#...........................................................#

# Please note that we do not include domed selectivity in our manuscript, but have left it in the code in case this shape is of use to you 

a<-c(1:8) # ages
ag<-c(2:7) 
p1<-(ag) # look at the documentation for these 6 variables. It is a segmented formula 
p2<-ag
p3<-rep(8,6)
p4<-rep(8,6)
p5<-rep(-10,6)
p6<-rep(-10,6)

# Create domed selectivity storage 
domeselec<-matrix(rep(NA), ncol = length(p1), nrow = 8)

# Loop through different values in the vectors to get multiple "shifted" domes 
for (i in 1:length(ag)){
  p<-c(p1[i],p2[i],p3[i],p4[i],p5[i],p6[i])
  domeselec[,i]<-domed(p,a)
}


#...........................................................#
#calculate knife edge selectivity - needs MQMF - HADDON 2011 
#...........................................................#

# Please note that we do not include knife edge selectivity in our manuscript, but have left it in the code in case this shape is of use to you 

# Storage matrix 
k_selec<-array(data=NA, dim = c(8,8))

# Fill it in 
for(i in 1:8){
  for(j in 1:8){
    if(i+j>=9){
      k_selec[i,j]<-max(s_log[,1])
    } else {
      k_selec[i,j]<-0.0001
    }
  }
}



# Matrix with all the selectivities together 
s_mean<-cbind(s_log, domeselec, k_selec)

# Option to scale so selectivity max out at the measure max
old_s_mean<- s_mean # This just allows you to compare to previous vector more easily 
for(i in 1:ncol(s_mean)){
  if (max(s_mean[,i]) > max(s_mean[,1])){
    s_mean[,i]<-s_mean[,i]*(max(s_mean[,1])/max(s_mean[,i])) #scale so that everything maxes at the max from lmesh prediction
  }
}

# Option to standardise the cumulative selectivity - this will 'flatten' a lot but may allow for interesting inferences
# Basically "sacrificing" catchability when selectivity is improved - keep catchability constant 

# old_s_mean<- s_mean #so stuff doesnt get overwritten
# for(i in 1:ncol(s_mean)){
#   for(j in 1:nrow(s_mean)){
#     if (sum(s_mean[,i]) > sum(s_mean[,1])){
#       s_mean[j,i]<-s_mean[j,i]*(sum(old_s_mean[,1])/sum(old_s_mean[,i])) #scale so that everything maxes at the max from lmesh prediction - #!little funky, def think about this 
#     }
#   }
# }



# Make variance matrix 
# Model to provide variance given selectivity 

# Organizing into a data frame 
ldf<-data.frame(var = predS_vari, selec = predS_mean) 

# Fit a model to the JAGS mean selectivity and variance 
l2<-lm(var ~ 0+ selec + I(selec^2), data = ldf) 

# summary(l2) # optional - check your summary 

lc<-coef(l2) # Fit is pretty good. Probably slightly overestimating low variance, which you may need to consider given your specific user case 

# Just putting a little fun(ction) for calculating variance
vari_from_mean<-function(x){ 
  vari<-(lc[1]*x)+(lc[2]*(x^2))  # Calculate variance based off the model coefficients. 
}

# # Alternatively this can be 'hard coded' as 
# vari_from_mean<-function(x){ 
#   vari<-(0.436*x)+(0.408*(x^2))  # Calculate variance based off the hard coded model coefficients. 
# }

vari_max<-function(x){ # Function for the max allowable variance for beta distribution
  vari<-x*(1-x)
}

# s_vari so that it fills in known vari, then pulls from l2, but if the variance is > mean(1-mean) selec then replace 
s_vari<-matrix(data = NA, nrow = nage, ncol = ncol(s_mean))

s_vari <- vari_from_mean(x = s_mean)

for (i in 1:nage){ #this replaces the impossible variances 
  for (j in 1:ncol(s_mean)){
    if (s_vari[i,j]>= vari_max(x=s_mean[i,j])){
      s_vari[i,j]<-vari_max(x = s_mean[i,j]) - (0.1*vari_max(x = s_mean[i,j]))
    }
  }
}

# vari_max.check<-vari_max(x = s_mean) # maximum variance allowed just to look and check

# Calculate ssa and ssb 

# Just setting up the  matrix structure 
ssa<- s_mean 
ssb<- s_mean #ditto above 

for (sscen in 1:ncol(s_mean)){
  ssa[,sscen]<-((1-s_mean[,sscen])/(s_vari[,sscen])-(1/s_mean[,sscen]))*(s_mean[,sscen]^2) #alpha for the selectivity scenario 
  ssb[,sscen]<-(ssa[,sscen]*(1/s_mean[,sscen]-1)) #beta for the selectivity scenario 
}
#! Check to make sure the alphas and betas aren't bad (impossible values) because of the variance calculations if applying to a new dataset 



##............................................
### Initialise effort and other variables####
##............................................

#Set up effort scenarios 
# Effort going from mean to max historic to 50x max historic by half max historic, plus a 100x max historic  
effort.scenario<- c(267,seq(341,3410, by = 170.5), seq(3751, 17050, by = 341), 34100) 

# How many years do you want to project the populations forward?
nysim<-50                                                    

# nyr<- 2                                                     # If using stable state 

# How many total years in storage vectors (simulated plus observed)
ny<-nysim+nyr                                                

# How many effort levels to test
n.scen.harv<-length(effort.scenario)                        

# Number of iterations
niter<-1000

# Age at maturity
## We pull an estimate from Wallsworth, Wallace, et al., (2023), but you can make this whatever age is relevant to you
matage<-matureage

# Read median catchability value from MCMC output
q<-exp(median(mcmcoutput[,colnames(mcmcoutput)=="lnqc"]))    

# Sigma for calculating recruitment anomalies 
sigrec<-1  

# Designate the recruitment scenario
## See below for a description of the different lake scenarios. 
recruit.scenario<-"stockLake" 

# Designate the lake level simulation scenario
## NOT means no lake change 
lake.scenario<-"periodic"  


# Function to simulate fluctuating lake level across simulation years
sim.lake.level<-function(ny){
  deg.to.rad<-pi/360
  lev<-cos(seq(1,720)*deg.to.rad*50)+rnorm(720,0,.25)
  start.yr<-sample(seq(1,720-ny),1)
  stop.yr<-start.yr-1+ny
  return(lev[start.yr:stop.yr])
}

# Read median annual mortality estimate from MCMC output
## You can change this mortality estimate if you have a known estimate for your species
m<-median(mcmcoutput[,colnames(mcmcoutput)=="m"]) 

# N at medians so that the model can initialise 
# If you have different population estimates, initialise them here 
Natqs<-apply(Natmcmc,MARGIN=2,median)
Natmed<-matrix(as.numeric(Natqs),nrow=8,ncol=nyr,byrow=T)
# sum(Natmed[,1]) # Check to make sure things look as expected for your population (NAs are a bad sign!)
# abundquants[3,]

# effort<-array(data=NA,dim=c(nysim+1,n.scen.harv,niter))    # Storage matrix for commercial effort by year
zlakedat<-(lakearea-mean(lakearea))/sd(lakearea)

# Set up array for mean selectivity by effort level 
selec.eff.u<-array(data=NA,dim=c(nage,n.scen.harv,ncol(s_mean)))   # Storage array for Numbers at age

# Fill in the values by selectivity and effort
for (sscen in 1:ncol(s_mean)){
  for (v in seq_along(effort.scenario)){
    for (k in 1:nage){
      selec.eff.u[k,v,sscen]<-sum(rbeta(effort.scenario[v],ssa[k,sscen], ssb[k,sscen]))
    }
  }
}


# Note that running sims forward can be slow - it takes me about an hour
## Note - Code CAN be run on a machine with 16GB RAM, the model output is "chunked" to not exceed this RAM

#### Run simulations forward ####
for (sscen in 1:17){   #! Index will need to change if you added additional selectivity scenarios 
  set.seed(10) # Set seed for reproducibility
  # reset arrays 
  NatPred<-array(data=NA,dim=c(nage,ny,n.scen.harv,niter))   # Storage array for Numbers at age
  catch.by.age<-NatPred                                        # Storage array for catch at age
  harvest.by.age<-catch.by.age                                 # storage array for harvest by age
  total.harvest<-array(data=NA,dim=c(ny,n.scen.harv,niter))  # Storage matrix for total harvest by year
  
  for(ww in 1:niter){
    print(ww) # Just to keep track of your run 
    mcuse<-sample(seq(1,1250),1,FALSE)
    
    #.........................................................................
    ## Create time series of recruitment anomalies
    ##  - USE "stockLake"!
    ##  - "random" recruitment has annual recruitment drawn from one
    ##      of two lognormal distributions. One for large recruitment years,
    ##      and one for low recruitment years
    ##     distributions
    ##  - "stock" recruitment scenario has recruitment driven by a stock 
    ##     recruitment function
    ##  - "stockLake: recruitment scenario has recruitment driven by
    ##     a stock-recruitment function including a lake area effect
    #...........................................................................
    if(recruit.scenario=="random"){
      
      rec.freq<-3         # How many years between big recruitment events (on average)?
      recruit.time<-rep(0,nysim)      # Storage for recruitment time series
      for(i in 1:nysim)           # For each year
      {
        rec.year<-runif(1,0,1)<=(1/rec.freq)       # Is there a big recruitment event?
        if(rec.year) {recruit.time[i]<-exp(rnorm(1,16.5,.5))       # If there is a big recruitment event, draw large recruitment
        } else{recruit.time[i]<-exp(rnorm(1,13.3,.5))}             # Otherwise, draw small recruitment
      }
      
    }
    
    if(recruit.scenario=="stock"){
      arec<-(mcmcoutput[mcuse,colnames(mcmcoutput)=="RecA"])
      brec<-(mcmcoutput[mcuse,colnames(mcmcoutput)=="RecB"])
      rec.anoms<-rnorm(nysim,-sigrec^2/2,sigrec)
      
    }
    
    
    if(recruit.scenario=="stockLake"){
      arec<-(mcmcoutput[mcuse,colnames(mcmcoutput)=="RecA"])  # Ricker alpha value drawn from the estimation model posterior distribution
      brec<-(mcmcoutput[mcuse,colnames(mcmcoutput)=="RecB"])  # Ricker beta value drawn from the estimation model posterior distribution
      lakeb<-(mcmcoutput[mcuse,colnames(mcmcoutput)=="LakeB"]) # Ricker lake level coefficient drawn from the estimation model posterior distribution
      rec.anoms<-rnorm(nysim,-sigrec^2/2,sigrec)  # Time series of random recruitment anomalies
      
      if(lake.scenario!="periodic") zlakeareasim<-rep(0,nysim)   # If lake level is held constant, set lake area to 0 (mean value of observed time series)
      if(lake.scenario=="periodic") zlakeareasim<-sim.lake.level(nysim) # If lake level is periodically changing, simulate lake level time series.
    }  
    
    for(v in 1:length(effort.scenario))                                        # For each harvest scenario
    {
      
      catch.by.age[,1:(nyr-1),v,ww]<-catch.out[,1:(nyr-1)]   # Set initial catch at age from model estimates
      
      #density.by.age[,1:nyr,v,ww]<- rep(0,8)                                 # Set initial density at age from model estimates
      NatPred[,1:nyr,v,ww]<-Natmed                                      # Set initial numbers at age from model estimates
      
      for(i in nyr:ny){                                                        # For each year being simulated...
        
        if(i>nyr){                                                             # If it is not the first year...
          
          
          
          for(j in 2:nage)                                                     # For all ages greater than 0
          {
            if(j<nage) NatPred[j,i,v,ww]<-(NatPred[j-1,i-1,v,ww]- #num individuals at age j is the num at j-1 from previous year minus harvest atage-1 previous year 
                                             catch.by.age[j-1,i-1,v,ww])*(1-m)    # Calculate the number of non-plus-age individuals
            
            if(j==nage) NatPred[j,i,v,ww]<-((NatPred[j-1,i-1,v,ww]-
                                               catch.by.age[j-1,i-1,v,ww])+
                                              (NatPred[j,i-1,v,ww]-
                                                 catch.by.age[j,i-1,v,ww]))*(1-m) # Calculate the number of plus-age individuals
          }
          
          #..........................................
          # Generate Age-0 recruits
          #..........................................
          if(recruit.scenario=="random") NatPred[1,i,v,ww]<-recruit.time[i-nyr]    
          if(recruit.scenario=="stock") NatPred[1,i,v,ww]<-exp(log(arec)-brec*sum(NatPred[c((matage):8),i,v,ww]*wgts[(matage):8])+rec.anoms[i-nyr])*sum(NatPred[c((matage):8),i,v,ww]*wgts[(matage):8])
          if(recruit.scenario=="stockLake") NatPred[1,i,v,ww]<-exp(log(arec)-brec*sum(NatPred[c((matage):8),i,v,ww]*wgts[(matage):8])+zlakeareasim[i-nyr]*lakeb+rec.anoms[i-nyr])*sum(NatPred[c((matage):8),i,v,ww]*wgts[(matage):8])
        }
        
        #......................................................
        # Adjust catchability (q) based on lake area using 
        #  parameter estimates from estimation model posterior
        #  distribution
        #......................................................
        if(i>nyr) q<-exp((mcmcoutput[mcuse,colnames(mcmcoutput)=="lnqc"]+zlakeareasim[i-nyr]*mcmcoutput[mcuse,colnames(mcmcoutput)=="qslope"])) # Calculate q for simulated future years
        if(i<=nyr) q<-exp((mcmcoutput[mcuse,colnames(mcmcoutput)=="lnqc"]+zlakedat[i]*mcmcoutput[mcuse,colnames(mcmcoutput)=="qslope"])) # Calculate q for previously observed years
        
        if(is.na(q)) stop("Q problem")  # Stop the simulations and throw an error message if q calculation above results in NA
        
        #.....................................................
        # Calculate catch at age, biomass harvested at age,
        #  and estimated density at age from simulated 
        #  standardized surveys
        #.....................................................
        for(k in 1:nage)                                           # For each age
        {
          catch.by.age[k,i,v,ww]<-NatPred[k,i,v,ww]*(1-exp(-1*q*selec.eff.u[k,v,sscen]))  # Calculate catch at age
          harvest.by.age[k,i,v,ww]<-catch.by.age[k,i,v,ww]*wgts[k]                               # Calculate biomass harvested at age
          if(is.na(harvest.by.age[k,i,v,ww])) stop(paste("harvest NA i =",i,"v=",v))                   # Throw an error if harvest at age is an NA
        }
        
        total.harvest[i,v,ww]<-sum(harvest.by.age[,i,v,ww])              # Calculate total harvest for the year (harvest summed across all ages)
        if(is.na(total.harvest[i,v,ww])) stop("TotHarv NA")                    # Throw an error if total harvest calculation results in an NA
        
      }
      
    } 
    
  } 
  # Store everything 
  # using saveRDS so I dont overwrite files when I read in :) 
  saveRDS(NatPred, file = paste0("outputs/sscen",sscen,"_NatPred.rds"))
  saveRDS(catch.by.age, file = paste0("outputs/sscen",sscen,"_CBA.rds"))
  saveRDS(harvest.by.age, file = paste0("outputs/sscen",sscen,"_HBA.rds"))
  saveRDS(total.harvest, file = paste0("outputs/sscen",sscen,"_TH.rds"))
  
}

# Afterwards save the whole shebang 
# save.image("outputs/230801_rs_pl_allcum_selecloop.rdata") # Not in GitHub due to size 

# Things to save for plotting - yes redundant with above but helpful for opening a smaller output file if you only want to plot and don't need the whole shebang 
save(effort.scenario,selec.eff.u,ny,niter,wgts,s_mean,s_log, domeselec, k_selec,file = "outputs/plot_info.RData")




# Figures (Manuscript) --------------------------------------------------------- 



# ##   # Plotting the manuscript figures 
# library(rcartocolor)
# library(paletteer) Make color palletes for cohesive figures ---- 

teals<-rev(carto_pal(7,"BluGrn")) # for "current selectivity" effort scenarios 
aro<-carto_pal(12, "Safe") # colorblind friendly pallette 
brgs<-(paletteer_c("grDevices::Burg", 11)) # for log selectivities 
grns<-paletteer_c("grDevices::Emrld", 8) # for knife selectivities
orng<-rev(paletteer_c("grDevices::Oranges", 7))

allcol<-c(brgs,orng[1:6],grns) # a single vector with all selectivityies when you have all the selectivities together 

## Load plot info files ----

load("outputs/plot_info.RData") # plot info from the simulations 

load("outputs/largesmall_jags_plotinfo.RData") # plot info for the JAGS selectivity outputs 

## Figures ----

### Figure 1. ----

# Plot the selectivity scenarios 

pdf("figures/Figure1.pdf", height = 5, width = 7) # PDF output for figure 

par(mar=c(5.6,6.6,3.1,1.25)) # Figure margins 

# For the manuscript, we only plot log selectivities which are indexed [,1-11] in s_mean

plot(c(0:7), s_mean[,1], type = 'l', ylim = c(-.01,1.01), xlim = c(-0.1,7.1), # Plot the baseline selectivity scenario to set up the plot
     xaxt = 'n', yaxt = 'n', yaxs ='i',xaxs='i', xlab = "Age", ylab = "Selectivity", 
     las =1, col=brgs[1],lwd=4, cex=1.2, cex.lab=2.0)
polygon(c(3,3,7,7,3),c(0,1,1,0,0), col = 'seashell', border = NA) # Add a box over all adult age classes for clarity when referenceing reproductive stages of age classes
axis(side=1, at=c(seq(0,7, by = 1)), # Add an x axis
     labels = c("0","1",'2','3','4','5','6','7+'), cex.axis = 2.0)
axis(side=2, at=c(seq(0,1, by = 1)), labels = c("0","1"), cex.axis = 2.0, las=1) # Add a y axis 
for (i in 1:11){ # Add all eleven log selectivity scenarios 
  lines(c(0:7),s_mean[,i], col=brgs[i], lwd=5.5)
}

dev.off() # Turn off PDF plotting window 

### Figure 2. ----

# Plot carp biomass with current selectivity and max historic projected forward

sscen1NatPred<-readRDS("outputs/sscen/sscen1_NatPred.rds") #Read in the simulations from the baseline selectivity scenario 

# Get biomass estimates by multiplying model output by weights 

abunds<-matrix(0, nrow=niter, ncol=ny) # Set up abundance matrix for storing output 


for(j in 1:niter){ # Loop through the NatPred (Predicted number of individuals in each age class at each time step)
  # The number in the sscen1NatPred[,,X,] will be the effort level (effort.scenario[X] gives the corresponding number of hauls)
  abunds[j,]<-as.numeric(colSums((sscen1NatPred[,,2,j]*wgts))) # abundance for each iteration for each year
}
abund_meds<-apply(abunds, MARGIN = 2, FUN = quantile, probs=.5) # Calculate the median abundance across simulation iterations

# Plotting 

pdf("figures/Figure2.pdf", height = 5, width = 5) # PDF for model outputs


par(mar=c(5.6,7.1,3.1,1.25)) # Set plotting margins 


plot(seq_along(abund_meds[14:64]),abund_meds[14:64], # Plot median biomass over time. Starts at 2022 ("0"; abund_meds[14]) and plots 50 years into future (abund_meds[64])
     type = 'l', ylim = c(0,170000000),
     xlim = c(1,51),
     yaxt='n',yaxs='i',xaxt="n", xaxs="i",
     lwd = 3,
     xlab = "Year",
     ylab = "Carp Biomass \nx 1Million Kg",
     col = brgs[1],
     cex.lab = 1.8
)
polygon(c(c(1:51),rev(c(1:51))),c(apply(abunds[,14:64],MARGIN=2,FUN=quantile,probs=0.025), # Add the 95% CI
                                  rev(apply(abunds[,14:64],MARGIN=2,FUN=quantile,probs=0.975))),
        col=adjustcolor(brgs[1], alpha.f = 0.25),
        border=adjustcolor(brgs[1], alpha.f = 0.25),lwd=1
)
polygon(c(c(1:51),rev(c(1:51))),c(apply(abunds[,14:64],MARGIN=2,FUN=quantile,probs=0.25), # Add the 50% CI
                                  rev(apply(abunds[,14:64],MARGIN=2,FUN=quantile,probs=0.75))),
        col=adjustcolor(brgs[1], alpha.f = 0.4),
        border=adjustcolor(brgs[1], alpha.f = 0.4),lwd=1
)
axis(side=2, at=seq(0,170000000,by=15000000),labels=seq(0,170,by=15), cex.axis=1, las =1) # Add the y-axis
axis(side=1,at=seq(1,51, by=5),labels=seq(0,50, by = 5), cex.axis = 1.2) # Add the x-axis


abline(h=abund_meds[1]*.25, lty=2, lwd = 3) # Suppression target - 25% of historic biomass 
abline(h=max(abund_meds), lty=3, lwd =3) # Historic maximum biomass 

legend("topleft", legend = c("Historic Max.","Target"), # Add figure legend 
       lty = c(3,2), lwd = c(3,3), bty = 'n',
       col = c(col = "black","black"),
       cex = 1.3)

dev.off() # Turn off PDF plotting window 



### Figure 3a. ----


# Plot current selectivity and max historic projected forward
sscen1NatPred<-readRDS("outputs/sscen/sscen1_NatPred.rds") # Do not need to reload if it is still loaded from Figure 2. 

# Calculate the biomass by effort level for every year 25-50 years in the future 
biom_by_eff_yr25to50<-array(0, c(niter,length(effort.scenario),length(c(25:50)))) # Storage array 
for(j in 1:niter){
  for(i in 1:26){
    biom_by_eff_yr25to50[j,,i]<-as.numeric(colSums((sscen1NatPred[,i+38,,j]*wgts))) # Biomass for each iteration for each year
  }
}

biom_by_eff<-apply(biom_by_eff_yr25to50, MARGIN = c(1,2), FUN = mean) # Average 
biom_by_eff_meds<-apply(biom_by_eff, MARGIN = 2, FUN = quantile, probs=.5) # Median 

#can calculate medians and then average across years, but I averaged across years and then took medians 
# biom_by_eff_meds<-apply(biom_by_eff, MARGIN = c(2,3), FUN = quantile, probs=.5)
# biom_by_eff_avg<-apply(biom_by_eff_meds, MARGIN = 1, FUN = mean)


rm(sscen1NatPred) # Remove data to save space 

# Plot

pdf("figures/Figure3.pdf", height = 5, width = 5) # PDF Output 

usemar<-c(5.1,7.1,4.1,1.25) # Margins 
par(mar=usemar) # Margins 

plot(effort.scenario/341,biom_by_eff_meds, type = 'l', col=brgs[1], lwd =3, # Plot median carp biomass by effort level 
     ylim =c(0,110000000), 
     cex.lab =1.8, cex.axis=1.2,
     ylab = "Carp Biomass \nx 1Million Kg",
     xlab = 'Effort Multiplier',
     yaxs='i', yaxt = 'n', xaxt = 'n'
)
abline(h=abund_meds[1]*.25, lty=2, lwd = 3) # Suppression target biomass
polygon(c(effort.scenario/341,rev(effort.scenario/341)), # 95% CI
        c(apply(biom_by_eff,MARGIN=2,FUN=quantile,probs=0.025),
          rev(apply(biom_by_eff,MARGIN=2,
                    FUN=quantile,probs=0.975))),
        col=adjustcolor(brgs[1], alpha.f = 0.25),
        border=adjustcolor(brgs[1], alpha.f = 0.25),lwd=1
)
polygon(c(effort.scenario/341,rev(effort.scenario/341)), # 50% CI
        c(apply(biom_by_eff,MARGIN=2,FUN=quantile,probs=0.25),
          rev(apply(biom_by_eff,MARGIN=2,
                    FUN=quantile,probs=0.75))),
        col=adjustcolor(brgs[1], alpha.f = 0.4),
        border=adjustcolor(brgs[1], alpha.f = 0.4),lwd=1
)
axis(side=2,at=seq(0,100000000, by = 10000000), labels = seq(0,100, by=10), las =1, cex.axis=1.3) # Add y-axis
eff.l<-effort.scenario/341 # Calculate effort level in relation to historic maximum effort 
eff.a<-c(1,seq(20,100,by=20)) # For effort axis  
axis(side=1,at=eff.a,labels = eff.a, las = 1, cex.axis=1.2) # Add effort axis 

# dev.off() # Turn off PDF plotting window if you just want 3a



### Figure 3b. ---- 

for (sscen in c(1:1)){ # Only plotting the 1 sscen for now, but the supplement will loop through same code  
  TH<-readRDS(paste0("outputs/sscen/sscen",sscen,"_TH.rds")) # Read in each selectivity scenario file 
  TH_y25to50_avg<-apply(TH[39:64,,], MARGIN = c(2,3), FUN = mean) # Mean 
  TH_meds_a<-apply(TH_y25to50_avg, MARGIN = 1, FUN = quantile, probs=.5) # Median 
  TH_.025_a<-apply(TH_y25to50_avg, MARGIN = 1, FUN = quantile, probs=.025) # CI
  TH_.975_a<-apply(TH_y25to50_avg, MARGIN = 1, FUN = quantile, probs=.975) # CI
  TH_.25_a<-apply(TH_y25to50_avg, MARGIN = 1, FUN = quantile, probs=.25) # CI
  TH_.75_a<-apply(TH_y25to50_avg, MARGIN = 1, FUN = quantile, probs=.75) # CI
  
  
  plot(effort.scenario/341,TH_meds_a, type = 'l', col=brgs[1], lwd =3, # Median 
       ylim =c(0,10000000), cex.lab =1.8, cex.axis=1.2,
       ylab = 'Harvest Biomass \nx 1Million Kg', xlab = 'Effort Multiplier',
       yaxs='i', yaxt = 'n', xaxt = 'n')
  polygon(c(effort.scenario/341,rev(effort.scenario/341)), # CI 
          c(TH_.025_a,rev(TH_.975_a)),
          col=adjustcolor(brgs[1], alpha.f = 0.25),
          border=adjustcolor(brgs[1], alpha.f = 0.25),lwd=1
  )
  polygon(c(effort.scenario/341,rev(effort.scenario/341)), # CI
          c(TH_.25_a,rev(TH_.75_a)),
          col=adjustcolor(brgs[1], alpha.f = 0.4),
          border=adjustcolor(brgs[1], alpha.f = 0.4),lwd=1
  )
  axis(side=2,at=seq(0,10000000, by = 1000000), labels = seq(0,10, by=1), las =1, cex.axis=1.2) # Y-axis
  axis(side=1,at=eff.a,labels = eff.a, las = 1, cex.axis=1.2) # X-axis
  
}


dev.off() # Turn off PDF plotting window 


### Figure 4. ----


# Get the effort level that achieves target for all selectivities

bbe_25to50<-biom_by_eff_yr25to50 # Temporary array that wont overwrite above
target_effort<-rep(NA,11) # Empty vector to store effort needed for each selectivity

# This will fill in the j, a , t selectivities, and in selec.ja[4] it will be a t/f if <50%
for(sscen in 1:ncol(s_mean[,1:11])){
  NatPred<-readRDS(paste0("outputs/sscen/sscen",sscen,"_NatPred.rds"))
  for(j in 1:niter){
    for(i in 1:26){
      bbe_25to50[j,,i]<-as.numeric(colSums((NatPred[,i+38,,j]*wgts)))
    }
  }
  bbe<-apply(bbe_25to50, MARGIN = c(1,2), FUN = mean)
  bbe_meds<-apply(bbe, MARGIN = 2, FUN = quantile, probs=.5)
  meets_target<-as.numeric(bbe_meds<0.25*max(abund_meds))
  target_effort[sscen]<-eff.l[which.max(meets_target)]
}




# Plot of the median target eff level and median msy eff level across sscen 

MSY_at<-c(1:11) # Empty vector to store MSY data - the effort level where msdy is reached 

for (sscen in c(1:11)){ #  For each log selectivity scenario 
  
  TH<-readRDS(paste0("outputs/sscen/sscen",sscen,"_TH.rds")) # Total harvest 
  TH_avg<-apply(TH[39:64,,], MARGIN = c(2,3), FUN = mean) # Average total harvest
  TH_med<-apply(TH_avg, MARGIN = 1, FUN = quantile, probs=.5) # Median total harvest
  TH_max<-which.max(TH_med) # Maximum total harvest 
  MSY_at[sscen]<-eff.l[TH_max] # MSY 
}

# Cumulative selectivity 
cum_selec<-apply(s_mean[,1:11],MARGIN = 2, FUN=sum)/8 # Cumulative selectivities 

pdf("figures/Figure4.pdf", height = 5, width = 5) # PDF for outputs 

par(mar=c(5.6,6.6,3.1,1.25)) # plot margins 

plot(cum_selec, target_effort, pch = 16, # Target 
     col = brgs,
     xlab = "Average Selectivity", ylab = 'Effort Multiplier',
     ylim =c(0,19), las=1, cex= 2.0, cex.axis =1.2,cex.lab=1.8,
     xlim = c(0.2,0.8))
points(cum_selec, MSY_at, pch = 2, col= brgs, cex=2.2, lwd =3.5) # MSY
legend("topright",pch=c(16,2), col=c(brgs[1],brgs[1]), legend = c("Target","MSY"), # Add a legend 
       cex=1.5, bty='n')

dev.off() # Turn off PDF plotting window 




# Additional Calculations ----

# What percent of simulations are below target year 
abund10.50<-abunds[,25:64]
a1<-apply(abund10.50, 1, mean)
a11<- which(a1 <= abund_meds[1]*.25)
length(a11)/1000


a1<-apply(abund10.50, 1, mean)

a15<-which(abund10.50 <= abund_meds[1]*.25, arr.ind = TRUE)
a15<-which(abund10.50 <= abund_meds[1]*.25)

length(a15)/41
195.297/1000 # 19 percent of iterations 




# Figures (Supplement) ---- 

# # Load packages
# library(rcartocolor)
# library(paletteer) 

## Load plot info files ----
  # These are the same files as above, so they are commented out. Uncomment them if you have not run "Figures"

# load("outputs/plot_info.RData") # plot info from the simulations 
# 
# load("outputs/largesmall_jags_plotinfo.RData") # plot info for the JAGS selectivity outputs 

### Figure S2. ----

pdf("figures/supplement/s2.pdf", height = 5, width = 7) # Save file as a PDF output in the JAE Figures folder

par(mar=c(5.6,6.6,3.1,1.25)) # Figure margins 

plot(seq(0,7,1),seq(0,1,length.out=8), col = scales::alpha ('white',1), # Plot an empty plot
     ylim = c(-.01,1.01), xlim = c(-0.1,7.1), 
     xaxt = 'n', yaxt = 'n', yaxs ='i',xaxs='i', xlab = "Age", ylab = "Selectivity", 
     las =1, lwd=4, cex=1.2, cex.lab=2.0)
axis(side=1, at=c(seq(0,7, by = 1)), # Add axis labels (year classes 0 to 7+) to the x axis 
     labels = c("0","1",'2','3','4','5','6','7+'), cex.axis = 2.0)
axis(side=2, at=c(seq(0,1, by = 1)), labels = c("0","1"), cex.axis = 2.0, las=1) # Add axis labels (0 to 1) to the y axis 
for (i in 1:length(l_lci_posth)){ # Add the large mesh age based CI
  lines(c((i-1.05),(i-1.05)),c(l_lci_posth[i],l_uci_posth[i]), col = brgs[1], lwd =1.5)
}
for (i in 1:length(s_lci_posth)){ # Add the small mesh age based CI
  lines(c((i-.95),(i-.95)),c(s_lci_posth[i],s_uci_posth[i]), col = aro[11], lwd =1.5)
}
polygon(c(age200, rev(age200)), c(s_lci_pred, rev(s_uci_pred)), # Add a white shape over the CI to make the model CI less awkward
        col = 'white', lty = 0) 
polygon(c(age200, rev(age200)), c(l_lci_pred, rev(l_uci_pred)), # Add a white shape over the CI to make the model CI less awkward
        col = 'white', lty = 0) 
polygon(c(age200, rev(age200)), c(l_lci_pred, rev(l_uci_pred)), # Large mesh model CI
        col = scales::alpha(brgs[1],.5), lty = 0)
polygon(c(age200, rev(age200)), c(s_lci_pred, rev(s_uci_pred)), # Small mesh model CI
        col = scales::alpha(aro[1],.4), lty = 0)
lines(age200,s_med_pred, type = "l", col = aro[11],lwd =3) # Small mesh model median
lines(age200,l_med_pred, type = "l", col = brgs[1],lwd =3) # Large mesh model median
legend("topleft",inset=c(.018,.02), # add a legend 
       legend= c("Large Mesh (Baseline)", "Small Mesh"),
       col = c(brgs[1], aro[11]),
       cex=1.2,lwd =c(4,4), bty = "n",
       lty = c(1,1))

dev.off() # Turn off the PDF plotting window 




### Figure S3. ----
## Note that this layout formatting is slightly different than the published supplemental figure (But the content is identical). The formatting takes up a lot of lines of code and has been ommited for clarity. 
# Adding in additional scenarios for supplemental material to address R1 feedback 

pdf("figures/supplement/s3.pdf", height = 11, width = 8.5) # Set up pdf to save all plots on one page 

layout(mat = matrix(c(1,6:14,2,15:23,3,24:32,4,33:41,5,42:50), ncol=5, byrow=F),
       heights = c(2,1,1,1,1,1,1,1,1,1),
       widths = c(1,1,1,1,1))
par(oma=c(10, 10, 3, 3), mar=c(6,1.5,2,1))

#plot the selectivity curves for the top five plots 
#which selectivity scenarios I am using


selec_examples<- c(1,3,5,7,11)
for (i in 1:length(selec_examples)){
  plot(c(0:7), s_mean[,1], type = 'l', ylim = c(-.01,1.01), xlim = c(-0.1,7.1),
       xaxt = 'n', yaxt = 'n', yaxs ='i',xaxs='i',
       las =1, col="white",lwd=4, cex=1.2, cex.lab=1.0) #using formatting from before but this is just a blank window 
  polygon(c(3,3,7,7,3),c(0,1,1,0,0), col = 'seashell', border = NA)
  axis(side=1, at=c(seq(0,7, by = 1)), 
       labels = c("0","1",'2','3','4','5','6','7+'), cex.axis = 1.0)
  if(i == 1){
    axis(side=2, at=c(seq(0,1, by = 1)), labels = c("0","1"), cex.axis = 1.0, las=1)
    ylab = "Selectivity"
  } else {
    axis(side=2, at=c(seq(0,1, by = 1)), labels = NA, cex.axis = 1.0, las=1)
  }
  lines(c(0:7),s_mean[,selec_examples[i]], col=brgs[selec_examples[i]], lwd=5.5)
}


#vector of effort levels to loop over for the supporting figure (1.5, 2.5,5,11,17,24,50,100)
effort_examples<-c(2,3,5,10,21,27,34,60,61)

# Loop through 
for(sscen in selec_examples){
  NatPred<-readRDS(paste0("outputs/sscen/sscen",sscen,"_NatPred.rds")) 

  
  abunds_effort_examples<-array(0, dim = c(niter, ny, length(effort_examples))) #storage for our abundance outputs 
  for(i in 1:length(effort_examples)){
    
    for(j in 1:niter){
      #the number in the sscen1NatPred will be the effort level
      #abunds[j,]<-as.numeric(colSums((sscen1NatPred[,,i,j]*wgts))) #abundance for each iteration for each year with 1.5x historic effort
      abunds_effort_examples[j,,i]<-as.numeric(colSums((NatPred[,,effort_examples[i],j]*wgts)))
    }
    
  }
  
  #medians
  abund_meds_effort_examples<-apply(abunds_effort_examples, MARGIN = c(2,3), FUN = quantile, probs=.5) #save that space 
  
  #upper 95th percentile
  abund_u95_effort_examples<-apply(abunds_effort_examples, MARGIN = c(2,3), FUN = quantile, probs=0.975)
  #lower 95th percentile 
  abund_l95_effort_examples<-apply(abunds_effort_examples, MARGIN = c(2,3), FUN = quantile, probs=.025)
  #upper 50th percentile 
  abund_u50_effort_examples<-apply(abunds_effort_examples, MARGIN = c(2,3), FUN = quantile, probs=0.75)
  #lower 50th percentile 
  abund_l50_effort_examples<-apply(abunds_effort_examples, MARGIN = c(2,3), FUN = quantile, probs=0.25)
  
  par( mar=c(1,1,1,1))
  for(i in 1:length(effort_examples)){
    
    plot(seq_along(abund_meds_effort_examples[14:64,i]),abund_meds_effort_examples[14:64,i],#this is 2022 plus 50 years into future 
         type = 'l', ylim = c(0,170000000),
         xlim = c(1,51),
         yaxt='n',yaxs='i',xaxt="n", xaxs="i",
         lwd = 3,
         xlab = "",
         ylab = "",
         col = brgs[sscen],
         #cex.lab = 1
    )
    polygon(c(c(1:51),rev(c(1:51))),c((abund_l95_effort_examples[14:64,i]),
                                      rev(abund_u95_effort_examples[14:64,i])),
            col=adjustcolor(brgs[sscen], alpha.f = 0.25),
            border=adjustcolor(brgs[sscen], alpha.f = 0.25),lwd=1
    )
    polygon(c(c(1:51),rev(c(1:51))),c((abund_l50_effort_examples[14:64,i]),
                                      rev(abund_u50_effort_examples[14:64,i])),
            col=adjustcolor(brgs[sscen], alpha.f = 0.4),
            border=adjustcolor(brgs[sscen], alpha.f = 0.4),lwd=1
    )
    axis(side=2, at=seq(0,170000000,by=15000000), labels=seq(0,170,by=15), cex.axis=1.2, las =1) #
    if(i==1){ylab = "Max Historic Effort"}
    if(i==2){ylab = "1.5X Effort"}
    if(i==3){ylab = "2.5x Effort"}
    if(i==4){ylab = "5x Effort"}
    if(i==5){ylab = "11x Effort"}
    if(i==6){ylab = "17x Effort"}
    if(i==7){ylab = "24x Effort"}
    if(i==8){ylab = "50x Effort"}
    if(i==9){ylab = "100x Effort"}
    
    if(i == 9){
      axis(side=1,at=seq(1,51, by=5),labels=seq(0,50, by = 5), cex.axis = 1)} else{
        axis(side=1,at=seq(1,51, by=5), labels = NA)  
      }
    #these relate to the abund of the real carp pop, not the simulated carp pop!
    abline(h=abund_meds[1]*.25, lty=2, lwd = 1)
    abline(h=max(abund_meds), lty=3, lwd =1)
    box()
    # legend("topleft", legend = c("Historic Max.","Target"),
    #        lty = c(3,2), lwd = c(3,3), bty = 'n',
    #        col = c(col = "black","black"),
    #        cex = 1.3) 
  }
  
}

dev.off()

### Figure S4 ----

pdf("figures/supplement/s4.pdf", height = 11, width = 8.5) # Set up pdf to save all plots on one page 
layout(mat = matrix(c(1:10,21:30,11:20), ncol=3, byrow=F), widths = c(3,1,3))
par(oma=c(20, 10, 5, 10), mar=c(0,1.5,2,1.5))

# Biomass over effort levels, 1 plot per sscen 

biom_by_eff_yr25to50<-array(0, c(niter,length(effort.scenario),length(c(25:50)))) # Storage array 


for (sscen in c(2:11)){ #this was looped for ESA, but i just want current ish for pub 
  NatPred<-readRDS(paste0("outputs/sscen/sscen",sscen,"_NatPred.rds")) 
  # Calculate the biomass by effort level for every year 25-50 years in the future 
  for(j in 1:niter){
    for(i in 1:26){
      biom_by_eff_yr25to50[j,,i]<-as.numeric(colSums((NatPred[,i+38,,j]*wgts))) # Biomass for each iteration for each year
    }
  }
  biom_by_eff<-apply(biom_by_eff_yr25to50, MARGIN = c(1,2), FUN = mean) # Average 
  biom_by_eff_meds<-apply(biom_by_eff, MARGIN = 2, FUN = quantile, probs=.5) # Median 
  
  # Plot
  plot(effort.scenario/341,biom_by_eff_meds, type = 'l', col=brgs[sscen], lwd =3, # Plot median carp biomass by effort level 
       ylim =c(0,110000000), 
       cex.lab =1.8, cex.axis=1.2,
       
       
       yaxs='i', yaxt = 'n', xaxt = 'n', log = "x"
  )
  abline(h=abund_meds[1]*.25, lty=2, lwd = 1.5) # Suppression target biomass
  polygon(c(effort.scenario/341,rev(effort.scenario/341)), # 95% CI
          c(apply(biom_by_eff,MARGIN=2,FUN=quantile,probs=0.025),
            rev(apply(biom_by_eff,MARGIN=2,
                      FUN=quantile,probs=0.975))),
          col=adjustcolor(brgs[sscen], alpha.f = 0.25),
          border=adjustcolor(brgs[sscen], alpha.f = 0.25),lwd=1
  )
  polygon(c(effort.scenario/341,rev(effort.scenario/341)), # 50% CI
          c(apply(biom_by_eff,MARGIN=2,FUN=quantile,probs=0.25),
            rev(apply(biom_by_eff,MARGIN=2,
                      FUN=quantile,probs=0.75))),
          col=adjustcolor(brgs[sscen], alpha.f = 0.4),
          border=adjustcolor(brgs[sscen], alpha.f = 0.4),lwd=1
  )
  abline(v=target_effort[sscen])
  #if(sscen == 2| sscen == 4| sscen == 6| sscen == 8| sscen == 10){
  axis(side=2,at=seq(0,100000000, by = 50000000), labels = seq(0,100, by=50), las =1, cex.axis=1) # Add y-axis 
  ylab = "Carp Biomass \nx 1Million Kg"
  #}
  if(sscen == 11){
    eff.l<-effort.scenario/341 # Calculate effort level in relation to historic maximum effort 
    eff.a<-c(1,seq(20,100,by=20)) # For effort axis  
    axis(side=1,at=eff.a,labels = eff.a, las = 1, cex.axis=1) # Add effort axis 
    xlab = 'Effort Multiplier'
  }
  
  
}


for (sscen in c(2:11)){ #looped  through all but the one in the pub 
  TH<-readRDS(paste0("outputs/sscen/sscen",sscen,"_TH.rds"))
  TH_y25to50_avg<-apply(TH[39:64,,], MARGIN = c(2,3), FUN = mean)
  TH_meds_a<-apply(TH_y25to50_avg, MARGIN = 1, FUN = quantile, probs=.5)
  TH_.025_a<-apply(TH_y25to50_avg, MARGIN = 1, FUN = quantile, probs=.025)
  TH_.975_a<-apply(TH_y25to50_avg, MARGIN = 1, FUN = quantile, probs=.975)
  TH_.25_a<-apply(TH_y25to50_avg, MARGIN = 1, FUN = quantile, probs=.25)
  TH_.75_a<-apply(TH_y25to50_avg, MARGIN = 1, FUN = quantile, probs=.75)
  
  
  plot(effort.scenario/341,TH_meds_a, type = 'l', col=brgs[sscen], lwd =3, 
       ylim =c(0,10000000), cex.lab =1, cex.axis=1,
       yaxs='i', yaxt = 'n', xaxt = 'n', log = "x")
  polygon(c(effort.scenario/341,rev(effort.scenario/341)),
          c(TH_.025_a,rev(TH_.975_a)),
          col=adjustcolor(brgs[sscen], alpha.f = 0.25),
          border=adjustcolor(brgs[sscen], alpha.f = 0.25),lwd=1
  )
  polygon(c(effort.scenario/341,rev(effort.scenario/341)),
          c(TH_.25_a,rev(TH_.75_a)),
          col=adjustcolor(brgs[sscen], alpha.f = 0.4),
          border=adjustcolor(brgs[sscen], alpha.f = 0.4),lwd=1
  )
  abline(v=MSY_at[sscen])
  #axis(side=1,at=eff.a,labels = eff.a, las = 1, cex.axis=1.2)
  # if(sscen == 2| sscen == 4| sscen == 6| sscen == 8| sscen == 10){
  axis(side=4,at=seq(0,10000000, by = 5000000), labels = seq(0,10, by=5), las =1, cex.axis=1)
  #}
  if(sscen == 11){
    eff.l<-effort.scenario/341 # Calculate effort level in relation to historic maximum effort 
    eff.a<-c(1,seq(20,100,by=20)) # For effort axis  
    axis(side=1,at=eff.a,labels = eff.a, las = 1, cex.axis=1) # Add effort axis 
  }
  
}

# Turn off plotting window 
dev.off()



# Thanks for reading through and/or using this code! 
# Please feel free to Rae Fadlovich via GitHub (rfadlovi) or email with any questions. 

