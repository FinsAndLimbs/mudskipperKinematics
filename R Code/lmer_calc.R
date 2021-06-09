###############################################################################
##                    MudskipperKine_Analysis.R lmer calc                    ##
##                                                                           ##
## Description:                                                              ##
##     1) Takes two dataframes as arguments (eg. pb_KneeAng_Combined,        ##
##        at_pec_KneeAng_Combined, max) and outputs a data frame containing  ##
##        lmer results calulated for either max or min joint angles          ##
##     2) Merges dataframes (keeping columns consistent) and adds a column   ##
##        containing individual ID's. Saves max or min joint angle for each  ##
##        trial (1 row = 1 trial) to addditional columns. Then calculates    ##
##        lmer with lmer(Max_Min ~ species + (1|Individual))                 ##
##                                                                           ##
###############################################################################                                                                           ##

## TODO: (@Zach) 1)Troubleshoot Tmax & Tmin singularity issues. Both are 
##                 commented out until the issue is fixed
##               2)Update function description to include fixedEffect argument

library(lme4)

lmer_calc <- function(appendage1, appendage2, fixedEffect = "species") {
    #merge dataframes by column
     jointAngMerged <- merge(appendage1, appendage2, all = T) 
  
    #save individual ID's (eg. at01, pb05) in new column 'Ind'
     jointAngMerged$Ind <- substring(jointAngMerged$filename, 1, 4)
     
    #save max/min of each trial (row) and save index (%stance)
    #of each max/min in new columns
      jointAngMerged$Max  <- apply(jointAngMerged[, 1:101], 1, max)
     jointAngMerged$Min  <- apply(jointAngMerged[, 1:101], 1, min)
    #jointAngMerged$Tmax <- apply(jointAngMerged[, 1:101], 1, which.max)
    #jointAngMerged$Tmin <- apply(jointAngMerged[, 1:101], 1, which.min)
     
     #run lmers for max and min joint angles, and the %stance of max and min
     if(fixedEffect == "species") {
       lmer_max <- lmer(Max ~ species + (1|Ind), data = jointAngMerged)
       lmer_min <- lmer(Min ~ species + (1|Ind), data = jointAngMerged)
      #lmer_Tmax <- lmer(Tmax ~ species + (1|Ind), data = jointAngMerged)
      #lmer_Tmin <- lmer(Tmin ~ species + (1|Ind), data = jointAngMerged)
     }
     
     else if(fixedEffect == "appendage"){
       lmer_max <- lmer(Max ~ appendage + (1|Ind), data = jointAngMerged)
       lmer_min <- lmer(Min ~ appendage + (1|Ind), data = jointAngMerged)
       #lmer_Tmax <- lmer(Tmax ~ appendage + (1|Ind), data = jointAngMerged)
       #lmer_Tmin <- lmer(Tmin ~ appendage + (1|Ind), data = jointAngMerged)
     }
     
     else {
       stop("Please specify the fixed effect: 'species' or 'appendage'")
     }
     
    #combine all lmer results into one list and return list
     lmer_results <- list(lmer_max, lmer_min)
    #lmer_results <- list(lmer_max, lmer_min, lmer_Tmax, lmer_Tmin)
     
     return(lmer_results)
}