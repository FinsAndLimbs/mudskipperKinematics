################################################################################
##                    Analyzing mudskipper kinematics on land                 ##
##                                    Sept. 9, 2020                           ##
##                                                                            ##
################################################################################

##TODO:
##1) @Zach Fix singularity issues with lmers (Especially Tmax/Tmin)
##2) @Zach Just leave in the joint angle stats section. Clean it up and make it 
##         pretty once the ms is written


## Note: terminology may still reflect salamander anatomy for convenience

## Clear Environment (workspace)
rm(list = ls())

## remove scientific notation from console
options(scipen = 999)

## collecting info about today's date 
today <- Sys.Date()
SaveDate <- format(today, format="%y%m%d")


#### STEP 1: LOAD LIBRARIES ####

## install libraries
#install.packages("devtools")
#install.packages("performance")

## load libraries
library(devtools)    # for install_github()
library(ggplot2)     # for ggplot()
library(sciplot)     # for se()
library(pspline)     # for smooth.Pspline()
library(signal)      # for interp1()
library(cowplot)     # for plot_grid()
library(grid)        # for textGrob()
library(gridExtra)   # for arrangeGrob() and grid.arrange()
library(lme4)        # for lmer()
library(performance) # for r2_xu()
library(emmeans)     # for emmeans() and pairs()

## use devtools to load the kraken repo from GitHub
?install_github # loads the help file for this function
install_github("MorphoFun/kraken", dependencies = TRUE)


## Now load the kraken repo as a library, so we can use the functions in my repo
library(kraken)

## list all functions and datasets within a package
ls("package:kraken")

percentStance <- seq(0,100,1)
Stance5 <- seq(1,101,5)


#### STEP 2: LOAD XY DATA ####
# Choose the directory containing the digitizing files (e.g., Step1_DigitizeVideos)
xy_path <- choose.dir(caption = "Choose the Step1_DigitizeVideos folder")
xy_list <- list.files(xy_path, pattern="xypts.csv", full=TRUE)


# Naming matrices in array with trial name
# e.g., pb01f01d
Trial <- list(substring(basename(xy_list), 1, 8))

# Identifying group
# e.g., Pec
Group <- list(substring(basename(xy_list), 10, 12))

# Getting data

Pb_RawFiles <- array(lapply(xy_list, read.csv, header=T), dimnames=Trial)
for (i in 1:length(Pb_RawFiles)) {
  colnames(Pb_RawFiles[[i]]) <- c("ToeTip_X", "ToeTip_Y", "KneeElbow_X", "KneeElbow_Y", "HipShoulder_X", "HipShoulder_Y",
                               "AnkleWrist_X", "AnkleWrist_Y", "MetaPhalJoint_X", "MetaPhalJoint_Y",
                               "Costal1_X", "Costal1_Y", "Costal2_X", "Costal2_Y")
}

## Separate dorsal from lateral views
Pb_RawFiles_d <- Pb_RawFiles[which(substring(names(Pb_RawFiles), 8,8)=="d")]
Pb_RawFiles_l <- Pb_RawFiles[which(substring(names(Pb_RawFiles), 8,8)=="l")]


## The order of the points in the xypts.csv mudskipper output from DTLdataviewer is: 
# 1. tip of pectoral fin
# 2. middle of "elbow" (intra-fin joint)
# 3. shoulder
# 4. middle of "wrist" 
# 5. "metacarpophalangeal" joint
# 6. midline point 1 (anterior); same as Costal1
# 7. midline point 2 (posterior); same as Costal2

## and the order of points in the MATLAB code is: 
# 1. tip of toe
# 2. Elbow
# 3. Shoulder
# 4. Wrist
# 5. Metacarpophalangeal joint
# 6. Midline 1 / Costal1
# 7. Midline 2 / Costal2

## so don't need to rearrange the order of the points!



#### STEP 3: SMOOTHING ####
## We'll use the kraken::smootheR() function to use a spline-based method with custom smoothing parameters
## The custom smoothing parameters are obtained by digitizing a file 3x and then taking the variance between these attempts
## That tells you how variable the digitizing was for each point, so the smoothing function can correct accordingly
?smootheR

## Calculating our custom smoothing parameters
# We'll calculate custom smoothing parameters from our repeatability efforts

# Retrieving repeatability files
RepeatPath <- paste(dirname(xy_path), "/Repeatability", sep = "")
RepeatFiles_list <- list.files(RepeatPath, pattern="xypts.csv", full=TRUE)

# Naming matrices in array with trial name
RepeatTrial <- list(substring(basename(RepeatFiles_list), 1, 12))

# Getting data
RepeatFiles <- array(lapply(RepeatFiles_list, read.csv, header=T, nrows=10), dimnames=RepeatTrial)
for (i in 1:length(RepeatFiles)) {
  colnames(RepeatFiles[[i]]) <- c("ToeTip_X", "ToeTip_Y", "KneeElbow_X", "KneeElbow_Y", "HipShoulder_X", "HipShoulder_Y", 
                                  "AnkleWrist_X", "AnkleWrist_Y", "MetaPhalJoint_X", "MetaPhalJoint_Y", 
                                  "Costal1_X", "Costal1_Y", "Costal2_X", "Costal2_Y")
}

VarNames <- list(c("ToeTip_X", "ToeTip_Y", "KneeElbow_X", "KneeElbow_Y", "HipShoulder_X", "HipShoulder_Y", 
                   "AnkleWrist_X", "AnkleWrist_Y", "MetaPhalJoint_X", "MetaPhalJoint_Y", 
                   "Costal1_X", "Costal1_Y", "Costal2_X", "Costal2_Y"))

RepeatFiles_PbPecDEx <- RepeatFiles[which(substring(names(RepeatFiles), 8,8)=="d")]
RepeatFiles_PbPecLEx <- RepeatFiles[which(substring(names(RepeatFiles), 8,8)=="l")]

NumFrames <- nrow(as.data.frame((RepeatFiles[1]))) # digitizers analyzed the same 10 frames
NumCoordCols <- ncol(RepeatFiles[[1]])
NumRepeats <- length(RepeatFiles_PbPecDEx)
Total <- NumFrames*NumCoordCols*NumRepeats

# Pb Pec Example dorsal view
for (j in 1:NumRepeats) {
  RepeatFiles_PbPecDEx[[j]] <- t(RepeatFiles_PbPecDEx[[j]])
}
PbPecDEx <- array(NA, dim=c(NumRepeats, NumFrames, NumCoordCols))

# calculating the variance bewteen the repeat digitizing attempts
PbPecDEx_Var <- matrix(NA, ncol=NumFrames, nrow=NumCoordCols)

for (i in 1:NumCoordCols) {
  PbPecDEx[,,i] <- rbind(RepeatFiles_PbPecDEx[[1]][i,], RepeatFiles_PbPecDEx[[2]][i,], RepeatFiles_PbPecDEx[[3]][i,])
  PbPecDEx_Var[i,] <- (apply(PbPecDEx[,,i], 2, var))
}

# calculating a mean variance for each variable based on the variances calculated for each frame
PbPecDEx_VarMean <- rowMeans(PbPecDEx_Var)
PbPecDEx_VarSE <- (apply(PbPecDEx_Var, 1, sd))/sqrt(3)
PbPecDEx_VarSave <- data.frame(rbind(PbPecDEx_VarMean, PbPecDEx_VarSE))
names(PbPecDEx_VarSave) <- c("ToeTip_X", "ToeTip_Y", "KneeElbow_X", "KneeElbow_Y", "HipShoulder_X", "HipShoulder_Y", 
                             "AnkleWrist_X", "AnkleWrist_Y", "MetaPhalJoint_X", "MetaPhalJoint_Y", 
                             "Costal1_X", "Costal1_Y", "Costal2_X", "Costal2_Y")
PbPecDEx_VarSave$Average <- apply(PbPecDEx_VarSave, 1, mean)



# Pb Pec Example lateral view
for (j in 1:NumRepeats) {
  RepeatFiles_PbPecLEx[[j]] <- t(RepeatFiles_PbPecLEx[[j]])
}
PbPecLEx <- array(NA, dim=c(NumRepeats, NumFrames, NumCoordCols))
PbPecLEx_Var <- matrix(NA, ncol=NumFrames, nrow=NumCoordCols)

for (i in 1:NumCoordCols) {
  PbPecLEx[,,i] <- rbind(RepeatFiles_PbPecLEx[[1]][i,], RepeatFiles_PbPecLEx[[2]][i,], RepeatFiles_PbPecLEx[[3]][i,])
  PbPecLEx_Var[i,] <- (apply(PbPecLEx[,,i], 2, var))
}

PbPecLEx_VarMean <- rowMeans(PbPecLEx_Var)
PbPecLEx_VarSE <- (apply(PbPecLEx_Var, 1, sd))/sqrt(3)
PbPecLEx_VarSave <- data.frame(rbind(PbPecLEx_VarMean, PbPecLEx_VarSE))
names(PbPecLEx_VarSave) <- c("ToeTip_X", "ToeTip_Y", "KneeElbow_X", "KneeElbow_Y", "HipShoulder_X", "HipShoulder_Y", 
                             "AnkleWrist_X", "AnkleWrist_Y", "MetaPhalJoint_X", "MetaPhalJoint_Y", 
                             "Costal1_X", "Costal1_Y", "Costal2_X", "Costal2_Y")

PbPecLEx_VarSave$Average <- apply(PbPecLEx_VarSave, 1, mean)

# We can now use the mean values in PbPecDEx_VarSave and PbPecLEx_VarSave
# as custom smoothing parameters

# we'll use method = 1 to use the value supplied in spar as our custom smoothing parameter
Pb_smoothed_d <- lapply(Pb_RawFiles_d, FUN = function(x) smootheR(x, method = 1, spar = PbPecDEx_VarSave[1,-length(PbPecDEx_VarSave)]))
Pb_smoothed_l <- lapply(Pb_RawFiles_l, FUN = function(x) smootheR(x, method = 1, spar = PbPecLEx_VarSave[1,-length(PbPecLEx_VarSave)]))




#### STEP 4: INTERPOLATION TO 101 POINTS ####
## Since the number of frames per trials differ, we want to standardize this value across the trials
## We can interpolate the frames to 101 points, so they represent 1%, from 0 to 100% of stance
## 
## we'll use the kraken::interpolater() function for this
?interpolateR

Pb_smooth_d_interp <- lapply(Pb_smoothed_d, FUN = function(x) interpolateR(x, 101))
Pb_smooth_l_interp <- lapply(Pb_smoothed_l, FUN = function(x) interpolateR(x, 101))

### Save these data so they can be analyzed with Rick's MATLAB code
# Saving each file as a separate .txt file
SmoothInterpPath <- paste(dirname(xy_path), "/Step2_SmoothInterpData", sep = "")
lapply(1:length(Pb_smooth_d_interp), FUN = function(x) write.table(data.frame(Pb_smooth_d_interp[[x]]), file = paste(SmoothInterpPath, "/", names(Pb_smooth_d_interp[x]), "_Smooth_101.txt", sep = ""), sep = "\t", row.names = FALSE))
lapply(1:length(Pb_smooth_l_interp), FUN = function(x) write.table(data.frame(Pb_smooth_l_interp[[x]]), file = paste(SmoothInterpPath, "/", names(Pb_smooth_l_interp[x]), "_Smooth_101.txt", sep = ""), sep = "\t", row.names = FALSE))


#### STEP 5a: CALCULATE KINEMATIC DATA - Pb ####

## Now read in the Periophthalmus barbarus kinematic data
## Select the *folder* containing the .KIN files (e.g., Step3_ProcessInMATLAB-->mudskipper fore-->forceanal results)
## Kinematics were calculated from MATLAB code to expedite the results, so loading files here
pb_kin_path <- choose.dir(caption = "select the mudskipper fore - forceanal results folder")
pb_kin_list <- list.files(pb_kin_path, pattern="KIN.txt", full=TRUE)

# Naming matrices in array with trial name (e.g., pb05f12)
# removing the directory info from the directory path
# this duplicates some of the steps taken before, but it helps to make sure that we correctly match the kinematic files with the correct trials
pb_kin_filenames <- basename(pb_kin_list)
pb_kin_trial <- list(substring(pb_kin_filenames, 1, 7))

# Identifying group
# e.g., Pec
pb_kin_group <- list(substring(pb_kin_filenames, 9, 11))


# Getting data
kinFiles_pb <- array(lapply(pb_kin_list, read.delim, header=T), dimnames=pb_kin_trial)

## pb - Yaw angle 
pb_Yaw <- rep(NA, length(kinFiles_pb))
for (i in 1:length(kinFiles_pb)) {pb_Yaw[i] <- kinFiles_pb[[i]]["Yaw"]}
pb_Yaw_Combined <- as.data.frame(do.call("rbind", pb_Yaw))
names(pb_Yaw_Combined) <- percentStance
pb_Yaw_Combined$filename <- unlist(pb_kin_trial)
pb_Yaw_Combined$appendage <- unlist(pb_kin_group)
pb_Yaw_Combined$species <- "pb"


## pb - Abduction / Adduction angle 
pb_AbdAdd <- rep(NA, length(kinFiles_pb))
for (i in 1:length(kinFiles_pb)) {pb_AbdAdd[i] <- kinFiles_pb[[i]]["FemHZAng"]}
pb_AbdAdd_Combined <- as.data.frame(do.call("rbind", pb_AbdAdd))
names(pb_AbdAdd_Combined) <- percentStance
pb_AbdAdd_Combined$filename <- unlist(pb_kin_trial)
pb_AbdAdd_Combined$appendage <- unlist(pb_kin_group)
pb_AbdAdd_Combined$species <- "pb"


## pb - Protaction / Retraction angle 
pb_ProRet <- rep(NA, length(kinFiles_pb))
for (i in 1:length(kinFiles_pb)) {pb_ProRet[i] <- kinFiles_pb[[i]]["FemTVAng"]}
pb_ProRet_Combined <- as.data.frame(do.call("rbind", pb_ProRet))

names(pb_ProRet_Combined) <- percentStance
pb_ProRet_Combined$filename <- unlist(pb_kin_trial)
pb_ProRet_Combined$appendage <- unlist(pb_kin_group)
pb_ProRet_Combined$species <- "pb"


## pb - Protraction / Retraction correction
#some values of pb_ProRet_Combined seem to have spontaneously changed sign 
#so I changed those values back to normal. I also subtracted 90 degrees
#to change 0 degrees to being perpendicular to the torso, which we were already doing

# Plotting the original data
orig_pbProRet_Combined <- reshape2::melt(pb_ProRet_Combined)
ggplot(data = orig_pbProRet_Combined, aes(x = factor(variable), y = value, color = filename)) +       
       geom_line(aes(group = filename))

pb_ProRet_Combined_fixed <- pb_ProRet_Combined

pb_ProRet_Combined_fixed[3,1:101]  <- pb_ProRet_Combined_fixed[3,1:101]*-1
pb_ProRet_Combined_fixed[7,1:101]  <- pb_ProRet_Combined_fixed[7,1:101]*-1
pb_ProRet_Combined_fixed[9,1:101]  <- pb_ProRet_Combined_fixed[9,1:101]*-1
pb_ProRet_Combined_fixed[12,1:101] <- pb_ProRet_Combined_fixed[12,1:101]*-1


pb_ProRet_Combined_fixed[18,72:81]  <- pb_ProRet_Combined_fixed[18,72:81]*-1
pb_ProRet_Combined_fixed[18,84:86]  <- pb_ProRet_Combined_fixed[18,84:86]*-1
pb_ProRet_Combined_fixed[18,89:100] <- pb_ProRet_Combined_fixed[18,89:100]*-1
pb_ProRet_Combined_fixed[21,45:101] <- pb_ProRet_Combined_fixed[21,45:101]*-1
pb_ProRet_Combined_fixed[25,68:101] <- pb_ProRet_Combined_fixed[25,68:101]*-1
pb_ProRet_Combined_fixed[40,70:101] <- pb_ProRet_Combined_fixed[40,70:101]*-1

pb_ProRet_Combined_fixed[,1:101] <- pb_ProRet_Combined_fixed[,1:101]-90

orig_pbProRet_Combined_fixed <- reshape2::melt(pb_ProRet_Combined_fixed)
ggplot(data = orig_pbProRet_Combined_fixed, aes(x = factor(variable), y = value, color = filename)) +       
  geom_line(aes(group = filename))

## pb - Knee / Elbow angle 
pb_KneeAng <- rep(NA, length(kinFiles_pb))
for (i in 1:length(kinFiles_pb)) {pb_KneeAng[i] <- kinFiles_pb[[i]][1]}
pb_KneeAng_Combined <- as.data.frame(do.call("rbind", pb_KneeAng))
names(pb_KneeAng_Combined) <- percentStance
pb_KneeAng_Combined$filename <- unlist(pb_kin_trial)
pb_KneeAng_Combined$appendage <- unlist(pb_kin_group)
pb_KneeAng_Combined$species <- "pb"


## pb - Ankle / Wrist angle 
pb_AnkAng <- rep(NA, length(kinFiles_pb))
for (i in 1:length(kinFiles_pb)) {pb_AnkAng[i] <- kinFiles_pb[[i]]["AnkAng"]}
pb_AnkAng_Combined <- as.data.frame(do.call("rbind", pb_AnkAng))
names(pb_AnkAng_Combined) <- percentStance
pb_AnkAng_Combined$filename <- unlist(pb_kin_trial)
pb_AnkAng_Combined$appendage <- unlist(pb_kin_group)
pb_AnkAng_Combined$species <- "pb"

## pb - Pitch angle 
pb_Pitch <- rep(NA, length(kinFiles_pb))
for (i in 1:length(kinFiles_pb)) {pb_Pitch[i] <- kinFiles_pb[[i]]["Pitch"]}
pb_Pitch_Combined <- as.data.frame(do.call("rbind", pb_Pitch))
names(pb_Pitch_Combined) <- percentStance
pb_Pitch_Combined$filename <- unlist(pb_kin_trial)
pb_Pitch_Combined$appendage <- unlist(pb_kin_group)
pb_Pitch_Combined$species <- "pb"

#Correct by 90 degrees to change how angles relate to body position
pb_Pitch_Combined[,1:101] <- pb_Pitch_Combined[,1:101]-90


#### STEP 5b: CALCULATE KINEMATIC DATA - At pec ####

## Now read in the Periophthalmus barbarus kinematic data
## Select the *folder* containing the .KIN files (e.g., AmbystomaKine_DataToAnalyze-->salamander fore-->forceanal results)
## Kinematics were calculated from MATLAB code to expedite the results, so loading files here
at_pec_kin_path <- choose.dir(caption = "select the salamander fore - forceanal results folder")
at_pec_kin_list <- list.files(at_pec_kin_path, pattern="KIN.txt", full=TRUE)

# Naming matrices in array with trial name (e.g., af05f12)
# removing the directory info from the directory path
# this duplicates some of the steps taken before, but it helps to make sure that we correctly match the kinematic files with the correct trials
at_pec_kin_filenames <- basename(at_pec_kin_list)
at_pec_kin_trial <- list(substring(at_pec_kin_filenames, 1, 7))

# Identifying group
# e.g., Pec
at_pec_kin_group <- list(substring(at_pec_kin_filenames, 9, 11))


# Getting data
kinFiles_at_pec <- array(lapply(at_pec_kin_list, read.delim, header=T), dimnames=at_pec_kin_trial)

## at_pec - Yaw angle 
at_pec_Yaw <- rep(NA, length(kinFiles_at_pec))
for (i in 1:length(kinFiles_at_pec)) {at_pec_Yaw[i] <- kinFiles_at_pec[[i]]["Yaw"]}
at_pec_Yaw_Combined <- as.data.frame(do.call("rbind", at_pec_Yaw))
names(at_pec_Yaw_Combined) <- percentStance
at_pec_Yaw_Combined$filename <- unlist(at_pec_kin_trial)
at_pec_Yaw_Combined$appendage <- unlist(at_pec_kin_group)
at_pec_Yaw_Combined$species <- "at"


## at_pec - Abduction / Adduction angle 
at_pec_AbdAdd <- rep(NA, length(kinFiles_at_pec))
for (i in 1:length(kinFiles_at_pec)) {at_pec_AbdAdd[i] <- kinFiles_at_pec[[i]]["FemHZAng"]}
at_pec_AbdAdd_Combined <- as.data.frame(do.call("rbind", at_pec_AbdAdd))
names(at_pec_AbdAdd_Combined) <- percentStance
at_pec_AbdAdd_Combined$filename <- unlist(at_pec_kin_trial)
at_pec_AbdAdd_Combined$appendage <- unlist(at_pec_kin_group)
at_pec_AbdAdd_Combined$species <- "at"


## at_pec - Protaction / Retraction angle 
at_pec_ProRet <- rep(NA, length(kinFiles_at_pec))
for (i in 1:length(kinFiles_at_pec)) {at_pec_ProRet[i] <- kinFiles_at_pec[[i]]["FemTVAng"]}
at_pec_ProRet_Combined <- as.data.frame(do.call("rbind", at_pec_ProRet))
names(at_pec_ProRet_Combined) <- percentStance
at_pec_ProRet_Combined$filename <- unlist(at_pec_kin_trial)
at_pec_ProRet_Combined$appendage <- unlist(at_pec_kin_group)
at_pec_ProRet_Combined$species <- "at"

## at_pec - Protaction / Retraction angle - corrected by 90 degrees
at_pec_ProRet_Combined_Corr <- at_pec_ProRet_Combined
at_pec_ProRet_Combined_Corr[,1:101] <- at_pec_ProRet_Combined[,1:101]-90


## at_pec - Knee / Elbow angle 
at_pec_KneeAng <- rep(NA, length(kinFiles_at_pec))
for (i in 1:length(kinFiles_at_pec)) {at_pec_KneeAng[i] <- kinFiles_at_pec[[i]][1]}
at_pec_KneeAng_Combined <- as.data.frame(do.call("rbind", at_pec_KneeAng))
names(at_pec_KneeAng_Combined) <- percentStance
at_pec_KneeAng_Combined$filename <- unlist(at_pec_kin_trial)
at_pec_KneeAng_Combined$appendage <- unlist(at_pec_kin_group)
at_pec_KneeAng_Combined$species <- "at"


## at_pec - Ankle / Wrist angle 
at_pec_AnkAng <- rep(NA, length(kinFiles_at_pec))
for (i in 1:length(kinFiles_at_pec)) {at_pec_AnkAng[i] <- kinFiles_at_pec[[i]]["AnkAng"]}
at_pec_AnkAng_Combined <- as.data.frame(do.call("rbind", at_pec_AnkAng))
names(at_pec_AnkAng_Combined) <- percentStance
at_pec_AnkAng_Combined$filename <- unlist(at_pec_kin_trial)
at_pec_AnkAng_Combined$appendage <- unlist(at_pec_kin_group)
at_pec_AnkAng_Combined$species <- "at"

## at_pec - Pitch angle 
at_pec_Pitch <- rep(NA, length(kinFiles_at_pec))
for (i in 1:length(kinFiles_at_pec)) {at_pec_Pitch[i] <- kinFiles_at_pec[[i]]["Pitch"]}
at_pec_Pitch_Combined <- as.data.frame(do.call("rbind", at_pec_Pitch))
names(at_pec_Pitch_Combined) <- percentStance
at_pec_Pitch_Combined$filename <- unlist(at_pec_kin_trial)
at_pec_Pitch_Combined$appendage <- unlist(at_pec_kin_group)
at_pec_Pitch_Combined$species <- "at"

#Correct by 90 degrees to change how angles relate to body position
at_pec_Pitch_Combined[,1:101] <- at_pec_Pitch_Combined[,1:101]-90




#### STEP 5c: CALCULATE KINEMATIC DATA - At pel ####

## Now read in the Periophthalmus barbarus kinematic data
## Select the *folder* containing the .KIN files (e.g., AmbystomaKine_DataToAnalyze-->salamander hind-->forceanal results)
## Kinematics were calculated from MATLAB code to expedite the results, so loading files here
at_pel_kin_path <- choose.dir(caption = "select the salamander hind - forceanal results folder")
at_pel_kin_list <- list.files(at_pel_kin_path, pattern="KIN.txt", full=TRUE)

# Naming matrices in array with trial name (e.g., af05f12)
# removing the directory info from the directory path
# this duplicates some of the steps taken before, but it helps to make sure that we correctly match the kinematic files with the correct trials
at_pel_kin_filenames <- basename(at_pel_kin_list)
at_pel_kin_trial <- list(substring(at_pel_kin_filenames, 1, 7))

# Identifying group
# e.g., pel
at_pel_kin_group <- list(substring(at_pel_kin_filenames, 9, 11))


# Getting data
kinFiles_at_pel <- array(lapply(at_pel_kin_list, read.delim, header=T), dimnames=at_pel_kin_trial)

## at_pel - Yaw angle 
at_pel_Yaw <- rep(NA, length(kinFiles_at_pel))
for (i in 1:length(kinFiles_at_pel)) {at_pel_Yaw[i] <- kinFiles_at_pel[[i]]["Yaw"]}
at_pel_Yaw_Combined <- as.data.frame(do.call("rbind", at_pel_Yaw))
names(at_pel_Yaw_Combined) <- percentStance
at_pel_Yaw_Combined$filename <- unlist(at_pel_kin_trial)
at_pel_Yaw_Combined$appendage <- unlist(at_pel_kin_group)
at_pel_Yaw_Combined$species <- "at"


## at_pel - Abduction / Adduction angle 
at_pel_AbdAdd <- rep(NA, length(kinFiles_at_pel))
for (i in 1:length(kinFiles_at_pel)) {at_pel_AbdAdd[i] <- kinFiles_at_pel[[i]]["FemHZAng"]}
at_pel_AbdAdd_Combined <- as.data.frame(do.call("rbind", at_pel_AbdAdd))
names(at_pel_AbdAdd_Combined) <- percentStance
at_pel_AbdAdd_Combined$filename <- unlist(at_pel_kin_trial)
at_pel_AbdAdd_Combined$appendage <- unlist(at_pel_kin_group)
at_pel_AbdAdd_Combined$species <- "at"


## at_pel - Protaction / Retraction angle 
at_pel_ProRet <- rep(NA, length(kinFiles_at_pel))
for (i in 1:length(kinFiles_at_pel)) {at_pel_ProRet[i] <- kinFiles_at_pel[[i]]["FemTVAng"]}
at_pel_ProRet_Combined <- as.data.frame(do.call("rbind", at_pel_ProRet))
names(at_pel_ProRet_Combined) <- percentStance
at_pel_ProRet_Combined$filename <- unlist(at_pel_kin_trial)
at_pel_ProRet_Combined$appendage <- unlist(at_pel_kin_group)
at_pel_ProRet_Combined$species <- "at"

## at_pel - Protaction / Retraction angle - corrected by 90 degrees
at_pel_ProRet_Combined_fixed<- at_pel_ProRet_Combined
at_pel_ProRet_Combined_fixed[,1:101] <- at_pel_ProRet_Combined[,1:101]-90
at_pel_ProRet_Combined_fixed <- at_pel_ProRet_Combined_fixed[-c(7,48),] # why are these being removed?
# this is missing two trials

## at_pel - Knee / Elbow angle 
at_pel_KneeAng <- rep(NA, length(kinFiles_at_pel))
for (i in 1:length(kinFiles_at_pel)) {at_pel_KneeAng[i] <- kinFiles_at_pel[[i]][1]}
at_pel_KneeAng_Combined <- as.data.frame(do.call("rbind", at_pel_KneeAng))
names(at_pel_KneeAng_Combined) <- percentStance
at_pel_KneeAng_Combined$filename <- unlist(at_pel_kin_trial)
at_pel_KneeAng_Combined$appendage <- unlist(at_pel_kin_group)
at_pel_KneeAng_Combined$species <- "at"


## at_pel - Ankle / Wrist angle 
at_pel_AnkAng <- rep(NA, length(kinFiles_at_pel))
for (i in 1:length(kinFiles_at_pel)) {at_pel_AnkAng[i] <- kinFiles_at_pel[[i]]["AnkAng"]}
at_pel_AnkAng_Combined <- as.data.frame(do.call("rbind", at_pel_AnkAng))
names(at_pel_AnkAng_Combined) <- percentStance
at_pel_AnkAng_Combined$filename <- unlist(at_pel_kin_trial)
at_pel_AnkAng_Combined$appendage <- unlist(at_pel_kin_group)
at_pel_AnkAng_Combined$species <- "at"


## at_pel - Pitch angle 
at_pel_Pitch <- rep(NA, length(kinFiles_at_pel))
for (i in 1:length(kinFiles_at_pel)) {at_pel_Pitch[i] <- kinFiles_at_pel[[i]]["Pitch"]}
at_pel_Pitch_Combined <- as.data.frame(do.call("rbind", at_pel_Pitch))
names(at_pel_Pitch_Combined) <- percentStance
at_pel_Pitch_Combined$filename <- unlist(at_pel_kin_trial)
at_pel_Pitch_Combined$appendage <- unlist(at_pel_kin_group)
at_pel_Pitch_Combined$species <- "at"

#Correct by 90 degrees to change how angles relate to body position
at_pel_Pitch_Combined[,1:101] <- at_pel_Pitch_Combined[,1:101]-90



#### STEP 6a: SUMMARIZING DATA - Pb ####

##  pb - Yaw 
pb_Yaw_Mean  <- sapply(pb_Yaw_Combined[,c(1:101)], FUN=function(x) mean(x, na.rm=TRUE))
pb_Yaw_SE  <- sapply(pb_Yaw_Combined[,c(1:101)], FUN=function(x) se(x, na.rm=TRUE))
pb_Yaw_Mean_SE <- data.frame(pb_Yaw_Mean[Stance5], pb_Yaw_SE[Stance5], seq(0,100,5))
names(pb_Yaw_Mean_SE) <- c("mean", "SE", "stance")
pb_Yaw_Mean_SE$type <- "pec"
pb_Yaw_Mean_SE$species <- "pb"

##  pb - AbdAdd 
pb_AbdAdd_Mean  <- sapply(pb_AbdAdd_Combined[,c(1:101)], FUN=function(x) mean(x, na.rm=TRUE))
pb_AbdAdd_SE  <- sapply(pb_AbdAdd_Combined[,c(1:101)], FUN=function(x) se(x, na.rm=TRUE))
pb_AbdAdd_Mean_SE <- data.frame(pb_AbdAdd_Mean[Stance5], pb_AbdAdd_SE[Stance5], seq(0,100,5))
names(pb_AbdAdd_Mean_SE) <- c("mean", "SE", "stance")
pb_AbdAdd_Mean_SE$type <- "pec"
pb_AbdAdd_Mean_SE$species <- "pb"

##  pb - ProRet
pb_ProRet_Mean  <- sapply(pb_ProRet_Combined[,c(1:101)], FUN=function(x) mean(x, na.rm=TRUE))
pb_ProRet_SE  <- sapply(pb_ProRet_Combined[,c(1:101)], FUN=function(x) se(x, na.rm=TRUE))
pb_ProRet_Mean_SE <- data.frame(pb_ProRet_Mean[Stance5], pb_ProRet_SE[Stance5], seq(0,100,5))
names(pb_ProRet_Mean_SE) <- c("mean", "SE", "stance")
pb_ProRet_Mean_SE$type <- "pec"
pb_ProRet_Mean_SE$species <- "pb"

##  pb - ProRet - corrected
pb_ProRet_Corr_Mean  <- sapply(pb_ProRet_Combined_fixed[,c(1:101)], FUN=function(x) mean(x, na.rm=TRUE))
pb_ProRet_Corr_SE  <- sapply(pb_ProRet_Combined_fixed[,c(1:101)], FUN=function(x) se(x, na.rm=TRUE))
pb_ProRet_Corr_Mean_SE <- data.frame(pb_ProRet_Corr_Mean[Stance5], pb_ProRet_Corr_SE[Stance5], seq(0,100,5))
names(pb_ProRet_Corr_Mean_SE) <- c("mean", "SE", "stance")
pb_ProRet_Corr_Mean_SE$type <- "pec"
pb_ProRet_Corr_Mean_SE$species <- "pb"

##  pb - KneeAng 
pb_KneeAng_Mean  <- sapply(pb_KneeAng_Combined[,c(1:101)], FUN=function(x) mean(x, na.rm=TRUE))
pb_KneeAng_SE  <- sapply(pb_KneeAng_Combined[,c(1:101)], FUN=function(x) se(x, na.rm=TRUE))
pb_KneeAng_Mean_SE <- data.frame(pb_KneeAng_Mean[Stance5], pb_KneeAng_SE[Stance5], seq(0,100,5))
names(pb_KneeAng_Mean_SE) <- c("mean", "SE", "stance")
pb_KneeAng_Mean_SE$type <- "pec"
pb_KneeAng_Mean_SE$species <- "pb"

##  pb - AnkAng 
pb_AnkAng_Mean  <- sapply(pb_AnkAng_Combined[,c(1:101)], FUN=function(x) mean(x, na.rm=TRUE))
pb_AnkAng_SE  <- sapply(pb_AnkAng_Combined[,c(1:101)], FUN=function(x) se(x, na.rm=TRUE))
pb_AnkAng_Mean_SE <- data.frame(pb_AnkAng_Mean[Stance5], pb_AnkAng_SE[Stance5], seq(0,100,5))
names(pb_AnkAng_Mean_SE) <- c("mean", "SE", "stance")
pb_AnkAng_Mean_SE$type <- "pec"
pb_AnkAng_Mean_SE$species <- "pb"

##  pb - Pitch
pb_Pitch_Mean  <- sapply(pb_Pitch_Combined[,c(1:101)], FUN=function(x) mean(x, na.rm=TRUE))
pb_Pitch_SE  <- sapply(pb_Pitch_Combined[,c(1:101)], FUN=function(x) se(x, na.rm=TRUE))
pb_Pitch_Mean_SE <- data.frame(pb_Pitch_Mean[Stance5], pb_Pitch_SE[Stance5], seq(0,100,5))
names(pb_Pitch_Mean_SE) <- c("mean", "SE", "stance")
pb_Pitch_Mean_SE$type <- "pec"
pb_Pitch_Mean_SE$species <- "pb"



#### STEP 6b: SUMMARIZING DATA - At pec ####

##  at_pec - Yaw 
at_pec_Yaw_Mean  <- sapply(at_pec_Yaw_Combined[,c(1:101)], FUN=function(x) mean(x, na.rm=TRUE))
at_pec_Yaw_SE  <- sapply(at_pec_Yaw_Combined[,c(1:101)], FUN=function(x) se(x, na.rm=TRUE))
at_pec_Yaw_Mean_SE <- data.frame(at_pec_Yaw_Mean[Stance5], at_pec_Yaw_SE[Stance5], seq(0,100,5))
names(at_pec_Yaw_Mean_SE) <- c("mean", "SE", "stance")
at_pec_Yaw_Mean_SE$type <- "pec"
at_pec_Yaw_Mean_SE$species <- "at"

##  at_pec - AbdAdd 
at_pec_AbdAdd_Mean  <- sapply(at_pec_AbdAdd_Combined[,c(1:101)], FUN=function(x) mean(x, na.rm=TRUE))
at_pec_AbdAdd_SE  <- sapply(at_pec_AbdAdd_Combined[,c(1:101)], FUN=function(x) se(x, na.rm=TRUE))
at_pec_AbdAdd_Mean_SE <- data.frame(at_pec_AbdAdd_Mean[Stance5], at_pec_AbdAdd_SE[Stance5], seq(0,100,5))
names(at_pec_AbdAdd_Mean_SE) <- c("mean", "SE", "stance")
at_pec_AbdAdd_Mean_SE$type <- "pec"
at_pec_AbdAdd_Mean_SE$species <- "at"

##  at_pec - ProRet
at_pec_ProRet_Mean  <- sapply(at_pec_ProRet_Combined[,c(1:101)], FUN=function(x) mean(x, na.rm=TRUE))
at_pec_ProRet_SE  <- sapply(at_pec_ProRet_Combined[,c(1:101)], FUN=function(x) se(x, na.rm=TRUE))
at_pec_ProRet_Mean_SE <- data.frame(at_pec_ProRet_Mean[Stance5], at_pec_ProRet_SE[Stance5], seq(0,100,5))
names(at_pec_ProRet_Mean_SE) <- c("mean", "SE", "stance")
at_pec_ProRet_Mean_SE$type <- "pec"
at_pec_ProRet_Mean_SE$species <- "at"

##  at_pec - ProRet - corrected
at_pec_ProRet_Corr_Mean  <- sapply(at_pec_ProRet_Combined_Corr[,c(1:101)], FUN=function(x) mean(x, na.rm=TRUE))
at_pec_ProRet_Corr_SE  <- sapply(at_pec_ProRet_Combined_Corr[,c(1:101)], FUN=function(x) se(x, na.rm=TRUE))
at_pec_ProRet_Corr_Mean_SE <- data.frame(at_pec_ProRet_Corr_Mean[Stance5], at_pec_ProRet_Corr_SE[Stance5], seq(0,100,5))
names(at_pec_ProRet_Corr_Mean_SE) <- c("mean", "SE", "stance")
at_pec_ProRet_Corr_Mean_SE$type <- "pec"
at_pec_ProRet_Corr_Mean_SE$species <- "at"

##  at_pec - KneeAng 
at_pec_KneeAng_Mean  <- sapply(at_pec_KneeAng_Combined[,c(1:101)], FUN=function(x) mean(x, na.rm=TRUE))
at_pec_KneeAng_SE  <- sapply(at_pec_KneeAng_Combined[,c(1:101)], FUN=function(x) se(x, na.rm=TRUE))
at_pec_KneeAng_Mean_SE <- data.frame(at_pec_KneeAng_Mean[Stance5], at_pec_KneeAng_SE[Stance5], seq(0,100,5))
names(at_pec_KneeAng_Mean_SE) <- c("mean", "SE", "stance")
at_pec_KneeAng_Mean_SE$type <- "pec"
at_pec_KneeAng_Mean_SE$species <- "at"

##  at_pec - AnkAng 
at_pec_AnkAng_Mean  <- sapply(at_pec_AnkAng_Combined[,c(1:101)], FUN=function(x) mean(x, na.rm=TRUE))
at_pec_AnkAng_SE  <- sapply(at_pec_AnkAng_Combined[,c(1:101)], FUN=function(x) se(x, na.rm=TRUE))
at_pec_AnkAng_Mean_SE <- data.frame(at_pec_AnkAng_Mean[Stance5], at_pec_AnkAng_SE[Stance5], seq(0,100,5))
names(at_pec_AnkAng_Mean_SE) <- c("mean", "SE", "stance")
at_pec_AnkAng_Mean_SE$type <- "pec"
at_pec_AnkAng_Mean_SE$species <- "at"

##  at_pec - Pitch
at_pec_Pitch_Mean  <- sapply(at_pec_Pitch_Combined[,c(1:101)], FUN=function(x) mean(x, na.rm=TRUE))
at_pec_Pitch_SE  <- sapply(at_pec_Pitch_Combined[,c(1:101)], FUN=function(x) se(x, na.rm=TRUE))
at_pec_Pitch_Mean_SE <- data.frame(at_pec_Pitch_Mean[Stance5], at_pec_Pitch_SE[Stance5], seq(0,100,5))
names(at_pec_Pitch_Mean_SE) <- c("mean", "SE", "stance")
at_pec_Pitch_Mean_SE$type <- "pec"
at_pec_Pitch_Mean_SE$species <- "at"


#### STEP 6c: SUMMARIZING DATA - at_pel ####

##  at_pel - Yaw 
at_pel_Yaw_Mean  <- sapply(at_pel_Yaw_Combined[,c(1:101)], FUN=function(x) mean(x, na.rm=TRUE))
at_pel_Yaw_SE  <- sapply(at_pel_Yaw_Combined[,c(1:101)], FUN=function(x) se(x, na.rm=TRUE))
at_pel_Yaw_Mean_SE <- data.frame(at_pel_Yaw_Mean[Stance5], at_pel_Yaw_SE[Stance5], seq(0,100,5))
names(at_pel_Yaw_Mean_SE) <- c("mean", "SE", "stance")
at_pel_Yaw_Mean_SE$type <- "pel"
at_pel_Yaw_Mean_SE$species <- "at"

##  at_pel - AbdAdd 
at_pel_AbdAdd_Mean  <- sapply(at_pel_AbdAdd_Combined[,c(1:101)], FUN=function(x) mean(x, na.rm=TRUE))
at_pel_AbdAdd_SE  <- sapply(at_pel_AbdAdd_Combined[,c(1:101)], FUN=function(x) se(x, na.rm=TRUE))
at_pel_AbdAdd_Mean_SE <- data.frame(at_pel_AbdAdd_Mean[Stance5], at_pel_AbdAdd_SE[Stance5], seq(0,100,5))
names(at_pel_AbdAdd_Mean_SE) <- c("mean", "SE", "stance")
at_pel_AbdAdd_Mean_SE$type <- "pel"
at_pel_AbdAdd_Mean_SE$species <- "at"

##  at_pel - ProRet
at_pel_ProRet_Mean  <- sapply(at_pel_ProRet_Combined[,c(1:101)], FUN=function(x) mean(x, na.rm=TRUE))
at_pel_ProRet_SE  <- sapply(at_pel_ProRet_Combined[,c(1:101)], FUN=function(x) se(x, na.rm=TRUE))
at_pel_ProRet_Mean_SE <- data.frame(at_pel_ProRet_Mean[Stance5], at_pel_ProRet_SE[Stance5], seq(0,100,5))
names(at_pel_ProRet_Mean_SE) <- c("mean", "SE", "stance")
at_pel_ProRet_Mean_SE$type <- "pel"
at_pel_ProRet_Mean_SE$species <- "at"

##  at_pel - ProRet - corrected
at_pel_ProRet_Corr_Mean  <- sapply(at_pel_ProRet_Combined_fixed[,c(1:101)], FUN=function(x) mean(x, na.rm=TRUE))
at_pel_ProRet_Corr_SE  <- sapply(at_pel_ProRet_Combined_fixed[,c(1:101)], FUN=function(x) se(x, na.rm=TRUE))
at_pel_ProRet_Corr_Mean_SE <- data.frame(at_pel_ProRet_Corr_Mean[Stance5], at_pel_ProRet_Corr_SE[Stance5], seq(0,100,5))
names(at_pel_ProRet_Corr_Mean_SE) <- c("mean", "SE", "stance")
at_pel_ProRet_Corr_Mean_SE$type <- "pel"
at_pel_ProRet_Corr_Mean_SE$species <- "at"

##  at_pel - KneeAng 
at_pel_KneeAng_Mean  <- sapply(at_pel_KneeAng_Combined[,c(1:101)], FUN=function(x) mean(x, na.rm=TRUE))
at_pel_KneeAng_SE  <- sapply(at_pel_KneeAng_Combined[,c(1:101)], FUN=function(x) se(x, na.rm=TRUE))
at_pel_KneeAng_Mean_SE <- data.frame(at_pel_KneeAng_Mean[Stance5], at_pel_KneeAng_SE[Stance5], seq(0,100,5))
names(at_pel_KneeAng_Mean_SE) <- c("mean", "SE", "stance")
at_pel_KneeAng_Mean_SE$type <- "pel"
at_pel_KneeAng_Mean_SE$species <- "at"

##  at_pel - AnkAng 
at_pel_AnkAng_Mean  <- sapply(at_pel_AnkAng_Combined[,c(1:101)], FUN=function(x) mean(x, na.rm=TRUE))
at_pel_AnkAng_SE  <- sapply(at_pel_AnkAng_Combined[,c(1:101)], FUN=function(x) se(x, na.rm=TRUE))
at_pel_AnkAng_Mean_SE <- data.frame(at_pel_AnkAng_Mean[Stance5], at_pel_AnkAng_SE[Stance5], seq(0,100,5))
names(at_pel_AnkAng_Mean_SE) <- c("mean", "SE", "stance")
at_pel_AnkAng_Mean_SE$type <- "pel"
at_pel_AnkAng_Mean_SE$species <- "at"

##  at_pel - Pitch
at_pel_Pitch_Mean  <- sapply(at_pel_Pitch_Combined[,c(1:101)], FUN=function(x) mean(x, na.rm=TRUE))
at_pel_Pitch_SE  <- sapply(at_pel_Pitch_Combined[,c(1:101)], FUN=function(x) se(x, na.rm=TRUE))
at_pel_Pitch_Mean_SE <- data.frame(at_pel_Pitch_Mean[Stance5], at_pel_Pitch_SE[Stance5], seq(0,100,5))
names(at_pel_Pitch_Mean_SE) <- c("mean", "SE", "stance")
at_pel_Pitch_Mean_SE$type <- "pel"
at_pel_Pitch_Mean_SE$species <- "at"






#### STEP 7a: PLOT DATA (Pec comparisons) ####

### Combining Ambystoma and Periophthalmus data - pectoral ###

## Yaw angle
pec_Yaw_Mean_SE <- rbind(pb_Yaw_Mean_SE, at_pec_Yaw_Mean_SE)
pec_Yaw_MaxMin <- aes(ymax=pec_Yaw_Mean_SE$mean + pec_Yaw_Mean_SE$SE, ymin=pec_Yaw_Mean_SE$mean-pec_Yaw_Mean_SE$SE)


pec_Yaw_Plot <- ggplot(data=pec_Yaw_Mean_SE, aes(x=stance, y=mean, fill=species, linetype=species))+
  scale_y_continuous("\nYaw (degrees)")+
  #scale_x_continuous("\nStance (%)\n")+
  scale_x_continuous(element_blank())+
  geom_line(size=1, alpha=0.75)+
  geom_ribbon(pec_Yaw_MaxMin, alpha=0.5)+
  scale_colour_manual(name="Species:", # changing legend title
                      labels=c("A. tigrinum  ", "P. barbarus  "), # Changing legend labels
                      values=c("ivory4", "ivory4"))+
  scale_fill_manual(name="Species:", 
                    labels=c("A. tigrinum  ", "P. barbarus  "),
                    values=c("red","blue"))+
  scale_linetype_manual(name="Species:", 
                        labels=c("A. tigrinum  ", "P. barbarus  "),
                        values=c("dashed", "solid"))+
  theme(legend.title = element_text(size = 15))+
  theme(legend.title = element_text(face = "bold"))+
  theme(legend.text = element_text(size = 15))+
  theme(axis.title.x=element_text(colour="black"))+ # vjust=0 puts a little more spacing btwn the axis text and label
  theme(axis.title.y=element_text(colour='black', size = 15))+
  theme(axis.text.x=element_text(colour='black', size = 15))+
  theme(axis.text.y=element_text(colour='black', size = 15))+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+ # get rid of gridlines
  theme(panel.background=element_blank())+ # make background white
  theme(axis.line=element_line(colour="black", linetype="solid")) # put black lines for axes
  #theme(legend.position="bottom", legend.direction="horizontal")+
  #theme(plot.title=element_text(size=8))+
  #annotate("text",  x=95, y = 160, label = "Abduction", size=4)+
  #annotate("text", label = "Adduction", x = 95, y = -60, size=4)+
  #ggtitle("A \n") + theme(plot.title=element_text(hjust=0, size=15, face="bold"))

## Abduction / Adduction angle
pec_AbdAdd_Mean_SE <- rbind(pb_AbdAdd_Mean_SE, at_pec_AbdAdd_Mean_SE)
pec_AbdAdd_MaxMin <- aes(ymax=pec_AbdAdd_Mean_SE$mean + pec_AbdAdd_Mean_SE$SE, ymin=pec_AbdAdd_Mean_SE$mean-pec_AbdAdd_Mean_SE$SE)


pec_AbdAdd_Plot <- ggplot(data=pec_AbdAdd_Mean_SE, aes(x=stance, y=mean, fill=species, linetype=species))+
  scale_y_continuous("\nAbduct / Adduct (degrees)")+
  #scale_x_continuous("\nStance (%)\n")+
  scale_x_continuous(element_blank())+
  geom_line(size=1, alpha=0.75)+
  geom_ribbon(pec_AbdAdd_MaxMin, alpha=0.5)+
  scale_colour_manual(name="Species:", # changing legend title
                      labels=c("A. tigrinum  ", "P. barbarus  "), # Changing legend labels
                      values=c("ivory4", "ivory4"))+
  scale_fill_manual(name="Species:", 
                    labels=c("A. tigrinum  ", "P. barbarus  "),
                    values=c("red","blue"))+
  scale_linetype_manual(name="Species:", 
                        labels=c("A. tigrinum  ", "P. barbarus  "),
                        values=c("dashed", "solid"))+
  theme(axis.title.x=element_text(colour="black"))+ # vjust=0 puts a little more spacing btwn the axis text and label
  theme(axis.title.y=element_text(colour='black', size = 15, vjust = 0.5))+
  theme(axis.text.x=element_text(colour='black', size = 15))+
  theme(axis.text.y=element_text(colour='black', size = 15))+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+ # get rid of gridlines
  theme(panel.background=element_blank())+ # make background white
  theme(axis.line=element_line(colour="black", linetype="solid"))+ # put black lines for axes
  #theme(legend.position="bottom", legend.direction="horizontal")+
  #theme(plot.title=element_text(size=8))+
  annotate("text",  x=95, y = 50, label = "Abduction", size=5)+
  annotate("text", label = "Adduction", x = 95, y = -70, size=5)
  #ggtitle("B\n") + theme(plot.title=element_text(hjust=0, size=15, face="bold"))




## Protraction / Retraction - Corrected
pec_ProRet_Corr_Mean_SE <- rbind(pb_ProRet_Corr_Mean_SE, at_pec_ProRet_Corr_Mean_SE)
pec_ProRet_Corr_MaxMin <- aes(ymax=pec_ProRet_Corr_Mean_SE$mean + pec_ProRet_Corr_Mean_SE$SE, ymin=pec_ProRet_Corr_Mean_SE$mean-pec_ProRet_Corr_Mean_SE$SE)


pec_ProRet_Corr_Plot <- ggplot(data=pec_ProRet_Corr_Mean_SE, aes(x=stance, y=mean, fill=species, linetype=species))+
  scale_y_continuous("\nProtract / Retract (degrees)")+
  #scale_x_continuous("\nStance (%)\n")+
  scale_x_continuous(element_blank())+
  geom_line(size=1, alpha=0.75)+
  geom_ribbon(pec_ProRet_Corr_MaxMin, alpha=0.5)+
  scale_colour_manual(name="Species:", # changing legend title
                      labels=c("A. tigrinum  ", "P. barbarus  "), # Changing legend labels
                      values=c("ivory4", "ivory4"))+
  scale_fill_manual(name="Species:", 
                    labels=c("A. tigrinum  ", "P. barbarus  "),
                    values=c("red","blue"))+
  scale_linetype_manual(name="Species:", 
                        labels=c("A. tigrinum  ", "P. barbarus  "),
                        values=c("dashed", "solid"))+
  theme(axis.title.x=element_text(colour="black"))+ # vjust=0 puts a little more spacing btwn the axis text and label
  theme(axis.title.y=element_text(colour='black', size = 15, vjust = 0.5))+
  theme(axis.text.x=element_text(colour='black', size = 15))+
  theme(axis.text.y=element_text(colour='black', size = 15))+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+ # get rid of gridlines
  theme(panel.background=element_blank())+ # make background white
  theme(axis.line=element_line(colour="black", linetype="solid"))+ # put black lines for axes
  #theme(legend.position="bottom", legend.direction="horizontal")+
  #theme(plot.title=element_text(size=8))+
  annotate("text",  x=95, y = 20, label = "Protraction", size=5)+
  annotate("text", label = "Retraction", x = 95, y = -60, size=5)
 # ggtitle("C\n") + theme(plot.title=element_text(hjust=0, size=15, face="bold"))



## Knee / elbow angle
pec_KneeAng_Mean_SE <- rbind(pb_KneeAng_Mean_SE, at_pec_KneeAng_Mean_SE)
pec_KneeAng_MaxMin <- aes(ymax=pec_KneeAng_Mean_SE$mean + pec_KneeAng_Mean_SE$SE, ymin=pec_KneeAng_Mean_SE$mean-pec_KneeAng_Mean_SE$SE)


pec_KneeAng_Plot <- ggplot(data=pec_KneeAng_Mean_SE, aes(x=stance, y=mean, fill=species, linetype=species))+
  scale_y_continuous("\nElbow angle (degrees)")+
  #scale_x_continuous("\nStance (%)\n")+
  scale_x_continuous(element_blank())+
  geom_line(size=1, alpha=0.75)+
  geom_ribbon(pec_KneeAng_MaxMin, alpha=0.5)+
  scale_colour_manual(name="Species:", # changing legend title
                      labels=c("A. tigrinum  ", "P. barbarus  "), # Changing legend labels
                      values=c("ivory4", "ivory4"))+
  scale_fill_manual(name="Species:", 
                    labels=c("A. tigrinum  ", "P. barbarus  "),
                    values=c("red","blue"))+
  scale_linetype_manual(name="Species:", 
                        labels=c("A. tigrinum  ", "P. barbarus  "),
                        values=c("dashed", "solid"))+
  theme(axis.title.x=element_text(colour="black"))+ # vjust=0 puts a little more spacing btwn the axis text and label
  theme(axis.title.y=element_text(colour='black', size = 15, vjust = 1.5))+
  theme(axis.text.x=element_text(colour='black', size = 15))+
  theme(axis.text.y=element_text(colour='black', size = 15))+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+ # get rid of gridlines
  theme(panel.background=element_blank())+ # make background white
  theme(axis.line=element_line(colour="black", linetype="solid"))+ # put black lines for axes
  #theme(legend.position="bottom", legend.direction="horizontal")+
  #theme(plot.title=element_text(size=8))+
  annotate("text",  x=95, y = 160, label = "Extension", size=5)+
  annotate("text", label = "Flexion", x = 95, y = 60, size=5)
  #ggtitle("D\n") + theme(plot.title=element_text(hjust=0, size=15, face="bold"))


## Ankle / Wrist angle
pec_AnkAng_Mean_SE <- rbind(pb_AnkAng_Mean_SE, at_pec_AnkAng_Mean_SE)
pec_AnkAng_MaxMin <- aes(ymax=pec_AnkAng_Mean_SE$mean + pec_AnkAng_Mean_SE$SE, ymin=pec_AnkAng_Mean_SE$mean-pec_AnkAng_Mean_SE$SE)


pec_AnkAng_Plot <- ggplot(data=pec_AnkAng_Mean_SE, aes(x=stance, y=mean, fill=species, linetype=species))+
  scale_y_continuous("\nWrist angle (degrees)")+
  #scale_x_continuous("\nStance (%)\n")+
  scale_x_continuous(element_blank())+
  geom_line(size=1, alpha=0.75)+
  geom_ribbon(pec_AnkAng_MaxMin, alpha=0.5)+
  scale_colour_manual(name="Species:", # changing legend title
                      labels=c("A. tigrinum  ", "P. barbarus  "), # Changing legend labels
                      values=c("ivory4", "ivory4"))+
  scale_fill_manual(name="Species:", 
                    labels=c("A. tigrinum  ", "P. barbarus  "),
                    values=c("red","blue"))+
  scale_linetype_manual(name="Species:", 
                        labels=c("A. tigrinum  ", "P. barbarus  "),
                        values=c("dashed", "solid"))+
  theme(axis.title.x=element_text(colour="black"))+ # vjust=0 puts a little more spacing btwn the axis text and label
  theme(axis.title.y=element_text(colour='black', size = 15, vjust = 1.5))+
  theme(axis.text.x=element_text(colour='black', size = 15))+
  theme(axis.text.y=element_text(colour='black', size = 15))+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+ # get rid of gridlines
  theme(panel.background=element_blank())+ # make background white
  theme(axis.line=element_line(colour="black", linetype="solid"))+ # put black lines for axes
  #theme(legend.position="bottom", legend.direction="horizontal")+
  #theme(plot.title=element_text(size=8))+
  annotate("text",  x=95, y = 165, label = "Extension", size=5)+
  annotate("text", label = "Flexion", x = 95, y = 100, size=5)
  #ggtitle("E\n") + theme(plot.title=element_text(hjust=0, size=15, face="bold"))


## Pitch angle
pec_Pitch_Mean_SE <- rbind(pb_Pitch_Mean_SE, at_pec_Pitch_Mean_SE)
pec_Pitch_MaxMin <- aes(ymax=pec_Pitch_Mean_SE$mean + pec_Pitch_Mean_SE$SE, ymin=pec_Pitch_Mean_SE$mean-pec_Pitch_Mean_SE$SE)


pec_Pitch_Plot <- ggplot(data=pec_Pitch_Mean_SE, aes(x=stance, y=mean, fill=species, linetype=species))+
  scale_y_continuous("\nPitch angle (degrees)")+
  #scale_x_continuous("\nStance (%)\n")+
  scale_x_continuous(element_blank())+
  geom_line(size=1, alpha=0.75)+
  geom_ribbon(pec_Pitch_MaxMin, alpha=0.5)+
  scale_colour_manual(name="Species:", # changing legend title
                      labels=c("A. tigrinum  ", "P. barbarus  "), # Changing legend labels
                      values=c("ivory4", "ivory4"))+
  scale_fill_manual(name="Species:", 
                    labels=c("A. tigrinum  ", "P. barbarus  "),
                    values=c("red","blue"))+
  scale_linetype_manual(name="Species:", 
                        labels=c("A. tigrinum  ", "P. barbarus  "),
                        values=c("dashed", "solid"))+
  theme(axis.title.x=element_text(colour="black"))+ # vjust=0 puts a little more spacing btwn the axis text and label
  theme(axis.title.y=element_text(colour='black', size = 15))+
  theme(axis.text.x=element_text(colour='black', size = 15))+
  theme(axis.text.y=element_text(colour='black', size = 15))+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+ # get rid of gridlines
  theme(panel.background=element_blank())+ # make background white
  theme(axis.line=element_line(colour="black", linetype="solid")) # put black lines for axes
  #theme(legend.position="bottom", legend.direction="horizontal")+
  #theme(plot.title=element_text(size=8))+
  #annotate("text",  x=95, y = 160, label = "Positive", size=4)+
  #annotate("text", label = "Negative", x = 95, y = 60, size=4)+
  #ggtitle("F\n") + theme(plot.title=element_text(hjust=0, size=15, face="bold"))



#### STEP 7b: PLOT DATA (Propulsor comparisons) ####

### Combining Ambystoma hind limb and Periophthalmus pectoral fin data  ###

## Yaw angle
propulsor_Yaw_Mean_SE <- rbind(pb_Yaw_Mean_SE, at_pel_Yaw_Mean_SE)
propulsor_Yaw_MaxMin <- aes(ymax=propulsor_Yaw_Mean_SE$mean + propulsor_Yaw_Mean_SE$SE, ymin=propulsor_Yaw_Mean_SE$mean-propulsor_Yaw_Mean_SE$SE)


propulsor_Yaw_Plot <- ggplot(data=propulsor_Yaw_Mean_SE, aes(x=stance, y=mean, fill=species, linetype=species))+
  scale_y_continuous("\nYaw (degrees)")+
  #scale_x_continuous("\nStance (%)\n")+
  scale_x_continuous(element_blank())+
  geom_line(size=1, alpha=0.75)+
  geom_ribbon(propulsor_Yaw_MaxMin, alpha=0.5)+
  scale_colour_manual(name="Species:", # changing legend title
                      labels=c("A. tigrinum  ", "P. barbarus  "), # Changing legend labels
                      values=c("ivory4", "ivory4"))+
  scale_fill_manual(name="Species:", 
                    labels=c("A. tigrinum  ", "P. barbarus  "),
                    values=c("red","blue"))+
  scale_linetype_manual(name="Species:", 
                        labels=c("A. tigrinum  ", "P. barbarus  "),
                        values=c("dashed", "solid"))+
  theme(legend.title = element_text(size = 15))+
  theme(legend.title = element_text(face = "bold"))+
  theme(legend.text = element_text(size = 15))+
  theme(axis.title.x=element_text(colour="black"))+ # vjust=0 puts a little more spacing btwn the axis text and label
  theme(axis.title.y=element_text(colour='black', size = 15))+
  theme(axis.text.x=element_text(colour='black', size = 15))+
  theme(axis.text.y=element_text(colour='black', size = 15))+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+ # get rid of gridlines
  theme(panel.background=element_blank())+ # make background white
  theme(axis.line=element_line(colour="black", linetype="solid")) # put black lines for axes
  #theme(legend.position="bottom", legend.direction="horizontal")+
  #theme(plot.title=element_text(size=8))+
  #ggtitle("A \n") + theme(plot.title=element_text(hjust=0, size=15, face="bold"))




## Abduction / Adduction angle
propulsor_AbdAdd_Mean_SE <- rbind(pb_AbdAdd_Mean_SE, at_pel_AbdAdd_Mean_SE)
propulsor_AbdAdd_MaxMin <- aes(ymax=propulsor_AbdAdd_Mean_SE$mean + propulsor_AbdAdd_Mean_SE$SE, ymin=propulsor_AbdAdd_Mean_SE$mean-propulsor_AbdAdd_Mean_SE$SE)


propulsor_AbdAdd_Plot <- ggplot(data=propulsor_AbdAdd_Mean_SE, aes(x=stance, y=mean, fill=species, linetype=species))+
  scale_y_continuous("\nAbduct / Adduct (degrees)")+
  #scale_x_continuous("\nStance (%)\n")+
  scale_x_continuous(element_blank())+
  geom_line(size=1, alpha=0.75)+
  geom_ribbon(propulsor_AbdAdd_MaxMin, alpha=0.5)+
  scale_colour_manual(name="Species:", # changing legend title
                      labels=c("A. tigrinum  ", "P. barbarus  "), # Changing legend labels
                      values=c("ivory4", "ivory4"))+
  scale_fill_manual(name="Species:", 
                    labels=c("A. tigrinum  ", "P. barbarus  "),
                    values=c("red","blue"))+
  scale_linetype_manual(name="Species:", 
                        labels=c("A. tigrinum  ", "P. barbarus  "),
                        values=c("dashed", "solid"))+
  theme(axis.title.x=element_text(colour="black"))+ # vjust=0 puts a little more spacing btwn the axis text and label
  theme(axis.title.y=element_text(colour='black', size = 15, vjust = 0.5))+
  theme(axis.text.x=element_text(colour='black', size = 15))+
  theme(axis.text.y=element_text(colour='black', size = 15))+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+ # get rid of gridlines
  theme(panel.background=element_blank())+ # make background white
  theme(axis.line=element_line(colour="black", linetype="solid"))+ # put black lines for axes
  #theme(legend.position="bottom", legend.direction="horizontal")+
  #theme(plot.title=element_text(size=8))+
  annotate("text",  x=95, y = 0, label = "Abduction", size=5)+
  annotate("text", label = "Adduction", x = 95, y = -60, size=5)
  #ggtitle("B\n") + theme(plot.title=element_text(hjust=0, size=15, face="bold"))




## Protraction / Retraction - Corrected
propulsor_ProRet_Corr_Mean_SE <- rbind(pb_ProRet_Corr_Mean_SE, at_pel_ProRet_Corr_Mean_SE)
propulsor_ProRet_Corr_MaxMin <- aes(ymax=propulsor_ProRet_Corr_Mean_SE$mean + propulsor_ProRet_Corr_Mean_SE$SE, ymin=propulsor_ProRet_Corr_Mean_SE$mean-propulsor_ProRet_Corr_Mean_SE$SE)


propulsor_ProRet_Corr_Plot <- ggplot(data=propulsor_ProRet_Corr_Mean_SE, aes(x=stance, y=mean, fill=species, linetype=species))+
  scale_y_continuous("\nProtract / Retract (degrees)")+
  #scale_x_continuous("\nStance (%)\n")+
  scale_x_continuous(element_blank())+
  geom_line(size=1, alpha=0.75)+
  geom_ribbon(propulsor_ProRet_Corr_MaxMin, alpha=0.5)+
  scale_colour_manual(name="Species:", # changing legend title
                      labels=c("A. tigrinum  ", "P. barbarus  "), # Changing legend labels
                      values=c("ivory4", "ivory4"))+
  scale_fill_manual(name="Species:", 
                    labels=c("A. tigrinum  ", "P. barbarus  "),
                    values=c("red","blue"))+
  scale_linetype_manual(name="Species:", 
                        labels=c("A. tigrinum  ", "P. barbarus  "),
                        values=c("dashed", "solid"))+
  theme(axis.title.x=element_text(colour="black"))+ # vjust=0 puts a little more spacing btwn the axis text and label
  theme(axis.title.y=element_text(colour='black', size = 15, vjust = 0.5))+
  theme(axis.text.x=element_text(colour='black', size = 15))+
  theme(axis.text.y=element_text(colour='black', size = 15))+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+ # get rid of gridlines
  theme(panel.background=element_blank())+ # make background white
  theme(axis.line=element_line(colour="black", linetype="solid"))+ # put black lines for axes
  #theme(legend.position="bottom", legend.direction="horizontal")+
  #theme(plot.title=element_text(size=8))+
  annotate("text",  x=95, y = 40, label = "Protraction", size=5)+
  annotate("text", label = "Retraction", x = 95, y = -60, size=5)
  #ggtitle("C\n") + theme(plot.title=element_text(hjust=0, size=15, face="bold"))



## Knee / elbow angle
propulsor_KneeAng_Mean_SE <- rbind(pb_KneeAng_Mean_SE, at_pel_KneeAng_Mean_SE)
propulsor_KneeAng_MaxMin <- aes(ymax=propulsor_KneeAng_Mean_SE$mean + propulsor_KneeAng_Mean_SE$SE, ymin=propulsor_KneeAng_Mean_SE$mean-propulsor_KneeAng_Mean_SE$SE)


propulsor_KneeAng_Plot <- ggplot(data=propulsor_KneeAng_Mean_SE, aes(x=stance, y=mean, fill=species, linetype=species))+
  scale_y_continuous("\nKnee angle (degrees)")+
  #scale_x_continuous("\nStance (%)\n")+
  scale_x_continuous(element_blank())+
  geom_line(size=1, alpha=0.75)+
  geom_ribbon(propulsor_KneeAng_MaxMin, alpha=0.5)+
  scale_colour_manual(name="Species:", # changing legend title
                      labels=c("A. tigrinum  ", "P. barbarus  "), # Changing legend labels
                      values=c("ivory4", "ivory4"))+
  scale_fill_manual(name="Species:", 
                    labels=c("A. tigrinum  ", "P. barbarus  "),
                    values=c("red","blue"))+
  scale_linetype_manual(name="Species:", 
                        labels=c("A. tigrinum  ", "P. barbarus  "),
                        values=c("dashed", "solid"))+
  theme(axis.title.x=element_text(colour="black"))+ # vjust=0 puts a little more spacing btwn the axis text and label
  theme(axis.title.y=element_text(colour='black', size = 15, vjust = 1.5))+
  theme(axis.text.x=element_text(colour='black', size = 15))+
  theme(axis.text.y=element_text(colour='black', size = 15))+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+ # get rid of gridlines
  theme(panel.background=element_blank())+ # make background white
  theme(axis.line=element_line(colour="black", linetype="solid"))+ # put black lines for axes
  #theme(legend.position="bottom", legend.direction="horizontal")+
  #theme(plot.title=element_text(size=8))+
  annotate("text",  x=95, y = 160, label = "Extension", size=5)+
  annotate("text", label = "Flexion", x = 95, y = 80, size=5)
  #ggtitle("D\n") + theme(plot.title=element_text(hjust=0, size=15, face="bold"))


## Ankle / Wrist angle
propulsor_AnkAng_Mean_SE <- rbind(pb_AnkAng_Mean_SE, at_pel_AnkAng_Mean_SE)
propulsor_AnkAng_MaxMin <- aes(ymax=propulsor_AnkAng_Mean_SE$mean + propulsor_AnkAng_Mean_SE$SE, ymin=propulsor_AnkAng_Mean_SE$mean-propulsor_AnkAng_Mean_SE$SE)


propulsor_AnkAng_Plot <- ggplot(data=propulsor_AnkAng_Mean_SE, aes(x=stance, y=mean, fill=species, linetype=species))+
  scale_y_continuous("\nAnkle angle (degrees)")+
  #scale_x_continuous("\nStance (%)\n")+
  scale_x_continuous(element_blank())+
  geom_line(size=1, alpha=0.75)+
  geom_ribbon(propulsor_AnkAng_MaxMin, alpha=0.5)+
  scale_colour_manual(name="Species:", # changing legend title
                      labels=c("A. tigrinum  ", "P. barbarus  "), # Changing legend labels
                      values=c("ivory4", "ivory4"))+
  scale_fill_manual(name="Species:", 
                    labels=c("A. tigrinum  ", "P. barbarus  "),
                    values=c("red","blue"))+
  scale_linetype_manual(name="Species:", 
                        labels=c("A. tigrinum  ", "P. barbarus  "),
                        values=c("dashed", "solid"))+
  theme(axis.title.x=element_text(colour="black"))+ # vjust=0 puts a little more spacing btwn the axis text and label
  theme(axis.title.y=element_text(colour='black', size = 15, vjust = 1.5))+
  theme(axis.text.x=element_text(colour='black', size = 15))+
  theme(axis.text.y=element_text(colour='black', size = 15))+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+ # get rid of gridlines
  theme(panel.background=element_blank())+ # make background white
  theme(axis.line=element_line(colour="black", linetype="solid"))+ # put black lines for axes
  #theme(legend.position="bottom", legend.direction="horizontal")+
  #theme(plot.title=element_text(size=8))+
  annotate("text",  x=95, y = 165, label = "Extension", size=5)+
  annotate("text", label = "Flexion", x = 95, y = 60, size=5)
  #ggtitle("E\n") + theme(plot.title=element_text(hjust=0, size=15, face="bold"))


## Pitch angle
propulsor_Pitch_Mean_SE <- rbind(pb_Pitch_Mean_SE, at_pel_Pitch_Mean_SE)
propulsor_Pitch_MaxMin <- aes(ymax=propulsor_Pitch_Mean_SE$mean + propulsor_Pitch_Mean_SE$SE, ymin=propulsor_Pitch_Mean_SE$mean-propulsor_Pitch_Mean_SE$SE)


propulsor_Pitch_Plot <- ggplot(data=propulsor_Pitch_Mean_SE, aes(x=stance, y=mean, fill=species, linetype=species))+
  scale_y_continuous("\nPitch angle (degrees)")+
  #scale_x_continuous("\nStance (%)\n")+
  scale_x_continuous(element_blank())+
  geom_line(size=1, alpha=0.75)+
  geom_ribbon(propulsor_Pitch_MaxMin, alpha=0.5)+
  scale_colour_manual(name="Species:", # changing legend title
                      labels=c("A. tigrinum  ", "P. barbarus  "), # Changing legend labels
                      values=c("ivory4", "ivory4"))+
  scale_fill_manual(name="Species:", 
                    labels=c("A. tigrinum  ", "P. barbarus  "),
                    values=c("red","blue"))+
  scale_linetype_manual(name="Species:", 
                        labels=c("A. tigrinum  ", "P. barbarus  "),
                        values=c("dashed", "solid"))+
  theme(axis.title.x=element_text(colour="black"))+ # vjust=0 puts a little more spacing btwn the axis text and label
  theme(axis.title.y=element_text(colour='black', size = 15))+
  theme(axis.text.x=element_text(colour='black', size = 15))+
  theme(axis.text.y=element_text(colour='black', size = 15))+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+ # get rid of gridlines
  theme(panel.background=element_blank())+ # make background white
  theme(axis.line=element_line(colour="black", linetype="solid")) # put black lines for axes
  #theme(legend.position="bottom", legend.direction="horizontal")+
  #theme(plot.title=element_text(size=8))+
  #annotate("text",  x=95, y = 160, label = "Positive", size=4)+
  #annotate("text", label = "Negative", x = 95, y = 60, size=4)+
  #ggtitle("F\n") + theme(plot.title=element_text(hjust=0, size=15, face="bold"))



#### STEP 7c: PLOT DATA (Ambystoma comparisons) ####

### Combining Ambystoma hind limb and forelimb data  ###

## Yaw angle
limb_Yaw_Mean_SE <- rbind(at_pec_Yaw_Mean_SE, at_pel_Yaw_Mean_SE)
limb_Yaw_MaxMin <- aes(ymax=limb_Yaw_Mean_SE$mean + limb_Yaw_Mean_SE$SE, ymin=limb_Yaw_Mean_SE$mean-limb_Yaw_Mean_SE$SE)


limb_Yaw_Plot <- ggplot(data=limb_Yaw_Mean_SE, aes(x=stance, y=mean, fill=type, linetype=type))+
  scale_y_continuous("\nYaw (degrees)")+
  #scale_x_continuous("\nStance (%)\n")+
  scale_x_continuous(element_blank())+
  geom_line(size=1, alpha=0.75)+
  geom_ribbon(limb_Yaw_MaxMin, alpha=0.5)+
  scale_colour_manual(name="limb:", # changing legend title
                      labels=c("forelimb  ", "hind limb  "), # Changing legend labels
                      values=c("ivory4", "ivory4"))+
  scale_fill_manual(name="limb:", 
                    labels=c("forelimb  ", "hind limb  "),
                    values=c("red","blue"))+
  scale_linetype_manual(name="limb:", 
                        labels=c("forelimb  ", "hind limb  "),
                        values=c("dashed", "solid"))+
  theme(legend.title = element_text(size = 15))+
  theme(legend.title = element_text(face = "bold"))+
  theme(legend.text = element_text(size = 15))+
  theme(axis.title.x=element_text(colour="black"))+ # vjust=0 puts a little more spacing btwn the axis text and label
  theme(axis.title.y=element_text(colour='black', size = 15))+
  theme(axis.text.x=element_text(colour='black', size = 15))+
  theme(axis.text.y=element_text(colour='black', size = 15))+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+ # get rid of gridlines
  theme(panel.background=element_blank())+ # make background white
  theme(axis.line=element_line(colour="black", linetype="solid")) # put black lines for axes
  #theme(legend.position="bottom", legend.direction="horizontal")+
  #theme(plot.title=element_text(size=8))+
  #annotate("text",  x=95, y = 160, label = "Abduction", size=4)+
  #annotate("text", label = "Adduction", x = 95, y = -60, size=4)+
  #ggtitle("A \n") + theme(plot.title=element_text(hjust=0, size=15, face="bold"))




## Pitch angle
limb_Pitch_Mean_SE <- rbind(at_pec_Pitch_Mean_SE, at_pel_Pitch_Mean_SE)
limb_Pitch_MaxMin <- aes(ymax=limb_Pitch_Mean_SE$mean + limb_Pitch_Mean_SE$SE, ymin=limb_Pitch_Mean_SE$mean-limb_Pitch_Mean_SE$SE)


limb_Pitch_Plot <- ggplot(data=limb_Pitch_Mean_SE, aes(x=stance, y=mean, fill=type, linetype=type))+
  scale_y_continuous("\nPitch angle (degrees)")+
  #scale_x_continuous("\nStance (%)\n")+
  scale_x_continuous(element_blank())+
  geom_line(size=1, alpha=0.75)+
  geom_ribbon(limb_Pitch_MaxMin, alpha=0.5)+
  scale_colour_manual(name="limb:", # changing legend title
                      labels=c("forelimb  ", "hind limb  "), # Changing legend labels
                      values=c("ivory4", "ivory4"))+
  scale_fill_manual(name="limb:", 
                    labels=c("forelimb  ", "hind limb  "),
                    values=c("red","blue"))+
  scale_linetype_manual(name="limb:", 
                        labels=c("forelimb  ", "hind limb  "),
                        values=c("dashed", "solid"))+
  theme(axis.title.x=element_text(colour="black"))+ # vjust=0 puts a little more spacing btwn the axis text and label
  theme(axis.title.y=element_text(colour='black', size = 15))+
  theme(axis.text.x=element_text(colour='black', size = 15))+
  theme(axis.text.y=element_text(colour='black', size = 15))+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+ # get rid of gridlines
  theme(panel.background=element_blank())+ # make background white
  theme(axis.line=element_line(colour="black", linetype="solid")) # put black lines for axes
  #theme(legend.position="bottom", legend.direction="horizontal")+
  #theme(plot.title=element_text(size=8))+
  #annotate("text",  x=95, y = 160, label = "Positive", size=4)+
  #annotate("text", label = "Negative", x = 95, y = 60, size=4)+
  #ggtitle("F\n") + theme(plot.title=element_text(hjust=0, size=15, face="bold"))


#### STEP 8a: GROUP PLOTS - Pec ####

xTitle <- textGrob("Percent of Stance", 
                   gp=gpar(fontface="bold", fontsize=14))


# extract the legend from one of the plots
pec_legend <- get_legend(
  pec_Yaw_Plot + guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "left") # align the legend options vertically
)


pec_PlotCompare <- plot_grid(pec_AbdAdd_Plot + theme(legend.position="none"), 
          pec_ProRet_Corr_Plot + theme(legend.position="none"), 
          pec_KneeAng_Plot + theme(legend.position="none"),
          pec_AnkAng_Plot + theme(legend.position="none"),
          pec_Pitch_Plot + theme(legend.position="none"),
          pec_Yaw_Plot + theme(legend.position="none"),
          
          #pec_legend,
          ncol = 2, align = "v",  labels = "AUTO",
          label_size = 20)


# Add common x-axis label to plots
grid.arrange(arrangeGrob(pec_PlotCompare, bottom = xTitle, right = pec_legend))


#### STEP 8b: GROUP PLOTS - Propulsors ####

xTitle <- textGrob("Percent of Stance", 
                   gp=gpar(fontface="bold", fontsize=14))


# extract the legend from one of the plots
propulsor_legend <- get_legend(
  propulsor_Yaw_Plot + guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "left") # align the legend options vertically
)


propulsor_PlotCompare <- plot_grid(propulsor_AbdAdd_Plot + theme(legend.position="none"), 
                             propulsor_ProRet_Corr_Plot + theme(legend.position="none"), 
                             propulsor_KneeAng_Plot + theme(legend.position="none"), 
                             propulsor_AnkAng_Plot + theme(legend.position="none"),
                             propulsor_Pitch_Plot + theme(legend.position="none"), 
                             propulsor_Yaw_Plot + theme(legend.position="none"),
                             
                             #propulsor_legend,
                             ncol = 2, align = "v",  labels = "AUTO", 
                             label_size = 20)


# Add common x-axis label to plots
grid.arrange(arrangeGrob(propulsor_PlotCompare, bottom = xTitle, right = propulsor_legend))


### lmer_prep() ####
lmer_prep <- function(appendage1, appendage2) {
  #merge dataframes by column. sort set to F since dataframes are already sorted
  jointAngMerged <- merge(appendage1, appendage2, all = T, sort = F) 
  #jointAngMerged <- rbind(do.call(rbind, appendage1), do.call(rbind, appendage2))
  
  #save individual ID's (eg. at01, pb05) in new column 'Ind'
  jointAngMerged$Ind <- substring(jointAngMerged$filename, 1, 4)
  
  #save max/min/mean of each trial (row) and save index (%stance)
  #of each max/min in new columns
  jointAngMerged$Max  <- apply(jointAngMerged[, 1:101], 1, max)
  jointAngMerged$Min  <- apply(jointAngMerged[, 1:101], 1, min)
  jointAngMerged$Mean <-  apply(jointAngMerged[, 1:101], 1, mean)
  jointAngMerged$Tmax <- apply(jointAngMerged[, 1:101], 1, which.max)
  jointAngMerged$Tmin <- apply(jointAngMerged[, 1:101], 1, which.min)

  return(jointAngMerged)
}


####LMER_CALC()####

#For Troubleshooting lmer_calc():
#appendage1 <-  pb_KneeAng_Combined
#appendage2 <-  at_pec_KneeAng_Combined

lmer_calc <- function(appendage1, appendage2, fixedEffect = "species") {
  #merge dataframes by column. sort set to F since dataframes are already sorted
  jointAngMerged <- merge(appendage1, appendage2, all = T, sort = F) 
  #jointAngMerged <- rbind(do.call(rbind, appendage1), do.call(rbind, appendage2))
  
  #save individual ID's (eg. at01, pb05) in new column 'Ind'
  jointAngMerged$Ind <- substring(jointAngMerged$filename, 1, 4)
  
  #save max/min/mean of each trial (row) and save index (%stance)
  #of each max/min in new columns
  jointAngMerged$Max  <- apply(jointAngMerged[, 1:101], 1, max)
  jointAngMerged$Min  <- apply(jointAngMerged[, 1:101], 1, min)
  jointAngMerged$Mean <-  apply(jointAngMerged[, 1:101], 1, mean)
  jointAngMerged$Tmax <- apply(jointAngMerged[, 1:101], 1, which.max)
  jointAngMerged$Tmin <- apply(jointAngMerged[, 1:101], 1, which.min)
  #TODO: Troubleshoot singularity issues in lmer's 
    
  #run lmers for max and min joint angles, and the %stance of max and min

  #compares between species
  if(fixedEffect == "species") {
    lmer_max   <- lmer(Max ~ species + (1|Ind), data = jointAngMerged)
    lmer_min   <- lmer(Min ~ species + (1|Ind), data = jointAngMerged)
    lmer_mean  <- lmer(Mean ~species + (1|Ind), data = jointAngMerged)
    lmer_Tmax <- lmer(Tmax ~ species + (1|Ind), data = jointAngMerged)
    lmer_Tmin <- lmer(Tmin ~ species + (1|Ind), data = jointAngMerged)

    #for stats troubleshooting
    #qqnorm(residuals(lmer_max), main = "max")
    #qqnorm(residuals(lmer_min), main = "min")
    #qqnorm(residuals(lmer_mean), main = "mean")
    #qqnorm(residuals(lmer_Tmax), main = "Tmax")
    #qqnorm(residuals(lmer_Tmin), main = "Tmin")
  }

  #compares within species
  else if(fixedEffect == "appendage"){
    lmer_max   <- lmer(Max ~ appendage + (1|Ind), data = jointAngMerged)
    lmer_min   <- lmer(Min ~ appendage + (1|Ind), data = jointAngMerged)
    lmer_mean  <- lmer(Mean ~appendage + (1|Ind), data = jointAngMerged)
    lmer_Tmax <- lmer(Tmax ~ appendage + (1|Ind), data = jointAngMerged)
    lmer_Tmin <- lmer(Tmin ~ appendage + (1|Ind), data = jointAngMerged)

    #for stats troubleshooting
    #qqnorm(residuals(lmer_max), main = "max")
    #qqnorm(residuals(lmer_min), main = "min")
    #qqnorm(residuals(lmer_mean), main = "mean")
    #qqnorm(residuals(lmer_Tmax), main = "Tmax")
    #qqnorm(residuals(lmer_Tmin), main = "Tmin")
    }

  #if neither appendage or species is specified
  else {
    stop("Please specify the fixed effect: 'species' or 'appendage'")
  }

  #combine all lmer results into one list and return list
  #lmer_results <- list(lmer_max, lmer_min, lmer_mean)
  lmer_results <- list(lmer_max, lmer_min, lmer_Tmax, lmer_Tmin, lmer_mean)

  return(lmer_results)
}


#### LMER PREP ####

pbVars      <- list(pb_AbdAdd_Combined, pb_AnkAng_Combined, pb_KneeAng_Combined, pb_Pitch_Combined, 
                    pb_ProRet_Combined_fixed, pb_Yaw_Combined)
at_pec_Vars <- list(at_pec_AbdAdd_Combined, at_pec_AnkAng_Combined, at_pec_KneeAng_Combined, 
                    at_pec_Pitch_Combined, at_pec_ProRet_Combined_Corr, at_pec_Yaw_Combined)
at_pel_Vars <- list(at_pel_AbdAdd_Combined, at_pel_AnkAng_Combined, at_pel_KneeAng_Combined, 
                    at_pel_Pitch_Combined, at_pel_ProRet_Combined_fixed, at_pel_Yaw_Combined)


pb_atpec_lmer <- matrix(NA,5, length(pbVars))
pb_atpel_lmer <- matrix(NA,5, length(pbVars))
at_pecpel_lmer <- matrix(NA,5, length(pbVars))

for(i in 1:length(pbVars)){
  
  pb_atpec_lmer[(5*i-4):(5*i)] <- lmer_calc(pbVars[i], at_pec_Vars[i], "species")
  pb_atpel_lmer[(5*i-4):(5*i)] <- lmer_calc(pbVars[i], at_pel_Vars[i], "species")
  at_pecpel_lmer[(5*i-4):(5*i)] <- lmer_calc(at_pec_Vars[i], at_pel_Vars[i], "appendage")
} 

### Prepping dataframes
pb_atpec <- list()
pb_atpel <- list()
at_pecpel <- list()

for(i in 1:length(pbVars)){
  pb_atpec[[i]] <- lmer_prep(pbVars[i], at_pec_Vars[i])
  pb_atpel[[i]] <- lmer_prep(pbVars[i], at_pel_Vars[i])
  at_pecpel[[i]] <- lmer_prep(at_pec_Vars[i], at_pel_Vars[i])

}

Vars <- c("AbdAdd_Combined", "AnkAng_Combined", "KneeAng_Combined", "Pitch_Combined", 
          "ProRet_Combined_fixed", "Yaw_Combined")

names(pb_atpec_lmer) <- paste(rep(Vars, each = 5),"_", rep(c("lmer_max", "lmer_min", "lmer_Tmax", "lmer_Tmin", "lmer_mean")), sep = "")
names(pb_atpel_lmer) <- paste(rep(Vars, each = 5),"_", rep(c("lmer_max", "lmer_min", "lmer_Tmax", "lmer_Tmin", "lmer_mean")), sep = "")
names(at_pecpel_lmer) <- paste(rep(Vars, each = 5),"_", rep(c("lmer_max", "lmer_min", "lmer_Tmax", "lmer_Tmin", "lmer_mean")), sep = "")

names(pb_atpec) <- Vars
names(pb_atpel) <- Vars
names(at_pecpel) <- Vars


### LMERs - with lmer_calc ####

## Commenting out because got a lot of warnings about "number of items to replace is not a multiple of replacement length"
# pb_atpec_lmers <- matrix(NA,5, length(pbVars))
# pb_atpel_lmer <- matrix(NA,5, length(pbVars))
# at_pecpel_lmer <- matrix(NA,5, length(pbVars))

# ## Creating empty variables to dump output from for loop
# pb_atpec_lmers <- list()
# pb_atpel_lmer <- list()
# at_pecpel_lmer <- list()
# 
# for(i in 1:length(pbVars)){
# ## Commenting the following three lines out because they led to the warnings about multiple of replacement length
#   # pb_atpec_lmers[(5*i-4):(5*i)] <- lmer_calc(pbVars[i], at_pec_Vars[i], "species")
#   # pb_atpel_lmer[(5*i-4):(5*i)] <- lmer_calc(pbVars[i], at_pel_Vars[i], "species") # generates error for yaw lmers
#   # at_pecpel_lmer[(5*i-4):(5*i)] <- lmer_calc(at_pec_Vars[i], at_pel_Vars[i], "appendage")
#   pb_atpec_lmers[[i]] <- lmer_calc(pbVars[i], at_pec_Vars[i], "species")
#   pb_atpel_lmer[[i]] <- lmer_calc(pbVars[i], at_pel_Vars[i], "species")
#   at_pecpel_lmer[[i]] <- lmer_calc(at_pec_Vars[i], at_pel_Vars[i], "species")
# }
# 
# names(pb_atpec_lmers) <- paste(rep(Vars, each = 5),"_", rep(c("lmer_max", "lmer_min", "lmer_Tmax", "lmer_Tmin", "lmer_mean")), sep = "")
# names(pb_atpel_lmer) <- paste(rep(Vars, each = 5),"_", rep(c("lmer_max", "lmer_min", "lmer_Tmax", "lmer_Tmin", "lmer_mean")), sep = "")
# names(at_pecpel_lmer) <- paste(rep(Vars, each = 5),"_", rep(c("lmer_max", "lmer_min", "lmer_Tmax", "lmer_Tmin", "lmer_mean")), sep = "")

#### LMERS - manual - Pec - pb vs. at ####
## Use the output from emmeans() to prepare the table 
## use the output from pairs() if you're interested in the pairwise difference between the groups

#sets emmeans() to return all sigfigs instead of only first 3 digits (rounded)
#(EX: returned 22.3 or 154 at first, now returns 22.331 and 153.61)
emm_options(opt.digits = F)

## AbdAdd_Combined
pb_atpec_lmer_AbdAdd_Combined_Max <- lmer(Max ~ species + (1|Ind), data = pb_atpec$AbdAdd_Combined)
pb_atpec_lmer_AbdAdd_Combined_Max_emm <- emmeans(pb_atpec_lmer_AbdAdd_Combined_Max, "species")
pb_atpec_lmer_AbdAdd_Combined_Max_emm
pairs(pb_atpec_lmer_AbdAdd_Combined_Max_emm)
performance::r2_xu(pb_atpec_lmer_AbdAdd_Combined_Max) # Xu's R2 = 0.984

pb_atpec_lmer_AbdAdd_Combined_Min <- lmer(Min ~ species + (1|Ind), data = pb_atpec$AbdAdd_Combined)
pb_atpec_lmer_AbdAdd_Combined_Min_emm <- emmeans(pb_atpec_lmer_AbdAdd_Combined_Min, "species")
pb_atpec_lmer_AbdAdd_Combined_Min_emm
pairs(pb_atpec_lmer_AbdAdd_Combined_Min_emm)
performance::r2_xu(pb_atpec_lmer_AbdAdd_Combined_Min) # Xu's R2 = 0.977

pb_atpec_lmer_AbdAdd_Combined_Mean <- lmer(Mean ~ species + (1|Ind), data = pb_atpec$AbdAdd_Combined)
pb_atpec_lmer_AbdAdd_Combined_Mean_emm <- emmeans(pb_atpec_lmer_AbdAdd_Combined_Mean, "species")
pb_atpec_lmer_AbdAdd_Combined_Mean_emm
pairs(pb_atpec_lmer_AbdAdd_Combined_Mean_emm)
performance::r2_xu(pb_atpec_lmer_AbdAdd_Combined_Mean) # Xu's R2 = 0.988

pb_atpec_lmer_AbdAdd_Combined_Tmax <- lmer(Tmax ~ species + (1|Ind), data = pb_atpec$AbdAdd_Combined)
pb_atpec_lmer_AbdAdd_Combined_Tmax_emm <- emmeans(pb_atpec_lmer_AbdAdd_Combined_Tmax, "species")
pb_atpec_lmer_AbdAdd_Combined_Tmax_emm
pairs(pb_atpec_lmer_AbdAdd_Combined_Tmax)
performance::r2_xu(pb_atpec_lmer_AbdAdd_Combined_Tmax) # Xu's R2 = 0.919

pb_atpec_lmer_AbdAdd_Combined_Tmin <- lmer(Tmin ~ species + (1|Ind), data = pb_atpec$AbdAdd_Combined)
pb_atpec_lmer_AbdAdd_Combined_Tmin_emm <- emmeans(pb_atpec_lmer_AbdAdd_Combined_Tmin , "species")
pb_atpec_lmer_AbdAdd_Combined_Tmin_emm
pairs(pb_atpec_lmer_AbdAdd_Combined_Tmin_emm)
performance::r2_xu(pb_atpec_lmer_AbdAdd_Combined_Tmin) # Xu's R2 = 0.938


## AnkAng_Combined
pb_atpec_lmer_AnkAng_Combined_Max <- lmer(Max ~ species + (1|Ind), data = pb_atpec$AnkAng_Combined)
performance::r2_xu(pb_atpec_lmer_AnkAng_Combined_Max) # Xu's R2 = 0.614

pb_atpec_lmer_AnkAng_Combined_Min <- lmer(Min ~ species + (1|Ind), data = pb_atpec$AnkAng_Combined)
performance::r2_xu(pb_atpec_lmer_AnkAng_Combined_Min) # Xu's R2 = 0.311

pb_atpec_lmer_AnkAng_Combined_Mean <- lmer(Mean ~ species + (1|Ind), data = pb_atpec$AnkAng_Combined)
performance::r2_xu(pb_atpec_lmer_AnkAng_Combined_Mean) # Xu's R2 = 0.610

pb_atpec_lmer_AnkAng_Combined_Tmax <- lmer(Tmax ~ species + (1|Ind), data = pb_atpec$AnkAng_Combined)
performance::r2_xu(pb_atpec_lmer_AnkAng_Combined_Tmax) # Xu's R2 = 0.287

pb_atpec_lmer_AnkAng_Combined_Tmin <- lmer(Tmin ~ species + (1|Ind), data = pb_atpec$AnkAng_Combined)
## singularity issue
performance::r2_xu(pb_atpec_lmer_AnkAng_Combined_Tmin) # Xu's R2 = 0.536


## KneeAng_Combined
pb_atpec_lmer_KneeAng_Combined_Max <- lmer(Max ~ species + (1|Ind), data = pb_atpec$KneeAng_Combined)
performance::r2_xu(pb_atpec_lmer_KneeAng_Combined_Max) # Xu's R2 = 0.577

pb_atpec_lmer_KneeAng_Combined_Min <- lmer(Min ~ species + (1|Ind), data = pb_atpec$KneeAng_Combined)
performance::r2_xu(pb_atpec_lmer_KneeAng_Combined_Min) # Xu's R2 = 0.859

pb_atpec_lmer_KneeAng_Combined_Mean <- lmer(Mean ~ species + (1|Ind), data = pb_atpec$KneeAng_Combined)
performance::r2_xu(pb_atpec_lmer_KneeAng_Combined_Mean) # Xu's R2 = 0.875

pb_atpec_lmer_KneeAng_Combined_Tmax <- lmer(Tmax ~ species + (1|Ind), data = pb_atpec$KneeAng_Combined)
# singularity issue
performance::r2_xu(pb_atpec_lmer_KneeAng_Combined_Tmax) # Xu's R2 = 0.749

pb_atpec_lmer_KneeAng_Combined_Tmin <- lmer(Tmin ~ species + (1|Ind), data = pb_atpec$KneeAng_Combined)
# sinularity issue
performance::r2_xu(pb_atpec_lmer_KneeAng_Combined_Tmin) # Xu's R2 = 0.021


## Pitch_Combined
pb_atpec_lmer_Pitch_Combined_Max <- lmer(Max ~ species + (1|Ind), data = pb_atpec$Pitch_Combined)
performance::r2_xu(pb_atpec_lmer_Pitch_Combined_Max) # Xu's R2 = 0.851

pb_atpec_lmer_Pitch_Combined_Min <- lmer(Min ~ species + (1|Ind), data = pb_atpec$Pitch_Combined)
performance::r2_xu(pb_atpec_lmer_Pitch_Combined_Min) # Xu's R2 = 0.878

pb_atpec_lmer_Pitch_Combined_Mean <- lmer(Mean ~ species + (1|Ind), data = pb_atpec$Pitch_Combined)
performance::r2_xu(pb_atpec_lmer_Pitch_Combined_Mean) # Xu's R2 = 0.877

pb_atpec_lmer_Pitch_Combined_Tmax <- lmer(Tmax ~ species + (1|Ind), data = pb_atpec$Pitch_Combined)
performance::r2_xu(pb_atpec_lmer_Pitch_Combined_Tmax) # Xu's R2 = 0.541

pb_atpec_lmer_Pitch_Combined_Tmin <- lmer(Tmin ~ species + (1|Ind), data = pb_atpec$Pitch_Combined)
performance::r2_xu(pb_atpec_lmer_Pitch_Combined_Tmin) # Xu's R2 = 0.300


## ProRet_Combined_fixed
pb_atpec_lmer_ProRet_Combined_fixed_Max <- lmer(Max ~ species + (1|Ind), data = pb_atpec$ProRet_Combined_fixed)
performance::r2_xu(pb_atpec_lmer_ProRet_Combined_fixed_Max) # Xu's R2 = 0.717

pb_atpec_lmer_ProRet_Combined_fixed_Min <- lmer(Min ~ species + (1|Ind), data = pb_atpec$ProRet_Combined_fixed)
performance::r2_xu(pb_atpec_lmer_ProRet_Combined_fixed_Min) # Xu's R2 = 0.554

pb_atpec_lmer_ProRet_Combined_fixed_Mean <- lmer(Mean ~ species + (1|Ind), data = pb_atpec$ProRet_Combined_fixed)
performance::r2_xu(pb_atpec_lmer_ProRet_Combined_fixed_Mean) # Xu's R2 = 0.435

pb_atpec_lmer_ProRet_Combined_fixed_Tmax <- lmer(Tmax ~ species + (1|Ind), data = pb_atpec$ProRet_Combined_fixed)
## singularity issue
performance::r2_xu(pb_atpec_lmer_ProRet_Combined_fixed_Tmax) # Xu's R2 = 0.087

pb_atpec_lmer_ProRet_Combined_fixed_Tmin <- lmer(Tmin ~ species + (1|Ind), data = pb_atpec$ProRet_Combined_fixed)
performance::r2_xu(pb_atpec_lmer_ProRet_Combined_fixed_Tmin) # Xu's R2 = 0.355


## Yaw_Combined
pb_atpec_lmer_Yaw_Combined_Max <- lmer(Max ~ species + (1|Ind), data = pb_atpec$Yaw_Combined)
performance::r2_xu(pb_atpec_lmer_Yaw_Combined_Max) # Xu's R2 = 0.752

pb_atpec_lmer_Yaw_Combined_Min <- lmer(Min ~ species + (1|Ind), data = pb_atpec$Yaw_Combined)
performance::r2_xu(pb_atpec_lmer_Yaw_Combined_Min) # Xu's R2 = 0.485

pb_atpec_lmer_Yaw_Combined_Mean <- lmer(Mean ~ species + (1|Ind), data = pb_atpec$Yaw_Combined)
performance::r2_xu(pb_atpec_lmer_Yaw_Combined_Mean) # Xu's R2 = 0.552

pb_atpec_lmer_Yaw_Combined_Tmax <- lmer(Tmax ~ species + (1|Ind), data = pb_atpec$Yaw_Combined)
performance::r2_xu(pb_atpec_lmer_Yaw_Combined_Tmax) # Xu's R2 = 0.424

pb_atpec_lmer_Yaw_Combined_Tmin <- lmer(Tmin ~ species + (1|Ind), data = pb_atpec$Yaw_Combined)
performance::r2_xu(pb_atpec_lmer_Yaw_Combined_Tmin) # Xu's R2 = 0.472



#### LMERS - manual - Propulsors - pb vs. at ####

## AbdAdd_Combined
pb_atpel_lmer_AbdAdd_Combined_Max <- lmer(Max ~ species + (1|Ind), data = pb_atpel$AbdAdd_Combined)
performance::r2_xu(pb_atpel_lmer_AbdAdd_Combined_Max) # Xu's R2 = 0.968

pb_atpel_lmer_AbdAdd_Combined_Min <- lmer(Min ~ species + (1|Ind), data = pb_atpel$AbdAdd_Combined)
performance::r2_xu(pb_atpel_lmer_AbdAdd_Combined_Min) # Xu's R2 = 0.955

pb_atpel_lmer_AbdAdd_Combined_Mean <- lmer(Mean ~ species + (1|Ind), data = pb_atpel$AbdAdd_Combined)
performance::r2_xu(pb_atpel_lmer_AbdAdd_Combined_Mean) # Xu's R2 = 0.978

pb_atpel_lmer_AbdAdd_Combined_Tmax <- lmer(Tmax ~ species + (1|Ind), data = pb_atpel$AbdAdd_Combined)
performance::r2_xu(pb_atpel_lmer_AbdAdd_Combined_Tmax) # Xu's R2 = 0.416

pb_atpel_lmer_AbdAdd_Combined_Tmin <- lmer(Tmin ~ species + (1|Ind), data = pb_atpel$AbdAdd_Combined)
performance::r2_xu(pb_atpel_lmer_AbdAdd_Combined_Tmin) # Xu's R2 = 0.387


## AnkAng_Combined
pb_atpel_lmer_AnkAng_Combined_Max <- lmer(Max ~ species + (1|Ind), data = pb_atpel$AnkAng_Combined)
performance::r2_xu(pb_atpel_lmer_AnkAng_Combined_Max) # Xu's R2 = 0.370

pb_atpel_lmer_AnkAng_Combined_Min <- lmer(Min ~ species + (1|Ind), data = pb_atpel$AnkAng_Combined)
performance::r2_xu(pb_atpel_lmer_AnkAng_Combined_Min) # Xu's R2 = 0.688

pb_atpel_lmer_AnkAng_Combined_Mean <- lmer(Mean ~ species + (1|Ind), data = pb_atpel$AnkAng_Combined)
performance::r2_xu(pb_atpel_lmer_AnkAng_Combined_Mean) # Xu's R2 = 0.899

pb_atpel_lmer_AnkAng_Combined_Tmax <- lmer(Tmax ~ species + (1|Ind), data = pb_atpel$AnkAng_Combined)
performance::r2_xu(pb_atpel_lmer_AnkAng_Combined_Tmax) # Xu's R2 = 0.430

pb_atpel_lmer_AnkAng_Combined_Tmin <- lmer(Tmin ~ species + (1|Ind), data = pb_atpel$AnkAng_Combined)
## singularity issue
performance::r2_xu(pb_atpel_lmer_AnkAng_Combined_Tmin) # Xu's R2 = 0.275


## KneeAng_Combined
pb_atpel_lmer_KneeAng_Combined_Max <- lmer(Max ~ species + (1|Ind), data = pb_atpel$KneeAng_Combined)
performance::r2_xu(pb_atpel_lmer_KneeAng_Combined_Max) # Xu's R2 = 0.368

pb_atpel_lmer_KneeAng_Combined_Min <- lmer(Min ~ species + (1|Ind), data = pb_atpel$KneeAng_Combined)
performance::r2_xu(pb_atpel_lmer_KneeAng_Combined_Min) # Xu's R2 = 0.809

pb_atpel_lmer_KneeAng_Combined_Mean <- lmer(Mean ~ species + (1|Ind), data = pb_atpel$KneeAng_Combined)
performance::r2_xu(pb_atpel_lmer_KneeAng_Combined_Mean) # Xu's R2 = 0.794

pb_atpel_lmer_KneeAng_Combined_Tmax <- lmer(Tmax ~ species + (1|Ind), data = pb_atpel$KneeAng_Combined)
performance::r2_xu(pb_atpel_lmer_KneeAng_Combined_Tmax) # Xu's R2 = 0.156

pb_atpel_lmer_KneeAng_Combined_Tmin <- lmer(Tmin ~ species + (1|Ind), data = pb_atpel$KneeAng_Combined)
performance::r2_xu(pb_atpel_lmer_KneeAng_Combined_Tmin) # Xu's R2 = 0.130


## Pitch_Combined
pb_atpel_lmer_Pitch_Combined_Max <- lmer(Max ~ species + (1|Ind), data = pb_atpel$Pitch_Combined)
performance::r2_xu(pb_atpel_lmer_Pitch_Combined_Max) # Xu's R2 = 0.967

pb_atpel_lmer_Pitch_Combined_Min <- lmer(Min ~ species + (1|Ind), data = pb_atpel$Pitch_Combined)
performance::r2_xu(pb_atpel_lmer_Pitch_Combined_Min) # Xu's R2 = 0.963

pb_atpel_lmer_Pitch_Combined_Mean <- lmer(Mean ~ species + (1|Ind), data = pb_atpel$Pitch_Combined)
performance::r2_xu(pb_atpel_lmer_Pitch_Combined_Mean) # Xu's R2 = 0.967

pb_atpel_lmer_Pitch_Combined_Tmax <- lmer(Tmax ~ species + (1|Ind), data = pb_atpel$Pitch_Combined)
performance::r2_xu(pb_atpel_lmer_Pitch_Combined_Tmax) # Xu's R2 = 0.447

pb_atpel_lmer_Pitch_Combined_Tmin <- lmer(Tmin ~ species + (1|Ind), data = pb_atpel$Pitch_Combined)
performance::r2_xu(pb_atpel_lmer_Pitch_Combined_Tmin) # Xu's R2 = 0.270

##I ran the Xu's R2 calculations and included the values that I got if they
# differ from what you had written 

## ProRet_Combined_fixed
pb_atpel_lmer_ProRet_Combined_fixed_Max <- lmer(Max ~ species + (1|Ind), data = pb_atpel$ProRet_Combined_fixed)
performance::r2_xu(pb_atpel_lmer_ProRet_Combined_fixed_Max) # Xu's R2 = 0.925

pb_atpel_lmer_ProRet_Combined_fixed_Min <- lmer(Min ~ species + (1|Ind), data = pb_atpel$ProRet_Combined_fixed)
performance::r2_xu(pb_atpel_lmer_ProRet_Combined_fixed_Min) # Xu's R2 = 0.713

pb_atpel_lmer_ProRet_Combined_fixed_Mean <- lmer(Mean ~ species + (1|Ind), data = pb_atpel$ProRet_Combined_fixed)
performance::r2_xu(pb_atpel_lmer_ProRet_Combined_fixed_Mean) # Xu's R2 = 0.845

pb_atpel_lmer_ProRet_Combined_fixed_Tmax <- lmer(Tmax ~ species + (1|Ind), data = pb_atpel$ProRet_Combined_fixed)
performance::r2_xu(pb_atpel_lmer_ProRet_Combined_fixed_Tmax) # Xu's R2 = 0.011

pb_atpel_lmer_ProRet_Combined_fixed_Tmin <- lmer(Tmin ~ species + (1|Ind), data = pb_atpel$ProRet_Combined_fixed)
# singularity issue
performance::r2_xu(pb_atpel_lmer_ProRet_Combined_fixed_Tmin) # Xu's R2 = 0.317


## Yaw_Combined
pb_atpel_lmer_Yaw_Combined_Max <- lmer(Max ~ species + (1|Ind), data = pb_atpel$Yaw_Combined)
performance::r2_xu(pb_atpel_lmer_Yaw_Combined_Max) # Xu's R2 = 0.719

pb_atpel_lmer_Yaw_Combined_Min <- lmer(Min ~ species + (1|Ind), data = pb_atpel$Yaw_Combined)
performance::r2_xu(pb_atpel_lmer_Yaw_Combined_Min) # Xu's R2 = 0.649

pb_atpel_lmer_Yaw_Combined_Mean <- lmer(Mean ~ species + (1|Ind), data = pb_atpel$Yaw_Combined)
performance::r2_xu(pb_atpel_lmer_Yaw_Combined_Mean) # Xu's R2 = 0.448

pb_atpel_lmer_Yaw_Combined_Tmax <- lmer(Tmax ~ species + (1|Ind), data = pb_atpel$Yaw_Combined)
performance::r2_xu(pb_atpel_lmer_Yaw_Combined_Tmax) # Xu's R2 = 0.218

pb_atpel_lmer_Yaw_Combined_Tmin <- lmer(Tmin ~ species + (1|Ind), data = pb_atpel$Yaw_Combined)
performance::r2_xu(pb_atpel_lmer_Yaw_Combined_Tmin) # Xu's R2 = 0.457


#### LMERS - manual - At - forelimb vs. hindlimb ####

## AbdAdd_Combined
at_pecpel_lmer_AbdAdd_Combined_Max <- lmer(Max ~ appendage + (1|Ind), data = at_pecpel$AbdAdd_Combined)
performance::r2_xu(at_pecpel_lmer_AbdAdd_Combined_Max) # Xu's R2 = 0.815

at_pecpel_lmer_AbdAdd_Combined_Min <- lmer(Min ~ appendage + (1|Ind), data = at_pecpel$AbdAdd_Combined)
performance::r2_xu(at_pecpel_lmer_AbdAdd_Combined_Min) # Xu's R2 = 0.621

at_pecpel_lmer_AbdAdd_Combined_Mean <- lmer(Mean ~ appendage + (1|Ind), data = at_pecpel$AbdAdd_Combined)
performance::r2_xu(at_pecpel_lmer_AbdAdd_Combined_Mean) # Xu's R2 = 0.851

at_pecpel_lmer_AbdAdd_Combined_Tmax <- lmer(Tmax ~ appendage + (1|Ind), data = at_pecpel$AbdAdd_Combined)
performance::r2_xu(at_pecpel_lmer_AbdAdd_Combined_Tmax) # Xu's R2 = 0.679

at_pecpel_lmer_AbdAdd_Combined_Tmin <- lmer(Tmin ~ appendage + (1|Ind), data = at_pecpel$AbdAdd_Combined)
performance::r2_xu(at_pecpel_lmer_AbdAdd_Combined_Tmin) # Xu's R2 = 0.587


## AnkAng_Combined
at_pecpel_lmer_AnkAng_Combined_Max <- lmer(Max ~ appendage + (1|Ind), data = at_pecpel$AnkAng_Combined)
performance::r2_xu(at_pecpel_lmer_AnkAng_Combined_Max) # Xu's R2 = 0.274

at_pecpel_lmer_AnkAng_Combined_Min <- lmer(Min ~ appendage + (1|Ind), data = at_pecpel$AnkAng_Combined)
performance::r2_xu(at_pecpel_lmer_AnkAng_Combined_Min) # Xu's R2 = 0.843

at_pecpel_lmer_AnkAng_Combined_Mean <- lmer(Mean ~ appendage + (1|Ind), data = at_pecpel$AnkAng_Combined)
performance::r2_xu(at_pecpel_lmer_AnkAng_Combined_Mean) # Xu's R2 = 0.859

at_pecpel_lmer_AnkAng_Combined_Tmax <- lmer(Tmax ~ appendage + (1|Ind), data = at_pecpel$AnkAng_Combined)
performance::r2_xu(at_pecpel_lmer_AnkAng_Combined_Tmax) # Xu's R2 = 0.068

at_pecpel_lmer_AnkAng_Combined_Tmin <- lmer(Tmin ~ appendage + (1|Ind), data = at_pecpel$AnkAng_Combined)
## singularity issue
performance::r2_xu(at_pecpel_lmer_AnkAng_Combined_Tmin) # Xu's R2 = 0.300


## KneeAng_Combined
at_pecpel_lmer_KneeAng_Combined_Max <- lmer(Max ~ appendage + (1|Ind), data = at_pecpel$KneeAng_Combined)
performance::r2_xu(at_pecpel_lmer_KneeAng_Combined_Max) # Xu's R2 = 0.509

at_pecpel_lmer_KneeAng_Combined_Min <- lmer(Min ~ appendage + (1|Ind), data = at_pecpel$KneeAng_Combined)
performance::r2_xu(at_pecpel_lmer_KneeAng_Combined_Min) # Xu's R2 = 0.425

at_pecpel_lmer_KneeAng_Combined_Mean <- lmer(Mean ~ appendage + (1|Ind), data = at_pecpel$KneeAng_Combined)
performance::r2_xu(at_pecpel_lmer_KneeAng_Combined_Mean) # Xu's R2 = 0.636

at_pecpel_lmer_KneeAng_Combined_Tmax <- lmer(Tmax ~ appendage + (1|Ind), data = at_pecpel$KneeAng_Combined)
performance::r2_xu(at_pecpel_lmer_KneeAng_Combined_Tmax) # Xu's R2 = 0.670

at_pecpel_lmer_KneeAng_Combined_Tmin <- lmer(Tmin ~ appendage + (1|Ind), data = at_pecpel$KneeAng_Combined)
performance::r2_xu(at_pecpel_lmer_KneeAng_Combined_Tmin) # Xu's R2 = 0.237


## Pitch_Combined
at_pecpel_lmer_Pitch_Combined_Max <- lmer(Max ~ appendage + (1|Ind), data = at_pecpel$Pitch_Combined)
performance::r2_xu(at_pecpel_lmer_Pitch_Combined_Max) # Xu's R2 = 0.867

at_pecpel_lmer_Pitch_Combined_Min <- lmer(Min ~ appendage + (1|Ind), data = at_pecpel$Pitch_Combined)
performance::r2_xu(at_pecpel_lmer_Pitch_Combined_Min) # Xu's R2 = 0.750

at_pecpel_lmer_Pitch_Combined_Mean <- lmer(Mean ~ appendage + (1|Ind), data = at_pecpel$Pitch_Combined)
performance::r2_xu(at_pecpel_lmer_Pitch_Combined_Mean) # Xu's R2 = 0.821

at_pecpel_lmer_Pitch_Combined_Tmax <- lmer(Tmax ~ appendage + (1|Ind), data = at_pecpel$Pitch_Combined)
performance::r2_xu(at_pecpel_lmer_Pitch_Combined_Tmax) # Xu's R2 = 0.116

at_pecpel_lmer_Pitch_Combined_Tmin <- lmer(Tmin ~ appendage + (1|Ind), data = at_pecpel$Pitch_Combined)
performance::r2_xu(at_pecpel_lmer_Pitch_Combined_Tmin) # Xu's R2 = 0.125


## ProRet_Combined_fixed
at_pecpel_lmer_ProRet_Combined_fixed_Max <- lmer(Max ~ appendage + (1|Ind), data = at_pecpel$ProRet_Combined_fixed)
performance::r2_xu(at_pecpel_lmer_ProRet_Combined_fixed_Max) # Xu's R2 = 0.811

at_pecpel_lmer_ProRet_Combined_fixed_Min <- lmer(Min ~ appendage + (1|Ind), data = at_pecpel$ProRet_Combined_fixed)
performance::r2_xu(at_pecpel_lmer_ProRet_Combined_fixed_Min) # Xu's R2 = 0.331

at_pecpel_lmer_ProRet_Combined_fixed_Mean <- lmer(Mean ~ appendage + (1|Ind), data = at_pecpel$ProRet_Combined_fixed)
performance::r2_xu(at_pecpel_lmer_ProRet_Combined_fixed_Mean) # Xu's R2 = 0.869

at_pecpel_lmer_ProRet_Combined_fixed_Tmax <- lmer(Tmax ~ appendage + (1|Ind), data = at_pecpel$ProRet_Combined_fixed)
performance::r2_xu(at_pecpel_lmer_ProRet_Combined_fixed_Tmax) # Xu's R2 = 0.399

at_pecpel_lmer_ProRet_Combined_fixed_Tmin <- lmer(Tmin ~ appendage + (1|Ind), data = at_pecpel$ProRet_Combined_fixed)
performance::r2_xu(at_pecpel_lmer_ProRet_Combined_fixed_Tmin) # Xu's R2 = 0.415


## Yaw_Combined
at_pecpel_lmer_Yaw_Combined_Max <- lmer(Max ~ appendage + (1|Ind), data = at_pecpel$Yaw_Combined)
performance::r2_xu(at_pecpel_lmer_Yaw_Combined_Max) # Xu's R2 = 0.193

at_pecpel_lmer_Yaw_Combined_Min <- lmer(Min ~ appendage + (1|Ind), data = at_pecpel$Yaw_Combined)
performance::r2_xu(at_pecpel_lmer_Yaw_Combined_Min) # Xu's R2 = 0.308

at_pecpel_lmer_Yaw_Combined_Mean <- lmer(Mean ~ appendage + (1|Ind), data = at_pecpel$Yaw_Combined)
performance::r2_xu(at_pecpel_lmer_Yaw_Combined_Mean) # Xu's R2 = 0.148

at_pecpel_lmer_Yaw_Combined_Tmax <- lmer(Tmax ~ appendage + (1|Ind), data = at_pecpel$Yaw_Combined)
performance::r2_xu(at_pecpel_lmer_Yaw_Combined_Tmax) # Xu's R2 = 0.673

at_pecpel_lmer_Yaw_Combined_Tmin <- lmer(Tmin ~ appendage + (1|Ind), data = at_pecpel$Yaw_Combined)
performance::r2_xu(at_pecpel_lmer_Yaw_Combined_Tmin) # Xu's R2 = 0.137


####Check for Singularity####
pb_atpec_Singular <- vector()
pb_atpel_Singular <- vector()
at_pecpel_Singular <- vector()

#Checks LMM results for singularity
for(i in 1:length(pb_atpec_lmer)) {
  
pb_atpec_Singular[i] <- isSingular(pb_atpec_lmer[[i]])
pb_atpel_Singular[i] <- isSingular(pb_atpel_lmer[[i]])
at_pecpel_Singular[i] <- isSingular(at_pecpel_lmer[[i]])
}

#Return index of singular results
if(any(pb_atpec_Singular)){
  pb_atpec_Singular <- which(pb_atpec_Singular)
  warning('Singular LMM Results in pb_atpec_lmer[...]')
  print(pb_atpec_Singular)
} else{pb_atpec_Singular <- NULL}
names(pb_atpec_lmer)[pb_atpec_Singular]

if(any(pb_atpel_Singular)){
  pb_atpel_Singular <- which(pb_atpel_Singular)
  warning('Singular LMM Results at pb_atpel_lmer[...]')
  print(pb_atpel_Singular)
} else{pb_atpel_Singular <- NULL}
names(pb_atpel_lmer)[pb_atpel_Singular]

if(any(at_pecpel_Singular)){
  at_pecpel_Singular <- which(at_pecpel_Singular)
  warning('Singular LMM Results at at_pecpel_lmer[...]')
  print(at_pecpel_Singular)
} else{at_pecpel_Singular <- NULL}
names(at_pecpel_lmer)[at_pecpel_Singular]
# Look at the lmer model that had the singularity
at_pecpel_lmer[at_pecpel_Singular]

# Manually running the lmer with the singularity issue to confirm it occurs
at_pel_AnkAng_Combined
at_pecpel_AnkAng_Combined <- rbind(at_pec_AnkAng_Combined, at_pel_AnkAng_Combined)
at_pecpel_AnkAng_Combined$Tmin <- apply(at_pecpel_AnkAng_Combined[, 1:101], 1, which.min)
at_pecpel_AnkAng_Combined$Ind <- substring(at_pecpel_AnkAng_Combined$filename, 1, 4)
lmer(Tmin ~ appendage + (1|Ind), data = at_pecpel_AnkAng_Combined)

# Checking what the mean and sd are for Tmin 
aggregate(at_pec$Tmin, list(at_pec$Ind),  function(x) c(mean = mean(x), sd = sd(x)))




#### UNUSED CODE ####


# ####### SAVING DATA FOR FIGURE 2 IN SHEFFIELD AND BLOB (2008) #########
# pb_kinSave <- rbind(	cbind(AfHL.FemHzAng.Combined, Variable = "AbductAdductAngle"),
#                    cbind(AfHL.AnkAng.Combined, Variable = "WristAnkleAngle"),
#                    cbind(AfHL.FemTVAng.Corr.Combined, Variable = "ProtractRetractAngle"),
#                    cbind(AfHL.AnkAng.Combined, Variable = "KneeElbowAngle"),
#                    cbind(AfFL.FemHzAng.Combined, Variable = "AbductAdductAngle"),
#                    cbind(AfFL.AnkAng.Combined, Variable = "WristAnkleAngle"),
#                    cbind(AfFL.FemTVAng.Corr.Combined, Variable = "ProtractRetractAngle"),
#                    cbind(pb_KneeAng_Combined, Variable = "KneeElbowAngle")
# 
# 
# 
# #setwd('/Users/SandyMKawano/Google Drive/Research/Fish and Salamander Bone Loading/')
# KineSaveName <- paste("Bone Load Kinematic Data_pb_", SaveDate, ".csv", sep="" )
# write.csv(KineSave, file = KineSaveName, row.names = FALSE)




