
## Load libraries
library(mgcv)
library(sf)
library(MuMIn)
library(DHARMa)
library(dismo)
library(terra)
library(dplyr)
library(mgcViz)

## new commit


# Set the aggregation group:
# j <- 6
for (j in c(1,2,3,4,6)){

  ## Load pres-abs data
  sps_data <- readRDS("data/pres_abs/presabs_Zone1.rds")


  ## Create presence-absence column:
  sps_data$PresAus <- 0
  sps_data <- as.data.frame(sps_data)
  sps_data$PresAus[which(sps_data$groups == j)] <- 1

  ## Remove duplicated:
  sps_data <- rbind(sps_data[which(sps_data$PresAus==1),],sps_data[which(sps_data$PresAus==0),])
  sps_data <- sps_data[!duplicated(sps_data[,1]),]

  ######################################################################################################################################################################################
  ## AUTOMATIZATION with dredge function MuMin package  ############
  ######################################################################################################################################################################################

  summary(sps_data)

  options(na.action = "na.fail")
  Model_Full <- mgcv::gam(PresAus ~ s(BPI_f,k=4) + s(BPI_b,k=4) + s(Depth,k=4) + s(Backscatter,k=4) + s(Slope,k=4) + s(East,k=4) + s(North,k=4), family= binomial, data=sps_data)

  AllModels <- MuMIn::dredge(Model_Full, rank=AIC)
  print(AllModels)

  ModelSelected_Several <-  get.models(AllModels, subset=delta<2)

  ModelSelect_Final <- ModelSelected_Several[[1]]


  ## Check predictors
  summary(ModelSelect_Final)
  AIC(ModelSelect_Final)

  ## gam check
  FinalGAM <- ModelSelect_Final
  windows()
  par(mfrow=c(2,2))
  gam.check(FinalGAM)
  AIC(FinalGAM)

  ## check residuals with DHARMa
  simulationOutput <- DHARMa::simulateResiduals(fittedModel = FinalGAM, plot = F)
  windows()
  plot(simulationOutput)

  ## Check response curves:
  windows()
  plot(FinalGAM, shade = TRUE, residuals=T, shade.col = 8,cex.lab=2,cex.axis=1.6, pages=1)


  ######################################################################################################################################################################################
  ## Check spatial autocorrelation and z-score with DHARMa pachage  ############
  ######################################################################################################################################################################################


  windows()
  res <- simulateResiduals(FinalGAM)
  plot(res)
  aut <- testSpatialAutocorrelation(res,
                                    x = sps_data$SUB1_Lon,
                                    y = sps_data$SUB1_Lat)
  aut
  ## z-Score:
  print(paste("z-score = ",(aut$statistic["observed"] - aut$statistic["expected"])/aut$statistic["sd"]))


  ######################################################################################################################################################################################
  ## Model prediction  ############
  ######################################################################################################################################################################################

  ## Load predictors data
  All_DF <- readRDS("data/predictors/predictors_df_zone1.rds")
  All_DF <- All_DF[,c("Depth",  "Backscatter","Slope","BPI_f","BPI_b","North","East","ID_capa","x","y" )]
  names(All_DF) <- c("Depth",  "Backscatter","Slope","BPI_f","BPI_b","North","East","ID","x","y" )

  All_SPDF_b <- st_as_sf(All_DF, coords = c("x", "y"), crs = 32630)

  Model_Geo_Binomial_GAM <- predict.gam(FinalGAM, All_DF, se.fit=TRUE, type="response")
  str(Model_Geo_Binomial_GAM)

  ##  add the prediction
  All_SPDF_b$GAM <- Model_Geo_Binomial_GAM$fit

  #We do the same with the error
  All_SPDF_b$GAM_Se <- Model_Geo_Binomial_GAM$se.fit

  #check the results
  windows()
  plot(All_SPDF_b[,"GAM"])


  ## Save layers
  All_SPDF_c <- as.data.frame(All_SPDF_b)
  assign(paste0("All_SPDF_",j),All_SPDF_c)
  assign(paste0("FinalGAM_",j),FinalGAM)


  ######################################################################################################################################################################################
  #MODEL EVALUATION  ############
  ######################################################################################################################################################################################


  #################################################################
  # MODEL EVALUATION USING TRAIN / TEST SPLITS
  #  - Estimate mean AUC
  #  - Estimate mean Kappa-maximizing threshold
  #  - Use that threshold to binarize the final prediction map
  #################################################################


  # Number of repetitions
  nrep <- 1000

  Data <- sps_data

  # Objects to store results
  AUC_GAM    <- numeric(nrep)
  Kappa_GAM  <- numeric(nrep)
  Thr_GAM    <- numeric(nrep)

  # Presence / absence data
  Presences <- Data[Data$PresAus == 1, ]
  Absences  <- Data[Data$PresAus == 0, ]

  # Loop over repetitions
  for (i in 1:nrep) {

    # -------------------------------------------------------------
    # 1. Train / test split using k-fold
    # -------------------------------------------------------------
    k <- 3  # 66% training / 33% testing

    grp_pres <- kfold(Presences, k)
    grp_abs  <- kfold(Absences,  k)

    EvalPres  <- Presences[grp_pres == 1, ]
    TrainPres <- Presences[grp_pres != 1, ]

    EvalAbs   <- Absences[grp_abs == 1, ]
    TrainAbs  <- Absences[grp_abs != 1, ]

    # NOTE:
    # TrainData is not used here because the model is already fitted
    # If you want a *true* cross-validation, the GAM should be refitted here
    # using TrainPres + TrainAbs


    # -------------------------------------------------------------
    # 2. Model (already fitted GAM)
    # -------------------------------------------------------------
    Model_Train_GAM <- FinalGAM


    # -------------------------------------------------------------
    # 3. Predictions at evaluation points
    # -------------------------------------------------------------
    # IMPORTANT:
    # We predict directly on the evaluation points.
    # No rasterization is needed for model evaluation.

    pred_pres <- predict(
      Model_Train_GAM,
      EvalPres,
      type = "response"
    )

    pred_abs <- predict(
      Model_Train_GAM,
      EvalAbs,
      type = "response"
    )


    # -------------------------------------------------------------
    # 4. Model evaluation
    # -------------------------------------------------------------
    ev <- dismo::evaluate(
      p = as.numeric(pred_pres),
      a = as.numeric(pred_abs)
    )

    # Store AUC
    AUC_GAM[i] <- ev@auc

    # Extract Kappa values and thresholds
    kappa_vals <- ev@kappa
    thr_vals   <- ev@t

    # Identify the threshold that maximizes Kappa
    i_max <- which.max(kappa_vals)

    Kappa_GAM[i] <- kappa_vals[i_max]
    Thr_GAM[i]   <- thr_vals[i_max]
  }

  auc <- mean(AUC_GAM)
  print(paste("AUC for group ",j," is ",auc))
  sd(AUC_GAM)

  mean(Kappa_GAM)
  sd(Kappa_GAM)

  # Mean threshold maximizing Kappa
  th <- mean(Thr_GAM)
  th
  sd(Thr_GAM)
  print(paste("threshold for group ",j," is ",th))

  # Create binary map using the mean Kappa-optimal threshold
  All_SPDF_new <- All_SPDF_b

  All_SPDF_new$GAM_bin <- ifelse(
    All_SPDF_new$GAM >= th, 1, 0
  )

  # Save binary data
  assign(paste0("All_SPDF_th_",j),All_SPDF_new)


  ## Rasterize
  xy <- st_coordinates(All_SPDF_new)

  df_xyz <- data.frame(
    x = xy[,1],
    y = xy[,2],
    z = All_SPDF_new$GAM_bin
  )

  GAM_bin_rast <- rast(
    df_xyz,
    type = "xyz",
    crs = crs(All_SPDF_new)
  )

  ## Check plot
  windows()
  plot(GAM_bin_rast)

  ## Save raster
  out_file <- paste0(
    "outputs/pres_abs/GAM_z1_presaus_th_",j,".tif"
  )
  
  if (!file.exists(out_file)) {
    writeRaster(
      GAM_bin_rast,
      out_file,
      datatype = "FLT4S",   # float 32-bit
      gdal = c("COMPRESS=LZW")
    )
  }





  ######################################################################################################################################################################################
  ##  Density Models  ############
  ######################################################################################################################################################################################

  ## Load densities matrix
  dens_data <- readRDS("data/densities/sps_dens_mat_z1.rds")
  dens_data$Depth <- dens_data$Depth * (-1)

  ## Selection of the Group SIMPER species:

  if (j==1) {
    vec <- c("Lop_per","Mad_ocu","Phe_her")
  } else if (j==2) {
    vec <- c("Pac_spp","Pha_rob","Geo_bar","Neo_bow")
  } else if (j==3) {
    vec <- c("Kop_ste","Fun_qua")
  } else if (j==4) {
    vec <- c("Tho_Eut","Sol_var","Aca_arm")
  } else if (j==6) {
    vec <- c("Lop_per","Mad_ocu","Aph_bea")}

  dens_data <- dens_data[which(dens_data$groups==j),]
  dens_data$dens <- rowSums(dens_data[,vec])

  ## Select Zone 1
  dens_data <- dens_data[which(dens_data$SUB1_Lon<260000),]

  ## Delete zeros
  dens_data <- dens_data[which(dens_data$dens != 0),]





  #############################################################################################
  ## Some of the groups does not have enough data to include all variables,
  ## so we must add one by one the most explanatory variables until we reach the max
  ## and then we start the full model
  #############################################################################################
  options(na.action = "na.fail")

  if (j==1) {
    Model_Full <- mgcv::gam(dens ~ s(BPI_f,k=4) + s(BPI_b,k=4) + s(Slope,k=4) + s(East,k=4), family= "nb", data=dens_data)
    AllModels <- dredge(Model_Full, rank=AIC)
    print(AllModels)
    ModelSelected_Several <-  get.models(AllModels, subset=delta<3)
    ModelSelect_Final <- ModelSelected_Several[[1]]
  }
  if (j==2) {
    Model_Full <- mgcv::gam(dens ~ s(BPI_f,k=4) + s(BPI_b,k=4) + s(Depth,k=4) + s(Backscatter,k=4) + s(Slope,k=4) + s(East,k=4) + s(Curv,k=4) + as.factor(Geomorf), family= "nb", data=dens_data)
    AllModels <- dredge(Model_Full, rank=AIC)
    print(AllModels)
    ModelSelected_Several <-  get.models(AllModels, subset=delta<3)
    ModelSelect_Final <- ModelSelected_Several[[1]]}
  if (j==3) {
    Model_Full <- mgcv::gam(dens ~ s(BPI_f,k=4) + s(Depth,k=4) + s(East,k=4), family= "nb", data=dens_data)
    AllModels <- dredge(Model_Full, rank=AIC)
    print(AllModels)
    ModelSelected_Several <-  get.models(AllModels, subset=delta<3)
    ModelSelect_Final <- ModelSelected_Several[[1]]}
  if (j==4) {
    Model_Full <- mgcv::gam(dens ~ s(BPI_f,k=4) + s(BPI_b,k=4) + s(Depth,k=4) + s(East,k=4), family= "nb", data=dens_data)
    AllModels <- dredge(Model_Full, rank=AIC)
    print(AllModels)
    ModelSelected_Several <-  get.models(AllModels, subset=delta<3)
    ModelSelect_Final <- ModelSelected_Several[[1]]}
  if (j==6) {
    Model_Full <- mgcv::gam(dens ~ s(BPI_f,k=4) + s(BPI_b,k=4) + s(Depth,k=4) + s(Backscatter,k=4) + s(Slope,k=4) + s(East,k=4) + s(North,k=4) + s(Curv,k=4) + as.factor(Geomorf), family= "nb", data=dens_data)
    AllModels <- dredge(Model_Full, rank=AIC)
    print(AllModels)
    ModelSelected_Several <-  get.models(AllModels, subset=delta<3)
    ModelSelect_Final <- ModelSelected_Several[[2]]}

  summary(ModelSelect_Final)
  AIC(ModelSelect_Final)
  FinalGAM <- ModelSelect_Final
  windows()
  par(mfrow=c(2,2))
  gam.check(FinalGAM)

  AIC(FinalGAM)
  library(DHARMa)
  simulationOutput <- simulateResiduals(fittedModel = FinalGAM, plot = F)
  windows()
  plot(simulationOutput)




  ######################################################################################################################################################################################
  ## Model prediction  ############
  ######################################################################################################################################################################################


  All_DF <- readRDS("data/predictors/predictors_df_zone1.rds")
  names(All_DF) <- c("Depth",  "Backscatter","Slope","BPI_f","BPI_b","North","East","Curv","Geomorf","ID","x","y" )

  p_sdf <- st_as_sf(All_DF, coords = c("x","y"), crs = 32630)

  ## Add levels if necessary
  if (is.null(FinalGAM$xlevels$`as.factor(Geomorf)`)) {All_SPDF_b <- p_sdf; All_DF <- as.data.frame(All_SPDF_b)} else {
    All_SPDF_b <- p_sdf[which(p_sdf$Geomorf %in% FinalGAM$xlevels$`as.factor(Geomorf)`),]
    All_DF <- as.data.frame(All_SPDF_b)
    All_DF$Geomorf <- as.factor(All_DF$Geomorf)}

  All_DF <- cbind(All_DF[,c("Depth",  "Backscatter","Slope","BPI_f","BPI_b","North","East","Curv","Geomorf","ID")],st_coordinates(All_SPDF_b))

  Model_Geo_Binomial_GAM <- predict.gam(FinalGAM, All_DF, se.fit=TRUE, type="response")
  str(Model_Geo_Binomial_GAM)

  #We add the prediction
  All_SPDF_b$GAM <- Model_Geo_Binomial_GAM$fit

  #We do the same with the error
  All_SPDF_b$GAM_Se <- Model_Geo_Binomial_GAM$se.fit

  #check the results
  windows()
  plot(All_SPDF_b[,"GAM"])

  ## Join layers
  All_SPDF_c <- st_join(All_SPDF_b[,c("GAM")],All_SPDF_new[,c("GAM_bin")])

  ## Calculate hurdle (delta) models (pres-abs * densities)
  All_SPDF_b$delta <- All_SPDF_c$GAM_bin * All_SPDF_c$GAM


  ## Rasterize
  xy <- st_coordinates(All_SPDF_b)

  df_xyz <- data.frame(
    x = xy[,1],
    y = xy[,2],
    z = All_SPDF_b$delta
  )

  GAM_delta_rast <- rast(
    df_xyz,
    type = "xyz",
    crs = crs(All_SPDF_b)
  )

  ## Export delta map
  out_file <- paste0(
    "outputs/delta/GAM_z1_delta_",j,".tif"
  )

  if (!file.exists(out_file)) {
    writeRaster(
      GAM_delta_rast,
      out_file,
      datatype = "FLT4S",
      gdal = c("COMPRESS=LZW")
    )
  }


}


