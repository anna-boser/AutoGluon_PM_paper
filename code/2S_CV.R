library(cvTools)
library(dplyr)
library(tidyr)
library(lme4)
library(spgwr)
library(leaps)
library(sp)

pod = TRUE # are you running this on pod or locally? 

if (pod){
  setwd("~/AutoGluon_PM_paper")
} else {
  setwd("~/Documents/GitHub/AutoGluon_PM_paper")
}

train_bw <- FALSE #train the bandwidth or use pre-trained bandwidth? 
spatial_cv <- FALSE #spatial CV or random 10-fold? 

################################################################################
########################### Load data ##########################################
################################################################################

Data <- read.csv(file = "Data/datasets/Train.csv")

Data$Day <- as.factor(Data$Day)

#split data between with and without AOD
w_AOD <- filter(Data, !is.na(AOD))

indfunc <- function(val){
  if (is.na(val)){
    return(0)
  } else {
    return(1)
  }
}
ind_AOD <- Data
ind_AOD$AOD <- unlist(lapply(ind_AOD$AOD, indfunc))  

################################################################################
########################### Function to get LME formula ########################
################################################################################

make_formula <- function(variables){
  time_varrying_variables <- c("Plumes_High", "Plumes_Med", "Plumes_Low", "Max_Temp", "Max_Wind", "Precip", "Rel_Humidity", "Wind_Dir", "BLH", "AOD")
  formula <- "PM ~"
  for (var in variables){
    formula <- paste(formula, var, "+")
  }
  variables <- variables[variables %in% time_varrying_variables]
  formula <- paste(formula, "(1")
  for (var in variables){
    formula <- paste(formula, "+", var)
  }
  formula <- paste(formula,  "| Day)")
}

get_r2 <- function(formula, Data){
  # print(paste0("training a model with formula: ", formula))
  model <- lmer(formula,
                data = Data)
  r2 <- 1 - (sum((predict(model) - Data$PM)^2)/sum((Data$PM - mean(Data$PM))^2))
  return(r2)
}

new_var <- function(previous_r2, current_variables, remaining_variables, Data){ #function that sees if adding another variable increases the r2, and if so chooses the best variable
  
  r2s <- c()
  
  for (new_var in remaining_variables){
    variables <- c(current_variables, new_var)
    formula <- make_formula(variables)
    r2s <- c(r2s, get_r2(formula, Data))
  }
  
  best_r2 <- max(r2s)
  best_variable <- remaining_variables[r2s == best_r2]
  
  if(best_r2 < previous_r2+0.005){#the r2 needs to increase by more than 005
    # print("should be done")
    return(NULL)
  } else {
    return(list("var"=best_variable, "r2"=best_r2))
  }
}

variables <- c("Elevation","Emissions","Forest", "Roads","Streets", "Plumes_High" ,"Plumes_Med" , "Plumes_Low", "Max_Temp"  , "Max_Wind" ,"Precip" , "Rel_Humidity","Wind_Dir" ,"BLH", "AOD")     

formula_search <- function(Data){
  current_variables = c()
  remaining_variables = variables
  previous_r2 = 0
  
  for (i in 1:length(variables)){
    newvar <- new_var(previous_r2, current_variables, remaining_variables, Data)
    if (is.null(newvar)){
      print("Done. Variables are:")
      print(current_variables)
      break
    } else {
      current_variables = c(current_variables, newvar$var)
      remaining_variables = remaining_variables[remaining_variables != newvar$var]
      previous_r2 = newvar$r2
      print(paste0("new r2 is ", previous_r2, " and the new variable is ", newvar$var))
    }
  }
  
  LME_formula <- make_formula(current_variables)
  
  return(LME_formula)
}

################################################################################
########################### Define model for cross-validation ##################
################################################################################

Two_stage_predictions <- function(train_AOD, train_wo_AOD, test_AOD, test_wo_AOD, bw = NA, bw2 = NA){
  # The way the two-stage model is broken down: 
  # 1: use the training data to fit a linear mixed effects model
  # 2: using the training data, predict the PM using the LME and calculate the residuals
  # 3: use the training data residuals to fit the GWR model on the residuals (actual - predicted) from the LME model
  # 4: use the test data to predict the PM with the LME model and the residuals from the LME model using the GWR model
  # 5: add the gwr-predicted residuals to the lme-predicted PM values to obtain the final prediction. 
  
  print("finding LME formula")
  lme_formula <- formula_search(train_AOD)
  gwr_formula <- "resid ~ AOD"
  
  test_AOD <- filter(test_AOD, Day %in% unique(train_AOD$Day)) #october 1, 13, and 31
  test_wo_AOD <- filter(test_wo_AOD, Day %in% unique(train_wo_AOD$Day))
  
  # 1:
  AOD_lme <- lmer(formula = lme_formula, 
                  data = train_AOD)
  wo_AOD_lme <- lmer(formula = lme_formula, 
                     data = train_wo_AOD)
  
  print("LME trained")
  
  # 2:
  train_AOD$resid <- train_AOD$PM - predict(AOD_lme) #training residuals
  AOD_lme_predictions <- predict(AOD_lme, newdata = test_AOD) #testing predictions
  
  train_wo_AOD$resid <- train_wo_AOD$PM - predict(wo_AOD_lme) #training residuals
  wo_AOD_lme_predictions <- predict(wo_AOD_lme, newdata = test_wo_AOD) #testing predictions
  
  print("LME predictions and residuals calculated")
  
  # 3 and 4: 
  
  #convert train and test data to sp for easier use in spgwr
  train_sp_AOD <- SpatialPointsDataFrame(coords = cbind(train_AOD$Lon, train_AOD$Lat), data = train_AOD, proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  test_sp_AOD <- SpatialPointsDataFrame(coords = cbind(test_AOD$Lon, test_AOD$Lat), data = test_AOD, proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  
  train_sp_wo_AOD <- SpatialPointsDataFrame(coords = cbind(train_wo_AOD$Lon, train_wo_AOD$Lat), data = train_wo_AOD, proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  test_sp_wo_AOD <- SpatialPointsDataFrame(coords = cbind(test_wo_AOD$Lon, test_wo_AOD$Lat), data = test_wo_AOD, proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  
  print("Data converted to sp objects")
  
  
  if (is.na(bw)){
    gwr.bw.AOD <- gwr.sel(gwr_formula, data = train_sp_AOD, longlat = TRUE)
    gwr_model_AOD <- gwr(gwr_formula, data = train_sp_AOD, bandwidth = gwr.bw.AOD, longlat = TRUE, fit.points = test_sp_AOD, predictions = TRUE)
  } else {
    gwr_model_AOD <- gwr(gwr_formula, data = train_sp_AOD, bandwidth = bw, longlat = TRUE, fit.points = test_sp_AOD, predictions = TRUE)
  }
  
  if (is.na(bw2)){
    gwr.bw.wo.AOD <- gwr.sel(gwr_formula, data = train_sp_wo_AOD, longlat = TRUE)
    gwr_model_wo_AOD <- gwr(gwr_formula, data = train_sp_wo_AOD, bandwidth = gwr.bw.wo.AOD, longlat = TRUE, fit.points = test_sp_wo_AOD, predictions = TRUE)
  } else {
    gwr_model_wo_AOD <- gwr(gwr_formula, data = train_sp_wo_AOD, bandwidth = bw, longlat = TRUE, fit.points = test_sp_wo_AOD, predictions = TRUE)
  }
  
  print("GWR model finished")
  
  gwr_resid_AOD <- gwr_model_AOD$SDF$pred
  gwr_resid_wo_AOD <- gwr_model_wo_AOD$SDF$pred
  
  AOD_test_predictions <- gwr_resid_AOD + AOD_lme_predictions
  wo_AOD_test_predictions <- gwr_resid_wo_AOD + wo_AOD_lme_predictions
  
  test_AOD$PM_pred <- AOD_test_predictions
  test_wo_AOD$PM_pred <- wo_AOD_test_predictions
  
  print("predictions calculated")
  
  full_predictions <- rbind(test_AOD, test_wo_AOD)
  
  full_predictions$station <- full_predictions$Id
  full_predictions$rmse <- abs(full_predictions$PM - full_predictions$PM_pred)
  full_predictions$bias <- full_predictions$PM_pred - full_predictions$PM
  
  full_predictions <- full_predictions %>% select('station', 'Day', 'PM', 'PM_pred', 'rmse', 'bias')
  
  return(full_predictions) #return grid of both final values and the modeled residuals from the second stage
}

################################################################################
########################### Define cross-validation function ###################
################################################################################

CV <- function(Data, indData, bw = NA, bw2 = NA, spatial_cv = TRUE){
  
  #initialize the dataframe
  fold_df <- data.frame(matrix(ncol = 6, nrow = 0))
  colnames(fold_df) <- c('station', 'Day', 'PM', 'PM_pred', 'rmse', 'bias')
  
  if (spatial_cv){ # if spatial cv, use the station as the fold
    Data$fold <- Data$Id
    indData$fold <- indData$Id
  } else { # if we are doing a random 10-fold CV, the fold should be a random number between 1 and 10
    Data$fold <- sample(1:10, nrow(Data), replace = TRUE)
    indData$fold <- sample(1:10, nrow(indData), replace = TRUE)
  }
  
  folds <- unique(Data$fold)
  
  #loop through the folds and fit/test the full model each time.
  for(i in 1:length(folds)){
    print(paste0("Fold ", i))
    # The way the two-stage model is broken down: 
    # 1: use the training data to fit a linear mixed effects model
    # 2: using the training data, predict the PM using the LME and calculate the residuals
    # 3: use the training data residuals to fit the GWR model on the residuals (actual - predicted) from the LME model
    # 4: use the test data to predict the PM with the LME model and the residuals from the LME model using the GWR model
    # 5: add the gwr-predicted residuals to the lme-predicted PM values to obtain the final prediction. 
    
    train_D <- Data[Data$fold != folds[i],]
    test_D <- Data[Data$fold == folds[i] & Data$Day %in% unique(train_D$Day),]
    # test_D <- Data[station_folds == i & Data$Day,] #some days are not trained on. 
    train_ind <- indData[indData$fold != folds[i],]
    test_ind <- indData[indData$fold == folds[i] 
                        & indData$Day %in% unique(train_ind$Day) 
                        & indData$AOD == 0,]
    fold_df <- rbind(fold_df, Two_stage_predictions(train_AOD = train_D, 
                                                    train_wo_AOD = train_ind, 
                                                    test_AOD = test_D, 
                                                    test_wo_AOD = test_ind, 
                                                    bw = bw, 
                                                    bw2 = bw2))
    print("tail(fold_df):")
    print(tail(fold_df))
    
  }
  
  return(fold_df)
}

################################################################################
########################### Run the cross-validation ###########################
################################################################################

start_time <- Sys.time()

if (train_bw){
  leave_one_out <- CV(w_AOD, ind_AOD, bw = NA, bw2 = NA, spatial_cv)
} else {
  leave_one_out <- CV(w_AOD, ind_AOD, bw = 15, bw2 = 13, spatial_cv)
}

if (spatial_cv){
  write.csv(leave_one_out, file = "Data/output/CV/2S.csv", row.names = FALSE)
} else {
  write.csv(leave_one_out, file = "Data/output/CV/2S_random_CV.csv", row.names = FALSE)
}

end_time <- Sys.time()

end_time - start_time



