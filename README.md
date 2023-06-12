# Accounting for spatio-temporal structure when validating environmental machine learning models
This repository holds all code necessary to conduct the analysis for *Accounting for spatio-temporal structure when validating environmental machine learning models*. 

Note: [*Full.csv*](https://drive.google.com/file/d/1eqrPDRcOFm9jOMIe5CWy8BxRTUBWGcf3/view?usp=sharing) and [*Small_Grid.csv*](https://drive.google.com/file/d/1M1mAIEVWL00jPFf2sNf22mitPDAUBPay/view?usp=sharing) are too large to be stored on GitHub, but can be downloaded at the provided links. 

Contents: 

**Data**
* *california*: A shapefile of California. 
* *Full.csv*: Dataset with all data on all 327 included days in 2017 for both training and target prediction areas. 
* *Small_Grid.csv*: Only the small grid values -- used to map and validate
* *Small_Ids.csv*: A list of the ids for the epa monitors found within the target prediction area (as opposed to the entire training area)
* *Train.csv*: The dataset used to train the models: only includes the pixels with monitors over the entire training area. 
* *Modeling_Grid*: Shapefile for the entire training area
* *small_grid*: Shapefile only for target prediction area. 
* *model_outputs*: a folder containing the model predictions over the *Small_Grid* for all the trained models
* *CV_outputs*: a folder containing tables of the location-specific leave-one-out crossvalidations performed for each of the models. Each file contains the model performances (with uncertainties) as measured by bias, mean squared error, and R^2 for each monitoring location, as well as the aggregate for all locations in the larger training area and the aggregate for only the small target prediction area. 

**code**
* *2S_CV.Rmd*: two stage model training and evaluation (crossvalidation). *2S_CV_train_bw.Rmd* can be used to train the bandwidth parameter. 
* *2S.Rmd*: two stage model application to prediction area to create map. 
* *ML_model_CV.ipynb*: autogluon model training and evaluation (crossvalidation)
* *ML_models_grid.ipynb*: autogluon model application to prediction area to create map
* *CV_figures.Rmd*: figures of model cross-validation results
* *Plot_grids.Rmd*: maps of models
