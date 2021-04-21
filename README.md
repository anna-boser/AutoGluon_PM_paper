# AutoGluon_PM_paper
This repository holds all code necessary to conduct the analysis found in SOME PAPER TITLE. 

Please be sure to add the necessary datasets (SI FOUND AT THIS LOCATION) to the data folder before running the code. 

Contents: 

**data**
* *Small_Grid.csv*: Only the small grid values -- used to map and validate
* *Small_Ids.csv*: A list of the ids for the epa monitors found within the target prediction area (as opposed to the entire training area)
* *Train.csv*: The dataset used to train the models: only includes the pixels with monitors over the entire training area. 
* *Full.csv*: Dataset with all data on all 327 included days in 2017 for both training and target prediction areas. 
* *Modeling_Grid*: Shapefile for the entire training area
* *small_grid*: Shapefile only for target prediction area. 
* *model_outputs*: a folder containing the model predictions over the *Small_Grid* for all the trained models
* *CV_outputs*: a folder containing tables of the location-specific leave-one-out crossvalidations performed for each of the models. Each file contains the model performances (with uncertainties) as measured by bias, mean squared error, and $R^2$ for each monitoring location, as well as the aggregate for all locations in the larger training area and the aggregate for only the small target prediction area. 

**code**
* *Two_stage.Rmd*: two stage model
* *autogluon_try.ipynb*: autogluon script
