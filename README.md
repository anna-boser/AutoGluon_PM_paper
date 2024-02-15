# Accounting for spatio-temporal structure when validating environmental machine learning models
This repository holds all code necessary to conduct the analysis for *Accounting for spatio-temporal structure when validating environmental machine learning models*. 

Note: [*Full.csv*](https://drive.google.com/file/d/1eqrPDRcOFm9jOMIe5CWy8BxRTUBWGcf3/view?usp=sharing) and [*Small_Grid.csv*](https://drive.google.com/file/d/1M1mAIEVWL00jPFf2sNf22mitPDAUBPay/view?usp=sharing) are too large to be stored on GitHub, but can be downloaded at the provided links. 

Setup: 
In order to run the portion of this code that is in python (the .ipynb files), install conda and run the following commands in the terminal: 
```{bash}
conda create -n test1 python=3.8 numpy==1.19.2 pandas==1.4.2  -y 
conda activate test1
# Ensure the most up to date pip version
# pip install -U pip  
# pip install -U setuptools
# pip install -U "mxnet<2.0.0"
# Install pre-release, frozen to a particual pre-release for stability
pip install --pre "autogluon==0.0.16b20201214" 
# pip install -U ipykernel
# Optional -- makes environment discoverable in VS Code
# python -m ipykernel install --user --name autogluon_e2 --display-name "Python (autogluon_e2)" 
```

What happened: 
- successfully made original conda env with `conda create -n test1 python=3.8 numpy==1.19.2 pandas==1.4.2  -y `
- then ran `pip install --pre "autogluon==0.0.16b20201214"`, which also ran successfully. it however uninstalled numpy 1.19.2
- opened pyhon and ran `from autogluon.tabular import TabularPrediction as task`. Got error: 
```
AttributeError: module 'numpy' has no attribute 'float'.
`np.float` was a deprecated alias for the builtin `float`. To avoid this error in existing code, use `float` by itself. Doing this will not modify any behavior and is safe. If you specifically wanted the numpy scalar type, use `np.float64` here.
The aliases was originally deprecated in NumPy 1.20; for more details and guidance see the original release note at:
    https://numpy.org/devdocs/release/1.20.0-notes.html#deprecations
```

Try with updated autogluon version: 
```
conda create -n updated_autogluon python=3.8  -y 
conda activate updated_autogluon
pip install autogluon
pip install -U ipykernel
python -m ipykernel install --user --name updated_autogluon --display-name "Python (updated_autogluon)" 
```

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
