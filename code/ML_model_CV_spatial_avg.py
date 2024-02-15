from autogluon.tabular import TabularPredictor, TabularDataset
import pandas as pd
import numpy as np
from sklearn.model_selection import GroupKFold
import gc
import os
import glob

# df = pd.read_csv("Data/datasets/Train.csv")

# def cv(df):
#     n_fold = len(set(df['Id']))
#     kf = GroupKFold(n_fold)
#     split = kf.split(df, groups=df['Id'])

#     metric_df = pd.DataFrame(columns=['station', 'Day', 'PM', 'PM_pred', 'rmse', 'bias'])

#     for i, (train_idx, test_idx) in enumerate(split):
#         print(f'Starting training fold {i}.')
#         _ = gc.collect()

#         # Calculate the daily mean PM for each day, excluding the current validation set
#         validation_station = df.loc[test_idx, 'Id'].iloc[0]
#         train_df = df[~df['Id'].isin([validation_station])]
#         daily_mean_PM_excl_validation = train_df.groupby('Day')['PM'].mean()

#         # Add the daily mean PM (excluding validation station) as a feature to the original df
#         df['Daily_Mean_PM_Excl_Val'] = df['Day'].map(daily_mean_PM_excl_validation)

#         features = ['Day', 'Lat', 'Lon', 'Elevation', 'Emissions', 'Forest',
#                     'Roads', 'Streets', 'Plumes_High', 'Plumes_Med', 'Plumes_Low',
#                     'Max_Temp', 'Max_Wind', 'Precip', 'Rel_Humidity', 'Wind_Dir', 'BLH',
#                     'AOD', 'Daily_Mean_PM_Excl_Val']  # Include the new feature

#         # Adjust y by subtracting the mean PM for each day (excluding validation station)
#         df['Adjusted_PM'] = df['PM'] - df['Daily_Mean_PM_Excl_Val']
#         label_column = 'Adjusted_PM'

#         X = df[features]
#         y = df[label_column]

#         # Ensure the index alignment between X and y when slicing
#         train_data = TabularDataset(pd.concat([X.loc[train_idx], y.loc[train_idx]], axis=1))
#         test_data = TabularDataset(pd.concat([X.loc[test_idx], y.loc[test_idx]], axis=1))

#         predictor = TabularPredictor(label=label_column, eval_metric='r2').fit(train_data=train_data)
        
#         # task.fit(train_data=train_data, 
#         #                     label=label_column, 
#         #                     hyperparameters=hyperparameters, 
#         #                     # time_limits = 60*60, # un-comment these for the tuned neural network. Train for 1 hour. 
#         #                     # hyperparameter_tune=True, 
#         #                     eval_metric='r2')

#         y_test = test_data[label_column]
#         test_data_nolab = test_data.drop(labels=[label_column], axis=1)
#         y_pred = predictor.predict(test_data_nolab)

#         # Calculate metrics
#         bias = y_pred - y_test
#         rmse = np.abs(bias)

#         station = np.repeat(df.loc[test_idx]['Id'].iloc[0], len(test_data['Day']))
#         df_to_append = pd.DataFrame({'station': station, 'Day': test_data['Day'], 'PM': y_test + df.loc[test_idx, 'Daily_Mean_PM_Excl_Val'], 'PM_pred': y_pred + df.loc[test_idx, 'Daily_Mean_PM_Excl_Val'], 'rmse': rmse, 'bias': bias})

#         metric_df = pd.concat([metric_df, df_to_append], ignore_index=True)

#     return metric_df

# AutoGluon = cv(df)
# AutoGluon.to_csv("Data/output/CV/AutoGluon_spatial_var.csv", index = False)

small_ids = pd.read_csv("Data/datasets/Small_Ids.csv")['x']
small_ids = [float(i) for i in small_ids]

full_df = pd.DataFrame(columns = ['Model', 'test_area', 'station', 'Day', 'PM', 'PM_pred'])

paths = ["Data/output/CV/AutoGluon_spatial_var.csv", "Data/output/CV/AutoGluon.csv"]

for path in paths:
  f = os.path.split(path)[1]
  Model = os.path.splitext(f)[0]
  Table = pd.read_csv(path)
  Table['test_area'] = Table['station'].isin(small_ids)
  Table['Model'] = Model

  full_df = pd.concat([full_df, Table], ignore_index = True)

full_df.to_csv("Data/output/Autogluon_cv_preds_w_spatial_var.csv", index = False)