# REGRESSION MODELLING ON MERFISH DATASET
#
# NOTE: use spearman

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn import metrics
from scipy import stats

plt.rcParams['figure.figsize'] = [8,5]
plt.rcParams['font.size'] =14
plt.rcParams['font.weight']= 'bold'
plt.style.use('seaborn-whitegrid')


df = pd.read_csv("data/merfish/merfishSpatial.csv", sep=" ")
df.drop(columns=['coord'], inplace=True)

df.info()
df.head()

reference = pd.read_csv("data/merfish/markerGene_for_merfish_data.csv")
reference.drop(columns=['p_value'], inplace=True)
reference.drop(range(148,168))
reference.drop(reference.index[reference.cell_type=='EpendymalInhibitory'].tolist(), inplace=True)
reference['cell_type'] = reference['cell_type'].astype(str).str.replace(" ", "")

markers = reference.groupby('cell_type').agg(list).marker_gene
celltypes = reference.cell_type.unique().tolist()
expressions = df.drop(columns=celltypes)

xcel = None #np.array([[1,1,1,1,1]])

for celltype in celltypes:

    row = [celltype, len(markers[celltype])]

    _df = pd.concat([expressions[markers[celltype]], df[celltype]], axis=1)
    _df.rename(columns={celltype:'cellcount'}, inplace=True)

    X = _df.drop(columns=['cellcount'])
    y = _df['cellcount']
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=0)

    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_test = scaler.transform(X_test)

    model = LinearRegression()
    model.fit(X_train, y_train)
    yhat = model.predict(X_test)
    y_pred = np.round(yhat)
    y_pred = yhat

    # np.vstack((y_pred,y_test)).T

    r2_lr = round(metrics.r2_score(y_test,y_pred), 2)
    sp_lr = round(stats.spearmanr(y_test,y_pred).correlation, 2)

    model = RandomForestRegressor()
    model.fit(X_train, y_train)
    yhat = model.predict(X_test)
    y_pred = np.round(yhat)
    y_pred = yhat

    # np.vstack((y_pred,y_test)).T

    r2_rf = round(metrics.r2_score(y_test,y_pred), 2)
    sp_rf = round(stats.spearmanr(y_test,y_pred).correlation, 2)

    # print("--------", celltype, "--------")
    # print('R^2 Score(LR):', r2_lr)
    # print('Correlation(LR):', sp_lr)
    # print('R^2 Score(rf):', r2_rf)
    # print('Correlation(rf):', sp_rf)

    row = np.append(row, [r2_lr, sp_lr, r2_rf, sp_rf])
    # print(row, type(row))
    if xcel is None:
        xcel = np.array([row])
    else:
        xcel = np.concatenate((xcel, [row]))

xcel = pd.DataFrame(xcel, columns=['celltype', 'n_markers', 'R2(LR)', 'Spearman(LR)', 'R2(RF)', 'Spearman(RF)'])
print(xcel)
xcel.to_csv("output/merfish/merfish_model.csv", sep=',')
