import pandas as pd
import numpy as np

import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import accuracy_score

from sklearn.preprocessing import StandardScaler, MinMaxScaler




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

scaler = StandardScaler()
expressionscaled = scaler.fit_transform(expressions)

# _df = []

# for celltype in celltypes:
#     _expressions = expressions
#     _expressions[[i for i in _expressions.columns if i not in markers[celltype]]] = 0
#     __df = pd.concat([_expressions,df[celltype]], axis=1)
#     __df.rename(columns={celltype:'cellcount'}, inplace=True)
#     _df.append(__df)

# _df = pd.concat(_df, ignore_index = True)

# _df.info()
# _df.describe()
# _df.groupby('cellcount').size()

# _df.hist()
# plt.show()


celltype = "Astrocyte"

_expressions = expressions[markers[celltype]]
_df = pd.concat([_expressions,df[celltype]], axis=1)

X = _df.drop(columns=[celltype])
y = _df[celltype]

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)


model = DecisionTreeClassifier()
model.fit(X_train, y_train)
yhat = model.predict(X_test)

accuracy_score(y_test, yhat)
