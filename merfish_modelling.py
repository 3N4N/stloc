import pandas as pd
import numpy as np
# import sklearn as skl

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

_X = []

for celltype in celltypes:
    _expressions = expressions
    _expressions[[i for i in _expressions.columns if i not in markers[celltype]]] = 0
    _df = pd.concat([_expressions,df[celltype]], axis=1)
    _df.rename(columns={celltype:'cellcount'}, inplace=True)
    _X.append(_df)

X = pd.concat(_X, ignore_index = True)
X.info()
