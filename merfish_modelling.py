import pandas as pd
import numpy as np
import sklearn as skl

df = pd.read_csv("data/merfish/merfishSpatial.csv", sep=" ")
df.drop(columns=['coord'], inplace=True)

df.info()
df.head()

reference = pd.read_csv("data/merfish/markerGene_for_merfish_data.csv")
reference.drop(columns=['p_value'], inplace=True)
reference['cell_type'] = reference['cell_type'].astype(str).str.replace(" ", "")
markers = reference.groupby('cell_type').agg(list).marker_gene
celltypes = list(reference.cell_type.unique())

_df = df[markers['Astrocyte'] + ['Astrocyte']]
