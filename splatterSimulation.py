import numpy as np
import pandas as pd


df = pd.read_csv("data/splatter/splatterCSV.csv")
df.info()

celltypes = df['cellType'].unique()

df_celltype = {}
for celltype in celltypes:
    df_celltype[celltype] = df[df["cellType"] == celltype]
