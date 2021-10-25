# Import pandas package

# Define a dictionary containing employee data
# data = {'Name':['Jai', 'Princi', 'Gaurav', 'Anuj'],
#         'Age':[27, 24, 22, 32],
#         'Address':['Delhi', 'Kanpur', 'Allahabad', 'Kannauj'],
#         'Qualification':['Msc', 'MA', 'MCA', 'Phd']}
 
# # Convert the dictionary into DataFrame 
# df = pd.DataFrame(data)
 
# # select two columns
# print(df[['Name', 'Qualification']])

import pandas as pd
from scipy.stats import mannwhitneyu

data = pd.read_csv("./data/moffitt/s7.csv")

# datatop = data.head()
# print(datatop)

data = data.drop(['Centroid_Y', 'Animal_ID', 'Animal_sex', 'Behavior', 'Bregma', 'Centroid_X', 'Centroid_Y', 'Neuron_cluster_ID'], axis = 1)


# var = data['Cell_class'] == "Astrocyte"
# x = data[var]
# var = data['Cell_class'] != "Astrocyte"
# y = data[var]

# x = x.iloc[:, [2]]
# y = y.iloc[:, [2]]

geneList = data.keys().tolist()
cell_types = ["Astrocyte", "Inhibitory", "Pericytes", "Ambiguous", "Endothelial 1",  "Excitatory", "OD Immature 1", "OD Immature 2", "Microglia", "OD Mature 2", "OD Mature 1", "Endothelial 3", "OD Mature 3", "OD Mature 4", "Endothelial 2", "Ependymal"]

markerGenes = dict()

for cell_type in cell_types:
  markerGenes[cell_type] = list()
  var = data['Cell_class'] == cell_type
  x = data[var]
  var = data['Cell_class'] != cell_type
  y = data[var]
  for gene in range(2, 157):
    newX = x.iloc[:, [gene]]
    newY = y.iloc[:, [gene]]
    U1, p = mannwhitneyu(newX, newY, alternative="greater")
    if p < .001:  
      markerGenes[cell_type].append(geneList[gene])

#print(markerGenes)

for cell_type in cell_types:
  print(cell_type)
  print(len(markerGenes[cell_type]))
  print(markerGenes[cell_type])