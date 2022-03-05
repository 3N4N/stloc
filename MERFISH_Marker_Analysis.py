#  ----------------------------------------------------------------------
#  This contains the code for determing the marker genes of merfish data
#  using the Mann-Whitneyu test.
#  ----------------------------------------------------------------------


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

data = pd.read_csv("./data/merfish/merfishVisium.csv")

# datatop = data.head()
# print(datatop)

data = data.drop(['Centroid_Y', 'Animal_ID', 'Animal_sex', 'Behavior',
                  'Bregma', 'Centroid_X', 'Centroid_Y', 'Neuron_cluster_ID'], axis=1)


# var = data['Cell_class'] == "Astrocyte"
# x = data[var]
# var = data['Cell_class'] != "Astrocyte"
# y = data[var]

# x = x.iloc[:, [2]]
# y = y.iloc[:, [2]]

geneList = data.keys().tolist()
cell_types = ["Astrocyte", "Inhibitory", "Pericytes", "Ambiguous", "Endothelial 1",
              "Excitatory", "OD Immature 1", "OD Immature 2", "Microglia", "OD Mature 2",
              "OD Mature 1", "Endothelial 3", "OD Mature 3", "OD Mature 4", "Endothelial 2",
              "Ependymal"]

markerGenes = dict()
markerGene_for_cell_types_top10 = pd.DataFrame(
    columns=['cell_type', 'marker_gene', 'p_value'])
tempo = pd.DataFrame(columns=['cell_type', 'marker_gene', 'p_value'])


cellTypeCSV = list()
markerGeneCSV = list()
pvalCSV = list()


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
            cellTypeCSV.append(cell_type)
            markerGeneCSV.append(geneList[gene])
            pvalCSV.append(p)
            markerGenes[cell_type].append(geneList[gene])
    temp = pd.DataFrame(list(zip(cellTypeCSV, markerGeneCSV, pvalCSV)), columns=[
                        'cell_type', 'marker_gene', 'p_value'])
    temp = temp.sort_values('p_value')
    # if(len(temp.index) > 10):
    #     temp = temp.iloc[0:10, :]
    cellTypeCSV.clear()
    markerGeneCSV.clear()
    pvalCSV.clear()
    markerGene_for_cell_types_top10 = pd.concat(
        [markerGene_for_cell_types_top10, temp], axis=0)

# print(markerGenes)

# ependymal-inhibitory
llen = list()
for x in range(20):
    llen.append("EpendymalInhibitory")

tempo = markerGene_for_cell_types_top10[markerGene_for_cell_types_top10['cell_type'] == "Inhibitory"]
tempo1 = markerGene_for_cell_types_top10[markerGene_for_cell_types_top10['cell_type'] == "Ependymal"]
tempo = pd.concat([tempo, tempo1], axis=0)
tempo["cell_type"] = llen
markerGene_for_cell_types_top10 = pd.concat(
    [markerGene_for_cell_types_top10, tempo], axis=0)

#markerGene_for_cell_types = pd.DataFrame(list(zip(cellTypeCSV, markerGeneCSV, pvalCSV)), columns=['cell_type', 'marker_gene', 'p_value'])

markerGene_for_cell_types_top10.to_csv(
    './data/merfish/markerGene_for_merfish_data.csv', index=False
)

# for cell_type in cell_types:
#   var = data[data['cell_type'] == cell_type]

# for cell_type in cell_types:
#   print(cell_type)
#   print(len(markerGenes[cell_type]))
#   print(markerGenes[cell_type])
