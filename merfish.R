library("readxl")

cr = read_xlsx("./datasets/moffitt/s7.xlsx")


z = c()
for (x in cr$Cell_class) {
    # print(x)
    if      (x=="Inhibitory")   z = append(z,2)
    # else if (x=="Excitatory")   z = append(z,2)
    # else if (x=="Mature OD")    z = append(z,3)
    # else if (x=="Immature OD")  z = append(z,4)
    # else if (x=="Astrocyte")    z = append(z,5)
    # else if (x=="Microglia")    z = append(z,6)
    # else if (x=="Ependymal")    z = append(z,7)
    # else if (x=="Endothelial")  z = append(z,8)
    # else if (x=="Mural")        z = append(z,9)
    else                        z = append(z,1)
}
cr$color = z
plot(cr$Centroid_X, cr$Centroid_Y, col=cr$color)
