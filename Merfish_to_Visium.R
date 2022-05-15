#  ----------------------------------------------------------------------
#  This file contains the code to generate Visium from the MERFISH data.
#  ----------------------------------------------------------------------


library("readxl")
require(dplyr)
library(ggplot2)


bregmaAnalysis <- FALSE
if (isTRUE(bregmaAnalysis)) {
    BregmaValue <- -0.19
    inputFile <- paste("./data/merfish/Bregma/BregmaExtracted/Bregma_", as.character(BregmaValue), ".xlsx", sep = "")
} else {
    inputFile <- "./data/merfish/s7.xlsx"
}
dataset <- read_xlsx(inputFile)

for (row in 1:nrow(dataset)) {
    dataset[row, "Cell_class"] <- gsub(" ", "", dataset[row, "Cell_class"], fixed = TRUE)
}

if (isTRUE(bregmaAnalysis)) {
    visiumFile <- paste("./data/merfish/Bregma/visiumData/BregmaVisium_", as.character(BregmaValue), ".csv")
} else {
    visiumFile <- "./data/merfish/merfishVisium.csv"
}
write.csv(dataset, visiumFile, row.names = FALSE)
cellType <- unique(dataset["Cell_class"])



# Merfish plot regenerate
df <- dataset[dataset$Cell_class == "Excitatory", ]
jpeg(file = "saving_plot1.jpeg")

sp <- ggplot(df, aes(Centroid_X, Centroid_Y, colour = Cell_class)) +
    geom_point(size = 0.8) +
    xlim(-3500, -2000)
sp <- sp + scale_color_manual(values = c("grey", "brown", "yellow", "purple", "#2398cf", "red", "pink", "black", "green", "#2d7d85"))
sp
dev.off()

for (cell in cellType)
{
    dataset[cell] <- 0
}

dataset <- dataset[-c(1:5)]
dataset$Centroid_X <- -1 * dataset$Centroid_X
dataset$Centroid_Y <- -1 * dataset$Centroid_Y

minX <- min(dataset$Centroid_X)
minY <- min(dataset$Centroid_Y)


tileWidth <- 50
dataset$Centroid_X <- dataset$Centroid_X - minX
dataset$Centroid_Y <- dataset$Centroid_Y - minY
dataset$new_X <- (floor(dataset$Centroid_X / tileWidth) * tileWidth) + tileWidth / 2
dataset$new_Y <- (floor(dataset$Centroid_Y / tileWidth) * tileWidth) + tileWidth / 2
dataset <- dataset %>% relocate(new_X, new_Y, .before = Centroid_X)
dataset <- dataset[order(dataset[, 1], dataset[, 2]), ]


# write.csv(dataset, "./changed.csv")
max_New_X <- max(dataset$new_X)
max_New_Y <- max(dataset$new_Y)

dataset <- as.data.frame(dataset)


startX <- 25
startY <- 25
tempDataFrame <- dataset[0, ]
tempRow <- dataset[1, ]
tempRow[dataset[1, ]$Cell_class] <- 1
checkOutside <- function(dataRow) {
    if (sqrt((dataRow$Centroid_X - dataRow$new_X)^2 + (dataRow$Centroid_Y - dataRow$new_Y)^2) >= tileWidth * 20) {
        return(TRUE)
    }
    return(FALSE)
}
for (k in 2:nrow(dataset))
{
    # if (sqrt((dataset[k, ]$Centroid_X - dataset[k, ]$new_X)^2 + (dataset[k, ]$Centroid_Y - dataset[k, ]$new_Y)^2) >= 20) {
    if (checkOutside(dataset[k, ])) {
        next
    }
    print(k)
    if (dataset[k, ]$new_X == startX && dataset[k, ]$new_Y == startY) {
        tempRow[dataset[k, ]$Cell_class] <- tempRow[dataset[k, ]$Cell_class] + 1
        for (i in 7:ncol(dataset))
        {
            tempRow[1, i] <- tempRow[1, i] + dataset[k, i]
        }
    } else {
        # print(tempRow[(ncol(dataset) - 16):ncol(dataset)])
        tempDataFrame <- rbind(tempDataFrame, tempRow)
        while (k <= nrow(dataset)) {
            if (checkOutside(dataset[k, ])) {
                k <- k + 1
            } else {
                break
            }
        }
        tempRow <- dataset[k, ]
        tempRow[dataset[k, ]$Cell_class] <- 1
        startX <- dataset[k, ]$new_X
        startY <- dataset[k, ]$new_Y
    }
}
tempDataFrame$coord <- paste(tempDataFrame$new_X, tempDataFrame$new_Y, sep = "x")
tempDataFrame <- tempDataFrame %>% select(coord, everything())
tempDataFrame <- select(tempDataFrame, -new_X, -new_Y, -Centroid_X, -Centroid_Y, -Cell_class, -Neuron_cluster_ID)

if (isTRUE(bregmaAnalysis)) {
    outputFile <- paste("./data/merfish/Bregma/spatialData/BregmaSpatial_", as.character((BregmaValue)), ".csv", sep = "")
} else {
    outputFile <- "./data/merfish/merfishSpatial.csv"
}
write.table(tempDataFrame, outputFile, row.names = FALSE, append = FALSE)