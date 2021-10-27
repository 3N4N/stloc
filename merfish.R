library("readxl")
require(dplyr)


dataset = read_xlsx("./data/merfish/s7.xlsx")

dataset <- dataset[ -c(1:5) ]

dataset$Centroid_X <- -1 * dataset$Centroid_X
dataset$Centroid_Y <- -1 * dataset$Centroid_Y
minX= min(dataset$Centroid_X)
minY=min(dataset$Centroid_Y)

dataset$Centroid_X <- dataset$Centroid_X - minX
dataset$Centroid_Y <- dataset$Centroid_Y - minY
dataset$new_X <- (floor(dataset$Centroid_X /50.0) *50 ) +25
dataset$new_Y <- (floor(dataset$Centroid_Y /50.0) *50) +25
dataset <- dataset %>% relocate(new_X,new_Y, .before = Centroid_X)
dataset <- dataset[order(dataset[,1], dataset[,2]),]



# write.table(dataset,"./data/temporary.csv",row.names = FALSE)


max_New_X = max(dataset$new_X)
max_New_Y = max(dataset$new_Y)


startX=25
startY=25
tempDataFrame <- dataset[0,]
tempRow <- dataset[1,]
for( k in 1: nrow(dataset))
{
    if( sqrt((dataset[k,]$Centroid_X -dataset[k,]$new_X)^2 +(dataset[k,]$Centroid_Y - dataset[k,]$new_Y )^2 ) >=20)
    {
        next
    }
    print(k)
    if(dataset[k,]$new_X == startX  && dataset[k,]$new_Y == startY)
    {

        for ( i in 7 : ncol(dataset))
        {
            tempRow[1,i] =  tempRow[1,i]+ dataset[k,i]
        } 
    }
    else{
        print(tempRow)
        tempDataFrame <- rbind(tempDataFrame,tempRow)
        tempRow <- dataset[k,]
        startX= dataset[k,]$new_X
        startY = dataset[k,]$new_Y
    }
    

}
tempDataFrame$coord<- paste( tempDataFrame$new_X ,tempDataFrame$new_Y,sep="x")
tempDataFrame <- tempDataFrame %>% select(coord,everything())
tempDataFrame = select(tempDataFrame,-new_X,-new_Y,-Centroid_X,-Centroid_Y,-Cell_class,-Neuron_cluster_ID)
write.table(tempDataFrame,"./data/moffitt/merfishSpatial.csv",row.names = FALSE,append= FALSE)
