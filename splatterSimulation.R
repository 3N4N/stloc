#  ----------------------------------------------------------------------
#  Generates Visium data from splatter-generated single-cell data.
#  ----------------------------------------------------------------------

library("splatter")
library("scater")
library("ggplot2")

set.seed(1)

sce <- mockSCE()
params <- newSplatParams()
# counts <- counts(sce)
# counts
sim <- splatSimulate(params, batchCells = c(300, 300, 300, 300, 300, 300),
                     verbose = FALSE,nGenes=200)


sim <- logNormCounts(sim)
simPCA <- runPCA(sim)
mat<- counts(sim)
mat = t(mat)
p<-colData(sim)
df_mtx <- as.data.frame(mat)
df_mtx["cellType"]<-as.list(p['Batch'])
# mat<-  cbind(mat,Class= as.list(p['Batch'])
r<-rownames(mat)
write.csv(df_mtx, file = "data/splatter/splatterCSV.csv", quote = TRUE)

plotPCA(simPCA, colour_by = "Batch")

cellTypes <- unique(df_mtx["cellType"])

#typeWiseDf
typeWiseDf <- vector(mode = "list", length = length(cellTypes))

# for ( type  in cellTypes)
for (i in 1: length(cellTypes[[1]]))
{
    # print(cellTypes[[1]][i])
    x<-df_mtx[ df_mtx["cellType"]==cellTypes[[1]][i] ,]
    # print(x)
    typeWiseDf[[i]]<-x
    # k=k+1
}


shape = 30
for (i in 1:shape)
{
    for( j in 1:shape)
    {
        if( 1<=i &&  i <=shape/2 && 1<=j && j<=shape/2)
        {
            "cellTypeA"
            s<-seq(2)
            tempDf <- typeWiseDf[[1]] #selecting subset of type A
            tempDf <- tempDf[1:2,] #selecting first two rows
            tempDf <- subset(tempDf,select=-c(ncol(tempDf))) # deleting last column of type name
            ro<-colSums(tempDf) # sum of the cell counts

            typeWiseDf[[1]]<-typeWiseDf[[1]][-s,] #deleting first two rows that were selected
            # typeWiseDf[[1]]<- subset(typeWiseDf[[1]],select = -c(ncol(typeWiseDf[[1]])))
            print(tempDf)
        }
        # else if (11<= i  &&  i <=20  && 1<=j && j<=10) {
        #    "cellTypeB"
        # }
        # else if (1<= i  && i <=10  &&  11<=j && j<=20) {
        #    "cellTypeC"
        # }
        # else if(11<=i && i<=20  &&  11<=j && j<=20){
        #     "cellType D"
        # }

    }
}
