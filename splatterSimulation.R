library("splatter")
library("scater")
library("ggplot2")
set.seed(1)
sce <- mockSCE()
params <- newSplatParams()
# counts <- counts(sce)
# counts
sim <- splatSimulate(params, batchCells = c(1000, 1000,1000,1000), verbose = FALSE,nGenes=1000)


sim <- logNormCounts(sim)
simPCA <- runPCA(sim)
mat<- counts(sim)
mat = t(mat)
p<-colData(sim)
df_mtx <- as.data.frame(mat)
df_mtx["cellType"]<-as.list(p['Batch'])
# mat<-  cbind(mat,Class= as.list(p['Batch'])
r<-rownames(mat)
write.csv(df_mtx, file = "splatterCSV.csv", append = FALSE, quote = TRUE, sep = ",")

plotPCA(simPCA, colour_by = "Batch")

cellTypes <-as.list(unique(df_mtx["cellType"]))

typeWiseDf<- c()
for ( type in cellTypes)
{
    x<-df_mtx[df_mtx["cellType"]==type,]
    typeWiseDf<- append(typeWiseDf,x)
}
for ( i in 1:20)
{
    for( j in 1:20)
    {
        if( 1<=i &&  i <=10 && 1<=j && j<=10)
        {

            "cellTypeA"
            
        }
        else if (11<= i  &&  i <=20  && 1<=j && j<=10) {
           "cellTypeB"
        }
        else if (1<= i  && i <=10  &&  11<=j && j<=20) {
           "cellTypeC"
        }
        else if(11<=i && i<=20  &&  11<=j && j<=20){
            "cellType D"
        }

    }
}