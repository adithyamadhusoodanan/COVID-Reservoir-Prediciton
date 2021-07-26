setwd("~/Users/adithya/COVID-Reservoir-Prediciton")

rm(list=ls())
setwd("") # Set local working directory where files are located

library(plyr)
library(h2o) # https://www.h2o.ai/products/h2o/
library(dplyr)
library(reshape2)
library(matrixStats)

`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))

localh20<-h2o.init(nthreads = -1)  # Start a local H2O cluster using nthreads = num available cores

# import and remove X and ID columns

f1<-read.csv(file="GEN_TFIDF_Fourier-3.csv",header=T)
f1 <- subset(f1, select = -c(X, ID)) # remove 

# renaming the reservoir to response and making it categorical
f1$response <-factor(f1$Reservoir)

# Remove orphans
f2<-subset(f1,f1$response!="Orphan")
f<-droplevels(f2)

# Group selection based on sample size thresholds
# Gropping any reservoir group with less than 15 sequences

t<-15 # Minimum sample size per group
s<-.7 # Proportion in the training set
host.counts<-table(f$response)
min.t <-host.counts[host.counts>=t] # minimum number of viruses per host group
f_st3<-f[f$response %in% c(names(min.t)),]
f_st3<-droplevels(f_st3) # remove empty spaces

rm(f,f2,f1)

# Evaluate patterns over many training sets
set.seed(1000000)
nloops<-1

vimps <-matrix(nrow=length(f_st3)-6-2, ncol=nloops) # create a m by n matrix where n number of rows (nloops) and m columns (number of features)
# -2 because two features are constant and will be dropped

for (i in 1:nloops){
    trains <-f_st3 %>% group_by(response) %>%  # create the training set from 70% of total data
    filter(Genbank.accession %in% sample(unique(Genbank.accession), ceiling(s*length(unique(Genbank.accession)))))
    trains<-droplevels(trains)
    f1_train <- subset(trains, select = -c(Genbank.accession,Reservoir,Viral.group,Vector.borne,Vector))
    train<-as.h2o(f1_train)
     # Identity the response column
    y <- "response"
    # Identify the predictor columns
    x <- setdiff(names(train), y)
    # Convert response to factor
    train[,y] <- as.factor(train[,y])
    # GBM with 5x cross validation of training set, test set is not used
    model1 <- h2o.gbm(x = x,
                      y = y,
                      training_frame = train,
                      ntrees = 150,
                      learn_rate = .1,
                      sample_rate = 1,
                      max_depth = 10,
                      col_sample_rate_per_tree = 1,
                      seed = 445,
                      nfolds = 5,
                      keep_cross_validation_predictions=T)
    # Retreive feature importance
    vi <- h2o.varimp(model1)
    data2  <- vi[order(vi[,1],decreasing=FALSE),] # order alphabetically
    vimps[,i]<-data2[,4] # "percentage" importance
    h2o.rm(model1)
    rm(trains,train,f1_train,vi)
}

# Average feature importance across all training sets
row.names(vimps)<-data2$variable
vimean<-rowMeans2(vimps)
visd<-rowSds(vimps)
vimps<-cbind(vimps,vimean,visd)
vistderr<-visd/sqrt(nloops)
vimps<-cbind(vimps,vimean,visd,vistderr)
vimps<- vimps[order(vimps[,nloops+1],decreasing=FALSE),] # sort by mean feature importance

# Write to file
write.csv(vimps,file="featureImportance_reservoir.csv",row.names = T)
