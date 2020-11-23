---
title: "COVID19_Mexico"
author: "Fahimeh"
date: "07/13/2020"
---

library(mclust)
library(EMCluster)
library(Rtsne)
library(dendextend)
library(clValid)
library(FRESA.CAD)
library(dplyr)

##### Sources ######

source('~/clusterStability.R', echo=TRUE)
source('~/getConsensusCluster.R', echo=TRUE)
source('~/jaccard.R', echo=TRUE)
source('~/nearestCentroid.R', echo=TRUE)
source('~/pamWrapper.R', echo=TRUE)

##### Load Data ######

COVdata <- read.csv("C:/Github/COVID19/200509COVID19MEXICO.csv", na.strings="99",
                    stringsAsFactors=FALSE)


UpdatedCOVdata <- read.csv("C:/Github/COVID19/200708COVID19MEXICO.csv", na.strings="99",
                           stringsAsFactors=FALSE)

COVID19Subset <- COVdata[,c("ID_REGISTRO","SEXO","TIPO_PACIENTE","FECHA_INGRESO","FECHA_SINTOMAS"
                              ,"FECHA_DEF","INTUBADO","NEUMONIA","EDAD","EMBARAZO","DIABETES"
                              ,"EPOC","ASMA","INMUSUPR","HIPERTENSION","OTRA_COM","CARDIOVASCULAR"
                              ,"OBESIDAD","RENAL_CRONICA","TABAQUISMO","RESULTADO","UCI")]

UpdatedCOVID19Subset <- UpdatedCOVdata[,c("ID_REGISTRO","SEXO","TIPO_PACIENTE","FECHA_INGRESO","FECHA_SINTOMAS"
                            ,"FECHA_DEF","INTUBADO","NEUMONIA","EDAD","EMBARAZO","DIABETES"
                            ,"EPOC","ASMA","INMUSUPR","HIPERTENSION","OTRA_COM","CARDIOVASCULAR"
                            ,"OBESIDAD","RENAL_CRONICA","TABAQUISMO","RESULTADO","UCI")]

COVID19MexicoComplete <- COVID19Subset[complete.cases(COVID19Subset),]
UpdatedCOVID19MexicoComplete <- UpdatedCOVID19Subset[complete.cases(UpdatedCOVID19Subset),]



varnames <- c("SEXO","TIPO_PACIENTE","FECHA_INGRESO","FECHA_SINTOMAS","FECHA_DEF","INTUBADO",
              "NEUMONIA","EDAD","EMBARAZO","DIABETES","EPOC","ASMA","INMUSUPR","HIPERTENSION",
              "OTRA_COM","CARDIOVASCULAR","OBESIDAD","RENAL_CRONICA","TABAQUISMO","RESULTADO","UCI")

for (vn in varnames)
{
  print(vn)
  tb <- table(COVID19MexicoComplete[,vn])
  if (length(tb)<5)
  {
    COVID19MexicoComplete[,vn] <- 1*(COVID19MexicoComplete[,vn] == 1)
    tb <- table(COVID19MexicoComplete[,vn])
  }
}

varforDiagnosisCovid <- c("ID_REGISTRO","SEXO","TIPO_PACIENTE","FECHA_DEF","EDAD","EMBARAZO","NEUMONIA","DIABETES",
                          "EPOC","ASMA","INMUSUPR","HIPERTENSION","OTRA_COM","CARDIOVASCULAR",
                          "OBESIDAD","RENAL_CRONICA","TABAQUISMO","RESULTADO","UCI")

covidDiagnosisData <- COVID19MexicoComplete[,varforDiagnosisCovid] 
covidDiagnosisData <- subset(covidDiagnosisData,(RESULTADO == 1 & TIPO_PACIENTE == 0))
covidDiagnosisData$RESULTADO <- NULL
covidDiagnosisData$outcome <- 1*(covidDiagnosisData$FECHA_DEF != "9999-99-99")
covidDiagnosisData$FECHA_DEF <- NULL
covidDiagnosisData$TIPO_PACIENTE <- NULL

######## Updated Dataset ######

for (vn in varnames)
{
  print(vn)
  tb <- table(UpdatedCOVID19MexicoComplete[,vn])
  if (length(tb)<5)
  {
    UpdatedCOVID19MexicoComplete[,vn] <- 1*(UpdatedCOVID19MexicoComplete[,vn] == 1)
    tb <- table(UpdatedCOVID19MexicoComplete[,vn])
  }
}

UpdatedcovidDiagnosisData <- UpdatedCOVID19MexicoComplete[,varforDiagnosisCovid] 
UpdatedcovidDiagnosisData <- subset(UpdatedcovidDiagnosisData,RESULTADO == 1 & TIPO_PACIENTE == 0)
UpdatedcovidDiagnosisData$RESULTADO <- NULL
UpdatedcovidDiagnosisData$outcome <- 1*(UpdatedcovidDiagnosisData$FECHA_DEF != "9999-99-99")
UpdatedcovidDiagnosisData$FECHA_DEF <- NULL
UpdatedcovidDiagnosisData$TIPO_PACIENTE <- NULL

id <- covidDiagnosisData$ID_REGISTRO
updatedCovidData <- subset(UpdatedcovidDiagnosisData,ID_REGISTRO %in% id)
id_updated <- updatedCovidData$ID_REGISTRO

diff_ids <- setdiff(id,id_updated)

restData <- subset(covidDiagnosisData,ID_REGISTRO %in% diff_ids)

PositiveCovidMexData <- rbind(updatedCovidData,restData)

covidDiagnosisData <- PositiveCovidMexData

table(covidDiagnosisData$outcome)

rownames(covidDiagnosisData) <- 1:nrow(covidDiagnosisData)


##### clustering features #####

varforSubtypesDiscovering <- c("SEXO","EDAD","EMBARAZO","DIABETES",
                          "EPOC","ASMA","INMUSUPR","HIPERTENSION","OTRA_COM","CARDIOVASCULAR",
                          "OBESIDAD","RENAL_CRONICA","TABAQUISMO")

covidSubsDiagnosisData <- covidDiagnosisData[,varforSubtypesDiscovering] 

########## Normalize age between 0 and 1 ####

range01 <- function(x){(x-min(x))/(max(x)-min(x))}
covidSubsDiagnosisData$EDAD <- range01(covidSubsDiagnosisData$EDAD)

############ PCA Reduce Dimension ####

 PCAcomps <- prcomp(covidSubsDiagnosisData)
 summary(PCAcomps)
 
 k <- 1:13
 plot(k, type='b', PCAcomps$sdev, xlab='Number of comps',ylab='PCAs', frame=FALSE)

 pcacovidSubsDiagnosisData <- prcomp(covidSubsDiagnosisData,rank. = 6)$x 
 
########## Consensus Clustering ######
 
### positive covid
rndSamples <- sample(nrow(pcacovidSubsDiagnosisData),7000)
ConsensusData <- as.data.frame(pcacovidSubsDiagnosisData[rndSamples,])

memory.limit(size=50000) 
 
start_time <- Sys.time()

pcs2 <- clusterStability(data=ConsensusData, clustermethod=pamCluster, randomTests = 70,
                         trainFraction = 0.7,k=2)

pcs3 <- clusterStability(data=ConsensusData, clustermethod=pamCluster, randomTests = 70,
                         trainFraction = 0.7,k=3)

pcs4 <- clusterStability(data=ConsensusData, clustermethod=pamCluster, randomTests = 70,
                         trainFraction = 0.7,k=4)

pcs5 <- clusterStability(data=ConsensusData, clustermethod=pamCluster, randomTests = 70,
                         trainFraction = 0.7,k=5)

pcs6 <- clusterStability(data=ConsensusData, clustermethod=pamCluster, randomTests = 70,
                         trainFraction = 0.7,k=6)

pcs7 <- clusterStability(data=ConsensusData, clustermethod=pamCluster, randomTests = 70,
                        trainFraction = 0.7,k=7)

end_time <- Sys.time()


ccluster <- getConsensusCluster(pcs7,who="training",thr=seq(0.80,0.25,-0.05))

#ccluster <- relaxConsensusCluster(pcs7,ccluster)

mycolors <- c("red","green","blue","yellow","cyan","orange","pink","violet")

ordermatrix <- pcs7$dataConcensus

heatmapsubsample <- sample(nrow(ordermatrix),500)

orderindex <- 10*ccluster + pcs7$trainJaccardpoint

orderindex <- orderindex[heatmapsubsample]
orderindex <- order(orderindex)
ordermatrix <- ordermatrix[heatmapsubsample,heatmapsubsample]
ordermatrix <- ordermatrix[orderindex,orderindex]
rowcolors <- mycolors[1+ccluster[heatmapsubsample]]
rowcolors <- rowcolors[orderindex]

#jpeg("pam4.jpeg", units="in", width=6, height=6, res=1200)

hplot <- gplots::heatmap.2(as.matrix(ordermatrix),Rowv=FALSE,Colv=FALSE,
                           RowSideColors = rowcolors,ColSideColors = rowcolors,dendrogram = "none",
                           trace="none",main="Cluster Co-Association \n (k=7)")
#dev.off()

######### PAC #########

max.temp <- c(0,pcs3$PAC,pcs4$PAC,pcs5$PAC,pcs6$PAC,pcs7$PAC)

color_range<- c(black="#FDFC74", orange="#76FF7A", skyblue="#B2EC5D", 
                bluegreen="#FFCF4A",yellow="#F0E442", blue="#0072B2", 
                reddish="#FFB653", purplish="#CC79A7",new="#FF5733")

jpeg("PAC.jpeg", units="in", width=6, height=6, res=1200)

barplot(max.temp,xlab = "Number of cluster",
        ylab = "PAC", names.arg = c( "2","3", "4","5","6","7"), col= color_range[1:length(c(1,6,2,6,1))])

dev.off()

############ Manipulate Train & Test Data ##########

prandomSamples <- sample(nrow(covidDiagnosisData),0.7*nrow(covidDiagnosisData))

covidDiagnosisData <- covidSubsDiagnosisData

trainCOVID19Medata <- as.data.frame(covidDiagnosisData[prandomSamples,])

pcatrainCOVID19Medata <- prcomp(trainCOVID19Medata,rank. = 6)

pcatestCOVID19Medata <- predict(pcatrainCOVID19Medata,covidDiagnosisData[-prandomSamples,])

##### PAM #####

mod1 <- pamCluster(pcatrainCOVID19Medata$x,k=6)

clusterLabels <- predict(mod1,pcatestCOVID19Medata)

plot(pcatestCOVID19Medata[,1:2],col = clusterLabels$classification,main="Mclust clustering(k=4): Testing set")
plot(pcatrainCOVID19Medata$x[,1:2],col = mod1$classification,main="Mclust clustering(k=4): Training set")


names(clusterLabels$classification) <- rownames(covidDiagnosisData[-prandomSamples,])
names(mod1$classification) <- rownames(covidDiagnosisData[prandomSamples,])


######## describe Subtypes #######
library(psych)

clab <- clusterLabels$classification  ## Test
#clab <- mod1$classification

sub1 <- clab[clab==1] 
sub2 <- clab[clab==2]
sub3 <- clab[clab==3]
sub4 <- clab[clab==4]
sub5 <- clab[clab==5]
sub6 <- clab[clab==6]

sub1Data <- covidDiagnosisData[as.numeric(names(sub1)),]
sub2Data <- covidDiagnosisData[as.numeric(names(sub2)),]
sub3Data <- covidDiagnosisData[as.numeric(names(sub3)),]
sub4Data <- covidDiagnosisData[as.numeric(names(sub4)),]
sub5Data <- covidDiagnosisData[as.numeric(names(sub5)),]
sub6Data <- covidDiagnosisData[as.numeric(names(sub6)),]

#AtRiskGroup <- rbind(sub2Data,sub3Data,sub5Data)


UCI1 <- length(sub1Data$UCI[sub1Data$UCI == 1 | sub1Data$outcome == 1]) 
UCI2 <- length(sub2Data$UCI[sub2Data$UCI == 1 | sub2Data$outcome == 1]) 
UCI3 <- length(sub3Data$UCI[sub3Data$UCI == 1 | sub3Data$outcome == 1]) 
UCI4 <- length(sub4Data$UCI[sub4Data$UCI == 1 | sub4Data$outcome == 1]) 
UCI5 <- length(sub5Data$UCI[sub5Data$UCI == 1 | sub5Data$outcome == 1]) 
UCI6 <- length(sub6Data$UCI[sub6Data$UCI == 1 | sub6Data$outcome == 1]) 

per1 <- (UCI1*100)/ nrow(sub1Data)
per2 <- (UCI2*100)/ nrow(sub2Data)
per3 <- (UCI3*100)/ nrow(sub3Data)
per4 <- (UCI4*100)/ nrow(sub4Data)
per5 <- (UCI5*100)/ nrow(sub5Data)
per6 <- (UCI6*100)/ nrow(sub6Data)

def1 <- length(sub1Data$outcome[sub1Data$outcome == 1])
def2 <- length(sub2Data$outcome[sub2Data$outcome == 1])
def3 <- length(sub3Data$outcome[sub3Data$outcome == 1])
def4 <- length(sub4Data$outcome[sub4Data$outcome == 1])
def5 <- length(sub5Data$outcome[sub5Data$outcome == 1])
def6 <- length(sub6Data$outcome[sub6Data$outcome == 1])
per1 <- (def1*100)/ nrow(sub1Data)
per2 <- (def2*100)/ nrow(sub2Data)
per3 <- (def3*100)/ nrow(sub3Data)
per4 <- (def4*100)/ nrow(sub4Data)
per5 <- (def5*100)/ nrow(sub5Data)
per6 <- (def6*100)/ nrow(sub6Data)


Obs1 <- length(sub1Data$OBESIDAD[sub1Data$OBESIDAD == 1])
Obs2 <- length(sub2Data$OBESIDAD[sub2Data$OBESIDAD == 1])
Obs3 <- length(sub3Data$OBESIDAD[sub3Data$OBESIDAD == 1])
Obs4 <- length(sub4Data$OBESIDAD[sub4Data$OBESIDAD == 1])
Obs5 <- length(sub5Data$OBESIDAD[sub5Data$OBESIDAD == 1])
Obs6 <- length(sub6Data$OBESIDAD[sub6Data$OBESIDAD == 1])
per1 <- (Obs1*100)/ nrow(sub1Data)
per2 <- (Obs2*100)/ nrow(sub2Data)
per3 <- (Obs3*100)/ nrow(sub3Data)
per4 <- (Obs4*100)/ nrow(sub4Data)
per5 <- (Obs5*100)/ nrow(sub5Data)
per6 <- (Obs6*100)/ nrow(sub6Data)

Sex1 <- length(sub1Data$SEXO[sub1Data$SEXO == 0])
Sex2 <- length(sub2Data$SEXO[sub2Data$SEXO == 0])
Sex3 <- length(sub3Data$SEXO[sub3Data$SEXO == 0])
Sex4 <- length(sub4Data$SEXO[sub4Data$SEXO == 0])
Sex5 <- length(sub5Data$SEXO[sub5Data$SEXO == 0])
Sex6 <- length(sub6Data$SEXO[sub6Data$SEXO == 0])
per1 <- (Sex1*100)/ nrow(sub1Data)
per2 <- (Sex2*100)/ nrow(sub2Data)
per3 <- (Sex3*100)/ nrow(sub3Data)
per4 <- (Sex4*100)/ nrow(sub4Data)
per5 <- (Sex5*100)/ nrow(sub5Data)
per6 <- (Sex6*100)/ nrow(sub6Data)


Tab1 <- length(sub1Data$TABAQUISMO[sub1Data$TABAQUISMO == 1])
Tab2 <- length(sub2Data$TABAQUISMO[sub2Data$TABAQUISMO == 1])
Tab3 <- length(sub3Data$TABAQUISMO[sub3Data$TABAQUISMO == 1])
Tab4 <- length(sub4Data$TABAQUISMO[sub4Data$TABAQUISMO == 1])
Tab5 <- length(sub5Data$TABAQUISMO[sub5Data$TABAQUISMO == 1])
Tab6 <- length(sub6Data$TABAQUISMO[sub6Data$TABAQUISMO == 1])
per1 <- (Tab1*100)/ nrow(sub1Data)
per2 <- (Tab2*100)/ nrow(sub2Data)
per3 <- (Tab3*100)/ nrow(sub3Data)
per4 <- (Tab4*100)/ nrow(sub4Data)
per5 <- (Tab5*100)/ nrow(sub5Data)
per6 <- (Tab6*100)/ nrow(sub6Data)


Inm1 <- length(sub1Data$INMUSUPR[sub1Data$INMUSUPR == 1])
Inm2 <- length(sub2Data$INMUSUPR[sub2Data$INMUSUPR == 1])
Inm3 <- length(sub3Data$INMUSUPR[sub3Data$INMUSUPR == 1])
Inm4 <- length(sub4Data$INMUSUPR[sub4Data$INMUSUPR == 1])
Inm5 <- length(sub5Data$INMUSUPR[sub5Data$INMUSUPR == 1])
Inm6 <- length(sub6Data$INMUSUPR[sub6Data$INMUSUPR == 1])
per1 <- (Inm1*100)/ nrow(sub1Data)
per2 <- (Inm2*100)/ nrow(sub2Data)
per3 <- (Inm3*100)/ nrow(sub3Data)
per4 <- (Inm4*100)/ nrow(sub4Data)
per5 <- (Inm5*100)/ nrow(sub5Data)
per6 <- (Inm6*100)/ nrow(sub6Data)

Neu1 <- length(sub1Data$NEUMONIA[sub1Data$NEUMONIA == 1])
Neu2 <- length(sub2Data$NEUMONIA[sub2Data$NEUMONIA == 1])
Neu3 <- length(sub3Data$NEUMONIA[sub3Data$NEUMONIA == 1])
Neu4 <- length(sub4Data$NEUMONIA[sub4Data$NEUMONIA == 1])
Neu5 <- length(sub5Data$NEUMONIA[sub5Data$NEUMONIA == 1])
Neu6 <- length(sub6Data$NEUMONIA[sub6Data$NEUMONIA == 1])
per1 <- (Neu1*100)/ nrow(sub1Data)
per2 <- (Neu2*100)/ nrow(sub2Data)
per3 <- (Neu3*100)/ nrow(sub3Data)
per4 <- (Neu4*100)/ nrow(sub4Data)
per5 <- (Neu5*100)/ nrow(sub5Data)
per6 <- (Neu6*100)/ nrow(sub6Data)

Em1 <- length(sub1Data$EMBARAZO[sub1Data$EMBARAZO == 1])
Em2 <- length(sub2Data$EMBARAZO[sub2Data$EMBARAZO == 1])
Em3 <- length(sub3Data$EMBARAZO[sub3Data$EMBARAZO == 1])
Em4 <- length(sub4Data$EMBARAZO[sub4Data$EMBARAZO == 1])
Em5 <- length(sub5Data$EMBARAZO[sub5Data$EMBARAZO == 1])
Em6 <- length(sub6Data$EMBARAZO[sub6Data$EMBARAZO == 1])
per1 <- (Em1*100)/ nrow(sub1Data)
per2 <- (Em2*100)/ nrow(sub2Data)
per3 <- (Em3*100)/ nrow(sub3Data)
per4 <- (Em4*100)/ nrow(sub4Data)
per5 <- (Em5*100)/ nrow(sub5Data)
per6 <- (Em6*100)/ nrow(sub6Data)


dia1 <- length(sub1Data$DIABETES[sub1Data$DIABETES == 1])
dia2 <- length(sub2Data$DIABETES[sub2Data$DIABETES == 1])
dia3 <- length(sub3Data$DIABETES[sub3Data$DIABETES == 1])
dia4 <- length(sub4Data$DIABETES[sub4Data$DIABETES == 1])
dia5 <- length(sub5Data$DIABETES[sub5Data$DIABETES == 1])
dia6 <- length(sub6Data$DIABETES[sub6Data$DIABETES == 1])
per1 <- (dia1*100)/ nrow(sub1Data)
per2 <- (dia2*100)/ nrow(sub2Data)
per3 <- (dia3*100)/ nrow(sub3Data)
per4 <- (dia4*100)/ nrow(sub4Data)
per5 <- (dia5*100)/ nrow(sub5Data)
per6 <- (dia6*100)/ nrow(sub6Data)


asma1 <- length(sub1Data$ASMA[sub1Data$ASMA == 1])
asma2 <- length(sub2Data$ASMA[sub2Data$ASMA == 1])
asma3 <- length(sub3Data$ASMA[sub3Data$ASMA == 1])
asma4 <- length(sub4Data$ASMA[sub4Data$ASMA == 1])
asma5 <- length(sub5Data$ASMA[sub5Data$ASMA == 1])
asma6 <- length(sub6Data$ASMA[sub6Data$ASMA == 1])
per1 <- (asma1*100)/ nrow(sub1Data)
per2 <- (asma2*100)/ nrow(sub2Data)
per3 <- (asma3*100)/ nrow(sub3Data)
per4 <- (asma4*100)/ nrow(sub4Data)
per5 <- (asma5*100)/ nrow(sub5Data)
per6 <- (asma6*100)/ nrow(sub6Data)


hip1 <- length(sub1Data$HIPERTENSION[sub1Data$HIPERTENSION == 1])
hip2 <- length(sub2Data$HIPERTENSION[sub2Data$HIPERTENSION == 1])
hip3 <- length(sub3Data$HIPERTENSION[sub3Data$HIPERTENSION == 1])
hip4 <- length(sub4Data$HIPERTENSION[sub4Data$HIPERTENSION == 1])
hip5 <- length(sub5Data$HIPERTENSION[sub5Data$HIPERTENSION == 1])
hip6 <- length(sub6Data$HIPERTENSION[sub6Data$HIPERTENSION == 1])
per1 <- (hip1*100)/ nrow(sub1Data)
per2 <- (hip2*100)/ nrow(sub2Data)
per3 <- (hip3*100)/ nrow(sub3Data)
per4 <- (hip4*100)/ nrow(sub4Data)
per5 <- (hip5*100)/ nrow(sub5Data)
per6 <- (hip6*100)/ nrow(sub5Data)

car1 <- length(sub1Data$CARDIOVASCULAR[sub1Data$CARDIOVASCULAR == 1])
car2 <- length(sub2Data$CARDIOVASCULAR[sub2Data$CARDIOVASCULAR == 1])
car3 <- length(sub3Data$CARDIOVASCULAR[sub3Data$CARDIOVASCULAR == 1])
car4 <- length(sub4Data$CARDIOVASCULAR[sub4Data$CARDIOVASCULAR == 1])
car5 <- length(sub5Data$CARDIOVASCULAR[sub5Data$CARDIOVASCULAR == 1])
car6 <- length(sub6Data$CARDIOVASCULAR[sub6Data$CARDIOVASCULAR == 1])
per1 <- (car1*100)/ nrow(sub1Data)
per2 <- (car2*100)/ nrow(sub2Data)
per3 <- (car3*100)/ nrow(sub3Data)
per4 <- (car4*100)/ nrow(sub4Data)
per5 <- (car5*100)/ nrow(sub5Data)
per6 <- (car6*100)/ nrow(sub6Data)


ren1 <- length(sub1Data$RENAL_CRONICA[sub1Data$RENAL_CRONICA == 1])
ren2 <- length(sub2Data$RENAL_CRONICA[sub2Data$RENAL_CRONICA == 1])
ren3 <- length(sub3Data$RENAL_CRONICA[sub3Data$RENAL_CRONICA == 1])
ren4 <- length(sub4Data$RENAL_CRONICA[sub4Data$RENAL_CRONICA == 1])
ren5 <- length(sub5Data$RENAL_CRONICA[sub5Data$RENAL_CRONICA == 1])
ren6 <- length(sub6Data$RENAL_CRONICA[sub6Data$RENAL_CRONICA == 1])
per1 <- (ren1*100)/ nrow(sub1Data)
per2 <- (ren2*100)/ nrow(sub2Data)
per3 <- (ren3*100)/ nrow(sub3Data)
per4 <- (ren4*100)/ nrow(sub4Data)
per5 <- (ren5*100)/ nrow(sub5Data)
per6 <- (ren6*100)/ nrow(sub6Data)

epo1 <- length(sub1Data$EPOC[sub1Data$EPOC == 1])
epo2 <- length(sub2Data$EPOC[sub2Data$EPOC == 1])
epo3 <- length(sub3Data$EPOC[sub3Data$EPOC == 1])
epo4 <- length(sub4Data$EPOC[sub4Data$EPOC == 1])
epo5 <- length(sub5Data$EPOC[sub5Data$EPOC == 1])
epo6 <- length(sub6Data$EPOC[sub6Data$EPOC == 1])
per1 <- (epo1*100)/ nrow(sub1Data)
per2 <- (epo2*100)/ nrow(sub2Data)
per3 <- (epo3*100)/ nrow(sub3Data)
per4 <- (epo4*100)/ nrow(sub4Data)
per5 <- (epo5*100)/ nrow(sub5Data)
per6 <- (epo6*100)/ nrow(sub6Data)


describe(sub1Data$EDAD)
describe(sub2Data$EDAD)
describe(sub3Data$EDAD)
describe(sub4Data$EDAD)
describe(sub5Data$EDAD)
describe(sub6Data$EDAD)

########### ANOVA Test #########

Data1 <- as.data.frame(cbind(sub1Data, class=1))
Data2 <- as.data.frame(cbind(sub2Data, class=2))
Data3 <- as.data.frame(cbind(sub3Data, class=3))
Data4 <- as.data.frame(cbind(sub4Data, class=4))
Data5 <- as.data.frame(cbind(sub5Data, class=5))
Data6 <- as.data.frame(cbind(sub6Data, class=6))
Data <- rbind(Data1,Data2,Data3,Data4,Data5,Data6) 

class <- as.numeric(Data$class)

aov_cont<- aov(class ~ Data$EDAD)
summary(aov_cont)[[1]][["Pr(>F)"]][[1]]

############# chi-square test #########

tbl_sex = table(Data$SEXO ,Data$class)
print(chisq.test(tbl_sex))

tbl_UCI = table(Data$UCI ,Data$class)
print(chisq.test(tbl_UCI))

tbl_DEF = table(Data$outcome ,Data$class)
print(chisq.test(tbl_DEF))

tbl_obs = table(Data$OBESIDAD ,Data$class)
print(chisq.test(tbl_obs))

tbl_tab = table(Data$TABAQUISMO ,Data$class)
print(chisq.test(tbl_tab))

tbl_inm = table(Data$INMUSUPR ,Data$class)
print(chisq.test(tbl_inm))

tbl_nue = table(Data$NEUMONIA ,Data$class)
print(chisq.test(tbl_nue))

tbl_emb = table(Data$EMBARAZO ,Data$class)
print(chisq.test(tbl_emb))

tbl_dia = table(Data$DIABETES ,Data$class)
print(chisq.test(tbl_dia))

tbl_asma = table(Data$ASMA ,Data$class)
print(chisq.test(tbl_asma))

tbl_hip = table(Data$HIPERTENSION ,Data$class)
print(chisq.test(tbl_hip))

tbl_car = table(Data$CARDIOVASCULAR ,Data$class)
print(chisq.test(tbl_car))

tbl_ren = table(Data$RENAL_CRONICA ,Data$class)
print(chisq.test(tbl_ren))

tbl_epo = table(Data$EPOC ,Data$class)
print(chisq.test(tbl_epo))


######## rpart algorithm ########
library(rpart)
library(rpart.plot)

# outcome classes
Data$class[Data$class == 1] <- "subtype 1"
Data$class[Data$class == 2] <- "subtype 5"
Data$class[Data$class == 3] <- "subtype 4"
Data$class[Data$class == 4] <- "subtype 2"
Data$class[Data$class == 5] <- "subtype 6"
Data$class[Data$class == 6] <- "subtype 3"

names(Data)[names(Data) == "UCI"] <- "ICU"
names(Data)[names(Data) == "EDAD"] <- "Age"
names(Data)[names(Data) == "NEUMONIA"] <- "Pneumonia"
names(Data)[names(Data) == "OBESIDAD"] <- "Obesity"
names(Data)[names(Data) == "SEXO"] <- "Sex"
names(Data)[names(Data) == "HIPERTENSION"] <- "Hypertension"
names(Data)[names(Data) == "DIABETES"] <- "Diabetes"
names(Data)[names(Data) == "TABAQUISMO"] <- "smoking"



fit <- rpart(class ~ . , method="class", data=Data,control=rpart.control(minsplit=10, minbucket=3, cp=0.02))#,control=rpart.control(minsplit=10, minbucket=3, cp=0.009)) #
             
rpart.plot(fit, main=" Positive COVID19 Subtypes ",yesno=2) 

