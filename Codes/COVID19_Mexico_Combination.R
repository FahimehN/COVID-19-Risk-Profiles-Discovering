---
title: "COVID19_Mexico_Combination"
author: "Fahimeh"
date: "02/08/2020"
---
  
knitr::opts_chunk$set(echo = TRUE)

#load(file="~/PositivCovidDiagnosisData.RDATA")
load(file="~/NegativeCovidDiagnosisData.RDATA")

covidDiagnosisData <- NegcovidDiagnosisData

varforSubtypesDiscovering <- c("SEXO","EDAD","EMBARAZO","DIABETES",
                               "EPOC","ASMA","INMUSUPR","HIPERTENSION","OTRA_COM","CARDIOVASCULAR",
                               "OBESIDAD","RENAL_CRONICA","TABAQUISMO")

PositiveCovidDiagnosisData <- covidDiagnosisData[,varforSubtypesDiscovering] 

COVID_Males_20_40 <-  subset(PositiveCovidDiagnosisData,SEXO == 0 & EDAD>20 & EDAD <=40)
COVID_Males_40_60 <-  subset(PositiveCovidDiagnosisData,SEXO == 0 & EDAD>40 & EDAD <=60)
COVID_Males_60_Older <-  subset(PositiveCovidDiagnosisData,SEXO == 0 & EDAD>60 )

COVID_Females_20_40 <-  subset(PositiveCovidDiagnosisData,SEXO == 1 & EDAD>20 & EDAD <=40)
COVID_Females_40_60 <-  subset(PositiveCovidDiagnosisData,SEXO == 1 & EDAD>40 & EDAD <=60)
COVID_Females_60_Older <-  subset(PositiveCovidDiagnosisData,SEXO == 1 & EDAD>60 )


  #COVIDData <- COVID_Males_20_40
  #COVIDData <- COVID_Males_40_60
  # COVIDData <- COVID_Males_60_Older

  #COVIDData <- COVID_Females_20_40
  #COVIDData <- COVID_Females_40_60
  #COVIDData <- COVID_Females_60_Older


Conditions <- c("DIABETES","EPOC","ASMA","INMUSUPR","HIPERTENSION","OTRA_COM","CARDIOVASCULAR","OBESIDAD","RENAL_CRONICA","TABAQUISMO","EMBARAZO")
Conditions1 <- c("DIAB","EPOC","ASMA","IMU","HYPR","OTR DIS","CARD","OBES","RENAL","SMOK","PREG")
 
 
COVIDData <- COVIDData[,Conditions]

rownames(COVIDData) <- 1:nrow(COVIDData)

# ConditionsEnglish <- c("Diabetes","EPOC","Asthma","Immunosuppression","Hypertension","Other diseases",
#                        "Cardiovascular","Obesity","Chronic kidney","Smoking","Pregnancy")


##### 1 combination of COVID19 Elements ####

countCOndition <- numeric(length(Conditions))

names(countCOndition) <- Conditions1 #ConditionsEnglish

m=1
for (m in 1:nrow(COVIDData))
{ 
  currentConditions <- COVIDData[m,Conditions]
  currentConditions[is.na(currentConditions)] <- 2
  currentConditions <- currentConditions == 1

  for (n in 1:length(Conditions))
  {
    whichones <- Conditions %in% Conditions[n]
    countCOndition[n] <- countCOndition[n] + 1*(sum(xor(whichones[n],currentConditions[n])) == 0)
  }
}

print(countCOndition)
which.max(countCOndition)
perMxConditions <- (countCOndition[order(-countCOndition)] / nrow(COVIDData))*100

##### 2 combination of COVID19 Elements #####

tableConditions <- combn(Conditions,2)
Conditionsnames2 <- combn(Conditions1,2)

countTwoCountCondition <- numeric(ncol(tableConditions))
combConditions <- numeric(2) == FALSE

m=1
for (m in 1:nrow(COVIDData))
{
  currentConditions <- COVIDData[m,Conditions]
  currentConditions[is.na(currentConditions)] <- 2
  currentConditions <- currentConditions == 1

  for (n in 1:ncol(tableConditions))
  {
    whichones <- tableConditions[,n]
    countTwoCountCondition[n] <- countTwoCountCondition[n] + 1*(sum(xor(combConditions,currentConditions[,whichones])) == 0)
  }
  
}

names(countTwoCountCondition) <- paste(Conditionsnames2[1,],Conditionsnames2[2,],sep="_")
print(countTwoCountCondition)
tableConditions[,which.max(countTwoCountCondition)]

perMxConditions <- (countTwoCountCondition[order(-countTwoCountCondition)] / nrow(COVIDData))*100


##### 3 combination of COVID19 Elements #####

tableConditions <-combn(Conditions,3)
Conditionsnames3 <- combn(Conditions1,3)

countThreeCountCondition <- numeric(ncol(tableConditions))
combConditions <- numeric(3) == FALSE

for (m in 1:nrow(COVIDData))
{
  currentConditions <- COVIDData[m,Conditions]
  currentConditions[is.na(currentConditions)] <- 2
  currentConditions <- currentConditions == 1
  
  for (n in 1:ncol(tableConditions))
  {
    whichones <- tableConditions[,n]
    countThreeCountCondition[n] <- countThreeCountCondition[n] + 1*(sum(xor(combConditions,currentConditions[,whichones])) == 0)
  }
}

names(countThreeCountCondition ) <- paste(Conditionsnames3[1,],Conditionsnames3[2,],Conditionsnames3[3,],sep="_")

print(countThreeCountCondition )

tableConditions[,which.max(countThreeCountCondition )]

perMxConditions <- (countThreeCountCondition[order(-countThreeCountCondition)] / nrow(COVIDData))*100
perMxConditions[1]
###### show ########

CountConsitions <- c(countCOndition,countTwoCountCondition,countThreeCountCondition)

CountConsitions <- CountConsitions[order(-CountConsitions)]

## min to max
color_range<- c(grange="#ffb74a",bluegreen="#FFCF4A",yellow="#F0E442",yellow2="#edf042",
                yellowgreen="#d3f042",green1="#9ff042",green2="#6ef042",green3="#42f065",
                green4="#51f595",green5="#76FF7A")

## max to min
# color_range<- c( green5="#76FF7A",green4="#51f595",green3="#42f065",green2="#6ef042",green1="#9ff042",
#                  yellow1="#d3f042",yellow2="#edf042",yellow3="#F0E442",yellow4="#FFCF4A",yellow5="#ffb74a")


barplot((CountConsitions[1:10]),#nrow(COVIDData))*100
        main = "Males of 60 to older aged group \n (N = 3989)", #older aged group
        xlab = "The tolal of participants",
        col= color_range,
        xlim = c( 0 , 2000 ),
        beside = TRUE,
        #cex.names=1.5, 
        cex.axis=1,
        plot = TRUE, axis.lty = 0, offset = 0,
        las=2,cex.names=0.5,horiz = TRUE) 


