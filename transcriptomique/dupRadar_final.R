#!/usr/bin/env Rscript
### Date : (creation) Wed Apr 14 10:27:31 2021 ; (last modified) Fri Apr 16 16:29:45 2021
### Author : Ombeline LAMER
### Purpose : Identification des BAM comportant un taux anormal de duplicats
### citation("dupRadar") #acceder aux crédits
### sessionInfo() # paramètre de la session R pour ce script

### Packages nécessaires au script
library(dupRadar)

### Recuperation des arguments
args = commandArgs(trailingOnly=TRUE)
setwd(dir = args[1])
CPU=as.integer(args[2])
fichierBam=list.files(pattern= "dupAllmarked.bam")
print(fichierBam)
### Creation d'un fichier de sauvegarde (header)

#modeleFit=file("../markedDupReport/duplicateFit.txt", open = "a")
#close(modeleFit)
write(paste("nomBase\tintercept\tcoefIntercept\tslope\tcoefSlope\tAIC"),file="../markedDupReport/duplicateFit.txt",append = FALSE)

###Analyses séquentielles des fichiers
for (BAM in fichierBam) {
  
  ### Nom racine du fichier = nom echantillon
  nomBase=strsplit(BAM,".dupAllmarked.bam",TRUE)[[1]]
  
  ### Comptage avec Rsubread / featureCounts
  dupMat<-analyzeDuprates(bam=BAM, gtf = "../BP_corr_nomChr.gff", stranded = 2, paired = TRUE, threads = CPU,
                        verbose =TRUE, GTF.featureType="gene", GTF.attrType= c("ID"))

  ### Plotting ###
  #Creation et ouverture d'un PDF pour l'echantillon pour les sorties graphiques
  pdf(paste("../markedDupReport/",nomBase,"duplicate.pdf"),height=10,width=10)
  par(mfrow=c(3,2))
  
  cumulativeDuprateBarplot(dupMat) #argument optionnel : stepSize=0.05 (default)
  title(main="",sub=paste("Read count cumulé de",nomBase))
  
  duprateExpBoxplot(dupMat)
  title(paste("Taux de duplication ~ RPK de",nomBase))
  #où RPK = total Read Count per Kilobase
  
  
  #An interesting quality metric are the fraction of reads taken up by groups of genes
  #binned by 5% expression levels.
  readcountExpBoxplot(dupMat)

  #help in identifying skewed distributions with unusual amount of lowly
  #expressed genes, or to detect no consensus between replicates
  expressionHist(dupMat)
  
  
  duprateExpDensPlot(dupMat)
  title(paste("Taux de duplication ~ comptage total des reads de",nomBase))
  
  #Localisation manuelle // via la souris
  #duprateExpIdentify(dupMat, idCol = "ID")
  
  duprateExpPlot(dupMat)
  title(paste("Taux de duplication ~ comptage total des reads de",nomBase))
  
  dev.off() #fermeture pdf
  
  ### Modelisation
  
  #Obtention du modèle utilisé par duprateExpDensPlot
  #duprateExpFit(dupMat) #renvoie le summary de la modelisation utilisée pour le scatterplot
  fit <- duprateExpFit(dupMat) #modèle logarithmique ajusté aux données
  
  #Interprétation rapide
  #duprate at low read counts: fit$intercept
  #progression of the duplication rate: fit$slope
  
  #Sauvegarde dans un fichier avec header
  write(paste(nomBase,"\t",fit$intercept,"\t",fit$glm$coefficients[[1]],"\t",fit$slope,"\t",fit$glm$coefficients[[2]],"\t",fit$glm$aic),
        file="../markedDupReport/duplicateFit.txt",append=TRUE)
  
} #fin de l'iteration sur les échantillons


### Fin du script