#!/usr/bin/env Rscript
### Date : Thu Jul 01 11:37:39 2021 ; (last modified)
### Author : Ombeline LAMER
### Purpose : Analyses statistiques des transcripts RNAseq par DESeq2
### citation("DESeq2") #acceder aux credits
### sessionInfo() # parametre de la session R pour ce script
### Usage : RScript DESeq2 -p cheminVersComptages -o cheminVersSortie

rm(list=ls())

### Packages necessaires au script

library(DESeq2)
library(magrittr)
library(apeglm)
library(genefilter)
library(pheatmap)

### Constantes d'environnement

args = commandArgs(trailingOnly=TRUE) #args[1] = 1er parametre apres nomFichier

#parametres : valeurs par defaut
chemin="./"

#valeurs specifiques si specifiees
index=which(regexpr("-p|--chemin",args)==1)+1
if (isEmpty(index)==F){
  chemin=args[index]
}

###############################################################################
# Quel.le.s sont les hypotheses et tests que l'on veut poser ?
# 1 ) Quels genes sont DE entre A2, B1 et C respectivement (dont A2 vs C )
# 2 ) Quels genes sont DE entre les souches et leur traitement AST2 AST17
# 3 ) Quels genes sont DE entre les souches et leur repiquage (chgmt ph?notype)
# 4 ) Quels genes sont DE entre AST2 et AST17 (avec et sans effet souche)
# 5 ) Quels genes sont DE entre les souches resistantes vs les non resistantes
# 5a) toutes souches confondues 5b) au sein des souches
###############################################################################

### Extraction des donnees d'interet 

## Recuperation de la totalite des donnees
sampleFiles <- regexpr("A|B|C",list.files(chemin))
sampleFiles <- list.files(chemin)[which(sampleFiles==1)]

sampleStrain <- c(rep("A",6*4),rep("B",6*4),rep("C",6))
sampleTreatment <- c( rep(
  c(rep("sans",6),rep("AST17",6),rep("AST2",6),rep("sans",6)), 2),
  rep("sans",6) )
sampleSubculture=c( rep(c(rep("0",6*3),rep("5",6)),2), rep("0",6) )

sampleTable <- data.frame(sampleName = as.character(strsplit(sampleFiles,".txt")),
                          fileName = sampleFiles,
                          strain = sampleStrain,
                          treatment = sampleTreatment,
                          subCulture = sampleSubculture)
sampleTable$strain <- factor(sampleTable$strain)
sampleTable$treatment <- factor(sampleTable$treatment)
sampleTable$subCulture <- factor(sampleTable$subCulture)

## Ordonner les facteurs pour le design en mettant les referents en premier
# sampleTable$strain # mettre "B" en referent
# sampleTable$treatment # mettre le "sans"
# sampleTable$subCulture # RAS : 0 en referent
sampleTable$strain %<>% relevel("B")
sampleTable$treatment %<>% relevel("sans")

## Adaptation de la liste des donnees suivant l'etude souhaitee

# etude de l'effet souche : (1)
ST_strainEffect <- sampleTable[c(1:6,25:30,49:54),1:3]
ST_strainEffect_comp <- sampleTable[c(1:6,49:54),1:3]
ST_strainEffect_comp$strain <-factor(ST_strainEffect_comp$strain)
ST_strainEffect_all <-sampleTable[,1:3]

# etude de l'effet traitement : (2) + (4?)
ST_TreatmentEffect_A2 <- sampleTable[c(1:18),c(1,2,4)]
ST_TreatmentEffect_B1 <- sampleTable[c(25:42),c(1,2,4)]
ST_TreatmentEffect_AB <- sampleTable[c(1:18,25:42),1:4]

# etude de l'effet repiquage : (3)
ST_subCultureEffect_A2 <- sampleTable[c(1:6,19:24),c(1,2,5)]
ST_subCultureEffect_B1 <- sampleTable[c(25:30,43:48),c(1,2,5)]
ST_subCultureEffect_AB <- sampleTable[c(1:6,19:24,25:30,43:48),c(1:3,5)]
# etude de la resistance des souches : (5)
  ### resistante vs. sensible (a faire si temps)

## Construction de DESeqDataSet (design simple)
# NB : 5807 genes pour 54 echantillons

# design= ~ strain + treatment + subCulture + # strain:treatment)
dds_strainEffect <- DESeqDataSetFromHTSeqCount(sampleTable = ST_strainEffect,
                                  directory = chemin,
                                  design = ~ strain)
ST_straintEffect_comp<-ST_strainEffect
ST_strainEffect_comp$strain = relevel(factor(ST_strainEffect_comp$strain),"C")
dds_strainEffect_comp <- DESeqDataSetFromHTSeqCount(sampleTable = ST_strainEffect_comp,
                                               directory = chemin,
                                               design = ~ strain)
dds_strainEffect_all <- DESeqDataSetFromHTSeqCount(sampleTable = ST_strainEffect_all,
                                               directory = chemin,
                                               design = ~ strain)

dds_ttA2 <- DESeqDataSetFromHTSeqCount(sampleTable = ST_TreatmentEffect_A2,
                                               directory = chemin,
                                               design = ~ treatment)
dds_ttB1 <- DESeqDataSetFromHTSeqCount(sampleTable = ST_TreatmentEffect_B1,
                                       directory = chemin,
                                       design = ~ treatment)
dds_ttAB <- DESeqDataSetFromHTSeqCount(sampleTable = ST_TreatmentEffect_AB,
                                       directory = chemin,
                                       design = ~ strain + treatment + strain:treatment)

dds_scA2 <- DESeqDataSetFromHTSeqCount(sampleTable = ST_subCultureEffect_A2,
                                       directory = chemin,
                                       design = ~ subCulture)
dds_scB1 <- DESeqDataSetFromHTSeqCount(sampleTable = ST_subCultureEffect_B1,
                                       directory = chemin,
                                       design = ~ subCulture)
dds_scAB <- DESeqDataSetFromHTSeqCount(sampleTable = ST_subCultureEffect_AB,
                                       directory = chemin,
                                       design = ~ strain + subCulture)

### Filtration des donnees (OUI !)
#au moins 6 echantillons (6 repliquats par ech) avec count >= 10

keep <- rowSums(counts(dds_strainEffect)>=10)>=6
dds_strainEffect<-dds_strainEffect[keep,] #4989
dds_strainEffect_comp<-dds_strainEffect_comp[keep,] #4989
keep <- rowSums(counts(dds_strainEffect_all)>=10)>=6
dds_strainEffect_all<-dds_strainEffect_all[keep,] #5141

keep <- rowSums(counts(dds_ttA2)>=10)>=6
dds_ttA2<-dds_ttA2[keep,] #4856
keep <- rowSums(counts(dds_ttB1)>=10)>=6
dds_ttB1<-dds_ttB1[keep,] #5006
keep <- rowSums(counts(dds_ttAB)>=10)>=6
dds_ttAB<-dds_ttAB[keep,] #5086


keep <- rowSums(counts(dds_scA2)>=10)>=6
dds_scA2<-dds_scA2[keep,] #4774
keep <- rowSums(counts(dds_scB1)>=10)>=6
dds_scB1<-dds_scB1[keep,] #4950
keep <- rowSums(counts(dds_scAB)>=10)>=6
dds_scAB<-dds_scAB[keep,] #5047

### Differential Expression Pipeline processing

de_strainEffect <- DESeq(dds_strainEffect) #cas : strain_A_vs_B et strain_C_vs_B
de_strainEffect_comp <- DESeq(dds_strainEffect_comp) #pour le strain_A_vs_C
#resultsNames(de_strainEffect)
#resultsNames(de_strainEffect_comp)



### Construction de la table des résultats (chaque expérience/DE)

################################################################################
# interpretation statistique et choix parametriques
# 
# log2 du Fold Change : ex. log2FC si log2FC=1 ou 1.5 alors expression genique
# resp. de 2**1=2 (doublee) et 2**1.5=2.82 (plus que doublee)
# On souhaite poser une limite a log2FC=1. Or probleme des incertitudes = 
# standard error (=res$lfcSE) qui bruite/masque les resultat significatifs.
# suivant la SE, il faudra peut-etre changer le FC seuil. 
# + si lfc>1 alors on ne pourra pas voir les genes qui ne sont pas double
# mais bouge qd meme !! => seuil a 0.5 ? soit 2**0.5 = 1.414
# 
# La p-value du test (les infos nous permettent de conclure que chmgt contre l'
# hypothese nulle = aucune DE du gene) represente la probabilite qu'un FC donne
# puisse etre cause par une situation ou la DE serait du a autre chose que la
# condition experimentale testee (ex. B vs B-AST2, la condition est AST2)
# Il nous faut donc une p-value seuil d'autant plus petite pour etre plus fiable
# la stringence sur la pval ajustee (+ petite), fait aussi baisser la FDR ( =
# False Dicovery Rate) i.e. le nombre de faux positifs.
################################################################################

# on pose seuilFC=0 (pour voir meme les petites variations,  quitte a avoir du bruit)
# et seuilPval=0.01 (ajuste selon Benjamini-Hochberg)
# PERMET d'avoir un nombre de gène DE relativement large

#Modele
res <- results(de_dds,
               contrast = c("condition","lvl_test","lvl_ref"),
               lfcThreshold = 0, # 0 by default
               alpha = 0.1, #0.1 by default stringence sur les padj
               parallel = T 
               )
## Table des résultats
res_SE_AvsB <- results(de_strainEffect,
                  contrast = c("strain","A","B"),
                  lfcThreshold = 0.5,
                  )#parallel = T)
res_SE_CvsB <- results(de_strainEffect,
                       lfcThreshold = 0.5,
                       parallel = T)
res_SE_AvsC <- results(de_strainEffect_comp,
                       lfcThreshold = 0.5,
                       parallel = T)
head(res_SE_AvsB) #verifier que les nom de l'acquision sont bons
head(res_SE_CvsB)
head(res_SE_AvsC)



### Obtenir la liste des genes les plus differentiellement exprimes

#NB : il faut faire un subset qd meme !!
res_SE_AvsC <- subset(res_SE_AvsC,padj<0.1)
down=head(res_SE_AvsC[order(res_SE_AvsC$log2FoldChange,decreasing = F),],50) #plus fort down-reg
up=head(res_SE_AvsC[order(res_SE_AvsC$log2FoldChange,decreasing = T),],50) #plus fort up-reg
res_SE_AvsC
up
down

################################################################################
###
### Plotter les resultats, CR en PDF et images individuelles
###
################################################################################

### ### A AUTOMATISER EN SORTIE !!! ### ###


### MA Plot

resultsNames(de_dds) #avec shrinkage
MA <- lfcShrink(de_dds, coef="strain_A_vs_B", type="apeglm")
plotMA(MA,ylim=c(-5,5))


##########
# NOTA : le MA plot mets en evidence le bruit "de Poisson". En effet, les genes
# aux comptes faibles et de DE significative ne peuvent etre distingue du bruit
# (les points gris au centre du plot ). 
##########


### Gene clustering -> A poursuivre pour voir si correspondance GO
topVarGenes <- head(order(rowVars(assay(de_dds)), decreasing = TRUE), 10)
mat  <- assay(de_dds)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(de_dds)[, c("strain","treatment","subCulture")])
pheatmap(mat, annotation_col = anno)

### Enrichment !! B. pseudomallei not in usual databasis... need to build our CSV

### Exportation des tables de resultats en csv !
# pour comparaison avec edgeR et construction de diagramme de Venn

write.csv(res,file="testvstest.csv") #POSER LES SEUILS !!
