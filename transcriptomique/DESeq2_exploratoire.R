#!/usr/bin/env Rscript
### Date : Tue Jun 29 16:48:48 2021 ; last modified : Wed Jul 21 11:20:47 2021
### Author : Ombeline LAMER
### Purpose : Visualisation pour echantillons RNAseq via DESeq2
### citation("DESeq2") #acceder aux credits
### sessionInfo() # parametre de la session R pour ce script
### Usage : RScript DESeq2_exploratoire -p cheminVersComptages -o cheminVersSortie

rm(list=ls())
args = commandArgs(trailingOnly=TRUE) #args[1] = 1er parametre apres nomFichier

### Affichage d'une aide
library(S4Vectors)
aide=regexpr("-h|--help",args)
if (isEmpty(aide)==F && isEmpty(which(aide==1))==F ){
  write("\n\t### AIDE ### \n\t# Analyses exploratoires automatisees pour le jeu de donnees ABC
  COMMANDE : R DESeq2_exporatoire.R -p ./ -o ./
  --path | -p\t: chemin vers les fichiers de comptage. Defaut : './'
  --output | -o\t: chemin vers le repertoire ou seront placees les sorties. Defaut : './' 
  --help | -h\t: affiche cette aide"
        , stdout())
  quit("no",1) #sortie
}

### Packages necessaires au script

library(DESeq2)
library(vsn)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(PoiClaClu)
library(glmpca)

### Constante d'environnement

#valeurs par defaut des parametres
chemin="./"
output="./"

#valeurs specifiees
  #index=grep("-p",args) #path
index=which(regexpr("-p|--path",args)==1)+1
if (isEmpty(index)==F){
  chemin=args[index]
}
index=which(regexpr("-o|--output",args)==1)+1
if (isEmpty(index)==F){
  output=args[index]
}

### test chemin
sampleFiles <- regexpr("A|B|C",list.files(chemin))
sampleFiles <- list.files(chemin)[which(sampleFiles==1)]
if (isEmpty(sampleFiles)==T){
  write("Erreur : Aucun comptage n'a ete trouve.
        Veuillez modifier le code R ou rectifier le chemin saisi :
        ex. R DESeq2_exploratoire.R -p ./chemin/vers/comptages ",stdout())
  quit("no",1) #sortie
}

### Sauvegarde dans un pdf
pdf(paste(output,"/DESeq2_exploratoire.pdf", sep = "/"),width = 10 , height = 10,
    title="Analyses exploratoires des donnees transcriptomiques, souches A2 B1 et C")

### Extraction des donnees d'interet

sampleFiles <- regexpr("A|B|C",list.files(chemin))
sampleFiles <- list.files(chemin)[which(sampleFiles==1)]

sampleStrain <- c(rep("A",6*4),rep("B",6*4),rep("C",6))
sampleTreatment <- c( rep(
  c(rep("sans",6),rep("AST17",6),rep("AST2",6),rep("sans",6)), 2),
  rep("sans",6) )
sampleSubculture=c(
  rep(c(rep("0",6*3),rep("5",6)),2),
  rep("0",6) )
sampleTable <- data.frame(sampleName = as.character(strsplit(sampleFiles,".txt")),
                          fileName = sampleFiles,
                          strain = sampleStrain,
                          treatment = sampleTreatment,
                          subCulture = sampleSubculture )
sampleTable$strain <- factor(sampleTable$strain)
sampleTable$treatment <- factor(sampleTable$treatment)
sampleTable$subCulture <- factor(sampleTable$subCulture)

#construction de DESeqDataSet (design simple)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = chemin,
                                       design = ~ strain + treatment +
                                         subCulture )


### Exploration des donnees et visualisation
  # NB1 : des transformations sont effectuees pour VISUALISER
  # NB2 : les RAW seront repris pour les STATISTIQUES des DE des g?nes

## Filtration des donn?es

#au moins 5 echantillons avec un count >= 15
keep <- rowSums(counts(ddsHTSeq)>=15)>=5
dds_filtered<-ddsHTSeq[keep,]
#nrow(dds_filtered) #5000 au lieu de 5807

## Transformation VST vs RLOG

  #VST : Variance Stabilizing Transformation
  #rlog : regularized-logarithm transformation

#From DESeq2 Manual, M. Love et al. 2019
#In RNAseq, the expected variance grows with the mean
#i.e. the higher the counts, the more weight they have in PCA.

# Strat 1 : log+1(pseudocount) des counts normalis?s
# biais : les genes de faibles comptes bruitent de mani?re plus importante !
#=> Autres strat : VST et regularized : ecrase la variance des faibles comptes
# a priori, vst qd n>30 et rlog qd n<30 ; ou n le nombre d'echantillon

# comparons leur performance sur notre ?chantillon
vsd<-vst(dds_filtered,blind = F)
rld<-rlog(dds_filtered,blind = F) #bien plus long ...
#Avec blind=F => prise en compte de l'effet souche et/ou l'effet traitement
meanSdPlot(assay(vsd),ranks=F)
meanSdPlot(assay(rld),ranks=F)
#CCL : les 2 methode sont a peu pres equivalente pour notre etude
#=> on gardera VST par la suite, car plus rapide a computater

#visualisation entre les differentes methodes de transformation

dds<-estimateSizeFactors(dds_filtered)
df<-bind_rows(
  as.data.frame(log2(counts(dds,normalized=T)[,1:2]+1)) %>% mutate(transformation="log2(x +1)"),
  as.data.frame(assay(vsd)[,1:2]) %>% mutate(transformation="vst"),
  as.data.frame(assay(rld)[,1:2]) %>% mutate(transformation="rlog"))
colnames(df)[1:2]<-c("x","y")
lvls<-c("log2(x +1)","vst","rlog")
df$transformation<-factor(df$transformation,levels=lvls)

ggplot(df,aes(x=x,y=y)) + geom_hex(bins=80) +
  coord_fixed() + facet_grid( . ~ transformation)
#CCL : On peut choisir indifferement entre VST et RLOG


## Distance entre echantillons

#construction de la matrice de distance
sampleDists_VSD <- dist(t(assay(vsd)))
sampleDists_rlog <- dist(t(assay(rld))) #pour prouver que ne change rien

#vers la visualisation
#VSD
cat_df = data.frame("Souche" = vsd$strain,
                    "Traitement" = vsd$treatment,
                    "Repiquage" = vsd$subCulture)
row.names(cat_df) = row.names(colData(vsd))
ann_colors = list(
  Souche = c(A="orangered1", B="deepskyblue",C="chartreuse2"),
  Traitement = c(AST17 = "orange", AST2 = "purple",sans = "grey" ),
  Repiquage = c("0" = "white", "5" = "darkred") )

sampleDistMatrix_vsd <- as.matrix( sampleDists_VSD )
colnames(sampleDistMatrix_vsd) <- colnames(vsd)
rownames(sampleDistMatrix_vsd) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix_vsd,
         clustering_distance_rows = sampleDists_VSD,
         clustering_distance_cols = sampleDists_VSD,
         col = colors, #les noms sont bien au bon endroit
         main = "Heatmap echantillons~souche, VST",
         cutree_cols = 8, cutree_rows= 3,
         annotation_col = cat_df,
         annotation_colors = ann_colors )

#Rlog
cat_df = data.frame("Souche" = rld$strain,
                    "Traitement" = rld$treatment,
                    "Repiquage" = rld$subCulture)
row.names(cat_df) = row.names(colData(rld))
ann_colors = list(
  Souche = c(A="orangered1", B="deepskyblue",C="chartreuse2"),
  Traitement = c(AST17 = "orange", AST2 = "purple",sans = "grey" ),
  Repiquage = c("0" = "white", "5" = "darkred") )

sampleDistMatrix_rlog <- as.matrix( sampleDists_rlog )
colnames(sampleDistMatrix_rlog) <- colnames(rld)
rownames(sampleDistMatrix_rlog) <- NULL
pheatmap(sampleDistMatrix_rlog,
         clustering_distance_rows = sampleDists_rlog,
         clustering_distance_cols = sampleDists_rlog,
         col = colors, #les noms sont bien au bon endroit
         main = "Heatmap echantillons~souche, rlog",
         cutree_cols = 8, cutree_rows= 3,
         annotation_col = cat_df,
         annotation_colors = ann_colors )

### DONNE LA MEME HEATMAP ! d'ou l'interet de ne travailler qu'avec 1 transfo

#############################################################
#derni?re strategie Poisson Distance : s?pare par dissimilarit? entre comptes
poisd <- PoissonDistance(t(counts(dds_filtered)))
samplePoisDistMatrix <- as.matrix( poisd$dd )

colnames(samplePoisDistMatrix) <- colnames(dds_filtered)
rownames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors, main="Heatmap echantillons~souche, dist de Poisson",
         cutree_cols = 8, cutree_rows= 3,
         annotation_col = cat_df,
         annotation_colors = ann_colors )
#si la strategie de la distance de poisson permet de trancher entre 2 groupes
#son utilisation ne permet pas une mise en evidence de fine variation intra-ech.
#############################################################


## PCA plot
plotPCA(vsd, intgroup = c("strain","treatment","subCulture"))+
  ggtitle("PCA, transformation VST")
plotPCA(rld, intgroup = c("strain","treatment","subCulture"))+
  ggtitle("PCA, transformation rlog")
#donne (quasi) la meme chose

##################################
#si double ou plus conditions
pcaData<-plotPCA(vsd, intgroup = c("strain","treatment","subCulture")
                 ,returnData=T)

grouplabel=c(rep(c(T, rep("",5)),9))
index=which(grouplabel=="TRUE")
hLabel=c(2,rep(0,5),  -1.5,rep(0,5),  0.5,rep(0,5),  3,rep(0,5),
         0,rep(0,5),  1,rep(0,5),  1,rep(0,5),  1.5,rep(0,5),
         2,rep(0,5) )
vLabel=c(5,rep(0,5),  2,rep(0,5),  2,rep(0,5),  1,rep(0,5),
         2,rep(0,5),  -3,rep(0,5),  2,rep(0,5),  4,rep(0,5),
         2,rep(0,5) )
grouplabel[index]=as.character(pcaData$strain[index])

ggplot(pcaData, aes(x = PC1, y = PC2, color = treatment, shape = subCulture)) +
  geom_point(size =3) +
  #xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  #ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +  ggtitle("PCA with VST data") + geom_text(aes(
    label = grouplabel, hjust = hLabel, vjust = vLabel,fontface = "bold"))

##################################

## PCA modelisation generalisee

gpca <- glmpca(counts(dds), L=2)
gpca.dat <- gpca$factors
gpca.dat$strain <- dds$strain
gpca.dat$subCulture <- dds$subCulture
gpca.dat$treatment <- dds$treatment

grouplabel=c(rep(c(T, rep("",5)),9))
index=which(grouplabel=="TRUE")
grouplabel[index]=as.character(gpca.dat$strain[index])

ggplot(gpca.dat, aes(x = dim1, y = dim2, color = treatment, shape = subCulture)
       ) + geom_point(size =3) + coord_fixed() + 
  ggtitle("glmpca - Generalized PCA") + geom_text(aes(label = grouplabel,
                                                      hjust = 0, vjust = -2))

## MDS plot
mds <- as.data.frame(colData(vsd))  %>%
  cbind(cmdscale(sampleDistMatrix_vsd))
ggplot(mds, aes(x = `1`, y = `2`, color = treatment, shape = subCulture)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with VST data") +
  geom_text(aes(label = grouplabel, hjust = -1.5, vjust = 0.5))

mdsPois <- as.data.frame(colData(dds)) %>%
  cbind(cmdscale(samplePoisDistMatrix))
ggplot(mdsPois, aes(x = `1`, y = `2`, color = treatment, shape = subCulture)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with PoissonDistances") +
  geom_text(aes(label = grouplabel, hjust = 2, vjust = -0.5))

###############################################################################
### Heatmap et PCA (VST) pour diff?rents sous-?chantillonage
###############################################################################

### Extraction des donn?es d'int?ret
  # On ne prend pas la souche C
sampleFiles <- regexpr("A|B",list.files(chemin))
sampleFiles <- list.files(chemin)[which(sampleFiles==1)]

sampleStrain <- c(rep("A",6*4),rep("B",6*4))
sampleTreatment <- c( rep(
  c(rep("sans",6),rep("AST17",6),rep("AST2",6),rep("sans",6)), 2))
sampleSubculture=c( rep(c(rep("0",6*3),rep("5",6)),2) )
sampleTable <- data.frame(sampleName = as.character(strsplit(sampleFiles,".txt")),
                          fileName = sampleFiles,
                          strain = sampleStrain,
                          treatment = sampleTreatment,
                          subCulture = sampleSubculture )
sampleTable$strain <- factor(sampleTable$strain)
sampleTable$treatment <- factor(sampleTable$treatment)
sampleTable$subCulture <- factor(sampleTable$subCulture)

ddsHTSeq_noC <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = chemin,
                                       design = ~ strain + treatment +
                                         subCulture )
  #Ni C ni repiquage

sampleFiles <- regexpr("A2-[1-9A]|B1-[1-9A]",list.files(chemin))
sampleFiles <- list.files(chemin)[which(sampleFiles==1)]

sampleStrain <- c(rep("A",6*3),rep("B",6*3))
sampleTreatment <- c( rep( c(rep("sans",6),rep("AST17",6),rep("AST2",6)), 2))
sampleTable <- data.frame(sampleName = as.character(strsplit(sampleFiles,".txt")),
                          fileName = sampleFiles,
                          strain = sampleStrain,
                          treatment = sampleTreatment )
sampleTable$strain <- factor(sampleTable$strain)
sampleTable$treatment <- factor(sampleTable$treatment)

ddsHTSeq_noCnoR <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                           directory = chemin,
                                           design = ~ strain + treatment )


### Exploration des donnees et visualisation

## Filtration des donn?es : au moins 5 echantillons avec un count >= 15
keep <- rowSums(counts(ddsHTSeq_noC)>=15)>=5
ddsnoC_filtered<-ddsHTSeq_noC[keep,] #4987 au lieu de 5807
keep <- rowSums(counts(ddsHTSeq_noCnoR)>=15)>=5
ddsnoCnoR_filtered<-ddsHTSeq_noCnoR[keep,] #4940 au lieu de 5807
  #Visualisation de la dispersion des donn?es => du bruit
vsd_noC<-vst(ddsnoC_filtered,blind = F)
vsd_noCnoR<-vst(ddsnoCnoR_filtered,blind = F)
meanSdPlot(assay(vsd_noC),ranks=F)
meanSdPlot(assay(vsd_noCnoR),ranks=F)


## Distance entre echantillons

#construction des matrices de distance
sampleDists_VSD_noC <- dist(t(assay(vsd_noC)))
sampleDists_VSD_noCnoR <- dist(t(assay(vsd_noCnoR)))

#vers la visualisation
  #noC
cat_df = data.frame("Souche" = vsd_noC$strain,
                    "Traitement" = vsd_noC$treatment,
                    "Repiquage" = vsd_noC$subCulture)
row.names(cat_df) = row.names(colData(vsd_noC))
ann_colors = list(
  Souche = c(A="orangered1", B="deepskyblue"),
  Traitement = c(AST17 = "orange", AST2 = "purple",sans = "grey" ),
  Repiquage = c("0" = "white", "5" = "darkred") )

sampleDistMatrix_vsd_noC <- as.matrix( sampleDists_VSD_noC )
colnames(sampleDistMatrix_vsd_noC) <- colnames(vsd_noC)
rownames(sampleDistMatrix_vsd_noC) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix_vsd_noC,
         clustering_distance_rows = sampleDists_VSD_noC,
         clustering_distance_cols = sampleDists_VSD_noC,
         col = colors,
         main = "Heatmap echantillons~souche, VST",
         cutree_cols = 8, cutree_rows= 2,
         annotation_col = cat_df,
         annotation_colors = ann_colors )

  #noCnoR
cat_df = data.frame("Souche" = vsd_noCnoR$strain,
                    "Traitement" = vsd_noCnoR$treatment)
row.names(cat_df) = row.names(colData(vsd_noCnoR))
ann_colors = list(
  Souche = c(A="orangered1", B="deepskyblue"),
  Traitement = c(AST17 = "orange", AST2 = "purple",sans = "white" ) )

sampleDistMatrix_vsd_noCnoR <- as.matrix( sampleDists_VSD_noCnoR )
colnames(sampleDistMatrix_vsd_noCnoR) <- colnames(vsd_noCnoR)
rownames(sampleDistMatrix_vsd_noCnoR) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix_vsd_noCnoR,
         clustering_distance_rows = sampleDists_VSD_noCnoR,
         clustering_distance_cols = sampleDists_VSD_noCnoR,
         col = colors,
         main = "Heatmap echantillons~souche, VST",
         cutree_cols = 6, cutree_rows= 2,
         annotation_col = cat_df,
         annotation_colors = ann_colors )

#Poisson Distance : s?pare par dissimilarit? entre comptes
poisd_noC <- PoissonDistance(t(counts(ddsnoC_filtered)))
poisd_noCnoR <- PoissonDistance(t(counts(ddsnoCnoR_filtered)))

samplePoisDistMatrix_noC <- as.matrix( poisd_noC$dd )
samplePoisDistMatrix_noCnoR <- as.matrix( poisd_noCnoR$dd )

#graphique heatmap
  #noC
cat_df = data.frame("Souche" = vsd_noC$strain,
                    "Traitement" = vsd_noC$treatment,
                    "Repiquage" = vsd_noC$subCulture)
row.names(cat_df) = row.names(colData(vsd_noC))
ann_colors = list( "Souche" = c(A="orangered1", B="deepskyblue"),
  "Traitement" = c(AST17 = "orange", AST2 = "purple",sans = "grey" ),
  "Repiquage" = c("0" = "white", "5" = "darkred") )

rownames(samplePoisDistMatrix_noC) <- NULL
colnames(samplePoisDistMatrix_noC) <- colnames(vsd_noC)
pheatmap(samplePoisDistMatrix_noC,
         clustering_distance_rows = poisd_noC$dd,
         clustering_distance_cols = poisd_noC$dd,
         col = colors, main="Heatmap echantillons~souche, dist. Poisson",
         cutree_cols = 8, cutree_rows= 2,
         annotation_col = cat_df,
         annotation_colors = ann_colors )

  #noCnoR
cat_df = data.frame("Souche" = vsd_noCnoR$strain,
                    "Traitement" = vsd_noCnoR$treatment)
row.names(cat_df) = row.names(colData(vsd_noCnoR))
ann_colors = list( Souche = c(A="orangered1", B="deepskyblue"),
  Traitement = c(AST17 = "orange", AST2 = "purple",sans = "white" ) )

colnames(samplePoisDistMatrix_noCnoR) <- colnames(ddsnoCnoR_filtered)
rownames(samplePoisDistMatrix_noCnoR) <- NULL
pheatmap(samplePoisDistMatrix_noCnoR,
         clustering_distance_rows = poisd_noCnoR$dd,
         clustering_distance_cols = poisd_noCnoR$dd,
         col = colors, main="Heatmap echantillons~souche, dist. Poisson",
         cutree_cols = 6, cutree_rows= 2,
         annotation_col = cat_df,
         annotation_colors = ann_colors )

## PCA plot
plotPCA(vsd_noC, intgroup = c("strain","treatment","subCulture"))+
  ggtitle("PCA - VST, w/o C")
plotPCA(vsd_noCnoR, intgroup = c("strain","treatment"))+
  ggtitle("PCA, PCA - VST, w/o C or Subculture")

##################################
pcaData_noC<-plotPCA(vsd_noC, intgroup = c("strain","treatment","subCulture")
                 ,returnData=T)
pcaData_noCnoR<-plotPCA(vsd_noCnoR, intgroup = c("strain","treatment")
                     ,returnData=T)

  #noC
grouplabel=c(rep(c(T, rep("",5)),8))
index=which(grouplabel=="TRUE")
hLabel=c(2,rep(0,5),  -1.5,rep(0,5),  0.5,rep(0,5),  3,rep(0,5),
         0,rep(0,5),  1,rep(0,5),  1,rep(0,5),  1.5,rep(0,5) )
vLabel=c(5,rep(0,5),  2,rep(0,5),  2,rep(0,5),  1,rep(0,5),
         2,rep(0,5),  -2,rep(0,5),  2,rep(0,5),  3,rep(0,5) )
grouplabel[index]=as.character(pcaData_noC$strain[index])
ggplot(pcaData_noC, aes(x = PC1, y = PC2, color = treatment, shape = subCulture)) +
  geom_point(size =3) +
  coord_fixed() +  ggtitle("PCA with VST data") + geom_text(aes(
    label = grouplabel, hjust = hLabel, vjust = vLabel,fontface = "bold"))

  #noCnoR - plutot inutile
grouplabel=c(rep(c(T, rep("",5)),6))
index=which(grouplabel=="TRUE")
hLabel=c(2,rep(0,5),  -1.5,rep(0,5),  0.5,rep(0,5),  3,rep(0,5),
         0,rep(0,5),  1,rep(0,5) )
vLabel=c(5,rep(0,5),  2,rep(0,5),  2,rep(0,5),  1,rep(0,5),
         2,rep(0,5),  -2,rep(0,5) )
grouplabel[index]=as.character(pcaData_noCnoR$strain[index])
ggplot(pcaData_noCnoR, aes(x = PC1, y = PC2, color = treatment, shape = strain 
                           )) + geom_point(size =3) +
  coord_fixed() +  ggtitle("PCA with VST data")


## PCA modelisation generalisee

  #noC
gpca_noC <- glmpca(counts(ddsnoC_filtered), L=2)
gpca.dat_noC <- gpca_noC$factors
gpca.dat_noC$strain <- ddsnoC_filtered$strain
gpca.dat_noC$subCulture <- ddsnoC_filtered$subCulture
gpca.dat_noC$treatment <- ddsnoC_filtered$treatment

grouplabel=c(rep(c(T, rep("",5)),8))
index=which(grouplabel=="TRUE")
grouplabel[index]=as.character(gpca.dat_noC$strain[index])

ggplot(gpca.dat_noC, aes(x = dim1, y = dim2, color = treatment,
                         shape = subCulture) ) +
  geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA") +
  geom_text(aes(label = grouplabel, hjust = 0, vjust = -1.5))

  #noCnoR
gpca_noCnoR <- glmpca(counts(ddsnoCnoR_filtered), L=2)
gpca.dat_noCnoR <- gpca_noCnoR$factors
gpca.dat_noCnoR$strain <- ddsnoCnoR_filtered$strain
gpca.dat_noCnoR$subCulture <- ddsnoCnoR_filtered$subCulture
gpca.dat_noCnoR$treatment <- ddsnoCnoR_filtered$treatment

grouplabel=c(rep(c(T, rep("",5)),6))
index=which(grouplabel=="TRUE")
grouplabel[index]=as.character(gpca.dat_noCnoR$strain[index])

ggplot(gpca.dat_noCnoR, aes(x = dim1, y = dim2, color = treatment,
                            shape = strain ) ) +
  geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA")

### AJOUTER LA SORTIE DES GRAPHIQUES !!
dev.off()
#Fin script