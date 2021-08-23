#!/bin/bash
# Auteur : Ombeline LAMER
# Date creation : 2021/08/12
# Description : Pipeline traitement RNAseq pour data "RNAseq_ABC" à partir de l'alignement à l'assemblage B1
# Usage : bash RNAseq_pipeline_B1.sh -r repertoireRaw -w repertoireSortie -m maxThread

###
### Lecture de la ligne de commande
###

# Options par défaut

#cheminRaw=/media/NextSeq/Fabienne/RNAseq_Fab/1-raw
#workDir=/media/NextSeq/RNASeq_Ombeline_2021
seed=2021

# Recuperation des options specifiees en ligne de commande
while [ "$1" != "" ]; do
    case $1 in
        -s | --seed )		shift
        				seed="$1"
                           		;;
        -r | --rawData )		shift
        				cheminRaw="$1"
					;;
        -w | --workDir )		shift
        				workDir="$1"
					;;
        -m | --maxThread )		shift
        				maxThread="$1"
					;;
        -t | --threadFastQC )		shift
        				tFastQC="$1"
					;;
        -c | --coreCutadapt )		shift
        				coreCut="$1"
					;;
        -i | --threadIndexBowtie )		shift
        				threadIndexBowtie="$1"
					;;
        -a | --threadAlignBowtie )		shift
        				threadAlignBowtie="$1"
					;;
		-R | --coreR )		shift
        				coreR="$1"
					;;
        -h | --help )		shift
        				aide="affichage"
					;;
    esac
    shift
done

# Verification - existence des paramètres
if [ "$cheminRaw" ]
then
	echo Chemin vers les raw : $cheminRaw
else
	echo Chemin vers les données brutes non trouvées. Préciser avec --rawData.
	exit
fi

if [ "$workDir" ]
then
	echo Dossier de sortie : $workDir
else
	echo Veuillez préciser le dossier où vous souhaitez vos sorties. Préciser avec --workDir.
	exit
fi

if [ "$maxThread" ]
then
	echo Nombre max. de CPU disponible : $maxThread
else
	echo Veuillez préciser le nombre maximum de CPU/cores disponibles. Préciser avec --maxThread.
	exit
fi

if [ "$aide" ]
then
	more help_RNAseq.txt #ATTENTION CHEMIN
	exit
fi

if [ "$tFastQC" ]
then
	echo CPU pour FastQC : $tFastQC
else
    tFastQC=$(($maxThread/4)) #10
fi

if [ "$coreCut" ]
then
	echo CPU pour Cutadapt : $coreCut
else
    coreCut=$(($maxThread/2)) #30 #nombre de coeur pour cutadapt (sans parallelisation, mais repartition de la tache sur les coeurs)
fi

if [ "$threadIndexBowtie" ]
then
	echo CPU pour indexation avec Bowtie2 : $threadIndexBowtie
else
    threadIndexBowtie=$maxThread #40
fi

if [ "$threadAlignBowtie" ]
then
	echo CPU pour alignement avec Bowtie2 : $threadAlignBowtie
else
    threadAlignBowtie=$(($maxThread/2)) #20
fi

if [ "$coreR" ]
then
	echo CPU pour statistiques avec R : $coreR
else
    coreR=$(($maxThread/2)) #20
fi


###
### Initialisation
###

# Creation du tableau des noms de fichiers à traiter
cd $cheminRaw
ls -f *.fastq | parallel gzip  #compresse rapidement tous les fastq qui ne le seraient pas
fichiers=(`ls -S *_R*.fastq.gz`)
#Nom racine des fichiers (i.e. echantillons) d'après la liste des raw
nomBase=(`ls -S *.fastq.gz | cut -d _ -f 1 | uniq`)
#ATTENTION le dossier des raw ne contient que les echantillons a traiter

cd $workDir

# Creation des dossiers de rangement des résultats, si nécessaire
[ ! -d "workDir/mapping_B1" ] && mkdir -p "workDir/mapping_B1"
[ ! -d "workDir/mapReport_B1" ] && mkdir -p "workDir/mapReport_B1"
[ ! -d "workDir/count_B1" ] && mkdir -p "workDir/count"
[ ! -d "workDir/analyseStat" ] && mkdir -p "workDir/analyseStat" #a voir si besoin ajout B1

###
### Traitement des fichiers
###

cd $workDir/mapping_B1

## Construction de l'index
#NB : On utilise le fichier fasta issu de l'assemblage de l'isolat B1
bowtie2-build $workDir/B1.fa ref-B1 \
	--threads $threadIndexBowtie \
	--seed $seed

## Alignement
parallel -k "bowtie2 --time \
	--very-sensitive \
	--threads $threadAlignBowtie \
	--seed $seed \
	--fr \
	-x ref-B1 \
	-1 $workDir/trimmed/{}_R1_trimmed.fastq.gz \
	-2 $workDir/trimmed/{}_R2_trimmed.fastq.gz \
	-S {}.sam &> $workDir/mapReport_B1/{}_time.txt" ::: ${nomBase[@]}

cd $workDir

#NOTA a supprimer conda genomic_bact

## QC alignement Bowtie2 par multiQC
multiqc -i B1_54PE_RNAseq_Samples_ \
	-b "108 fichiers de reads ont été trimmés avec cutadapt. Ce sont donc 54 couples de reads (PE) qui sont alignés avec Bowtie2 sur la séquence de l'isolat B1" \
	-n B1_allAligned \
	-o ${workDir}/mapReport_B1 \
	--interactive ${workDir}/mapReport_B1

#newList
cd $workDir/mapping_B1
nomBase=(`ls -S *.sam | cut -d . -f 1`) #actualisation
#cd $workDir

## Conversion SAM en BAM
parallel -k "samtools view {}.sam -o {}.bam" ::: ${nomBase[@]}

## Sorting des reads alignés dans les bam
parallel -k "samtools sort {}.bam -o {}.bam" ::: ${nomBase[@]}

## Indexage des BAM ordonnés
parallel -k "samtools index {}.bam" ::: ${nomBase[@]}

## Statistiques sur les BAM
parallel -k "samtools flagstat {}.bam > $workDir/mapReport_B1/{}.flagstat.txt && \
	samtools idxstat {}.bam > $workDir/mapReport_B1/{}.idxstat.txt && \
	samtools stats {}.bam > $workDir/mapReport_B1/${}.stats.txt" ::: ${nomBase[@]}

## multiQC sur les stats
multiqc -i B1_BAM_54*PE_RNAseq_Samples \
	-b "108 fichiers de reads ont été trimmés avec cutadapt. Ce sont donc 54 couples de reads (PE) qui sont alignés avec Bowtie2 sur la séquence de l'isolat B1. Ils sont ensuite convertis, triés et analysés avec SAMtools" \
	-n B1_sortedBAM \
	-o ${workDir}/mapReport_B1 \
	--interactive ${workDir}/mapReport_B1

### Observation du taux de duplication avec script R 
# meme situation que pour K96243. ressort le meme resultat (meme conclusion)

### Comptage des reads
cd $workDir/count_B1

parallel -k -j27 '"tseq-count -f bam \
	--order pos \
	--stranded reverse \
	--type gene \
	--idattr locus_tag \
	-m intersection-nonempty \
	{}.bam \
	${workDir}/BP_corr_nomChr.gff > {}.txt" ::: ${nomBase[@]}

# Fin du script
