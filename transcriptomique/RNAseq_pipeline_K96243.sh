#!/bin/bash
# Auteur : Ombeline LAMER
# Date creation : 2021/03/29
# Description : Pipeline traitement RNAseq pour data "RNAseq_ABC"
# Usage : bash RNAseq_pipeline_K96243.sh -r repertoireRaw -w repertoireSortie -m maxThread

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
#ATTENTION mettre dans le repertoire les données d'interet seulement (pas de undetermined)

cd $workDir

# Creation des dossiers de rangement des résultats, si nécessaire
[ ! -d "workDir/rawReport" ] && mkdir -p "workDir/rawReport"
[ ! -d "workDir/trimmed" ] && mkdir -p "workDir/trimmed"
[ ! -d "workDir/trimmedReport" ] && mkdir -p "workDir/trimmedReport"
[ ! -d "workDir/markedDup" ] && mkdir -p "workDir/markedDup"
[ ! -d "workDir/markedDupReport" ] && mkdir -p "workDir/markedDupReport"
[ ! -d "workDir/mapping" ] && mkdir -p "workDir/mapping"
[ ! -d "workDir/mapReport" ] && mkdir -p "workDir/mapReport"
[ ! -d "workDir/count" ] && mkdir -p "workDir/count"
[ ! -d "workDir/analyseStat" ] && mkdir -p "workDir/analyseStat"

###
### Traitement des fichiers
###

### QC report sur les fichiers brut (raw)

## QC report par fastQC
parallel -X fastqc ${cheminRaw}/{} -o ${workDir}/rawReport/ -t $tFastQC ::: ${fichiers[@]} #voir si le -X marche

## QC global report par multiQC
multiqc -i 108_RNAseq_Samples \
	-n allSamplesAnalysis \
	-o ${workDir}/rawReport/globalAnalysis \
	-p ${workDir}/rawReport #obtenir les graphiques figés en png et svg
multiqc -i 108_RNAseq_Samples \
	-n allSamplesInteractive \
	-o ${workDir}/rawReport/globalAnalysis \
	--interactive ${workDir}/rawReport #version graphique interactive

### Trimming par Cutadapt

#OPTIONS CHOISIES :
# -a 	: 	R1_3end ANCHORED et Illumina ScriptSeq adapter  
# -A 	: 	R2_3end ANCHORED et Illumina SriptSeq adapter
# --times : nombre de fois où adapter peut être retire (1x en 3' + si jamais accroché une 2e fois)
# -q	: 	qualité min des fin en 5' et 3' // FOLLOWING OPTIONS APPLIED AFTER TRIMMING
# --trim-n :trimme les N terminaux
# --max-n : abandonne les reads ayant plus de 10% de N dans leur séquence
# -m	:	longueur minimale post-traitement = 20 (R1 et R2)

# -k : la parallelisation conserve l'ordre de liste arguments
parallel -k "cutadapt --cores=$coreCut \ 
	-a 64nt=AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG$ \
	-a R1-IAC=AGATCGGAAGAGCACACGTCTGAAC \
	-A 58nt=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT$ \
	-A R2-IAC=AGATCGGAAGAGCGTCGTGTAGGGA \
	--times=2 \
	-q 20,20 \
	--trim-n \
	--max-n=0.1 \
	-m 20:20 \
	-o $workDir/trimmed/${nomBase}_R1_trimmed.fastq.gz \
	--paired-output=$workDir/trimmed/${nomBase}_R2_trimmed.fastq.gz \
	$cheminRaw/${nomBase}_R1.fastq.gz \
	$cheminRaw/${nomBase}_R2.fastq.gz > $workDir/trimmedReport/${nomBase}.txt" ::: ${nomBase[@]}


#constitution de la liste des fichiers trimmés
cd $workDir/trimmed
trimFile=(`ls -S *_R*_trimmed.fastq.gz | uniq`)
#actualisation des nomBase au cas où ordre change
nomBase=(`ls -S *_R*_trimmed.fastq.gz | cut -d _ -f 1 | uniq`)
cd $workDir

### QC report sur les fichiers trimmés

## QC report par fastQC
parallel -X fastqc ${workDir}/trimmed/{} -o ${workDir}/trimmedReport/ -t $tFastQC ::: ${trimFile[@]} #voir si le -X marche #peut etre oter le -t

## QC global report par multiQC
multiqc -i 108_trimmed-RNAseq_Samples \
	-n allSamplesAnalysis \
	-o ${workDir}/trimmedReport/globalAnalysis \
	-p ${workDir}/trimmedReport
multiqc -i 108_trimmed-RNAseq_Samples \
	-n allSamplesInteractive \
	-o ${workDir}/trimmedReport/globalAnalysis \
	--interactive ${workDir}/trimmedReport


### Alignement des données
cd $workDir/mapping

## Construction de l'index 
#NB : on a au préalable concaténé BX571965.fna et BX571966.fna dans un fichier K96243.fna
bowtie2-build $workDir/ref_K96243/K96243.fna ref-K96243 \
	--threads $threadIndexBowtie \
	--seed $seed

## Alignement
parallel -k "bowtie2 --time \
	--very-sensitive \
	--threads $threadAlignBowtie \
	--seed $seed \
	--fr \
	-x ref-K96243 \
	-1 $workDir/trimmed/{}_R1_trimmed.fastq.gz \
	-2 $workDir/trimmed/{}_R2_trimmed.fastq.gz \
	-S {}.sam &> $workDir/mapReport/{}_time.txt" ::: ${nomBase[@]}

cd $workDir

#normalement inutile
#conda deactivate #on retombe sur (genomic_bact)

## QC alignement Bowtie2 par multiQC
multiqc -i 54PE_RNAseq_Samples \
	-b "108 fichiers de reads ont été trimmés avec cutadapt. Ce sont donc 54 couples de reads (PE) qui sont alignés avec Bowtie2." \
	-n allAligned \
	-o ${workDir}/mapReport \
	--interactive ${workDir}/mapReport

#newList
cd $workDir/mapping
nomBase=(`ls -S *.sam | cut -d . -f 1`) #actualisation
#cd $workDir

## Conversion SAM en BAM
parallel -k "samtools view {}.sam -o {}.bam" ::: ${nomBase[@]}

## Sorting des reads alignés dans les bam
parallel -k "samtools sort {}.bam -o {}.bam" ::: ${nomBase[@]}

## Indexage des BAM ordonnés
parallel -k "samtools index {}.bam" ::: ${nomBase[@]}

## Statistiques sur les BAM
parallel -k "samtools flagstat {}.bam > $workDir/mapReport/{}.flagstat.txt && \
	samtools idxstat {}.bam > $workDir/mapReport/{}.idxstat.txt && \
	samtools stats {}.bam > $workDir/mapReport/${}.stats.txt" ::: ${nomBase[@]}

## multiQC sur les stats
multiqc -i BAM_54*PE_RNAseq_Samples \
	-b "108 fichiers de reads ont été trimmés avec cutadapt. Ce sont donc 54 couples de reads (PE) qui sont alignés avec Bowtie2. Ils sont ensuite convertis, triés et analysés avec SAMtools" \
	-n sortedBAM \
	-o ${workDir}/mapReport \
	--interactive ${workDir}/mapReport


cd $workDir
parallel -k "picard MarkDuplicates I=./mapping/{}.bam \
	O=./markedDup/{}.dupAllmarked.bam \
	M=./markedDup/{}.dupAllmarkedMetrics.txt \
	MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1004 \
	TAGGING_POLICY=All \
	ASSUME_SORT_ORDER=coordinate" ::: ${nomBase[@]}

## Statistiques pour markedDup
cd $workDir/markedDup
parallel -k "samtools flagstat {}.dupAllmarked.bam > {}.dupMarked.flagstat.txt && \
	samtools stats {}.dupAllmarked.bam  > {}.dupMarked.stats.txt" ::: ${nomBase[@]}

# MultiQC pour markedDup flagstat
multiqc -i BAM_54*PE_RNAseq_Samples \
	-b "108 fichiers de reads ont été trimmés avec cutadapt. Ce sont donc 54 couples de reads (PE) qui sont alignés avec Bowtie2. Ils sont ensuite convertis, triés et analysés avec SAMtools. Enfin les duplicats sont taggués par MarkDuplicates" \
	-n markedBAM \
	-o ${workDir}/markedDup \
	--interactive ${workDir}/markedDup


### Observation du taux de duplication avec script R 
# i.e. : Faut-il faire une deduplication ?
cd $workDir
Rscript $workDir/dupRadar.R $workDir/markedDup $coreR #nombre de CPU
#NOTA : taux de duplication tolerable (modelisation controle OK)
# Donc on n'effectue pas de deduplication.


### Comptage des reads

cd $workDir/count
parallel -k -j27 "htseq-count -f bam \
	--order pos \
	--stranded reverse \
	--type gene \
	--idattr locus_tag \
	-m intersection-nonempty \
	{}.bam \
	${workDir}/BP_corr_nomChr.gff > {}.txt" ::: ${nomBase[@]}

# Fin du script
