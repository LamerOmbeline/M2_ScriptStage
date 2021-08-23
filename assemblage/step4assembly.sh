#######################################################################
### Obtenir les séquences génomiques des isolats par assemblage hybride
### Obtenir les annotations structurelles des isolats
#######################################################################

### Commandes Unicycler

## Obtenir des assemblages "brut" - env. 6h par assemblage
## /!\ Utiliser les scripts java et d'appel de Pilon modifié dans le module de conda où se trouve Unicycler 0.4.8

# isolat A2
unicycler --threads 48 \
  --mode conservative \
  -1 A2_Illumina_reads_1.fastq.gz -2 A2_Illumina_reads_2.fastq.gz \
  -l A2_MinION_reads.fastq.gz \
  -o UnicyclerAssembly_A2

# isolat B1
unicycler --threads 48 \
  --mode conservative \
  -1 B1_Illumina_reads_1.fastq.gz -2 B1_Illumina_reads_2.fastq.gz \
  -l B1_MinION_reads.fastq.gz \
  -o UnicyclerAssembly_B1

### Commandes QUAST - etudier les recombinaisons

## Comparaison globale entre la reference K96243, A2 et B1 (REFvsA2vsB1)
quast A2_assembly.fasta B1_assembly.fasta _
  -r K96243.fna \
  --feature BP.gff \
  --threads 45 \
  --labels "A2, B1" \
  -o  REFvsA2vsB1 \
  --circos \
  --split-scaffolds

## Construction des upperbounds des isolats sur A2 et B1
#on donne les reads courts et long pour mieux comprendre les alignements et missassembly

# isolat A2
quast A2_assembly.fasta \
  -r K96243.fna \
  --feature gene:BP.gff \
  --threads 45 \
  --labels "A2" \
  -o  REFvsA2_upperBoundAssembly \
  --circos \
  --split-scaffolds \
  --upper-bound-assembly \
  --pe1 A2_Illumina_R1.fastq.gz \
  --pe2 A2_Illumina_R2.fastq.gz \
  --nanopore A2_MinION.fastq.gz \
  --glimmer

# isolat B1
quast B1_assembly.fasta \
  -r K96243.fna \
  --feature gene:BP.gff \
  --threads 45 \
  --labels "B1" \
  -o  REFvsB1_upperBoundAssembly \
  --circos \
  --split-scaffolds \
  --upper-bound-assembly \
  --pe1 B1_Illumina_R1.fastq.gz \
  --pe2 B1_Illumina_R2.fastq.gz \
  --nanopore B1_MinION.fastq.gz \
  --glimmer

# NB : QUAST et bedtools ne sont pas compatibles sur le format des séquences, l'un considère la première base à 1 l'autre à 0.

## Comparaison A2 et B1 à la référence K96243 
quast A2_curated_assembly.fasta B1_curated_assembly.fasta _
  -r K96243.fna \
  --feature BP.gff \
  --threads 45 \
  --labels "A2, B1" \
  -o  REFvsA2vsB1 \
  --circos \
  --glimmer

### Commandes Prokka - annotation structurelle des assemblages des isolats

# construction du fichier d'entrainements pour Prokka par Prodigal v2.6.3 (Feb 2016)
prodigal -i K96243.fna -t tf_K96243

# isolat A2
### pour A2_full
prokka --outdir ./A2_full_clean -\
  -prefix A2 \
  --addgenes --addmrna \
  --locustag BPA2 \
  --gffver 3 \
  --compliant \
  --genus Burkholderia --species pseudomallei \
  --strain A2 \
  --prodigaltf tf_K96243  \
  --proteins K96243.gb \
  --cpus 0 \
  --rfam A2_curated_assembly.fa

# isolat B1
prokka --outdir ./B1_full_clean \
  --prefix B1 \
  --addgenes --addmrna \
  --locustag BPB1 \
  --gffver 3 \
  --compliant \
  --genus Burkholderia --species pseudomallei \
  --strain B1 \
  --prodigaltf tf_K96243  \
  --proteins K96243.gb \
  --cpus 0 \
  --rfam B1_curated_full.fa
