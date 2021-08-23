# Transcriptomique

1. Script R utilisé pour les analyses exploratoires
2. DESeq2_explo_graph : dossier regroupant les sorties graphiques du script précédent
3. DE : dossier regroupant les sorties des analyses calculatoires suivant les différentes conditions considérées. Ce sont les listes des gènes (de l'isolat B1, utilisé pour l'alignement) dont l'expression varie le plus entre 2 situations données.
4. RNAseq_pipeline_K96243.sh : script de traitement des données raw, l'alignement et le comptage sur K96243.
5. RNAseq_pipeline_B1.sh : script de traitement des données trimmés (par le script précédent) avec l'alignement et le comptage des reads sur B1.
6. dupRadar.R : script R permettant de déterminer si les librairies RNAseq présentent un taux anormale de duplication des reads par ses sorties graphiques.
7. RNAseq.yaml : specification de l'environnement conda nécessaire aux bon déroulement des scripts utilisés pour les différentes analyses.
