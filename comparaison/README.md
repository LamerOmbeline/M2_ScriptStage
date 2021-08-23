# Comparaison de génomes

Vous trouverez ici les scripts python ayant servi à comparer les génomes des isolats et K96243 comme présenté dans le rapport associé. L'ordre d'explication des scripts correspond à l'ordre dans lequels ils doivent être appelés. Des exemples d'appel et la définition des options sont donnés au sein des scripts.


1. blastdb.py : construction des bases de données protéiques.
2. blastp.py : BLAST protéines-protéines d'un organisme sur un autre, pour tous les organismes considérés
3. matches.py : Identification du meilleur hit significatif d'une protéine contre une autre de 2 organismes différents pour tous les organismes considérés.
4. reciprocite.py : Identification des BBH (Bi-directionnal Best Hit) pour une paire d'organisme, pour tous les organismes considérés. Arrêt ici si seulement 2 organismes considérés au total.
5. clique.py : Recherche des cliques (3 organismes) pour constituer le core génome et les génomes accessoires des organismes. Visualisation avec UpsetPlot. Les diagrammes de Venn sont constitués manuellement à partir des listes écrites en sorties de ce script.
6. Upset_XXX.png : Exemples de sortie graphique de clique.py. Permet une visualisation 2D synthétique sans limite du nombre d'organismes considérés (contrairement aux diagrammes de Venn dont la limite est 5 organismes). minID = pourcentage minimum d'AA identiques entre les 2 protéines alignées. eval = evalue maximale. NR = non restrictif. R = Restrictif.   
