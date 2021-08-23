#!/usr/bin/env python3
# Auteur : Ombeline LAMER
# But : Blastp des protéines contre protéines par paire d'organismes
# Commande : python3 blastp.py -org ../org.txt -faDir ../1_fastaProt -dbDir ../2_BlastDB -outDir ../3_BlastpOutput

import string, sys, os

#Menu d'aide
def usage():
    print ("""

     Obligatoire:
     ===========

	-org	--> liste des noms d'organisme dans un fichier texte
	-dbDir	--> répertoire contenant la database
	-faDir	--> répertoire des fichiers .fa de séquences protéiques
	-outDir	--> répertoire de sortie pour les résultats des BLASTp

  """)

#Interception des arguments en ligne de commande
try:
	org = sys.argv[sys.argv.index("-org")+1]
except:
	usage()
	sys.exit()

try:
	dbDir = sys.argv[sys.argv.index("-dbDir")+1]
except:
	usage()
	sys.exit()

try:
	faDir = sys.argv[sys.argv.index("-faDir")+1]
except:
	usage()
	sys.exit()

try:
	outDir = sys.argv[sys.argv.index("-outDir")+1]
except:
	usage()
	sys.exit()

#Récupération des noms de fichier à traiter (organismes listés)
f = open(org, "r")
lines = f.readlines()
f.close()

org = []
for line in lines :
	org.append(line.rstrip("\n")) #Oter \n en fin de ligne


f = open("../BLASTp_log.txt", "w") #Fichier de log des BLASTp
for o1 in org :
	print(o1+"vs. :")
	for o2 in org :
		print("\t"+o2)
		if (o1 == o2) :
            print("\t\t=> IGNORE")
		if (o1 != o2) :
			print("\t\t=> BLASTp")
			
			#Sauvegarde des paires d'organismes dont les protéines sont BLASTées
			f.write(o1 + "_" + o2 + "\n")
			
			#BLAST o1 vs o2 et sauvegarde des résultats dans le répertoire de sortie
			os.system("blastp -query {}/{}.fa -db {}/{} -out {}/{}_{} -outfmt '7 qseqid sseqid pident qlen length slen mismatch evalue bitscore gaps'".format(faDir, o1, dbDir, o2, outDir, o1, o2))

print("\n\tBLASTp tous effectués.\n")
f.close() #Fermeture du fichier en cours d'ecriture

#Fin script
