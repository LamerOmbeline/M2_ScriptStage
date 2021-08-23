#!/usr/bin/env python3
# Auteur : Ombeline LAMER
# But : Création de la base de données : protéines par organisme
# Commande : python3 blastdb.py -org ../org.txt -inDir ../1_fastaProt -outDir ../2_BlastDB

import string, sys, os

#Menu d'aide
def usage():
    print ("""

     Obligatoire:
     ===========
     
	-org	--> liste des noms d'organisme dans un fichier texte sans extension
	-inDir	--> répertoire contenant les fichiers de séquences protéiques au format fasta .fa explicite
	-outDir	--> répertoire de dépôt pour la base de données à construire

  """)

#Interception des arguments en ligne de commande
try:
	org = sys.argv[sys.argv.index("-org")+1]
except:
	usage()
	sys.exit()

try:
	inDir = sys.argv[sys.argv.index("-inDir")+1]
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

#Construction de la database de protéines
for o in org :
	os.system("makeblastdb -in {}/{}.fa -dbtype \"prot\" -out {}/{}".format(inDir, o, outDir, o))

#Fin script
