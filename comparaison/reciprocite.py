#!/usr/bin/env python3
# Auteur : Ombeline LAMER
# But : Création des listes de Bidirectionnal Best Hits BBH
# Commande : python3 reciprocite.py -org ../org.txt -inDir ../4_BestMatch -outDir ../5_BBH

import string, sys

#Menu d'aide
def usage():
    print ("""

     Obligatoire:
     ===========
	 -org		--> liste des noms d'organismes
	 -inDir		--> répertoire contenant les listes de best match par paire d'organismes
	 -outDir	--> répertoire de sorties pour les BBH
	 
  """)

#Interception des arguments en ligne de commande
try:
	inDir = sys.argv[sys.argv.index("-inDir")+1]
except:
	usage()
	sys.exit()
try:
	org = sys.argv[sys.argv.index("-org")+1]
except:
	usage()
	sys.exit()
try:
	outDir = sys.argv[sys.argv.index("-outDir")+1]
except:
	usage()
	sys.exit()

### Initialisation

f = open(org, "r")
organism = f.readlines() #contient les noms des organismes
f.close()

organismes=[] #stockage des noms d'organisme
for name in organism :
	organismes.append(name.rstrip('\n')) #oter le \n
#ok
allPairs=[]

#constituer les paires non redondantes (NR) de fichiers
for file1 in organismes :
	for file2 in organismes :
		#si les noms de fichiers sont differents
		#et que la paire n'est pas déjà notée dans un autre sens
		if file1!=file2 and [file2, file1] not in allPairs :
			#alors on stocke cette nouvelle paire
			allPairs.append([file1, file2])

del organism #destruction car inutile
del organismes

# enregistrement des paires NR dans un fichier (besoin pour clique.py)
f=open("../pairOrg.txt","w")
for pair in allPairs :
#NB : ne pas rassembler avec la boucle suivante sinon conflit d'ecriture
	f.write(pair[0]+"\t"+pair[1]+"\n")
f.close()

### Traitement

#Pour chaque paire de fichiers NR :
for pair in allPairs :
	
	## Recupération des Best Hits(BH)
	#entre les 2 organismes dans un sens
	f1 = open(inDir + "/" + pair[0] + "_" + pair[1], "r")
	BH1 = f1.readlines()
	f1.close()
	#et dans l'autre sens
	f2 = open(inDir + "/" + pair[1] + "_" + pair[0], "r")
	BH2 = f2.readlines()
	f2.close()
	
	#Initialisation des dictionnaires des BH
	org1={}
	org2={}
	#Ecriture des dictionnaires de BH
	for line in BH1 :
		#formatage
		line = line.rstrip("\n")
		line = line.split("\t")
		#entrée : key = prot from Org1, value = prot from Org2
		org1[line[0]]=line[1]
	for line in BH2 :
		line = line.rstrip("\n")
		line = line.split("\t")
		#entrée : key = prot from Org2, value = prot from Org1
		org2[line[0]]=line[1]
	
    #Liberation de la RAM
	del BH1
	del BH2
    
	## Constitution des Bidirectionnal BH (BBH)
    
	#Initialisation du fichier de sortie
	f=open(outDir + "/" + pair[0]+"_"+pair[1],"w")
	
	#Pour toutes les protéines de O1
	for protO1 in org1.keys():
		# pas besoin de chercher egalement pour toutes les prot de O2 
		# car si pas BH dans O1 ne pourra pas etre BBH (i.e. reciproque)
		
		#comme on a filtré les hits par un seuil, l'organisme reciproque peut ne pas avoir conservé de hit pour la proteine O2 matche par protO1 => c'est le premier test
		if org1[protO1] in org2.keys() :
			# 2e test : la protO2 pointée par le dico org1 pointe-t-elle sur protO1 ? i.e. retrouve-t-on protO1 ? => BBH => stockage
			if (org2[org1[protO1]]==protO1) :
				f.write(protO1+"\t"+org1[protO1]+"\n")
	f.close()

#Fin script
