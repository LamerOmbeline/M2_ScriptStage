#!/usr/bin/env python3
# Auteur : Ombeline LAMER
# But : Constitution des cliques complètes (core genome) et genomes accessoires
# Commande : python3 clique.py -files ../org.txt ../pairOrg.txt -faDir ../1_fastaProt -inDir ../5_BBH -outDir ../6_cliques -outProt ../7_proteines -workDir ../

### ATTENTION ! A DEBUGGUER ==> mauvais comptage de protéines pour Upset !

import string, sys, os, networkx as nx
#installer networkx dans le yaml conda!

#Menu d'aide
def usage():
    print ("""

     obligatoire :
     ===========
     
     -faDir		--> répertoire contenant les fichiers multi-fasta de protéines
	 -files		--> 2 fichiers attendus dans cet ordre (1) liste de noms des organismes (2) liste des paires d'organisme
     -outProt	--> repertoire allant contenir les fichiers listant les protéines par organisme
	 -outDir	--> répertoire contenant les cliques
	 -inDir		--> répertoire contenant les BBH formattés pour networkx (pour chaque paire d'organismes, non redondant)
	 -workDir	--> répertoire de travail global
	 
  """)

#Interception des arguments en ligne de commande
try:
	faDir = sys.argv[sys.argv.index("-faDir")+1]
except:
	usage()
	sys.exit()

try:
	nameOrg = sys.argv[sys.argv.index("-files")+1]
	pairOrg = sys.argv[sys.argv.index("-files")+2]
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

try:
	outProt = sys.argv[sys.argv.index("-outProt")+1]
except:
	usage()
	sys.exit()
try:
	workDir = sys.argv[sys.argv.index("-workDir")+1]
except:
	usage()
	sys.exit()

### Initialiation

#Construction des listes des noms d'organismes
f=open(nameOrg,"r")
nomOrg=f.readlines()
f.close()

listOrg=[] #liste de noms des organismes
for nom in nomOrg :
	listOrg.append(nom.rstrip('\n')) #oter le saut de ligne
del nomOrg

#dico des noms de protéine par organisme KEY=Org : VALUE:proteinName
dProtParOrg={}
# pour tous les fichiers de prot, recuperer les noms de protéines
for organisme in listOrg :
	#creer le fichier en dur
	os.system("grep '>' {}/{}.fa | sed 's/>//g' > {}/{} ".format(faDir, organisme, outProt, organisme))
	
	# charger la liste protéines par organismes
	dProtParOrg[organisme]={ "allProt" : [], "protIsolees" : [] }
	
	f = open(outProt+"/"+organisme, "r")
	lines = f.readlines()
	f.close()
	
	for prot in lines :
		dProtParOrg[organisme]["allProt"].append(prot.rstrip("\n"))
del lines

## init listes :
# 1) dico pour les protéines isolees par organismes
#=> Fait plus haut, dans le dictionnaire avec protIsolees
# 2) liste globales des protéines traitees :
classifiedProt=[]

## initialisation graph
G = nx.Graph()

# ecriture des nodes avec les listes de protéines
for org in dProtParOrg.keys() :
	G.add_nodes_from(dProtParOrg[org]["allProt"], organisme=org)

# recuperation donnees : BBH

#nom des paires d'organisme non redondant (NR)
f=open(pairOrg,"r")
paireOrg=f.readlines()
f.close()

pairBBH=[] #liste des paires de noms d organisme
for line in paireOrg :
	couple=line.rstrip('\n') #oter le saut de ligne
	couple=couple.split('\t')
	pairBBH.append([couple[0],couple[1]])
del(paireOrg,pairOrg)

BBH4Graph=[]
for couple in pairBBH :
#pour tous les fichiers de BBH
	BBHlist=str(couple[0]+'_'+couple[1])
	#lire le fichier
	f=open(inDir+'/'+BBHlist,"r")
	BBHlines=f.readlines()
	f.close()
	
	for BBH in BBHlines :
		edge=BBH.rstrip("\n")
		edge=edge.split("\t")
		#NB : 1ligne = 1tuple =>format des edges
		BBH4Graph.append(tuple(edge))

# ecriture des arêtes avec les BBH
G.add_edges_from(BBH4Graph)
del(BBH4Graph)#degagement RAM


## constituer les listes non redondantes de toutes les cliques possibles
		#NB : aide au stokage pour les comptages/Venn
#automatiser après pour plus de cliques...
#TO DO upgrade : separer la mise en forme sortie. automatiser pour plus de 3 organismes.

print("Initialisation terminée. Graph créée.")
print("Lancement Recherche cliques completes. Constitution G. core et accessoires")

coreGenome=[] #liste des cliques completes (taille 3)

dPaire={'KB' : {}, 'KA' : {}, 'BA' : {} } #cliques de taille 2
dPaire['KB']={'org' : ['K96243','B1'], 'clique' : []}
dPaire['KA']={'org' : ['K96243','A2'], 'clique' : []}
dPaire['BA']={'org' : ['B1','A2'], 'clique' : []}

paire=[['K96243','B1'],['K96243','A2'],['B1','A2']] #ordre des types de cliques

#pour les cliques de taille 1 :
	#seront rangés dans dProtParOrg[ORG]['protIsolees']=[...]


### Traitement ### 


allCliques=list(nx.find_cliques(G))

###svg temporaire au cas ou :
#f=open(outDir+'allCliquesList','w')
#for clique in allCliques :
	#f.write(str(clique)+'\n')
#f.close()

for clique in allCliques :
	
	if len(clique) == 3 : #clique complete !
		coreGenome.append(clique)
		#pas besoin de stockes dans classifiedProt parce que find_cliques les exclues d'autre cliques commme 1) BBH et cas de 3 org
		
	elif len(clique) == 1 : #proteine isolee
		organisme=G.nodes[clique[0]]['organisme']
		dProtParOrg[organisme]['protIsolees'].append(clique[0])
		#pas la peine de mettre dans classifiedProt parce que exclue de ailleurs par def des BBH sinon.
		
	else : #couple len(clique==2)
		prot1=clique[0] #nom des prot
		prot2=clique[1]
		o1=G.nodes[clique[0]]['organisme'] #nom des organismes des prots
		o2=G.nodes[clique[1]]['organisme']
		
		#quelle est la protéine problematique pour la stringence ?
		trouble=0
		if G.degree(prot1) >= 2 :
			trouble=1
			classifiedProt.append(prot1)
			
			#dans ce cas, est à exclure :
			if prot1 not in dProtParOrg[o1]['protIsolees'] :
				dProtParOrg[o1]['protIsolees'].append(prot1)
		
		elif G.degree(prot2) >= 2 :
			trouble=2
			classifiedProt.append(prot2)
			
			if prot2 not in dProtParOrg[o2]['protIsolees'] :
				dProtParOrg[o2]['protIsolees'].append(prot2)
			
		elif G.degree(prot1) ==2 and G.degree(prot2) == 2 :
			trouble=3
			classifiedProt.append(prot1)
			classifiedProt.append(prot2)
			
			if prot1 not in dProtParOrg[o1]['protIsolees'] :
				dProtParOrg[o1]['protIsolees'].append(prot1)
			if prot2 not in dProtParOrg[o2]['protIsolees'] :
				dProtParOrg[o2]['protIsolees'].append(prot2)
		
		else : # G.degree(prot1) == 1 and G.degree(prot2) == 1 :
			#aucun probleme c'est bien une paire !
			
			#maintenant il va falloir la stocker au bon endroit...
			coupleOO1=[o1,o2]
			coupleOO2=[o2,o1]
		
			#retrouver le bon endroit
			initialPair=list(dPaire.keys())
			
			notmarked=0
			for ip in initialPair :
				really=False
				if coupleOO1 == dPaire[ip]['org'] :
					orgs=dPaire[ip]['org'][0][0]+dPaire[ip]['org'][1][0] #selection 1e lettre
					dPaire[orgs]['clique'].append(clique) # ATTENTION pas dans l'ordre
				elif coupleOO2 == dPaire[ip]['org'] :
					orgs=dPaire[ip]['org'][0][0]+dPaire[ip]['org'][1][0]
					dPaire[orgs]['clique'].append(clique)
				else :
					really=True
			if really==False :
				notmarked+=1

#EndFor, sur les cliques identifiées

### Stockage des cliques !!


### Sorties
#listes cliques completes (nb de lignes = nb de cliques)
f=open(outDir+'/CoreGenome_K96-B1-A2','w')
core=0
for clique in coreGenome :
	f.write(clique[0]+'\t'+clique[1]+'\t'+clique[2]+'\n')
	core+=1
f.close()

#listes cliques pseudo-completes (1 par ensemble d'organismes différents)
dCount=[ 0, 0, 0 ]
	# KB KA BA
	# NS en 1er dans la liste, S ensuite

### ATTENTION constat les organisme sont stocké au lieu des protéines en clique !
for cliqueT2 in dPaire.keys() : #verifier que les listes-clés sont autorisées
	f=open(outDir+'/'+dPaire[cliqueT2]['org'][0]+'_'+
		dPaire[cliqueT2]['org'][1],"w")
	
	for clique in dPaire[cliqueT2]['clique'] :
		f.write(clique[0]+clique[1]+'\n')
		
		if cliqueT2 == 'KB' :
			dCount[0]+=1
		elif cliqueT2 == 'KA' :
			dCount[1]+=1
		elif cliqueT2 == 'BA' :
			dCount[2]+=1
		
	f.close()


#listes prot isolees (une par organisme)
A=0;B=0;K=0

f1=open(outDir+'/K96243','w')
for prot in dProtParOrg['K96243']['protIsolees'] :
	f1.write(prot+'\n')
	K+=1
f1.close()

f2=open(outDir+'/B1','w')
for prot in dProtParOrg['B1']['protIsolees'] :
	f2.write(prot+'\n')
	B+=1
f2.close()

f3=open(outDir+'/A2','w')
for prot in dProtParOrg['A2']['protIsolees'] :
	f3.write(prot+'\n')
	A+=1
f3.close()

#ecriture du diagramme de venn => "a la main" /OU/ sortir un fichier avec les comptes :
    #si le nombre est compris entre 2 et 5 organismes
    #1e colonne liste les noms d organismes impliques
    #2e colonne donne le nb de protéine concernees

#ecriture du doc avec le package UpSetPlot (v 0.6.0 channel conda-forge)
#appel d'un script myUpset.py
import upsetplot

mesDataList=[#[], # TEST sans ensemble vide car inutile en situation biologique
			 ['B1'],
			 ['A2'],
			 ['A2','B1'],
			 ['K96243'],
			 ['K96243','B1'],
			 ['K96243','A2'],
			 ['K96243','A2','B1']
			]
mesCounts=[B,A,dCount[2],K,dCount[0],dCount[1],core]

#formattage upset
unsaneData=upsetplot.from_memberships(mesDataList,data=mesCounts)

#impression des plots
import matplotlib.pyplot as plt
imageNS = upsetplot.UpSet(unsaneData, subset_size='sum', show_counts=True).plot()
plt.savefig(outDir+'/UpSet.png')
print("Terminee")
#Fin Script
