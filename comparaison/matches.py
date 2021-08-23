#!/usr/bin/env python3
# Auteur : Ombeline LAMER
# But : Sélection du meilleur match protéine-protéine parmi les matchs significatifs des sorties de BLASTp
# Commande : python3 matches.py -paire ../BLASTp_log.txt -inDir ../3_BlastpOutput -outDir ../4_BestMatch

import string, sys

#Menu d'aide
def usage():
    print ("""

     Obligatoire:
     ===========

	 -paire	    --> liste des paires d organismes comparés en BLASTp (hors extension)
	 -inDir	    --> répertoire contenant les fichier de résultats des BLASTp
	 -outDir	--> répertoire de sortie listant les meilleurs matchs proteine-protéine par paire d'organisme

     Optionnel:
     ===========
     
    Seuils de significativité des matchs
	   -minID		--> pourcentage minimal d'AA identique entre la query et le subject
     -minLEN	--> valeur minimale d'alignement entre les deux protéines
     -minQCOV   --> pourcentage minimale de la couverture de la query par l'alignement
     -maxEval	--> valeur maximale d'evalue
     -minBITS	--> valeur minimale du bitscore
    
    Equation - définition d'un meilleur match :
        A changer au sein du script : lignes 165 & 174.
     
  """)

#Interception des arguments en ligne de commande

try:
	outDir = sys.argv[sys.argv.index("-outDir")+1]
except:
	usage()
	sys.exit()

try:
	inDir = sys.argv[sys.argv.index("-inDir")+1]
except:
	usage()
	sys.exit()

try:
	files = sys.argv[sys.argv.index("-paire")+1]
except:
	usage()
	sys.exit()

try:
	seuil_id = float(sys.argv[sys.argv.index("-minID")+1])
except:
	seuil_id = 80 #minimum de 80% de AA identique entre query et subject
try:
	seuil_align = float(sys.argv[sys.argv.index("-minLEN")+1])
except:
	seuil_align = 30 #minimum de 30 AA aligné entre query et subject
try:
	seuil_qcov = float(sys.argv[sys.argv.index("-minQCOV")+1])
except:
	seuil_qcov = 0.8 #minimum de 80% de la query aligné
try:
	seuil_evalue = float(sys.argv[sys.argv.index("-maxEval")+1])
except:
	seuil_evalue = 0.05 #maximum de 5% de probabilité que le match soit attribuable au hasard
try:
	seuil_bitscore = float(sys.argv[sys.argv.index("-minBITS")+1])
except:
	seuil_bitscore = 40 #minimum de score global (basé sur la similarité des séquences)

print("""

    RECAPITULATIF - Valeurs seuil utilisées pour la significativité des matchs :
    % min. AA identiques {}
    nb AA min. aligné {}
    % min. de la query alignée {}
    evalue max. {}
    bitscore min. {}
"""
.format(seuil_id,seuil_align,seuil_qcov,seuil_evalue,seuil_bitscore)
)

### Initialisation

f = open(files, "r")
lines = f.readlines() #lines contient tous les noms de fichiers de paires d'organismes de BLASTp
f.close()

blasts = [] #Liste permanente de tous les noms de fichiers
for line in lines :
	blasts.append(line.rstrip("\n")) #Oter le "\n" terminal

### Traitement des fichiers :

    #1) identification des matchs significatifs
    #2) identification des meilleurs matchs parmi les significatifs
    
for file in blasts : #pour tous les fichiers de BLASTp
	
	f = open(inDir+"/" + file, "r") #avec chemin vers le fichier d'intérêt
	lines = f.readlines() #contient toutes les infos du blastp-file en cours
	f.close()
	
	#dictionnaire temporaire : infos des BLASTp 1 key par protéineQuery
	dHits = {}

	for line in lines : #pour chaque ligne du fichier BLASTp
		
		if line[0:7] == "# Query" : #si on repère une protéineQuery
			
			line = line.rstrip("\n")
			query = line[9:] #récupération du nom de la protéine
			
			dHits[query] = {} #creation d'un dictionnaire par query
			dHits[query]["matchs"]=[] #initialisation de la liste des hits significatifs
		
        
		elif line[0] != "#" : #traitement des lignes informatives (hits) suivant la mention # Query
			line = line.rstrip("\n") #oter le "\n terminal. Utile si recuperation de la derniere valeur de la ligne
			line_split = line.split("\t") #separation des différents champs
			#TEST : print(line_split)
			#Attribution des valeurs de la ligne :
			#qseqid = line_split[0] #correspond à la valeur dans query
			subject = line_split[1] #nom de la proteine blastee contre la query
			pident = float(line_split[2]) #pourcentage d'AA identique
			#qlen = float(line_split[3]) #longueur query
			length = float(line_split[4]) #longueur de l'alignement
			#slen = float(line_split[5]) #longueur subject
			qcov = float(line_split[4])/float(line_split[3]) #pourcentage de couverture de la Query par l'alignement
			#mismatch = float(line_split[6]) #nb de mismatch
			evalue = float(line_split[7])
			bitscore = float(line_split[8]) #de formule inconnue (interne au BLAST)
			#gaps = float(line_split[9]) #nb de gap
			
			#Le hit est-il significatif ? => Si oui, c'est un match.
			if pident >= seuil_id and evalue < seuil_evalue and qcov >= seuil_qcov and length >= seuil_align :
				
				#Chaque match est retenu : KEY=nomDuMatch, VALUE=[%id,length,%qcov,eval,bitscore]
				dHits[query][subject]=[pident,length,qcov,evalue,bitscore]
				dHits[query]["matchs"].append(subject)
                
    #INFO : Toutes les lignes du fichier ont été lues !
    #Des matches de chaque query, on selectionne le meilleur
	for protQ in dHits.keys(): #pour chaque protéineQuery
		
        #si 0 match, la protéineQuery n'apparaitra pas dans la liste de sortie.
        #si 1 seul match => le meilleur par défaut
		if len(dHits[protQ]["matchs"]) == 1 and dHits[protQ]["matchs"] == [] :	
			dHits[protQ]["best"]=dHits[protQ]["matchs"][0]
		
		elif len(dHits[protQ]["matchs"]) > 1 : #Le meilleur match est choisi par comparaison de critère
			#print(dHits[protQ]["matchs"])
			best=dHits[protQ]["matchs"][0] #stockage du meilleur temporaire
			#stockage de ces valeurs [pident,length,qcov,evalue,bitscore]
			best_pident = dHits[protQ][best][0]
			best_length = dHits[protQ][best][1]
			best_qcov = dHits[protQ][best][2]
			best_evalue = dHits[protQ][best][3]
			best_bitscore = dHits[protQ][best][4]
			best_formulaScore = ( (best_pident/100)+best_qcov+(1-best_evalue) ) * best_length
            
			for match in dHits[protQ]["matchs"][1:] : #pour tous les matchs de la protéineQuery excepté le premier best temporaire
				
				match_pident = dHits[protQ][match][0]
				match_length = dHits[protQ][match][1]
				match_qcov = dHits[protQ][match][2]
				match_evalue = dHits[protQ][match][3]
				match_bitscore = dHits[protQ][match][4]
				match_formulaScore = ( (match_pident/100)+match_qcov+(1-match_evalue) ) * match_length
                
				#Si tous les critères sont les meilleurs : alors le match est meilleur (absolu)
				if match_pident >= best_pident and match_length >= best_length and match_qcov >= best_qcov and match_evalue <= best_evalue and match_bitscore >= best_bitscore :
					best_pident = match_pident
					best_length = match_length
					best_qcov = match_qcov
					best_evalue = match_evalue
					best_bitscore = match_bitscore
					best_formulaScore = match_formulaScore
					best = match
				
                #Que faire si l'un des critères est meilleur et un autre pire ?
                #On définit un meilleur match de manière relative par une équation.
				#=> on résume arbitrairement les critères retenus en un score par la formulaScore				
				elif match_formulaScore >= best_formulaScore :
					best_pident = match_pident
					best_length = match_length
					best_qcov = match_qcov
					best_evalue = match_evalue
					best_bitscore = match_bitscore
					best_formulaScore = match_formulaScore
					best = match
			
            #Fin Boucle multiMatch : meilleur définitif stocké dans best
			dHits[protQ]["best"]=best
			
	#Sauvegarde du meilleur match par protéine
	f = open(outDir + "/" + file, "w")
	for protQ in dHits.keys():
		if "best" in dHits[protQ].keys() and dHits[protQ]["best"]!='':
			f.write(protQ+"\t"+dHits[protQ]["best"]+"\n")
	f.close
#End script
