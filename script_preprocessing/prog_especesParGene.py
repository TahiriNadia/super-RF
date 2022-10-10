
import os
from sys import path_hooks

proteines = [	
                'Mevalonate kinase',
                'Chorismate synthase',
                'Phosphoribosyl-AMP cyclohydrolase',
                'Adenylosuccinate lyase',
                'Argininosuccinate lyase',
                'Argininosuccinate synthase',
                'Histidinol dehydrogenase',
                'Glutamyl-tRNA reductase',
                'Diphthine synthase',
                '3-phosphoshikimate 1-carboxyvinyltransferase',
                'Prephenate dehydratase',
                'Phosphoribosylformylglycinamidine cyclo-ligase',
                '6,7-dimethyl-8-ribityllumazine synthase',
                'CTP synthase',
                'Aspartate carbamoyltransferase',
                'Diaminopimelate epimerase',
                'Valyl-tRNA synthetase',
                'Aspartate-semialdehyde dehydrogenase',
                'Acetylglutamate kinase',
                'Tyrosyl-tRNA synthetase',
                'Signal recognition particle-docking protein FtsY',
                'Leucyl-tRNA synthetase',
                'Arginyl-tRNA synthetase',
                'Adenylosuccinate synthetase',
                'UbiD family decarboxylase',
                'Glutamyl-tRNA synthetase',
                'Methionyl-tRNA synthetase',
                'Aspartyl-tRNA synthetase',
                'Tryptophanyl-tRNA synthetase',
                'Ribonuclease Z',
                "2'-5' RNA ligase",
                'Prolyl-tRNA synthetase',
                'Porphobilinogen synthase',
                'Glycyl-tRNA synthetase',
                'Protein-L-isoaspartate O-methyltransferase'
            ]

nomFichierEspeces = os.listdir('PATH' + '/genesParEspece')
nomFichierProteines = os.listdir('PATH' + 'genes/1_raw')

#Fonction permettant de liste l'ensemble des protéines présentes chez les différentes espcèces
def extraireLesProteines():
    
    print("Début de l'exctraction")
    setProteines = {""}
    
    #Dresser la liste de protéines présentes dans l'ensemble des génomes
    for nomFichierEspece in nomFichierEspeces:
        if (nomFichierEspece != ".DS_Store"):
            #Ouverture et lecture du fichier contenant les gènes des espèces
            fichierEspece = open('PATH/' + nomFichierEspece)
            lignes = fichierEspece.readlines()

            #Lecture de chacune des ligne du fichier
            for ligne in lignes:
                if ligne[0] == '>':
                    #Déterminer le nom de la protéine
                    idxDepart = int(ligne.find("[protein=") + 9)
                    idxFin = int(ligne.find("]", idxDepart))
                    nomProteine = ligne[idxDepart:idxFin]
                    nomProteine = nomProteine.replace('/', '#')

                    #Remplacer la première lettre minuscule par une majuscule
                    if ord(nomProteine[0]) >= 97 and ord(nomProteine[0]) <= 122:
                        nouvLettre = chr(ord(nomProteine[0]) - 32)
                        proteine = nouvLettre
                        idxLettre = 1
                        while idxLettre < len(nomProteine):
                            proteine += nomProteine[idxLettre]
                            idxLettre+=1
                    else:
                        proteine = nomProteine
                    

                    #Ajout de la proterine au répertoire (set)
                    setProteines.add(proteine)
        fichierEspece.close()

    #Noter le résultat
    listeProteines = open([PATH], "w")
    listeProteines.write("proteines = [")
    for proteine in setProteines:
        listeProteines.write('\t"' + proteine + '",\n')
    listeProteines.write("]")
    listeProteines.close()
    return 0

#Cette fonction permet d'écrire un fichier .fasta par protéine, où la séquence de cette protéine
#pour chacune des espèces qui en portent le gène est indiquée. 
#Accessoirement, la fonction permet de calculer la fréquence d'une protéine chez les différentes espèces 
#en plus de détecter si la protéine est présente plus d'une fois chez dans un même génome.
def ecrireFichiersProteines():
    print("Début du calcul de la fréquence")
    
    #Variable pour le suivi du fonctionnement du programme
    stats = []

    #Vérification de la présence de la protéine d'intérêt dans le génome des espèces, pour toutes les protéines
    for proteine in proteines:
        print(proteine)
        #Création du fichier pour noter les génomes
        nomFichierGenes = proteine.replace(" ", "_") + ".fasta"
        fichierGenes = open('PATH' + nomFichierGenes, "w")

        #Ajuster les caractères problématiques 
        proteine = proteine.replace('#', '/')

        #Déterminer le nom de la protéine avec une lettre majuscule et sans lettre majuscule
        #À noter que dans le répertoire, elles ont toujours une lettre majuscule
        #Il faut donc vérifier si la protéine commence par une lettre majuscule, pour trouver 
        #la version sans lettre majuscule également. 
        if (ord(proteine[0]) >= 65 and ord(proteine[0]) <=90):
            appellationProMaj = "[protein=" + proteine  + "]"
            nouvLettre = chr(ord(proteine[0]) + 32)
            protMin = nouvLettre
            idxLettre = 1
            while idxLettre < len(proteine):
                protMin += proteine[idxLettre]
                idxLettre+=1
            appellationProtMin = "[protein=" + protMin  + "]"
        #Si elle commence pas autre chose qu'une lettre, les deux appellations sont identiques
        else:
            appellationProtMin = appellationProMaj = "[protein=" + proteine  + "]"

        

        #Vérification de la présence de la protéine dans le génone de toutes les espèces
        frequence = 0
        doublonDeLaProteineChezUneEspece = False
        for nomFichierEspece in nomFichierEspeces:
            proteinePresenteChezLespece = False
            if (nomFichierEspece != ".DS_Store"):
                #Ouverture et lecture du fichier contenant les gènes des espèces
                fichierEspece = open('PATH/' + nomFichierEspece)
                lignes = fichierEspece.readlines()

                #Lecture de chacune des ligne du fichier
                for idxLigne in range(len(lignes)):
                    #Vérifier si la proétine d'intérêt se trouve sur la ligne.
                    #Si la protéine se trouve sur la liste, elle est ajoutée au fichiers en sortie
                    if appellationProtMin in lignes[idxLigne] or appellationProMaj in lignes[idxLigne]:
                        #Mise à jour de stats
                        frequence += 1
                        if proteinePresenteChezLespece:
                            doublonDeLaProteineChezUneEspece = True
                        proteinePresenteChezLespece = True
                        
                        #Nom de l'espèce
                        nomEspece = nomFichierEspece[0:-4]
                        fichierGenes.write(">" + nomEspece + "\n")
                        #Gène de la protéine
                        idxLigne += 1
                        while lignes[idxLigne][0] != '>' and idxLigne < len(lignes) - 1:
                            fichierGenes.write(lignes[idxLigne])
                            idxLigne+=1

        #Update des stats
        stats.append((proteine, frequence, doublonDeLaProteineChezUneEspece))

        #Fermeture du fichier des gènes d'une espèce    
        fichierGenes.close()

    #Noter des statistiques
    statistiques = open('PATH' + stats.txt', "w")
    statistiques.write(str("Fréquence de chacune des protéines présentes dans le génome de l'ensemble des 43 organismes : \n"))
    for stat in stats:
        statistiques.write("%s\tfrequence: \t%d\tdoublon: \t%i\n" % stat)
    statistiques.close()

    print("Stats : ")
    for stat in stats:
        print (" - %s \t frequence: %d, doublon: %i" % stat)
    
    return 0

#Cette fonction permet de calculer la fréquence d'une protéine chez les différentes espèces 
#en plus de détecter si la protéine est présente plus d'une fois chez dans un même génome.
def calculerFrequenceProteine():
    print("Début du calcul de la fréquence")
    
    #Variable pour le suivi du fonctionnement du programme
    stats = []

    #Vérification de la présence de la protéine d'intérêt dans le génome des espèces, pour toutes les protéines
    for proteine in proteines:
        print(proteine)

        #Ajuster les caractères problématiques 
        proteine = proteine.replace('#', '/')
        
        #Déterminer le nom de la protéine avec une lettre majuscule et sans lettre majuscule
        #À noter que dans le répertoire, elles ont toujours une lettre majuscule
        #Il faut donc vérifier si la protéine commence par une lettre majuscule, pour trouver 
        #la version sans lettre majuscule également. 
        if (ord(proteine[0]) >= 65 and ord(proteine[0]) <=90):
            appellationProMaj = "[protein=" + proteine  + "]"
            nouvLettre = chr(ord(proteine[0]) + 32)
            protMin = nouvLettre
            idxLettre = 1
            while idxLettre < len(proteine):
                protMin += proteine[idxLettre]
                idxLettre+=1
            appellationProtMin = "[protein=" + protMin  + "]"
        #Si elle commence pas autre chose qu'une lettre, les deux appellations sont identiques
        else:
            appellationProtMin = appellationProMaj = "[protein=" + proteine  + "]"

        #Vérification de la présence de la protéine dans le génone de toutes les espèces
        frequence = 0
        doublonDeLaProteineChezUneEspece = False
        for nomFichierEspece in nomFichierEspeces:
            proteinePresenteChezLespece = False
            if (nomFichierEspece != ".DS_Store"):
                #Ouverture et lecture du fichier contenant les gènes des espèces
                fichierEspece = open('PATH' + '/genesParEspece/' + nomFichierEspece)
                lignes = fichierEspece.readlines()

                #Lecture de chacune des ligne du fichier
                for idxLigne in range(len(lignes)):
                    #Vérifier si la proétine d'intérêt se trouve sur la ligne.
                    #Si la protéine se trouve sur la liste, elle est ajoutée au fichiers en sortie
                    if appellationProtMin in lignes[idxLigne] or appellationProMaj in lignes[idxLigne]:
                        #Mise à jour de stats
                        frequence += 1
                        if proteinePresenteChezLespece:
                            doublonDeLaProteineChezUneEspece = True
                        proteinePresenteChezLespece = True
                        

        #Update des stats
        stats.append((proteine, frequence, doublonDeLaProteineChezUneEspece))

    #Noter des statistiques
    statistiques = open('PATH' + 'gestionGenes/stats.txt', "w")
    statistiques.write(str("Fréquence de chacune des protéines présentes dans le génome de l'ensemble des 43 organismes : \n"))
    for stat in stats:
        statistiques.write("%s\tfrequence: \t%d\tdoublon: \t%i\n" % stat)
    statistiques.close()
    
    return 0

#Fonction permettant de lister les espèces présentes dans les fichiers de protéines
def listerEspeces():
    
    setEspeces = {""}
    
    #Dresser la liste de protéines présentes dans l'ensemble des génomes
    for nomFichierProteine in nomFichierProteines:
        if (nomFichierProteine != ".DS_Store"):
            #Ouverture et lecture du fichier contenant la séquence de gène chez chacune des espèces
            fichierEspece = open('PATH' + 'genes/1_raw/' + nomFichierProteine)
            lignes = fichierEspece.readlines()

            #Lecture de chacune des ligne du fichier
            for ligne in lignes:
                if ligne[0] == '>':
                    #Déterminer le nom de l'espèce
                    idxDepart = 1
                    idxFin = len(ligne) - 1
                    nomEspece = ligne[idxDepart:idxFin]

                    #Ajout de la proterine au répertoire (set)
                    setEspeces.add(nomEspece)
        fichierEspece.close()

    #Afficher en sortie le résultat
    for espece in setEspeces:
        print(espece)

    return 0

def __main__():

    #extraireLesProteines()
    #calculerFrequenceProteine()
    ecrireFichiersProteines()
    #listerEspeces()

    return 0

__main__()
