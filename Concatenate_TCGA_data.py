import gzip
import os
import numpy as np

os.chdir('DATA/')
#%%
###### GET FILE'S PATH IN THE MANIFEST #######

#MANIFEST contient le nom de tout les fichier, avec le nom de chemin relatif
manifest_brest = open('TGCA_breast_control/MANIFEST.txt')
list_paths_brest = ['TGCA_breast_control/'+line.split('\t')[1] for line in manifest_brest][1:] #Récupération de tous les paths des fichiers brest


#%%
####### OPEN AND PUT IN LIST ALL THE FILES #######

def get_files_from_paths(list_paths):
#Renvoie la liste des fichier ouvert associé à la liste des paths. annotation pour le moment pas utilisé mais quand meme stoqué si besoin.
    files_list = []
    annotation_files = []
    files_missing = 0
    for path in list_paths: 
        try:
#            if path[-11:-7] == 'FPKM':    #For FPKM normalized expression lvl
            if path[-9:-3] == 'counts':    #Il semblerais que le clustering est plus pertinants sur les raw counts mais le choix est libre
                file = gzip.open(path, 'r')
                files_list.append(file)
            if path[-15:-4] == 'annotations':
                file = open(path, 'r')
                annotation_files.append(file)
        except:
            files_missing +=1
            
    print('missing', files_missing, '/', files_missing+len(files_list))
    return files_list, annotation_files

files_list_brest, annotation_liste_brest = get_files_from_paths(list_paths_brest)

#%%   
####### GET DATA IN FILES (LONG TIME) #######

def get_data_from_RNAseq_files(files_list, name):
#Récuparation sous la forme d'un tableau de tout les raw counts de chacun des fichiers dans "files_list"
    list_genes = []
    expression_tot = []
    expression = []
    
    for line in files_list[0]:
        line = line.decode("utf-8")
        list_genes.append(line.split('\t')[0]) 
        expression.append(float(line.split('\t')[1]))
    
    i = 0
    for file in files_list:
        expression_tot.append(expression)
        expression = []
        try:
            for line in file:
                line = line.decode("utf-8")
                expression.append(float(line.split('\t')[1]))
        except:
                print('bug:')
        i+=1
        print("", end=f"\r{name} {i} / {len(files_list)} complete")
        
    del expression_tot[1] #Correspond à un élément vide car le premier fichier de la liste à déja été ajouté pour comparé les gènes avec les autres fichiers.
    print("")
    return list_genes, expression_tot

genes_brest, expression_brest = get_data_from_RNAseq_files(files_list_brest, 'brest:')

#%%

import pandas as pd
expression_brest = pd.DataFrame(expression_brest)

#%%
# CONVERSION STABLE ID TO SYMBOL via un dictionnaire(plus facile à lire)

dico_conversion = dict()
file = open("BiomartData/human_symbol_ID.txt")
dup = 0

for line in file:
    line = line.strip("\n")
    if line.split(",")[0] not in dico_conversion:
        dico_conversion[line.split(",")[0]] = line.split(",")[1]
    else:
        dico_conversion[line.split(",")[0]] = [dico_conversion[line.split(",")[1]], line.split(",")[1]]
        dup +=1
        
file.close()

#%%
#Conversion des genes ID -> STABLE
converted_genes_names = list()
nop = 0
for gene in genes_brest:
    gene = gene.split(".")[0]
    try:
        converted_genes_names.append(dico_conversion[gene])
    except:
        converted_genes_names.append(gene)
        print(gene)    
        nop +=1
        
print(nop)

#%%
   
expression_brest.columns = converted_genes_names

#%%
expression_brest.T
export_csv = expression_brest.to_csv ('breast_dataframe_OMICS.csv', index = True) 

#%%

