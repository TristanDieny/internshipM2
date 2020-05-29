# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 10:57:12 2020

@author: tristan
"""
import pandas as pd 
from scipy import stats
import numpy as np
import functions as fc
#%%

def mice_to_human(liste_names): #Fonction pour convertir les noms des genes de souris en humains (pour les orthologues..)
    file = open('mouse_to_human.txt')
    
    translateMouse_symbol = dict()
    for line in file:
        line = line.split('\t')
        translateMouse_symbol[line[1]] = [line[2], line[3].strip('\n')]
        
    converted = []
    missing = []
    keeped = []
    for gene in liste_names:
        try:
            converted.append(translateMouse_symbol[gene][1])
            keeped.append(gene)
        except:
            missing.append(gene)
            print("", end=f"\r{len(missing)} / {len(liste_names)} genes missing for translation")
            
    return converted, missing, keeped
        

#%%
######     GET HUMAN CONTROL DATA      #######
control = pd.read_csv("breast_dataframe_control.csv")
######     GET MICE DATA      #######
mice = pd.read_csv("MICE.csv")
######     GET HUMAN TUMORS DATA      #######
breast = pd.read_csv("BREAST_tumors.csv")

#%%
control = control.T
mice = mice.T
breast = breast.T

mice.columns = mice.iloc[0]
mice = mice.drop("Unnamed: 0", axis = 0)
#%%

# =============================================================================
# Conversion RPM !! (fonction dans le fichier fonction.py)
# =============================================================================

RPM_mice = fc.rawToRPM(mice)
RPM_control = fc.rawToRPM(control)
RPM_breast = fc.rawToRPM(breast)


#%%
#Vérification des meme genes chez l'humain
for i in range(len(breast.index)):
    if RPM_breast.index[i] != RPM_control.index[i]:
        print(i)
#Conversion des gènes de souris en humain

mousse_to_human_genes, missing, newlist = mice_to_human(RPM_mice.index)

#%%
# =============================================================================
# Nouveaux DF souris
# =============================================================================

RPM_mice = RPM_mice.loc[newlist]
RPM_mice.index = mousse_to_human_genes

#%%
# =============================================================================
# CONCATENATE: Common genes (pour utiliser les meme genes chez les humains et les souris
# =============================================================================

yes = 0
no = 0
gene_in_common = []
gene_exclude = []
for gene in breast.index: #control à les meme genes
    gene = gene.split('.')[0]
    if gene in RPM_mice.index:
        yes +=1
        gene_in_common.append(gene)
    else:
        gene_exclude.append(gene)
        no += 1
    print("", end=f"\r{yes+ no} / {len(breast.columns)} complete, {yes} find, {no} not found")     
    
#%%

# =============================================================================
# DOUBLONS!!! Certains gènes sont présent en double chez la souris, pour éviter toute confusion ils sont retité
# =============================================================================
RPM_mice = RPM_mice.loc[gene_in_common] 

gene_in_common_no_DUP = []
already_drop = []   
for gene in RPM_mice.index:
    if gene in already_drop:
        continue
    elif  gene in gene_in_common_no_DUP:
        gene_in_common_no_DUP.remove(gene)
        already_drop.append(gene)
    elif gene not in gene_in_common_no_DUP:
        gene_in_common_no_DUP.append(gene)
    
    
#%%
# =============================================================================
# New DF
# =============================================================================

RPM_mice = RPM_mice.loc[gene_in_common_no_DUP] 
RPM_breast = RPM_breast.loc[gene_in_common_no_DUP]
RPM_control = RPM_control.loc[gene_in_common_no_DUP]

#%%

# =============================================================================
# SELECTION TNBC # Filtré au préalable en observant les expression des 3 récepteurs HER2, ESR1 et PGR
# =============================================================================

TNBC_indices = [0, 3, 7, 22, 31, 35, 41, 49, 51, 53, 54, 55, 56, 59, 61, 63, 65, 71, 76, 77, 78, 86, 88, 89, 91, 95, 98, 100, 102, 106, 116, 120, 129, 131, 132, 140, 142, 144, 161, 162, 166, 168, 169, 179, 188, 189, 194, 208, 209, 213, 221, 223, 234, 235, 237, 241, 243, 249, 252, 256, 278, 286, 300, 302, 311, 315, 323, 324, 327, 331, 335, 336, 338, 344, 345, 346, 348, 352, 360, 362, 396, 399, 402, 406, 412, 425, 431, 448, 450, 451, 454, 457, 460, 476, 477, 482, 492, 494, 495, 496, 499, 502, 503, 504, 517, 521, 522, 526, 532, 534, 540, 544, 550, 561, 568, 571, 573, 580, 581, 584, 590, 608, 620, 625, 627, 642, 644, 649, 661, 676, 678, 691, 693, 694, 706, 723, 725, 729, 738, 741, 743, 764, 766, 768, 776, 778, 781, 786, 797, 798, 802, 805, 814, 818, 819, 820, 821, 822, 828, 833, 851, 857, 888, 890, 892, 893, 894, 895, 900, 904, 911, 913, 916, 919, 925, 927, 928, 933, 936, 940, 942, 943, 944, 951, 953, 956, 957, 959, 965, 967, 973]
RPM_TNBC = RPM_breast[TNBC_indices]
RPM_breast = RPM_breast.drop(TNBC_indices, axis = 1)

#%%

# =============================================================================
# CONCATENATE: colnames
# =============================================================================
RPM_control_names = []
RPM_breast_names = []
RPM_TNBC_names = [] 

for i in range(len(RPM_control.columns)):
    RPM_control_names.append("control_"+str(i))
    
breast_names = []

for i in range(len(RPM_breast.columns)):
    RPM_breast_names.append("breast_"+str(i))
    
TNBC_names = []

for i in range(len(RPM_TNBC.columns)):
    RPM_TNBC_names.append("TNBC_"+str(i))

#%%

RPM_breast.columns = RPM_breast_names
RPM_control.columns = RPM_control_names
RPM_TNBC.columns = RPM_TNBC_names
RPM_mice.columns.name = ""

#%%

RPM_mice.to_csv("RPM_FINAL_mice.csv", index = True)
RPM_breast.to_csv("RPM_FINAL_breast.csv", index = True)
RPM_control.to_csv("RPM_FINAL_control.csv", index = True)
RPM_TNBC.to_csv("RPM_FINAL_tnbc.csv", index = True)
