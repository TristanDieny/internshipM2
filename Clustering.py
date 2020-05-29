# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 10:57:12 2020

@author: tristan
"""
import pandas as pd 
from scipy import stats
import numpy as np
from termcolor import colored
import functions as fc

mice = pd.read_csv("RPM_FINAL_mice.csv", index_col = 0 )
breast = pd.read_csv("RPM_FINAL_breast.csv", index_col = 0 )
control = pd.read_csv("RPM_FINAL_control.csv", index_col = 0 )
TNBC = pd.read_csv("RPM_FINAL_tnbc.csv", index_col = 0 )

#%%

# =============================================================================
# Genes selection based on stats: sélection des genes les moins différencié entre la souris et l'homme
# =============================================================================
allGenes2 = [] 

for i in range(len(control.index)):
    souris = list(mice[["Control_1", "Control_2"]].iloc[i])
    ctrl = list(control.iloc[i])
    pval = stats.ttest_ind(ctrl, souris)[1]
    if pval<= 0.05:
        allGenes2.append(control.index[i])

#%%
####################################################################
#######################  CLUSTERIZATION   ##########################
####################################################################

##### Structuration des données: récupération des genes fournis par Lehmann sous forme de dictionnaire.
#Chaque cluster possède une liste de genes UP régulé et une liste de gènes DOWN régulé
    
file = open("genes_subtypes.txt", "r", encoding="utf8")

index = 0

genes_subtypes = dict()
clusters = list()

for line in file:
    if index == 0:
        index = 1
        continue
    elif index == 1:
        line = line.strip('\n')
        line = line.split('\t')
        for cluster in line:
            if cluster != '':
                genes_subtypes[cluster] = [[],[]]       #pour chaque cluster on associe une liste de genes 'UP' et 'DOWN'
                clusters.append(cluster)
        index = 2
    elif index == 2 :
        index = 3
        continue
    else:
        line = line.strip('\n')
        line = line.split('\t')     
        for i in range(0,28,4):       
            if line[i] != '':
                if line[i+1] == "UP":
                    genes_subtypes[clusters[int(i/4)]][0].append(line[i])
                if line[i+1] == "DOWN":
                    genes_subtypes[clusters[int(i/4)]][1].append(line[i])
    

#%%
        
##### Vérification de la présence des gènes dans mes données, on enlève ceux qui n'y sont pas
        
for cluster in clusters:
    for i in range(2):
        present = 0
        absent = 0
        genes_temp = []
        for gene in genes_subtypes[cluster][i]:
            if gene in mice.index:
                present += 1
                genes_temp.append(gene)
            else:
                absent += 1
        genes_subtypes[cluster][i] = genes_temp
        print(cluster, '\t\t', present, "/", present + absent, '\t', len(genes_subtypes[cluster][i]))
#        print(cluster, "\t\t", len(genes_subtypes_unique[cluster][i]), "\t", len(genes_subtypes[cluster][i]))

    
    
#%%   

##### Creation d'un dictionnaire contenant les genes uniquement associé à un seul et unique cluster (correspond à la list "filtered" dans mon mémoire)
genes_subtypes_unique = dict()

for cluster in clusters:
    genes_subtypes_unique[cluster] = [[],[]]       
 
subtype_temp = list()

for i in range(2):
    for cluster1 in clusters:
        subtype_temp = list()
        genes_temp = list()
        for cluster2 in clusters:
            if cluster1 == cluster2:
                continue
            for gene in genes_subtypes[cluster2][i]:
                genes_temp.append(gene)
        genes_temp = list(dict.fromkeys(genes_temp))
        for gene in genes_subtypes[cluster1][i]:
            if gene not in genes_temp:
                subtype_temp.append(gene)
        genes_subtypes_unique[cluster1][i] = subtype_temp

for cluster in clusters:
    for i in range(2):
        print(cluster, "\t", i , "\t",  len(genes_subtypes_unique[cluster][i]))
#%%%

controls = ["Control_1", "Control_2"]
earlys = ["Early_1", "Early_2", "Early_3"]
advanced = ["Advanced_1", "Advanced_2", "Advanced_3", "Advanced_4"]
#Selection des TNBC samples choisis dans mon mémoire suite aux résultats obtenus par TNBCtype
TNBC_samples = ["TNBC_177", "TNBC_126", "TNBC_29", "TNBC_9", "TNBC_139", "TNBC_88", "TNBC_24", "TNBC_115", "TNBC_42", "TNBC_181", "TNBC_12", "TNBC_2",]

#%%

##### récuparation DE LA MOYENNE DES CONTROLS des souris et des humains
    
mean_control = mice[controls].mean(axis = 1)
mean_control_human = control.mean(axis = 1)

#%%   

#######################  COMPARISON OF MEANS  ########################## MICE

alpha = 0.05

meanC = open("mean_comparaison.txt", "w+")  #Les résultats sont aussi écrit dans un fichier

print("Cluster \t\t\t UP(0)/DOWN(1)\tPvalue \t\tMean(ctrl)/Mean(test) \t\tscore")
for advance in advanced: 
    print("\n"+advance)
    for i in range(2):   
        for cluster in clusters:
            test = list(mice.loc[genes_subtypes[cluster][i]][advance])  #expression des genes pour le sample en question
            ctrl = list(mean_control.loc[genes_subtypes[cluster][i]])   #expression des genes des controls (moyenne déja calculé au dessus)
            pval = stats.wilcoxon(ctrl, test)[1]
            if pval <= alpha:                                     #alpha peut etre égal à 0.05 ou 1 en fonction de l'envie de restreindre au résultat significatif ou non
#                if i == 0 and np.mean(ctrl) < np.mean(test):     #Dans ce cas on affiche que les scores positif
                 if i == 0:                                       #Dans ce cas on affiche tout les résultats
                    print("{0:40} {1} \t".format(cluster, i), end = "")
                    print("{:2.2}\t\t".format(pval), end = "")
                    print("{:2.2f}\t/ {:2.2f}\t\t".format(np.mean(ctrl), np.mean(test)), end = "")    
                    print("{:2.2f}\t ".format(np.mean(test) - np.mean(ctrl)))
                    meanC.write(cluster + "\t" + advance + "\t" + str(i) + "\t" + str(np.mean(ctrl)) + "\t" + str(np.mean(test)) + "\t"  + "\n")
#                elif i == 1 and np.mean(ctrl) > np.mean(test):
                 if i == 1:
                    print("{0:40} {1} \t".format(cluster, i), end = "")
                    print("{:2.2}\t\t".format(pval), end = "")
                    print("{:2.2f}\t/ {:2.2f}\t\t".format(np.mean(ctrl), np.mean(test)), end = "")    
                    print("{:2.2f}\t".format(np.mean(ctrl) - np.mean(test))) 
                    meanC.write(cluster + "\t" + advance + "\t" + str(i) + "\t" + str(np.mean(ctrl)) + "\t" + str(np.mean(test)) + "\t" +  "\n")
    print("\n=======================\n")
    
meanC.close()

#%%   

#######################  COMPARISON OF MEANS  ########################## HUMAN
# Meme fonctionnement que pour les souris
alpha = 1

meanC = open("mean_comparaison_human.txt", "w+")

print("Cluster \t\t\t UP(0)/DOWN(1) \t Sample \tPvalue \t\tMean(ctrl)/Mean(test) \t\tscore")
for TNBC_sample in TNBC_samples: 
    print("\n"+TNBC_sample)
    for i in range(2):   
        for cluster in clusters:
            test = list(TNBC.loc[genes_subtypes[cluster][i]][TNBC_sample])
            ctrl = list(mean_control_human.loc[genes_subtypes[cluster][i]])
            pval = stats.wilcoxon(ctrl, test)[1]
            if pval <= alpha:
#                if i == 0 and np.mean(ctrl) < np.mean(test):
                 if i == 0:
                    print("{0:40} {1} \t".format(cluster, i), end = "")
                    print("{:2.2}\t\t".format(pval), end = "")
                    print("{:2.2f}\t/ {:2.2f}\t\t".format(np.mean(ctrl), np.mean(test)), end = "")    
                    print("{:2.2f}\t ".format(np.mean(test) - np.mean(ctrl)))
                    meanC.write(cluster + "\t" + TNBC_sample + "\t" + str(i) + "\t" + str(np.mean(ctrl)) + "\t" + str(np.mean(test)) + "\t"  + str(np.mean(test) - np.mean(ctrl)) +  "\n")
#                elif i == 1 and np.mean(ctrl) > np.mean(test):
                 if i == 1:
                    print("{0:40} {1} \t".format(cluster, i), end = "")
                    print("{:2.2}\t\t".format(pval), end = "")
                    print("{:2.2f}\t/ {:2.2f}\t\t".format(np.mean(ctrl), np.mean(test)), end = "")    
                    print("{:2.2f}\t".format(np.mean(ctrl) - np.mean(test))) 
                    meanC.write(cluster + "\t" + TNBC_sample + "\t" + str(i) + "\t" + str(np.mean(ctrl)) + "\t" + str(np.mean(test)) + "\t" +  str(np.mean(ctrl) - np.mean(test))+ "\n")
    print("\n=======================\n")
    
meanC.close()


#%%
        
####################################################################
######################  DEGs comparison   ##########################
####################################################################       

##### DICO CONVERSION SYMBOL MOUSSE TO HUMAN (pour les genes récuperer par DESEQ2)

symbol_mouse_to_human = dict()

file = open("mouse_to_human.txt")

for line in file:
    line = line.split("\t")
    symbol_mouse_to_human[line[1]] = line[3].strip("\n")

file.close

#%%

#Récupération des genes obtenus avec DESeq2
files_DOWN = [open("DESEQ2/DOWN_genes_control_vs_tumor1.csv"),
         open("DESEQ2/DOWN_genes_control_vs_tumor2.csv"),
         open("DESEQ2/DOWN_genes_control_vs_tumor3.csv"),
         open("DESEQ2/DOWN_genes_control_vs_tumor4.csv")]
files_UP = [open("DESEQ2/UP_genes_control_vs_tumor1.csv"),
         open("DESEQ2/UP_genes_control_vs_tumor2.csv"),
         open("DESEQ2/UP_genes_control_vs_tumor3.csv"),
         open("DESEQ2/UP_genes_control_vs_tumor4.csv")]

DESEQ2 = open("deseq2_analysis.txt", "w+")  #Ecriture des résultats dans un fichier

for i in range(4):    #Pour les 4 samples "Advanced"
    #Get donw genes in the file
    DE_genes_DOWN = list()
    index = 0
    for line in files_DOWN[i]:
        if index == 0:
            index = 1 
            continue
        line = line.split(',')
        DE_genes_DOWN.append(line[1].strip('"')) 
    print(files_DOWN[i].name, "\t", len(DE_genes_DOWN), "DOWN regulated genes")
    files_DOWN[i].close()
    #Get UP genes in the file
    DE_genes_UP = list()
    index = 0
    for line in files_UP[i]:
        if index == 0:
            index = 1 
            continue
        line = line.split(',')
        DE_genes_UP.append(line[1].strip('"')) 
    print(files_UP[i].name, "\t", len(DE_genes_UP), "UP regulated genes")
    files_UP[i].close()
    
    print("cluster\t\t\t\t\t    regulation\t  % of present genes    % of present genes for unique genes" )

# =============================================================================
#     DE_genes_UP = X.index
#     DE_genes_DOWN = X.index
# =============================================================================
    #Ici on regarde pour chcaun des cluster si il y a des genes en commun entre la liste de Lehmann et la liste obtenus avec DESEQ2
    for cluster in clusters:
        okUp, okDown = 0, 0
        okUniqueUp, okUniqueDown = 0, 0
        not_converted = 0
        for gene in DE_genes_UP:
            try:
                if symbol_mouse_to_human[gene] in genes_subtypes[cluster][0]:
#                if gene in genes_subtypes[cluster][0]:
                    okUp += 1
                if symbol_mouse_to_human[gene] in genes_subtypes_unique[cluster][0]:
#                if gene in genes_subtypes_unique[cluster][0]:
                    okUniqueUp += 1
            except:
                not_converted +=1
        for gene in DE_genes_DOWN:
            try:
                if symbol_mouse_to_human[gene] in genes_subtypes[cluster][1]:
#                if gene in genes_subtypes[cluster][1]:
                    okDown += 1
                if symbol_mouse_to_human[gene] in genes_subtypes_unique[cluster][1]:
#                if gene in genes_subtypes_unique[cluster][1]:
                    okUniqueDown += 1
            except:
                not_converted +=1
#        print("{0:40}".format(cluster), "\tUP\t", okUp , "/" , len(genes_subtypes[cluster][0]), "\t" , okUniqueUp , "/", len(genes_subtypes_unique[cluster][0]))
#        print("{0:40}".format(cluster), "\tDOWN\t", okDown , "/" , len(genes_subtypes[cluster][1]), "\t" , okUniqueDown , "/", len(genes_subtypes_unique[cluster][1]))7
        outputUp = "{0:40}".format(cluster)+ "\tUP\t\t"+ "{:2.2f}".format(okUp/len(genes_subtypes[cluster][0])*100)+ "\t\t" + "{:2.2f}".format(okUniqueUp/len(genes_subtypes_unique[cluster][0])*100)
        outputDown = "{0:40}".format(cluster)+ "\tDOWN\t\t"+ "{:2.2f}".format(okDown/len(genes_subtypes[cluster][1])*100)+ "\t\t" + "{:2.2f}".format(okUniqueDown/len(genes_subtypes_unique[cluster][1])*100)
        
        W_outputUp = cluster+ "\tUP\t"+ str(okUp/len(genes_subtypes[cluster][0])*100)+ "\t" + str(okUniqueUp/len(genes_subtypes_unique[cluster][0])*100) + "\n"
        W_outputDown = cluster + "\tDOWN\t"+ str(okDown/len(genes_subtypes[cluster][1])*100) + "\t" + str(okUniqueDown/len(genes_subtypes_unique[cluster][1])*100) + "\n"

#        DESEQ2.write(W_outputUp)
#        DESEQ2.write(W_outputDown)
        if okUp/len(genes_subtypes[cluster][0])*100 >= 20 or okUniqueUp/len(genes_subtypes_unique[cluster][0])*100 >= 20:   #Ecriture des score élevé en rouge
            print(colored(str('\033[1m'+outputUp), 'red'))
        else:
            print(outputUp)       
        if okDown/len(genes_subtypes[cluster][1])*100 >= 20 or okUniqueDown/len(genes_subtypes_unique[cluster][1])*100 >= 20:   #Ecriture des score élevé en rouge
            print(colored(str('\033[1m'+outputDown), 'red'))
        else:
            print(outputDown)
#    print("conversions are missing for:", not_converted, "genes")
    print() 
    file.close
    DESEQ2.close()
  
#%%

#Exactement pareil que précedement mais pour les humains

files_DOWN = [open("DESEQ2/Human/DOWN_genes_control_vs_tumor177.csv"),
         open("DESEQ2/Human/DOWN_genes_control_vs_tumor126.csv"),
         open("DESEQ2/Human/DOWN_genes_control_vs_tumor29.csv"),
         open("DESEQ2/Human/DOWN_genes_control_vs_tumor9.csv"),
         open("DESEQ2/Human/DOWN_genes_control_vs_tumor139.csv"),
         open("DESEQ2/Human/DOWN_genes_control_vs_tumor88.csv"),
         open("DESEQ2/Human/DOWN_genes_control_vs_tumor24.csv"),
         open("DESEQ2/Human/DOWN_genes_control_vs_tumor115.csv"),
         open("DESEQ2/Human/DOWN_genes_control_vs_tumor42.csv"),
         open("DESEQ2/Human/DOWN_genes_control_vs_tumor181.csv"),
         open("DESEQ2/Human/DOWN_genes_control_vs_tumor12.csv"),
         open("DESEQ2/Human/DOWN_genes_control_vs_tumor2.csv")]
files_UP = [open("DESEQ2/Human/UP_genes_control_vs_tumor177.csv"),
         open("DESEQ2/Human/UP_genes_control_vs_tumor126.csv"),
         open("DESEQ2/Human/UP_genes_control_vs_tumor29.csv"),
         open("DESEQ2/Human/UP_genes_control_vs_tumor9.csv"),
         open("DESEQ2/Human/UP_genes_control_vs_tumor139.csv"),
         open("DESEQ2/Human/UP_genes_control_vs_tumor88.csv"),
         open("DESEQ2/Human/UP_genes_control_vs_tumor24.csv"),
         open("DESEQ2/Human/UP_genes_control_vs_tumor115.csv"),
         open("DESEQ2/Human/UP_genes_control_vs_tumor42.csv"),
         open("DESEQ2/Human/UP_genes_control_vs_tumor181.csv"),
         open("DESEQ2/Human/UP_genes_control_vs_tumor12.csv"),
         open("DESEQ2/Human/UP_genes_control_vs_tumor2.csv")]

DESEQ2 = open("deseq2_analysis.txt", "w+")

for i in range(12):
    #Get donw genes
    DE_genes_DOWN = list()
    index = 0
    for line in files_DOWN[i]:
        if index == 0:
            index = 1 
            continue
        line = line.split(',')
        DE_genes_DOWN.append(line[1].strip('"')) 
    print(files_DOWN[i].name, "\t", len(DE_genes_DOWN), "DOWN regulated genes")
    files_DOWN[i].close()
    #Get UP genes
    DE_genes_UP = list()
    index = 0
    for line in files_UP[i]:
        if index == 0:
            index = 1 
            continue
        line = line.split(',')
        DE_genes_UP.append(line[1].strip('"')) 
    print(files_UP[i].name, "\t", len(DE_genes_UP), "UP regulated genes")
    files_UP[i].close()
    
    print("cluster\t\t\t\t\t    regulation\t  % of present genes    % of present genes for unique genes" )

# =============================================================================
#     DE_genes_UP = X.index
#     DE_genes_DOWN = X.index
# =============================================================================
    for cluster in clusters:
        okUp, okDown = 0, 0
        okUniqueUp, okUniqueDown = 0, 0
        not_converted = 0
        for gene in DE_genes_UP:
            try:
                if gene in genes_subtypes[cluster][0]:
#                if gene in genes_subtypes[cluster][0]:
                    okUp += 1
                if symbol_mouse_to_human[gene] in genes_subtypes_unique[cluster][0]:
#                if gene in genes_subtypes_unique[cluster][0]:
                    okUniqueUp += 1
            except:
                not_converted +=1
        for gene in DE_genes_DOWN:
            try:
                if gene in genes_subtypes[cluster][1]:
#                if gene in genes_subtypes[cluster][1]:
                    okDown += 1
                if gene in genes_subtypes_unique[cluster][1]:
#                if gene in genes_subtypes_unique[cluster][1]:
                    okUniqueDown += 1
            except:
                not_converted +=1
#        print("{0:40}".format(cluster), "\tUP\t", okUp , "/" , len(genes_subtypes[cluster][0]), "\t" , okUniqueUp , "/", len(genes_subtypes_unique[cluster][0]))
#        print("{0:40}".format(cluster), "\tDOWN\t", okDown , "/" , len(genes_subtypes[cluster][1]), "\t" , okUniqueDown , "/", len(genes_subtypes_unique[cluster][1]))7
        outputUp = "{0:40}".format(cluster)+ "\tUP\t\t"+ "{:2.2f}".format(okUp/len(genes_subtypes[cluster][0])*100)+ "\t\t" + "{:2.2f}".format(okUniqueUp/len(genes_subtypes_unique[cluster][0])*100)
        outputDown = "{0:40}".format(cluster)+ "\tDOWN\t\t"+ "{:2.2f}".format(okDown/len(genes_subtypes[cluster][1])*100)+ "\t\t" + "{:2.2f}".format(okUniqueDown/len(genes_subtypes_unique[cluster][1])*100)
        
        W_outputUp = cluster+ "\tUP\t"+ str(okUp/len(genes_subtypes[cluster][0])*100)+ "\t" + str(okUniqueUp/len(genes_subtypes_unique[cluster][0])*100) + "\n"
        W_outputDown = cluster + "\tDOWN\t"+ str(okDown/len(genes_subtypes[cluster][1])*100) + "\t" + str(okUniqueDown/len(genes_subtypes_unique[cluster][1])*100) + "\n"

#        DESEQ2.write(W_outputUp)
#        DESEQ2.write(W_outputDown)
        if okUp/len(genes_subtypes[cluster][0])*100 >= 20 or okUniqueUp/len(genes_subtypes_unique[cluster][0])*100 >= 20:
            print(colored(str('\033[1m'+outputUp), 'red'))
        else:
            print(outputUp)       
        if okDown/len(genes_subtypes[cluster][1])*100 >= 20 or okUniqueDown/len(genes_subtypes_unique[cluster][1])*100 >= 20:
            print(colored(str('\033[1m'+outputDown), 'red'))
        else:
            print(outputDown)
#    print("conversions are missing for:", not_converted, "genes")
    print() 
    file.close
    DESEQ2.close()
       
#%%
# =============================================================================
# DEEP LEARNING : TRAINNING
# =============================================================================

import numpy as np
import tensorflow as tf
import tensorflow.keras as keras
import numpy
from sklearn.model_selection import StratifiedKFold   #cross validation
from tensorflow.keras.callbacks import EarlyStopping
from tensorflow.keras.callbacks import ModelCheckpoint
from tensorflow.keras.models import load_model
from sklearn.preprocessing import OneHotEncoder
from keras.utils import to_categorical
import os
     
##### Récupération des données
        
# GET Y VALUES: obtenus depuis TNBCtype pour les humains
file = open("V2_GOOD_results.csv")
Y = []
for line in file:
    Y.append(line.split(",")[1])
    
Y =Y[1:181]

for i in range(len(Y)):
    if Y[i] == "UNS" :
        Y[i] = 0
    elif Y[i] == "LAR":
        Y[i] = 1
    elif Y[i] == "BL1":
        Y[i] = 2
    elif Y[i] == "BL2":
        Y[i] = 3
    elif Y[i] == "M":
        Y[i] = 4
    elif Y[i] == "IM":
        Y[i] = 5
    elif Y[i] == "MSL":
        Y[i] = 6 

for i in range(len(control.columns)):
    Y.append(7)
for i in range(30):   #On ne récupere que 30 sur les 1041 samples de breast cancer puisqu'il ne sont pas des triple négatif est permette uniquement de servir de test pour les souris 
    Y.append(8)
 
#Unification de tous les genes UP ou DOWN de tous les cluster en une seul liste (le modele ne vas s'entrainner qu'avec ces genes pertinant plutot que sur les 15000 genes qui compléxifi énormément le réseau de neuronnes)
allGenes1 = []
for i in  genes_subtypes_unique.values():
    allGenes1 += i
allGenes = []
for i in  allGenes1:
    allGenes += i

#Soit on utilise les 3000 genes de Lehmann, soit on utilise les 5000 genes peut variable entre les souris et les huimains pour une bonne prédiction
#Liste de gènes correspondant à Lehmanns
final_DF = pd.concat([TNBC.loc[allGenes], control.loc[allGenes], breast.loc[allGenes].iloc[:,0:30], mice.loc[allGenes]], axis=1)

#Liste de gènes peut différencié entre les controles souris et humains
#final_DF = pd.concat([TNBC.loc[allGenes2], control.loc[allGenes2], breast.loc[allGenes2].iloc[:,0:30], mice.loc[allGenes2]], axis=1)

X = np.asarray(final_DF.T)
Y = np.asarray(Y)

norm = np.linalg.norm(X)
X = X/norm
Y = to_categorical(Y)

#%%
x = X[0:-10]        #Retrait des souris qui ne sont pas labelé

activation = 'relu'
optimizer = keras.optimizers.Adam(lr=1e-3)
nb_noeud = len(final_DF.index)
nb_couche = 2
batch_size = 50
dropout_rate = 0.2
nb_cv = 2
epoch = 110000000   #On utilise le early stop de toute façon

es = EarlyStopping(monitor='categorical_accuracy', mode='max', patience=500, verbose=1)
kfold = StratifiedKFold(n_splits=nb_cv, shuffle=True)

mc = ModelCheckpoint('Model/best_model.h5', monitor='categorical_accuracy', mode='max', verbose=0, save_best_only=True) #Pour sauvegarder le meilleur modele
                                                            
model = keras.models.Sequential()
model.add(keras.layers.Dropout(dropout_rate, input_shape=(len(final_DF.index),)))
model.add(keras.layers.Dense(nb_noeud, kernel_initializer = "uniform", activation=activation, input_dim=len(final_DF.index)))   #, use_bias=False

for i in range(nb_couche-1):
     model.add(keras.layers.Dense(int(nb_noeud/2), kernel_initializer = "uniform", activation=activation)) #
          
model.add(keras.layers.Dense(9, activation=tf.nn.softmax))
model.compile(optimizer = optimizer,
loss=keras.losses.categorical_crossentropy,
metrics=[keras.metrics.categorical_accuracy])
            
 #model.summary()
model.fit(x, Y,  callbacks=[es, mc], epochs=epoch, verbose=1, batch_size = batch_size)
            

#%%
# =============================================================================
# DEEP LEARNING : PREDICT
# =============================================================================

best_model = load_model('Model/best_model.h5')
miceDF = X[-10:]

predictions = best_model.predict(miceDF)

for i in range(len(predictions)):
    print(mice.columns[i], "->", np.argmax(predictions[i])) #Affichage de la prédiction pour chacun des samples

