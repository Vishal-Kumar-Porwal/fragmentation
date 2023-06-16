
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from pandas import DataFrame
import math as m
import numpy as np
import os
import Bio
from Bio import *
from Bio.PDB.PDBParser import PDBParser
import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)
from mpl_toolkits import mplot3d
import ast
import pylab as pl
import random
from random import randint
#import babel
from os import listdir
from chemml.chem import Molecule
#from rdkit import Chem 




file_xyz = input('path of the xyz file:')
#file_xyz = str('A234.xyz')

file_connection= input('path of the file with only connection information:')
#file_connection= str('connection-A234.txt')

#rotors_output= input('path of the file where output of rotors will be stored:')
#rotors_output= str('rotors.txt')



connection_dict={}
n_connection=1
with open(file_connection,'r') as something:
    
    for line in something:
        # print(line)
        line = line.split()
        
        if line!= []:
            line_1 = []
            for j in line:
                j = float(j)
                line_1.append(j)
            connection_dict[n_connection] = line_1
            n_connection=n_connection+1
            # print(structure_dict)
        else:
            break
        
structure_dict={}
n_structure=1
with open(file_xyz,'r') as something:
    
    for line in something:
        # print(line)
        line = line.split()
        if line!= []:
            structure_dict[n_structure] = line
            n_structure=n_structure+1
            # print(structure_dict)
        else:
            break

#print('structure_dict=',structure_dict)# 
#print('connection_dict=',connection_dict)

bonds_dataframe = pd.DataFrame(index=connection_dict.keys(),columns=connection_dict.keys())
#print(bonds_dataframe)

for i in connection_dict.keys():
    for j in range(len(connection_dict[i])):
        
        if j%2 == 0 and j > 0:
            bonds_dataframe.at[i,connection_dict[i][j-1]] = connection_dict[i][j]
            bonds_dataframe.at[connection_dict[i][j-1],i] = connection_dict[i][j]
        else:
            bonds_dataframe.at[i,connection_dict[i][j]] = 0
            bonds_dataframe.at[connection_dict[i][j],i] = 0


bonds_dataframe = bonds_dataframe.replace(np.nan,0)

#print(bonds_dataframe)        

def find_bond(dataframe,atom1):
    list_atoms=[] 
    for atom2 in dataframe.columns:
        #if dataframe[atom1][atom2] >= 1:
            #print('still trying','atom1=',atom1,'atom2=',atom2)
          
        #    list_atoms.append(atom2)
            #return atom2
        #print(dataframe[atom1][atom2])
        if dataframe[atom1][atom2] > 0:
            #print('still trying','atom1=',atom1,'atom2=',atom2)
            #print('bond found with',atom2)
            list_atoms.append(atom2)
            
        elif dataframe[atom1][atom2] == 0:
             #print('no bond found with',atom2) 
             continue
             #print('still trying','atom1=',atom1,'atom2=',atom2)
                            
  
        else:
            continue
            #print('error: atoms indices are wrong, correctng and still trying')
    
    return list_atoms
    
#print(set([1,2])==set([2,1]))
    
for i in bonds_dataframe.columns:
    bondds = find_bond(bonds_dataframe,i)
#    print('bonds of:',i,structure_dict[i][0],bondds)


def linked_atoms(dataframe,atom1,atom2):
    list_linked_atoms_atom1 = find_bond(dataframe,atom1) 
    
    
    list_linked_atoms_atom2 = find_bond(dataframe,atom2)
    
    
    dict_linked_atoms={}
    dict_linked_atoms[atom1] = list_linked_atoms_atom1
    dict_linked_atoms[atom2] = list_linked_atoms_atom2
    total_length = len(dict_linked_atoms[atom1]) + len(dict_linked_atoms[atom2]) 
  
    while total_length != len(dataframe.columns):
        #print(total_length) 
        #for atom in dataframe.columns:
        
        atoms = random.choice(dataframe.columns)
                #print(atoms)
        if atoms == atom2:
           #print(atoms,"same atom; still trying")
            continue 
           
        elif atoms == atom1:
           #print(atoms,"same atom; still trying")
            continue
           
        elif atoms in list_linked_atoms_atom1:
           #print(atoms,"already considered; still trying")
            continue
           
        elif atoms in list_linked_atoms_atom2:
           #print(atoms,"already considered; still trying")
           
            continue
        else:
            #atoms = random.choice(dataframe.columns)
        
            for linked_atoms_atom1 in list_linked_atoms_atom1:
               
                if dataframe[linked_atoms_atom1][atoms] >= 1:
                    #print(atoms,"atoms linked to", atom1)
                    list_linked_atoms_atom1.append(atoms)
                else:
                    continue
                    
            
            for linked_atoms_atom2 in list_linked_atoms_atom2:
                if dataframe[linked_atoms_atom2][atoms] >= 1:
                    #print(atoms,"atoms linked to", atom2)
                    list_linked_atoms_atom2.append(atoms)
                else:
                    continue
        total_length = len(dict_linked_atoms[atom1]) + len(dict_linked_atoms[atom2])    
            #atoms = random.choice(dataframe.columns) 
    
             
        
    list_linked_atoms_atom1 = list_linked_atoms_atom1.remove(atom2)
    list_linked_atoms_atom2 = list_linked_atoms_atom2.remove(atom1)     
    #print(list_linked_atoms_atom1)
    #print(list_linked_atoms_atom2)  
    return dict_linked_atoms

#print(linked_atoms(bonds_dataframe,13,14))
                

    

hydrogen = 'H'
tuples=[]
corner_atoms_i=[]
corner_atoms_j=[]
for i in bonds_dataframe.columns:
    #print("i=",i)
    no_bonds_i=find_bond(bonds_dataframe,i)
 
     
        #print(i,structure_dict[i],connection_dict[i][2])
    
    for j in bonds_dataframe.columns:
        #print("j=",j)
        no_bonds_j=find_bond(bonds_dataframe,j)
        #if len(no_bonds_j) <= 1 :
        #    no_bonds_i=find_bond(bonds_dataframe,i)
        #    corner_atoms_j.append(j)
        if set([i,j]) in tuples:
            #print('already used this combination',i,j) 
            continue
        elif bonds_dataframe[i][j] == 0 :
            tuples.append(set([i,j]))
            continue
        elif bonds_dataframe[i][j] == 1 :
            tuples.append(set([i,j])) 
            dict_1 = linked_atoms(bonds_dataframe,i,j)
#            print(dict_1)
            multiplicity = int(1 + 2*1*(1/2))
#            list_of_files = listdir('*.xyz')
#            for files in list_of_files:
                 
            with open('XYZ/%s%s_%s%s_%s.xyz'% (structure_dict[i][0],i,structure_dict[j][0],j,multiplicity),'w') as file:
                #file.write(str(str(0)+" "+str(multiplicity)))
                file.write(str(len(dict_1[i])+1)) 
                file.write("\n")
                file.write("\n")
                file.write(str(" "+structure_dict[i][0]+"                 "+structure_dict[i][1]+"   "+structure_dict[i][2]+"   "+structure_dict[i][3]))
                file.write("\n")
            for elements in dict_1[i]:
                with open('XYZ/%s%s_%s%s_%s.xyz'% (structure_dict[i][0],i,structure_dict[j][0],j,multiplicity),'a') as file:
                   
                    file.write(str(" "+structure_dict[elements][0]+"                 "+structure_dict[elements][1]+"   "+structure_dict[elements][2]+"   "+structure_dict[elements][3]))
                    file.write("\n")

            with open('XYZ/%s%s_%s%s_%s.xyz'% (structure_dict[j][0],j,structure_dict[i][0],i,multiplicity),'w') as file:
                #file.write(str(str(0)+" "+str(multiplicity)))
                file.write(str(len(dict_1[j])+1))
                file.write("\n")
                file.write("\n")
                file.write(str(" "+structure_dict[j][0]+"                 "+structure_dict[j][1]+"   "+structure_dict[j][2]+"   "+structure_dict[j][3]))
                file.write("\n")
            for elements in dict_1[j]:
                with open('XYZ/%s%s_%s%s_%s.xyz'% (structure_dict[j][0],j,structure_dict[i][0],i,multiplicity),'a') as file:
                   
                    file.write(str(" "+structure_dict[elements][0]+"                 "+structure_dict[elements][1]+"   "+structure_dict[elements][2]+"   "+structure_dict[elements][3]))
                    file.write("\n")
                
            
        elif bonds_dataframe[i][j] == 2 :                                  
             dict_2 = linked_atoms(bonds_dataframe,i,j)
#             print(dict_2)
             multiplicity = int(1 + 2*2*(1/2))
             with open('XYZ/%s%s_%s%s_%s.xyz'% (structure_dict[i][0],i,structure_dict[j][0],j,multiplicity),'w') as file:
                 #file.write(str(str(0)+" "+str(multiplicity))) 
                 file.write(str(len(dict_2[i])+1))
                 file.write("\n")
                 file.write("\n")
                 file.write(str(" "+structure_dict[i][0]+"                 "+structure_dict[i][1]+"   "+structure_dict[i][2]+"   "+structure_dict[i][3]))
                 file.write("\n")
             for elements in dict_2[i]:
                 with open('XYZ/%s%s_%s%s_%s.xyz'% (structure_dict[i][0],i,structure_dict[j][0],j,multiplicity),'a') as file:
                    
                     file.write(str(" "+structure_dict[elements][0]+"                 "+structure_dict[elements][1]+"   "+structure_dict[elements][2]+"   "+structure_dict[elements][3]))
                     file.write("\n")
   
             with open('XYZ/%s%s_%s%s_%s.xyz'% (structure_dict[j][0],j,structure_dict[i][0],i,multiplicity),'w') as file:
                 #file.write(str(str(0)+" "+str(multiplicity)))
                 file.write(str(len(dict_2[j])+1))
                 file.write("\n")
                 file.write("\n")
                 file.write(str(" "+structure_dict[j][0]+"                 "+structure_dict[j][1]+"   "+structure_dict[j][2]+"   "+structure_dict[j][3]))
                 file.write("\n")
             for elements in dict_2[j]:
                 with open('XYZ/%s%s_%s%s_%s.xyz'% (structure_dict[j][0],j,structure_dict[i][0],i,multiplicity),'a') as file:
                    
                     file.write(str(" "+structure_dict[elements][0]+"                 "+structure_dict[elements][1]+"   "+structure_dict[elements][2]+"   "+structure_dict[elements][3]))
                     file.write("\n")
           
             tuples.append(set([i,j]))
        else:
            continue 


list_inchi = []
list_of_files = listdir('XYZ/')
#print(list_of_files)
for files in list_of_files:
    print(files)
    
    mol = Molecule('XYZ/%s' % files, 'xyz')
    
    mol.to_smiles()
    mol.to_smarts()
    mol.to_inchi()
    inchi_form = mol.inchi
    
    if inchi_form in list_inchi:
        os.remove('XYZ/%s' % files)
        
    else:
        list_inchi.append(inchi_form)
        
