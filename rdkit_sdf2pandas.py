"""
  PandasTools is an efficient way to load and analyze
  large SDF file.  
  CK 9/9/25
"""

import os
import sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import Descriptors
import pandas as pd
from rdkit.Chem import PandasTools
import matplotlib.pyplot as plt
import numpy as np

# Import SDF file directly as a DataFrame


"""
Read file in SDF format and return as Pandas data frame. If embedProps=True all properties also get embedded in Mol objects in the molecule column. If molColName=None molecules would not be present in resulting DataFrame (only properties would be read).

Sanitize boolean is passed on to Chem.ForwardSDMolSupplier sanitize. If neither molColName nor smilesName are set, sanitize=false.

https://rdkit.org/docs/source/rdkit.Chem.PandasTools.html#rdkit.Chem.PandasTools.LoadSDF
"""

df = PandasTools.LoadSDF("/Users/ck/Downloads/Div_set_50k.sdf", idName = 'IDNUMBER', sanitize=True, embedProps=True, removeHs=True)

print(df.head())
print(df.info())

mol = df['ROMol']  # this is the structure
print(Chem.MolToSmiles(mol[0]))

MWs = [Descriptors.ExactMolWt(mol) for mol in mol]
HBAs = [Descriptors.NumHAcceptors(mol) for mol in mol]
HBDs = [Descriptors.NumHDonors(mol) for mol in mol]
LogPs = [Descriptors.MolLogP(mol) for mol in mol]
title = [mol.GetProp("_Name") for mol in mol]
print(title[0], MWs[0], HBAs[ 0 ], HBDs[ 0 ], LogPs[ 0 ])

mwt = df['MW'].values.tolist()
new_mwt = [float(i) for i in mwt]  # convert dataframe object str format to float 
#print(type(new_mwt[50]))
clogp = df['cLogP'].values.tolist()
new_clogp = [float(i) for i in clogp]  # convert dataframe object str format to float 

figure, axis = plt.subplots(2,2)
axis[0,0].hist(MWs, bins=5)
#axis[0,0].set_xticks(np.arange(100,501,100))
#axis[0,0].set_yticks(np.arange(0,25001,5000))
axis[0,0].set_title("MWt",fontsize=7)
axis[0,1].hist(HBAs, bins=10)
#axis[0,1].set_xticks(np.arange(1,11,2))
#axis[0,1].set_yticks(np.arange(0,15000,5000))
axis[0,1].set_title("HB-Acceptors", fontsize=7)
axis[1,0].hist(HBDs, bins=10)
#axis[1,0].set_xticks(np.arange(1,6,1))
#axis[1,0].set_yticks(np.arange(0,25000,5000))
axis[1,0].set_title("HB-Donors", fontsize=7)
axis[1,1].hist(LogPs, bins=10)
#axis[1,1].set_xticks(np.arange(-6,6,2))
#axis[1,1].set_yticks(np.arange(0,15000,5000))
axis[1,1].set_title("cLogP", fontsize=7)
axis[0,0].tick_params(labelsize=5)
axis[0,1].tick_params(labelsize=5)
axis[1,0].tick_params(labelsize=5)
axis[1,1].tick_params(labelsize=5)
#plt.hist(new_mwt,bins=5)
#plt.hist(new_clogp,bins=10)
#plt.hist(LogPs,bins=10)
#plt.xticks(np.arange(min(new_mwt), max(new_mwt)+1, 100))
#plt.xticks(np.arange(100, 501, 100))
#plt.yticks(np.arange(0, 25001, 5000))
plt.xticks(np.arange(-5, 5.01, 2))
plt.yticks(np.arange(0, 25001, 5000))
#plt.xlabel('Molecular Weight')
#plt.xlabel('cLogP')
#plt.ylabel('Frequency')
#plt.title('ChemDiv Div Set 50K')
figure.suptitle('ChemDiv Div Set 50K', fontsize=20)
#plt.savefig('mwt.png')
figure.savefig("fig.png")
#plt.savefig('clogp.png')
#plt.show()
