"""
  This script reads a commercial library ChemDiv's Div_set_50k.sdf
  and performs analysis, structure conversion etc.,

  This script was tested with: python3 (3.13.5), numpy (2.1.3), 
  matplotlib (3.10.0), pandas (2.2.3) and rdkit (2025.03.6) 

  In the main function below, all relevant functions are commented.
  Depending on your needs, you can un-comment and use.
  CK, 9/29/25
"""

import os
import sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
import pandas as pd
from rdkit.Chem import PandasTools
import matplotlib.pyplot as plt
import numpy as np

def make_png(structure):
    img = Draw.MolToImage(structure)
    img.save('img.png')

def to_smiles(structures, titles, MWs, LogPs, tpsa, rotb, HBAs, HBDs):
    smi_out = open('test.csv', 'w')
    print("Smiles,IDNUMBER, MWt, cLogP, TPSA, RotB, HBA, HBD", file=smi_out)
    for i in range(len(structures)):
        print(Chem.MolToSmiles(structures[ i ]), titles[ i ], MWs[ i ], LogPs[ i ],
                  tpsa[ i ], rotb[ i ], HBAs[ i ], HBDs[ i ], sep=',', file=smi_out) 

def plot_properties(structures):
    MWs = [Descriptors.ExactMolWt(structure) for structure in structures]
    HBAs = [Descriptors.NumHAcceptors(structure) for structure in structures]
    HBDs = [Descriptors.NumHDonors(structure) for structure in structures]
    LogPs = [Descriptors.MolLogP(structure) for structure in structures]
    tpsa = [Descriptors.TPSA(structure) for structure in structures]
    rotb = [Descriptors.NumRotatableBonds(structure) for structure in structures]
    figure, axis = plt.subplots(2,3, figsize=(15,10))
    axis[0,0].hist(MWs, bins=5)
    #axis[0,0].set_title('MWt', fontsize=7)
    axis[0,0].set_xlabel('MWt', fontsize=7)
    axis[0,0].tick_params(labelsize=10)
    axis[0,1].hist(LogPs, bins=10)
    #axis[0,1].set_title('cLogP', fontsize=7)
    axis[0,1].set_xlabel('cLogP', fontsize=7)
    axis[0,1].tick_params(labelsize=10)
    axis[0,2].hist(tpsa, bins=5)
    #axis[0,2].set_title('TPSA', fontsize=7)
    axis[0,2].set_xlabel('TPSA', fontsize=7)
    axis[0,2].tick_params(labelsize=10)
    axis[1,0].hist(HBAs, bins=10)
    axis[1,0].tick_params(labelsize=10)
    #axis[1,0].set_title('Hbond Acceptors', fontsize=7)
    axis[1,0].set_xlabel('Hbond Acceptors', fontsize=7)
    axis[1,1].hist(HBDs, bins=10)
    #axis[1,1].set_title('Hbond Donors', fontsize=7)
    axis[1,1].set_xlabel('Hbond Donors', fontsize=7)
    axis[1,1].tick_params(labelsize=10)
    axis[1,2].hist(rotb, bins=10)
    #axis[1,2].set_title('Num. Rotatable Bonds', fontsize=7)
    axis[1,2].set_xlabel('Num. Rotatable Bonds', fontsize=7)
    axis[1,2].tick_params(labelsize=10)
    #plt.show()
    figure.savefig('props.png') 

def substr_search(structures):
    matches = Chem.SDWriter('substr_matches.sdf')
    substr = input("SDF or SMILES filename [with .sdf/.smi extension]  of a substructure: ")
    if "smi" in substr:
       template = Chem.MolFromSmiles(substr) 
    else:
       template = Chem.MolFromMolFile(substr) 
    substr_match = []
    for structure in structures:
        if structure.HasSubstructMatch(template):
           substr_match.append(structure)
           matches.write(structure)

def minim_str(structure):
    strH = Chem.AddHs(structure)
    AllChem.EmbedMolecule(strH)
    return strH

def confgen_str(structure, num_conf):
    confs = Chem.SDWriter('confs.sdf')
    strH = Chem.AddHs(structure)
    cids = AllChem.EmbedMultipleConfs(strH, numConfs=num_conf)
    print('Num Conformers = ', len(cids))
    rmslist=[]
    AllChem.AlignMolConformers(strH, RMSlist=rmslist)
    print('Num conformers aligned to first str = ',len(rmslist))
    #minimize conformers
    minConf = AllChem.MMFFOptimizeMoleculeConfs(strH)
    print('Num minimized conformers = ',len(minConf)) 
    for i in range(len(minConf)):
        confs.write(strH, i) 
        #print(minConf[ i ])  # Tuple showing convergence (0 = not converged, 1 = converged) and energy

def main():
    sd_file = input("Give library file in SDF format: ")
    df = PandasTools.LoadSDF(sd_file)
    print(df.head())
    print(df.info())
# Following properties are given in the SDF file.  Therefore, I am using extracting and using them.
# If not give, use plot_properties function to calculate them.
    ids = df['IDNUMBER'].values.tolist()
    HBAs = df['HBA'].values.tolist()
    HBDs = df['HBD'].values.tolist()
    MWs = df['MW'].values.tolist()
    cLogPs = df['cLogP'].values.tolist()
    tpsa = df['TPSA'].values.tolist()
    mol = df['ROMol']
    rotb = [Descriptors.NumRotatableBonds(structure) for structure in mol]
    mol = df['ROMol']
    to_smiles(mol, ids, MWs, cLogPs, tpsa, rotb, HBAs, HBDs)
    #plot_properties(mol)  
    #substr_search(mol)
    #min_str = minim_str(mol[0])
    #confgen_str(mol[ 0 ],10)  # 10 conformers requested
    #make_png(min_str)
 
if __name__ == '__main__':
   main()
