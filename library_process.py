"""
  This script reads a commercial library ChemDiv's Div_set_50k.sdf
  and performs analysis, structure conversion etc.,
  CK, 9/11/25
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

def make_png(structure):
    img = Draw.MolToImage(structure)
    img.save('img.png')

def to_smiles(structures, titles):
    smi_out = open('test.smi', 'w')
    print("Smiles,IDNUMBER", file=smi_out)
    for i in range(len(structures)):
        print(Chem.MolToSmiles(structures[ i ]), titles[ i ], sep=',', file=smi_out) 

def plot_properties(structures):
    MWs = [Descriptors.ExactMolWt(structure) for structure in structures]
    HBAs = [Descriptors.NumHAcceptors(structure) for structure in structures]
    HBDs = [Descriptors.NumHDonors(structure) for structure in structures]
    LogPs = [Descriptors.MolLogP(structure) for structure in structures]
    figure, axis = plt.subplots(2,2)
    axis[0,0].hist(MWs, bins=5)
    axis[0,0].set_title('MWt', fontsize=7)
    axis[0,0].tick_params(labelsize=5)
    axis[0,1].hist(LogPs, bins=10)
    axis[0,1].set_title('cLogP', fontsize=7)
    axis[0,1].tick_params(labelsize=5)
    axis[1,0].hist(HBAs, bins=10)
    axis[1,0].tick_params(labelsize=5)
    axis[1,0].set_title('Hbond Acceptors', fontsize=7)
    axis[1,1].hist(HBDs, bins=10)
    axis[1,1].set_title('Hbond Donors', fontsize=7)
    axis[1,1].tick_params(labelsize=5)
    #plt.show()
    figure.savefig('props.png') 

def substr_search(structures):
    matches = Chem.SDWriter('substr_matches.sdf')
    substr = input("Smiles strings of a substructure: ")
    template = Chem.MolFromSmiles(substr) 
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
    ids = df['IDNUMBER'].values.tolist()
    mol = df['ROMol']
    #to_smiles(mol,ids)
    #plot_properties(mol)  
    #substr_search(mol)
    min_str = minim_str(mol[0])
    #confgen_str(mol[ 0 ],10)  # 10 conformers requested
    #make_png(min_str)
 
if __name__ == '__main__':
   main()
