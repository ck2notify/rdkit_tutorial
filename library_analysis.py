"""
  discoveryTools.py script can do: 
                          1) Small molecule library analysis
                          2) Minimize one or more structure (small molecule)
                          3) Generate given number of conformers of small molecules
                          4) Substructure search of a fragment from a library
  Please download lib_analysis.py, minimize.py, confgen.py and substr_search.py.
  You need to have python3 and rdkit installed on your device for these scripts to work.  
  This script was tested using python3 (3.13.5), numpy (2.1.3), 
  matplotlib (3.10.0), pandas (2.2.3) and rdkit (2025.03.6) 

  CK, 10/02/2025
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


def lib_properties():
    inFile = input("SDF or SMILES structure file name [with .sdf or .smi extension]: ")
    if 'smi' in inFile:
       smiFile = pd.read_csv(inFile)
       mol = [Chem.MolFromSmiles(smi) for smi in smiFile['Smiles']]
       print(len(mol))
    if 'sdf' in inFile:
       molFile = PandasTools.LoadSDF(inFile)
       mol = molFile['ROMol']
    mwt = [Descriptors.ExactMolWt(structure) for structure in mol]
    clogp = [Descriptors.MolLogP(structure) for structure in mol]
    tpsa = [Descriptors.TPSA(structure) for structure in mol]
    hba = [Descriptors.NumHAcceptors(structure) for structure in mol]
    hbd = [Descriptors.NumHDonors(structure) for structure in mol]
    rotb = [Descriptors.NumRotatableBonds(structure) for structure in mol]
# plotting the data
    figure, axis = plt.subplots(2,3, figsize=(13,7))
    axis[0,0].hist(mwt, bins=5)
    axis[0,0].set_xlabel('MWt', fontsize=7)
    axis[0,0].set_ylabel('Count', fontsize=7)
    axis[0,0].tick_params(labelsize=10)
    axis[0,1].hist(clogp, bins=10)
    axis[0,1].set_xlabel('cLogP', fontsize=7)
    axis[0,1].tick_params(labelsize=10)
    axis[0,2].hist(tpsa, bins=5)
    axis[0,2].set_xlabel('TPSA', fontsize=7)
    axis[0,2].tick_params(labelsize=10)
    axis[1,0].hist(hba, bins=10)
    axis[1,0].tick_params(labelsize=10)
    axis[1,0].set_xlabel('Hbond Acceptors', fontsize=7)
    axis[1,0].set_ylabel('Count', fontsize=7)
    axis[1,1].hist(hbd, bins=10)
    axis[1,1].set_xlabel('Hbond Donors', fontsize=7)
    axis[1,1].tick_params(labelsize=10)
    axis[1,2].hist(rotb, bins=10)
    axis[1,2].set_xlabel('Num. Rotatable Bonds', fontsize=7)
    axis[1,2].tick_params(labelsize=10)
    plt.show()
    figure.savefig('props.png') 

def main():
    lib_properties()

if __name__ == '__main__':
   main()
