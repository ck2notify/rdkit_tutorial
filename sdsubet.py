"""
  This script can be used to extract a range of structures from 
  a large sdf file.  Results can be stored as test.sdf file.
  CK, 9/30/2025
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

def main():
    df = PandasTools.LoadSDF('/Users/ck/Downloads/Div_set_50k.sdf')
    print(df.head())
    print(df.info())
    mol = df['ROMol']
    nstart, nend = eval(input("Give starting and ending range of structures (eg 1,20): "))
    min_sdf = mol[nstart:nend]
    print(len(min_sdf))
    test = Chem.SDWriter('test.sdf')
    for structure in min_sdf:
        test.write(structure)

if __name__ == '__main__':
   main()
