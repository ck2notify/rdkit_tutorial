"""
  First RDKit code to write SMILES string into SDF file
  CK, 8/25/25
"""

import os
import sys
from rdkit import Chem

infile = input("Give Input SMILES file Name (with.smi extension): ")
oufile = input("Give Output SDF file Name (with .sdf extension)': ")

smiles = open(infile, 'r')  # for writing sdf structure as a block
sdfs = Chem.SDWriter(oufile)  # An SDWriter object is created

def main():
    count = 0
    for smi in smiles:
        if count == 10:
           break
        else:
           str_split = smi.split(',')
           sdf = Chem.MolFromSmiles(str_split[0])   # converting SMILES to 2d SDF
           title = str_split[1].strip('\n')
           sdf.SetProp("_Name", title)       # adding to Title to SDF file block
           sdfs.write(sdf)  # this writes with $$$$ with proper title
        count += 1

if __name__ == '__main__':
   main() 
