"""
  RDKit tool to convert a library of SDF structures to SMILES strings
  CK, 8/25/25
"""

import os
import sys
from rdkit import Chem


infile = input("Give Input SDF file name (with .sdf extension): ")
oufile = input("Give Output SMILES file name (with .smi extension): ")

lib_sdf = Chem.SDMolSupplier(infile)
smi_file = open(oufile, 'w')  
#print("Smiles,Title", file=smi_file)  # you can open using Datawarrior

def main():
    for mol in lib_sdf:
        x = Chem.MolToSmiles(mol)  # SDF to SMILES
        t = mol.GetProp('_Name') # Get Vendor ID (it may change with Vendor
        smi_file.write("{0},{1}\n".format(x,t)) # csv style
       # for kekulization of smiles strings
       # y = Chem.MolToSmiles(mol, kekuleSmiles=True)
       # smi_file.write("{0},{1}\n".format(y,t))

if __name__ == '__main__':
   main()
