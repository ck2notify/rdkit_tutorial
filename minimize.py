import os
import sys
from rdkit import Chem
from rdkit.Chem import AllChem

def minim_str(structure):
    strH = Chem.AddHs(structure)
    AllChem.EmbedMolecule(strH)
    AllChem.MMFFOptimizeMolecule(strH)
    return strH

def main():
    inFile = input("Give an input structure file name with .smi or .sdf extension: ")
    strFile = open(inFile, "r")
    strOut = Chem.SDWriter("minimized_str.sdf")
# str conversion
    if "smi" in inFile:
       fread = strFile.read()
       mol = Chem.MolFromSmiles(fread)
    if "sdf" in inFile:
       molFile = Chem.SDMolSupplier(inFile)
       mol = molFile[0]
    minimStr = minim_str(mol)
    strOut.write(minimStr) 

if __name__ == '__main__':
   main()
