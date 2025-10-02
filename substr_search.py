import os
import sys
from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd
from rdkit.Chem import PandasTools

def substrSearch(structures,template):
    matches = Chem.SDWriter('substr_matched.sdf')
    substr_match = []
    for structure in structures:
        if structure.HasSubstructMatch(template):
           substr_match.append(structure)
           matches.write(structure)
    print("Number of matches found: ",len(substr_match))

def main():
    inFile = input("Library file name [with .sdf or .smi extension]: ")
    queryFile = input("Substructure file name [with .sdf or .smi extension]: ")
    queryStr = open(queryFile, 'r')

    if "smi" in queryFile:
       fread = queryStr.read()
       substr = Chem.MolFromSmiles(fread)
    if "sdf" in queryFile:
       molFile = Chem.SDMolSupplier(queryFile)
       substr = molFile[0]
    
    if 'smi' in inFile:
       smiFile = pd.read_csv(inFile)
       mol = [Chem.MolFromSmiles(smi) for smi in smiFile['Smiles']]
       print('Number of structures in library:', len(mol))
    if 'sdf' in inFile:
       molFile = PandasTools.LoadSDF(inFile)
       mol = molFile['ROMol']

    substrSearch(mol, substr)

if __name__ == '__main__':
   main()
