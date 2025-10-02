import os
import sys
from rdkit import Chem
from rdkit.Chem import AllChem

def confgen_str(structure, num_conf):
    confs = Chem.SDWriter('minimized_conformers.sdf')
    strH = Chem.AddHs(structure)
    cids = AllChem.EmbedMultipleConfs(strH, numConfs=num_conf)
    print('Num Conformers = ', len(cids))
    rmslist=[]
    AllChem.AlignMolConformers(strH, RMSlist=rmslist)
    print('Num conformers aligned to first str = ',len(rmslist))
#    minimize conformers
    minConf = AllChem.MMFFOptimizeMoleculeConfs(strH)
    print('Num minimized conformers = ',len(minConf)) 
    for i in range(len(minConf)):
        confs.write(strH, i) 
        #print(minConf[ i ])  # Tuple showing convergence (0 = not converged, 1 = converged) and energy

def main():
    inFile = input("Give an input structure file name with .smi or .sdf extension: ")
    nconf = eval(input("Number of conformers required: "))
    strFile = open(inFile, "r")
    strOut = Chem.SDWriter("minimized_conformers.sdf")
# str conversion
    if "smi" in inFile:
       fread = strFile.read()
       mol = Chem.MolFromSmiles(fread)
    if "sdf" in inFile:
       molFile = Chem.SDMolSupplier(inFile)
       mol = molFile[0]
    confgen_str(mol,nconf)

if __name__ == '__main__':
   main()
