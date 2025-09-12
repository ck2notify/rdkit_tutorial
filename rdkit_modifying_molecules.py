import os
import sys
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem

# RDKit normally treats hydrogens implicitly.  However, 
# if they need to be explicitly present for optimizing 
# 3D geometry, we can use following method to addHs.

def make_png(structure):
    img = Draw.MolToImage(structure)
    img.save('m.png')

m = Chem.MolFromSmiles('CCO')
print(m.GetNumAtoms())

mH = Chem.AddHs(m)
print(mH.GetNumAtoms())

# Hydrogens can be removed using RemoveHs
mNoH = Chem.RemoveHs(mH)
print(mNoH.GetNumAtoms())

# RDKit stores bonds in aromatic rings having aromatic bond types
# can be changed with kekulize (only bond type aromatic is stored)
# can be converted to double bond

mAr = Chem.MolFromSmiles('c1ccccc1')
print(mAr.GetBondWithIdx(0).GetBondType())
Chem.Kekulize(mAr)
print(mAr.GetBondWithIdx(0).GetBondType())
print(mAr.GetBondWithIdx(1).GetBondType())

# However, bonds are still marked as aromatic, though it is still 
# identified as double or single bond.

print(mAr.GetBondWithIdx(1).GetIsAromatic())

# Bond type aromatic flag can be changed using the following steps

mAr = Chem.MolFromSmiles('c1ccccc1')
print(mAr.GetBondWithIdx(0).GetIsAromatic())

Chem.Kekulize(mAr, clearAromaticFlags=True)
print(mAr.GetBondWithIdx(0).GetIsAromatic())

# Bonds can be restored to the aromatic bond type using 
# rdkit.Chem.rdmolops.SanitizeMol() function

Chem.SanitizeMol(mAr)  
print(mAr.GetBondWithIdx(0).GetBondType())

# Generating 2D coordinates and aligning structures to a 
# common scaffold example

template = Chem.MolFromSmiles('c1nccc2n1ccc2')
AllChem.Compute2DCoords(template) # generating 2D coords of template scaffold

# Align these structures to the template

m1 = Chem.MolFromSmiles('OCCc1ccn2cnccc12')
m2 = Chem.MolFromSmiles('C1CC1Oc1cc2ccncn2c1')
m3 = Chem.MolFromSmiles('CNC(=O)c1nccc2cccn12')

AllChem.GenerateDepictionMatching2DStructure(m1, template)
AllChem.GenerateDepictionMatching2DStructure(m2, template)
AllChem.GenerateDepictionMatching2DStructure(m3, template)

f_align = Chem.SDWriter('aligned.sdf') # I used chimerax to view these strs
f_align.write(template)
f_align.write(m1)
f_align.write(m2)
f_align.write(m3)


# minmization MMFF94


f = Chem.SDWriter('m2.sdf')
m = Chem.MolFromSmiles('C1CCC1OC')
m2 = Chem.AddHs(m)
#f.write(m2)  # 2D str
AllChem.EmbedMolecule(m2)  # Embed step assigns 3d coordinates to each atom
                           # using Cambridge Structural Database (CDD)
                           # torsion angle preferences.
  
#f.write(m2)  # 3D str, but not minimized

AllChem.MMFFOptimizeMolecule(m2)  #MMFF94 is the default minimization FF after RDKit 2024.03
#f.write(m2)  # Minimized str
#make_png(m2)


# Generating multiple conformers
"""
To generate multiple conformers, call embed multiple times to generate 
multiple starting points.  The option 'numConfs' allows to set number 
of conformers to be generated.  
"""

m = Chem.MolFromSmiles('N(C)(C)S(=O)(=O)c1ccc(S(=O)(=O)Nc2c(N3CCCCC3)cccc2)cc1')
m2 = Chem.AddHs(m)
# run ETKDG 10 times (generate 10 set of starting coordinates)
cids = AllChem.EmbedMultipleConfs(m2, numConfs=25)
print(len(cids)) # check how many generated
# To align them all together
rmslist = []
AllChem.AlignMolConformers(m2, RMSlist=rmslist)  # all conformers are aligned to first str
print(len(rmslist))  # you'll see one less than 10

# rmslist contains rms values of all 9 conformers to first str
# if we want rms value of any conformers use this (1 & 9)
rms_19 = AllChem.GetConformerRMS(m2, 1, 9, prealigned=True)
rms_09 = AllChem.GetConformerRMS(m2, 0, 9, prealigned=True)
print(rmslist)
print('RMS_1and_9 = {:0.4f}'.format(rms_19))
print('RMS_0and_9 = {:0.4f}'.format(rms_09))

# To minimize all conformers (although not needed) use this:
res = AllChem.MMFFOptimizeMoleculeConfs(m2)
# res is a 2-tuples (not_converged, energy) for each conformer
# if not_converged = 0, minimization for that conformer is converged

print(type(res))
print(len(res))
for i in range(len(res)):
    f.write(m2,i)
    print(res[i])

"""
By default AllChem.EmbedMultipleConfs and AllChem.MMFFOptimizeMoleculeConfs()
run single thread, but this can be changed to use multiple threads  by this:

params = AllChem.ETKDGv3()
params.numThreads = 0
cids = AllChem.EmbedMultipleConfs(m2, 25, params)
res = AllChem.MMFFOptimizeMoleculeConfs(m2, numThreads=0)

Setting numThreads to 0 causes the software to use maximum number of threads
available 

Although ETKDG conformer generation is useful for most purposes, there is a 
separate conformer analysis tool exists.
"""
