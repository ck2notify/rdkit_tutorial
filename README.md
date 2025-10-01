The script library_process_v2.py can read a chemical library 
(commercial or internal) containing structures in SDF format.
It can convert structures into SMILES format, plot properties
like molecular weight, cLogP, TPSA, H-Bond donor, H-Bond acceptor,
and Number of rotatable bonds (see prop.png file).  These properties 
are also written along with structures (SMILES) in a csv file (test.csv)
so that Datawarrior program can be used to filter structures.  A sample
library of (test.sdf) 24 compounds used as an example.

In addition substructure search can be performed on a library using 
either a SMILES query or an SDF query.  Also chemical structures can 
be minimized and conformations can be generated.  For convenience
one can also create image file in PNG format of structures.  


Hope you may find it useful.

CK
Oct 01, 2025
ck2notify@gmail.com
