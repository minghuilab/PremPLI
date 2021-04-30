import chimera
from DockPrep import prep
import Midas
import sys
import os

PDB_file = sys.argv[1]
lig_name = sys.argv[2]
num = sys.argv[3]
ligand=chimera.openModels.open("%s_%s_CH%s.pdb"%(PDB_file,lig_name,num))
prep(ligand,addCharges=False)
#prep(ligand,addCharges=True)
Midas.write(ligand,None,PDB_file+"_"+lig_name+"_ADDH_CH" + num +".pdb")
print ("success!")
