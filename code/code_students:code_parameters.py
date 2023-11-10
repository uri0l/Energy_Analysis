#!/usr/bin/env python
# coding: utf-8

# In[1]:


import argparse
import sys
import os
import math

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.NACCESS import NACCESS_atomic
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.PDBIO import PDBIO, Select


# This are functions that you will need to import the parameters for VanderWaals or the residue library:

# In[2]:


class ResiduesDataLib():
    def __init__(self, fname):
        self.residue_data = {}
        try:
            fh = open(fname, "r")
        except OSError:
            print("#ERROR while loading library file (", fname, ")")
            sys.exit(2)
        for line in fh:
            if line[0] == '#':
                continue
            data = line.split()
            r = Residue(data)
            self.residue_data[r.id] = r
        self.nres = len(self.residue_data)

    def get_params(self, resid, atid):
        atom_id = resid + ':' + atid
        if atom_id in self.residue_data:
            return self.residue_data[atom_id]
        else:
            print("WARNING: atom not found in library (", atom_id, ')')
            return None

class Residue():
    def __init__(self,data):
        self.id     = data[0]+':'+data[1]
        self.at_type = data[2]
        self.charge  = float(data[3])
        
class AtomType():
    def __init__(self, data):
        self.id   = data[0]
        self.eps  = float(data[1])
        self.sig  = float(data[2])
        self.mass = float(data[3])
        self.fsrf = float(data[4])
        self.rvdw = self.sig * 0.5612
        
class VdwParamset(): #extracted from GELPI's github
    #parameters for the VdW
    def __init__ (self, file_name):
        self.at_types = {}
        try:
            fh = open(file_name, "r")
        except OSError:
            print ("#ERROR while loading parameter file (", file_name, ")")
            sys.exit(2)
        for line in fh:
            if line[0] == '#':
                continue
            data = line.split()
            self.at_types[data[0]] = AtomType(data)
        self.ntypes = len(self.at_types)
        fh.close()


# We need to load the residue library. 

# In[3]:


# loading residue library from data/aaLib.lib
residue_library = ResiduesDataLib('/Users/urileal/ESCI/2nd_year_BDBI/Biophysics/GROUP_PROJECT/aalLib.lib')


# We also need to load the Vdw parameters.

# In[ ]:


# loading VdW parameters
ff_params = VdwParamset('/Users/urileal/ESCI/2nd_year_BDBI/Biophysics/GROUP_PROJECT/VdW_parameters.txt')


# Now, let's load the strucutre of the PDB

# In[ ]:


# set the pdb_path and load the structure
pdb_path = "/Users/urileal/ESCI/2nd_year_BDBI/Biophysics/GROUP_PROJECT/6m0j.pdb"
# Setting the Bio.PDB.Parser object
parser = PDBParser(PERMISSIVE=1)
# Loading structure
st = parser.get_structure('st', pdb_path)

