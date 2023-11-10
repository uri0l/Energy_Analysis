#!/usr/bin/env python
# coding: utf-8

# ## Introduction
# 
# This preparation steps of our protein complex is essential for the subsequent energy analysis. When analyzing a protein-protein interaction, it's crucial to focus solely on the protein components that make up the complex. Removing heteroatoms effectively trims away all non-protein elements, leaving behind only the amino acids that are part of the protein complex itself. This step helps in working with the most natural and biologically relevant representation of the protein complex, so there are no other elements that do not belong, naturally, in the molecule
# 
# - Basic Structure Manipulations: Selecting specific models, chains, and alternative locations within a protein structure.
# 
# - Amide Assignment Correction: Detecting and correcting issues with amide assignments and chirality.
# 
# - Protein Backbone Inspection and Repair: Identifying and fixing problems with the protein backbone, including missing fragments, atoms, and capping.
# 
# - Missing Side-Chain Atom Detection and Correction: Detecting and addressing missing side-chain atoms.
# 
# - Hydrogen Atom Addition: Adding hydrogen atoms to the structure based on specific criteria.
# 
# - Clash Detection: Identifying clashes (steric hindrances) between atoms within the structure.
# 
# - Disulfide Bond Detection: Detecting and classifying possible disulfide bonds in the protein structure.

# ## Installation and importing requirements

# In[2]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')


# ### Basic imports and initialization

# In[3]:


import biobb_structure_checking
import biobb_structure_checking.constants as cts
from biobb_structure_checking.structure_checking import StructureChecking
base_dir_path=biobb_structure_checking.__path__[0]
args = cts.set_defaults(base_dir_path,{'notebook':True})


# ### General Help
# 

# In[4]:


with open(args['commands_help_path']) as help_file:
    print(help_file.read())
#TODO: prepare a specific help method
# print_help(command)


# Set input (PDB or local file, pdb or mmCif formats allowed) and output (local file, pdb format).
# Use pdb:pdbid for downloading structure from PDB (RCSB)

# In[5]:


base_path = "/Users/urileal/ESCI/2nd_year_BDBI/Biophysics/GROUP_PROJECT/"
args['output_format'] = "pdb"
args['keep_canonical'] = False
args['input_structure_path'] = base_path + '6m0j.cif'
args['output_structure_path'] = base_path + '6m0j_fixed.pdb'
args['output_structure_path_charges'] = base_path + '6m0j_fixed.pdbqt'
args['time_limit'] = False
args['nocache'] = False
args['copy_input'] = False
args['build_warnings'] = False
args['debug'] = False
args['verbose'] = False
args['coords_only'] = False
args['overwrite'] = False


# Initializing checking engine, loading structure and showing statistics

# In[7]:


st_c = StructureChecking(base_dir_path, args)


# ## Protein complex structural preparation

# In this series of code blocks with the respective brief explanations, we are preparing the protein complex in order to study it in the most natural form possible. To fully understand the whole series of operations done in the following lines, here's a concise explanation, in a general manner, of how it works:
# 
# - st_c is the object, in our case it reffers to a the inital structure of our protein complex
# 
# Then the code calls various methods on the st_c object, which are the ones that actually perform the assigned operation and return the output.
# 
# - st_c.models()
# - st_c.chains()
# - st_c.altloc()
# - st_c.metals()
# - st_c.ligands()
# - st_c.rem_hydrogen()
# - st_c.water()
# - st_c.amide()
# - st_c.chiral()
# 
# We use arguments like "all", "yes" or "auto" as options within the method call to perform a specific operation.
# 

# ### models 
# Checks afor the presence of models in the structure. MD simulations require a single structure, although some structures (e.g. biounits) may be defined as a series of models, in such case all of them are usually required.
# Use models('--select N') to select model num N for further analysis
# 

# In[19]:


st_c.models()


# ### chains
# Checks for chains (also obtained from print_stats), and allow to select one or more.
# MD simulations are usually performed with complete structures. However input structure may contain several copies of the system, or contains additional chains like peptides or nucleic acids that may be removed. Use chains('X,Y') to select chain(s) X and Y to proceed

# In[20]:


st_c.chains()


# #### altloc
# Checks for the presence of residues with alternative locations. Atoms with alternative coordinates and their occupancy are reported.
# MD simulations requires a single position for each atom.
# Use altloc('occupancy | alt_ids | list of res:id) to select the alternative

# In[21]:


st_c.altloc()


# In[22]:


# We need to choose one of the alternative forms for each residue

st_c.altloc('occupancy')


# In[23]:


st_c.altloc()


# #### metals
# Detects HETATM being metal ions allow to selectively remove them.
# To remove use metals (' All | None | metal_type list | residue list ')

# In[24]:


st_c.metals()


# #### ligands
# Detects HETATM (excluding Water molecules) to selectively remove them.
# To remove use ligands('All | None | Residue List (by id, by num)')

# In[25]:


st_c.ligands()


# In[26]:


st_c.ligands('All')


# In[27]:


st_c.ligands()


# #### remove hydrogens
# Detects and remove hydrogen atoms. MD setup can be done with the original H atoms, however to prevent from non standard labelling, remove them is safer.
# To remove use rem_hydrogen('yes')

# In[28]:


st_c.rem_hydrogen()


# #### water
# Detects water molecules and allows to remove them Crystallographic water molecules may be relevant for keeping the structure, however in most cases only some of them are required. These can be later added using other methods (titration) or manually.
# 
# To remove water molecules use water('yes')

# In[29]:


st_c.water()


# In[30]:


st_c.water("yes")


# #### amide
# Amide terminal atoms in Asn ang Gln residues can be labelled incorrectly.
# amide suggests possible fixes by checking the sourrounding environent.
# 
# To fix use amide ('All | None | residue_list')
# 
# Note that the inversion of amide atoms may trigger additional contacts.
# 

# In[31]:


st_c.amide()


# In[32]:


st_c.amide('all')


# In[33]:


st_c.amide('A42,A103')


# In[34]:


st_c.amide('E394')


# #### chiral
# Side chains of Thr and Ile are chiral, incorrect atom labelling lead to the wrong chirality.
# To fix use chiral('All | None | residue_list')

# In[35]:


st_c.chiral()


# ### Backbone
# Detects and fixes several problems with the backbone use any of --fix_atoms All|None|Residue List --fix_chain All|None|Break list --add_caps All|None|Terms|Breaks|Residue list --no_recheck --no_check_clashes

# In[36]:


st_c.backbone()


# In[37]:


st_c.backbone('--fix_atoms All --fix_chain none --add_caps none')


# #### fixside
# Detects and re-built missing protein side chains.
# To fix use fixside('All | None | residue_list')

# In[38]:


st_c.fixside()


# #### getss
# Detects possible -S-S- bonds based on distance criteria. Proper simulation requires those bonds to be correctly set. Use All|None|residueList to mark them

# In[39]:


st_c.getss()


# In[40]:


st_c.getss('all')


# ### Add_hydrogens
# Add Hydrogen Atoms. Auto: std changes at pH 7.0. His->Hie. pH: set pH value list: Explicit list as [*:]HisXXHid, Interactive[_his]: Prompts for all selectable residues Fixes missing side chain atoms unless --no_fix_side is set Existing hydrogen atoms are removed before adding new ones unless --keep_h set.

# In[41]:


st_c.add_hydrogen()


# In[42]:


st_c.add_hydrogen('auto')


# #### clashes
# Detects steric clashes based on distance criteria.
# Contacts are classified in:
# 
# Severe: Too close atoms, usually indicating superimposed structures or badly modelled regions. Should be fixed.
# Apolar: Vdw colissions.Usually fixed during the simulation.
# Polar and ionic. Usually indicate wrong side chain conformations. Usually fixed during the simulation
# 

# In[43]:


st_c.clashes()


# In[44]:


st_c.checkall()


# ### Save the output file in pdb format
# 

# Specify the location of the output file

# In[55]:


st_c._save_structure(args['output_structure_path'])


# In[46]:


st_c.rem_hydrogen('yes')


# ### Alternative way of calling through the command line

# In[52]:


#st_c.add_hydrogen('--add_charges --add_mode auto')
#Alternative way calling through command line
import os
os.system('check_structure -i ' + args['output_structure_path'] + ' -o ' + args['output_structure_path_charges'] + ' add_hydrogen --add_charges --add_mode auto')


# ### To revert the changes if we want to start again the process

# In[53]:


#st_c._save_structure(args['output_structure_path_charges'])


# In[54]:


#st_c.revert_changes()

