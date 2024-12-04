from Bio.PDB.PDBParser import PDBParser 	## Module to parse the structure of the protein
from Bio.PDB.NeighborSearch import NeighborSearch	## Module to search atoms within a distance

def get_interface_residues(st, id1, id2, distance):
	interface_chain1 = set()
	interface_chain2 = set()

	chain1 = st[0][id1]	## Get the desired chains from the protein
	chain2 = st[0][id2]

	NeighborSearch_chain2 = NeighborSearch(list(chain2.get_atoms())) 	## Prepare the Neighbour search for chain E

	for res_chain1 in chain1:
		for atom_chain1 in res_chain1:	## Iterate over all the atoms from chain A
		
			Neighbor_atom_chain2 = NeighborSearch_chain2.search(atom_chain1.coord, distance)	## Look for atoms in chain E within the distance from the atom in chain A
			
			for atom_chain2 in Neighbor_atom_chain2:	## Itarate over the atoms we know they are within the distance
			
				res_chain2 = atom_chain2.get_parent()		## Get the residues to which the atom belongs to
				
				interface_chain1.add(res_chain1)		## Use a set to only save each residue once
				interface_chain2.add(res_chain2)		

	return list(interface_chain1), list(interface_chain2)	## Return the list for each chain

pdb_path = "/Users/urileal/ESCI/2nd_year_BDBI/Biophysics/GROUP_PROJECT/6m0j_fixed.pdb"	## Path to the proteins's structure file
parser = PDBParser()
st = parser.get_structure('st', pdb_path)	# Get the model of the protein

interface_residues_chainA, interface_residues_chainE = get_interface_residues(st, 'A', 'E', 4)	## Call the function with the desired arguments

print(f"Interface residues in chain A: {interface_residues_chainA}") 	## List of residues in chain A within the distance from chain E
print(f"Interface residues in chain E: {interface_residues_chainE}")		## List of residues in chain E within the distance from chain A
