# PremPLI-1.0
## About
<font size=4> 
  
Before running PremPLI, you need to create a folder including three input files. 

(see the example of 2021031002563493343426069)
  
</font>

## Three input files in the folder of 2021031002563493343426069
<font size=4> 

1. The 3D structure of a protein (1OSS.pdb1), which can be obtained from the Protein Data Bank (PDB) or created by the user.

2. The SDF file of a ligand (BEN_ideal.sdf), which can be obtained from the Protein Data Bank (PDB) or created by the user.

3. The file containing mutation information (2021031002563493343426069.input), who's name must be consistent with the input folder name.

- PDBfile: coordinate file of a protein structure.
- Partner1: the selected protein chains that will be taken into account during the calculation.
- Partner2: the selected ligand chains that will be taken into account during the calculation.
- MutChain: the protein chain where the mutation occurs.
- Mutation_PDB: the first character is one letter amino acid code for the wild-type residue, the second to penultimate characters indicate residue number in MutChain, and the final character indicates the mutant amino acid.
- Result_Id: a number defined by the user.
- LIG_ID: the identifier of the small molecular.
- LIG_POS: the position in ligand chains.

  The columns are separated by tabs.

</font>
