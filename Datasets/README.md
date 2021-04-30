# Datasets

## About

<font size=4>

S796.txt: The training dataset for parameterizing PremPLI model and it contains 796 single point mutations.

S144.txt: The independent testing dataset, which contains 144 single point mutations.

S129.txt: The independent testing dataset, which contains 129 single point mutations.

S99.txt: The subset of dataset S796, which contains 99 single point mutations.

</font> 

## Terms in the text files

<font size=4>

PDB ID: The PDB entry of the protein-ligand complex.

Ligand ID: The small molecule identifier in the protein-ligand complex.

TKI: The tyrosine kinase inhibitors.

Mutated Chain: Protein chain with mutation.

Mutation_PDB: The mutation corresponding to the residue numbering found in the protein databank. The first character is one letter amino acid code for the wild-type residue, the second to penultimate characters indicate residue number, and the final character indicates the mutant amino acid.

MUTATION: The mutation corresponding to the residue numbering in the 'cleaned' pdb files from the study of Hauser et al.[PMID: 30159405].

DDGexp: Experimental changes of binding affinity upon a single mutation (in kcal/mol).

Interface: Whether the mutation occurs at the protein-ligand binding interface (yes/no). The interface residues were defined if any heavy atoms of proteins are within 5 Å distance from any heavy atoms of ligands.

PremPLI: Predicted binding affinity change (in kcal/mol), and positive and negative signs correspond to the mutations decreasing and increasing binding affinity respectively.

PremPLI(CV3): "Leave-one-complex-out" cross-validation results (in kcal/mol).

PremPLI_C: PremPLI was retrained after removing all mutations in the overlapped complexes with S99 from the training dataset.

mCSM-lig: The binding affinity changes predicted by mCSM-lig, taken from [PMID: 27384129].

ML1: The ML1 model trained on 484 single mutations from the Platinum database, taken from [PMID: 31482130]. 

Prime: MD simulations combined with the solution of the Generalized Born equation calculated by Prime, taken from [PMID: 30159405].

FEP+: Alchemical free-energy perturbation calculations using FEP+, taken from [PMID: 30159405].

R15: Rosetta using the flex_ddg protocol and REF2015 scoring function, taken from [PMID: 31482130].

R14: Rosetta using the flex_ddg protocol and talaris2014 scoring function, taken from [PMID: 30648154].

A14: Amber14sb and GAFF(v2.1)/AM1-BCC force fields were used for proteins and ligands respectively, taken from [PMID: 30648154].

RMD: The combination of R14 and A14, taken from [PMID: 30648154].

</font> 
