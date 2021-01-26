# ligand_superimposition
Small python script for superimposing a ligand into a protein crystal structure.

This script is an quick and easy to use tool for inserting a new ligand over a similar ligand in a crystal structure.  
The packages used are MDAnalysis (1) version 1.0.0 and RDKit (2) version 2020.09.4.
).  

A new ligand needs to be supplied as a RDKit Mol object, which can be generated using SMILEs for example. Hydrogens should be present and the molecule needs to be embedded in 3D. The crystal structure can be supplied by a pdb file and the ligand on which superimposition is to be performed needs to be declared.  
The maximum common substructure is then determined and the Tanimoto distance is returned to give a sense of compound similarity. Using the maximum common substructure, which can be also manually supplied or manipulated to include essential atoms, the rotational matrix of the two substructures is computed. The new ligand is then moved to the center of mass of the old ligand and the rotation matrix is applied. The final position of the new ligand is contained in a MDAnalysis Universe object and can then be merged with the old system, also stored as a Universe object, by selecting all atoms except for the old ligand.  
Before usage in molecular dynamics simulations equalibration with constraints is recommended to prevent atom clashes and subsequent unexpected behaviour.



(1) N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein. MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations. J. Comput. Chem. 32 (2011), 2319–2327. doi:10.1002/jcc.21787  
(2) RDKit: Open-source cheminformatics; http://www.rdkit.org  
