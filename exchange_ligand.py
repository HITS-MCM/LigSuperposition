from rdkit import Chem
from rdkit.Chem import rdFMCS
from rdkit import DataStructs
from rdkit.Chem import AllChem
import MDAnalysis as mda
from MDAnalysis.analysis import align
import os


class Superimpose:
    """Used to superimpose a new ligand into a given crystal structure.

    Old System is loaded from specified file. New ligand needs to be supplied as rdkit mol object with Hydrogens and 3D-
    coordinates. The maximum match between the two ligands is found and used to extract corresponding coordinates of
    the ligands. Atom indices can be manually added or can be completely supplied manually. Rotation matrix is computed
    and used to rotate the centered new ligand. The new ligand is translated by the center of mass of the old ligand.
    This results in a mda Universe containing the superimposed new ligand.

    """
    def __init__(self, topology, ligand_old, ligand_new, data_path="./data/", coordinates=None):
        """Inits superimpose class.

        Args:
            system (str): File name of the file containing the equilibrated old coordinate file. Any file that can be
                          read by MDAnalysis can be used.
            ligand_old (str): Name of the ligand so that it can be extracted.
            ligand_new (rdkit.Chem.rdchem.Mol): New ligand with added hydrogens and 3D structure.
            data_path (str): Path to where the data is stored.
        """
        self.path = data_path
        if coordinates is None:
            self.system = mda.Universe(self.path + topology)
        else:
            self.system = mda.Universe(self.path+topology, self.path+coordinates)
        self.ligand_new = ligand_new
        self.ligand_old = self.system.select_atoms("resname "+ligand_old)
        self.max_common_struc = None
        self.old_match = None
        self.new_match = None

    def find_similar(self):
        """Finds the maximum common substructure.

        Returns:
            Nothing.
        """
        #os.chdir(self.path)
        sim_path = os.path.join(self.path, "similar_temp")
        old_path = os.path.join(sim_path, "old_temp.pdb")
        if not os.path.isdir(sim_path):
            os.mkdir(sim_path, mode=0o777)
        self.ligand_old.write(old_path)
        ligand_old = Chem.MolFromPDBFile(old_path, removeHs=False)
        ligand_old_fp = Chem.RDKFingerprint(ligand_old)
        ligand_new_fp = Chem.RDKFingerprint(self.ligand_new)
        similarity = DataStructs.FingerprintSimilarity(ligand_old_fp, ligand_new_fp)
        print("Calculated tanimoto similarity between ligands:   ", similarity)
        max_common_struc = rdFMCS.FindMCS([ligand_old, self.ligand_new]).smartsString
        print("Maximum common substructure found to be: ", max_common_struc)
        self.max_common_struc = Chem.MolFromSmarts(max_common_struc)

    def get_match(self, add_index_old=None, add_index_new=None, overwrite_index_old=None, overwrite_index_new=None):
        """Creates mol object of matching atoms between ligands and their coordinates.

        Note that the indices of the ligands need to have the same shape.

        Args:
            add_index_old (list): Indices that should be additionally included in the old ligand.
            add_index_new (list): Indices that should be additionally included in the new ligand.
            overwrite_index_old (list): Overwrites automatically found indices.
            overwrite_index_new (list): Overwrites automatically found indices.

        Returns:
            Nothing.

        """
        ligand_old = Chem.MolFromPDBFile(os.path.join(self.path, "similar_temp", "old_temp.pdb"), removeHs=False)
        old_indx = list(ligand_old.GetSubstructMatch(self.max_common_struc))
        new_indx = list(self.ligand_new.GetSubstructMatch(self.max_common_struc))
        if add_index_old is not None:
            old_indx = old_indx + add_index_old
        if add_index_new is not None:
            new_indx = new_indx + add_index_new
        if overwrite_index_old is not None:
            old_indx = overwrite_index_old
        if overwrite_index_new is not None:
            new_indx = overwrite_index_new
        old = (old_indx, ligand_old, "old")
        new = (new_indx, self.ligand_new, "new")
        for molecule in [old, new]:
            index, mol, kind = molecule
            match = Chem.EditableMol(Chem.Mol())
            conformer = Chem.rdchem.Conformer()
            i = 0
            for entry in index:
                pos = mol.GetConformer().GetAtomPosition(entry)
                match.AddAtom(mol.GetAtoms()[entry])
                conformer.SetAtomPosition(i, pos)
                i += 1
            if kind == "old":
                self.old_match = match.GetMol()
                self.old_match.AddConformer(conformer)
            elif kind == "new":
                self.new_match = match.GetMol()
                self.new_match.AddConformer(conformer)

    def rotation_matrix(self):
        if not os.path.isdir(os.path.join(self.path, "rotation_temp")):
            os.mkdir(os.path.join(self.path, "rotation_temp"), mode=0o777)
        path = os.path.join(self.path, "rotation_temp/")
        Chem.MolToPDBFile(self.old_match, path + "old.pdb")
        Chem.MolToPDBFile(self.new_match, path + "new.pdb")
        old = mda.Universe(path + "old.pdb")
        new = mda.Universe(path + "new.pdb")
        weights = (old.atoms.masses + new.atoms.masses)/2
        old_pos = old.select_atoms("all").positions - old.atoms.center_of_mass()
        new_pos = new.select_atoms("all").positions - new.atoms.center_of_mass()
        rot_matrix, rmsd = align.rotation_matrix(new_pos, old_pos, weights)
        return rot_matrix, rmsd

    def rotate_translate_ligand(self, lig_name="NEW"):
        if not os.path.isdir(os.path.join(self.path, "rotation_temp")):
            os.mkdir(os.path.join(self.path, "rotation_temp"), mode=0o777)
        path = os.path.join(self.path, "rotation_temp/")
        rot_matr, rmsd = self.rotation_matrix()
        Chem.MolToPDBFile(self.ligand_new, path + "new_lig.pdb")
        new_lig = mda.Universe(path + "new_lig.pdb")
        new = mda.Universe(path + "new.pdb")
        old = mda.Universe(path + "old.pdb")
        new_lig.atoms.translate(-new.atoms.center_of_mass())
        new_lig.atoms.rotate(rot_matr)
        new_lig.atoms.translate(old.atoms.center_of_mass())
        new_lig.atoms.residues.resnames = lig_name
        new_lig.atoms.record_types = ["HETATM"] * len(new_lig.atoms)
        self.ligand_new = new_lig

    def remove_temp_data(self):
        """Removes all generated temporary data.

        Returns:
            Nothing.
        """
        rot_path = os.path.join(self.path, "rotation_temp")
        rot_files = ["old.pdb", "new.pdb", "new_lig.pdb"]
        sim_files = ["old_temp.pdb"]
        sim_path = os.path.join(self.path, "similar_temp")
        try:
            for file in rot_files:
                os.remove(os.path.join(rot_path, file))
            for file in sim_files:
                os.remove(os.path.join(sim_path, file))
            os.rmdir(rot_path)
            os.rmdir(sim_path)
        except OSError:
            print("There was an error in deleting all temporary files. Please delete them manually.")


if __name__ == "__main__":
    lig_new = Chem.MolFromSmiles("C1CCN(C1)CC#CCN2CCCC2=O")
    lig_new = Chem.AddHs(lig_new)
    AllChem.EmbedMolecule(lig_new)
    superimp = Superimpose(topology="struct.pdb", ligand_old="IXO", ligand_new=lig_new)
    superimp.find_similar()
    superimp.get_match(add_index_old=[6], add_index_new=[14])
    superimp.rotate_translate_ligand("OXO")
    superimp.remove_temp_data()
    output = mda.Merge(superimp.system.atoms, superimp.ligand_new.atoms)
    output.atoms.write("output.pdb")
