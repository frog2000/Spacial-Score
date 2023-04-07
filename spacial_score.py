# Version 1.0

from rdkit import Chem
import rdkit.Chem.Descriptors as Desc
import numpy as np
import argparse
import sys
import csv


class SpacialScore:
    """Class intended for calculating spacial score (SPS) and size-normalised SPS (nSPS) for small organic molecules"""
    def __init__(self, smiles, mol, verbose=False):
        self.smiles = smiles
        self.mol = mol
        self.verbose = verbose

        self.hyb_score = {}
        self.stereo_score = {}
        self.ring_score = {}
        self.bond_score = {}
        self.chiral_idxs = self.find_stereo_atom_idxs()
        self.doublebonds_stereo = self.find_doublebonds_stereo()
        self.score = self.calculate_spacial_score()
        self.per_atom_score = self.score/Desc.HeavyAtomCount(self.mol)

        if self.verbose:
            self.display_scores()


    def display_scores(self):
        """Displays the individual scores for each molecule atom"""

        print("SMILES:", self.smiles)
        print("Atom Idx".ljust(10, " "), end="")
        print("Element".ljust(10, " "), end="")
        print("Hybrid".ljust(10, " "), end="")
        print("Stereo".ljust(10, " "), end="")
        print("Ring".ljust(10, " "), end="")
        print("Neighbs".ljust(10, " "))
        print("".ljust(60, "-"))

        for atom in self.mol.GetAtoms():
            atom_idx = atom.GetIdx()
            print(str(atom_idx).ljust(10, " "), end="")
            print(str(Chem.rdchem.Atom.GetSymbol(atom)).ljust(10, " "), end="")
            print(str(self.hyb_score[atom_idx]).ljust(10, " "), end="")
            print(str(self.stereo_score[atom_idx]).ljust(10, " "), end="")
            print(str(self.ring_score[atom_idx]).ljust(10, " "), end="")
            print(str(self.bond_score[atom_idx]).ljust(10, " "))
       
        print("".ljust(60, "-"))
        print("Total Spacial Score:", self.score)
        print("Per-Atom Score:", self.per_atom_score.__round__(2), "\n")


    def find_stereo_atom_idxs(self, includeUnassigned=True):
        """Finds indeces of atoms that are (pseudo)stereo/chiralcentres, in repsect to the attached groups (does not account for double bond isomers)"""
        stereo_centers = Chem.FindMolChiralCenters(self.mol, includeUnassigned=includeUnassigned, includeCIP=False, useLegacyImplementation=False)
        stereo_idxs = [atom_idx for atom_idx, _ in stereo_centers]
        return stereo_idxs
    

    def find_doublebonds_stereo(self):
        """Finds indeces of stereo double bond atoms (E/Z)"""
        db_stereo = {}
        for bond in self.mol.GetBonds():
            if str(bond.GetBondType()) == "DOUBLE":
                db_stereo[(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())] = bond.GetStereo()
        return db_stereo


    def calculate_spacial_score(self):
        """Calculates the total spacial score for a molecule"""
        score = 0
        for atom in self.mol.GetAtoms():
            atom_idx = atom.GetIdx()
            self.hyb_score[atom_idx] = self._account_for_hybridisation(atom)
            self.stereo_score[atom_idx] = self._account_for_stereo(atom_idx)
            self.ring_score[atom_idx] = self._account_for_ring(atom)
            self.bond_score[atom_idx] = self._account_for_neighbours(atom)
            score += self._calculate_score_for_atom(atom_idx)
        return score


    def _calculate_score_for_atom(self, atom_idx):
        """Calculates the total score for a single atom in a molecule"""
        atom_score = self.hyb_score[atom_idx] * self.stereo_score[atom_idx] * self.ring_score[atom_idx] * self.bond_score[atom_idx]
        return atom_score


    def _account_for_hybridisation(self, atom):
        """Calculates the hybridisation score for a single atom in a molecule"""
        hybridisations = {"SP": 1, "SP2": 2, "SP3": 3}
        hyb_type = str(atom.GetHybridization())

        if hyb_type in hybridisations.keys():
            return hybridisations[hyb_type]
        return 4  # h score for any other hybridisation than sp, sp2 or sp3 
    

    def _account_for_stereo(self, atom_idx):
        """Calculates the stereo score for a single atom in a molecule"""
        if atom_idx in self.chiral_idxs:
            return 2
        for bond_atom_idxs, stereo in self.doublebonds_stereo.items():
            if atom_idx in bond_atom_idxs and not(str(stereo).endswith("NONE")):
                return 2
        return 1


    def _account_for_ring(self, atom):
        """Calculates the ring score for a single atom in a molecule"""
        if atom.GetIsAromatic():  # aromatic rings are not promoted 
            return 1
        if atom.IsInRing():
            return 2
        return 1


    def _account_for_neighbours(self, atom):
        """Calculates the neighbour score for a single atom in a molecule
        The second power allows to account for branching in the molecular structure"""
        return (len(atom.GetNeighbors()))**2


def smiles_to_mol(smiles: str):
    """ Generate a RDKit Molecule from SMILES.

    Parameters:
    ===========
    smiles: the input string

    Returns:
    ========
    The RDKit Molecule. If the Smiles parsing failed, NAN is returned instead.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            return mol
        return np.nan
    except:
        return np.nan


def close_files(open_files:tuple):
    """Closes open files"""
    for file in open_files:
        file.close()


def process_input(smiles:str, filename:str, output_name:str, total_score:bool, verbose:bool, confirmation:False):
    """Processes the command line input to print out the resulting score or create a file with added results"""

    if smiles:  # process a directly provided SMILES string
        result = calculate_score_from_smiles(smiles, per_atom=(not total_score), verbose=verbose)
        if not verbose:
            score_type = "SPS" if total_score else "nSPS" 
            print(f"\nNormalisation Applied: {not total_score}\nSMILES: {smiles}\nCalculated {score_type}: {result}")
        if result is np.nan:
            print("\nPlease double-check your input SMILES string...")
    
    elif filename:  # process provided file
        provided_filename_base = filename.split(".")[0]
        output_filename = output_name if output_name else provided_filename_base + "_SPS.csv"
        outfile = open(output_filename, "w")
        
        # read the input .csv or .tsv file
        if filename.endswith("csv"):
            infile = open(filename, "r")
            reader = csv.DictReader(infile, dialect="excel")
        elif filename.endswith("tsv"):
            infile = open(filename, "r")
            reader = csv.DictReader(infile, dialect="excel-tab")
        else:
            raise ValueError(f"Unknown input file format: {filename}")

        print("\nProcessing, please wait...")
        for idx, row in enumerate(reader):
            if idx == 0:
                header = [column_name for column_name in row]  # read existing headers
                # add SPS or nSPS column to the file 
                if total_score:
                    header.append("SPS")  
                header.append("nSPS")
                outfile.write(",".join(header) + "\n")
            
            try:
                if total_score:
                    row["SPS"] = calculate_score_from_smiles(row["Smiles"], per_atom=False, verbose=verbose)
                row["nSPS"] = calculate_score_from_smiles(row["Smiles"], per_atom=True, verbose=verbose)
            except KeyError:
                close_files((outfile, infile))
                raise KeyError("Please make sure that your file contains column called 'Smiles' with SMILES strings")
            
            line = [str(row[x]) for x in row]  # reconstruct the row 
            outfile.write(",".join(line) + "\n")
            if confirmation:
                print("Finished calculations for:", row["Smiles"])

        close_files((outfile, infile))
        print(f"Finished. {output_filename} was saved.")
    else:
        raise ValueError(f"No input was provided")


def calculate_score_from_smiles(smiles: str, per_atom=False, verbose=False) -> float:
    """ Calculates the spacial score as a total SPS or size-normalised, per-atom nSPS for a molecule.

    Parameters:
    ===========
    smiles: valid SMILES string
    per_atom: flag to denote if the normalised per-atom result (nSPS) should be returned
    verbose: flag to denote if the detailed scores for each atom should be printed

    Returns:
    ========
    Total or per-atom numeric spacial score for the provided molecule.
    """
    mol = smiles_to_mol(smiles)
    if mol is np.nan:
        return np.nan
    sps = SpacialScore(smiles, mol, verbose)
    if per_atom:
        return sps.per_atom_score
    return sps.score


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=
    'Script for calculating Spacial Score (SPS) or normalised SPS (nSPS) for small molecules.\
    \nThe script can calculate the scores for a direct SMILES input or for a .csv or .tsv file containing a list of SMILES.\
    \nnSPS is calculated by deafult.', 
    usage=None, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-s', action="store",
                        metavar="SMILES string",
                        help='Your input SMILES string for which to calculate the score', default=None)
    parser.add_argument('-i', action="store",
                        help='Your .csv or .tsv file containing column called "Smiles" which contains SMILES strings. Resutls will be saved in a new .csv file', 
                        metavar='filename.ext', 
                        default=None)
    parser.add_argument('-o', action="store",
                        help='You can specify name of the output .csv file. Not required.', 
                        metavar='filename.csv', 
                        default=None)
    parser.add_argument('-t', action="store_true",
                        help='Option to calculate total SPS (no normalisation).',
                        default=False)
    parser.add_argument('-v', action="store_true",
                        help='Option to print verbose results, with information for each atom index.', 
                        default=False)
    parser.add_argument('-p', action="store_true",
                        help='Option to print confirmation after processing of each SMILES string in a file.',
                        default=False)
                                            
    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(1)

    ARGS = parser.parse_args()
    process_input(smiles=ARGS.s, filename=ARGS.i, output_name=ARGS.o, total_score=ARGS.t, verbose=ARGS.v, confirmation=ARGS.p) 
