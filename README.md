# Spacial-Score
A Comprehensive Topological Indicator for Small Molecule Complexity
(Created by Waldmann Lab at the Max Planck Institute of Molecular Physiology, Dortmund)

The score is intended for assessing molecular topology of organic molecules and improve upon the idea of the fraction of stereo and sp<sup>3</sup> carbons.

The script requires [RDKit package](https://www.rdkit.org/) and numpy.

The script can be used directly from a command line, reading either a directly provided SMILES string (-s) or a .csv/.tsv file (-i):
```
usage: spacial_score_terminal.py [-h] [-s SMILES string] [-i filename.ext] [-o filename.ext] [-t] [-v] [-p]

Script for calculating Spacial Score (SPS) or normalised SPS (nSPS) for small molecules.
The script can calculate the scores for a direct SMILES input or for a .csv or .tsv file containing a list of SMILES.
nSPS is calculated by deafult.

optional arguments:
  -h, --help        show this help message and exit
  -s SMILES string  Your input SMILES string for which to calculate the score
  -i filename.ext   Your .csv or .tsv file containing column called "Smiles" which contains SMILES strings. Resutls will be
                    saved in a new .csv file
  -o filename.ext   You can specify name of the output .csv file. Not required.
  -t                Option to calculate total SPS (no normalisation).
  -v                Option to print verbose results, with information for each atom index.
  -p                Option to print confirmation after procession of each SMILES string in a file.
```

The scores can also be calculated by using function:
```
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
```

