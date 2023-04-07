# Script for testing of spacial_score.py

import numpy as np
import rdkit
import spacial_score as sps


valid_smiles = r"C/C=C\C1C=CCC(Br)C1C2=CC(C#C)=CC=C2"
invalid_smiles = r"abcx--a"


def test_smiles_to_mol():
    valid_mol = sps.smiles_to_mol(valid_smiles)
    invalid_mol =  sps.smiles_to_mol(invalid_smiles)

    assert isinstance(valid_mol, rdkit.Chem.rdchem.Mol)
    assert invalid_mol is np.nan


def create_sps_object():
    mol = sps.smiles_to_mol(valid_smiles)
    sps_object = sps.SpacialScore(valid_smiles, mol)
    return sps_object


def calc_sum_score(score_type):
    sum_score = 0
    for atom_score in score_type.values():
        sum_score += atom_score
    return sum_score


def test_hybridisation_score():
    hyb_score = calc_sum_score(create_sps_object().hyb_score)
    assert hyb_score == 40


def test_stereo_score():
    stereo_score = calc_sum_score(create_sps_object().stereo_score)
    assert stereo_score == 23


def test_ring_score():
    ring_score = calc_sum_score(create_sps_object().ring_score)
    assert ring_score == 24


def test_bond_score():
    bond_score = calc_sum_score(create_sps_object().bond_score)
    assert bond_score == 88


def test_SPS_from_smiles():
    sps_score = sps.calculate_score_from_smiles(valid_smiles)
    assert sps_score == 491


def test_nSPS_from_smiles():
    sps_score = sps.calculate_score_from_smiles(valid_smiles, per_atom=True).__round__(2)
    assert sps_score == 27.28
