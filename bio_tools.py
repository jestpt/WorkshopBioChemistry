import Bio
from Bio.PDB import *
from Bio.PDB import DSSP
from Bio.PDB.Polypeptide import PPBuilder
from Bio.PDB.Polypeptide import three_to_one
from Bio.PDB.Polypeptide import is_aa
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os



###------------------basic PDB input processor for feature extraction------------------###

def pdb_parser(input_pdb):

    ###Parses pdb from .pdb file
    parser = PDBParser()
    pdb_name = input_pdb[0:-4]
    structure = parser.get_structure(pdb_name, input_pdb)
    return structure, pdb_name

def amino_sequence(input_pdb):
        
    ###Retrieves amino acid sequence from bioython structure object
    structure = pdb_parser(input_pdb)[0]
    chain_sequences = []
    for model in structure:
        for chain in model:
            seq = []
            for residue in chain:
                ## The test below checks if the amino acid
                ## is one of the 20 standard amino acids
                ## Some proteins have "UNK" or "XXX", or other symbols
                ## for missing or unknown residues
                if is_aa(residue.get_resname(), standard=True):
                    seq.append(three_to_one(residue.get_resname()))
                else:
                    seq.append("X")
            chain_sequences.append(str("".join(seq)))
    return  chain_sequences

def amino_acid_by_secondary_structure(sequence, secondary_sequence, secondary_structure_tag):
        
    ###Fetches the amino acid counts by TM:
    all_chains = []
    for chain in sequence:
        PC_dict = {'G': 0, 'A': 0,'V': 0,'L': 0,'M': 0,'I': 0,'F': 0,
               'Y': 0,'W': 0,'S': 0,'T': 0,'C': 0,'P': 0,'N': 0,
               'Q': 0,'K': 0,'R': 0,'H': 0,'D': 0,'E': 0, 'X': 0}
        for residue, residue_TM in zip(sequence, secondary_sequence):
            if residue_TM == secondary_structure_tag:
                PC_dict[residue] = PC_dict[residue] + 1
            else:
                continue
        all_chains.append(PC_dict)
    return all_chains

def b_factor(input_pdb, input_atom = "CA"):

    ###Returns b-factor for the input .pdb atoms, alpha carbon as default
    structure = pdb_parser(input_pdb)[0]
    chain_tagger = []
    for model in structure:
        for chain in model:
            count = 0
            b_factors = {}
            for residue in chain:
                count = count + 1
                for atom in residue:
                    if atom.get_name() == input_atom:
                        B = atom.get_bfactor()
                        b_factors[count] = B
            chain_tagger.append(b_factors)
    return chain_tagger



main = r'C:/Users/A.J. Preto/Desktop/workshop_ANB/biopython'
os.chdir(main)

write_overall_csv("1h2s.pdb")
