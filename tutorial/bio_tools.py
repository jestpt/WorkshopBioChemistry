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

print amino_sequence("1h2s.pdb")
print b_factor("1h2s.pdb")