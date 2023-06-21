import tarfile
import gzip
import io
from Bio.PDB import PDBParser

def extract_pdb_sequences():
    # Create a PDBParser object
    parser = PDBParser()
    # Open the tar archive containing the PDB files
    with tarfile.open("UP000005640_9606_HUMAN_v4.tar", "r") as tar:
        # Iterate over the file names in the tar archive
        for member in tar.getmembers():
            # Check if the file is a PDB file and model is F1 model
            if member.name.endswith(".pdb.gz") and member.name.split("-")[2] == "F1":
                # Extract the protein ID from the file name
                protein_id = member.name.split("-")[1]
                # Extract the contents of the PDB file from the tar archive
                contents = tar.extractfile(member).read()
                # Decompress the contents of the PDB file using the gzip module
                contents = gzip.decompress(contents)
                # Parse the PDB file using Biopython's PDBParser
                structure = parser.get_structure(member.name, io.StringIO(contents.decode('utf-8')))
                # conversion between each amino acid 3 letter representation to 1 letter representation
                d3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
                         'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
                         'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
                         'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
                # Extract the FASTA sequence from the structure
                fasta_sequence = ""
                for chain in structure[0]:
                    for residue in chain:
                        if residue.get_resname() in (
                                "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS",
                                "MET",
                                "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"):
                            fasta_sequence += d3to1[residue.get_resname()]

                # write fasta sequence of protein to fasta sequences file
                with open("sequences_without_repeats.txt", "a") as sequences_file:
                    sequences_file.write(f">{protein_id}\n{fasta_sequence}\n")
