# Final Project – Computational Biology
Authors: Jonathan Mandl and Dina Englard

## Table of Contents
1. [Extract_pdb_sequences](#extract_pdb_sequences)
2. [Blast_cmd.txt](#blast_cmdtxt)
3. [Parse_blast](#parse_blast)
4. [PDB_Locations](#pdb_locations)
5. [Extract_Locations](#extract_locations)
6. [Find_triplets](#find_triplets)

---

## Extract_pdb_sequences <a name="extract_pdb_sequences"></a>
- This code extracts the PDB files from a TAR archive, decompresses them, parses them using Biopython's PDBParser, converts the amino acid residues to their one-letter representation, and writes the protein ID and FASTA sequence to a file.

## Blast_cmd.txt <a name="blast_cmdtxt"></a>
- These commands, to be run in cmd, create a protein sequence database using "makeblastdb" and then perform a protein sequence comparison using "blastp" between human protein sequences and the sequences in the Pan Troglodytes NR (non-redundant) database which we downloaded from the NCBI website. The top five matches are saved in the "nr_blast_output.txt" file.

## Parse_blast <a name="parse_blast"></a>
- This code reads the BLAST output file which contains the best Pan Troglodytes proteins alignments for each human protein sequence. The code parses the file using NCBIStandalone.BlastParser, and iterates over the blast records and their alignments. It extracts relevant information such as query and hit IDs, alignment statistics, and calls the align() method with the human protein sequence and its closest Pan Troglodyte homolg sequence. The align() method returns the indexes  where there are differences between the 2 given protein sequences. The code then prints the alignment statistics and writes the different indexes (where there are differences between human and Pan Troglodyte sequences) to a file for further analysis.
- The output file looks like this: human ID, Pan Troglodytes ID, different indexes.

## PDB_Locations <a name="pdb_locations"></a>
- This code reads a file containing different amino acid indexes, extracts the indexes for each protein ID, and then parses a collection of PDB files. For each PDB file, it retrieves the atom locations for the specified amino acid indexes and stores them in a dictionary. Finally, it writes the protein ID and atom locations to a file in a specific format.
- The format looks like this: Human ID (A0A024RCN7), index of amino acid (44), the 3D location of each atom in the amino acid (x, y, z values.) The atoms in each amino acid are separated by semicolons.

## Extract_Locations <a name="extract_locations"></a>
- This code reads a file that contains protein IDs and their corresponding amino acid locations. It extracts the locations and indexes for each protein ID and stores them in dictionaries (locations_per_id and indexes_per_id). The locations_per_id dictionary maps each protein ID to a dictionary of different amino acid indexes and their corresponding atom locations. The indexes_per_id dictionary maps each protein ID to a list of different amino acid indexes. The function then returns these dictionaries as the output.

## Find_triplets <a name="find_triplets"></a>
- find_triplets(): This function calls the extract_locations() function to obtain a dictionary of atom locations and indexes per protein ID. It also retrieves identity and coverage data using the get_identities_coverage() function. The function then iterates over the protein IDs and their corresponding atom locations. It performs a series of steps to find triplets of amino acids where one amino acid is within a certain distance of the other two. The resulting triplets are stored in the triplets_per_id dictionary
