import tarfile
import gzip
import io
from Bio.PDB import PDBParser

def pull_locations():
    # create dictionary of different aa indexes for each protein id.
    diff_indexes = {}
    # create dictionary of different aa locations for each protein id.
    location_per_id = {}

    # open diff_indexes file to extract the different indexes in each protein
    with open('diff_indexes_updated.txt', 'r') as f:
        # Iterate over each line in the file
        for line in f:
            cols = line.strip().split('\t')  # split line by tab because file is tab separated
            human_id = cols[0]  # human id is first column

            if len(cols) >= 3 and cols[2]:  # ensure 2nd column which contains indexes exists
                indexes = [int(i) for i in cols[2].split(',')]  # split indexes column by , because the indexes are separated by commas
                diff_indexes[human_id] = indexes  # store different indexes for the human id

    # open tar file which contains all pdb files of human alphafold proteins
    with tarfile.open("UP000005640_9606_HUMAN_v4.tar", "r") as tar:
        # Iterate over the file names in the tar archive
        for member in tar.getmembers():
            # Check if the file is a PDB file and only select F1 model
            if member.name.endswith(".pdb.gz") and member.name.split("-")[2] == "F1":
                # Extract the contents of the PDB file from the tar archive
                contents = tar.extractfile(member).read()
                # Decompress the contents of the PDB file using the gzip module
                contents = gzip.decompress(contents)
                # extract human id from file name
                human_id = member.name.split("-")[1]
                if human_id in diff_indexes:
                    differences = diff_indexes[human_id]
                else:
                    continue
                # Create a PDBParser object
                parser = PDBParser()
                # Parse the PDB file using Biopython's PDBParser
                structure = parser.get_structure(member.name, io.StringIO(contents.decode('utf-8')))
                # Get the first model (there may be multiple models in a PDB file)
                model = structure[0]
                locations = {}  # create dictionary of different locations for specific protein
                if len(differences):  # check that there are more than 0 different indexes
                    # Iterate over all residues in the model
                    for chain in model:
                        for residue in chain:
                            # Get the residue ID and index
                            res_id = residue.get_id()[1]
                            # check if residue is in the list of different indexes
                            if res_id in differences:
                                atom_locations = []
                                # Iterate over all atoms in the residue
                                for atom in residue:
                                    # atom_name = atom.get_name()
                                    # Get the coordinates of the atom
                                    x, y, z = atom.get_coord()
                                    # Create a unique identifier for the atom using residue ID and atom name
                                    # Add the atom location to the locations dictionary
                                    atom_locations.append((x, y, z))
                                # add list of atom locations to specific res_id in the dictionary
                                locations[res_id] = atom_locations
                    # add locations dictionary to location_per_id dictionary at the entry of human id
                    location_per_id[human_id] = locations

    # add atom locations and indexes of different amino acids for each protein to the 3d_locations file
    with open("3D_locations_every_atom", "a") as output_file:
        for human_id, locations in location_per_id.items():
            # Create the row string with the human ID
            row = f"{human_id}\t"
            # Iterate over each key and value pair: aa index and atom location list in the locations dictionary
            for index, atom_locations in locations.items():
                row += f"{index}:"  # add index in the beginning
                for location in atom_locations:
                    # add each atom location of specific aa index
                    row += f"{location[0]:.3f},{location[1]:.3f},{location[2]:.3f};"
                row += "\t"
            # Add a newline character to the end of the row
            row += "\n"
            # Write the row to the output file
            output_file.write(row)
