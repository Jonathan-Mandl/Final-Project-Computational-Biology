import itertools
import numpy as np
from Bio.SearchIO._legacy import NCBIStandalone
from extract_locations import extract_locations
import tarfile
import gzip


def get_identities_coverage(id_list):
    """
    Retrieves identity and coverage data for protein IDs.

    Args:
        id_list (list): List of protein IDs.

    Returns:
        tuple: A tuple containing two dictionaries:
            - identities_per_id: A dictionary mapping protein IDs to their identity information.
            - coverage_per_id: A dictionary mapping protein IDs to their coverage information.
    """
    seq_lengths = {}
    with open('human_sequences_no_repeats.fasta', 'r') as file:
        for line in file:
            if line.startswith(">"):
                id = line.strip()[1:]
                if id in id_list:
                    seq = next(file).strip()
                    length = len(seq)
                    seq_lengths[id] = length

    result_handle = open("nr_blast_output.txt", 'r')
    blast_parser = NCBIStandalone.BlastParser()
    blast_iterator = NCBIStandalone.Iterator(result_handle, blast_parser)

    identities_per_id = {}
    coverage_per_id = {}
    for blast_record in blast_iterator:
        query_id = blast_record.query.split()[0]
        if query_id in id_list:
            for alignment in blast_record.alignments:
                first_hsp = alignment.hsps[0]
                identities = first_hsp.identities
                identities_per_id[query_id] = identities
                alignment_length = first_hsp.query_end - first_hsp.query_start + 1
                query_length = seq_lengths[query_id]
                coverage = alignment_length / query_length
                coverage_per_id[query_id] = coverage
                break

    return identities_per_id, coverage_per_id


def get_title(id_list):
    """
    Retrieves the titles of PDB files for protein IDs.

    Args:
        id_list (list): List of protein IDs.

    Returns:
        dict: A dictionary mapping protein IDs to their corresponding PDB titles.
    """
    title_per_id = {}

    with tarfile.open("UP000005640_9606_HUMAN_v4.tar", "r") as tar:
        pdb_files = [member for member in tar.getmembers() if
                     member.name.endswith(".pdb.gz") and member.name.split("-")[2] == "F1" and
                     member.name.split("-")[1] in id_list]

        for member in pdb_files:
            protein_id = member.name.split("-")[1]

            with tar.extractfile(member) as f:
                contents = gzip.decompress(f.read())

            pdb_title = None

            for line_idx, line in enumerate(contents.decode('utf-8').split("\n")):
                if line.startswith("TITLE"):
                    next_line_idx = line_idx + 1
                    if next_line_idx < len(contents.decode('utf-8').split("\n")):
                        next_line = contents.decode('utf-8').split("\n")[next_line_idx]
                        if next_line.startswith("TITLE"):
                            pdb_title = line.split("FOR")[1].strip() + next_line.split("TITLE")[1].strip()
                        else:
                            pdb_title = line.split("FOR")[1].strip()
                        break

            title_per_id[protein_id] = pdb_title

    return title_per_id


def find_triplets():
    """
    Finds protein triplets based on atom locations, identity, and coverage data.
    Writes the results to a file.
    """
    # use extract locations method to get a dictionary of locations per id for every protein and dictionary index list for every protein iD
    locations_per_id, indexes_per_id = extract_locations()
    triplets_per_id = {}  # dictionary that contains triplets list for every protein ID
    id_list = locations_per_id.keys()  # list of protein IDs
    # extract identity and coverage per id dictionaries from blast output file using get_identities_coverage() method
    identity_per_id, coverage_per_id = get_identities_coverage(id_list)

    for human_id, locations in locations_per_id.items():
        identity = identity_per_id[human_id]
        identity_percentage = identity[0] / identity[1]
        coverage = coverage_per_id[human_id]
        # check that the identity and coverage meet our threshold
        if identity_percentage > 0.94 and coverage > 0.94:
            index_list = indexes_per_id[human_id]  # get index list of protein from indexes_per_id dictionary
            triplets = []
            neighbor_dict = {}
            # Create a dictionary of neighboring amino acids within 5 angstrom of distance
            num_aa = len(index_list)  # number of different indexes in protein
            for i in range(num_aa):
                # iterate over every atom location in index i and in every other index j and find indexes where there is euclidian distance less than 5
                neighbor_dict[i] = [j for j in range(num_aa) for point1 in locations[index_list[j]] for point2 in
                                    locations[index_list[i]] if (np.linalg.norm(point1 - point2) < 5 and i != j)]

            # Find triplets where one amino acid is within 5 angstrom of distance from the other two
            # iterate over every combination of indexes
            for i, j in itertools.combinations(range(num_aa), 2):
                # if indexes are equal, skip
                if i == j:
                    continue
                # define sets of neighbors of each amino acid index using neighbor_dict we created before
                neighbors_i = set(neighbor_dict[i])
                neighbors_j = set(neighbor_dict[j])
                # find common neighbors using intersection between 2 sets of neighbors
                common_neighbors = neighbors_i.intersection(neighbors_j)
                # iterate over common neighbors
                for k in common_neighbors:
                    # ensure common neighbors are not one of the 2 indexes we are checking
                    if k != i and k != j:
                        # triplets is the sorted tuple of the 2 indexes plus the common neighbors
                        triplet = sorted([index_list[i], index_list[j], index_list[k]])
                        # ensure triplets are more than 4 units apart in the sequence itself
                        if triplet[1] - triplet[0] > 4 and triplet[2] - triplet[1] > 4:
                            triplets.append(triplet)

            # filter only the unique triplets
            if len(triplets):
                unique_triplets = []
                for triplet in triplets:
                    if triplet not in unique_triplets:
                        unique_triplets.append(triplet)
                # add the unique triplets for our protein to the triplets_per_id dictionary in the human_id entry
                triplets_per_id[human_id] = unique_triplets

    id_list = triplets_per_id.keys() # id list of protein which contain triplets
    # get the titles for our proteins that contain triplets using title_per_id() method
    title_per_id = get_title(id_list)

    # add each protein and all of its triplets to the 3d triplets file
    with open("3d_triplets_every_atom.txt", "a") as output_file:
        for human_id, triplets in triplets_per_id.items():
            title = title_per_id[human_id]
            identity = identity_per_id[human_id]
            coverage = coverage_per_id[human_id]
            row = f"{human_id}\t{title}\tIdentity:{identity[0]}/{identity[1]}\tCoverage:{coverage}\t"
            for triplet in triplets:
                row += f"{triplet[0]},{triplet[1]},{triplet[2]}"
                row += "\t"
            row += "\n"
            output_file.write(row)
