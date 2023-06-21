import numpy as np

def extract_locations():
    """
    Extracts atom locations and indexes per protein ID from a file.

    Returns:
        tuple: A tuple containing two dictionaries:
            - locations_per_id: A dictionary mapping protein IDs to dictionaries of different amino acid indexes
                               and their corresponding atom locations.
            - indexes_per_id: A dictionary mapping protein IDs to lists of different amino acid indexes.
    """

    with open('3D_locations_every_atom', 'r') as file: # open 3d locations file
        locations_per_id = {}  # mapping between different aa locations and id
        indexes_per_id = {}  # mapping between different aa indexes and id
        # parse locations file line by line. each line is for different protein
        for line in file:
            cols = line.strip().split('\t') # split column by tabs because file is tab separated
            human_id = cols[0] # human id is first column
            locations = {}  # mapping of different indexes to a list of their atom's locations
            index_list = []  # list of different indexes
            for column in cols[1:]: # iterate over rest of columns
                index, atom_locations = column.strip().split(":")
                atom_locations=atom_locations[:-1]
                atom_locations=atom_locations.strip().split(";")
                atom_coordinates_list=[]
                for atom in atom_locations:
                    x,y,z=atom.strip().split(",")
                    x = float(x)
                    y = float(y)
                    z = float(z)
                    atom_coordinates_list.append(np.array([x,y,z]))
                index = int(index)
                locations[index] = atom_coordinates_list
                index_list.append(index)

            locations_per_id[human_id] = locations # add locations dictionary to specific human id entry in locations_per_id dictionary
            indexes_per_id[human_id] = index_list # add index_list dictionary to specific human id entry in indexes_per_id dictionary

        return locations_per_id, indexes_per_id
