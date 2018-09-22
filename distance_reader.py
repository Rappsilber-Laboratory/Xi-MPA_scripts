from structure.StructureContainer import StructureContainer
import pandas as pd


def get_distances_from_pdb(filename, uniprot_id, chain_id=" "):
    tmp_struct = StructureContainer(contact_atom="CA")
    tmp_struct.load_structure("xxxx", chain_id, filename, seqsep=1)
    distances = tmp_struct.get_contact_map().get_mapped_distances()
    residue_type_data = add_residue_types(distances, tmp_struct)
    distance_data_frame = transform_to_data_frame(residue_type_data, uniprot_id)
    return distance_data_frame


def add_residue_types(distances, structure):
    residue_type_data = []
    for index_i, index_j, distance in distances:
        if index_i <= index_j:
            res_type_i = structure.get_residue_type(index_i)
            res_type_j = structure.get_residue_type(index_j)
            residue_type_data.append((index_i, res_type_i, index_j, res_type_j, distance))
    return residue_type_data


def transform_to_data_frame(distances, uniprot_id):
    distance_data = pd.DataFrame()
    transpose_distance_data = zip(*distances)
    name_column = [uniprot_id] * len(distances)
    distance_data['uniprot_id'] = name_column
    distance_data['residue_1'] = transpose_distance_data[0]
    distance_data['residue_type_1'] = transpose_distance_data[1]
    distance_data['residue_2'] = transpose_distance_data[2]
    distance_data['residue_type_2'] = transpose_distance_data[3]
    distance_data['distance'] = transpose_distance_data[4]
    return distance_data
