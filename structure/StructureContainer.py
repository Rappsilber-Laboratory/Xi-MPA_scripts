import ContactMap
import Bio.PDB


class StructureContainer:
    
    def __init__(self, contact_atom="CB"):
        
        self.id = None 
        self.chain_id = None 
        self.pdb_filename = None
        self.bio_pdb_pose = None 
        self.pose = None
        self.contact_map = None 
        self.centroid = None
        self.contact_atom=contact_atom
           
    def load_structure(self, input_id, input_chain_id, input_pdb_filename , seqsep = 12):
        """Loads a structure and initializes respective objects"""
        
        self.id = input_id 
        self.chain_id = input_chain_id 
        self.pdb_filename = input_pdb_filename
        
        structure = Bio.PDB.PDBParser().get_structure(self.id, self.pdb_filename)
        model = structure[0]
        self.bio_pdb_pose = model[ self.chain_id ]
        
        self.contact_map = ContactMap.ContactMap(contact_atom=self.contact_atom)
        self.contact_map.calculate_cm_from_struct( self.bio_pdb_pose , seqsep)
        self.centroid = self.contact_map.calc_centroid(self.bio_pdb_pose)
        
    def get_centroid(self):
        return self.centroid

    def get_id(self):
        return self.id

    def get_chain_id(self):
        return self.id

    def get_residue_type(self, residue_id):
        return (self.bio_pdb_pose[residue_id]).get_resname()

    def get_residue_coord(self, residue_id):
        residue = self.bio_pdb_pose[residue_id]
        
        if residue.get_resname() == "GLY":
            atom_res1 = "CA"
        else:
            atom_res1 = "CB" # CB here for CB definition!
    
        return residue[atom_res1].coord 

    def get_pdb_filename(self):
        return self.pdb_filename

    def get_bioPBD_pose(self):
        return self.bio_pdb_pose

    def get_rosetta_pose(self):
        return self.pose

    def get_contact_map(self):
        return self.contact_map

    def get_number_of_residues(self):
        return len(self.bio_pdb_pose)

