import numpy
import Bio.PDB

"""It might be a little confusing that ContactMap class is currently used to calculate contact maps AND distance maps. It is mostly done for efficiency, since
   one needs to compute the distance map to derive the contact map anyway.... Maybe it should be split into two classes that inherent from a base"""


class ContactMap:
    def __init__(self, contact_atom="CB"):
        self.contact_map = None
        self.distance_map = None
        self.cutoff = 8
        self.contact_atom = contact_atom
        self.resid_map = {}

        self.tp_predictions = None
        self.fp_predictions = None

    def load_cm_from_pdb(self, id, chain_id, pdb_file):
        structure = Bio.PDB.PDBParser().get_structure(id, pdb_file)
        model = structure[0]
        self.calc_contact_matrix(model[chain_id], model[chain_id])

    def calculate_cm_from_struct(self, bio_pdb_struct, input_seqsep):
        self.calc_contact_matrix(bio_pdb_struct, bio_pdb_struct, input_seqsep)

    def load_tp_predictions(self, contact_list):
        if len(contact_list) >= 1:
            self.tp_predictions = contact_list

    def load_fp_predictions(self, contact_list):
        if len(contact_list) >= 1:
            self.fp_predictions = contact_list

    def calc_residue_dist(self, residue_one, residue_two):
        """Returns the C-alpha distance between two residues"""
        atom_res1 = ''
        atom_res2 = ''
        if residue_one.get_resname() == "GLY":
            # print residue_one.get_resname()
            atom_res1 = "CA"
        else:
            atom_res1 = self.contact_atom  # CB here for CB definition!

        if residue_two.get_resname() == "GLY":
            # print residue_one.get_resname()
            atom_res2 = "CA"
        else:
            atom_res2 = self.contact_atom  # CB here for CB definition!

        try:
            # print residue_one.get_resname(), residue_two.get_resname()
            diff_vector = residue_one[atom_res1].coord - residue_two[atom_res2].coord
        except:
            diff_vector = residue_one["CA"].coord - residue_two["CA"].coord

        return numpy.sqrt(numpy.sum(diff_vector * diff_vector))

    def calc_contact_matrix(self, chain_one, chain_two, seqsep=12):
        """Returns a matrix of C-beta distances (C-alpha wheren residue is glycine) between two chains"""
        counter = 0
        row_count = 0
        col_count = 0
        self.contact_map = numpy.zeros((len(chain_one), len(chain_two)), numpy.float)
        self.distance_map = numpy.zeros((len(chain_one), len(chain_two)), numpy.float)
        for row, residue_one in enumerate(chain_one):

            for col, residue_two in enumerate(chain_two):

                dist = self.calc_residue_dist(residue_one, residue_two)

                r_id_1 = int(residue_one.__repr__().split()[3].split('=')[1])
                r_id_2 = int(residue_two.__repr__().split()[3].split('=')[1])
                # self.resid_map
                # print r_id_1, r_id_2
                self.resid_map[(r_id_1, r_id_2)] = (row, col)
                self.distance_map[row, col] = dist
                # print "GAAA"

                # if dist <= float(self.cutoff) and abs(col_count - row_count) < 24 and abs(col_count - row_count) >= 12 :
                if dist <= float(self.cutoff) and abs(row - col) >= seqsep:

                    self.contact_map[row, col] = 1.0
                    counter += 1
                else:

                    self.contact_map[row, col] = 0
                col_count += 1
            col_count = 0
            row_count += 1

    def number_of_contacts(self, seqsep_min=6, seqsep_max=999, resnum=None):

        num_contacts = 0

        for resnum_lower, residue_lower in enumerate(self.contact_map):

            for resnum_upper, residue_upper in enumerate(residue_lower):

                if resnum_lower < resnum_upper and self.contact_map[resnum_lower][resnum_upper] == 1 and abs(
                                resnum_upper - resnum_lower) >= seqsep_min and abs(
                                resnum_upper - resnum_lower) <= seqsep_max:
                    if resnum == None:
                        num_contacts += 1
                    else:
                        if (resnum_lower + 1 == resnum):
                            num_contacts += 1

                            # print resnum_lower, resnum_upper
        return num_contacts

    def calc_centroid(self, chain_one):
        x = []
        y = []
        z = []
        len_chain = float(len(chain_one))
        for residue in chain_one:
            atom_res1 = ''
            if residue.get_resname() == "GLY":
                # print residue_one.get_resname()
                atom_res1 = "CA"
            else:
                atom_res1 = "CA"  # CB here for CB definition!
            x.append(residue[atom_res1].coord[0])
            y.append(residue[atom_res1].coord[1])
            z.append(residue[atom_res1].coord[2])
        return [numpy.mean(x), numpy.mean(y), numpy.mean(z)]

    def print_cm(self):
        print self.contact_map

    def print_dm(self):
        print self.distance_map

    def print_res_format(self):
        # residue_lower = 0
        # residue_upper = 0
        res_format = []
        for resnum_lower, residue_lower in enumerate(self.contact_map):
            for resnum_upper, residue_upper in enumerate(residue_lower):
                if resnum_lower <= resnum_upper and self.contact_map[resnum_lower][resnum_upper] == 1:
                    # print resnum_lower, resnum_upper
                    res_format.append((resnum_lower, resnum_upper))
        return res_format

    def distances(self, seq_sep, res):
        d_vector = []
        for resnum_lower, residue_lower in enumerate(self.distance_map):
            for resnum_upper, residue_upper in enumerate(residue_lower):
                if resnum_lower <= resnum_upper and abs(resnum_lower - resnum_upper) >= 12:
                    if res == 1:
                        d_vector.append(
                            (resnum_lower + 1, resnum_upper + 1, self.distance_map[resnum_lower][resnum_upper]))
                    else:
                        d_vector.append(self.distance_map[resnum_lower][resnum_upper])
        return d_vector

    def draw_cm(self, len_protein, name="cm"):
        pylab.rcParams.update({'font.size': 9,
                               'font.name': 'Arial',
                               'xtick.major.size': 1,
                               'xtick.major.width': 0.1,
                               'ytick.major.size': 1,
                               'ytick.major.width': 0.1,
                               'xtick.minor.size': 0,
                               'xtick.minor.width': 0.0,
                               'ytick.minor.size': 0,
                               'ytick.minor.width': 0.0})

        draw_map = numpy.zeros((len_protein, len_protein), numpy.float)
        native_map_x = []
        native_map_y = []
        for i in xrange(0, self.contact_map.shape[0]):
            for j in xrange(0, self.contact_map.shape[1]):
                if self.contact_map[i, j] == 1:
                    # print i+1,j+1
                    native_map_x.append(i + 1)
                    native_map_y.append(j + 1)
                draw_map[i, j] = self.contact_map[i, j]

        # draw_map = self.contact_map
        draw_map = draw_map * 0.2
        # draw_map[0][0] = 1

        # m = pylab.imshow(draw_map,cmap="binary", interpolation="nearest")
        # m = pylab.scatter(native_map_x, native_map_y,c="k", s="40")

        # print native_map_x, native_map_y
        fig = pylab.figure(linewidth=0.5)
        ax = fig.add_subplot(111)
        m = pylab.scatter(native_map_x, native_map_y, 8, [0.5, 0.5, 0.5], marker='o', linewidth=0.0)

        if self.fp_predictions != None:
            x_list = []
            y_list = []
            for x, y, p in self.fp_predictions:
                x_list.append(x)
                y_list.append(y)
            # m = pylab.scatter(x_list,y_list, c='#B40426', s=100,  linewidths=0.0, marker = (5,1))#3B4CC0
            # UNCOMMENT HERE
            m = pylab.scatter(x_list, y_list, c='#B40426', s=10, linewidths=0.0, marker=(5, 1))
        if self.tp_predictions != None:
            x_list = []
            y_list = []
            for x, y, p in self.tp_predictions:
                x_list.append(x)
                y_list.append(y)
            # m = pylab.scatter(x_list,y_list, c='#B40426', s=100,  linewidths=0.0, marker = (5,1))#3B4CC0
            # m = pylab.scatter(x_list,y_list, c='#3B4CC0', s=60,  linewidths=0.0, marker = (5,1))
            # UNCOMMENT HERE
            m = pylab.scatter(x_list, y_list, c='#3B4CC0', s=10, linewidths=0.0, marker=(5, 1))

            # print ax.get_yticks()
        ax.set_yticks(ax.get_yticks()[2:-1])
        ax.set_xticks(ax.get_xticks()[2:-1])
        # leg = pylab.legend(('true positives', "false positives", "native"), #'Rosetta Top5',# 'Cluster Top5'  ),
        # loc='upper center', bbox_to_anchor=(0.5, -0.05), fancybox=False, shadow=False, ncol=5,scatterpoints = 1,prop={'size':6})
        # leg.draw_frame(False)
        # for t in leg.get_texts():
        #    t.set_fontsize('small')
        m = pylab.plot([0, len_protein], [0, len_protein], '-k', linewidth=0.25)

        pylab.xlim((0, len_protein))
        pylab.ylim((0, len_protein))

        # fig = matplotlib.pyplot.gcf()
        # fig.legend()
        [i.set_linewidth(0.25) for i in ax.spines.itervalues()]
        fig.set_size_inches(2.2, 2.2)
        # print fig.gca().get_lines()
        pylab.setp(fig.gca().get_lines()[0], markeredgewidth=0.5)
        # rec = fig.patch
        # print rec
        pylab.savefig("%s.svg" % (name), dpi=300)
        pylab.savefig("%s.png" % (name), dpi=300)

    def get_mapped_distance(self, i, j):
        if self.resid_map.has_key((i, j)):
            mapped_tuple = self.resid_map[(i, j)]
            return self.distance_map[mapped_tuple[0]][mapped_tuple[1]]
        else:
            return 999.9

    def get_mapped_distances(self):
        distances = []
        for key, value in self.resid_map.iteritems():
            residue_tuple = key
            distance = self.distance_map[value[0]][value[1]]
            distances.append((residue_tuple[0], residue_tuple[1], distance))
        return distances

    def get_cm(self):
        return self.contact_map

    def get_dm(self):
        return self.distance_map

#
