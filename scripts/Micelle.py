import numpy as np
from utils.utils import custom_distance_to_point

mass_dictionary = {"O": 15.999, "C": 12.0096, "H": 1.00784}


class Micelle:

    def __init__(self, hash_id, frame, indexes):
        self.hash_id = hash_id
        self.indexes = indexes
        self.positions = frame.positions[indexes]
        self.size = len(self.positions)
        self.gyradius = 0
        self.mass_center = None
        self.masses = np.array([mass_dictionary[frame.topology.atoms[i].name[0]] for i in self.indexes])
        self.updated = True

    def update_positions(self, new_frame, indexes):
        self.updated = True
        self.positions = new_frame.positions[indexes]
        self.masses = np.array([mass_dictionary[new_frame.topology.atoms[i].name[0]] for i in self.indexes])

    def drop_updated(self):
        self.updated = False

    def calculate_mass_center(self, cell):
        '''
        Calculate the mass center for given coordinates with respect to the atom-mass dictionary.
        :param frame: chemfiles frame with topology
        :param indexes_in_obj: indexes of the object on the frame
        :param cell: cell lengths
        :return: mass center position
        '''
        self.mass_center = np.zeros(3, dtype=float)
        theta = self.positions/cell * 2 * np.pi
        ksi = np.cos(theta)
        dzeta = np.sin(theta)
        ksi = np.average(ksi, weights=self.masses, axis=0)
        dzeta = np.average(dzeta, weights=self.masses, axis=0)
        theta = np.arctan2(-ksi, -dzeta) + np.pi
        self.mass_center = theta / (2*np.pi) * cell

    def calculate_gyradius(self, cell):
        self.gyradius = np.average(custom_distance_to_point(self.positions, self.mass_center, cell),
                                   weights=self.masses, axis=0)

    def ask_size(self):
        return self.size