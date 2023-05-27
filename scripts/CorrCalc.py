import numpy as np


class CorrCalcer:

    def __init__(self, hash_id, aggregate_size):
        self.hash_id = hash_id
        self.aggregate_size = aggregate_size
        self.coordinates = np.zeros((2, 3), dtype=float)
        self.cells = np.zeros((2, 3), dtype=float)
        self.boxes = np.zeros((2, 3), dtype=int)
        self.dr = None
        self.calculated_times = 0
        self.step_iterator = 0
        self.ready_to_calculate = False

    def update_coordinates(self, new_coordinates, cell):
        # update mass center
        self.coordinates[0] = self.coordinates[1]
        self.coordinates[1] = new_coordinates
        # update cell
        self.cells[0] = self.cells[1]
        self.cells[1] = np.array(cell)
        # box inheritance
        if self.step_iterator != 0:
            dr = self.coordinates[1] - self.coordinates[0]
            self.boxes[1] = self.boxes[0] - np.round(dr / self.cells[0])
        self.step_iterator += 1

    def calc_shift(self):
        if self.check_size():
            self.calculated_times += 1
            d_r = self.coordinates[1] - self.coordinates[0]
            d_box = (self.boxes[1] - self.boxes[0]) * 0.5 * (self.cells[0] + self.cells[1])
            self.dr = d_r + d_box
            self.ready_to_calculate = True

    def check_size(self):
        if (self.aggregate_size >= 402) and (self.aggregate_size <= 670):
            return True
        else:
            return False
