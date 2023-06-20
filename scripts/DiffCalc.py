import numpy as np


class DiffCalcer:

    def __init__(self, hash_id, window, shift, aggregate_size):
        self.window = window
        self.shift = shift
        self.hash_id = hash_id
        self.aggregate_size = aggregate_size

        self.coordinates = np.zeros((window, 3), dtype=float)
        self.cells = np.zeros((window, 3), dtype=float)
        self.boxes = np.zeros((window, 3), dtype=int)
        self.dr = np.zeros(window, dtype=float)
        self.step_iterator = 0
        self.till_last_time_calculated = shift
        self.ready_to_calculate = False
        self.shifts = np.zeros(window, dtype=float)
        self.calculated_times = 0

    def update_coordinates(self, new_coordinates, cell):
        # update mass center
        self.coordinates[self.step_iterator % self.window] = new_coordinates
        # update cell
        self.cells[self.step_iterator % self.window] = np.array(cell)
        # box inheritance
        if self.step_iterator != 0:
            dr = self.coordinates[self.step_iterator % self.window] - self.coordinates[(self.step_iterator -1) % self.window]
            # if any(elt != 0 for elt in np.round(dr / self.cells[self.step_iterator % self.window])):
            #     print("PBC")
            self.boxes[self.step_iterator % self.window] = self.boxes[(self.step_iterator - 1) % self.window] - np.round(dr / self.cells[self.step_iterator % self.window])

        if (self.step_iterator >= self.window - 1) and (self.till_last_time_calculated >= self.shift):
            self.calculate_shifts()

        self.step_iterator += 1
        self.till_last_time_calculated += 1

    def check_size(self):
        # if (self.aggregate_size >= 536) and (self.aggregate_size <= 536):
        if (self.aggregate_size >= 402) and (self.aggregate_size <= 670):
            return True
        else:
            return False

    def calculate_shifts(self):
        if self.check_size():
            self.calculated_times += 1
            self.till_last_time_calculated = 0

            d_r = self.coordinates - self.coordinates[(self.step_iterator + 1) % self.window]
            d_box = (self.boxes - self.boxes[(self.step_iterator + 1) % self.window])\
                        * 0.5 * (self.cells[(self.step_iterator + 1) % self.window] + self.cells)
            distance = ((d_r + d_box) ** 2).sum(axis=-1)
            for i in range(self.window):
                self.shifts[i] += distance[(self.step_iterator + 1 + i) % self.window]

    def collect_data(self, shifts, shifts_calculated_times):
        shifts += self.shifts
        shifts_calculated_times += self.calculated_times

class DiffCalcerMol:

    def __init__(self, hash_id, window, shift, aggregate_size):
        self.window = window
        self.shift = shift
        self.hash_id = hash_id
        self.aggregate_size = aggregate_size

        self.coordinates = np.zeros((window, 3), dtype=float)
        self.cells = np.zeros((window, 3), dtype=float)
        self.boxes = np.zeros((window, 3), dtype=int)
        self.dr = np.zeros(window, dtype=float)
        self.step_iterator = 0
        self.till_last_time_calculated = shift
        self.ready_to_calculate = False
        self.shifts = np.zeros(window, dtype=float)
        self.calculated_times = 0

    def update_coordinates(self, new_coordinates, cell):
        # update mass center
        self.coordinates[self.step_iterator % self.window] = new_coordinates
        # update cell
        self.cells[self.step_iterator % self.window] = np.array(cell)
        # box inheritance
        if self.step_iterator != 0:
            dr = self.coordinates[self.step_iterator % self.window] - self.coordinates[(self.step_iterator -1) % self.window]
            # if any(elt != 0 for elt in np.round(dr / self.cells[self.step_iterator % self.window])):
            #     print("PBC")
            self.boxes[self.step_iterator % self.window] = self.boxes[(self.step_iterator - 1) % self.window] - np.round(dr / self.cells[self.step_iterator % self.window])

        if (self.step_iterator >= self.window - 1) and (self.till_last_time_calculated >= self.shift):
            self.calculate_shifts()

        self.step_iterator += 1
        self.till_last_time_calculated += 1

    def calculate_shifts(self):

        self.calculated_times += 1
        self.till_last_time_calculated = 0

        d_r = self.coordinates - self.coordinates[(self.step_iterator + 1) % self.window]
        d_box = (self.boxes - self.boxes[(self.step_iterator + 1) % self.window])\
                    * 0.5 * (self.cells[(self.step_iterator + 1) % self.window] + self.cells)
        distance = ((d_r + d_box) ** 2).sum(axis=-1)
        for i in range(self.window):
            self.shifts[i] += distance[(self.step_iterator + 1 + i) % self.window]

    def collect_data(self, shifts, shifts_calculated_times):
        shifts += self.shifts
        shifts_calculated_times += self.calculated_times
