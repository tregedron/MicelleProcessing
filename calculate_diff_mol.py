from chemfiles import Trajectory, Selection
from tqdm import tqdm
from scripts.DiffCalc import DiffCalcerMol
from scripts.CorrCalc import CorrCalcer
from utils.utils import calc_list_correlation
import numpy as np
import pandas as pd
import os
import time
import argparse

def calculate_diff_coefficient_molecule(trj_path: str,
                                        topol_path: str,
                                        window: int,
                                        shift: int,
                                        atom="C1",
                                        trj_step_time=500):

    '''
    usage example:
    python main.py -trj data/100ns_NPT_1_8_ext.xtc -top data/100ns_NPT_1_8micelles.gro -w 5000 -s 1 > logs/log_1_8_1_ext

    :param trj_path:path to trajectory. The xtc format is usually used here.
    :param topol_path: path to topology in gro format. The atoms names are provided by it.
    :param window: number of steps in dr array calculations. The big enough window is necessary to be chosen for dr
    to pass balistic regime.
    :param shift: number of steps between dr array calculations.
    :param trj_step_time: time of 1 step in trajectory. The time that elapse between writing position coordinates
    :param atom: str, string of atom to calculate diffusion coefficient
    :return:
    '''
    print("Working on trajectory", trj_path)
    print("topology: ", topol_path)
    print("Window: ", window, " shift: ", shift)

    shifts = np.zeros(window, dtype=float)
    shifts_calculated_times = np.zeros(1, dtype=int)

    diff_calcer_list = []

    with Trajectory(topol_path) as trajectory:
        trajectory.set_topology(topol_path)
        for frame in trajectory:
            selection = Selection("name " + atom +" and resname C7H16")
            list_selection_ids = selection.evaluate(frame)
            positions = frame.positions[list_selection_ids]
            for i in range(len(positions)):
                diff_calcer_list.append(DiffCalcerMol(i, window, shift, 1))

            print(len(diff_calcer_list))

    with Trajectory(trj_path) as trajectory:
        trajectory.set_topology(topol_path)

        for frame in tqdm(trajectory):
            positions = frame.positions[list_selection_ids]
            for ind, position in enumerate(positions):
                diff_calcer_list[ind].update_coordinates(position, frame.cell.lengths)

            # if frame.step == 4500:
            #     for diff_calcer in diff_calcer_list:
            #         diff_calcer.collect_data(shifts, shifts_calculated_times)

                break
    for diff_calcer in diff_calcer_list:
        diff_calcer.collect_data(shifts, shifts_calculated_times)

    print(f"there were {shifts_calculated_times} shifts")
    shifts = shifts/shifts_calculated_times
    times = np.array([i * trj_step_time for i in range(0, window)])
    # ddr - numerical deviation d(dr^2)/d(t), dt is expected to be in ps.
    ddr = np.zeros(len(shifts))
    for i in range(1, len(shifts)):
        ddr[i] = (shifts[i] - shifts[i - 1]) / (times[i] - times[i - 1]) * 1000

    dr = pd.DataFrame({'time, fs': times, 'dr': shifts, "ddr": ddr})

    path_out_dir = os.path.join("results_heptane", trj_path.split("/")[-1].split(".")[0])
    os.makedirs(path_out_dir, exist_ok=True)

    dr.to_csv(os.path.join(path_out_dir, f"dr_{atom}_{window}_{shift}"), sep='\t')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Diffusion calculator')
    parser.add_argument('-trj', '--trajectory', default=None, type=str,
                        help='Trajectory file in xtc format (or not in xtc)')
    parser.add_argument('-top', '--topology', default=None, type=str,
                        help='Topology file in gro format')
    parser.add_argument('-w', '--window', default=None, type=int,
                        help='Number of frames in dr window')
    parser.add_argument('-s', '--shift', default=None, type=int,
                        help='Shift size between dr calculations')
    parser.add_argument('-a', '--atom-mass-center', default=None, type=str,
                        help='Atoms to consider as a mass center. For example, in heptane it can be C1, or O in water')
    args = parser.parse_args()

    start_time = time.time()
    # trj_temp = os.path.join("data", "100ns_NPT_8_8micelles.xtc")
    # top_temp = os.path.join("data", "100ns_NPT_8_8micelles.gro")
    # window_temp = 10
    # shift_temp = 1

    calculate_diff_coefficient_molecule(args.trajectory, args.topology, args.window, args.shift)
    # calculate_corr(args.trajectory, args.topology)
    # calculate_diff_coefficients(trj_temp, top_temp, window_temp, shift_temp)
    # calculate_corr(trj_temp, top_temp)
    print("--- %s seconds ---" % (time.time() - start_time))