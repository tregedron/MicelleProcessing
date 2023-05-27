from chemfiles import Trajectory, Selection
from tqdm import tqdm
from scripts.Cluster import Cluster
from scripts.Micelle import Micelle
from scripts.DiffCalc import DiffCalcer
from scripts.CorrCalc import CorrCalcer
from utils.utils import calc_list_correlation
import numpy as np
import pandas as pd
import os
import time


def calculate_diff_coefficients(trj_path: str, topol_path: str, window: int, shift: int):
    print("Working on trajectory", trj_path)
    print("topology: ", topol_path)
    print("Window: ", window, " shift: ", shift)

    shifts = np.zeros(window, dtype=float)
    shifts_calculated_times = np.zeros(1, dtype=int)

    with Trajectory(trj_path) as trajectory:
        trajectory.set_topology(topol_path)

        selection = Selection("name O1 or name O2 or name O3 or name O4 or name O5")

        micelles_list = []
        diff_calcer_list = []
        ids = {}

        for frame in tqdm(trajectory):

            for micelle in micelles_list:
                micelle.drop_updated()

            # cluster the frame
            cluster = Cluster(frame, selection)
            cluster.cluster_frame()

            number_of_clusters = len(set(cluster.clustering))
            cluster_indexes = []

            # collect clusters
            for cluster_id in range(number_of_clusters):
                cluster_indexes.append([])
            for atom_ind, cluster in enumerate(cluster.clustering):
                cluster_indexes[cluster].append(atom_ind)

            for cluster_id in range(number_of_clusters):
                micelle_hash_id = hash(''.join(str(x) for x in cluster_indexes[cluster_id]))
                if micelle_hash_id in ids.keys():
                    micelles_list[ids[micelle_hash_id]].update_positions(frame, cluster_indexes[cluster_id])
                else:
                    ids[micelle_hash_id] = len(micelles_list)
                    micelles_list.append(Micelle(micelle_hash_id, frame, cluster_indexes[cluster_id]))
                    diff_calcer_list.append(DiffCalcer(micelle_hash_id, window, shift, len(cluster_indexes[cluster_id])))

            positions_to_delete = []
            for ind, micelle in enumerate(micelles_list):
                if micelle.updated:
                    micelle.calculate_mass_center(frame.cell.lengths)
                    micelle.calculate_gyradius(frame.cell.lengths)
                    diff_calcer_list[ind].update_coordinates(micelle.mass_center, frame.cell.lengths)
                else:
                    positions_to_delete.append(ids.pop(micelle.hash_id))
            if positions_to_delete:
                for position in reversed(sorted(positions_to_delete)):
                    for hash_id in ids.keys():
                        if ids[hash_id] > position:
                            ids[hash_id] -= 1

                    diff_calcer_list[position].collect_data(shifts, shifts_calculated_times)
                    micelles_list.pop(position)
                    diff_calcer_list.pop(position)

            print("Step = ", frame.step, len(micelles_list), " micelles ", ids)
            # if frame.step == 5000:
            #     for diff_calcer in diff_calcer_list:
            #         diff_calcer.collect_data(shifts, shifts_calculated_times)
            #     break

    print(f"there were {shifts_calculated_times} shifts")
    shifts = shifts/shifts_calculated_times

    df = pd.DataFrame({'dr': shifts})
    df["ddr"] = df["dr"]
    for i in range(1, window):
        df["ddr"][i] = df["dr"][i] - df["dr"][i-1]
    df['time, ps'] = np.array([i*500 for i in range(0, window)])
    df = df[['time, ps', 'dr', 'ddr']]

    path_out_dir = os.path.join("results", trj_path.split("/")[-1].split(".")[0])
    os.makedirs(path_out_dir, exist_ok=True)

    df.to_csv(os.path.join(path_out_dir, f"dr_{window}_{shift}"), sep='\t')


def calculate_corr(trj_path: str, topol_path: str):
    print("Working on trajectory", trj_path)
    print("topology: ", topol_path)
    print("Calculating corr coefficient")

    corr = np.zeros(1, dtype=float)
    corr_calculated_times = np.zeros(1, dtype=int)

    df = pd.DataFrame(columns=['time', 'corr'])

    with Trajectory(trj_path) as trajectory:
        trajectory.set_topology(topol_path)

        selection = Selection("name O1 or name O2 or name O3 or name O4 or name O5")

        micelles_list = []
        ids = {}
        corr_calcer_list = []

        for frame in tqdm(trajectory):
            corr_frame = np.zeros(1, dtype=float)
            corr_calculated_times_frame = 0
            for micelle in micelles_list:
                micelle.drop_updated()

            # cluster the frame
            cluster = Cluster(frame, selection)
            cluster.cluster_frame()

            number_of_clusters = len(set(cluster.clustering))
            cluster_indexes = []

            # collect clusters
            for cluster_id in range(number_of_clusters):
                cluster_indexes.append([])
            for atom_ind, cluster in enumerate(cluster.clustering):
                cluster_indexes[cluster].append(atom_ind)

            for cluster_id in range(number_of_clusters):
                micelle_hash_id = hash(''.join(str(x) for x in cluster_indexes[cluster_id]))
                if micelle_hash_id in ids.keys():
                    micelles_list[ids[micelle_hash_id]].update_positions(frame, cluster_indexes[cluster_id])
                else:
                    ids[micelle_hash_id] = len(micelles_list)
                    micelles_list.append(Micelle(micelle_hash_id, frame, cluster_indexes[cluster_id]))
                    corr_calcer_list.append(CorrCalcer(micelle_hash_id, len(cluster_indexes[cluster_id])))

            positions_to_delete = []
            for ind, micelle in enumerate(micelles_list):
                if micelle.updated:
                    micelle.calculate_mass_center(frame.cell.lengths)
                    micelle.calculate_gyradius(frame.cell.lengths)
                    corr_calcer_list[ind].update_coordinates(micelle.mass_center, frame.cell.lengths)
                    corr_calcer_list[ind].calc_shift()
                else:
                    positions_to_delete.append(ids.pop(micelle.hash_id))
            if positions_to_delete:
                for position in reversed(sorted(positions_to_delete)):
                    for hash_id in ids.keys():
                        if ids[hash_id] > position:
                            ids[hash_id] -= 1
                    micelles_list.pop(position)
                    corr_calcer_list.pop(position)

            corr_frame, corr_calculated_times_frame = calc_list_correlation(corr_calcer_list)

            corr += corr_frame
            corr_calculated_times += corr_calculated_times_frame
            if corr_calculated_times != 0:
                df.loc[len(df.index)] = [frame.step, (corr_frame/corr_calculated_times_frame).item()]

            # if frame.step == 5000:
            #     break

    path_out_dir = os.path.join("results", trj_path.split("/")[-1].split(".")[0])
    os.makedirs(path_out_dir, exist_ok=True)

    df.to_csv(os.path.join(path_out_dir, f"corr_full"), sep='\t')

    print(corr/corr_calculated_times)
    print(corr_calculated_times)


if __name__ == '__main__':
    # parser = argparse.ArgumentParser(description='PyTorch Template')
    # parser.add_argument('-trj', '--trajectory', default=None, type=str,
    #                     help='trajectory file in xtc format (or not in xtc)')
    # parser.add_argument('-top', '--topology', default=None, type=str,
    #                     help='topology file in gro format')
    # parser.add_argument('-w', '--window', default=None, type=str,
    #                     help='number of frames in dr window')
    # parser.add_argument('-s', '--shift', default=None, type=str,
    #                     help='shift size between dr calculations')
    # args = parser.parse_args()

    start_time = time.time()
    trj_temp = os.path.join("data", "100ns_NPT_8_8micelles.xtc")
    top_temp = os.path.join("data", "100ns_NPT_8_8micelles.gro")
    window_temp = 5000
    shift_temp = 1

    # calculate_diff_coefficients(args.trajectory, args.topology, args.window, args.shift)
    calculate_diff_coefficients(trj_temp, top_temp, window_temp, shift_temp)
    calculate_corr(trj_temp, top_temp)
    print("--- %s seconds ---" % (time.time() - start_time))
