from utils import linear
from utils import constant_func

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import argparse
import os


def func(x, a, b):
    return a*x+b


def fit_dr(dr):
    pass


def process_dr(dr_path, start):
    # TODO Add the parameter for conversion of dr.index to dt. 0.5 fix! maybe convert the second column to the time.
    '''
    The function gets path to the set of time-mean square dispalecement (MSD) dependencies. The approximation is the
    following: 6*D = lim d<\delta r^2>/dt. Be careful the time is expected to be in ps, dr^2 in A^2. dr.index has 0.5
    because in my trajectories I wrote coordinates in file every 0.5 ps.

    :param dr_path:
    :param start:
    :return:
    '''
    print(dr_path)
    out_dir = dr_path.split("/")[:-1]
    out_dir = os.path.join(*out_dir)
    dr = pd.read_csv(dr_path, sep='\t', index_col=0)

    # # this code will be deleted in the future.
    # # ddr - numerical deviation d(dr^2)/d(t), dt is expected to be in ps.
    # for i in range(1, len(dr["dr"])):
    #     dr["ddr"].loc[i] = (dr["dr"][i] - dr["dr"][i - 1]) / (dr['time, fs'][i] - dr['time, fs'][i - 1]) * 1000

    if not dr.isnull().any().any():

        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot(111)
        ax.plot(0.5*dr.index, dr["dr"], 'b-', label='dr', markersize=3)
        ax.set_title(f'dr of time {dr_path.split("/")[-1]}', fontsize=22, pad=8)
        plt.legend()
        plt.savefig(os.path.join(out_dir, f'{dr_path.split("/")[-1]}_dr_t.png'), bbox_inches='tight')
        fig.tight_layout()
        plt.close()

        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot(111)
        ax.plot(0.5*dr.index[1:], dr["ddr"][1:], "o", color="red", label='ddr', markersize=3)
        popt, pcov = curve_fit(linear, 0.5*dr.index[start:], dr["ddr"][start:])
        print(popt)
        ax.plot(0.5*dr.index[start:], linear(0.5*dr.index[start:], *popt), 'b-', label='fit: a=%f, b=%f' % tuple(popt))
        ax.set_title(f'ddr of time {dr_path.split("/")[-1]}', fontsize=22, pad=8)
        plt.legend()
        plt.savefig(os.path.join(out_dir, f'{dr_path.split("/")[-1]}_ddr_fit_linear.png'), bbox_inches='tight')
        fig.tight_layout()
        plt.close()


        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot(111)
        ax.plot(0.5*dr.index[1:], dr["ddr"][1:], "o", color="red", label='ddr', markersize=3)
        popt, pcov = curve_fit(constant_func, 0.5*dr.index[start:], dr["ddr"][start:])
        print(popt)
        ax.plot(0.5*dr.index[start:], constant_func(0.5*dr.index[start:], *popt), 'b-', label='fit: const=%f' % tuple(popt))
        ax.set_title(f'ddr of time {dr_path.split("/")[-1]}', fontsize=22, pad=8)
        plt.legend()
        plt.savefig(os.path.join(out_dir, f'{dr_path.split("/")[-1]}_ddr_fit_const.png'), bbox_inches='tight')
        fig.tight_layout()
        plt.close()

def process_size_distribution(dist_path):
    monomer_size = 67
    out_dir = dist_path.split("/")[:-1]
    out_dir = os.path.join(*out_dir)
    dist = pd.read_csv(dist_path, sep='\t', index_col=0)

    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111)
    ax.plot(dist.index[0:20*monomer_size]/monomer_size, dist["size"][0:20*monomer_size], 'b-', label='dr', markersize=3)
    ax.set_title(f'Micelle size distribution in {dist_path.split("/")[-1]}', fontsize=22, pad=8)
    plt.legend()
    plt.savefig(os.path.join(out_dir, f'{dist_path.split("/")[-1]}_size_dist.png'), bbox_inches='tight')
    fig.tight_layout()
    plt.close()

def process_number_distribution(dist_path):
    out_dir = dist_path.split("/")[:-1]
    out_dir = os.path.join(*out_dir)
    dist = pd.read_csv(dist_path, sep='\t', index_col=0)

    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111)
    ax.plot(dist.index, dist["number"], 'b-', label='dr', markersize=3)
    ax.set_title(f'Number of micelles distribution in {dist_path.split("/")[-1]}', fontsize=22, pad=8)
    plt.legend()
    plt.savefig(os.path.join(out_dir, f'{dist_path.split("/")[-1]}_number_dist.png'), bbox_inches='tight')
    fig.tight_layout()
    plt.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Use dr plotter to plot dr, fit it and find diffusion coefficient')
    parser.add_argument('-start', '--start', default=1000, type=int,
                        help='start fitting from')
    parser.add_argument('-path', '--path-to-results', default="../results", type=str,
                        help='path to results directory')
    args = parser.parse_args()

    start = args.start
    path_to_results = args.path_to_results

    for dirpath, dirnames, filenames in os.walk(os.path.join(path_to_results)):
        for filename in [f for f in filenames if "dr" in f and ".png" not in f]:
            file = os.path.join(dirpath, filename)
            process_dr(file, start)

    for dirpath, dirnames, filenames in os.walk(os.path.join(path_to_results)):
        for filename in [f for f in filenames if "size_distribution" in f and ".png" not in f]:
            file = os.path.join(dirpath, filename)
            process_size_distribution(file)

    for dirpath, dirnames, filenames in os.walk(os.path.join(path_to_results)):
        for filename in [f for f in filenames if "number_distribution" in f and ".png" not in f]:
            file = os.path.join(dirpath, filename)
            process_number_distribution(file)

