from utils import linear

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
    print(dr_path)
    out_dir = dr_path.split("/")[:-1]
    out_dir = os.path.join(*out_dir)
    dr = pd.read_csv(dr_path, sep='\t', index_col=0)

    print(dr)

    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111)
    ax.plot(dr.index, dr["dr"], 'b-', label='dr', markersize=3)
    ax.set_title(f'dr of time {dr_path.split("/")[-1]}', fontsize=22, pad=8)
    plt.legend()
    plt.savefig(os.path.join(out_dir, f'{dr_path.split("/")[-1]}_dr_t.png'), bbox_inches='tight')
    fig.tight_layout()
    plt.close()

    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111)
    ax.plot(dr.index, dr["ddr"], "o", color="red", label='ddr', markersize=3)
    popt, pcov = curve_fit(linear, dr.index[start:], dr["ddr"][start:])
    print(popt)
    ax.plot(dr.index[start:], linear(dr.index[start:], *popt), 'b-', label='fit: a=%f, b=%f' % tuple(popt))
    ax.set_title(f'ddr of time {dr_path.split("/")[-1]}', fontsize=22, pad=8)
    plt.legend()
    plt.savefig(os.path.join(out_dir, f'{dr_path.split("/")[-1]}_ddr_fit.png'), bbox_inches='tight')
    fig.tight_layout()
    plt.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-start', '--start', default=1000, type=int,
                        help='start fitting from')
    args = parser.parse_args()

    start = args.start

    for dirpath, dirnames, filenames in os.walk(os.path.join("../results")):
        for filename in [f for f in filenames if "dr" in f and ".png" not in f]:
            file = os.path.join(dirpath, filename)
            process_dr(file, start)
