from utils import linear
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import argparse
import os
import math


def calc_viscosity(alpha, T=223):
    '''
    D = D_0 - kT/(6 * pi * eta) * 1/L
    [D] = A^2/ps = 10^(-20)/10^(-12) m^2/s
    [1/L] = 1/nm = 10^(9) m
    k = 1.38*10^(-23) J/K

    1.38*10^(-23) * T / (6 * pi * alpha * 10^(-17) ) = eta, ms/kg
    [eta] = mPa*s = 10^(-3) ms/kg

    :param alpha: parameter from fitting D of 1/L
    :return: viscosity of solution in mPa*s
    '''

    viscosity = 1.38 * 10**(-23) * T / (6 * math.pi * alpha * 10**(-17))
    viscosity = -1 * viscosity*1000
    if viscosity<0:
        print("You are currently breaking the laws of hydrodynamics... viscosity can't be negative")
    return viscosity


def process_DofL(DofL_path):
    print(DofL_path)
    out_dir = DofL_path.split("/")[:-1]
    out_dir = os.path.join(*out_dir)
    DofL = pd.read_csv(DofL_path, sep='\t')

    L_inv = [1/DofL["L"].values]

    DofL["D"] = DofL["6D"]/6

    print(DofL["L"].values)
    print(1 / DofL["L"].values)
    print(DofL["D"].values)

    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111)
    ax.plot(1/DofL["L"], DofL["D"], "o", color="red", label='D in data', markersize=3)
    popt, pcov = curve_fit(linear, 1/DofL["L"], DofL["D"])
    print(popt)
    viscosity = calc_viscosity(popt[0])

    ax.plot(1/DofL["L"], linear(1/DofL["L"], *popt), 'b-', label='fit: a=%f, b=%f' % tuple(popt))
    font_size = 22
    ax.set_title(f'D of 1/L', fontsize=30, pad=8)
    ax.set_ylabel("D, A^2/ps", fontsize=font_size)
    ax.set_xlabel("1/L, 1/nm", fontsize=font_size)

    text_visc = f'Viscosity = {viscosity:.3f} mPa*s'
    text_D = f"$D_0$ = {popt[1]*10**(-4):.3E} cm^2/s"
    plt.scatter([], [], color="w", alpha=0, label=text_visc)
    plt.scatter([], [], color="w", alpha=0, label=text_D)
    plt.legend()

    plt.savefig(os.path.join(out_dir, f'{DofL_path.split("/")[-1]}_fit.png'), bbox_inches='tight')
    fig.tight_layout()
    plt.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-path', '--path', default="../results/DofL", type=str,
                        help='path to D of L data')
    args = parser.parse_args()

    process_DofL(args.path)
