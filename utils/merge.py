import argparse
import numpy as np
import pandas as pd
import os


def merge_two_drs(dr_path_1: str, dr_path_2: str, out_name="merged_dr"):
    '''
    usage example: python merge.py -dr1 ../results/100ns_NPT_1_8micelles/dr_5000_1 -dr2 ../results/100ns_NPT_1_8micelles/dr_5000_1_ext
    :param dr_path_1:
    :param dr_path_2:
    :param out_name:
    :return:
    '''
    dr1 = pd.read_csv(dr_path_1, sep='\t')
    dr2 = pd.read_csv(dr_path_2, sep='\t')
    print(dr1)
    print(dr2)
    dr1["dr"] = (dr1["dr"] + dr2["dr"])*0.5
    dr1["ddr"] = (dr1["ddr"] + dr2["ddr"])*0.5
    out_dir = dr_path_1.split("/")[:-1]
    out_dir = os.path.join(*out_dir, out_name)
    print(out_dir)
    print(dr1)
    dr1.to_csv(out_dir, sep='\t', index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='use merger to merge two drs of the same system')
    parser.add_argument('-dr1', '--dr_path_1', default=None, type=str,
                        help='dr1 path')
    parser.add_argument('-dr2', '--dr_path_2', default=None, type=str,
                        help='dr2 path')
    parser.add_argument('-out', '--out_path', default=None, type=str,
                        help='out path')
    args = parser.parse_args()

    merge_two_drs(args.dr_path_1, args.dr_path_2, args.out_path)

