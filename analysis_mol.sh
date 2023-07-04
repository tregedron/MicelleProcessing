xtc_file=$1
echo $xtc_file
window=20000
shift=1

xtc_file=$(basename "${xtc_file%.xtc}")

source venv/bin/activate

python calculate_diff_mol.py -trj "data_heptane/${xtc_file}.xtc" -top "data_heptane/${xtc_file}.gro" -w "${window}" -s "${shift}" > "logs/heptane_${xtc_file}_${window}_${shift}"

