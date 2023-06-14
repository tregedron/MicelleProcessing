xtc_file=$1

window=10000
shift=1

xtc_file=$(basename "${xtc_file%.xtc}")

source venv/bin/activate

python main.py -trj "data/${xtc_file}.xtc" -top "data/${xtc_file}.gro" -w "${window}" -s "${shift}" > "logs/${xtc_file}_${window}_${shift}"

