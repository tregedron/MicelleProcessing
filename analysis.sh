xtc_file=$1
results=$2
data_path=$3

window=15000
shift=1

xtc_file=$(basename "${xtc_file%.xtc}")

source venv/bin/activate

python main.py -trj "${data_path}/${xtc_file}.xtc" -top "${data_path}/${xtc_file}.gro" -w "${window}" -s "${shift}" -res "${results}" > "logs/${xtc_file}_${window}_${shift}"

