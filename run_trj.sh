JOBS_COUNTER=0
MAX_CHILDREN=11
MY_PID=$$
results="results_8_15k"
data_path="data_8"

for xtc_file in $data_path/*.xtc;
do
    JOBS_COUNTER=$((`ps ax -Ao ppid | grep $MY_PID | wc -l`))
    while [ $JOBS_COUNTER -ge $MAX_CHILDREN ]
    do
        JOBS_COUNTER=$((`ps ax -Ao ppid | grep $MY_PID | wc -l`))
    done

    ./analysis.sh $xtc_file $results $data_path &
done

echo Finishing children ...
# wait for children here
while [ $JOBS_COUNTER -gt 1 ]
do
    JOBS_COUNTER=$((`ps ax -Ao ppid | grep $MY_PID | wc -l`))
    sleep 1
done

echo Done

source venv/bin/activate

python utils/dr_plotter.py -path "${results}"
