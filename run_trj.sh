JOBS_COUNTER=0
MAX_CHILDREN=11
MY_PID=$$

source venv/bin/activate

for xtc_file in data/*.xtc;
do
    JOBS_COUNTER=$((`ps ax -Ao ppid | grep $MY_PID | wc -l`))
    while [ $JOBS_COUNTER -ge $MAX_CHILDREN ]
    do
        JOBS_COUNTER=$((`ps ax -Ao ppid | grep $MY_PID | wc -l`))
    done

    ./analysis.sh $xtc_file &
done

echo Finishing children ...
# wait for children here
while [ $JOBS_COUNTER -gt 1 ]
do
    JOBS_COUNTER=$((`ps ax -Ao ppid | grep $MY_PID | wc -l`))
    sleep 1
done

python utils/dr_plotter.py

echo Done
