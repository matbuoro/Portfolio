#!/bin/bash

#DIR="results/" # directory where to save results
#cd $DIR

#nSIMUL=10 # Nb simulations 100

# SIMULATION
echo "BEGINNING OF SIMULATIONS"
echo "PID du processus courant : $$"
STARTTIME=$(date +%s);

scn=( 103 )
#scn=( 103 303 503 703 903 1103)

for s in "${scn[@]}"
do

#STARTTIME=$(date +%s);
Rscript --vanilla DEMOGRAPHY.R $s &
Rscript --vanilla PHENOGENOTYPE.R $s &

echo "Extraction of scenario $s started!" 
#sleep 2 # wait x seconds before starting the next extraction
done # end loop 


#ENDTIME=$(date +%s);
#MINUTES=$(( ($ENDTIME - $STARTTIME) / 60 ));
#echo "Simulations successful! Duration: $MINUTES minutes" 
#echo "END OF SIMULATIONS"
