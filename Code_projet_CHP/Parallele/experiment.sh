#!/bin/bash
Np_start=2
Np_max=8
Nb_exp=10
charge_totale=25.3416

echo "1 $charge_totale 1 1" >> temps.txt

for ((Np = Np_start; Np <= Np_max; Np++)); do
    moy=0
    val=0
    for ((k = 0; k < Nb_exp; k++)); do
        val=$(mpiexec -n "$Np" ./run Datafiles/datafile.dat)

        printf "\n"
        printf "\n"
        echo "temps $Np $val"
        moy=$(python -c "print($val+$moy)")
        sleep 0.5
    done
    moy=$(python -c "print($moy/$Nb_exp)")
    speedup=$(python -c "print($charge_totale/$moy)")
    eff=$(python -c "print($charge_totale/$moy/$Np)")
    echo "$Np $moy $speedup $eff" >> temps.txt
done
printf "\n"
