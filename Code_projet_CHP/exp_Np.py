import subprocess
import time
import numpy as np
Nx = 1200
Ny = 1200

file_name = "Temps/exp_Np_1200.dat"

Np_start = 1
Np_max = 24
Nb_exp = 3
fact_elarg =4#2.3

moy = []
incert = []  # equart type
speedup = []
i_speedup = []  # equart type
eff = []
i_eff = []  # equart type
val = []

# bloucle sur le nb de procs
for Np in range(Np_start, Np_max+1):
    val = []
    # on repete plusieurs fois l'exp
    for k in range(1, Nb_exp+1):
        exc = subprocess.run(["mpiexec", "-n", str(Np), "./run", "datafile.dat",str(Nx), str(Ny)], stdout=subprocess.PIPE, text=True)
        val.append(float(exc.stdout))
        print("Np %d/%d, it %d/%d, val %f" % (Np, Np_max, k, Nb_exp, val[-1]))
        time.sleep(0.1)
    if (Np == 1):
        charge_totale = np.array(val).mean()
    moy.append(np.array(val).mean())
    incert.append(np.array(val).std()*fact_elarg)
    speedup.append(charge_totale/moy[-1])
    i_speedup.append(speedup[-1] * np.sqrt((incert[-1]/moy[-1])**2+(incert[0]/moy[0])**2))
    eff.append(charge_totale/moy[-1]/Np)
    i_eff.append(eff[-1] * np.sqrt((incert[-1]/moy[-1])**2+(incert[0]/moy[0])**2))


    # ecriture fichier
    file = open(file_name, "w+")
    for Np in range(len(moy)):
        file.write("%d %f %f %f %f %f %f \n" %(Np+1, moy[Np], incert[Np], speedup[Np], i_speedup[Np], eff[Np], i_eff[Np]))
    file.close()
