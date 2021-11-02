import subprocess
import time
import numpy as np
file_name = "Temps/exp_matsize_log.dat"

N_start = 30
N_max = 20000
Np = 20
Nb_exp = 5

# Déclarations

moy = []
incert = []  # equart typemake 
val = []
eff = []
i_eff = []  # equart type

Ns = np.arange(N_start,N_max,100)
Ns = np.linspace(N_start,N_max,200)
Ns = np.concatenate((np.linspace(30, 1000,50),np.logspace(3,5,50)))
Ns = Ns
# bloucle sur le nb de procs
for N in Ns:
    Nx = N
    Ny = N
    val = []
    # on repete plusieurs fois l'exp
    for k in range(1, Nb_exp+1):
        exc = subprocess.run(["mpiexec", "-n", str(Np), "./run", "datafile.dat",
                             str(Nx), str(Ny)], stdout=subprocess.PIPE, text=True)
        val.append(float(exc.stdout))
        print("N %d/%d, it %d/%d, val %f" % (N, N_max, k, Nb_exp, val[-1]))
        time.sleep(0.1)
    moy.append(np.array(val).mean())
    incert.append(np.array(val).std()*2/np.sqrt(Nb_exp))
    # eff.append(moy[0]/moy[-1]*Nx**2)
    # i_eff.append(eff[-1] * np.sqrt((incert[-1]/moy[-1])**2+(incert[0]/moy[0])**2))


    # ecriture fichier
    file = open(file_name, "w+")
    for k in range(len(moy)):
        file.write("%d %f %f \n" %
                   (Ns[k], moy[k], incert[k]))
    file.close()
