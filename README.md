# Projet Calcul parallele --- ENSEIRB-MATMECA
# Auteurs : FR Hammer, A Kien et M Praud

------------------------------
* 'main.cpp' : programme pour la resolution de l'equation de diffusion par decompostion de domaine
	* pour compiler, lancer 'make'
	* pour compiler en mode debug, lancer 'make DEBUG=1'
	* pour executer en sequentiel, lancer './exec 1'
	* pour executer en parallele, lancer 'mpiexec -n X ./exec 0' avec X le nombre de procs
* 'parameters.dat' contient les informations pour le systeme et la resolution, dans l'ordre suivant:
(valeurs par defaut)
	* Nx (6), Ny (6) : nombre de noeuds du maillage en x et en y
	* Lx (1), Ly (1) : taille du domaine rectangulaire
	* D (1) : coefficient de diffusion
	* dt (0.01) : pas de temps
	* mode (1) : choix pour les fonctions f,g,h (3 modes actuellement)
	* h_part (2) : parametre de recouvrement, ie nombre de lignes communes aux deux domaines
	* alpha (1) : coefficient pour la condition de Neumann
	* beta (1) : coefficient pour la condition de Dirichlet 
* Arguments facultatifs pour l'executable
	* h_part, parametre de recouvrement: nombre de lignes communes aux deux domaines
	* mode, pour le choix des fonctions f,g,h (mode=1,2,3)
	* Nt, nombre d'iterations en temps
------------------------------
Commande pour visualiser avec gnuplot les fichiers solutionMe.txt ecrits dans le dossier Results,
avec Me le numero d'un proc

* gnuplot
* splot for [i=0:X-1] 'solution'.i.'.txt'
------------------------------
