# Compilateur
CXX = mpicxx

# Pour compiler en mode debug, lancez:
# make DEBUG=1
#####

DEBUG ?= 0
ifeq ($(DEBUG), 1)
	# Options en mode debug - La variable est DEBUG est definie comme vraie
	CXXFLAGS = -g -Wall -DDEBUG
else
	# Options en mode optimise - La variable DEBUG est definie comme fausse
	CXXFLAGS = -O2 -DNDEBUG
endif

# Nom de l'executable
main_PROG = exec

# Fichiers sources a compiler
SRC = charge.cpp update.cpp solveur.cpp matrix.cpp fonctions.cpp
main_SRC = $(SRC) main.cpp
main_OBJECTS = $(main_SRC:.cpp=.o)

# La commande complete : rassemble tous les .o pour l'edition de liens
$(main_PROG): message $(main_OBJECTS)
	$(CXX)  $(CXXFLAGS) $(main_OBJECTS) -o $(main_PROG)

all: $(main_PROG)

# Supprime l'executable, les fichiers binaires (.o)
# les fichiers temporaires de sauvegarde (~),
# et les fichiers textes (.txt)
clean :
	rm -f $(main_OBJECTS) $(main_PROG) *~ *.txt

# Regle commune pour compiler chaque .cpp en .o
# $@ signifie "la destination" (a gauche de la regle : le fichier.o)
# $< signifie "la source" (a droite de la regle : le fichier.cpp)
%.o: %.cpp %.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

.PHONY: message
message:
	@echo "******************************"
	@if [ $(DEBUG) -eq "1" ]; then echo "** Building in debug mode **";  else echo "** Building in release mode **" ; fi
	@echo "******************************"


#exec :
#	 g++ -std=c++11 -g -o exec main.cpp charge.cpp solveur.cpp matrix.cpp fonctions.cpp
