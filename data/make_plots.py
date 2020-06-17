#!/usr/bin/python3
"""
 Copyright (C) 2020, Commissariat à l'Energie Atomique
 
 Contributeur : Hunter Belanger (hunter.belanger@cea.fr)
 
 Ce logiciel est un programme informatique servant à faire des comparaisons
 entre les méthodes de transport qui sont capable de traiter les milieux
 continus avec la méthode Monte Carlo. Il résoud l'équation de Boltzmann
 pour les particules neutres, à une vitesse et dans une dimension.
 
 Ce logiciel est régi par la licence CeCILL soumise au droit français et
 respectant les principes de diffusion des logiciels libres. Vous pouvez
 utiliser, modifier et/ou redistribuer ce programme sous les conditions
 de la licence CeCILL telle que diffusée par le CEA, le CNRS et l'INRIA
 sur le site "http://www.cecill.info".
 
 En contrepartie de l'accessibilité au code source et des droits de copie,
 de modification et de redistribution accordés par cette licence, il n'est
 offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,
 seule une responsabilité restreinte pèse sur l'auteur du programme,  le
 titulaire des droits patrimoniaux et les concédants successifs.
 
 A cet égard  l'attention de l'utilisateur est attirée sur les risques
 associés au chargement,  à l'utilisation,  à la modification et/ou au
 développement et à la reproduction du logiciel par l'utilisateur étant
 donné sa spécificité de logiciel libre, qui peut le rendre complexe à
 manipuler et qui le réserve donc à des développeurs et des professionnels
 avertis possédant  des  connaissances  informatiques approfondies.  Les
 utilisateurs sont donc invités à charger  et  tester  l'adéquation  du
 logiciel à leurs besoins dans des conditions permettant d'assurer la
 sécurité de leurs systèmes et ou de leurs données et, plus généralement,
 à l'utiliser et l'exploiter dans les mêmes conditions de sécurité.
 
 Le fait que vous puissiez accéder à cet en-tête signifie que vous avez
 pris connaissance de la licence CeCILL, et que vous en avez accepté les
 termes.
"""
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({'font.size': 14})

file = open("collision_density.txt")

x = []
L = 2.0
Nx = 100
dx = L/Nx
for i in range(Nx):
    x.append(i*dx + 0.5*dx)

data = {}
cnt = 0
for line in file:
    if("##" in line): # New tracker type
        tracker = line.replace("##","")
        tracker = tracker.strip()
        data[xsname][tracker] = {}
        data[xsname][tracker]["real"] = {}
        data[xsname][tracker]["total"] = {}

        data[xsname][tracker]["real"]["coll"] = []
        data[xsname][tracker]["real"]["err"] = []
        data[xsname][tracker]["real"]["FOM"] = []

        data[xsname][tracker]["total"]["coll"] = []
        data[xsname][tracker]["total"]["err"] = []
        data[xsname][tracker]["total"]["FOM"] = []

    elif("#" in line): # Make dictionary for new XS type
        xsname = line.replace("#","")
        xsname = xsname.strip()
        data[xsname] = {}
    
    else:
        line = line.split(",")

        if (cnt < 3):
            estimator = "real"
            indx = cnt
            
        else:
            estimator = "total"
            indx = cnt - 3

        if(indx == 0):
            val = "coll"
        elif(indx == 1):
            val = "err"
        elif(indx == 2):
            val = "FOM"

        for elem in line:
            data[xsname][tracker][estimator][val].append(float(elem))

        cnt += 1
        if(cnt == 6):
            cnt = 0

file.close()

# All possible sets
#xss = ["LI", "LD", "EI", "ED", "SG", "BG"]
#trns = ["SS", "DS", "DT", "RDT", "NWDT", "RNWDT", "CCTT", "RCCTT", "IRCCTT"]

# Sets plotted in paper
xss = ["EI", "ED", "SG", "BG"]
trns = ["DT", "RDT", "NWDT", "RNWDT", "CCTT", "RCCTT"]


# Plot real estimator data
for xs in xss:
    # Plots collision density, dont care but leave for cure IDing
    for trn in trns:
        plt.errorbar(x, data[xs][trn]["real"]["coll"], yerr=data[xs][trn]["real"]["err"], label=trn)

    plt.legend()
    plt.show()
    
    # Plots FOM
    for trn in trns:
        if trn == "DT":
            plt.plot(x, data[xs][trn]["real"]["FOM"], color='blue', linestyle='-', label=trn)
        elif trn == "RDT":
            plt.plot(x, data[xs][trn]["real"]["FOM"], color='blue', linestyle='--', label=trn)
        elif trn == "NWDT":
            plt.plot(x, data[xs][trn]["real"]["FOM"], color='green', linestyle='-', label=trn)
        elif trn == "RNWDT":
            plt.plot(x, data[xs][trn]["real"]["FOM"], color='green', linestyle='--', label=trn)
        elif trn == "CCTT":
            plt.plot(x, data[xs][trn]["real"]["FOM"], color='red', linestyle='-', label=trn)
        else:
            plt.plot(x, data[xs][trn]["real"]["FOM"], color='red', linestyle='--', label=trn)
        
    
    # Only put key in exponentially decreasing plot
    if(xs == "ED"):
        leg = plt.legend(ncol = 2)
        leg.draggable()

    plt.xlabel("Position [cm]")
    plt.ylabel("FOM [Arb. Units]")
    plt.xlim([0.,2.])
    #plt.title(xs + " FOM")
    plt.tight_layout()
    plt.show()

# Plot total estimator data
for xs in xss:
    for trn in trns:
        plt.errorbar(x, data[xs][trn]["total"]["coll"], yerr=data[xs][trn]["total"]["err"], label=trn)

    plt.legend()
    plt.show()

    # Plots FOM
    for trn in trns:
        if trn == "DT":
            plt.plot(x, data[xs][trn]["total"]["FOM"], color='blue', linestyle='-', label=trn)
        elif trn == "RDT":
            plt.plot(x, data[xs][trn]["total"]["FOM"], color='blue', linestyle='--', label=trn)
        elif trn == "NWDT":
            plt.plot(x, data[xs][trn]["total"]["FOM"], color='green', linestyle='-', label=trn)
        elif trn == "RNWDT":
            plt.plot(x, data[xs][trn]["total"]["FOM"], color='green', linestyle='--', label=trn)
        elif trn == "CCTT":
            plt.plot(x, data[xs][trn]["total"]["FOM"], color='red', linestyle='-', label=trn)
        else:
            plt.plot(x, data[xs][trn]["total"]["FOM"], color='red', linestyle='--', label=trn)
    
    if(xs == "ED"):
        leg = plt.legend()
        leg = plt.legend(ncol = 2)
        leg.draggable()

    plt.xlabel("Position [cm]")
    plt.ylabel("FOM [Arb. Units]")
    plt.xlim([0.,2.])
    #plt.title(xs + " FOM")
    plt.tight_layout()
    plt.show()


