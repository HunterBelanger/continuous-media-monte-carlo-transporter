/*=============================================================================*
 * Copyright (C) 2020, Commissariat à l'Energie Atomique
 *
 * Contributeur : Hunter Belanger (hunter.belanger@cea.fr)
 *
 * Ce logiciel est un programme informatique servant à faire des comparaisons
 * entre les méthodes de transport qui sont capable de traiter les milieux
 * continus avec la méthode Monte Carlo. Il résoud l'équation de Boltzmann
 * pour les particules neutres, à une vitesse et dans une dimension.
 *
 * Ce logiciel est régi par la licence CeCILL soumise au droit français et
 * respectant les principes de diffusion des logiciels libres. Vous pouvez
 * utiliser, modifier et/ou redistribuer ce programme sous les conditions
 * de la licence CeCILL telle que diffusée par le CEA, le CNRS et l'INRIA
 * sur le site "http://www.cecill.info".
 *
 * En contrepartie de l'accessibilité au code source et des droits de copie,
 * de modification et de redistribution accordés par cette licence, il n'est
 * offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,
 * seule une responsabilité restreinte pèse sur l'auteur du programme,  le
 * titulaire des droits patrimoniaux et les concédants successifs.
 *
 * A cet égard  l'attention de l'utilisateur est attirée sur les risques
 * associés au chargement,  à l'utilisation,  à la modification et/ou au
 * développement et à la reproduction du logiciel par l'utilisateur étant
 * donné sa spécificité de logiciel libre, qui peut le rendre complexe à
 * manipuler et qui le réserve donc à des développeurs et des professionnels
 * avertis possédant  des  connaissances  informatiques approfondies.  Les
 * utilisateurs sont donc invités à charger  et  tester  l'adéquation  du
 * logiciel à leurs besoins dans des conditions permettant d'assurer la
 * sécurité de leurs systèmes et ou de leurs données et, plus généralement,
 * à l'utiliser et l'exploiter dans les mêmes conditions de sécurité.
 *
 * Le fait que vous puissiez accéder à cet en-tête signifie que vous avez
 * pris connaissance de la licence CeCILL, et que vous en avez accepté les
 * termes.
 *============================================================================*/
#include "sub_programs.hpp"

#include <cmath>
#include "docopt.h"

static const char USAGE[] =
    R"(Continuous Media Monte Carlo Transporter

Usage:
  cmmct all [-p=<num> -b=<num> --output=<file>]
  cmmct trials <xs> <transporter> [-t=<num> -p=<num> -b=<num>]

Options:
  -h --help            Show this message
  -v --version         Show version
  -o --output=<file>   Filename for collision density (default output.txt)
  -t=<num>             Number of independent trials to run
  -p=<num>             Number of particles
  -b=<num>             Number of batches

Valid arguments for the cross section (xs) are:
  C : Constant
  LI: Linearly Increasing
  LD: Linearly Decreasing
  EI: Exponentially Increasing
  ED: Exponentially Decreasing
  SG: Sharp Gaussian
  BG: Broad Gaussian

Valid arguments for transporter are:
  SS    : Substeping
  DS    : Direct Sampling
  DT    : Delta Tracking
  RDT   : Meshed Delta Tracking
  NWDT  : Negative Weighted Delta Tracking
  RNWDT : Meshed Negative Weighted Delta Tracking
  CCTT  : Carter Tracking
  RCCTT : Meshed Carter Tracking
  IRCCTT: Improving Meshed Carter Tracking)";

static const char VERSION[] =
    R"(Continuous Media Monte Carlo Transporter 1.0.0
Written by Hunter Belanger (hunter.belanger@cea.fr)

Copyright (C) 2020, Commissariat à l'Energie Atomique
Released under the CeCILL license.)";

int main(int argc, const char** argv) {
  std::map<std::string, docopt::value> args =
      docopt::docopt(USAGE, {argv + 1, argv + argc}, true, VERSION);

  // Constants for XS functions and transport methods
  double A = std::sqrt(2.0 / M_PI);  // Amplitude fo Gaussians
  double z = 1.23;                   // Center of Gaussian XSs
  double s_b = 1.0;                  // Width of Broad Gaussian
  double s_s = 0.05;                 // Width of Sharp Gaussian
  double b = 0.1;                    // Background for Gaussian XSs
  double P = 0.85;                   // Esmp / Emaj ratio for NWDT
  double q = 0.3;                    // Probability of real collision for NWDT
  double P_c = 0.85;                 // Esmp / Emaj ratio for CCTT

  // If needed, open output file
  if (args["all"].asBool()) {
    if (args["--output"])
      file.open(args["--output"].asString());
    else
      file.open("output.txt");
  }

  if (args["-p"]) {
    int NP = std::stoi(args["-p"].asString());
    if (NP >= 1000)
      NPARTICLES = NP;
    else {
      std::cout << " Invalid number of particles. Using " << NPARTICLES
                << " instead.\n";
    }
  }
  if (args["-b"]) {
    int NB = std::stoi(args["-b"].asString());
    if (NB >= 10)
      NBATCHES = NB;
    else {
      std::cout << " Invalid number of batches. Using " << NBATCHES
                << " instead.\n";
    }
  }

  // Determine which program to run
  if (args["trials"].asBool()) {
    // Determine which XS to use
    std::unique_ptr<XS> crs;
    if (args["<xs>"].asString() == "C")
      crs = std::make_unique<Constant>(1.0);
    else if (args["<xs>"].asString() == "LI")
      crs = std::make_unique<Linear>(0.0, 2.0);
    else if (args["<xs>"].asString() == "LD")
      crs = std::make_unique<Linear>(2.0, 0.0);
    else if (args["<xs>"].asString() == "EI")
      crs = std::make_unique<Exponential>(0.1, 2.0);
    else if (args["<xs>"].asString() == "ED")
      crs = std::make_unique<Exponential>(1., -3.0);
    else if (args["<xs>"].asString() == "SG")
      crs = std::make_unique<Gaussian>(A, z, s_s, b);
    else if (args["<xs>"].asString() == "BG")
      crs = std::make_unique<Gaussian>(A, z, s_b, b);
    else {
      std::cout << " ERROR: Invalid cross section. Program closing...\n";
      file.close();
      return 1;
    }

    // Determine which transport method to use
    std::unique_ptr<Transporter> trnsprt;
    if (args["<transporter>"].asString() == "SS")
      trnsprt = std::make_unique<Substeper>(50);
    else if (args["<transporter>"].asString() == "DS")
      trnsprt = std::make_unique<Direct_Sampler>();
    else if (args["<transporter>"].asString() == "DT")
      trnsprt = std::make_unique<Delta_Tracker>();
    else if (args["<transporter>"].asString() == "RDT")
      trnsprt = std::make_unique<Regional_Delta_Tracker>();
    else if (args["<transporter>"].asString() == "NWDT")
      trnsprt = std::make_unique<Negative_Weighted_Delta_Tracker>(P, q);
    else if (args["<transporter>"].asString() == "RNWDT")
      trnsprt =
          std::make_unique<Regional_Negative_Weighted_Delta_Tracker>(P, q);
    else if (args["<transporter>"].asString() == "CCTT")
      trnsprt = std::make_unique<Carter_Transporter>(P_c);
    else if (args["<transporter>"].asString() == "RCCTT")
      trnsprt = std::make_unique<Regional_Carter_Transporter>(P_c);
    else if (args["<transporter>"].asString() == "IRCCTT")
      trnsprt = std::make_unique<Improving_Regional_Carter_Transporter>(P_c);
    else {
      std::cout << " ERROR: Invalid transporter. Program closing...\n";
      file.close();
      return 1;
    }

    // Get number of trials to conduct if set
    if (args["-t"]) {
      int provided_NTRIALS = std::stoi(args["-t"].asString());
      if (provided_NTRIALS > 0)
        NTRIALS = provided_NTRIALS;
      else {
        std::cout << " ERROR: Invalid number of trials. Will perform "
                  << NTRIALS << " trial...\n\n";
      }
    }

    std::cout << " NBATCHES = " << NBATCHES << "\n";

    independent_trials(trnsprt, crs);

  } else if (args["all"].asBool()) {
    if (args["--output"])
      file.open(args["--output"].asString());
    else
      file.open("output.txt");

    std::cout << " NBATCHES = " << NBATCHES << "\n";

    fixed_source_all();
  } else {
    std::cout << " ERROR: Invalid Command. Program closing...\n";
    file.close();
    return 1;
  }

  // Close output file
  file.close();

  return 0;
}
