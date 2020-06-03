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
#include "cross_sections.hpp"
#include "transporter.hpp"

#include <cmath>
#include <memory>

void fixed_source_all() {
  // Initialize RNG seeds
  Seed_RNGs();

  // Constants for Gaussian XSs
  double A = std::sqrt(2.0 / M_PI);
  double z = 1.23;
  double s_b = 1.0;
  double s_s = 0.05;
  double b = 0.1;

  // Constants for transporters
  double P = 0.85;
  double q = 0.3;
  double P_c = 0.85;

  // Make Source
  std::vector<Particle> bank;
  for (int i = 0; i < NPARTICLES; i++) bank.push_back(Particle(0.0, 1.0, 1.0));

  // Iterate through all cross section types
  for (int x = 0; x < 6; x++) {
    std::unique_ptr<XS> crs;
    if (x == 0) {
      crs = std::make_unique<Linear>(0.0, 2.0);
      file << "#LI\n";
    } else if (x == 1) {
      crs = std::make_unique<Linear>(2.0, 0.0);
      file << "#LD\n";
    } else if (x == 2) {
      crs = std::make_unique<Exponential>(0.1, 2.0);
      file << "#EI\n";
    } else if (x == 3) {
      crs = std::make_unique<Exponential>(1.0, -3.0);
      file << "#ED\n";
    } else if (x == 4) {
      crs = std::make_unique<Gaussian>(A, z, s_s, b);
      file << "#SG\n";
    } else if (x == 5) {
      crs = std::make_unique<Gaussian>(A, z, s_b, b);
      file << "#BG\n";
    }

    // Iterate through all tracking methods
    for (int t = 0; t < 9; t++) {
      std::unique_ptr<Transporter> trnprt;
      if (t == 0) {
        trnprt = std::make_unique<Substeper>(100);
        file << "##SS\n";
      } else if (t == 1) {
        trnprt = std::make_unique<Direct_Sampler>();
        file << "##DS\n";
      } else if (t == 2) {
        trnprt = std::make_unique<Delta_Tracker>();
        file << "##DT\n";
      } else if (t == 3) {
        trnprt = std::make_unique<Regional_Delta_Tracker>();
        file << "##RDT\n";
      } else if (t == 4) {
        trnprt = std::make_unique<Negative_Weighted_Delta_Tracker>(P, q);
        file << "##NWDT\n";
      } else if (t == 5) {
        trnprt =
            std::make_unique<Regional_Negative_Weighted_Delta_Tracker>(P, q);
        file << "##RNDWDT\n";
      } else if (t == 6) {
        trnprt = std::make_unique<Carter_Transporter>(P_c);
        file << "##CCTT\n";
      } else if (t == 7) {
        trnprt = std::make_unique<Regional_Carter_Transporter>(P_c);
        file << "##RCCTT\n";
      } else if (t == 8) {
        trnprt = std::make_unique<Improving_Regional_Carter_Transporter>(P_c);
        file << "##IRCCTT\n";
      }

      // Do all batches in trial
      for (int b = 1; b <= NBATCHES; b++) {
        trnprt->Transport(bank, crs);
        Accumulate_Batch(b);
        Zero_Batch_Scores();
      }
      Output();
      Zero_All_Scores();
      crs->Reset_XSs();

    }  // All Transporters
    std::cout << " =======================================================\n";
  }  // All Cross Sections
}
