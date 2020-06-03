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
#include "tallies.hpp"

// Tallies for current batch
double batch_escaped;          // Number of particles which escaped
double batch_real_collisions;  // Number of collisions with real estimator
double batch_all_collisions;   // Number of collisions with all estimator

double batch_real_coll_dens[N_FLUX_BINS];
double batch_all_coll_dens[N_FLUX_BINS];

// Problem tallies
double avg_leakage_rate;
double avg_leakage_rate_variance;

double avg_coll_per_part_real;
double avg_coll_per_part_real_variance;

double avg_coll_per_part_all;
double avg_coll_per_part_all_variance;

double avg_coll_dens_real[N_FLUX_BINS][2];  // [0] is avg, [1] is variance
double avg_coll_dens_all[N_FLUX_BINS][2];

// Other tally variables
int n_particles_transported;
int n_sign_changes_in_weight;
double n_xs_evals;

void Score_Escape(double w) {
#pragma omp atomic
  batch_escaped += w;
}

void Score_Collision_Real(double w, double x) {
#pragma omp atomic
  batch_real_collisions += w;

  int bin = std::floor(x / dx_flux_bins);
  if (bin >= N_FLUX_BINS) bin = N_FLUX_BINS - 1;

#pragma omp atomic
  batch_real_coll_dens[bin] += w;
}

void Score_Collision_All(double w, double x) {
#pragma omp atomic
  batch_all_collisions += w;

  int bin = std::floor(x / dx_flux_bins);
  if (bin >= N_FLUX_BINS) bin = N_FLUX_BINS - 1;

#pragma omp atomic
  batch_all_coll_dens[bin] += w;
}

void Accumulate_Batch(int batch) {
  double n = static_cast<double>(batch);
  double batch_leak_rate = batch_escaped / static_cast<double>(NPARTICLES);
  double batch_coll_per_part_real =
      batch_real_collisions / static_cast<double>(NPARTICLES);
  double batch_coll_per_part_all =
      batch_all_collisions / static_cast<double>(NPARTICLES);

  double M_leakage = (batch_leak_rate - avg_leakage_rate);
  double M_n_coll_real = (batch_coll_per_part_real - avg_coll_per_part_real);
  double M_n_coll_all = (batch_coll_per_part_all - avg_coll_per_part_all);

  avg_leakage_rate += (batch_leak_rate - avg_leakage_rate) / n;
  avg_coll_per_part_real +=
      (batch_coll_per_part_real - avg_coll_per_part_real) / n;
  avg_coll_per_part_all +=
      (batch_coll_per_part_all - avg_coll_per_part_all) / n;

  avg_leakage_rate_variance += (batch_leak_rate - avg_leakage_rate) * M_leakage;
  avg_coll_per_part_real_variance +=
      (batch_coll_per_part_real - avg_coll_per_part_real) * M_n_coll_real;
  avg_coll_per_part_all_variance +=
      (batch_coll_per_part_all - avg_coll_per_part_all) * M_n_coll_all;

  // Collision densities
  for (int b = 0; b < N_FLUX_BINS; b++) {
    double M_real = (batch_real_coll_dens[b] - avg_coll_dens_real[b][0]);
    double M_all = (batch_all_coll_dens[b] - avg_coll_dens_all[b][0]);

    avg_coll_dens_real[b][0] +=
        (batch_real_coll_dens[b] - avg_coll_dens_real[b][0]) / n;
    avg_coll_dens_all[b][0] +=
        (batch_all_coll_dens[b] - avg_coll_dens_all[b][0]) / n;

    avg_coll_dens_real[b][1] +=
        (batch_real_coll_dens[b] - avg_coll_dens_real[b][0]) * M_real;
    avg_coll_dens_all[b][1] +=
        (batch_all_coll_dens[b] - avg_coll_dens_all[b][0]) * M_all;
  }
}

void Zero_Batch_Scores() {
  batch_escaped = 0.0;
  batch_real_collisions = 0.0;
  batch_all_collisions = 0.0;

  for (int b = 0; b < N_FLUX_BINS; b++) {
    batch_real_coll_dens[b] = 0.0;
    batch_all_coll_dens[b] = 0.0;
  }
}

void Zero_All_Scores() {
  avg_leakage_rate = 0.0;
  avg_leakage_rate_variance = 0.0;

  avg_coll_per_part_real = 0.0;
  avg_coll_per_part_real_variance = 0.0;

  avg_coll_per_part_all = 0.0;
  avg_coll_per_part_all_variance = 0.0;

  for (int i = 0; i < N_FLUX_BINS; i++) {
    avg_coll_dens_all[i][0] = 0.0;
    avg_coll_dens_all[i][1] = 0.0;
    avg_coll_dens_real[i][0] = 0.0;
    avg_coll_dens_real[i][1] = 0.0;
  }

  Zero_Batch_Scores();

  // Zero other tallies
  n_particles_transported = 0;
  n_sign_changes_in_weight = 0;
  n_xs_evals = 0.0;
}
