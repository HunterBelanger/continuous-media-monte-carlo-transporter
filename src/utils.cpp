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
#include "utils.hpp"

// Global Variable Declarations
std::vector<uint64_t> pcg_seeds;
std::ofstream file;
std::ofstream source_dists;

// Function Definitions
void Roulette(Particle& p, PCG& rng) {
  if (std::abs(p.wgt) < WGT_CUTOFF) {
    double P_kill = 1.0 - (std::abs(p.wgt) / WGT_SURVIVE);
    if (rng.rand() < P_kill)
      p.kill();
    else {
      if (p.wgt > 0.0)
        p.wgt = WGT_SURVIVE;
      else
        p.wgt = -WGT_SURVIVE;
    }
  }
}

void Output() {
  double n_batches = static_cast<double>(NBATCHES);

  avg_leakage_rate_variance /= (n_batches - 1.0);
  double avg_leakage_rate_error =
      std::sqrt(avg_leakage_rate_variance / n_batches);

  avg_coll_per_part_real_variance /= (n_batches - 1.0);
  double avg_c_pp_real_error =
      std::sqrt(avg_coll_per_part_real_variance / n_batches);

  avg_coll_per_part_all_variance /= (n_batches - 1.0);
  double avg_c_pp_all_error =
      std::sqrt(avg_coll_per_part_all_variance / n_batches);

  double leak_rel_err = avg_leakage_rate_error / avg_leakage_rate;
  double coll_pp_rel_err = avg_c_pp_real_error / avg_coll_per_part_real;
  double coll_pp_all_err = avg_c_pp_all_error / avg_coll_per_part_all;

  double FOM_leak = 1.0 / (n_xs_evals * leak_rel_err * leak_rel_err);
  double FOM_c_rel = 1.0 / (n_xs_evals * coll_pp_rel_err * coll_pp_rel_err);
  double FOM_c_all = 1.0 / (n_xs_evals * coll_pp_all_err * coll_pp_all_err);

  std::cout << std::fixed << std::setprecision(6);
  std::cout << " Real Coll Rate: " << avg_coll_per_part_real << " +/- "
            << avg_c_pp_real_error;
  std::cout << std::scientific << ", FOM_c_rel = " << FOM_c_rel << "\n"
            << std::fixed << std::setprecision(6);
  std::cout << "  All Coll Rate: " << avg_coll_per_part_all << " +/- "
            << avg_c_pp_all_error;
  std::cout << std::scientific << ", FOM_c_all = " << FOM_c_all << "\n"
            << std::fixed << std::setprecision(6);
  std::cout << "   Leakage Rate: " << avg_leakage_rate << " +/- "
            << avg_leakage_rate_error;
  std::cout << std::scientific << ", FOM_leak  = " << FOM_leak << "\n";
  std::cout << " NParticle Transported: " << n_particles_transported;
  std::cout << " , XS Evals: " << n_xs_evals << "\n\n";

  // -------------------------------------------------------------------------
  // File Output
  // -------------------------------------------------------------------------

  // First Write real flux
  for (int i = 0; i < N_FLUX_BINS - 1; i++) {
    file << avg_coll_dens_real[i][0] / dx_flux_bins << ",";
  }
  file << avg_coll_dens_real[N_FLUX_BINS - 1][0] / dx_flux_bins << "\n";
  // Write real flux error
  for (int i = 0; i < N_FLUX_BINS; i++) {
    double err = avg_coll_dens_real[i][1] / (n_batches - 1.0);
    err = std::sqrt(err / n_batches);
    err /= dx_flux_bins;
    if (i != N_FLUX_BINS - 1)
      file << err << ",";
    else
      file << err << "\n";
  }
  // Write real flux FOM
  for (int i = 0; i < N_FLUX_BINS; i++) {
    double err = avg_coll_dens_real[i][1] / (n_batches - 1.0);
    err = std::sqrt(err / n_batches);
    double rel_err = err / avg_coll_dens_real[i][0];
    double FOM = 1.0 / (n_xs_evals * rel_err * rel_err);
    if (i != N_FLUX_BINS - 1)
      file << FOM << ",";
    else
      file << FOM << "\n";
  }

  // First Write all flux
  for (int i = 0; i < N_FLUX_BINS; i++) {
    if (i != N_FLUX_BINS - 1)
      file << avg_coll_dens_all[i][0] / dx_flux_bins << ",";
    else
      file << avg_coll_dens_all[i][0] / dx_flux_bins << "\n";
  }
  // Write all flux error
  for (int i = 0; i < N_FLUX_BINS; i++) {
    double err = avg_coll_dens_all[i][1] / (n_batches - 1.0);
    err = std::sqrt(err / n_batches);
    err /= dx_flux_bins;
    if (i != N_FLUX_BINS - 1)
      file << err << ",";
    else
      file << err << "\n";
  }
  // Write all flux FOM
  for (int i = 0; i < N_FLUX_BINS; i++) {
    double err = avg_coll_dens_all[i][1] / (n_batches - 1.0);
    err = std::sqrt(err / n_batches);
    double rel_err = err / avg_coll_dens_all[i][0];
    double FOM = 1.0 / (n_xs_evals * rel_err * rel_err);
    if (i != N_FLUX_BINS - 1)
      file << FOM << ",";
    else
      file << FOM << "\n";
  }
  // file << "\n";
  file << std::flush;
}

void Seed_RNGs() {
  size_t nthreads;
#ifdef _OPENMP
  nthreads = omp_get_max_threads();
#else
  nthreads = 1;
#endif
  pcg_seeds.resize(nthreads);
  for (size_t i = 0; i < nthreads; i++) {
    uint64_t seed = i + 1;
    pcg_seeds[i] = seed;
  }
}

double Dist_to_Bin(const Particle& p) {
  double x_low = static_cast<double>(p.bin) * dx_xs_bins + X_MIN;
  double x_hi = static_cast<double>(p.bin + 1) * dx_xs_bins + X_MIN;

  if (p.u == -1.0)
    return p.x - x_low;
  else
    return x_hi - p.x;
}
