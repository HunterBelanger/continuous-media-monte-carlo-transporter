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
#include "transporter.hpp"

void Delta_Tracker::Transport(const std::vector<Particle>& bank,
                              const std::unique_ptr<XS>& xs) {
  int cnts_sum = 0;

#pragma omp parallel
  {
    PCG rng;
    int thread_id;
#ifdef _OPENMP
    thread_id = omp_get_thread_num();
#else
    thread_id = 0;
#endif
    uint64_t pcg_seed;
#pragma omp atomic read
    pcg_seed = pcg_seeds[thread_id];
    rng.seed(pcg_seed);

#pragma omp for schedule(dynamic)
    for (int n = 0; n < static_cast<int>(bank.size()); n++) {
#pragma omp atomic
      n_particles_transported++;

      bool virtual_collision;
      Particle p = bank[n];
      while (p.alive) {
        virtual_collision = true;
        while (virtual_collision) {
          double d = -std::log(rng.rand()) / (xs->Emax);
          p.move(d);
          if ((p.x >= X_MAX) or (p.x <= X_MIN)) {
            Score_Escape(p.wgt);
            virtual_collision = false;
            p.kill();
#pragma omp atomic
            cnts_sum += p.xs_evals;
          } else {
            p.xs_eval();
            if (rng.rand() < (xs->E(p.x) / xs->Emax)) {
              virtual_collision = false;
            }
            // Score every collision
            double score = p.wgt * xs->E(p.x) / xs->Emax;
            Score_Collision_All(score, p.x);
          }
        }

        if (p.alive) {
          // Score real collision
          Score_Collision_Real(p.wgt, p.x);

          // Implicit capture
          p.wgt *= 1.0 - P_abs;  // Equivalent to 1.0 - (Ea/Et)

          // Scatter
          if (rng.rand() > P_straight_ahead) p.turn();

          // Russian Roulette
          Roulette(p, rng);
          if (p.alive == false) {
#pragma omp atomic
            cnts_sum += p.xs_evals;
          }
        }  // If alive for real collision
      }    // While alive
    }      // For all particles
#pragma omp critical
    { pcg_seeds[thread_id] = rng.get_seed(); }

  }  // Parallel
  n_xs_evals += cnts_sum;
}

void Regional_Delta_Tracker::Transport(const std::vector<Particle>& bank,
                                       const std::unique_ptr<XS>& xs) {
  int cnts_sum = 0;

#pragma omp parallel
  {
    PCG rng;
    int thread_id;
#ifdef _OPENMP
    thread_id = omp_get_thread_num();
#else
    thread_id = 0;
#endif
    uint64_t pcg_seed;
#pragma omp atomic read
    pcg_seed = pcg_seeds[thread_id];
    rng.seed(pcg_seed);

#pragma omp for schedule(dynamic)
    for (int n = 0; n < static_cast<int>(bank.size()); n++) {
#pragma omp atomic
      n_particles_transported++;

      bool virtual_collision = true;
      Particle p = bank[n];
      double Emax = xs->E_bin_max[p.bin];

      double d, d_bin;
      while (p.alive) {
        virtual_collision = true;
        while (virtual_collision) {
          d_bin = Dist_to_Bin(p);
          d = -std::log(rng.rand()) / Emax;
          if (d_bin < d) {
            d = d_bin + 1e-6;
            p.move(d);
            if ((p.x >= X_MAX) or (p.x <= X_MIN)) {
              Score_Escape(p.wgt);
              virtual_collision = false;
              p.kill();

// Score other tallies
#pragma omp atomic
              cnts_sum += p.xs_evals;
            } else {
              Emax = xs->E_bin_max[p.bin];
            }
          } else {
            p.xs_eval();
            p.move(d);
            double xi = rng.rand();
            double Pr = xs->E(p.x) / Emax;
            if (xi < Pr) {
              virtual_collision = false;
            }

            double score = p.wgt * xs->E(p.x) / Emax;
            Score_Collision_All(score, p.x);
          }
        }  // While virtual

        if (p.alive) {
          // Score real collision
          Score_Collision_Real(p.wgt, p.x);

          // Implicit capture
          p.wgt *= 1.0 - P_abs;

          // Scatter
          if (rng.rand() > P_straight_ahead) p.turn();

          // Russian Roulette
          Roulette(p, rng);
          if (p.alive == false) {
#pragma omp atomic
            cnts_sum += p.xs_evals;
          }
        }
      }  // While alive
    }    // For all particles
#pragma omp atomic write
    pcg_seeds[thread_id] = rng.get_seed();

  }  // Parallel
  n_xs_evals += cnts_sum;
}

void Negative_Weighted_Delta_Tracker::Transport(
    const std::vector<Particle>& bank, const std::unique_ptr<XS>& xs) {
  int cnts_sum = 0;
  double sign_change = 0.0;

  std::vector<Particle> Bank = bank;  // Copy for re-writing after splits
  int n_particles = static_cast<int>(Bank.size());

  while (n_particles > 0) {
#pragma omp parallel
    {
      PCG rng;
      int thread_id;
#ifdef _OPENMP
      thread_id = omp_get_thread_num();
#else
      thread_id = 0;
#endif
      uint64_t pcg_seed;
#pragma omp atomic read
      pcg_seed = pcg_seeds[thread_id];
      rng.seed(pcg_seed);

      std::vector<Particle> this_thread_Splits;

#pragma omp for schedule(dynamic)
      for (int n = 0; n < n_particles; n++) {
#pragma omp atomic
        n_particles_transported++;

        double Esamp = P * (xs->Emax);
        Particle p = Bank[n];
        while (p.alive) {
          double d = -std::log(rng.rand()) / Esamp;
          p.move(d);

          if ((p.x >= X_MAX) or (p.x <= X_MIN)) {
            Score_Escape(p.wgt);

            p.kill();
#pragma omp atomic
            cnts_sum += p.xs_evals;
          } else if (rng.rand() < q) {
            // Collision is real, addjust weight
            p.wgt *= xs->E(p.x) / (Esamp * q);
            double score = p.wgt * q;
            Score_Collision_All(score, p.x);
            Score_Collision_Real(p.wgt, p.x);
            p.xs_eval();

            // Implicit capture
            p.wgt *= 1.0 - P_abs;

            // Scatter
            if (rng.rand() > P_straight_ahead) p.turn();

            // Russian Roulette
            Roulette(p, rng);
            if (p.alive == false) {
#pragma omp atomic
              cnts_sum += p.xs_evals;
            }
          } else {  // Collision is virtual
            double dw = ((1.0 - (xs->E(p.x) / Esamp)) / (1.0 - q));
            if (dw < 0.0) {
#pragma omp atomic
              sign_change += 1.0;
            }
            double score = p.wgt * (xs->E(p.x) / Esamp);
            Score_Collision_All(score, p.x);
            p.wgt = p.wgt * dw;
            p.xs_eval();
          }

          // Split if needed
          if (p.alive and (std::abs(p.wgt) >= WGT_SPLIT)) {
            double n_new = std::round(std::abs(p.wgt));
            p.wgt /= n_new;
            for (int j = 0; j < static_cast<int>(n_new - 1); j++) {
              this_thread_Splits.push_back(Particle(p.x, p.u, p.wgt));
            }
          }  // split
        }    // While alive
      }      // For all particles
#pragma omp atomic write
      pcg_seeds[thread_id] = rng.get_seed();

// Clear Bank to accept splits
#pragma omp single
      { Bank.clear(); }
// Add thread splits to Bank
#pragma omp barrier
#pragma omp critical
      {
        Bank.insert(std::end(Bank), std::begin(this_thread_Splits),
                    std::end(this_thread_Splits));
      }
    }  // Parallel
    n_particles = static_cast<int>(Bank.size());

  }  // While still split particles

  n_xs_evals += cnts_sum;
}

void Regional_Negative_Weighted_Delta_Tracker::Transport(
    const std::vector<Particle>& bank, const std::unique_ptr<XS>& xs) {
  // Adjust sampling XSs in xs object by factor P
  xs->Adust_Smp_XSs(P);

  int cnts_sum = 0;
  double sign_change = 0.0;

  std::vector<Particle> Bank = bank;
  int n_particles = static_cast<int>(Bank.size());

  while (n_particles > 0) {
#pragma omp parallel
    {
      PCG rng;
      int thread_id;
#ifdef _OPENMP
      thread_id = omp_get_thread_num();
#else
      thread_id = 0;
#endif
      uint64_t pcg_seed;
#pragma omp atomic read
      pcg_seed = pcg_seeds[thread_id];
      rng.seed(pcg_seed);

      std::vector<Particle> this_thread_Splits;

#pragma omp for schedule(dynamic)
      for (int n = 0; n < n_particles; n++) {
#pragma omp atomic
        n_particles_transported++;

        Particle p = Bank[n];
        double Esamp = xs->E_bin_smp[p.bin];
        double d, d_bin;
        while (p.alive) {
          d_bin = Dist_to_Bin(p);
          d = -std::log(rng.rand()) / Esamp;
          if (d_bin < d) {
            d = d_bin + 1e-6;
            p.move(d);
            if ((p.x >= X_MAX) or (p.x <= X_MIN)) {
              Score_Escape(p.wgt);
              p.kill();

#pragma omp atomic
              cnts_sum += p.xs_evals;
            } else {
              Esamp = xs->E_bin_smp[p.bin];
            }
          } else {
            p.move(d);
            if (rng.rand() < q) {  // Real collision
              // update weight
              p.wgt *= (xs->E(p.x) / (Esamp * q));
              double score = p.wgt * (q);
              Score_Collision_All(score, p.x);
              Score_Collision_Real(p.wgt, p.x);
              p.xs_eval();

              // Implicit capture
              p.wgt *= 1.0 - P_abs;

              // Scatter
              if (rng.rand() > P_straight_ahead) p.turn();

              // Russian Roulette
              Roulette(p, rng);
              if (p.alive == false) {
#pragma omp atomic
                cnts_sum += p.xs_evals;
              }
            } else {
              p.xs_eval();
              double dw = (1.0 - (xs->E(p.x) / Esamp)) / (1.0 - q);
              if (dw < 0.0) {
#pragma omp atomic
                sign_change += 1.0;
              }
              double score = p.wgt * (xs->E(p.x) / Esamp);
              Score_Collision_All(score, p.x);
              p.wgt = p.wgt * dw;
            }

            // Split if needed
            if (p.alive and (std::abs(p.wgt) >= WGT_SPLIT)) {
              double n_new = std::round(std::abs(p.wgt));
              p.wgt /= n_new;
              for (int j = 0; j < static_cast<int>(n_new - 1); j++) {
                Particle p_daughter(p.x, p.u, p.wgt);
                p_daughter.bin = p.bin;
                this_thread_Splits.push_back(p_daughter);
              }
            }  // split
          }
        }  // While alive
      }    // For all particles
#pragma omp atomic write
      pcg_seeds[thread_id] = rng.get_seed();

// Clear bank for new particles
#pragma omp single
      { Bank.clear(); }

// Add thread splits to Bank
#pragma omp barrier
#pragma omp critical
      {
        Bank.insert(std::end(Bank), std::begin(this_thread_Splits),
                    std::end(this_thread_Splits));
      }
    }  // Parallel

    n_particles = static_cast<int>(Bank.size());

  }  // While still split particles

  n_xs_evals += cnts_sum;
}

void Carter_Transporter::Transport(const std::vector<Particle>& bank,
                                   const std::unique_ptr<XS>& xs) {
  int cnts_sum = 0;
  double sign_change = 0.0;

  // Particle bank vectors
  std::vector<Particle> Bank = bank;
  int n_particles = static_cast<int>(Bank.size());

  while (n_particles > 0) {
// std::cout << " nparticles = " << n_particles << "\n";
#pragma omp parallel
    {
      PCG rng;
      int thread_id;
#ifdef _OPENMP
      thread_id = omp_get_thread_num();
#else
      thread_id = 0;
#endif
      uint64_t pcg_seed;
#pragma omp atomic read
      pcg_seed = pcg_seeds[thread_id];
      rng.seed(pcg_seed);

      std::vector<Particle> this_thread_Splits;

#pragma omp for schedule(dynamic)
      for (int n = 0; n < n_particles; n++) {
#pragma omp atomic
        n_particles_transported++;

        Particle p = Bank[n];
        double Esmp = P * xs->Emax;
        bool real_collision = false;

        while (p.alive) {
          double d = -std::log(rng.rand()) / Esmp;
          p.move(d);
          real_collision = false;

          // Fist check for leak
          if ((p.x >= X_MAX) or (p.x <= X_MIN)) {
            Score_Escape(p.wgt);
            p.kill();
#pragma omp atomic
            cnts_sum += p.xs_evals;
          } else {
            double E_tot = xs->E(p.x);
            p.xs_eval();
            if (E_tot > Esmp) {  // First negative branch
              double D = E_tot / (2 * E_tot - Esmp);
              double F = (E_tot / (D * Esmp));
              double score = p.wgt * E_tot / Esmp;
              Score_Collision_All(score, p.x);

              // if(rand(rng) < D_alpha) {
              if (rng.rand() < D) {
                real_collision = true;
                p.wgt *= F;
                // wgt *=  F*(D/D_alpha);
              } else {
                p.wgt *= -F;
// wgt *= -F*((1. - D)/(1. - D_alpha));
#pragma omp atomic
                sign_change += 1.0;
              }

            } else {  // Delta tracking branch
              double P_real = E_tot / Esmp;
              if (rng.rand() < P_real) {
                real_collision = true;
              }
              double score = p.wgt * E_tot / Esmp;
              Score_Collision_All(score, p.x);
            }

            if (real_collision) {
              // Score real collision
              Score_Collision_Real(p.wgt, p.x);

              // Implicit caputure
              p.wgt *= 1.0 - P_abs;

              // Scatter
              if (rng.rand() > P_straight_ahead) p.turn();

              // Russian Roulette
              Roulette(p, rng);
              if (p.alive == false) {
#pragma omp atomic
                cnts_sum += p.xs_evals;
              }

            }  // End real coll.
          }

          // Split if needed
          if (p.alive and (std::abs(p.wgt) >= WGT_SPLIT)) {
            double n_new = std::floor(std::abs(p.wgt));
            p.wgt /= n_new;
            for (int j = 0; j < static_cast<int>(n_new - 1); j++) {
              Particle p_daughter(p.x, p.u, p.wgt);
              this_thread_Splits.push_back(p_daughter);
            }
          }  // split

        }  // While alive
      }    // For all particles
#pragma omp atomic write
      pcg_seeds[thread_id] = rng.get_seed();

#pragma omp single
      { Bank.clear(); }

#pragma omp barrier
#pragma omp critical
      {
        Bank.insert(std::end(Bank), std::begin(this_thread_Splits),
                    std::end(this_thread_Splits));
      }
    }  // Parallel

    n_particles = static_cast<int>(Bank.size());

  }  // While still particles

  n_xs_evals += cnts_sum;
}

void Regional_Carter_Transporter::Transport(const std::vector<Particle>& bank,
                                            const std::unique_ptr<XS>& xs) {
  int cnts_sum = 0;
  double sign_change = 0.0;

  // Particle bank vectors
  std::vector<Particle> Bank = bank;
  int n_particles = static_cast<int>(Bank.size());

  while (n_particles > 0) {
#pragma omp parallel
    {
      PCG rng;
      int thread_id;
#ifdef _OPENMP
      thread_id = omp_get_thread_num();
#else
      thread_id = 0;
#endif
      uint64_t pcg_seed;
#pragma omp atomic read
      pcg_seed = pcg_seeds[thread_id];
      rng.seed(pcg_seed);

      std::vector<Particle> this_thread_Splits;

#pragma omp for schedule(dynamic)
      for (int n = 0; n < n_particles; n++) {
#pragma omp atomic
        n_particles_transported++;

        Particle p = Bank[n];
        double d_bin;
        double Esmp = P * xs->E_bin_smp[p.bin];
        bool real_collision = false;

        while (p.alive) {
          double d = -std::log(rng.rand()) / Esmp;
          d_bin = Dist_to_Bin(p);
          real_collision = false;

          if (d_bin < d) {
            p.move(d_bin + 1e-6);
            if ((p.x >= X_MAX) or (p.x <= X_MIN)) {
              p.kill();
              Score_Escape(p.wgt);
#pragma omp atomic
              cnts_sum += p.xs_evals;
            } else {
              Esmp = P * xs->E_bin_smp[p.bin];
            }

          } else {
            p.move(d);
            // Fist check for leak
            if ((p.x >= X_MAX) or (p.x <= X_MIN)) {
              p.kill();
              Score_Escape(p.wgt);
#pragma omp atomic
              cnts_sum += p.xs_evals;
            } else {
              double E_tot = xs->E(p.x);
              p.xs_eval();
              if (E_tot > Esmp) {  // First negative branch
                double D = E_tot / (2 * E_tot - Esmp);
                double F = E_tot / (D * Esmp);
                double score = p.wgt * E_tot / Esmp;
                Score_Collision_All(score, p.x);
                p.wgt *= F;
                if (rng.rand() < D) {
                  real_collision = true;
                } else {
                  p.wgt *= -1.0;
#pragma omp atomic
                  sign_change += 1.0;
                }

              } else {  // Delta tracking branch
                double P_real = E_tot / Esmp;
                double score = p.wgt * E_tot / Esmp;
                Score_Collision_All(score, p.x);
                if (rng.rand() < P_real) {
                  real_collision = true;
                }
              }

              if (real_collision) {
                // Score real collision
                Score_Collision_Real(p.wgt, p.x);

                // Implicit capture
                p.wgt *= 1.0 - P_abs;

                // Scatter
                if (rng.rand() > P_straight_ahead) p.turn();

                // Russian Roulette
                Roulette(p, rng);
                if (p.alive == false) {
#pragma omp atomic
                  cnts_sum += p.xs_evals;
                }
              }

              // Split if needed
              if (p.alive and (std::abs(p.wgt) >= WGT_SPLIT)) {
                double n_new = std::round(std::abs(p.wgt));
                p.wgt /= n_new;
                for (int j = 0; j < static_cast<int>(n_new - 1); j++) {
                  Particle p_daughter(p.x, p.u, p.wgt);
                  p_daughter.bin = p.bin;
                  this_thread_Splits.push_back(p_daughter);
                }
              }  // split
            }
          }
        }  // While alive
      }    // For all particles
#pragma omp atomic write
      pcg_seeds[thread_id] = rng.get_seed();

#pragma omp single
      { Bank.clear(); }

#pragma omp barrier
#pragma omp critical
      {
        Bank.insert(std::end(Bank), std::begin(this_thread_Splits),
                    std::end(this_thread_Splits));
      }
    }  // Parallel

    n_particles = static_cast<int>(Bank.size());

  }  // While split particles

  n_xs_evals += cnts_sum;
}

void Improving_Regional_Carter_Transporter::Transport(
    const std::vector<Particle>& bank, const std::unique_ptr<XS>& xs) {
  // Adjust XSs (only works first time, so updates are saved which is good)
  xs->Adust_Smp_XSs(P);

  int cnts_sum = 0;
  double sign_change = 0.0;

  // Particle bank vectors
  std::vector<Particle> Bank = bank;
  int n_particles = static_cast<int>(Bank.size());

  while (n_particles > 0) {
#pragma omp parallel
    {
      PCG rng;
      int thread_id;
#ifdef _OPENMP
      thread_id = omp_get_thread_num();
#else
      thread_id = 0;
#endif
      uint64_t pcg_seed;
#pragma omp atomic read
      pcg_seed = pcg_seeds[thread_id];
      rng.seed(pcg_seed);

      std::vector<Particle> this_thread_Splits;

#pragma omp for schedule(dynamic)
      for (int n = 0; n < n_particles; n++) {
#pragma omp atomic
        n_particles_transported++;

        Particle p = Bank[n];
        double d_bin;
        double Esmp = xs->E_bin_smp[p.bin];
        bool real_collision = false;

        while (p.alive) {
          double d = -std::log(rng.rand()) / Esmp;
          d_bin = Dist_to_Bin(p);
          real_collision = false;

          if (d_bin < d) {
            p.move(d_bin + 1e-6);
            if ((p.x >= X_MAX) or (p.x <= X_MIN)) {
              p.kill();
              Score_Escape(p.wgt);
#pragma omp atomic
              cnts_sum += p.xs_evals;
            } else {
              Esmp = xs->E_bin_smp[p.bin];
            }

          } else {
            p.move(d);
            // Fist check for leak
            if ((p.x >= X_MAX) or (p.x <= X_MIN)) {
              p.kill();
              Score_Escape(p.wgt);
#pragma omp atomic
              cnts_sum += p.xs_evals;
            } else {
              double E_tot = xs->E(p.x);
              p.xs_eval();
              if (E_tot > Esmp) {  // First negative branch
                double D = E_tot / (2 * E_tot - Esmp);
                double F = E_tot / (D * Esmp);
                double score = p.wgt * E_tot / Esmp;
                Score_Collision_All(score, p.x);
                p.wgt *= F;
                if (rng.rand() < D) {
                  real_collision = true;
                } else {
                  p.wgt *= -1.0;
#pragma omp atomic
                  sign_change += 1.0;
                }

              } else {  // Delta tracking branch
                double P_real = E_tot / Esmp;
                double score = p.wgt * E_tot / Esmp;
                Score_Collision_All(score, p.x);
                if (rng.rand() < P_real) {
                  real_collision = true;
                }
              }

              if (real_collision) {
                // Score real collision
                Score_Collision_Real(p.wgt, p.x);

                // Implicit capture
                p.wgt *= 1.0 - P_abs;

                // Scatter
                if (rng.rand() > P_straight_ahead) p.turn();

                // Russian Roulette
                Roulette(p, rng);
                if (p.alive == false) {
#pragma omp atomic
                  cnts_sum += p.xs_evals;
                }
              }

              if (E_tot > Esmp) {
                // Improve Esmp
                Esmp = E_tot;
#pragma omp atomic write
                xs->E_bin_smp[p.bin] = Esmp;
              }

              // Split if needed
              if (p.alive and (std::abs(p.wgt) >= WGT_SPLIT)) {
                double n_new = std::round(std::abs(p.wgt));
                p.wgt /= n_new;
                for (int j = 0; j < static_cast<int>(n_new - 1); j++) {
                  Particle p_daughter(p.x, p.u, p.wgt);
                  p_daughter.bin = p.bin;
                  this_thread_Splits.push_back(p_daughter);
                }
              }  // split
            }
          }
        }  // While alive
      }    // For all particles
#pragma omp atomic write
      pcg_seeds[thread_id] = rng.get_seed();

#pragma omp single
      { Bank.clear(); }

#pragma omp barrier
#pragma omp critical
      {
        Bank.insert(std::end(Bank), std::begin(this_thread_Splits),
                    std::end(this_thread_Splits));
      }
    }  // Parallel

    n_particles = static_cast<int>(Bank.size());

  }  // While split particles

  n_xs_evals += cnts_sum;
}
