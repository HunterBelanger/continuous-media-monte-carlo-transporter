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

void Substeper::fill_xs_table(const std::unique_ptr<XS>& xs) {
  E.clear();

  double dx = LENGTH / static_cast<double>(Nsteps);
  for (int i = 0; i < Nsteps; i++) {
    double x = static_cast<double>(i) * dx + 0.5 * dx;
    E.push_back(xs->E(x));
  }
}

int Substeper::get_xs_bin(double x) {
  double dx = LENGTH / static_cast<double>(Nsteps);
  int bin = std::floor(x / dx);
  return bin;
}

double Substeper::dist_to_bin(const Particle& p) {
  double dx = LENGTH / static_cast<double>(Nsteps);
  int bin = std::floor(p.x / dx);

  if (p.u == -1) {
    double x_bin = static_cast<double>(bin) * dx + X_MIN;
    return p.x - x_bin;
  } else {
    double x_bin = static_cast<double>(bin + 1) * dx + X_MIN;
    return x_bin - p.x;
  }
}

void Substeper::Transport(const std::vector<Particle>& bank,
                          const std::unique_ptr<XS>& xs) {
  fill_xs_table(xs);
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

      bool collision;
      Particle p = bank[n];
      double Esmp = E[get_xs_bin(p.x)];
      p.xs_eval();
      while (p.alive) {
        collision = false;
        while (collision == false) {
          double d = -std::log(rng.rand()) / (Esmp);
          double d_xs_bin = dist_to_bin(p);
          if (d_xs_bin < d) {
            p.move(d_xs_bin + 1e-6);
            if ((p.x >= X_MAX) or (p.x <= X_MIN)) {
              Score_Escape(p.wgt);
              collision = true;
              p.kill();
#pragma omp atomic
              cnts_sum += p.xs_evals;
            } else {
              Esmp = E[get_xs_bin(p.x)];
              p.xs_eval();
            }
          } else {
            p.move(d);
            collision = true;
          }
        }

        if (p.alive) {
          // Score real collision
          Score_Collision_Real(p.wgt, p.x);
          Score_Collision_All(p.wgt, p.x);

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
