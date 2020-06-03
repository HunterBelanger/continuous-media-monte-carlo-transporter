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

double Direct_Sampler::get_collision_distance(Particle& p, double P_nc,
                                              const std::unique_ptr<XS>& xs,
                                              PCG& rng) {
  double x_0 = p.x;
  double x_max;
  double G = 1.0 - P_nc;
  if (p.u == -1)
    x_max = X_MIN;
  else
    x_max = X_MAX;

  double xi = rng.rand();
  double tau_hat = -std::log(1.0 - G * xi);

  // Initial guess is half way between current location and surface

  double s_1 = 0.5 * (x_0 + x_max);
  double s_0 = x_max;
  double d = std::abs(s_1 - s_0);
  double d_old = d + 1.0;
  while (d > 1e-6) {
    s_0 = s_1;
    double g = tau_hat - xs->T(x_0, s_0);
    p.xs_eval();
    double dg = -xs->E(s_0);
    p.xs_eval();

    s_1 = s_0 - (g / dg);

    d_old = d;
    d = std::abs(s_1 - s_0);

    // Call Bisection Method if Needed
    if (p.u == -1.0 and (s_1 > p.x or s_1 < X_MIN)) {
      s_1 = bisection(p, tau_hat, xs);
      s_0 = s_1;
    } else if (p.u == 1.0 and (s_1 > X_MAX or s_1 < p.x)) {
      s_1 = bisection(p, tau_hat, xs);
      s_0 = s_1;
    } else if (d_old - d >= 0.0) {  // Needed due to circular cases occuring
      s_1 = bisection(p, tau_hat, xs);
      s_0 = s_1;
    }
  }

  double x_coll;
  if (p.u == -1)
    x_coll = std::max(s_1, s_0);
  else
    x_coll = std::min(s_1, s_0);

  return std::abs(x_0 - x_coll);
};

double Direct_Sampler::bisection(Particle& p, double tau_hat,
                                 const std::unique_ptr<XS>& xs) {
  // Determine limits
  double x_low, x_hi;
  if (p.u == -1) {
    x_hi = p.x;
    x_low = X_MIN;
  } else {
    x_low = p.x;
    x_hi = X_MAX;
  }

  // Do Bisection Method
  while (x_hi - x_low > 1e-6) {
    double x = 0.5 * (x_hi + x_low);
    double f = tau_hat - xs->T(p.x, x);

    p.xs_eval();

    if (p.u == -1) {
      if (f < 0.0)
        x_low = x;
      else if (f > 0.0)
        x_hi = x;
    } else {
      if (f < 0.0)
        x_hi = x;
      else if (f > 0.0)
        x_low = x;
    }
  }

  return 0.5 * (x_hi + x_low);
}

double Direct_Sampler::get_P_no_collision(Particle& p,
                                          const std::unique_ptr<XS>& xs) {
  double x = p.x;
  double x_surf;
  if (p.u == -1)
    x_surf = X_MIN;
  else
    x_surf = X_MAX;

  double tau = xs->T(x, x_surf);
  p.xs_eval();
  return std::exp(-tau);
}

void Direct_Sampler::Transport(const std::vector<Particle>& bank,
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
    for (int n = 0; n < static_cast<double>(bank.size()); n++) {
#pragma omp atomic
      n_particles_transported++;

      Particle p = bank[n];
      while (p.alive) {
        double P_leak = get_P_no_collision(p, xs);
        if (rng.rand() < P_leak) {  // Particle leaks
          Score_Escape(p.wgt);
          p.kill();
#pragma omp atomic
          cnts_sum += p.xs_evals;
        } else {  // Collision
          // Get collision distance
          double d = get_collision_distance(p, P_leak, xs, rng);
          p.move(d);

          // Score real collision
          Score_Collision_Real(p.wgt, p.x);
          Score_Collision_All(p.wgt, p.x);

          // Implicit capture
          p.wgt *= 1.0 - P_abs;

          // Scatter
          if (rng.rand() < P_straight_ahead) p.turn();

          // Russin Roulette
          Roulette(p, rng);
          if (p.alive == false) {
#pragma omp atomic
            cnts_sum += p.xs_evals;
          }
        }
      }  // While alive
    }    // For all particles in bank
#pragma omp critical
    { pcg_seeds[thread_id] = rng.get_seed(); }
  }  // Paralllel
  n_xs_evals += cnts_sum;
}
