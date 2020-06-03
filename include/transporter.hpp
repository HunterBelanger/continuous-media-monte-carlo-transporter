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
#ifndef TRANSPORTERS_H
#define TRANSPORTERS_H

#include <omp.h>
#include <memory>
#include <vector>

#include "cross_sections.hpp"
#include "particle.hpp"
#include "settings.hpp"
#include "utils.hpp"

class Transporter {
 public:
  virtual ~Transporter(){};
  virtual void Transport(const std::vector<Particle>& bank,
                         const std::unique_ptr<XS>& xs) = 0;
};

class Delta_Tracker : public Transporter {
 public:
  Delta_Tracker() { std::cout << " Using Delta Tracking\n\n"; };
  ~Delta_Tracker() = default;

  void Transport(const std::vector<Particle>& bank,
                 const std::unique_ptr<XS>& xs);
};

class Regional_Delta_Tracker : public Transporter {
 public:
  Regional_Delta_Tracker() {
    std::cout << " Using Regional Delta Tracking\n\n";
  };
  ~Regional_Delta_Tracker() = default;

  void Transport(const std::vector<Particle>& bank,
                 const std::unique_ptr<XS>& xs);
};

class Negative_Weighted_Delta_Tracker : public Transporter {
 public:
  Negative_Weighted_Delta_Tracker(double P_ = 0.8, double q_ = 0.3)
      : P{P_}, q{q_} {
    std::cout << " Using Negative Weighted Delta Tracking\n\n";
  };
  ~Negative_Weighted_Delta_Tracker() = default;

  void Transport(const std::vector<Particle>& bank,
                 const std::unique_ptr<XS>& xs);

 private:
  double P;
  double q;
};

class Regional_Negative_Weighted_Delta_Tracker : public Transporter {
 public:
  Regional_Negative_Weighted_Delta_Tracker(double P_ = 0.8, double q_ = 0.3)
      : P{P_}, q{q_} {
    std::cout << " Using Regional Netative Weighted Delta Tracking\n\n";
  };
  ~Regional_Negative_Weighted_Delta_Tracker() = default;

  void Transport(const std::vector<Particle>& bank,
                 const std::unique_ptr<XS>& xs);

 private:
  double P;
  double q;
};

class Carter_Transporter : public Transporter {
 public:
  Carter_Transporter(double P_ = 0.85) : P{P_} {
    std::cout << " Using Carter Transport\n\n";
  };
  ~Carter_Transporter() = default;

  void Transport(const std::vector<Particle>& bank,
                 const std::unique_ptr<XS>& xs);

 private:
  double P;
};

class Regional_Carter_Transporter : public Transporter {
 public:
  Regional_Carter_Transporter(double P_ = 0.85) : P{P_} {
    std::cout << " Using Regional Carter Transport\n\n";
  };
  ~Regional_Carter_Transporter() = default;

  void Transport(const std::vector<Particle>& bank,
                 const std::unique_ptr<XS>& xs);

 private:
  double P;
};

class Improving_Regional_Carter_Transporter : public Transporter {
 public:
  Improving_Regional_Carter_Transporter(double P_ = 0.85) : P{P_} {
    std::cout << " Using Improving Regional Carter Transport\n\n";
  };
  ~Improving_Regional_Carter_Transporter() = default;

  void Transport(const std::vector<Particle>& bank,
                 const std::unique_ptr<XS>& xs);

 private:
  double P;
};

class Substeper : public Transporter {
 public:
  Substeper(int n_steps) : Nsteps{n_steps} {
    std::cout << " Using Sub-Stepping Transport\n\n";
  };

  ~Substeper() = default;

  void Transport(const std::vector<Particle>& bank,
                 const std::unique_ptr<XS>& xs);

 private:
  int Nsteps;             // Number of steps or bins
  std::vector<double> E;  // Value of XS used in each step bin

  void fill_xs_table(const std::unique_ptr<XS>& xs);
  int get_xs_bin(double x);
  double dist_to_bin(const Particle& p);
};

class Direct_Sampler : public Transporter {
 public:
  Direct_Sampler() { std::cout << " Using Direct-Sampling\n\n"; };
  ~Direct_Sampler() = default;

  void Transport(const std::vector<Particle>& bank,
                 const std::unique_ptr<XS>& xs);

 private:
  double get_collision_distance(Particle& p, double P_nc,
                                const std::unique_ptr<XS>& xs, PCG& rng);
  double get_P_no_collision(Particle& p, const std::unique_ptr<XS>& xs);
  double bisection(Particle& p, double tau_hat, const std::unique_ptr<XS>& xs);
};
#endif
