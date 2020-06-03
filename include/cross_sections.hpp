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
#ifndef XS_H
#define XS_H

#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>

#include "settings.hpp"

// -----------------------------------------------------------------------------
// XS Class
//    Virtual base class for all XS types
// -----------------------------------------------------------------------------
class XS {
 public:
  XS();
  virtual ~XS() = default;

  virtual double E(double x) = 0;

  virtual double T(double x_0, double x) = 0;

  virtual void Adust_Smp_XSs(double P);
  virtual void Reset_XSs() = 0;

  // Data
  double Emax;                    // max over entire spacial domain
  double Esmp;                    // Sampling xs for entire domain
  std::vector<double> E_bin_max;  // Max in each bin
  std::vector<double> E_bin_smp;  // Sampling XS in each bin
  bool adjusted;
};

// -----------------------------------------------------------------------------
// Constant Class
//    Constant xs extends XS (E(x) = XS)
// -----------------------------------------------------------------------------
class Constant : public XS {
 public:
  Constant(double xs);
  ~Constant() = default;

  double E(double x);
  double T(double x_0, double x);

  void Reset_XSs();

  double XS;
};

// -----------------------------------------------------------------------------
// Linear Class
//    E(x) = mx + b  and must be positive over entire domain !
// -----------------------------------------------------------------------------
class Linear : public XS {
 public:
  Linear(double E_x_min, double E_x_max);
  ~Linear() = default;

  double E(double x);
  double T(double x_0, double x);

  void Reset_XSs();

  double E_x_min;
  double E_x_max;
  double m;
  double b;
};

// -----------------------------------------------------------------------------
// Exponential Class
//    E(x) = Coeff * exp(Exp_coeff * x)  and Coeff must be positive
// -----------------------------------------------------------------------------
class Exponential : public XS {
 public:
  Exponential(double Coeff, double Exp_coeff);
  ~Exponential() = default;

  double E(double x);
  double T(double x_0, double x);

  void Reset_XSs();

  double Coeff;
  double Exp_coeff;
};

// -----------------------------------------------------------------------------
// Gaussian Class
//    E(x) = Coeff * exp(-(x - z)^2/sig^2) + b  and Coeff must be positive
// -----------------------------------------------------------------------------
class Gaussian : public XS {
 public:
  Gaussian(double Coeff, double z, double sig, double b);
  ~Gaussian() = default;

  double E(double x);
  double T(double x_0, double x);

  void Reset_XSs();

  double Coeff;
  double z;
  double sig;
  double b;
};
#endif
