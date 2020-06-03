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

XS::XS() {
  for (int i = 0; i < N_XS_BINS; i++) {
    E_bin_max.push_back(0.0);
    E_bin_smp.push_back(0.0);
  }
  adjusted = false;
}

void XS::Adust_Smp_XSs(double P) {
  if (adjusted == false) {
    for (int i = 0; i < N_XS_BINS; i++) {
      E_bin_smp[i] *= P;
    }
    adjusted = true;
  }
}

// Constant XS Implementation
Constant::Constant(double xs) : XS() {
  this->XS = xs;
  this->Emax = xs;
  this->Esmp = xs;

  for (auto& bin : E_bin_max) bin = xs;

  for (auto& bin : E_bin_smp) bin = xs;

  std::cout << std::fixed << std::setprecision(2);
  std::cout << " Using constant XS Σ(x) = " << xs << " ,  ∀ x ∈ [";
  std::cout << X_MIN << "," << X_MAX << "]\n\n";
}

double Constant::E(double x) { return this->XS; }

double Constant::T(double x_0, double x) {
  double x_min = std::min(x_0, x);
  double x_max = std::max(x_0, x);

  // Can only integrate within problem bounds
  if (x_min < X_MIN) x_min = X_MIN;
  if (x_max > X_MAX) x_max = X_MAX;

  return this->XS * (x_max - x_min);  // Always positive
}

void Constant::Reset_XSs() {
  for (auto& bin : E_bin_max) bin = this->XS;
  for (auto& bin : E_bin_smp) bin = this->XS;
  adjusted = false;
}

// Linear XS Implementation
Linear::Linear(double E_x_min, double E_x_max) : XS() {
  assert(E_x_max >= 0.0);
  assert(E_x_min >= 0.0);

  this->E_x_min = E_x_min;
  this->E_x_max = E_x_max;

  this->m = (E_x_max - E_x_min) / (X_MAX - X_MIN);
  this->b = E_x_max - this->m * X_MAX;

  this->Emax = std::max(E_x_max, E_x_min);
  this->Esmp = this->Emax;

  // Get max and sample for all bins
  if (this->m > 0) {
    for (size_t b = 1; b <= N_XS_BINS; b++) {
      double x_bin_max = static_cast<double>(b) * dx_xs_bins + X_MIN;
      this->E_bin_max[b - 1] = this->E(x_bin_max);
      this->E_bin_smp[b - 1] = this->E(x_bin_max);
    }
  } else {
    for (size_t b = 0; b < N_XS_BINS; b++) {
      double x_bin_max = static_cast<double>(b) * dx_xs_bins + X_MIN;
      this->E_bin_max[b] = this->E(x_bin_max);
      this->E_bin_smp[b] = this->E(x_bin_max);
    }
  }

  std::cout << std::fixed << std::setprecision(2);
  std::cout << "\n Using linear XS Σ(x) = " << this->m << " * x + ";
  std::cout << this->b << " ,  ∀ x ∈ [" << X_MIN << "," << X_MAX << "]\n\n";
}

double Linear::E(double x) { return (this->m * x + this->b); }

double Linear::T(double x_0, double x) {
  double x_min = std::min(x_0, x);
  double x_max = std::max(x_0, x);

  // Can only integrate within problem bounds
  if (x_min < X_MIN) x_min = X_MIN;
  if (x_max > X_MAX) x_max = X_MAX;

  double integral = this->m * 0.5 * x_max * x_max + this->b * x_max;
  integral -= this->m * 0.5 * x_min * x_min + this->b * x_min;

  return integral;
}

void Linear::Reset_XSs() {
  // Get max and sample for all bins
  if (this->m > 0) {
    for (size_t b = 1; b <= N_XS_BINS; b++) {
      double x_bin_max = static_cast<double>(b) * dx_xs_bins + X_MIN;
      this->E_bin_max[b - 1] = this->E(x_bin_max);
      this->E_bin_smp[b - 1] = this->E(x_bin_max);
    }
  } else {
    for (size_t b = 0; b < N_XS_BINS; b++) {
      double x_bin_max = static_cast<double>(b) * dx_xs_bins + X_MIN;
      this->E_bin_max[b] = this->E(x_bin_max);
      this->E_bin_smp[b] = this->E(x_bin_max);
    }
  }

  adjusted = false;
}

// Exponential XS Implementation
Exponential::Exponential(double Coeff, double Exp_coeff) : XS() {
  assert(Coeff > 0.0);

  this->Coeff = Coeff;
  this->Exp_coeff = Exp_coeff;

  if (Exp_coeff > 0.0) {
    this->Emax = this->E(X_MAX);
    for (size_t b = 1; b <= N_XS_BINS; b++) {
      double x_bin_max = static_cast<double>(b) * dx_xs_bins + X_MIN;
      this->E_bin_max[b - 1] = this->E(x_bin_max);
      this->E_bin_smp[b - 1] = this->E(x_bin_max);
    }
  } else {
    this->Emax = this->E(X_MIN);
    for (size_t b = 0; b < N_XS_BINS; b++) {
      double x_bin_max = static_cast<double>(b) * dx_xs_bins + X_MIN;
      this->E_bin_max[b] = this->E(x_bin_max);
      this->E_bin_smp[b] = this->E(x_bin_max);
    }
  }
  this->Esmp = this->Emax;

  std::cout << std::fixed << std::setprecision(2);
  std::cout << "\n Using Exponential XS Σ(x) = " << this->Coeff << " * Exp[";
  std::cout << this->Exp_coeff << " * x] ,  ∀ x ∈ [";
  std::cout << X_MIN << "," << X_MAX << "]\n\n";
}

double Exponential::E(double x) {
  return this->Coeff * std::exp(this->Exp_coeff * x);
}

double Exponential::T(double x_0, double x) {
  double x_min = std::min(x_0, x);
  double x_max = std::max(x_0, x);

  // Can only integrate within problem bounds
  if (x_min < X_MIN) x_min = X_MIN;
  if (x_max > X_MAX) x_max = X_MAX;

  double integral = this->E(x_max) / this->Exp_coeff;
  integral -= this->E(x_min) / this->Exp_coeff;

  return integral;
}

void Exponential::Reset_XSs() {
  if (Exp_coeff > 0.0) {
    this->Emax = this->E(X_MAX);
    for (size_t b = 1; b <= N_XS_BINS; b++) {
      double x_bin_max = static_cast<double>(b) * dx_xs_bins + X_MIN;
      this->E_bin_max[b - 1] = this->E(x_bin_max);
      this->E_bin_smp[b - 1] = this->E(x_bin_max);
    }
  } else {
    this->Emax = this->E(X_MIN);
    for (size_t b = 0; b < N_XS_BINS; b++) {
      double x_bin_max = static_cast<double>(b) * dx_xs_bins + X_MIN;
      this->E_bin_max[b] = this->E(x_bin_max);
      this->E_bin_smp[b] = this->E(x_bin_max);
    }
  }
  this->Esmp = this->Emax;

  adjusted = false;
}

// Gaussian XS Implementation
Gaussian::Gaussian(double Coeff, double z, double sig, double b) : XS() {
  assert(Coeff > 0.0);
  assert(b >= 0.0);
  assert(z > X_MIN);
  assert(z < X_MAX);

  this->Coeff = Coeff;
  this->z = z;
  this->sig = sig;
  this->b = b;

  this->Emax = Coeff + b;
  this->Esmp = Coeff + b;

  // Get xs bin maxes
  for (size_t b = 0; b < N_XS_BINS; b++) {
    double x_low_bin = static_cast<double>(b) * dx_xs_bins + X_MIN;
    double x_hi_bin = static_cast<double>(b + 1) * dx_xs_bins + X_MIN;

    if (x_hi_bin < z) {
      this->E_bin_max[b] = this->E(x_hi_bin);
      this->E_bin_smp[b] = this->E(x_hi_bin);
    } else if (x_low_bin > z) {
      this->E_bin_max[b] = this->E(x_low_bin);
      this->E_bin_smp[b] = this->E(x_low_bin);
    } else if (x_low_bin < z and z < x_hi_bin) {
      this->E_bin_max[b] = this->E(z);
      this->E_bin_smp[b] = this->E(z);
    } else if (x_low_bin == z) {
      this->E_bin_max[b] = this->E(z);
      this->E_bin_smp[b] = this->E(z);
    } else if (x_hi_bin == z) {
      this->E_bin_max[b] = this->E(z);
      this->E_bin_smp[b] = this->E(z);
    }
  }

  std::cout << std::fixed << std::setprecision(2);
  std::cout << " Using Gaussian XS Σ(x) = " << this->Coeff << " * Exp[";
  std::cout << "-(x-" << this->z << ")^2/(" << this->sig << ")^2]";
  std::cout << " + " << this->b << " ,  ∀ x ∈ [" << X_MIN << ",";
  std::cout << X_MAX << "]\n\n";
}

double Gaussian::E(double x) {
  double exponent = -((x - this->z) / this->sig) * ((x - this->z) / this->sig);
  return this->Coeff * std::exp(exponent) + this->b;
}

double Gaussian::T(double x_0, double x) {
  double x_min = std::min(x_0, x);
  double x_max = std::max(x_0, x);

  // Can only integrate within problem bounds
  if (x_min < X_MIN) x_min = X_MIN;
  if (x_max > X_MAX) x_max = X_MAX;

  double integral =
      0.5 * std::sqrt(M_PI) * Coeff * sig * std::erf((x_max - z) / sig) +
      b * x_max;
  integral -=
      0.5 * std::sqrt(M_PI) * Coeff * sig * std::erf((x_min - z) / sig) +
      b * x_min;

  return integral;
}

void Gaussian::Reset_XSs() {
  // Get xs bin maxes
  for (size_t b = 0; b < N_XS_BINS; b++) {
    double x_low_bin = static_cast<double>(b) * dx_xs_bins + X_MIN;
    double x_hi_bin = static_cast<double>(b + 1) * dx_xs_bins + X_MIN;

    if (x_hi_bin < z) {
      this->E_bin_max[b] = this->E(x_hi_bin);
      this->E_bin_smp[b] = this->E(x_hi_bin);
    } else if (x_low_bin > z) {
      this->E_bin_max[b] = this->E(x_low_bin);
      this->E_bin_smp[b] = this->E(x_low_bin);
    } else if (x_low_bin < z and z < x_hi_bin) {
      this->E_bin_max[b] = this->E(z);
      this->E_bin_smp[b] = this->E(z);
    } else if (x_low_bin == z) {
      this->E_bin_max[b] = this->E(z);
      this->E_bin_smp[b] = this->E(z);
    } else if (x_hi_bin == z) {
      this->E_bin_max[b] = this->E(z);
      this->E_bin_smp[b] = this->E(z);
    }
  }

  adjusted = false;
}
