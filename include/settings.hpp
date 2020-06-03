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
#ifndef SETTINGS_H
#define SETTINGS_H
// Problem Domain
constexpr double X_MIN = 0.0;
constexpr double X_MAX = 2.0;
constexpr double LENGTH = X_MAX - X_MIN;

// Problem Divisions
constexpr int N_XS_BINS = 5;
constexpr int N_FLUX_BINS = 100;
constexpr double dx_xs_bins = LENGTH / static_cast<double>(N_XS_BINS);
constexpr double dx_flux_bins = LENGTH / static_cast<double>(N_FLUX_BINS);

// Number of Particle and batches
extern int NPARTICLES;  // Particles per batch
extern int NBATCHES;    // Number of batches
extern int NTRIALS;     // Number of independent trials to run

// Weight Settings
constexpr double WGT_CUTOFF = 0.6;
constexpr double WGT_SURVIVE = 1.0;
constexpr double WGT_SPLIT = 2.0;

// System Properties
constexpr double P_abs = 0.3;
constexpr double P_sct = 0.7;
constexpr double P_straight_ahead = 0.5;

#endif
