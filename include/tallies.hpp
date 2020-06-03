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
#ifndef TALLIES_H
#define TALLIES_H

#include <cmath>

#include "settings.hpp"

// Tallies for current batch
extern double batch_escaped;  // Number of particles which escaped
extern double
    batch_real_collisions;           // Number of collisions with real estimator
extern double batch_all_collisions;  // Number of collisions with all estimator

extern double batch_real_coll_dens[N_FLUX_BINS];
extern double batch_all_coll_dens[N_FLUX_BINS];

// Problem tallies
extern double avg_leakage_rate;
extern double avg_leakage_rate_variance;

extern double avg_coll_per_part_real;
extern double avg_coll_per_part_real_variance;

extern double avg_coll_per_part_all;
extern double avg_coll_per_part_all_variance;

extern double avg_coll_dens_real[N_FLUX_BINS]
                                [2];  // [0] is avg, [1] is variance
extern double avg_coll_dens_all[N_FLUX_BINS][2];

// Other tally variables
extern int n_particles_transported;
extern int n_sign_changes_in_weight;
extern double n_xs_evals;

// Scoring Functions
void Score_Escape(double w);
void Score_Collision_Real(double w, double x);
void Score_Collision_All(double w, double x);
void Accumulate_Batch(int batch);
void Zero_Batch_Scores();
void Zero_All_Scores();
#endif
