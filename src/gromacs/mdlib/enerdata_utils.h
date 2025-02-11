/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017,2018,2019, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */

#ifndef GMX_MDLIB_ENERDATA_UTILS_H
#define GMX_MDLIB_ENERDATA_UTILS_H

#include "gromacs/math/arrayrefwithpadding.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/arrayref.h"

struct gmx_enerdata_t;
struct gmx_grppairener_t;
struct t_lambda;

void reset_foreign_enerdata(gmx_enerdata_t *enerd);
/* Resets only the foreign energy data */

void reset_enerdata(gmx_enerdata_t *enerd);
/* Resets the energy data */

void sum_epot(gmx_grppairener_t *grpp, real *epot);
/* Locally sum the non-bonded potential energy terms */

void sum_dhdl(gmx_enerdata_t *enerd, gmx::ArrayRef<const real> lambda, const t_lambda &fepvals);
/* Sum the free energy contributions */

void write_dhdl_to_log(FILE *fplog, gmx_enerdata_t *enerd, const t_lambda &fepvals);

void write_energies_to_file(FILE *fplog, gmx_enerdata_t *enerd, gmx::ArrayRef<const real> lambda, t_lambda *fepvals, gmx_bool bDoDHDL);

void write_dhdl_to_log(FILE *fplog, gmx_enerdata_t *enerd);

real lambda_var_morph(real *lambda);

real lambda_var_morph(gmx::ArrayRef<const real> lambda);

real lambda_var_morph(real lam_coul, real lam_vdw);

real delta_lambda_var_morph(t_lambda *fepvals);

real get_energy_var_morph(int state, real lam, gmx_enerdata_t *enerd, t_lambda *fepvals);

void reduce_AB_energies(gmx_enerdata_t *enerd, real *lambda, t_lambda *fepvals);


#endif
