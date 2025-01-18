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
#include "gmxpre.h"

#include "enerdata_utils.h"
#include <iostream>
#include <fstream>

#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

gmx_enerdata_t::gmx_enerdata_t(int numEnergyGroups,
                               int numFepLambdas) :
    grpp(numEnergyGroups),
    enerpart_lambda(numFepLambdas == 0 ? 0 : numFepLambdas + 1),
    enerpart_lambdaA(numFepLambdas == 0 ? 0 : numFepLambdas + 1),
    enerpart_lambdaB(numFepLambdas == 0 ? 0 : numFepLambdas + 1),
    foreign_grpp(numEnergyGroups)
{
}

static real sum_v(int n, gmx::ArrayRef<const real> v)
{
    real t;
    int  i;

    t = 0.0;
    for (i = 0; (i < n); i++)
    {
        t = t + v[i];
    }

    return t;
}

void sum_epot(gmx_grppairener_t *grpp, real *epot)
{
    int i;

    /* Accumulate energies */
    epot[F_COUL_SR]  = sum_v(grpp->nener, grpp->ener[egCOULSR]);
    epot[F_LJ]       = sum_v(grpp->nener, grpp->ener[egLJSR]);
    epot[F_LJ14]     = sum_v(grpp->nener, grpp->ener[egLJ14]);
    epot[F_COUL14]   = sum_v(grpp->nener, grpp->ener[egCOUL14]);
    epot[F_EPOTA]    = sum_v(grpp->nener, grpp->enerA);
    epot[F_EPOTB]    = sum_v(grpp->nener, grpp->enerB);


/* lattice part of LR doesnt belong to any group
 * and has been added earlier
 */
    epot[F_BHAM]     = sum_v(grpp->nener, grpp->ener[egBHAMSR]);

    epot[F_EPOT] = 0;
    for (i = 0; (i < F_EPOT); i++)
    {
        if (i != F_DISRESVIOL && i != F_ORIRESDEV)
        {
            epot[F_EPOT] += epot[i];
        }
    }
}

void sum_dhdl(gmx_enerdata_t *enerd, gmx::ArrayRef<const real> lambda, const t_lambda &fepvals)
{
    int    index;
    int    this_state = fepvals.init_fep_state;
    double E_this = 0.;

    enerd->dvdl_lin[efptVDW] += enerd->term[F_DVDL_VDW];  /* include dispersion correction */
    for (int i = 0; i < efptNR; i++)
    {
        if (fepvals.separate_dvdl[i])
        {
            /* could this be done more readably/compactly? */
            switch (i)
            {
                case (efptMASS):
                    index = F_DKDL;
                    break;
                case (efptCOUL):
                    index = F_DVDL_COUL;
                    break;
                case (efptVDW):
                    index = F_DVDL_VDW;
                    break;
                case (efptBONDED):
                    index = F_DVDL_BONDED;
                    break;
                case (efptRESTRAINT):
                    index = F_DVDL_RESTRAINT;
                    break;
                default:
                    index = F_DVDL;
                    break;
            }
            enerd->term[index] = enerd->dvdl_lin[i] + enerd->dvdl_nonlin[i];
            if (debug)
            {
                fprintf(debug, "dvdl-%s[%2d]: %f: non-linear %f + linear %f\n",
                        efpt_names[i], i, enerd->term[index], enerd->dvdl_nonlin[i], enerd->dvdl_lin[i]);
            }
        }
        else
        {
            enerd->term[F_DVDL] += enerd->dvdl_lin[i] + enerd->dvdl_nonlin[i];

            if (debug)
            {
                fprintf(debug, "dvd-%sl[%2d]: %f: non-linear %f + linear %f\n",
                        efpt_names[0], i, enerd->term[F_DVDL], enerd->dvdl_nonlin[i], enerd->dvdl_lin[i]);
            }
        }
    }

    // calculate dhdl variable for TI
    if (fepvals.var_morph){

    	real Va, Vb;
    	real dvAdl = 0.;
    	real dvBdl = 0.;

        if (fepvals.sc_alpha == 0){
        	Va = enerd->term[F_EPOTA];
        	Vb = enerd->term[F_EPOTB] + enerd->term[F_VM_LIN_CORR] + enerd->term[F_VM_COUL_CORR];
        }
        else{

         	real *Va_get = enerd->grpp.enerA.data();
    		real *Vb_get = enerd->grpp.enerB.data();

    		dvAdl  = enerd->dvdl_AB[0];
    		dvBdl  = enerd->dvdl_AB[1];

        	// VM_LIN_CORR: long range dispersion correction term for energy
        	// in case the end states are lambda dependent

    		Va = Va_get[0];
    		Vb = Vb_get[0] + enerd->term[F_VM_LIN_CORR] + enerd->term[F_VM_COUL_CORR];
        }

        /* linear test: the Va, Vb, dvAdla and dvBdl variables were tracked separately,
           but are combined here in the (regular) linear fashion.
           Test that the same results (and, therefore, not contribution remains
           untracked, as when var_morph is disabled
        */
    	if (fepvals.linear_test) {
    		real   lam          = lambda[0];
    		// if efptVDW == efptCOUL, index is the same
	        enerd->term[F_DVDL] = Vb - Va + (1.-lam) * dvAdl + lam * dvBdl;
    	}
    	else{
    		real   lam          = lambda[0];
    		const real kbT      = 2.479; // kJ/mol
    		real   s            = fepvals.smoothing;
    		real   deltag_est   = fepvals.deltag_est;
    		double expA         = exp( -s * Va / kbT );
    		double expB         = exp( -s * (Vb - deltag_est) / kbT );
    		double denom        = s * ( (1.-lam) * expA + lam * expB );
    		double gA           = ( s * (1.-lam) * dvAdl / kbT + 1. ) * expA;
    		double gB           = ( s * lam      * dvBdl / kbT - 1. ) * expB;
    		enerd->term[F_DVDL] = kbT * ( gA + gB ) / denom;
    	}
    }


    if (fepvals.separate_dvdl[efptBONDED])
    {
        enerd->term[F_DVDL_BONDED] += enerd->term[F_DVDL_CONSTR];
    }
    else
    {
        enerd->term[F_DVDL] += enerd->term[F_DVDL_CONSTR];
    }

    for (int i = 0; i < fepvals.n_lambda; i++)
    {

    	/* Parts treated by var_morph have been tracked separately
    	 * in A and B and are not included in  enerpart_lambda
    	 * Therefore, collecting these terms is conducted no matter
    	 * if var_morph is used or not.
    	 * In the next step, if enabled, all var_morph contributions
    	 * are added up.
    	 */

        double &enerpart_lambda  = enerd->enerpart_lambda[i + 1];
        double &enerpart_lambdaA = enerd->enerpart_lambdaA[i + 1];
        double &enerpart_lambdaB = enerd->enerpart_lambdaB[i + 1];

        /* for linear interpolation:
         * note we are iterating over fepvals here!
           For the current lam, dlam = 0 automatically,
           so we don't need to add anything to the
           enerd->enerpart_lambda[0] */

        /* we don't need to worry about dvdl_lin contributions to dE at
           current lambda, because the contributions to the current
           lambda are automatically zeroed */

        for (gmx::index j = 0; j < lambda.ssize(); j++)
        {
            /* Note that this loop is over all dhdl components, not just the separated ones */
            const double dlam  = fepvals.all_lambda[j][i] - lambda[j];

            if (! fepvals.var_morph)
            {
                enerpart_lambda   += dlam*enerd->dvdl_lin[j];
            }

            /* Constraints can not be evaluated at foreign lambdas, so we add
             * a linear extrapolation. This is an approximation, but usually
             * quite accurate since constraints change little between lambdas.
             */
            if ((j == efptBONDED && fepvals.separate_dvdl[efptBONDED]) ||
                (j == efptFEP && !fepvals.separate_dvdl[efptBONDED]))
            {
                enerpart_lambda += dlam*enerd->term[F_DVDL_CONSTR];
            }

            if (j == efptMASS && !fepvals.separate_dvdl[j])
            {
                enerpart_lambda += dlam*enerd->term[F_DKDL];
            }

            if (debug)
            {
                fprintf(debug, "enerdiff lam %g: (%15s), non-linear %f linear %f*%f\n",
                        fepvals.all_lambda[j][i], efpt_names[j],
                        enerpart_lambda - enerd->enerpart_lambda[0],
                        dlam, enerd->dvdl_lin[j]);
            }
        }

        if (fepvals.var_morph){

        	real Va, Vb;

            if (fepvals.sc_alpha == 0){
            	Va = enerd->term[F_EPOTA];
            	Vb = enerd->term[F_EPOTB] + enerd->term[F_VM_LIN_CORR] + enerd->term[F_VM_COUL_CORR];
            }
            else{
            	Va = enerpart_lambdaA;
            	Vb = enerpart_lambdaB + enerd->term[F_VM_LIN_CORR] + enerd->term[F_VM_COUL_CORR] ;
            }

        	if (fepvals.linear_test) {
        		// if efptVDW == efptCOUL, index is the same
        		real lam = lambda_var_morph(fepvals.all_lambda[efptCOUL][i], fepvals.all_lambda[efptVDW][i]);
                enerpart_lambda    += (1. - lam) * Va + lam * Vb;
        	}
        	else{
        		real   lam         = lambda_var_morph(fepvals.all_lambda[efptCOUL][i], fepvals.all_lambda[efptVDW][i]);
        		const real kbT     = 2.479; // kJ/mol
        		real   s           = fepvals.smoothing;
        		real   deltag_est  = fepvals.deltag_est;
        		double gA          = (1.-lam) * exp( -s * Va / kbT );
        		double gB          = lam      * exp( -s * (Vb - deltag_est) / kbT );
        		enerpart_lambda    += -kbT / s * log(gA + gB);
        	}
        }

        if (i == this_state )
        	E_this = enerpart_lambda;
    }
    if (fepvals.var_morph)
    {
		for (int i = 0; i < fepvals.n_lambda; i++)
		{
			double &enerpart_lambda  = enerd->enerpart_lambda[i + 1];
			enerpart_lambda -= E_this;
		}
    }


    /* The constrain contribution is now included in other terms, so clear it */
    enerd->term[F_DVDL_CONSTR] = 0;
}


real get_energy_var_morph(int state, real lam, gmx_enerdata_t *enerd, t_lambda *fepvals){

    const double    kbT = 2.479; // kJ/mol

    real s           = fepvals->smoothing;
	real Vout        = 0.;
	real deltag_est  = fepvals->deltag_est;
    real Va          = 0.;
    real Vb          = 0.;

    if (fepvals->sc_alpha == 0){
    	Va = enerd->term[F_EPOTA];
    	Vb = enerd->term[F_EPOTB]+ enerd->term[F_VM_LIN_CORR] + enerd->term[F_VM_COUL_CORR];
    }
    else{
    	Va = enerd->enerpart_lambdaA[state+1];
    	Vb = enerd->enerpart_lambdaB[state+1] + enerd->term[F_VM_LIN_CORR] + enerd->term[F_VM_COUL_CORR];
    }

    Vb -= deltag_est;

	if ( lam == 0 ){
		Vout = Va;
	}
	else if ( lam == 1 ){
		Vout = Vb;
	}
	else if (fepvals->var_morph && !fepvals->linear_test)
	{

        double expA = exp( -s * Va  / kbT );
	    double expB = exp( -s * Vb  / kbT );

	    Vout = -kbT / s * log( (1.-lam) * expA + lam * expB );
	}
	else if (fepvals->var_morph && fepvals->linear_test == 1)
	{
		Vout = (1. -lam) * Va + lam * Vb;
	}

	return Vout ;

}

void write_energies_to_file(FILE *fplog, gmx_enerdata_t *enerd, gmx::ArrayRef<const real> lambda, t_lambda *fepvals, gmx_bool bDoDHDL){

	real lam         = lambda_var_morph(lambda);
	int  state       = fepvals->init_fep_state;
    real delta_lam   = delta_lambda_var_morph(fepvals);

	if ( bDoDHDL ){

		if (lam == 0){
			fprintf(fplog, "eXx  %f  %f \n", get_energy_var_morph(state, lam, enerd, fepvals), get_energy_var_morph(state+1, lam+delta_lam, enerd, fepvals));
		}
		else if (lam == 1){
			fprintf(fplog, "eXx  %f  %f \n", get_energy_var_morph(state, lam, enerd, fepvals), get_energy_var_morph(state-1, lam-delta_lam, enerd, fepvals));
		}
		else{
			fprintf(fplog, "eXx  %f  %f  %f \n", get_energy_var_morph(state, lam, enerd, fepvals), get_energy_var_morph(state-1, lam-delta_lam, enerd, fepvals), get_energy_var_morph(state+1, lam+delta_lam, enerd, fepvals));
		}

		//energyfile.close();
	}
}


real lambda_var_morph(real *lambda){
	real lam_out = 0.;

	if ( lambda[efptCOUL] == lambda[efptVDW] )
		lam_out = lambda[efptCOUL];
	else if( lambda[efptCOUL] == 0 && lambda[efptVDW] != 0)
		lam_out = lambda[efptVDW];
	else if( lambda[efptCOUL] != 0 && lambda[efptVDW] == 0)
		lam_out = lambda[efptCOUL];

	// in case of non-eq. case lambda is only in lambda[0]
	lam_out = std::max(lam_out, lambda[0]);
    return lam_out;
}

real lambda_var_morph(real lam_coul, real lam_vdw){
	if (lam_coul == lam_vdw)
		return lam_coul;
	else if (lam_coul == 0 && lam_vdw != 0 )
		return lam_vdw;
	else if (lam_coul != 0 && lam_vdw == 0 )
        return lam_coul;
    else
    {
    	std::cout << "Problem: cannot determine a var_morph lambda " << std::endl;
	    return 0.;
    }
}

real lambda_var_morph(gmx::ArrayRef<const real> lambda){
	real lam_out = 0.;

	if ( lambda[efptCOUL] == lambda[efptVDW] )
		lam_out = lambda[efptCOUL];
	else if( lambda[efptCOUL] == 0 && lambda[efptVDW] != 0)
		lam_out = lambda[efptVDW];
	else if( lambda[efptCOUL] != 0 && lambda[efptVDW] == 0)
		lam_out = lambda[efptCOUL];
    return lam_out;
}

real delta_lambda_var_morph(t_lambda *fepvals){
	real delta_lam = 0.;

	if ( fepvals->all_lambda[efptCOUL][1] == fepvals->all_lambda[efptVDW][1] )
		delta_lam = fepvals->all_lambda[efptCOUL][1] - fepvals->all_lambda[efptCOUL][0];
	else if( fepvals->all_lambda[efptCOUL][1] == 0 &&  fepvals->all_lambda[efptVDW][1] != 0)
		delta_lam = fepvals->all_lambda[efptVDW][1] -  fepvals->all_lambda[efptVDW][0];
	else if( fepvals->all_lambda[efptCOUL][1] != 0 &&  fepvals->all_lambda[efptVDW][1] == 0)
		delta_lam = fepvals->all_lambda[efptCOUL][1] - fepvals->all_lambda[efptCOUL][0];
    return delta_lam;
}





void write_dhdl_to_log(FILE *fplog, gmx_enerdata_t *enerd, const t_lambda &fepvals)
{
	int index;
    for (int i = 0; i < efptNR; i++)
    {
        if (fepvals.separate_dvdl[i])
        {
            /* could this be done more readably/compactly? */
            switch (i)
            {
                case (efptMASS):
                    index = F_DKDL;
                    break;
                case (efptCOUL):
                    index = F_DVDL_COUL;
                    break;
                case (efptVDW):
                    index = F_DVDL_VDW;
                    break;
                case (efptBONDED):
                    index = F_DVDL_BONDED;
                    break;
                case (efptRESTRAINT):
                    index = F_DVDL_RESTRAINT;
                    break;
                default:
                    index = F_DVDL;
                    break;
            }

            fprintf(fplog, "dvdl-%s[%2d]: %f: non-linear %f + linear %f\n",
                    efpt_names[i], i, enerd->term[index], enerd->dvdl_nonlin[i], enerd->dvdl_lin[i]);

        }
        else
        {
           fprintf(fplog, "dvd-%sl[%2d]: %f: non-linear %f + linear %f\n",
                     efpt_names[0], i, enerd->term[F_DVDL], enerd->dvdl_nonlin[i], enerd->dvdl_lin[i]);

        }
    }
    fprintf(fplog, "enerd->term[F_DVDL_VDW]:  %f  \n", enerd->term[F_DVDL_VDW]);

}

void reset_foreign_enerdata(gmx_enerdata_t *enerd)
{
    int  i, j;

    /* First reset all foreign energy components.  Foreign energies always called on
       neighbor search steps */
    for (i = 0; (i < egNR); i++)
    {
        for (j = 0; (j < enerd->grpp.nener); j++)
        {
            enerd->foreign_grpp.ener[i][j] = 0.0;
        }
    }

    for (j = 0; (j < enerd->grpp.nener); j++)
    {
        enerd->foreign_grpp.enerA[j] = 0.0;
        enerd->foreign_grpp.enerB[j] = 0.0;
    }

    /* potential energy components */
    for (i = 0; (i <= F_EPOT); i++)
    {
        enerd->foreign_term[i] = 0.0;
    }
}

void reset_enerdata(gmx_enerdata_t *enerd)
{
    int      i, j;

    /* First reset all energy components. */
    for (i = 0; (i < egNR); i++)
    {
        for (j = 0; (j < enerd->grpp.nener); j++)
        {
            enerd->grpp.ener[i][j] = 0.0;
        }
    }
    for (i = 0; i < efptNR; i++)
    {
        enerd->dvdl_lin[i]    = 0.0;
        enerd->dvdl_nonlin[i] = 0.0;
        enerd->dvdl_AB[0]     = 0.0;
        enerd->dvdl_AB[1]     = 0.0;
    }

    for (j = 0; (j < enerd->grpp.nener); j++)
    {
        enerd->grpp.enerA[j] = 0.0;
        enerd->grpp.enerB[j] = 0.0;
    }

    /* Normal potential energy components */
    for (i = 0; (i <= F_EPOT); i++)
    {
        enerd->term[i] = 0.0;
    }
    enerd->term[F_DVDL]            = 0.0;
    enerd->term[F_DVDL_COUL]       = 0.0;
    enerd->term[F_DVDL_VDW]        = 0.0;
    enerd->term[F_DVDL_BONDED]     = 0.0;
    enerd->term[F_DVDL_RESTRAINT]  = 0.0;
    enerd->term[F_DKDL]            = 0.0;
    enerd->term[F_VM_COUL_CORR]    = 0.0;
    std::fill(enerd->enerpart_lambda.begin(), enerd->enerpart_lambda.end(), 0);
    std::fill(enerd->enerpart_lambdaA.begin(), enerd->enerpart_lambdaA.end(), 0);
    std::fill(enerd->enerpart_lambdaB.begin(), enerd->enerpart_lambdaB.end(), 0);

    /* reset foreign energy data - separate function since we also call it elsewhere */
    reset_foreign_enerdata(enerd);
}
