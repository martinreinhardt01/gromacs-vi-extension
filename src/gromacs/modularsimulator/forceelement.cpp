/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
/*! \libinternal
 * \brief Defines the force element for the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 */

#include "gmxpre.h"

#include "forceelement.h"

#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/force.h"
#include "gromacs/mdlib/force_flags.h"
#include "gromacs/mdlib/mdatoms.h"

#include "energyelement.h"
#include "freeenergyperturbationelement.h"
#include "statepropagatordata.h"

struct gmx_edsam;
struct gmx_enfrot;
struct gmx_multisim_t;
class history_t;
struct t_graph;

namespace gmx
{
ForceElement::ForceElement(
        StatePropagatorData           *statePropagatorData,
        EnergyElement                 *energyElement,
        FreeEnergyPerturbationElement *freeEnergyPerturbationElement,
        bool                           isDynamicBox,
        FILE                          *fplog,
        const t_commrec               *cr,
        const t_inputrec              *inputrec,
        const MDAtoms                 *mdAtoms,
        t_nrnb                        *nrnb,
        t_forcerec                    *fr,
        t_fcdata                      *fcd,
        gmx_wallcycle                 *wcycle,
        MdrunScheduleWorkload         *runScheduleWork,
        gmx_vsite_t                   *vsite,
        ImdSession                    *imdSession,
        pull_t                        *pull_work) :
    nextNSStep_(-1),
    nextEnergyCalculationStep_(-1),
    nextVirialCalculationStep_(-1),
    nextFreeEnergyCalculationStep_(-1),
    statePropagatorData_(statePropagatorData),
    energyElement_(energyElement),
    freeEnergyPerturbationElement_(freeEnergyPerturbationElement),
    localTopology_(nullptr),
    isDynamicBox_(isDynamicBox),
    ddBalanceRegionHandler_(cr),
    lambda_(),
    fplog_(fplog),
    cr_(cr),
    inputrec_(inputrec),
    mdAtoms_(mdAtoms),
    nrnb_(nrnb),
    wcycle_(wcycle),
    fr_(fr),
    vsite_(vsite),
    imdSession_(imdSession),
    pull_work_(pull_work),
    fcd_(fcd),
    runScheduleWork_(runScheduleWork)
{
    lambda_.fill(0);
}

void ForceElement::scheduleTask(
        Step step, Time time,
        const RegisterRunFunctionPtr &registerRunFunction)
{
    unsigned int flags = (
            GMX_FORCE_STATECHANGED |
            GMX_FORCE_ALLFORCES |
            (isDynamicBox_ ? GMX_FORCE_DYNAMICBOX : 0) |
            (nextVirialCalculationStep_ == step ? GMX_FORCE_VIRIAL : 0) |
            (nextEnergyCalculationStep_ == step ? GMX_FORCE_ENERGY : 0) |
            (nextFreeEnergyCalculationStep_ == step ? GMX_FORCE_DHDL : 0) |
            (nextNSStep_ == step ? GMX_FORCE_NS : 0));

    (*registerRunFunction)(
            std::make_unique<SimulatorRunFunction>(
                    [this, step, time, flags]()
                    {run(step, time, flags); }));
}

void ForceElement::elementSetup()
{
    GMX_ASSERT(localTopology_, "Setup called before local topology was set.");

}

void ForceElement::run(Step step, Time time, unsigned int flags)
{
    // Disabled functionality
    Awh            *awh              = nullptr;
    gmx_edsam      *ed               = nullptr;
    gmx_multisim_t *ms               = nullptr;
    gmx_enfrot     *enforcedRotation = nullptr;
    t_graph        *graph            = nullptr;

    /* The coordinates (x) are shifted (to get whole molecules)
     * in do_force.
     * This is parallelized as well, and does communication too.
     * Check comments in sim_util.c
     */
    auto           x       = statePropagatorData_->positionsView();
    auto           forces  = statePropagatorData_->forcesView();
    auto           forcesA = statePropagatorData_->forcesAView();
    auto           forcesB = statePropagatorData_->forcesBView();
    auto           box     = statePropagatorData_->constBox();
    history_t     *hist    = nullptr; // disabled

    tensor         force_vir = {{0}};
    // TODO: Make lambda const (needs some adjustments in lower force routines)
    ArrayRef<real> lambda = freeEnergyPerturbationElement_ ?
        freeEnergyPerturbationElement_->lambdaView() : lambda_;

    do_force(fplog_, cr_, ms, inputrec_, awh, enforcedRotation, imdSession_,
             pull_work_,
             step, nrnb_, wcycle_, localTopology_,
             box, x, hist,
             forces, forcesA, forcesB, force_vir, mdAtoms_->mdatoms(), energyElement_->enerdata(), fcd_,
             lambda_, graph,
             fr_, runScheduleWork_, vsite_, energyElement_->muTot(), time, ed,
             static_cast<int>(flags), ddBalanceRegionHandler_);
    energyElement_->addToForceVirial(force_vir, step);
}

void ForceElement::setTopology(const gmx_localtop_t *top)
{
    localTopology_ = top;
}

SignallerCallbackPtr ForceElement::registerNSCallback()
{
    return std::make_unique<SignallerCallback>(
            [this](Step step, Time gmx_unused time)
            {this->nextNSStep_ = step; });
}

SignallerCallbackPtr ForceElement::
    registerEnergyCallback(EnergySignallerEvent event)
{
    if (event == EnergySignallerEvent::energyCalculationStep)
    {
        return std::make_unique<SignallerCallback>(
                [this](Step step, Time)
                {nextEnergyCalculationStep_ = step; });
    }
    if (event == EnergySignallerEvent::virialCalculationStep)
    {
        return std::make_unique<SignallerCallback>(
                [this](Step step, Time){nextVirialCalculationStep_ = step; });
    }
    if (event == EnergySignallerEvent::freeEnergyCalculationStep)
    {
        return std::make_unique<SignallerCallback>(
                [this](Step step, Time){nextFreeEnergyCalculationStep_ = step; });
    }
    return nullptr;
}
}  // namespace std
