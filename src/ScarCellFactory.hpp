/*

Copyright (c) 2005-2019, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef SCARCELLFACTORY_HPP_
#define SCARCELLFACTORY_HPP_

#ifdef CHASTE_CVODE

#include <boost/shared_ptr.hpp>
#include <boost/pointer_cast.hpp>
#include <utility> // for std::pair

#include "AbstractCardiacCellFactory.hpp"
#include "FakeBathCell.hpp"
#include "RegularStimulus.hpp"
#include "SteadyStateRunner.hpp"

#include "li_mouse_2010Cvode.hpp"
#include "sachse_moreno_abildskov_2008_bCvode.hpp"

/**
 * Structure encapsulating the enumeration of different scar shapes
 */
typedef enum ScarShape_
{
    CIRCLE = 0,
    SQUARE
} ScarShape;

/**
 * This class creates the right kind of cells for each node in the computational mesh, and gives them back to the
 * tissue simulation as it is being set up.
 */
template<unsigned DIM>
class ScarCellFactory : public AbstractCardiacCellFactory<DIM>
{
private:
    boost::shared_ptr<RegularStimulus> mpLeftStimulus;
    boost::shared_ptr<RegularStimulus> mpLesionStimulus;
    double mSuctionElectrodeRadius;
    double mScarRadius;
    bool mUseFakeBathCell;
    bool mUseFibroblastModel;
    ScarShape mScarShape;
    bool mLesionPacing;
    bool mCreateCut;
    double mRegionWidth;
    double mRegionHeight;
    double mRegionDepth;

    /** Initial conditions for the cell model */
    std::vector<double> mCellModelICs;

    /** Helper method */
    void StoreNearestNodeInfo(Node<DIM>* pNode,
                              const double& rX,
                              const double& rY,
                              const double& rZ,
                              const double& rX_target,
                              const double& rY_target,
                              const double& rZ_target,
                              std::pair<unsigned, double>& rNearestNodeInfo)
    {
        double temp = (rX-rX_target)*(rX-rX_target) + (rY-rY_target)*(rY-rY_target) + (rZ-rZ_target)*(rZ-rZ_target);
        if ( temp < rNearestNodeInfo.second )
        {
            rNearestNodeInfo.first = pNode->GetIndex();
            rNearestNodeInfo.second = temp;
        }
    }

    /**
     * Calculate the nearest node to each of the following locations and also the distance,
     * store it for retrieval later on.
     */
    std::pair<unsigned, double> mCentreOfScar;      // Nearest to (0.25, 0.25, 0.025)
    std::pair<unsigned, double> mTopOfScar;         // Nearest to (0.25, 0.25, 0.05) 3D
    std::pair<unsigned, double> mLeftOfScar;        // Nearest to (0.18, 0.25, 0.025)
    std::pair<unsigned, double> mRightOfScar;       // Nearest to (0.32, 0.25, 0.025)
    std::pair<unsigned, double> mSideOfScar;        // Nearest to (0.25, 0.18, 0.025)
    std::pair<unsigned, double> mCentralTissue;     // Nearest to (0.25, 0.1, 0.025)
    std::pair<unsigned, double> mTopOfCentralTissue;// Nearest to (0.25, 0.1, 0.05) 3D
    std::pair<unsigned, double> mRightTissue;       // Nearest to (0.4, 0.25, 0.025)
    std::pair<unsigned, double> mTopOfRightTissue;  // Nearest to (0.4, 0.25, 0.05) 3D
    std::pair<unsigned, double> mLeftTissue;        // Nearest to (0.1, 0.25, 0.025)
    std::pair<unsigned, double> mTopOfLeftTissue;   // Nearest to (0.1, 0.25, 0.05) 3D
    std::pair<unsigned, double> mBaseTissue;        // Nearest to (0.25, 0, 0.025)

public:
    ScarCellFactory(const double& rRegionWidth,
                    const double& rScarRadius,
                    const double& rPeriod,
                    const double& rScarChiScaling,
                    bool neutralCellModel,
                    bool useFibroblastModel,
                    ScarShape shape,
                    bool lesionPacing,
                    bool createCut)
    : AbstractCardiacCellFactory<DIM>(),
      mpLeftStimulus(new RegularStimulus(-40000.0, 2, rPeriod, 10)),
      // We do a special scaling on the next line since I_stim must be scaled to have a larger/smaller effect with membrane density
      mpLesionStimulus(new RegularStimulus(-10000.0*rScarChiScaling, 4, rPeriod, 10)),
      //mpLesionStimulus(new RegularStimulus(-200000.0, 4, rPeriod, 10)),
      mSuctionElectrodeRadius(0.015),
      mScarRadius(rScarRadius),
      mUseFakeBathCell(neutralCellModel),
      mUseFibroblastModel(useFibroblastModel),
      mScarShape(shape),
      mLesionPacing(lesionPacing),
      mCreateCut(createCut),
      mCentreOfScar(UNSIGNED_UNSET, DBL_MAX), // These are all to track sensible node locations.
      mTopOfScar(UNSIGNED_UNSET, DBL_MAX),
      mLeftOfScar(UNSIGNED_UNSET, DBL_MAX),
      mRightOfScar(UNSIGNED_UNSET, DBL_MAX),
      mSideOfScar(UNSIGNED_UNSET, DBL_MAX),
      mCentralTissue(UNSIGNED_UNSET, DBL_MAX),
      mTopOfCentralTissue(UNSIGNED_UNSET, DBL_MAX),
      mRightTissue(UNSIGNED_UNSET, DBL_MAX),
      mTopOfRightTissue(UNSIGNED_UNSET, DBL_MAX),
      mLeftTissue(UNSIGNED_UNSET, DBL_MAX),
      mTopOfLeftTissue(UNSIGNED_UNSET, DBL_MAX),
      mBaseTissue(UNSIGNED_UNSET, DBL_MAX)
    {
        // Get the cell model to a sensible steady state before we start.
        // (NB We need a stimulus that is much smaller magnitude to stop it blowing up the cells)
        boost::shared_ptr<AbstractIvpOdeSolver> p_empty_solver;
        boost::shared_ptr<AbstractStimulusFunction> p_single_cell_stim(new RegularStimulus(-50.0, 2, rPeriod, 10));

        boost::shared_ptr<AbstractCvodeCell> p_cell(new Cellli_mouse_2010FromCellMLCvode(p_empty_solver, p_single_cell_stim));
        p_cell->SetParameter("membrane_non_inactivating_steady_state_potassium_current_conductance", 0.0);

        if (mLesionPacing)
        {
            p_cell->SetParameter("membrane_fast_sodium_current_conductance", 0.0);
        }

        std::cout << "Running AP model to initial steady state... ";
        SteadyStateRunner steady_runner(p_cell);
        steady_runner.RunToSteadyState();
        CopyToStdVector(p_cell->rGetStateVariables(),mCellModelICs);
        std::cout << "done!" << std::endl << std::flush;

        if (mCreateCut)
        {
            mRegionWidth = 0.7;  // Hardcoded from the mesh geometry in 2D_with_cuts_7_by_9_mm.msh.
            mRegionHeight = 0.9; // Hardcoded from the mesh geometry in 2D_with_cuts_7_by_9_mm.msh.
        }
        else
        {
            mRegionWidth = rRegionWidth;
            mRegionHeight = rRegionWidth; // Square so same width as height.
        }
        mRegionDepth = 0.05; // Hardcoded from the mesh geometry, only really used in 3D.
    }

    AbstractCardiacCellInterface* CreateCardiacCellForTissueNode(Node<DIM>* pNode)
    {
        AbstractCardiacCellInterface* p_cell;
        boost::shared_ptr<AbstractIvpOdeSolver> p_empty_solver;

        const double x = pNode->rGetLocation()[0];
        double y = 0.25; // Calculations should still work if we pretend this is the location in 1D (central line of 2D).
        if (DIM>=2u)
        {
            y = pNode->rGetLocation()[1];
        }

        double z = 0.0; // Calculations should still work if this is constant in 2D
        if (DIM==3)
        {
            z = pNode->rGetLocation()[2];
        }

        // We can work these out in any dimension, even if they are only used in certain ones.
        StoreNearestNodeInfo(pNode,x,y,z, 0.5*mRegionWidth, 0.5*mRegionHeight,   0.5*mRegionDepth, mCentreOfScar);
        StoreNearestNodeInfo(pNode,x,y,z, 0.5*mRegionWidth, 0.5*mRegionHeight,       mRegionDepth, mTopOfScar);
        StoreNearestNodeInfo(pNode,x,y,z, 0.36*mRegionWidth, 0.5*mRegionHeight,   0.5*mRegionDepth, mLeftOfScar);
        StoreNearestNodeInfo(pNode,x,y,z, 0.64*mRegionWidth, 0.5*mRegionHeight,   0.5*mRegionDepth, mRightOfScar);
        StoreNearestNodeInfo(pNode,x,y,z, 0.5*mRegionWidth, 0.36*mRegionHeight,   0.5*mRegionDepth, mSideOfScar);
        StoreNearestNodeInfo(pNode,x,y,z, 0.5*mRegionWidth, 0.8*mRegionHeight,   0.5*mRegionDepth, mCentralTissue);
        StoreNearestNodeInfo(pNode,x,y,z, 0.5*mRegionWidth, 0.8*mRegionHeight,       mRegionDepth, mTopOfCentralTissue);
        StoreNearestNodeInfo(pNode,x,y,z, 0.8*mRegionWidth, 0.5*mRegionHeight,   0.5*mRegionDepth, mRightTissue);
        StoreNearestNodeInfo(pNode,x,y,z, 0.8*mRegionWidth, 0.5*mRegionHeight,       mRegionDepth, mTopOfRightTissue);
        StoreNearestNodeInfo(pNode,x,y,z, 0.2*mRegionWidth, 0.5*mRegionHeight,   0.5*mRegionDepth, mLeftTissue);
        StoreNearestNodeInfo(pNode,x,y,z, 0.2*mRegionWidth, 0.5*mRegionHeight,       mRegionDepth, mTopOfLeftTissue);
        StoreNearestNodeInfo(pNode,x,y,z, 0.5*mRegionWidth,                0,   0.5*mRegionDepth, mBaseTissue);

        if ( ( mScarShape==CIRCLE
                  && (x-mRegionWidth/2.0)*(x-mRegionWidth/2.0) +(y-mRegionHeight/2.0)*(y-mRegionHeight/2.0)
                     < mScarRadius*mScarRadius) // Central region radius 0.1cm.
             || (mScarShape==SQUARE
                  && (x-mRegionWidth/2.0)*(x-mRegionWidth/2.0) < mScarRadius*mScarRadius
                  && (y-mRegionHeight/2.0)*(y-mRegionHeight/2.0) < mScarRadius*mScarRadius)
                )
        {
            // We are in the lesion region
            if (mUseFakeBathCell)
            {
                // First thing to try is a scar region that is still conductive,
                // but has a completely neutral cell model
                // (i.e. sum of ionic currents across membrane is always zero)
                // and voltage change is due to diffusion of ions from neighbours.
                p_cell = new FakeBathCell(p_empty_solver, this->mpZeroStimulus);

                // Give it a sensible starting voltage (otherwise it acts as a stimulus!).
                p_cell->SetVoltage(-80);
            }
            else if (mUseFibroblastModel)
            {
                p_cell = new Cellsachse_moreno_abildskov_2008_bFromCellMLCvode(p_empty_solver, this->mpZeroStimulus);
                // Give it a sensible starting voltage too.
                p_cell->SetVoltage(-80);

                const double scaling_factor = 10.0;

                const double leak_conductance = p_cell->GetParameter("membrane_leakage_current_conductance");
                p_cell->SetParameter("membrane_leakage_current_conductance",leak_conductance*scaling_factor);

                const double iK1_conductance = p_cell->GetParameter("membrane_inward_rectifier_potassium_current_conductance");
                p_cell->SetParameter("membrane_inward_rectifier_potassium_current_conductance",iK1_conductance*scaling_factor);

                const double shaker_conductance = p_cell->GetParameter("membrane_delayed_rectifier_potassium_current_conductance");
                p_cell->SetParameter("membrane_delayed_rectifier_potassium_current_conductance",shaker_conductance*scaling_factor);
            }
            else
            {
                p_cell = new Cellli_mouse_2010FromCellMLCvode(p_empty_solver, this->mpZeroStimulus);
            }

            if (mLesionPacing && (x-mRegionWidth/2.0)*(x-mRegionWidth/2.0) +(y-mRegionHeight/2.0)*(y-mRegionHeight/2.0)
                    < mSuctionElectrodeRadius*mSuctionElectrodeRadius)
            {
                // Applies to all cell models near the centre (will also be ones in the middle if we run this in 3D).
                p_cell->SetStimulusFunction(mpLesionStimulus);
            }
        }
        else // We are not in the lesion region
        {
            if (!mLesionPacing && (x<0.05+1e-6)  ) // Left hand edge
            {
                p_cell = new Cellli_mouse_2010FromCellMLCvode(p_empty_solver, mpLeftStimulus);
            }
            else
            {
                p_cell = new Cellli_mouse_2010FromCellMLCvode(p_empty_solver, this->mpZeroStimulus);
            }
        }

        // If this is a mouse myocyte cell model, then we want to make the APD a bit longer to match experiment
        // in control conditions. So we are going to block this current.
        //
        // Unfortunately AbstractCvodeCell doesn't have a HasParameter method itself and we need to tell
        // the compiler that it is really an AbstractUntemplatedParameterisedSystem to use this.
        if (boost::dynamic_pointer_cast<AbstractUntemplatedParameterisedSystem>(p_cell)
                ->HasParameter("membrane_non_inactivating_steady_state_potassium_current_conductance"))
        {
            p_cell->SetParameter("membrane_non_inactivating_steady_state_potassium_current_conductance", 0.0);
            // Use sensible steady state ICs.
            p_cell->SetStateVariables(mCellModelICs);
        }

        // If we are doing a current injection experiment there is Flecainide present to block fast sodium
        // and prevent activation waves - i.e. we mostly see diffusion of voltage and can see space constant.
        if (mLesionPacing && boost::dynamic_pointer_cast<AbstractUntemplatedParameterisedSystem>(p_cell)
                ->HasParameter("membrane_fast_sodium_current_conductance"))
        {
            p_cell->SetParameter("membrane_fast_sodium_current_conductance", 0.0);
            // Use sensible steady state ICs.
            p_cell->SetStateVariables(mCellModelICs);
        }

        return p_cell;
    }

    unsigned GetNodeIndexCentreOfScar()
    {
        return DoAnAllReduceOnThis(mCentreOfScar);
    }

    unsigned GetNodeIndexTopOfScar()
    {
        return DoAnAllReduceOnThis(mTopOfScar);
    }

    unsigned GetNodeIndexLeftOfScar()
    {
        return DoAnAllReduceOnThis(mLeftOfScar);
    }

    unsigned GetNodeIndexRightOfScar()
    {
        return DoAnAllReduceOnThis(mRightOfScar);
    }

    unsigned GetNodeIndexSideOfScar()
    {
        return DoAnAllReduceOnThis(mSideOfScar);
    }

    unsigned GetNodeIndexCentralTissue()
    {
        return DoAnAllReduceOnThis(mCentralTissue);
    }

    unsigned GetNodeIndexTopOfCentralTissue()
    {
        return DoAnAllReduceOnThis(mTopOfCentralTissue);
    }

    unsigned GetNodeIndexRightTissue()
    {
        return DoAnAllReduceOnThis(mRightTissue);
    }

    unsigned GetNodeIndexTopOfRightTissue()
    {
        return DoAnAllReduceOnThis(mTopOfRightTissue);
    }

    unsigned GetNodeIndexLeftTissue()
    {
        return DoAnAllReduceOnThis(mLeftTissue);
    }

    unsigned GetNodeIndexTopOfLeftTissue()
    {
        return DoAnAllReduceOnThis(mTopOfLeftTissue);
    }

    unsigned GetNodeIndexBaseOfTissue()
    {
        return DoAnAllReduceOnThis(mBaseTissue);
    }

    /**
     * Get the global node (across all processes) that has the smallest distance from a point of interest.
     *
     * @param rNodeIndexAndDistance A pair including a global node index and its distance.
     *
     * @return The node index nearest this point, across all processes.
     */
    unsigned DoAnAllReduceOnThis(std::pair<unsigned, double>& rNodeIndexAndDistance)
    {
        // This is a handy data structure that will work with MPI_DOUBLE_INT data type.
        struct
        {
            double distance;
            int node_index;
        } value, minval;

        value.node_index = rNodeIndexAndDistance.first;
        value.distance = rNodeIndexAndDistance.second;

        //std::cout << "Min distance = " << value.distance << " on nodes I own, global index = " << value.node_index << "\n";

        /* global minloc */
        MPI_Allreduce( &value, &minval, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD );

        //std::cout << "Min distance = " << minval.distance << " on all nodes, global index = " << minval.node_index << "\n";

        return minval.node_index;
    }


};

#endif // CHASTE_CVODE
#endif // SCARCELLFACTORY_HPP_
