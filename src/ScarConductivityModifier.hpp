/*

Copyright (c) 2005-2015, University of Oxford.
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

#ifndef SCARCONDUCTIVITYMODIFIER_HPP_
#define SCARCONDUCTIVITYMODIFIER_HPP_

#include "AbstractConductivityModifier.hpp"
#include "ScarCellFactory.hpp"

/**
 * This is a special conductivity modifier class, which is the easiest way to get
 * complicated variation in conductivity throughout the tissue.
 *
 * This class handles conductivity for the scar (circular or square) and its boundary
 * together with a 'cut' introduced half way along the domain.
 */
template<unsigned DIM>
class ScarConductivityModifier : public AbstractConductivityModifier<DIM,DIM>
{
private:
    /** Pointer to mesh */
    AbstractTetrahedralMesh<DIM,DIM>* mpMesh;
    /** Special shape enumeration defined in the ScarCellFactory.hpp */
    ScarShape mShape;
    /** Working memory for returning a zero tensor */
    c_matrix<double,DIM,DIM> mZeroTensor;
    /** Working memory for returning modified tensors */
    c_matrix<double,DIM,DIM> mReturnTensor;

    /** Width of the domain we're simulating. */
    double mRegionWidth;
    /** Radius of the scar in the centre. */
    double mScarRadius;
    /** Scaling factor to multiply the intracellular conductivity tensor entries by */
    double mScarScaling;
    /** Scaling factor for 'chi' to divide the conductivity tensor entries by (NB both intra and extra need dividing here.*/
    double mScarChiScaling;
    /** Width of the boundary region (myocytes but different conductivity) */
    double mBoundaryWidth;
    /** Scaling factor for conductivity tensor entries */
    double mBoundaryScaling;
    /** Whether a lesion area should be applied */
    bool mLesionApplied;
    /** Whether a cut should also be applied*/
    bool mCutApplied;
    /** The width of the cut */
    double mCutWidth;

public:
    /**
     * Constructor
     *
     * @param pMesh  A pointer to the mesh we're going to use.
     * @param shape  CIRCLE or SQUARE enumeration (see ScarCellFactory)
     * @param regionWidth  Width of the entire domain.
     * @param scarRadius  Radius (or now 1/2 width for square) of the scar
     * @param scarScaling  Scaling factor for condictivity in the scar
     * @param chiScaling  Scaling factor for chi - surface area to volume tissue ratio.
     * @param boundaryWidth  Extra radius/width of border around the scar
     * @param boundaryScaling  Scaling factor for condictivity in the boundary
     * @param scarActive  Whether to introduce a scar
     * @param cutActive  Whether to introduce a cut
     */
    ScarConductivityModifier(AbstractTetrahedralMesh<DIM,DIM>* pMesh,
                             ScarShape shape,
                             double regionWidth,
                             double scarRadius,
                             double scarScaling,
                             double chiScaling,
                             double boundaryWidth,
                             double boundaryScaling,
                             bool scarActive,
                             bool cutActive)
        : AbstractConductivityModifier<DIM,DIM>(),
          mpMesh(pMesh),
          mShape(shape),
          mRegionWidth(regionWidth),
          mScarRadius(scarRadius),
          mScarScaling(scarScaling),
          mScarChiScaling(chiScaling),
          mBoundaryWidth(boundaryWidth),
          mBoundaryScaling(boundaryScaling),
          mLesionApplied(scarActive),
          mCutApplied(cutActive),
          mCutWidth(0.01)
    {
        assert(pMesh != NULL);
        mZeroTensor = zero_matrix<double>(DIM,DIM);
    }

    /**
     * This is the pure virtual method that this class provides to the CardiacTissue to modify conductivities.
     *
     * @param elementIndex  Global element index
     * @param rOriginalConductivity  The original tensor for this element and domain.
     * @param domainIndex  0 for intracellular, 1 for extracellular.
     *
     * @return Modified conductivity tensor for the scar/border/cut.
     */
    c_matrix<double,DIM,DIM>& rCalculateModifiedConductivityTensor(unsigned elementIndex,
                                                               const c_matrix<double,DIM,DIM>& rOriginalConductivity,
                                                               unsigned domainIndex)
    {
        c_vector<double, DIM> loc = mpMesh->GetElement(elementIndex)->CalculateCentroid();
        const double x = loc[0];
        const double y = loc[1];

        bool intracellular = (domainIndex==0u);
        bool inside_boundary = false;
        bool inside_scar = false;

        // Just work out whether we are within the scar or boundary regions
        if (mLesionApplied)
        {
            if (mShape==CIRCLE)
            {
                double distance_from_origin_squared = (x-mRegionWidth/2.0)*(x-mRegionWidth/2.0)+(y-mRegionWidth/2.0)*(y-mRegionWidth/2.0);
                inside_boundary = (distance_from_origin_squared < (mScarRadius+mBoundaryWidth)*(mScarRadius+mBoundaryWidth));
                inside_scar = (distance_from_origin_squared < mScarRadius*mScarRadius);
            }
            else
            {
                assert(mShape==SQUARE);
                double x_distance_from_origin_squared = (x-mRegionWidth/2.0)*(x-mRegionWidth/2.0);
                double y_distance_from_origin_squared = (x-mRegionWidth/2.0)*(x-mRegionWidth/2.0);
                inside_boundary = ((x_distance_from_origin_squared < (mScarRadius+mBoundaryWidth)*(mScarRadius+mBoundaryWidth)) &&
                                   (y_distance_from_origin_squared < (mScarRadius+mBoundaryWidth)*(mScarRadius+mBoundaryWidth))   );
                inside_scar = ((x_distance_from_origin_squared < mScarRadius*mScarRadius) &&
                               (y_distance_from_origin_squared < mScarRadius*mScarRadius)   );
            }
        }

        if (inside_scar && intracellular)
        {
            mReturnTensor = rOriginalConductivity;
            if (intracellular)
            {
                mReturnTensor *= mScarScaling;
            }
            mReturnTensor /= mScarChiScaling;
            return mReturnTensor;
        }
        else if ( DIM==3u && mCutApplied && ((x-mRegionWidth/2.0)*(x-mRegionWidth/2.0) < mCutWidth*mCutWidth) )
        {
            // (If DIM==2 then the mesh takes care of the cut).
            // In the cut region (fresh air, not extracellular conductivity either)
            return mZeroTensor;
        }
        else if (inside_boundary && intracellular)
        {
            mReturnTensor = rOriginalConductivity*mBoundaryScaling;
            return mReturnTensor;
        }
        else
        {
            // everywhere else.
            mReturnTensor = rOriginalConductivity;
            return mReturnTensor;
        }
    }
};


#endif // SCARCONDUCTIVTYMODIFIER_HPP_
