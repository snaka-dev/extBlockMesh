/*---------------------------------------------------------------------------*\
  extBlockMesh
  Copyright (C) 2014 Etudes-NG
  Copyright (C) 2020 OpenCFD Ltd.
  ---------------------------------
License
    This file is part of extBlockMesh.

    extBlockMesh is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    extBlockMesh is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with extBlockMesh.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#ifndef MESHSMOOTHER_H
#define MESHSMOOTHER_H

#include "polyMesh.H"
#include "scalarField.H"
#include "PtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class blockMesh;
class SmootherControl;
class SmootherParameter;
class SmootherBoundary;

/*---------------------------------------------------------------------------*\
                        Class MeshSmoother Declaration
\*---------------------------------------------------------------------------*/

class MeshSmoother
{
    // Private Data

        // Pointer of parents
        polyMesh *_polyMesh;

        // Current cell quality. See SmootherCell
        scalarField _cellQuality;

        // MeshSmoother child items
        SmootherControl* _ctrl;
        SmootherParameter* _param;
        SmootherBoundary* _bnd;


    // Private Member Functions

        // Quality analysis
        void analyseMeshQuality();
        void analyseMeshQuality(const labelHashSet& cell);
        void qualityStats();

        // GETMe smoothing
        labelHashSet addTransformedElementNodeWeight();
        void addUnTransformedElementNodeWeight(labelHashSet &tp);
        bool untransformedAndhavePointTransformed
        (
            const label cellI,
            const labelHashSet& tp
        );

        void iterativeNodeRelaxation(labelHashSet &tP, const scalarList &r);

        // Run one iteration
        bool runIteration();
        pointField getMovedPoints() const;

        // Smoothing algo
        void GETMeSmoothing();
        void snapSmoothing();

public:

    //- Constructors

        //- Construct from polyMesh and smoother dictionary
        MeshSmoother(polyMesh *mesh, const dictionary& smootherDict);


    //- Destructor
    ~MeshSmoother();


    // Member Functions

        //- Smooth the mesh acording to smoothDict
        void update();

        //- Smooth mesh acording to smoothDict,
        //- writing intermediate mesh, optionally write the quality field
        //  \return the number of time-steps written
        label updateAndWrite(Time& runTime, const bool withQuality=false);

        //- Get tranformation treshold
        scalar getTransformationTreshold() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // MESHSMOOTHER_H

// ************************************************************************* //
