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

#ifndef MESHSMOOTHERCELL_H
#define MESHSMOOTHERCELL_H

#include "cellShape.H"
#include "polyMesh.H"

#include "SmootherBoundary.H"
#include "SmootherPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class MeshSmootherCell Declaration
\*---------------------------------------------------------------------------*/

class SmootherCell
{
    // Static Data

        // Pointer to smoother boundary
        static SmootherBoundary* _bnd;

        // Transformation parameter
        static scalar _transParam;


    // Private Data

        // Reference of cell shape
        const cellShape& _cellShape;

        //- Cell quality (mean ratio)
        scalar _quality;


    // Private Member Functions

        // Tetrahedral quality computation
        scalar tetCellQuality(const label ref) const;

        // get point from cell position
        inline const point& initPt(const label p) const;
        inline const point& relaxPt(const label p) const;


public:

    //- Constructors

        //- Construct from cellShape
        explicit SmootherCell(const cellShape& shape);


    //- Destructor
    ~SmootherCell() = default;


    // Member Functions

        // Set static values
        static void setStaticItems(SmootherBoundary* bnd, const scalar t);

        // Copy clone
        autoPtr<SmootherCell> clone() const
        {
            return autoPtr<SmootherCell>(new SmootherCell(*this));
        }

        // The cell quality
        constexpr scalar quality() const noexcept { return _quality; }

        void computeQuality();

        // Transform cell
        pointField geometricTransform();
};

const point &SmootherCell::initPt(const label p) const
{
    return _bnd->pt(_cellShape[p])->getInitialPoint();
}

const point &SmootherCell::relaxPt(const label p) const
{
    return _bnd->pt(_cellShape[p])->getRelaxedPoint();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // MESHSMOOTHERCELL_H

// ************************************************************************* //
