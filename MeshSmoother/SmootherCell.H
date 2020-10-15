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
        static scalar tetCellQuality
        (
            const cellShape& shape,
            const label ref
        );


        inline static const point& initPt(const label pointi)
        {
            return _bnd->pt(pointi)->getInitialPoint();
        }

        inline static const point& relaxPt(const label pointi)
        {
            return _bnd->pt(pointi)->getRelaxedPoint();
        }


public:

    //- Constructors

        //- Construct from cellShape
        explicit SmootherCell(const cellShape& shape);


    //- Destructor
    ~SmootherCell() = default;


    // Static Functions

        //- Set static values
        static void setStaticItems(SmootherBoundary* bnd, const scalar t);

        //- The cell quality for (hex) cell shape
        static scalar quality(const cellShape& shape);

        //- Transform (hex) cell
        static pointField geometricTransform(const cellShape& shape);


    // Member Functions

        //- Copy clone
        autoPtr<SmootherCell> clone() const
        {
            return autoPtr<SmootherCell>(new SmootherCell(*this));
        }

        //- The current cell quality
        constexpr scalar quality() const noexcept { return _quality; }

        //- Calculate the cell quality
        inline void computeQuality()
        {
            _quality = SmootherCell::quality(_cellShape);
        }

        //- Transform (hex) cell
        inline pointField geometricTransform() const
        {
            return SmootherCell::geometricTransform(_cellShape);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // MESHSMOOTHERCELL_H

// ************************************************************************* //
