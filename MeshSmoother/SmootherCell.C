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

#include "SmootherCell.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::scalar Foam::SmootherCell::_transParam = 1.0;
Foam::SmootherBoundary* Foam::SmootherCell::_bnd = nullptr;


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

void Foam::SmootherCell::setStaticItems(SmootherBoundary* bnd, const scalar t)
{
    _transParam = t;
    _bnd = bnd;
}


// * * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * //

Foam::scalar Foam::SmootherCell::tetCellQuality
(
    const cellShape& shape,
    const label ref
)
{
    const label v1[] = {3, 0, 1, 2, 7, 4, 5, 6};
    const label v2[] = {4, 5, 6, 7, 5, 6, 7, 4};
    const label v3[] = {1, 2, 3, 0, 0, 1, 2, 3};

    const Tensor<scalar> mA
    (
        relaxPt(shape[v1[ref]]) - relaxPt(shape[ref]),
        relaxPt(shape[v2[ref]]) - relaxPt(shape[ref]),
        relaxPt(shape[v3[ref]]) - relaxPt(shape[ref])
    );
    const scalar sigma(det(mA));

    if (sigma < ROOTVSMALL)
    {
        return 0;
    }

    // 3 * pow(sigma, 2.0/3.0)/magSqr(mA)
    return scalar(3) * Foam::cbrt(Foam::sqr(sigma)) / magSqr(mA);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SmootherCell::SmootherCell(const cellShape& shape)
:
    _cellShape(shape),
    _quality(0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::SmootherCell::quality(const cellShape& shape)
{
    scalar avgQuality = 0;

    forAll(shape, pointi)
    {
        const scalar tetQuality(tetCellQuality(shape, pointi));

        if (tetQuality < ROOTVSMALL)
        {
            return 0;
        }

        avgQuality += tetQuality;
    }

    return (avgQuality / scalar(8));
}


Foam::pointField Foam::SmootherCell::geometricTransform(const cellShape& shape)
{
    // TODO if found a way to have the face label from cellShape,
    // store face centre in MeshSmoother and avoid use of pointFields dO
    // faceList fList = cell.faces();
    // fList[0].centre(pts);

    // Faces centres
    pointField fc(6);
    const label f[] = {0, 0, 1, 2, 0, 4};
    const label g[] = {1, 4, 5, 6, 3, 7};
    const label h[] = {2, 5, 6, 7, 7, 6};
    const label i[] = {3, 1, 2, 3, 4, 5};

    for (label j = 0; j < 6; ++j)
    {
        fc[j] =
        (
            initPt(shape[f[j]])
          + initPt(shape[g[j]])
          + initPt(shape[h[j]])
          + initPt(shape[i[j]])
        ) / 4.0;
    }

    // Transformed points
    pointField H(8);
    const label a[] = {0, 0, 0, 0, 5, 5, 5, 5};
    const label b[] = {1, 2, 3, 4, 4, 1, 2, 3};
    const label d[] = {4, 1, 2, 3, 1, 2, 3, 4};

    for (label j = 0; j < 8; ++j)
    {
        const point c = (fc[a[j]] + fc[b[j]] + fc[d[j]])/3.0;
        const point n = (fc[b[j]] - fc[a[j]]) ^ (fc[d[j]] - fc[a[j]]);
        H[j] = c + _transParam/std::sqrt(mag(n))*n;
    }

    // Length scale
    scalar length = 0.0, lengthN = 0.0;
    const label k[] = {0, 1, 2, 3, 0, 1, 2, 3, 4, 5, 6, 7};
    const label l[] = {1, 2, 3, 0, 4, 5, 6, 7, 5, 6, 7, 4};

    for (label j = 0; j < 12; ++j)
    {
        length += mag
        (
            initPt(shape[k[j]])
          - initPt(shape[l[j]])
        );
        lengthN += mag(H[k[j]] - H[l[j]]);
    }
    length /= lengthN;

    const point c = (H[0] + H[1] + H[2] + H[3] + H[4] + H[5] + H[6] + H[7])/8.0;
    const pointField C(8, c);

    return pointField(C + length*(H - C));
}


// ************************************************************************* //
