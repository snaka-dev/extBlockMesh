/*---------------------------------------------------------------------------*\
  extBlockMesh
  Copyright (C) 2014 Etudes-NG
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

#include "argList.H"
#include "Time.H"
#include "IOdictionary.H"
#include "columnFvMesh.H"

#include "MeshSmoother.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::addBoolOption
    (
        "writeStep",
        "write mesh at different smoothing step"
    );

    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // Smoothing
    {
        #include "createSmoother.H"
        MeshSmoother& smoother = smootherPtr();

        if (args.optionFound("writeStep"))
        {
            smoother.updateAndWrite(runTime);
        }
        else
        {
            smoother.update();

            // Reset mesh directory to constant/polyMesh !!! Hard to find :P
            mesh.setInstance(runTime.constant());
        }
    }


    // Set the precision of the points data to 10
    IOstream::defaultPrecision(max(10u, IOstream::defaultPrecision()));

    Info<< nl << "Writing polyMesh" << endl;

    if (!args.optionFound("writeStep"))
    {
        mesh.removeFiles();
        if (!mesh.write())
        {
            FatalErrorIn(args.executable())
                << "Failed writing polyMesh."
                << exit(FatalError);
        }
    }

    //
    // write some information
    //
    {
        const polyPatchList& patches = mesh.boundaryMesh();

        Info<< "----------------" << nl
            << "Mesh Information" << nl
            << "----------------" << nl
            << "  " << "boundingBox: " << boundBox(mesh.points()) << nl
            << "  " << "nPoints: " << mesh.nPoints() << nl
            << "  " << "nCells: " << mesh.nCells() << nl
            << "  " << "nFaces: " << mesh.nFaces() << nl
            << "  " << "nInternalFaces: " << mesh.nInternalFaces() << nl;

        Info<< "----------------" << nl
            << "Patches" << nl
            << "----------------" << nl;

        forAll(patches, patchI)
        {
            const polyPatch& p = patches[patchI];

            Info<< "  " << "patch " << patchI
                << " (start: " << p.start()
                << " size: " << p.size()
                << ") name: " << p.name()
                << nl;
        }
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
