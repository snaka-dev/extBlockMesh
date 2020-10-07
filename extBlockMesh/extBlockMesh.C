/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    blockMesh

Description
    A multi-block mesh generator.

    Uses the block mesh description found in
    \a constant/polyMesh/blockMeshDict
    (or \a constant/\<region\>/polyMesh/blockMeshDict).

Usage

    - blockMesh [OPTION]

    \param -blockTopology \n
    Write the topology as a set of edges in OBJ format.

    \param -region \<name\> \n
    Specify an alternative mesh region.

    \param -dict \<filename\> \n
    Specify alternative dictionary for the block mesh description.

\*---------------------------------------------------------------------------*/
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

#include "Time.H"
#include "IOdictionary.H"
#include "IOPtrList.H"

#include "blockMesh.H"
#include "attachPolyTopoChanger.H"
#include "emptyPolyPatch.H"

#include "cellSet.H"

#include "argList.H"
#include "OSspecific.H"
#include "OFstream.H"

#include "Pair.H"
#include "slidingInterface.H"

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
    argList::addOption("dict", "file", "Alternative blockMeshDict");

    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"

    word regionName(polyMesh::defaultRegion);
    word regionPath;

    // Check if the region is specified otherwise mesh the default region
    if (args.optionReadIfPresent("region", regionName))
    {
        Info<< nl << "Generating mesh for region " << regionName << endl;
        regionPath = regionName;
    }

    // Locate appropriate blockMeshDict
    #include "findBlockMeshDict.H"

    blockMesh blocks(meshDict, regionName);
    blocks.verbose(false);

    word defaultFacesName = "defaultFaces";
    word defaultFacesType = emptyPolyPatch::typeName;

    Info<< nl << "Creating polyMesh from blockMesh" << endl;

    polyMesh mesh
    (
        IOobject
        (
            regionName,
            runTime.constant(),
            runTime
        ),
        pointField(blocks.points()),
        blocks.cells(),
        blocks.patches(),
        blocks.patchNames(),
        blocks.patchDicts(),
        defaultFacesName,
        defaultFacesType
    );

    // ########################################################################
    // Smoothing

    Info<< nl << "Initialize smoother algorithm" << nl;

    const word smootherDictName("smootherDict");
    IOobject smootherDictIO
    (
        smootherDictName,
        runTime.system(),
        "",
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    if (!smootherDictIO.typeHeaderOk<IOdictionary>(true))
    {
        FatalErrorIn(args.executable())
            << "Cannot open mesh smoothing file\n    "
            << smootherDictIO.objectPath()
            << nl
            << exit(FatalError);
    }

    IOdictionary smootherDict(smootherDictIO);
    MeshSmoother meshSmoother(&mesh, &smootherDict, &blocks);

    if (args.optionFound("writeStep"))
    {
        meshSmoother.updateAndWrite
        (
            regionName,
            defaultFacesName,
            defaultFacesType,
            runTime
        );
    }
    else
    {
        meshSmoother.update();

        // Reset mesh directory to constant/polyMesh !!! Hard to find :P
        mesh.setInstance(runTime.constant());
    }

    // End of smoothing
    //##########################################################################

    // Merge patch pairs (dictionary entry "mergePatchPairs")
    #include "mergePatchPairs.H"

    // Set any cellZones
    #include "addCellZones.H"


    // Set the precision of the points data to 10
    IOstream::defaultPrecision(max(10u, IOstream::defaultPrecision()));

    Info<< nl << "Writing polyMesh" << endl;

    // #########################################################################

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

    // #########################################################################

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
