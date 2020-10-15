// Search for the appropriate blockMesh dictionary....
const word dictName("blockMeshDict");

autoPtr<IOdictionary> meshDictPtr;

{
    fileName dictPath;

    if
    (
        #if (OPENFOAM < 1812)
        args.optionReadIfPresent("dict", dictPath)
        #else
        args.readIfPresent("dict", dictPath)
        #endif
    )

    {
        // Dictionary specified on the command-line ...

        if (isDir(dictPath))
        {
            //Newer:  dictPath /= dictName;
            //Older:
            dictPath = dictPath / dictName;
        }
    }
    else if
    (
        exists
        (
            runTime.path()/runTime.constant()
            /regionPath/polyMesh::meshSubDir/dictName
        )
    )
    {
        // Dictionary present in constant polyMesh directory (old-style)

        dictPath =
            runTime.constant()
                /regionPath/polyMesh::meshSubDir/dictName;


        // Warn that constant/polyMesh/blockMeshDict was used
        // instead of system/blockMeshDict
        WarningIn(args.executable())
            << "Using the old blockMeshDict location: "
            << dictPath << nl
            << "    instead of the default location:  "
            << runTime.system()/regionPath/dictName << nl
            << endl;
    }
    else
    {
        // Assume dictionary is to be found in the system directory

        dictPath = runTime.system()/regionPath/dictName;
    }

    IOobject meshDictIO
    (
        dictPath,
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    if (!meshDictIO.typeHeaderOk<IOdictionary>(true))
    {
        FatalErrorInFunction
            << meshDictIO.objectPath() << nl
            << exit(FatalError);
    }

    Info<< "Creating block mesh from "
        //Newer: << runTime.relativePath(meshDictIO.objectPath()) << endl;
        //Older:
        << meshDictIO.objectPath() << endl;

    meshDictPtr.reset(new IOdictionary(meshDictIO));
}

//Newer:  const IOdictionary& meshDict = *meshDictPtr;

//Older:
const IOdictionary& meshDict = meshDictPtr();
