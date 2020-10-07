autoPtr<IOdictionary> smootherDictPtr;
autoPtr<MeshSmoother> smootherPtr;

{
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
        FatalErrorInFunction
            << smootherDictIO.objectPath() << nl
            << exit(FatalError);
    }

    Info<< nl << "Initialize smoother algorithm" << nl;

    smootherDictPtr.reset(new IOdictionary(smootherDictIO));
    smootherPtr.reset(new MeshSmoother(&mesh, &(smootherDictPtr())));
}
