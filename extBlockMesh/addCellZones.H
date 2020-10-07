// Set any cellZones
// Note cell labelling unaffected by previous mergePatchPairs

{
    const label nZones = blocks.numZonedBlocks();
    if (nZones)
    {
        Info<< nl << "Adding cell zones" << endl;

        // Map from zoneName to cellZone index
        HashTable<label> zoneMap(nZones);

        // Cells per zone.
        List<DynamicList<label> > zoneCells(nZones);

        // Running cell counter
        label celli = 0;

        // Largest zone so far
        label freeZonei = 0;

        forAll(blocks, blocki)
        {
            const block& b = blocks[blocki];
            const word& zoneName = b.zoneName();
            const label nCellsInBlock = b.cells().size();

            if (zoneName.size())
            {
                HashTable<label>::const_iterator iter = zoneMap.find(zoneName);

                label zonei = freeZonei;

                if (iter != zoneMap.end())
                {
                    zonei = iter();
                }
                else
                {
                    zoneMap.insert(zoneName, zonei);
                    ++freeZonei;

                    Info<< "    " << zonei << '\t' << zoneName << endl;
                }


                // Fill with cell ids

                zoneCells[zonei].reserve
                (
                    zoneCells[zonei].size() + nCellsInBlock
                );

                const label endOfFill = celli + nCellsInBlock;

                for (; celli < endOfFill; ++celli)
                {
                    zoneCells[zonei].append(celli);
                }
            }
            else
            {
                celli += nCellsInBlock;
            }
        }

        List<cellZone*> cz(zoneMap.size());

        forAllConstIter(HashTable<label>, zoneMap, iter)
        {
            label zoneI = iter();

            cz[zoneI] = new cellZone
            (
                iter.key(),
                zoneCells[zoneI].shrink(),
                zoneI,
                mesh.cellZones()
            );

            // Write as cellSet for ease of processing
            cellSet cset(mesh, iter.key(), zoneCells[zoneI].shrink());
            cset.write();
        }

        mesh.pointZones().resize(0);
        mesh.faceZones().resize(0);
        mesh.cellZones().resize(0);
        mesh.addZones(List<pointZone*>(), List<faceZone*>(), cz);
    }
}


// ************************************************************************* //
