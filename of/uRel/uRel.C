/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    uRel

Description
    For each time: calculate the relative Velocity in MRFZones.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "MRFZone.H"
#include "MRFZoneList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();

#   include "setRootCase.H"
#   include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

#   include "createMesh.H"

    IOdictionary mrfDict
            (
                IOobject
                (
                    "MRFProperties",
                    mesh.time().constant(),
                    mesh.time(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );

    MRFZoneList mrfZones(mesh, mrfDict);

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time = " << runTime.timeName() << endl;


        IOobject Uheader
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );


        // Check U exist
        if (Uheader.typeHeaderOk<volVectorField>(true,true,false))
        {
            mesh.readUpdate();


            Info<< "    Reading U" << endl;
            volVectorField U(Uheader, mesh);

            Info<< "    Calculating uRel" << endl;
            if (U.dimensions() == dimensionSet(0, 1, -1, 0, 0))
            {
                volVectorField uRel
                (
                    IOobject
                    (
                        "uRel",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ
                    ),
                    U
		);

		mrfZones.makeRelative(uRel);
		uRel.write();
            }
        }
        else
        {
            Info<< "    No U input	" << endl;
        }

        Info<< endl;
    }

    return 0;
}


// ************************************************************************* //
