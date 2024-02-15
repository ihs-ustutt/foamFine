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
 eulerHead

Description
    Calculates the integral of r x U dphi.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "cyclicPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
#   include "addRegionOption.H"
    timeSelector::addOptions();
    argList::validArgs.append("patchName");
    argList::validOptions.insert("origin", "origin");
    argList::validOptions.insert("rotationAxis", "rotationAxis");
    argList::validOptions.insert("omega", "omega");
#   include "setRootCase.H"
#   include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
#   include "createNamedMesh.H"

    word patchName( args.additionalArgs()[0] );
    
    //
    // options
    //
    point origin(0,0,0);
    vector z_(0,0,1);
    vector omega_(0,0,0);
    if (args.optionFound("origin")) {
      origin = point( args.optionLookup("origin")() );
    }
    if (args.optionFound("rotationAxis")) {
      z_ = vector( args.optionLookup("rotationAxis")() );
    }    
    if (args.optionFound("omega")) {
      omega_ = vector( args.optionLookup("omega")() );
    }    
            
    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

        IOobject UHeader
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        IOobject phiHeader
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );
        
        // Check field exists
        if (UHeader.headerOk() && phiHeader.headerOk() )
        {
            mesh.readUpdate();

            label patchi = mesh.boundaryMesh().findPatchID(patchName);
            if (patchi < 0)
            {
                FatalError
                    << "Unable to find patch " << patchName << nl
                    << exit(FatalError);
            }


            Info << "   origin = " << origin << endl;
            Info << "   rotationAxis = " << z_ << endl;
            Info << "   omega = " << omega_ << endl;
            
            // norm z
            z_ = z_ / mag(z_);
            
            // defining vector origin to face center
            vectorField d_ = (mesh.Cf().boundaryField()[patchi] - origin);
            
            //
            // z coordinate
            //
            Info<< "    Calculate z " << nl;
            scalarField z = (z_ & d_);
            
            // radius vector
            Info<< "    Calculate r_ " << nl;
            vectorField r_ =  d_ - z * z_;
          
            
            //
            // read fields
            //
            Info 
              << "    Reading " << volVectorField::typeName 
              << " U" << endl;
            volVectorField UField(UHeader, mesh);
            Info 
              << "    Reading " << surfaceScalarField::typeName 
              << " phi" << endl; 
            surfaceScalarField phiField(phiHeader, mesh);

            //
            // output
            //
            Info
                << "    min( mag(r_) ) = " << min( mag(r_) ) << nl
                << "    max( mag(r_) ) = " << max( mag(r_) ) << nl
                << "    min( mag(z) ) = " << min( z ) << nl
                << "    max( mag(z) ) = " << max( z ) << nl
                << "    Integral of phi of patch "
                << patchName << '[' << patchi << ']' << " = "
                << gSum
                   (
                       phiField.boundaryField()[patchi]
                   )
                << nl
                << "    Integral of magSf of patch "
                << patchName << '[' << patchi << ']' << " = "
                << gSum
                  (
                      mesh.magSf().boundaryField()[patchi]
                  )
                << nl
                << "    Integral of r ^ (r ^ omega ) over phi of patch "
                << patchName << '[' << patchi << ']' << " = "
                << gSum
                  (
                      phiField.boundaryField()[patchi]
                     * ( r_ ^ (r_ ^ omega_) )
                  )
                << nl
                << "    Integral of r ^ U over phi of patch "
                << patchName << '[' << patchi << ']' << " = "
                << gSum
                  (
                      phiField.boundaryField()[patchi]
                     * (r_ ^ UField.boundaryField()[patchi])
                  )
                << nl
                << "    Integral of d_ ^ U over phi of patch "
                << patchName << '[' << patchi << ']' << " = "
                << gSum
                  (
                      phiField.boundaryField()[patchi]
                     * (d_ ^ UField.boundaryField()[patchi])
                  )
               << nl;
        }
        else
        {
            Info<< "    Header not ok." << endl;
        }

        Info<< endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
