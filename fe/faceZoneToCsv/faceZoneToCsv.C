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
 faceZoneToCSV

Description
    faceZonesToCSV.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "cyclicPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "addRegionOption.H"
    Foam::timeSelector::addOptions();
    Foam::argList::validArgs.append("faceSetName");
    Foam::argList::validArgs.append("fieldName");
#   include "setRootCase.H"
#   include "createTime.H"
    instantList timeDirs = Foam::timeSelector::select0(runTime, args);
#   include "createNamedMesh.H"

    Foam::word faceZoneName( args.additionalArgs()[0] );
    Foam::word fieldName( args.additionalArgs()[1] );

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

        IOobject fieldHeader(
          fieldName,
          runTime.timeName(),
          mesh,
          IOobject::MUST_READ
        );

        IOobject phiHeader(
          "phi",
          runTime.timeName(),
          mesh,
          IOobject::MUST_READ
        );
        
        // Check wheather field exists
        if (fieldHeader.headerOk() && phiHeader.headerOk() ) {
          mesh.readUpdate();
          
          //get faceZoneID of faceZone 'faceZoneName'
          label faceZoneID = mesh.faceZones().findZoneID(faceZoneName);
          if (faceZoneID < 0)
          {
            FatalError
            << "Unable to find faceZone " << faceZoneName << nl
            << exit(FatalError);
          }
            
          //collect FaceIDs contained in faceZone 'faceZoneID'
          const labelList& zoneFaces = 
          mesh.faceZones()[faceZoneID].zoneMesh()[faceZoneID];   

          //determine if the fieldvalue is a 
          //    scalar -> volScalarField
          //    vector -> volVectorField
          if (fieldHeader.headerClassName() == volScalarField::typeName) {
          //
          // read fields
          //
          volScalarField fieldField(fieldHeader, mesh);
          //interpolate fieldValues from cell to the cellfaces
          surfaceScalarField fieldFaceField =
            fvc::interpolate(fieldField, fieldName) ;
          surfaceScalarField phiField(phiHeader, mesh);
                
          //set output-filename 
          std::string fstr(faceZoneName+"_"+fieldName+"_"+runTime.timeName()+".csv");
          
          
          //Write csv-header
          Info<< "    Write " << fstr << endl;              
          std::ofstream of;
          of.open( fstr.c_str(), std::ios::out | std::ios::trunc );
                
          of 
            << "# 1 x" << std::endl
            << "# 2 y" << std::endl
            << "# 3 z" << std::endl
            << "# 4 value" << std::endl
            << "# 5 sfX " << std::endl          
            << "# 6 sfY " << std::endl
            << "# 7 sfZ " << std::endl
            << "# 8 phi" << std::endl;
                               
          //write fieldvalue for each face in faceZone      
          forAll(zoneFaces, ii){
            point pp = mesh.Cf().internalField()[zoneFaces[ii]];
            scalar ss = fieldFaceField.internalField()[zoneFaces[ii]];
            vector sf = mesh.Sf().internalField()[zoneFaces[ii]];
            scalar phi = phiField.internalField()[zoneFaces[ii]];

                                                     
            of
              << pp.x() << ", "
              << pp.y() << ", "
              << pp.z() << ", "
              << ss << ", "
              << sf.x() << ", "
              << sf.y() << ", "
              << sf.z() << ", "
              << phi << std::endl;
          }

          of.close();
                
        }
        else if (fieldHeader.headerClassName() == volVectorField::typeName) {
          //
          // read fields
          //          
          volVectorField fieldField(fieldHeader, mesh);

          //interpolate fieldValues from cell to the cellfaces
          surfaceVectorField fieldFaceField = 
            fvc::interpolate(fieldField, fieldName) ;
          surfaceScalarField phiField(phiHeader, mesh);

          //set output-filename
          std::string fstr( 
            faceZoneName+"_"+fieldName+"_"+runTime.timeName()+".csv" 
          );
          
          //Write csv-header
          Info<< "    Write " << fstr << endl;              
          std::ofstream of;
          of.open( fstr.c_str(), std::ios::out | std::ios::trunc );

          of 
            << "# 1  x" << std::endl
            << "# 2  y" << std::endl
            << "# 3  z" << std::endl
            << "# 4  valueX" << std::endl
            << "# 5  valueY" << std::endl
            << "# 6  valueZ" << std::endl
            << "# 7  sfX " << std::endl          
            << "# 8  sfY " << std::endl
            << "# 9  sfZ " << std::endl
            << "# 10 phi" << std::endl;
              
          //write fieldvalue for each face in faceZone                    
          forAll(zoneFaces, ii){
            point pp = mesh.Cf().internalField()[zoneFaces[ii]];
            vector ss = fieldFaceField.internalField()[zoneFaces[ii]]; 

            vector sf = mesh.Sf().internalField()[zoneFaces[ii]];
            scalar phi = phiField.internalField()[zoneFaces[ii]];

                
            of
              << pp.x() << ", "
              << pp.y() << ", "
              << pp.z() << ", "
              << ss.x() << ", "
              << ss.y() << ", "
              << ss.z() << ", "
              << sf.x() << ", "
              << sf.y() << ", "
              << sf.z() << ", "
              << phi << std::endl;
          }
              
          of.close();
        }
        else if (fieldHeader.headerClassName() == volSymmTensorField::typeName) {
          //
          // read fields
          //          
          volSymmTensorField fieldField(fieldHeader, mesh);

          //interpolate fieldValues from cell to the cellfaces
          surfaceSymmTensorField fieldFaceField = 
            fvc::interpolate(fieldField, fieldName) ;
          surfaceScalarField phiField(phiHeader, mesh);

          //set output-filename
          std::string fstr( 
            faceZoneName+"_"+fieldName+"_"+runTime.timeName()+".csv" 
          );
          
          //Write csv-header
          Info<< "    Write " << fstr << endl;              
          std::ofstream of;
          of.open( fstr.c_str(), std::ios::out | std::ios::trunc );

          of 
            << "# 1  x" << std::endl
            << "# 2  y" << std::endl
            << "# 3  z" << std::endl
            << "# 4  valueXX" << std::endl
            << "# 5  valueXY" << std::endl
            << "# 6  valueXZ" << std::endl
            << "# 7  valueYY" << std::endl
            << "# 8  valueYZ" << std::endl
            << "# 9  valueZZ" << std::endl
            << "# 10 sfX " << std::endl          
            << "# 11 sfY " << std::endl
            << "# 12 sfZ " << std::endl
            << "# 13 phi" << std::endl;
              
          //write fieldvalue for each face in faceZone                    
          forAll(zoneFaces, ii){
            point pp = mesh.Cf().internalField()[zoneFaces[ii]];
            symmTensor ss = fieldFaceField.internalField()[zoneFaces[ii]]; 

            vector sf = mesh.Sf().internalField()[zoneFaces[ii]];
            scalar phi = phiField.internalField()[zoneFaces[ii]];

                
            of
              << pp.x() << ", "
              << pp.y() << ", "
              << pp.z() << ", "
              << ss.component(symmTensor::XX) << ", "
              << ss.component(symmTensor::XY) << ", "
              << ss.component(symmTensor::XZ) << ", "
              << ss.component(symmTensor::YY) << ", "
              << ss.component(symmTensor::YZ) << ", "
              << ss.component(symmTensor::ZZ) << ", "
              << sf.x() << ", "
              << sf.y() << ", "
              << sf.z() << ", "
              << phi
              << std::endl;
          }
              
          of.close();
        }
        else {
          FatalError
            << "Only possible to write "
            << volScalarField::typeName 
            << "and " << volVectorField::typeName
            << " and " <<volSymmTensorField::typeName
            << nl << exit(FatalError);
        }
      }
      else {
        Info
        << "    Header not ok." << endl
        << "    fieldHeader [ " 
        << fieldName 
        << " ]: headerOk = " 
        << fieldHeader.headerOk() << endl
        << "    phiHeader [ phi ]: headerOk = " 
        << phiHeader.headerOk() << endl;
    }
        
    Info<< endl;
  }
    
  Info<< "End\n" << endl;

  return 0;
}



// ************************************************************************* //
