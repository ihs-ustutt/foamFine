/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

Class
    Foam::mappedFieldFixedValueFvPatchField

Group
    grpInletBoundaryConditions grpCoupledBoundaryConditions

Description
    This boundary conditions interpolates the values from a set of supplied
    points in space. If the volume flow needs to be scaled after the 
    interpolation, this can be done by defining the subdict "scaleVolumeFlow"
    in 0/U.

    Supplied data should be specified in constant/boundaryData/\<patchname\>/
      - points             : pointField of locations
      - \<time\>/\<field\> : field of values at time \<time\>

    The default mode of operation (mapMethod planarInterpolation) is to project
    the points onto a plane (constructed from the first three points) and
    construct a 2D triangulation and finds for the face centres the triangle it
    is in and the weights to the 3 vertices.

    The optional mapMethod nearest will avoid all projection and triangulation
    and just use the value at the nearest vertex.

    

Usage
    \table
        Property     | Description                      | Required | Default
        setAverage   | Use average value                | no    | false
        perturb      | Perturb points for regular geometries | no | 1e-5
        points       | Name of points file              | no    | points
        fieldTable   | Alternative field name to sample | no | this field name
        mapMethod    | Type of mapping              | no | planarInterpolation
        offset       | Offset to mapped values      | no | Zero
    \endtable

    \subDict
    // Optional scaling of the volume flow. Define either requested
    // normal velocity or requested volume flow.
    scaleVolumeFlow
    {
        normalVelocity  1.0;
        volumeFlow      1.0;
    }
    \endsubDict
    Additionally, a path to the data can be defined in the controlDict. The
    data is then copied to the case at the beginning of runtime. This requires
    the following entry in the controlDict:
    \controlDict
    mappedFieldFixedValue
    {
        pathToData      pathToData;
        absolutePath    1; //Set if absolute path is used, default 0
    }
    \endcontrolDict
    \verbatim
    <patchName>
    {
        type            mappedFieldFixedValue;
    }
    \endverbatim

See also
    Foam::fixedValueFvPatchField
    Foam::Function1Types

SourceFiles
    mappedFieldFixedValueFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef mappedFieldFixedValueFvPatchField_H
#define mappedFieldFixedValueFvPatchField_H

#include "fixedValueFvPatchFields.H"
#include "MappedFile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class mappedFieldFixedValueFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class mappedFieldFixedValueFvPatchField
:
    public fixedValueFvPatchField<Type>
{
    // Private data

        autoPtr<PatchFunction1Types::MappedFile<Type>> uniformValue_;


public:

    //- Runtime type information
    TypeName("mappedFieldFixedValue");


    // Constructors

        //- Construct from patch and internal field
        mappedFieldFixedValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        mappedFieldFixedValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given mappedFieldFixedValueFvPatchField
        //  onto a new patch
        mappedFieldFixedValueFvPatchField
        (
            const mappedFieldFixedValueFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        mappedFieldFixedValueFvPatchField
        (
            const mappedFieldFixedValueFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new mappedFieldFixedValueFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        mappedFieldFixedValueFvPatchField
        (
            const mappedFieldFixedValueFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new mappedFieldFixedValueFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<Type>&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            //virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "mappedFieldFixedValueFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
