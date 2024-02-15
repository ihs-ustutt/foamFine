/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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
    Foam::functionObjects::pressureTotCarnot

Group
    grpFieldFunctionObjects

Description
    Calculates a total pressure field which can be used to calculate the
    total pressure at an outlet taking into account the Carnot losses.
    Therefore, the total pressure is the sum of static pressure and kinetic
    energy if backflow occurs, otherwise only the static pressure is used.  

    - total pressure Carnot
        \f[
            if backflow: p_0 = p + 0.5 |\vec U|^2
            else: p_0 = p 
        \f]

    where
    \vartable
      \vec U      | Velocity [m/s]
      p           | Static pressure (p/rho)[m2/s2]
      p_0         | Total pressureTotCarnot (p/rho)[m2/s2]
    \endvartable

    Operands:
    \table
      Operand        | Type           | Location
      input          | volScalarField | $FOAM_CASE/\<time\>/\<inpField\>
      output file    | -              | -
      output field   | volScalarField | $FOAM_CASE/\<time\>/\<outField\>
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    pressureTotCarnot1
    {
        // Mandatory entries (unmodifiable)
        type            pressureTotCarnot;
        libs            (fieldFunctionObjects);

        // Mandatory entries (runtime modifiable)
        mode            \<option\>;

        // Optional entries (runtime modifiable)
        flowDirection   z;
        field           \<pName\>;
        U               \<UName\>;
        rho             <rhoName>;


        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Req'd | Dflt
      type         | Type name: pressureTotCarnot       | word |  yes  | -
      libs         | Library name: fieldFunctionObjects | word |  yes  | -
      mode         | Calculation mode (see below)       | word |  yes  | -
      flowDirection| Name of the flow direction         | word |  no   | z
      field        | Name of the pressure field         | word |  no   | p
      U            | Name of the velocity field         | word |  no   | U
      rho          | Name of the density field          | word |  no   | rho
    \endtable

    Options for the \c flowDirection entry:
    \verbatim
        x      | positiv x-direction
        -x     | negativ x-direction
        y      | positiv y-direction
        -y     | negativ y-direction
        z      | positiv z-direction
        z      | negativ z-direction
    \endverbatim

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link fieldExpression.H \endlink

    Usage by the \c postProcess utility is not available.

See also
    - Foam::functionObject
    - Foam::functionObjects::fieldExpression
    - Foam::functionObjects::fvMeshFunctionObject
    - ExtendedCodeGuide::functionObjects::field::pressureTotCarnot

SourceFiles
    pressureTotCarnot.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_pressureTotCarnot_H
#define functionObjects_pressureTotCarnot_H

#include "fieldExpression.H"
#include "volFieldsFwd.H"
#include "dimensionedVector.H"
#include "dimensionedScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                        Class pressureTotCarnot Declaration
\*---------------------------------------------------------------------------*/

class pressureTotCarnot
:
    public fieldExpression
{
public:

    // Public Data Types
        
        enum flowDirection
        {
            diXPos,             
            diXNeg,
            diYPos,
            diYNeg,
            diZPos,
            diZNeg,
        };
        
        //- Operation type names
        static const Enum<flowDirection> flowDirectionNames;



private:

    // Private Data
        
        flowDirection flowDirection_;

        //- Hydrostatic constribution mode
        //hydrostaticMode hydrostaticMode_;

        //- Name of velocity field
        word UName_;

        //- Name of density field
        word rhoName_;


    // Private Member Functions

        //- Return the name of the derived pressureTotCarnot field
        word resultName() const;

        //- Calculate and return the pressureTotCarnot
        tmp<volScalarField> calcPressure
        (
            const volScalarField& p
        ) const;

        //- Calculate the derived pressureTotCarnot field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("pressureTotCarnot");


    // Constructors

        //- Construct from Time and dictionary
        pressureTotCarnot
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );

        //- No copy construct
        pressureTotCarnot(const pressureTotCarnot&) = delete;

        //- No copy assignment
        void operator=(const pressureTotCarnot&) = delete;


    //- Destructor
    virtual ~pressureTotCarnot() = default;


    // Member Functions

        //- Read the pressureTotCarnot data
        virtual bool read(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
