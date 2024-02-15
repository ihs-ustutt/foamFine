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

\*---------------------------------------------------------------------------*/

#include "pressureTotCarnot.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(pressureTotCarnot, 0);
    addToRunTimeSelectionTable(functionObject, pressureTotCarnot, dictionary);
}
}


const Foam::Enum
<
    Foam::functionObjects::pressureTotCarnot::flowDirection
>
Foam::functionObjects::pressureTotCarnot::flowDirectionNames
({
    // Normal operations
    { flowDirection::diXPos, "x" },
    { flowDirection::diXNeg, "-x" },
    { flowDirection::diYPos, "y" },
    { flowDirection::diYNeg, "-y" },
    { flowDirection::diZPos, "z" },
    { flowDirection::diZNeg, "-z" },

});



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::word Foam::functionObjects::pressureTotCarnot::resultName() const
{
    word rName;

    rName = "totalCarnot(" + fieldName_ + ")";

    return rName;
}



Foam::tmp<Foam::volScalarField>
Foam::functionObjects::pressureTotCarnot::calcPressure
(
    const volScalarField& p
) const
{
    // Initialise 
    auto tresult =
        tmp<volScalarField>::New
        (
            IOobject
            (
                scopedName("p"),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ
            ),
            mesh_,
            dimensionedScalar("p", dimPressure/dimDensity, 0)
        );

    volScalarField& result = tresult.ref();
    
    switch (flowDirection_)
    {
        // Calc total pressure, use kinetic energy only if backflow occurs
        case diXPos:
        {
            result += 
                p + neg(lookupObject<volVectorField>(UName_).component(0)) *
                0.5*magSqr(lookupObject<volVectorField>(UName_)); 
        }
        break;

        case diXNeg:
        {
            result += 
                p + pos(lookupObject<volVectorField>(UName_).component(0)) *
                0.5*magSqr(lookupObject<volVectorField>(UName_)); 
        }
        break;

        case diYPos:
        {
            result += 
                p + neg(lookupObject<volVectorField>(UName_).component(1)) *
                0.5*magSqr(lookupObject<volVectorField>(UName_));;
        }
        break;

        case diYNeg:
        {
            result += 
                p + pos(lookupObject<volVectorField>(UName_).component(1)) *
                0.5*magSqr(lookupObject<volVectorField>(UName_)); 
        }
        break;

        case diZPos:
        {
            result += 
                p + neg(lookupObject<volVectorField>(UName_).component(2)) *
                0.5*magSqr(lookupObject<volVectorField>(UName_));;
        }
        break;

        case diZNeg:
        {
            result += 
                p + pos(lookupObject<volVectorField>(UName_).component(2)) *
                0.5*magSqr(lookupObject<volVectorField>(UName_));
        }
        break;

        default:
        {
             FatalErrorInFunction
                << "Entry flowDirection not found which is "
                << "needed for patchAverage" << nl
                << abort(FatalError);
        }
        break;

    }

    return tresult;
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::pressureTotCarnot::calc()
{
    if (foundObject<volScalarField>(fieldName_))
    {
        const volScalarField& p = lookupObject<volScalarField>(fieldName_);

        auto tp = tmp<volScalarField>::New
        (
            IOobject
            (
                resultName_,
                p.mesh().time().timeName(),
                p.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::REGISTER
            ),
            std::move(calcPressure(p))
        );

        return store(resultName_, tp);
    }

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::pressureTotCarnot::pressureTotCarnot
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict, "p"),
    flowDirection_(diZPos),
    UName_("U")
{
    
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::pressureTotCarnot::read(const dictionary& dict)
{

    Info<< type() << " " << name() << ":" << nl;

    fieldExpression::read(dict);

    UName_   = dict.getOrDefault<word>("U", "U");
    
    if (dict.found("flowDirection"))
    {
        flowDirection_ = flowDirectionNames.read(dict.lookup("flowDirection"));   
    }
    Info << "   Using flow direction: " 
         << flowDirectionNames.get(flowDirection_) << nl;

    resultName_ = dict.getOrDefault<word>("result", resultName());

    Info<< endl;

    return true;
}


// ************************************************************************* //
