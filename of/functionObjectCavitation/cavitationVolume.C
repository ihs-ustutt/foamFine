/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd.
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

#include "cavitationVolume.H"
#include "fvMesh.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace fieldValues
{
    defineTypeNameAndDebug(cavitationVolume, 0);
    addToRunTimeSelectionTable(fieldValue, cavitationVolume, runTime);
    addToRunTimeSelectionTable(functionObject, cavitationVolume, dictionary);
}
}
}

const Foam::Enum
<
    Foam::functionObjects::fieldValues::cavitationVolume::operationType
>
Foam::functionObjects::fieldValues::cavitationVolume::operationTypeNames_
({
    // Normal operations
    { operationType::opNone, "none" },
    { operationType::opSum, "sum" },

});

const Foam::Enum
<
    Foam::functionObjects::fieldValues::cavitationVolume::heightAxis
>
Foam::functionObjects::fieldValues::cavitationVolume::heightAxisNames_
({
    // Normal operations
    { heightAxis::hXPos, "x" },
    { heightAxis::hXNeg, "-x" },
    { heightAxis::hYPos, "y" },
    { heightAxis::hYNeg, "-y" },
    { heightAxis::hZPos, "z" },
    { heightAxis::hZPos, "-z" },

});

const Foam::Enum
<
    Foam::functionObjects::fieldValues::cavitationVolume::flowDirection
>
Foam::functionObjects::fieldValues::cavitationVolume::flowDirectionNames_
({
    // Normal operations
    { flowDirection::diXPos, "x" },
    { flowDirection::diXNeg, "-x" },
    { flowDirection::diYPos, "y" },
    { flowDirection::diYNeg, "-y" },
    { flowDirection::diZPos, "z" },
    { flowDirection::diZNeg, "-z" },

});



// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


void Foam::functionObjects::fieldValues::cavitationVolume::writeFileHeader
(
    Ostream& os
) const
{
    volRegion::writeFileHeader(*this, os);

    writeCommented(os, "Time");

    // TBD: add in postOperation information?

    /*for (const word& fieldName : fields_)
    {
        os  << tab << operationTypeNames_[operation_]
            << "(" << fieldName << ")";
    }*/

    os  << endl;
}


Foam::label Foam::functionObjects::fieldValues::cavitationVolume::writeAll
(
    const scalarField& V
)
{
    label nProcessed = 0;
    
    if (fields_[0] != "p" || fields_.size() != 1)
    {
        FatalErrorInFunction
            << "Only supported for Field p" << nl
            << abort(FatalError);
    }

    for (const word& fieldName : fields_)
    {
        if
        (
            writeValues<scalar>(fieldName, V)
        )
        {
            ++nProcessed;
        }
        else
        {
            WarningInFunction
                << "Requested field " << fieldName
                << " not found in database and not processed"
                << endl;
        }
    }

    return nProcessed;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldValues::cavitationVolume::cavitationVolume
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldValue(name, runTime, dict, typeName),
    volRegion(fieldValue::mesh_, dict),
    operation_(operationTypeNames_.get("operation", dict)),
    useHeightAxis_(false),
    installDepth_(0)
    
{
    read(dict);
    writeFileHeader(file());
}


Foam::functionObjects::fieldValues::cavitationVolume::cavitationVolume
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    fieldValue(name, obr, dict, typeName),
    volRegion(fieldValue::mesh_, dict),
    operation_(operationTypeNames_.get("operation", dict)),
    useHeightAxis_(false),
    installDepth_(0)
    
{
    
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fieldValues::cavitationVolume::read
(
    const dictionary& dict
)
{
    fieldValue::read(dict);
    Info<< type() << ' ' << name() << ':' << nl;
    
    if (dict.readIfPresent("patchAverageName", patchAverageName_))
    {
        patchAverage_ = true;
        Info << "   Using patch average for patch: " 
            << patchAverageName_ << endl;
    }
    else patchAverage_ = false; 
    
    dict.readIfPresent("installDepth", installDepth_);
    Info << "   Using installDepth: " << installDepth_ << nl;
    
    if (dict.found("heightAxis"))
    {
        heightAxis_ = heightAxisNames_.read(dict.lookup("heightAxis"));
        useHeightAxis_ = true;
        Info << "   Using height axis: " << dict.lookup("heightAxis") << nl;
    }
    if (dict.found("flowDirection"))
    {
        flowDirection_ = flowDirectionNames_.read(dict.lookup("flowDirection"));
        Info << "   Using flow direction: " 
            << dict.lookup("flowDirection") << nl;
    }
    
    Info << nl << endl;
    
    return true;
}


bool Foam::functionObjects::fieldValues::cavitationVolume::write()
{
    volRegion::update();        // Ensure cached values are valid

    fieldValue::write();
    
    if (Pstream::master())
    {
        writeCurrentTime(file());
    }

    scalarField V;
    
    V = filterField(fieldValue::mesh_.V());

    // Process the fields
    writeAll(V);

    if (Pstream::master())
    {
        file()<< endl;
    }

    Log << endl;

    return true;
}


// ************************************************************************* //
