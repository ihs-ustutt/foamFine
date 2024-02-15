/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2018 OpenCFD Ltd.
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

#include "mappedFieldFixedValueFvPatchField.H"
#include "Time.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::mappedFieldFixedValueFvPatchField<Type>::
mappedFieldFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    uniformValue_()
{}


template<class Type>
Foam::mappedFieldFixedValueFvPatchField<Type>::
mappedFieldFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF, dict, false),
    uniformValue_
    (
        new PatchFunction1Types::MappedFile<Type>
        (
            p.patch(),
            "uniformValue",
            dict,
            iF.name(),          // field table name
            true                // face values
        )
    )
{
    if (this->db().time().controlDict().found("mappedFieldFixedValue"))
    {
        const string pathToCopy = 
            this->db().time().globalPath()+"/"+this->db().time().constant();
     
        if (!isDir(pathToCopy+"/boundaryData"))
        {
            string pathToData;
            bool absolutePath = 0;
            
            this->db().time().controlDict().subDict("mappedFieldFixedValue").
                     readEntry
                    (
                        "pathToData", pathToData
                    );  
                    
            this->db().time().controlDict().subDict("mappedFieldFixedValue").
                     readIfPresent
                    (
                        "absolutePath", absolutePath
                    );
                    
            if (!absolutePath)
            {
                pathToData = this->db().time().globalPath() + "/" + pathToData;
            }

            Info << "MappedFieldFixedValue: Copy data from " << pathToData
                    << " to " << pathToCopy << endl;
            
            cp(pathToData, pathToCopy);
        }
        else
        {
            Info << "MappedFieldFixedValue: Folder " 
                    << pathToCopy+"/boundaryData" 
                    << " already exists, skipping copy" << endl;
        }
    }
    if (dict.found("scaleVolumeFlow"))
    {

        const scalar t = this->db().time().timeOutputValue();
        Field<Type > mappedVelocity = uniformValue_->value(t);
        
        const scalar mappedPhi = mag(
            gSum(
                mappedVelocity.component(0) * p.Sf().component(0) + 
                mappedVelocity.component(1) * p.Sf().component(1) + 
                mappedVelocity.component(2) * p.Sf().component(2)
            )
        );
        
        scalar scaleFactor = 1.0;
        
        if (dict.subDict("scaleVolumeFlow").found("normalVelocity"))
        {
            scalar normalVelocity;
            dict.subDict("scaleVolumeFlow").readEntry
                (
                    "normalVelocity", normalVelocity
                );
            
            scaleFactor = gSum(normalVelocity*p.magSf())/mappedPhi;
        }
        else if (dict.subDict("scaleVolumeFlow").found("volumeFlow"))
        {
            scalar volumeFlow;
            dict.subDict("scaleVolumeFlow").readEntry
                (
                    "volumeFlow", volumeFlow
                );
            
            scaleFactor = volumeFlow/mappedPhi;
        }
        else
        {
            FatalIOErrorIn
            (
                "Foam::mappedFieldFixedValueFvPatchField<Type>::\n"
                "mappedFieldFixedValueFvPatchField\n"
                "(\n"
                "    const fvPatch& p,\n"
                "    const DimensionedField<Type, volMesh>& iF,\n"
                "    const dictionary& dict\n"
                ")\n",
                dict.subDict("scaleVolumeFlow")
            )   << "Cannot find entry normalVelocity or volumeFlow." << nl
                << "Please specify one of them in subDict scaleVolumeFlow"
                << abort(FatalIOError);
        }
        
        fvPatchField<Type>::operator==(mappedVelocity * scaleFactor);	
    }
    else if (dict.found("value"))
    {
        fvPatchField<Type>::operator==(Field<Type>("value", dict, p.size()));
    }
    else
    {
        const scalar t = this->db().time().timeOutputValue();
        fvPatchField<Type>::operator==(uniformValue_->value(t));
    }
}


template<class Type>
Foam::mappedFieldFixedValueFvPatchField<Type>::
mappedFieldFixedValueFvPatchField
(
    const mappedFieldFixedValueFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    uniformValue_
    (
        new PatchFunction1Types::MappedFile<Type>
        (
            ptf.uniformValue_(),
            p.patch()
        )
    )
{}


template<class Type>
Foam::mappedFieldFixedValueFvPatchField<Type>::
mappedFieldFixedValueFvPatchField
(
    const mappedFieldFixedValueFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    uniformValue_
    (
        new PatchFunction1Types::MappedFile<Type>
        (
            ptf.uniformValue_(),
            this->patch().patch()
        )
    )
{}


template<class Type>
Foam::mappedFieldFixedValueFvPatchField<Type>::
mappedFieldFixedValueFvPatchField
(
    const mappedFieldFixedValueFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    uniformValue_
    (
        new PatchFunction1Types::MappedFile<Type>
        (
            ptf.uniformValue_(),
            this->patch().patch()
        )
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::mappedFieldFixedValueFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchField<Type>::autoMap(m);
    uniformValue_().autoMap(m);
}


template<class Type>
void Foam::mappedFieldFixedValueFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField<Type>::rmap(ptf, addr);

    const mappedFieldFixedValueFvPatchField<Type>& tiptf =
        refCast<const mappedFieldFixedValueFvPatchField<Type>>(ptf);

    uniformValue_().rmap(tiptf.uniformValue_(), addr);
}


template<class Type>
void Foam::mappedFieldFixedValueFvPatchField<Type>::write
(
    Ostream& os
) const
{
    fvPatchField<Type>::write(os);
    uniformValue_->writeData(os);
    this->writeEntry("value", os);
}


// ************************************************************************* //
