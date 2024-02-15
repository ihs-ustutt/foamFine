/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2021 OpenCFD Ltd.
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
#include "volFields.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
bool Foam::functionObjects::fieldValues::cavitationVolume::validField
(
    const word& fieldName
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;
    typedef typename VolFieldType::Internal IntVolFieldType;

    return
    (
        obr_.foundObject<VolFieldType>(fieldName)
     || obr_.foundObject<IntVolFieldType>(fieldName)
    );
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::functionObjects::fieldValues::cavitationVolume::getFieldValues
(
    const word& fieldName,
    const bool mandatory
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;
    typedef typename VolFieldType::Internal IntVolFieldType;

    if (obr_.foundObject<VolFieldType>(fieldName))
    {
        return filterField(obr_.lookupObject<VolFieldType>(fieldName));
    }
    else if (obr_.foundObject<IntVolFieldType>(fieldName))
    {
        return filterField(obr_.lookupObject<IntVolFieldType>(fieldName));
    }

    if (mandatory)
    {
        FatalErrorInFunction
            << "Field " << fieldName << " not found in database" << nl
            << abort(FatalError);
    }

    return tmp<Field<Type>>::New();
}


template<class Type>
Type Foam::functionObjects::fieldValues::cavitationVolume::processValues
(
    const Field<Type>& values,
    const scalarField& V
) const
{
    Type result = Zero;
    
    switch (operation_)
    {
        case opNone:
        {
            break;
        }
        
        case opSum:
        {
            const volVectorField& cellCenter = mesh_.C();
            
            // Get values of the axis defined as height axis
            scalarField heightAxisValue(values.size(), pTraits<Type>::zero);
            if (useHeightAxis_)
            {
                switch (heightAxis_)
                {
                    case hXPos:
                    {
                        heightAxisValue = cellCenter.component(0);
                    }
                    break;
                    
                    case hXNeg:
                    {
                        heightAxisValue = -1*cellCenter.component(0);
                    }
                    break;
                    
                    case hYPos:
                    {
                        heightAxisValue = cellCenter.component(1);
                    }
                    break;
                    
                    case hYNeg:
                    {
                        heightAxisValue = -1*cellCenter.component(1);
                    }
                    break;
                    
                    case hZPos:
                    {
                        heightAxisValue = cellCenter.component(2);
                    }
                    break;
                    
                    case hZNeg:
                    {
                        heightAxisValue = -1*cellCenter.component(2);
                    }
                    break;
                }
            }
            
            scalar patchAverage = 0;
            scalar maxDistance = 1E+30;
            scalarField distanceValue(values.size(), -1E+30);
            label sign = 1;
            
            // Calc area averaged pressure at outlet patch
            if (patchAverage_)
            {
                // Get coordinate values in flow direction
                label flowDirectionLabel = 0;
                switch (flowDirection_)
                {
                    case diXPos:
                    {
                        flowDirectionLabel = 0;
                        distanceValue = cellCenter.component(0);
                    }
                    break;
                    
                    case diXNeg:
                    {
                        flowDirectionLabel = 0;
                        distanceValue = cellCenter.component(0);
                        sign = -1;
                    }
                    break;
                    
                    case diYPos:
                    {
                        flowDirectionLabel = 1;
                        distanceValue = cellCenter.component(1);
                    }
                    break;
                    
                    case diYNeg:
                    {
                        flowDirectionLabel = 1;
                        distanceValue = cellCenter.component(1);
                        sign = -1;
                    }
                    break;
                    
                    case diZPos:
                    {
                        flowDirectionLabel = 2;
                        distanceValue = cellCenter.component(2);
                    }
                    break;
                    
                    case diZNeg:
                    {
                        flowDirectionLabel = 2;
                        distanceValue = cellCenter.component(2);
                        sign = -1;
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
                
                // Get pressure field at outlet patch and calc averaged pressure
                const label patchi = 
                    mesh_.boundaryMesh().findPatchID(patchAverageName_);
                
                if (patchi == -1)
                {
                    FatalErrorInFunction
                        << "Patch " << patchAverageName_ << " not found" << nl
                        << abort(FatalError);
                }
                
                const polyPatch& pp = mesh_.boundaryMesh()[patchi];              

                maxDistance = 
                    gMin( pp.faceCentres().component(flowDirectionLabel) );
                
                typedef GeometricField<Type, fvPatchField, volMesh> vf;
                const vf& fld = lookupObject<vf>("p");
                
                const scalarField& pValue = fld.boundaryField()[patchi];
                
                patchAverage = gSum( pValue * pp.magFaceAreas() )
                        / gSum( pp.magFaceAreas() );
                
            }

            
            Field<Type> weights(values.size(), pTraits<Type>::zero);
            
            // Correct the pressure
            scalarField pressure = values.component(0) - patchAverage
                + 9.81 * (installDepth_ - heightAxisValue);
            
            // Get cavitating cells
            forAll(weights,i)
            {
                if
                (
                    pressure[i] < -97.7 
                 && sign * distanceValue[i] < sign * maxDistance
                )
                {
                    weights[i] = pTraits<Type>::one;
                }
            }
            // Calc cavitation volume
            result = gSum(V*weights);
        }
    }

    return result;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::functionObjects::fieldValues::cavitationVolume::writeValues
(
    const word& fieldName,
    const scalarField& V
)
{
    const bool ok = validField<Type>(fieldName);

    if (ok)
    {
        Field<Type> values(getFieldValues<Type>(fieldName));

        if (writeFields_)
        {
            word outName = fieldName + '_' + regionTypeNames_[regionType_];
            if (this->volRegion::regionName_ != polyMesh::defaultRegion)
            {
                outName = outName + '-' + this->volRegion::regionName_;
            }

            IOField<Type>
            (
                IOobject
                (
                    outName,
                    obr_.time().timeName(),
                    obr_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                values
            ).write();
        }

        if (operation_ != opNone)
        {

            Type result = processValues(values, V);


            // Write state/results information
            word prefix, suffix;
            {
                prefix += operationTypeNames_[operation_];
                prefix += '(';
                suffix += ')';
            }

            word regionPrefix;
            if (this->volRegion::regionName_ != polyMesh::defaultRegion)
            {
                regionPrefix = this->volRegion::regionName_ + ',';
            }

            word resultName = prefix + regionPrefix + fieldName + suffix;

            Log << "    " << prefix << this->volRegion::regionName_ << suffix
                << " = ";

            file()<< tab << result;

            Log << result << endl;

            this->setResult(resultName, result);
        }
    }

    return ok;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::functionObjects::fieldValues::cavitationVolume::filterField
(
    const Field<Type>& field
) const
{
    if (this->volRegion::useAllCells())
    {
        return field;
    }

    return tmp<Field<Type>>::New(field, cellIDs());
}


// ************************************************************************* //
