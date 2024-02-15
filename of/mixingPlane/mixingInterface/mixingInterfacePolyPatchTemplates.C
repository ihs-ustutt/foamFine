/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 CFD support s.r.o. 
     \\/     M anipulation  |
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

Note
    Developed by IHS, University of Stuttgart, in 2021.
 
 \*---------------------------------------------------------------------------*/


#include "mixingInterfacePolyPatch.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::mixingInterfacePolyPatch::toProfile
(
    const Field<Type>& pf
) const
{
    if (not divided_)
    {
        patchDivision();
    }
    
    tmp<Field<Type> > trField
    (
        new Field<Type> (numPlanes(), pTraits<Type>::zero)
    );
    
    Field<Type>& rField = trField.ref();
    
    // Transform into cylinidrical coordinate system
    Field<Type> pfCyl = transformToCylindrical(pf);
    
    const Field<scalar>& area = this->magFaceAreas();
   
    for (label iplane = 0; iplane < numPlanes(); iplane++)
    {
        List<label> const & faces = mxpFaces_[iplane];
        List<scalar> const & weights = mxpWeights_[iplane];
        
        scalar loc_wgt = 0;
        
        forAll (faces, i)
        {
            loc_wgt += area[faces[i]] * weights[i];
        }
        
        reduce(loc_wgt, sumOp<scalar>());

        forAll (faces, i)
        {
            rField[iplane] += 
                (area[faces[i]] * weights[i] * pfCyl[faces[i]])/loc_wgt;
        }      
    }
    
    reduce(rField, sumOp<Field<Type> >()); 
        
    return trField;
}

template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::mixingInterfacePolyPatch::circumferentialAverage
(
    const Field<Type>& tpf
)const
{
    if (not divided_)
    {
        patchDivision();
    }
    
    tmp<Field<Type> > trField
    (
        new Field<Type> (tpf.size(), pTraits<Type>::zero)
    );
    
    Field<Type>& rField = trField.ref();
    
    // Interpolate from patch to ribbon
    Field<Type> ribbonField = this->toProfile(tpf);
    
    // Interpolate from ribbon to patch
    for (label iplane = 0; iplane < numPlanes(); iplane++)
    {
        List<scalar> const & weights = mxpWeights_[iplane];
        
        List<label> const & faces = mxpFaces_[iplane];
        
        forAll(faces, i)
        {
        
            rField[faces[i]] += ribbonField[iplane] * weights[i]; 
        }      
    }
    // Transform back to cartesion coordinates
    rField = transformToCartesian(rField);
    
    return trField;
}

template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::mixingInterfacePolyPatch::interpolate
(
    const Field<Type>& spf
)const
{
    if (not divided_)
    {
        patchDivision();
    }

    tmp<Field<Type> > trField
    (
            new Field<Type> (this->size(), pTraits<Type>::zero)
    );
    
    Field<Type>& rField = trField.ref();
    
    // Interpolate from shadow patch to ribbon
    Field<Type> ribbonField = 
        shadowPatch().toProfile(spf);
    
    // Interpolate from ribbon to patch
    for (label iplane = 0; iplane < numPlanes(); iplane++)
    {       
        List<scalar> const & weights = mxpWeights_[iplane];
        
        List<label> const & faces = mxpFaces_[iplane];
        
        forAll(faces, i)
        {
        
            rField[faces[i]] += ribbonField[iplane] * weights[i]; 
        }
    }
    
    // Transform back to cartesion coordinates
    rField = transformToCartesian(rField);
    
    return trField;
}


template<class Type>
Foam::tmp<Foam::Field<Type> > 
Foam::mixingInterfacePolyPatch::transformToCylindrical
(
    const Field<Type>& tpf
) const
{
    if (not rotationTensorsUpdated_)
    {
        updateRotationTensors();
    }
    
    return transform(cartesianToCylindrical_, tpf);
}


template<class Type>
Foam::tmp<Foam::Field<Type> > 
Foam::mixingInterfacePolyPatch::transformToCartesian
(
    const Field<Type>& tpf
) const
{
    if (not rotationTensorsUpdated_)
    {
        updateRotationTensors();
    }
    
    return transform(cylindricalToCartesian_, tpf);

}




