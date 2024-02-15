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

Class
    Foam::mixingInterfaceFvPatch

Description
    Foam::mixingInterfaceFvPatch

SourceFiles
    mixingInterfaceFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef mixingInterfaceFvPatch_H
#define mixingInterfaceFvPatch_H

#include "fvPatch.H"
#include "mixingInterfacePolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class mixingInterfaceFvPatch Declaration
\*---------------------------------------------------------------------------*/

class mixingInterfaceFvPatch
:
    public fvPatch
{
    // Private Data
    
    //- Reference to polyPatch
    const mixingInterfacePolyPatch& mixingPlanePolyPatch_;

public:

    //- Runtime type information
    TypeName(mixingInterfacePolyPatch::typeName_());


    // Constructors

        //- Construct from components
        mixingInterfaceFvPatch(const polyPatch& patch, const fvBoundaryMesh& bm)
        :
            fvPatch(patch, bm),
            mixingPlanePolyPatch_(refCast<const mixingInterfacePolyPatch>(patch))
        {}
        
    // Member functions
        
        // Access
        
            //- Return shadow patch name
            word shadowPatchName() const
            {
                return mixingPlanePolyPatch_.shadowPatchName();
            }
            
            //- Return shadow patch id
            label shadowPatchID() const
            {
                return mixingPlanePolyPatch_.shadowPatchID();
            }
            
            //- Return rotation axis
            vector axis() const
            {
                return mixingPlanePolyPatch_.axis();
            }
            
            //- Return shadow patch
            const mixingInterfacePolyPatch& shadowPatch() const
            {
                return mixingPlanePolyPatch_.shadowPatch();
            };
            
            //- Return number of profile bands
            label numPlanes() const
            {
                return mixingPlanePolyPatch_.numPlanes();
            }
                        
        // Interpolation
        
            template<class Type>
            tmp<Field<Type> > toProfile(const Field<Type>& pf) const
            {
                return 
                    mixingPlanePolyPatch_.toProfile(pf);
            }
            
            template<class Type>
            tmp<Field<Type> > circumferentialAverage
            (
                const Field<Type>& tpf
            ) const
            {
                return 
                    mixingPlanePolyPatch_.circumferentialAverage(tpf);
            }
            
            template<class Type>
            tmp<Field<Type> > interpolate(const Field<Type>& spf) const
            {
                return 
                    mixingPlanePolyPatch_.interpolate(spf);
            }
            
        // Transformation
            
            template<class Type>
            tmp<Field<Type> > transformToCylindrical
            (   const Field<Type>& tpf
            ) const
            {
                return 
                    mixingPlanePolyPatch_.transformToCylindrical(tpf);
            }
            
        
            template< class Type>
            tmp<Field<Type> > transformToCartesian 
            (
                const Field<Type>& tpf
            ) const
            {
                return 
                    mixingPlanePolyPatch_.transformToCartesian(tpf);
            }
            
            tensorField calcRotationTensor
            (
                const vector& rotationAxis,
                const vectorField& vi,
                const vectorField& vf
            ) const
            {
                return 
                    mixingPlanePolyPatch_.calcRotationTensor
                    (
                        rotationAxis,
                        vi,
                        vf
                    );
            }
                
            
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
