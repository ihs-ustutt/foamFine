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
    Foam::mixingInterfaceFvPatchField

Group
    grpGenericBoundaryConditions grpCoupledBoundaryConditions

Description
    Patch for the mixingInterface boundary condition based on the theory of an
    explicit mixing plane interface.
    
    Example usage:
    
        rotor_stator_interface
        {
            type    mixingInterface;
            value   $internalField;
        }
    
    Includes the averaging methods area averaging, flux averaging and the
    consistent averaging method, which provides the conservation of the mass
    flow as well as the conservation of momentum and angular momentum. The
    averaging method must be set in system/fvSchmes, for example:
    
    mixingInterface
    {
        U       consistentAveraging;
        p       consistentAveraging;
        epsilon fluxAveraging;
        k       fluxAveraging;
    }
    
Note
    Developed by IHS, University of Stuttgart, in 2021.

SourceFiles
    mixingInterfaceFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef mixingInterfaceFvPatchField_H
#define mixingInterfaceFvPatchField_H

#include "fixedValueFvPatchFields.H"
#include "interpolation.H"
#include "mixingInterfaceFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    
/*---------------------------------------------------------------------------*\
                  Class mixingInterfaceFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template <class Type>
class mixingInterfaceFvPatchField
:
    public fixedValueFvPatchField<Type>
{
    
    // Private data
    
        const mixingInterfaceFvPatch& mixingPlanePatch_;

        mutable word mixing_;

        mutable scalarField fluxMask_;

        mutable scalarField fluxWeights_;

    // Private member functions

        //- Read mixing Type
        void readMixingType() const;
        
        //- Return flux mask.  
        const scalarField& fluxMask() const;

        //- Return flux weights.  
        const scalarField& fluxWeights() const;

        //- Return velocity normal to the patch
        const scalarField calcNormalU() const;
        
        //- Return velocity vector consistent averaged
        tmp<Field<Type> > consistentAveU() const;
        
        //- Return pressure consistent averaged
        tmp<Field<Type> > consistentAveP() const;
        
        //- Return the values on the patch
        tmp<Field<Type> > patchValues() const;

protected:
    
    //- Transformation tensors
    
    mutable tensorField forwardTransform_;
        
    mutable tensorField backwardTransform_;
        
    mutable bool transformTensorsUpdated_;
    
public:

    //- Runtime type information
    TypeName("mixingInterface");


    // Constructors

        //- Construct from patch and internal field
        mixingInterfaceFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        mixingInterfaceFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  mixingInterfaceFvPatchField
        //  onto a new patch
        mixingInterfaceFvPatchField
        (
            const mixingInterfaceFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        mixingInterfaceFvPatchField( const mixingInterfaceFvPatchField&) = delete;


        //- Construct as copy setting internal field reference
        mixingInterfaceFvPatchField
        (
            const mixingInterfaceFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type> > clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type> >
            (
                new mixingInterfaceFvPatchField<Type>
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions
        
        // Access
        
            //- Return shadow patch field
            const mixingInterfaceFvPatchField<Type>&
            shadowPatchField() const;
        
        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();
            
            //- Initialise the evaluation of the patch field
            virtual void initEvaluate(const Pstream::commsTypes commsType);

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType
            );
            

        //- Write
            
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "mixingInterfaceFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
