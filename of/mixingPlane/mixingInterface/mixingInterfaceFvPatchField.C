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

\*---------------------------------------------------------------------------*/

#include "mixingInterfaceFvPatchField.H"
#include "mixingInterfacePolyPatch.H"

#include "volFields.H"
#include "surfaceFields.H"
#include "interpolationCell.H"
#include "cylindricalCS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
    
template<class Type>
void mixingInterfaceFvPatchField<Type>::readMixingType() const
{
    const dictionary& dict =
        this->patch().boundaryMesh().mesh().schemesDict().subDict
        (
            "mixingInterface"
        );
    
    word fieldName = this->internalField().name();
 
    if (dict.found(fieldName))
    {
        mixing_ = dict.getWord(fieldName);

    }
    
    else if (dict.found("default"))
    {
        mixing_ = dict.getWord("default");
    }
    
    else
    {
        FatalIOErrorIn
        (
            "void mixingInterfaceFvPatchField<Type>"
            "::readMixingType() const",
            dict
        )   << "Cannot find mixing type for field "
            <<  this->internalField().name() << nl
            << "Please specify in fvSchemes in mixingPlane, "
            << "under field name " << nl
            << "Available types are:" << nl
            << "areaAveraging" << nl << "consistentAveraging" << nl
            << "fluxAveraging"
            << abort(FatalIOError);
    }
    
    if (mixing_ == "consistentAveraging")
    {
        if 
        (
            this->internalField().name() != "p"
         && this->internalField().name() != "U"
        )
        {
            FatalIOErrorIn
            (   
                "void mixingInterfaceFvPatchField<Type>"
                "::readMixingType() const",
                dict
            )   << "Disallowed mixing type consistentAveraging for field "
                << this->internalField().name() << nl
                << "consistentAveraging is only a valid mixing type for p and U"
                << abort(FatalIOError);
        }
        else if
        (
            (
                this->internalField().name() == "p" 
             && dict.getWord("U") != "consistentAveraging"  
            )
         ||
            (
                this->internalField().name() == "U" 
             && dict.getWord("p") != "consistentAveraging"   
            )
        )
        {
            FatalIOErrorIn
            (   
                "void mixingInterfaceFvPatchField<Type>"
                "::readMixingType() const",
                dict
            )   << "Disallowed mixing type consistentAveraging for field "
                << this->internalField().name() << nl
                << "consistentAveraging has to be the mixing type"
                << " for both p and U"
                << abort(FatalIOError);   
        }
    }

    
}


template<class Type>
const scalarField& mixingInterfaceFvPatchField<Type>::fluxMask() const
{
    if (!this->db().objectRegistry::found("phi"))
    {
        InfoIn
            (
                "void mixingInterfaceFvPatchField<Type>::calcFluxMask()"
                ""
            )   << "Flux not found for flux averaging mixing plane on "
                << this->patch().name() << " for field "
                << this->internalField().name()
                << ".  Flux field: phi" 
                << endl;

            fluxMask_.setSize(this->size(), 0);
    }
    else
    {  
        const scalarField& phip = this->patch().lookupPatchField
        (
            "phi",
            reinterpret_cast<const surfaceScalarField*>(0),
            reinterpret_cast<const scalar*>(0)
        );
        
        // Use flux mask to decide wheter to use zeroGradient or fixedValue
        // If flux mask = 1 use interpolated value, otherwise use zeroGradient
        if (this->internalField().name() == "p")
        {
            fluxMask_ =  pos(mixingPlanePatch_.circumferentialAverage(phip));
        }
        else
        {
            fluxMask_ =  neg(mixingPlanePatch_.circumferentialAverage(phip));
        }
    }

    return fluxMask_;
}


template<class Type>
const scalarField& mixingInterfaceFvPatchField<Type>::fluxWeights() const
{
    // Calculate flux weights needed for flux averaging. 
    // Normal velocity is needed for flux weights, not phi
    scalar maxUNormal = 0;

    const scalarField& uNormal = calcNormalU();

    if (!uNormal.empty())
    {
        maxUNormal = max(uNormal);
    }

    reduce(maxUNormal, maxOp<scalar>());

    if (maxUNormal > SMALL)
    {
        fluxWeights_ = Foam::max(uNormal, scalar(0));
        fluxWeights_ /=
            mixingPlanePatch_.circumferentialAverage(fluxWeights_) + SMALL;
   }
    else
    {
        fluxWeights_ = 1;
    }

    return fluxWeights_;
}


template<class Type>
const scalarField mixingInterfaceFvPatchField<Type>::calcNormalU() const
{
    const fvPatchField<vector>& patchField =
        mixingPlanePatch_.lookupPatchField<volVectorField, vector>("U");

    tmp<Field<vector> > patchInternalFieldU = patchField.patchInternalField();

    const vectorField patchFieldU = patchInternalFieldU;

    const scalarField uNormal = patchFieldU & mixingPlanePatch_.nf();
    
    return uNormal;
}


template<class Type>
tmp<Field<Type> > mixingInterfaceFvPatchField<Type>::consistentAveU() const
{ 
    // Radial, tangential, normal dircetion change depending on rotation axis
    // 1. rotation axis (1,0,0) -> (normal, radial, tangential)
    // 2. rotation axis (0,1,0) -> (radial, normal, tangential)
    // 3. rotation axis (0,0,1) -> (radial, tangential, normal)

    label radial = 0, tangential = 1, normal = 2;
    
    if (mixingPlanePatch_.axis().component(0) == 1)
    {
        radial = 1, tangential = 2, normal = 0;   
    }
    else if (mixingPlanePatch_.axis().component(1) == 1)
    {
        radial = 0, tangential = 2, normal = 1;
    }
    
    if (not transformTensorsUpdated_)
    {
        // Transformation of the velocity vector is needed in such a way, that  
        // we have a velocity component normal to the interface. We go from the 
        // shadow patch to this patch, so use shadow patch side to build
        // forward transformation tensor and this patch side to build backward
        // transformation tensor. Angle for transformation is defined by the 
        // angle between rotation axis and normal vector
        
        // Get rotation axis and change sign if negative
        const vector& axis = cmptMultiply
            (
                mixingPlanePatch_.axis(), mixingPlanePatch_.axis()
            );
        
        // Take the opposite sign for the rotation axis for the shadow side
        // since the sign of the normal vector is reversed. Which direction we 
        // take does not matter, since this only changes the rotation angle by 
        // 180 degrees. This only flips the sign of the velocity vector, which 
        // makes no difference for the averaging. The backward transformation
        // will correct the sign of the velocity vector in this case.
        vectorField shRotationAxis(this->shadowPatchField().size());
        shRotationAxis = -1 * axis;
        
        vectorField rotationAxis(this->patch().size());
        rotationAxis = axis;

        // Get normal vectors 
        // Use no reference here, otherwise you will get garbage
        const vectorField shNormalVec = this->shadowPatchField().patch().nf();
        const vectorField normalVec = mixingPlanePatch_.nf();
        
        // Transform normal vectors into cylindrical coordinate system and 
        // eliminate phi direction
        vectorField shNormalVecCyl = 
            mixingPlanePatch_.shadowPatch().transformToCylindrical(shNormalVec);
        shNormalVecCyl.replace(tangential,0);

        vectorField normalVecCyl = 
                mixingPlanePatch_.transformToCylindrical(normalVec);
        normalVecCyl.replace(tangential,0);
        
        // The axis around the velocity vector rotates is the phi axis
        vector transformAxis(vector::zero);    
        transformAxis.replace(tangential,1);
        
        // Build the transformation tensors
        forwardTransform_ = 
            mixingPlanePatch_.calcRotationTensor
                (
                    transformAxis, shNormalVecCyl, shRotationAxis
                );

        backwardTransform_ = mixingPlanePatch_.calcRotationTensor
            (
                transformAxis, rotationAxis, normalVecCyl
            );
        
        transformTensorsUpdated_ = true;
    }
    
    tmp<Field<Type> > trField
    (
        new Field<Type> (this->patch().size(), pTraits<Type>::zero)
    );
    
    Field<Type>& rField = trField.ref();
    
    // Get shadow patch field
    Field<Type> shField = this->shadowPatchField().patchInternalField();
    
    // Get radius of faces needed for averaging of the tangential velocity
    const vectorField& shCellCenter = this->shadowPatchField().patch().Cf();
    
    const vectorField shCellCenterCyl = 
        mixingPlanePatch_.shadowPatch().transformToCylindrical(shCellCenter);
    
    const scalarField radius = shCellCenterCyl.component(radial);
    
    //Build an average radius needed for averaging of the tangential veloctiy
    const scalarField radiusAve = radius / 
        (
            mixingPlanePatch_.shadowPatch().circumferentialAverage(radius)
            + SMALL
        );
    
    // Transform the velocity Field
    Field<Type> shFieldCyl = 
        mixingPlanePatch_.shadowPatch().transformToCylindrical(shField);
    
    Field<Type> shFieldLoc = transform(forwardTransform_, shFieldCyl);
    
    // Averaging of the velocity vector
    const scalarField uNormal = shFieldLoc.component(normal);
    const scalarField uNormalAveraged = mixingPlanePatch_.interpolate(uNormal);
    
    // Use flux averaging for the radial component
    const scalarField& shFluxWeights = this->shadowPatchField().fluxWeights();
    
    const scalarField uWRadial = shFieldLoc.component(radial) * shFluxWeights;
    
    const scalarField uRadialAveraged = mixingPlanePatch_.interpolate(uWRadial);
    
    const scalarField uWTangential = shFieldLoc.component(normal) * 
        shFieldLoc.component(tangential) * radiusAve;
    
    const scalarField uTangentialAveraged = mixingPlanePatch_.interpolate
        (
            uWTangential
        ) / uNormalAveraged;
    
    Field<Type> hsFieldLoc(this->patch().size());
    hsFieldLoc.replace(radial,uRadialAveraged);
    hsFieldLoc.replace(tangential,uTangentialAveraged);
    hsFieldLoc.replace(normal,uNormalAveraged);
    
    // Transform the velocity vector back into cartesian coordinates
    Field<Type> uAveragedCyl = transform(backwardTransform_, hsFieldLoc);
    
    rField = mixingPlanePatch_.transformToCartesian(uAveragedCyl);
    
    return trField;
}


template<class Type>
tmp<Field<Type> > mixingInterfaceFvPatchField<Type>::consistentAveP() const
{
    tmp<Field<Type> > trField
    (
        new Field<Type> (this->patch().size(), pTraits<Type>::zero)
    );
    
    Field<Type>& rField = trField.ref();
    
    // Get shadow pressure field
    Field<Type> shField = this->shadowPatchField().patchInternalField();
    
    // Normal velocity is needed for the correction of the pressure
    const scalarField& uNormal = calcNormalU();
    
    const scalarField uNormalPow = pow(uNormal, 2);
    
    const scalarField uNormalAveraged =
        mixingPlanePatch_.circumferentialAverage(uNormal);
    
    const scalarField pressure = shField.component(0);
    
    // Average the pressure
    scalarField pAveraged = mixingPlanePatch_.interpolate(pressure)
        - mixingPlanePatch_.circumferentialAverage(uNormalPow)
        + pow(uNormalAveraged, 2);
    
    rField.replace(0, pAveraged);
    
    return trField;
}


template<class Type>
tmp<Field<Type> > mixingInterfaceFvPatchField<Type>::patchValues() const
{
    // Read mixing type
    readMixingType();
    
    tmp<Field<Type> > ptrField
    (
        new Field<Type> (this->size(), pTraits<Type>::zero)
    );
    
    Field<Type>& prField = ptrField.ref();
    
    Field<Type> shField = this->shadowPatchField().patchInternalField();
    
    const scalarField& mask = fluxMask();
    
    if (mixing_ == "areaAveraging")
    {
        prField = mask * mixingPlanePatch_.interpolate(shField)
            + (1-mask) * this->patchInternalField();
    }
    
    else if (mixing_ == "fluxAveraging")
    {
        const scalarField& shadowFluxWeights = shadowPatchField().fluxWeights();
        
       Field<Type> wShField = shadowFluxWeights.component(0) * shField;
        
        prField = mask * mixingPlanePatch_.interpolate(wShField)
            + (1-mask) * this->patchInternalField();
    }
    
    else if (mixing_ == "consistentAveraging")
    {
        if (this->internalField().name() == "p")
        {
            Field<Type> pAveraged = this->consistentAveP();
            
            prField = mask * pAveraged + (1-mask) * this->patchInternalField();
        }
        else
        {
            Field<Type> uAveraged = this->consistentAveU();
            
            prField = mask * uAveraged + (1-mask) * this->patchInternalField();
        }
    }
    
    else
    {
       FatalErrorIn
        (
            "tmp<Field<Type> > mixingInterfaceFvPatchField<Type>::"
            "patchValues() const"
        )   << "Unknown mixing type for patch " << this->patch().name()
            << " for field "
            <<  this->internalField().name()
            << abort(FatalError); 
    }

    
    return ptrField;
}    
  
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class Type>
mixingInterfaceFvPatchField<Type>::mixingInterfaceFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    mixingPlanePatch_(refCast<const mixingInterfaceFvPatch>(p)),
    mixing_("mixing_unknown"),
    fluxMask_(),
    fluxWeights_(p.size(), 0),
    forwardTransform_(),
    backwardTransform_(),
    transformTensorsUpdated_(false)
        
{}
    

template <class Type>
mixingInterfaceFvPatchField<Type>::mixingInterfaceFvPatchField
(
    const mixingInterfaceFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    mixingPlanePatch_(refCast<const mixingInterfaceFvPatch>(p)),
    mixing_(ptf.mixing_),
    fluxMask_(),
    fluxWeights_(ptf.fluxWeights_, mapper),
    forwardTransform_(ptf.forwardTransform_),
    backwardTransform_(ptf.backwardTransform_),
    transformTensorsUpdated_(ptf.transformTensorsUpdated_)
{}
    

template <class Type>
mixingInterfaceFvPatchField<Type>::mixingInterfaceFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF, dict),
    mixingPlanePatch_(refCast<const mixingInterfaceFvPatch>(p)),
    mixing_("mixing_unknown"),
    fluxMask_(),
    fluxWeights_(p.size(), 0),
    forwardTransform_(),
    backwardTransform_(),
    transformTensorsUpdated_(false)
{}
    

template <class Type>
mixingInterfaceFvPatchField<Type>::mixingInterfaceFvPatchField
(
    const mixingInterfaceFvPatchField& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    mixingPlanePatch_(refCast<const mixingInterfaceFvPatch>(ptf.patch())),
    mixing_(ptf.mixing_),
    fluxMask_(),
    fluxWeights_(ptf.fluxWeights_),
    forwardTransform_(ptf.forwardTransform_),
    backwardTransform_(ptf.backwardTransform_),
    transformTensorsUpdated_(ptf.transformTensorsUpdated_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const mixingInterfaceFvPatchField<Type>&
mixingInterfaceFvPatchField<Type>::shadowPatchField() const
{
    const GeometricField<Type, fvPatchField, volMesh>& fld =
        static_cast<const GeometricField<Type, fvPatchField, volMesh>&>
        (
            this->internalField()
        );

    return refCast<const mixingInterfaceFvPatchField<Type> >
    (
        fld.boundaryField()[mixingPlanePatch_.shadowPatchID()]
    );
}


template<class Type>
void mixingInterfaceFvPatchField<Type>::updateCoeffs()
{   
    // Read mixing type
    readMixingType();

    fixedValueFvPatchField<Type>::updateCoeffs();  
}


template<class Type>
void mixingInterfaceFvPatchField<Type>::initEvaluate
(
    const Pstream::commsTypes commsType
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }
    
    // Get values on the patch
    Field<Type> patchFieldValues = this->patchValues();
    
    Field<Type>::operator=(patchFieldValues);
}


template<class Type>
void mixingInterfaceFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes
)
{
    fixedValueFvPatchField<Type>::evaluate();
}


template <class Type>
void mixingInterfaceFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
