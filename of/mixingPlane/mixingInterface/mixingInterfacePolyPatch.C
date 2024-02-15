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

#include "mixingInterfacePolyPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "polyBoundaryMesh.H"
#include "polyMesh.H"
#include "vectorList.H"
#include "OFstream.H"


const Foam::scalar Inf = 1E+30;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mixingInterfacePolyPatch, 0);
 
    addToRunTimeSelectionTable(polyPatch, mixingInterfacePolyPatch, word);
    addToRunTimeSelectionTable(polyPatch, mixingInterfacePolyPatch, dictionary);
}

const Foam::Enum<Foam::mixingInterfacePolyPatch::division>
    Foam::mixingInterfacePolyPatch::divisionNames_
    ({
        { division::UNIFORM, "uniform" },
        { division::MANUAL, "manual" },
        { division::MESH_DEPENDENT, "meshDependent" },
        { division::DIVISION_UNKNOWN, "divisionUnknown" },
    });

const Foam::Enum<Foam::mixingInterfacePolyPatch::stackAxis>
    Foam::mixingInterfacePolyPatch::stackAxisNames_
    ({
        { stackAxis::STACK_R, "R" },
        { stackAxis::STACK_Z, "Z" },
        { stackAxis::STACK_UNKNOWN, "stackUnknown" },
     });

const Foam::Enum<Foam::mixingInterfacePolyPatch::discretisation>
    Foam::mixingInterfacePolyPatch::discretisationNames_
    ({
        { discretisation::BOTH_PATCHES, "bothPatches" },
        { discretisation::THIS_PATCH, "thisPatch" },
        { discretisation::SHADOW_PATCH, "shadowPatch" },
        { discretisation::USER_DEFINED, "userDefined" },
        { discretisation::DISCRETISATION_UNKNOWN, "discretisationUnknown" },
    });


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::mixingInterfacePolyPatch::updateRotationTensors() const
{
    // Define radial-axis, use (1, 0, 0) except for the case, when the 
    // rotation axis is (1, 0, 0)
    vector radialAxis(1,0,0);
            
    if (axis().component(0) == 1)
    {
        radialAxis.replace(0, 0);
        radialAxis.replace(1, 1);
    }

    vectorField globalVectors = this->faceCentres() - origin();
     
    // Transform globalVectors in such a way, that phi is zero. The angle
    // between the transformed and non-transformed vector corresponds to the
    // angle phi needed to calculate the transformation tensor.
    vectorField transformedVectors(globalVectors.size(), vector::zero);
    
    forAll (globalVectors, i)
    {
        vector er = globalVectors[i]/mag(globalVectors[i]);
        
        vector en = axis();          en /= mag(en);
        vector et = axis() ^ er;   et /= mag(et);
        vector em = en ^ et;     em /= mag(em);
        
        // Positive values are needed for the radial component, otherwise
        // transformation normal to the interface in
        // mixingInterfaceFvPatchField<Type>::consistentAveU() is wrong
        transformedVectors[i] = 
            (
                (globalVectors[i] & en) * en + 
                mag((globalVectors[i] & em)) * radialAxis
            );
    }
    
    cartesianToCylindrical_ = 
            calcRotationTensor(axis(), globalVectors, transformedVectors);
    
    cylindricalToCartesian_ = 
            calcRotationTensor(axis(), transformedVectors, globalVectors);
    
    rotationTensorsUpdated_ = true;
}

void Foam::mixingInterfacePolyPatch::patchDivision() const
{
    for (label iplane = 0; iplane < numPlanes(); iplane++)
    {
        List<label> thisFaces;
        List<scalar> thisWeights;  

        planeFaces(iplane, thisFaces, thisWeights);

        mxpFaces_.append(thisFaces);
        mxpWeights_.append(thisWeights);
    }

    divided_ = true;
}

Foam::direction Foam::mixingInterfacePolyPatch::getStackDirection() const
{
    switch (stackAxis_)
    {
        case STACK_R:
        {
            return vector::X;
        }
        break;
    
        case STACK_Z:
        {
            return vector::Z;
        }
        break;

        default:
        {
            FatalErrorIn
            (
                "direction mixingInterfacePolyPatch::getStackDirection() const"
            )   << "Bad stackAxis type "  << stackAxisNames_[stackAxis_] 
                << " Available types: R, Z"
                << exit(FatalError);
            
            // Dummy return
            return vector::X;
        }
        break;
    }
    
}

Foam::scalar Foam::mixingInterfacePolyPatch::calcGradingDist
(
    const scalar& varG,
    const scalar& varX
) const
{
    return 0.5 + 0.5 * tanh(varG * (-1 + 2 * varX)) / tanh(varG);
}

Foam::pointField Foam::mixingInterfacePolyPatch::calcUserDefinedRibbons
(
    const pointField tempProfile
) const
{   
    // Current number of planes after discretisation bothPatches
    const label curNplanes = tempProfile.size() - 1;
    
    // Check if user definitions are appropriate
    if (2*nPlanesBl_ >= nPlanes_ || nPlanes_ > curNplanes)
    {
        FatalErrorIn
        (
            "void mixingInterfacePolyPatch::calcUserDefinedRibbons() const"
        )   << "Something is wrong with the number of Planes specified"
            << " for patch " << shadowPatchName() << nl 
            << "Check if the number of planes in the boundary layer is greater " 
            << "or equal than the total number of planes: " << nl 
            << 2*nPlanesBl_ << " >= " << nPlanes_ << nl << "Check if the " 
            << "number of planes is greater than that resulting from "
            << "discretisation bothPatches: " << nl << nPlanes_ << " > "
            << curNplanes
            << exit(FatalError);
    }
    
    // Total number of ribbons requested
    const label nRibbons = nPlanes_ + 1;
    
    // Number of Ribbons requested for the boundary layer
    const label nRibbonsBl = nPlanesBl_ + 1;
    
    // New ribbon profile
    pointField profile(nRibbons);
    
    // Number of ribbons between boundary layer
    const label nRibbonsIf = nRibbons - 2 * nRibbonsBl;
    
    // Calculate the distance of the ribbons between the boundary layer. In
    // case of a conical interface, the distance between the last ribbons of the
    // boundary layer is not the patch distance. Therefore, build the distance
    // as a sum of the points
    scalar distIf = 0;
    for (label i=nPlanesBl_; i< curNplanes-nPlanesBl_; i++)
    {
        distIf += mag(tempProfile[i+1] - tempProfile[i]);
    }
    
    // Build the boundary layer of the new ribbon profile
    for (label i=0; i<nRibbonsBl; i++)
    {
        profile[i] = tempProfile[i];
        profile[nPlanes_-i] = tempProfile[curNplanes-i];
    }
    
    // Calculate the distance of the ribbons between the boundary layer. 
    // If gradingIF_ = false use equidistant ribbon distribution
    scalarField distRibbons(nRibbonsIf, distIf/(nRibbonsIf+1));
    
    if (gradingIf_ && nRibbonsIf > 2)
    {
        // Calculate distance of the first internal ribbon depending on
        // the distance between last boundary layer ribbons 
        const scalar distFr = 0.5 * 1.3 * 
            (
                mag(tempProfile[nPlanesBl_] - tempProfile[nPlanesBl_-1]) 
                + mag(tempProfile[curNplanes-nPlanesBl_+1]
                - tempProfile[curNplanes-nPlanesBl_])
            );
        
        if (distFr >= distRibbons[0])
        {
            FatalErrorIn
            (
                "void mixingInterfacePolyPatch::calcUserDefinedRibbons() const"
            )   << "Calculated distance of the first internal ribbon with a "
                << "grading is >= distance with an equidistant ribbon "
                << "distribution" << nl << "Please set gradingIf to false or "
                << "decrease number of planes."
                << exit(FatalError);
        }

        const scalar dist = 1.0/(nRibbonsIf + 1);
        
        // Calculate varG according to distance of the first ribbon
        scalar varG = 10;
        
        // Make sure g is great enough
        if (calcGradingDist(varG, dist) * distIf - distFr > 0)
        {
            do
            {
                varG += 0.1;
            } while (calcGradingDist(varG, dist) * distIf - distFr > 0);
        }
        // Reduce g to get a good first guess
        do
        {
            varG -= 0.1;
        } while (calcGradingDist(varG, dist) * distIf - distFr < 0);
        
        // Increase g until tolerance is reached
        for (label i=0; i<1000; i++)
        {
            varG = varG + i * 1e-4;
            // Tolerance reached or difference between calculated distance and
            // needed distance for the first ribbon switched sign? In both cases
            // it makes no sense to continue the loop
            if (calcGradingDist(varG, dist) * distIf - distFr < distFr*1e-3) 
                break;
        }
        
        // Calculate the new distances of the ribbons between the boundary layer
        for (label i=0; i<round(nRibbonsIf/2.0); i++)
        {
            // First side
            distRibbons[i] = distIf *
                (
                    calcGradingDist(varG, (i+1)*dist) 
                    - calcGradingDist(varG, i*dist)
                );
            
            // Second side. If the number of ribbons between the boundary layer
            // is an odd number, the last ribbon is taken from the first side 
            if (i < nRibbonsIf/2)
            {
                distRibbons[nRibbonsIf-1-i] = distIf *
                (
                    calcGradingDist(varG, (i+1)*dist) 
                    - calcGradingDist(varG, i*dist)
                );
            }
        }
        
        if (debug)
        {
            Info << "Calculated grading for internal ribbon field" << nl
            << "\t" << "distFirstRibbon :" << distFr << nl
            << "\t" << "finalDistFirstRibbon :" << distRibbons[0] << nl
            << "\t" 
            << "Percentage deviation distFirstRibbon, finalDistFirstRibbon :"
            << distFr/distRibbons[0] - 1 << nl
            << "\t" << "calcRibbonDistances :" << distRibbons << endl;
        }  
    }
    
    // Build the ribbon patch between the boundary layer. Build it from both 
    // sides to get a smoother distribution
    label indexFirstSide = nRibbonsBl;
    label indexSecondSide = curNplanes - nRibbonsBl;
    
    for (label i=0; i<round(nRibbonsIf/2.0); i++)
    {
        // First side
        scalar previousDevFirstSide = GREAT;

        // Iterate over all ribbons resulting from discretisation bothPatches  
        // between the boundary layer and current second index to find the next 
        // ribbon with the requested distance
        for (label k=indexFirstSide; k<indexSecondSide; k++)
        {
            scalar ribbonDistFirstSide = 
                mag(profile[i+nPlanesBl_] - tempProfile[k]);

            scalar deviationFirstSide = 
                mag(ribbonDistFirstSide - distRibbons[i]);
            
            if (deviationFirstSide < previousDevFirstSide)
            {
                indexFirstSide = k;
            }
            
            previousDevFirstSide = deviationFirstSide;
        }
        
        profile[i+nRibbonsBl] = tempProfile[indexFirstSide];
        indexFirstSide += 1;
        
        // Do the same for the second side. If the number of ribbons between 
        // the boundary layer is an odd number, the last ribbon is taken from  
        // the first side. In this case, we don´t have to go in the loop
        // for the last ribbon.
        if (i < nRibbonsIf/2)
        {
            scalar previousDevSecondSide = GREAT;
            
            for (label k=indexSecondSide; k>indexFirstSide; k--)
            {
                scalar ribbonDistSecondSide = 
                    mag(profile[nPlanes_-nPlanesBl_-i] - tempProfile[k]);

                scalar deviationSecondSide = 
                    mag(ribbonDistSecondSide - distRibbons[nRibbonsIf-1-i]);

                if (deviationSecondSide < previousDevSecondSide)
                {
                    indexSecondSide = k;
                }

                previousDevSecondSide = deviationSecondSide;
            }

            profile[nPlanes_-nRibbonsBl-i] = tempProfile[indexSecondSide];
            indexSecondSide -= 1;
        }
    }
    
    return profile;
}

Foam::pointField Foam::mixingInterfacePolyPatch::computeProfileFromHistograms
(
    const profileHistogram& patchHisto,
    const profileHistogram& shadowHisto,
    const scalar halfSizeBin
) const
{
    // Find min, max bounds
    scalar histoMinValue =
        Foam::min
        (
            patchHisto.begin()->first,
            shadowHisto.begin()->first
        );

    scalar histoMaxValue =
        Foam::max
        (
            (--patchHisto.end())->first,
            (--shadowHisto.end())->first
        );
    
    // Next, we compare both histograms, leap-frogging from one histogram
    // to the other, everytime jumping to the next largest value from
    // a given position
    
    scalar curRvalue = histoMinValue;
    
    List<scalar> sep;
    sep.append(curRvalue);
    
    do
    {
        profileHistogram::const_iterator nextPatch =
            patchHisto.lower_bound(curRvalue + halfSizeBin - SMALL);

        profileHistogram::const_iterator nextShadow  =
            shadowHisto.lower_bound(curRvalue + halfSizeBin - SMALL);

        if
        (
            nextPatch == patchHisto.end()
         || nextShadow  == shadowHisto.end()
        )
        {
            // We are done
            if (curRvalue != histoMaxValue)
            {
                sep.append(histoMaxValue);
            }
            break;
        }

        // Leap frog to the next largest delta

        curRvalue = max(nextPatch->first, nextShadow->first);
        
        sep.append(curRvalue);
        
    } while (curRvalue != histoMaxValue);
    
    std::list<point> leapFrogProfile(0);
        
    forAll(sep, sI)
    {
        scalar curRvalue = sep[sI];
        if (patchHisto.find(curRvalue) != patchHisto.end())
        {
            leapFrogProfile.push_back
            (
                *(patchHisto.find(curRvalue)->second.begin())
            );
        }
        else
        {
            leapFrogProfile.push_back
            (
                *(shadowHisto.find(curRvalue)->second.begin())
            );

        }    
    }

    pointField profile(sep.size());

    label pI = 0;

    forAllIter (std::list<point>, leapFrogProfile, lI)
    {
        profile[pI++] = *lI;
    }
    
    return profile;    
}

void Foam::mixingInterfacePolyPatch::updateProfileHistogram
(
    profileHistogram& histo,
    const point& profileCoord,  
    const direction dir,        
    const scalar halfSizeBin    
) const
{
    bool foundNewBin = true;

    scalar keyValue = profileCoord.component(dir);

    forAllIter (profileHistogram, histo, histoI)
    {
        if
        (
            keyValue >= histoI->first - halfSizeBin
         && keyValue < histoI->first + halfSizeBin
        )
        {
            foundNewBin = false;
            histoI->second.push_back(profileCoord);
            break;
        }
    }

    if (foundNewBin)
    {
        std::list<point> initValue;

        initValue.push_back(profileCoord);
        histo.insert
        (
            std::pair<scalar, std::list<point> >(keyValue, initValue)
        );
    }
}

Foam::pointField Foam::mixingInterfacePolyPatch::calcRibbons() const
{
    List<scalar> sep;
    
    const direction stackAxis = getStackDirection();
    
    switch(patchDivision_)
    {
        case UNIFORM:
        {
            // Patch points.
            vectorField const & patchPoints = this->localPoints();
            
            // get point closest to and furthest from the axis
            scalar rmin = +Inf, rmax = -Inf;
            forAll (patchPoints, I)
            {
                scalar distance = 
                    mag((patchPoints[I] - origin_) ^ axis_);
                if (stackAxis_ == STACK_Z)
                {
                    distance = (patchPoints[I] - origin_) & axis_;
                }
                rmin = min(distance, rmin);
                rmax = max(distance, rmax);
            }

            // synchronize across processes
            reduce(rmin, minOp<scalar>());
            reduce(rmax, maxOp<scalar>());

            // calculate radial separation distances
            for (int i = 0; i <= nPlanes_; i++)
            {
                // replace lower bound by zero in case stackAxis = R
                if (i == 0 && stackAxis_ == STACK_R) sep.append(0.);
                // replace lower bound by -Infinity in case stackAxis = Z
                else if (i==0 && stackAxis_ == STACK_Z) sep.append(-Inf);
                // replace upper bound by Infinity
                else if (i == nPlanes_) sep.append(Inf);
                // others keep as requested
                else sep.append(rmin + i * (rmax - rmin) / nPlanes_);
            }            
            pointField profile(sep.size());
            
            profile.replace(stackAxis, sep);
            
            return profile;
        }
        break;

        case MANUAL:
        {   
            sep.append(orgsep_);
            
            pointField profile(sep.size());
            
            profile.replace(stackAxis, sep);
            
            // replace lower bound by zero in case stackAxis = R
            if (stackAxis_ == STACK_R)
                profile[0].replace(stackAxis, 0);
            // replace lower bound by -Infinity in case stackAxis = Z
            else if (stackAxis_ == STACK_Z)
                profile[0].replace(stackAxis, -Inf);
            // replace upper bound by Infinity
            profile[profile.size() - 1].replace(stackAxis, Inf);
            
            return profile;
        }
        break;

        case MESH_DEPENDENT:
        {
            // Define a local cylindrical coordinate system and transform mesh
            // points to the local coordinate system
            
            vector dirn = vector(1,0,0);
            if (axis().component(0) == 1) dirn = vector(0,1,0);
                    
            autoPtr<coordinateSystem> cs_ (
                new coordSystem::cylindrical
                (
                    "mixingCS",
                    vector::zero,
                    axis(),
                    dirn
                )
            );

            pointField patchGlobalProfile = 
                cs_().localPosition(this->localPoints());

            pointField shadowGlobalProfile = 
                cs_().localPosition(shadowPatch().localPoints());
            
            // Set sweep axis to zero. I don´t think we need a definition of a 
            // sweep axis, because averaging is always done in circumferential
            // direction. 
            patchGlobalProfile.replace(1,0);
            shadowGlobalProfile.replace(1,0);
            
            // Find minimum edge length for this patch and shadow patch
            scalar patchMinEdgeLength = GREAT;
            const edgeList& patchEdgeList = this->edges();

            forAll (patchEdgeList, pEi)
            {
                patchMinEdgeLength = 
                    Foam::min
                    (
                        patchMinEdgeLength,
                        patchEdgeList[pEi].mag(this->localPoints())
                    );
            }

            scalar shadowMinEdgeLength = GREAT;
            const edgeList& shadowEdgeList = shadowPatch().edges();

            forAll (shadowEdgeList, sEi)
            {
                shadowMinEdgeLength = 
                    Foam::min
                    (
                        shadowMinEdgeLength,
                        shadowEdgeList[sEi].mag(shadowPatch().localPoints())
                    );     
            }

            scalar halfMinSizeBin = 
                    Foam::max(patchMinEdgeLength, shadowMinEdgeLength)/2.0;

            // Build the histogramms
            profileHistogram patchHistogram;
            profileHistogram shadowHistogram;

            forAll (patchGlobalProfile, pI)
            {
                updateProfileHistogram
                (
                    patchHistogram,
                    patchGlobalProfile[pI],
                    stackAxis,
                    halfMinSizeBin
                );
            }

            forAll (shadowGlobalProfile, sI)
            {
                updateProfileHistogram
                (
                    shadowHistogram,
                    shadowGlobalProfile[sI],
                    stackAxis,
                    halfMinSizeBin   
                );
            }
            
            // Ribbon profile
            pointField profile;

            switch (discretisation_)
            {
                case BOTH_PATCHES:
                {
                    profile = computeProfileFromHistograms
                    (
                        patchHistogram,
                        shadowHistogram,
                        halfMinSizeBin
                    );
                }
                break;

                // We are on the second patch and properties for the ribbon  
                // patch are defined at first patch. Therefore, use  
                // shadowHistogram (histogram first patch) for THIS_PATCH and 
                //  patchHistogram (histogram second patch) for SHADOW_PATCH
                case THIS_PATCH:
                {
                    profile = computeProfileFromHistograms
                    (
                        shadowHistogram,
                        shadowHistogram,
                        halfMinSizeBin
                    );
                }
                break;

                case SHADOW_PATCH:
                 {
                    profile = computeProfileFromHistograms
                    (
                        patchHistogram,
                        patchHistogram,
                        halfMinSizeBin
                    );
                }
                break;
                
                // Use discretisation bothPatches for initial calculation of the
                // ribbons.
                case USER_DEFINED:
                {
                    pointField tempProfile;
                    
                    tempProfile = computeProfileFromHistograms
                    (
                        patchHistogram,
                        shadowHistogram,
                        halfMinSizeBin
                    );
                    
                    // Calculate the ribbons according to the user definition
                    profile = calcUserDefinedRibbons(tempProfile);                    
                }
                break;
                
                default:
                {
                    FatalErrorIn
                    (
                        "pointField mixingInterfacePolyPatch::"
                        "calcRibbons() const"
                    )   << "Bad discretisation type: "  
                        << discretisationNames_[discretisation_] 
                        << " Available types: bothPatches, thisPatch, "
                        << "shadowPatch, userDefined"
                        << exit(FatalError);
                }
                break;
            }
            
            // replace lower bound by zero in case stackAxis = R
            if (stackAxis_ == STACK_R)
                profile[0].replace(stackAxis, 0);
            // replace lower bound by -Infinity in case stackAxis = Z
            else if (stackAxis_ == STACK_Z)
                profile[0].replace(stackAxis, -Inf);
            // replace upper bound by Infinity
            profile[profile.size() - 1].replace(stackAxis, Inf);
            
            return profile;
        }
        break;
        
        default:
        {
            FatalErrorIn
            (
                "pointField mixingInterfacePolyPatch::calcRibbons() const"
            )   << "Bad division type: "  << divisionNames_[patchDivision_] 
                << " Available types : uniform, manual, meshDependent"
                << exit(FatalError);
            
            // Dummy return
            pointField profile;
            return profile;
        }
        break;
    } 
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::mixingInterfacePolyPatch::mixingInterfacePolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    polyPatch(name, size, start, index, bm, patchType),
    shadowPatch_(""),
    axis_(vector::zero),
    origin_(vector::zero),
    orgsep_(),
    patchDivision_(DIVISION_UNKNOWN),
    stackAxis_(STACK_UNKNOWN),
    discretisation_(DISCRETISATION_UNKNOWN),
    nPlanes_(-1),
    nPlanesBl_(-1),
    gradingIf_(false),
    sep_(),
    mxpFaces_(),
    mxpWeights_(),
    divided_(false),
    cartesianToCylindrical_(),
    cylindricalToCartesian_(),
    rotationTensorsUpdated_(false),
    master_(false)
{
    // mixingInterface is not constraint type so add mixingInterface group
    // explicitly
    if (inGroups().find(typeName) == -1)
    {
        inGroups().append(typeName);
    }
}


Foam::mixingInterfacePolyPatch::mixingInterfacePolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    polyPatch(name, dict, index, bm, patchType),
    shadowPatch_(dict.lookup("shadowPatch")),
    axis_(dict.lookupOrDefault<vector>("axis", vector::zero)),
    origin_(dict.lookupOrDefault<vector>("origin", vector::zero)),
    orgsep_(),
    patchDivision_(DIVISION_UNKNOWN),
    stackAxis_(STACK_UNKNOWN),
    discretisation_(DISCRETISATION_UNKNOWN),
    nPlanes_(-1),
    nPlanesBl_(-1),
    gradingIf_(false),
    sep_(),
    mxpFaces_(),
    mxpWeights_(),
    divided_(false),
    cartesianToCylindrical_(),
    cylindricalToCartesian_(),
    rotationTensorsUpdated_(false),
    master_(false)
{
    // Normalize axis.
    if (mag(axis_) == 0)
    {
        FatalErrorIn
        (
            "Foam::mixingInterfacePolyPatch::mixingInterfacePolyPatch"
        )   << "\n    " << name 
            << ": Axis needs to be nonzero vector." 
            << exit(FatalError);
    }
    axis_ /= mag(axis_);
    
    // Get axis_ with positiv sign in case axis_ is negative. We don`t need a
    // negative axis definition in case of a negative rotation axis
    axis_ = cmptMultiply(axis_, axis_);
    
    // Check if axis_ points in x-, y-, or z-direction. At the moment an
    // arbitrarily oriented rotation axis is not possible
    if 
    (
        axis_ != vector(1,0,0)
     && axis_ != vector(0,1,0)
     && axis_ != vector(0,0,1)
    )
    {
        FatalErrorIn
        (
            "Foam::mixingInterfacePolyPatch::mixingInterfacePolyPatch"
        )   << "\n    " << name
            << ": axis must point in x-, y-, or z-direction "
            << exit(FatalError);
    }
    
    if
    (
        shadowPatchID() == -1
     || 
        (
            shadowPatchID() != -1 
         && !boundaryMesh()[shadowPatchID()].inGroup(type())
        )
    )
    {
        // We are on the master side
        master_ = true;
    }

    
    // Read ribbon patch definition at first patch
    if (master())
    {
        // Stack axis must be read
        stackAxis_ = stackAxisNames_.read
        (
            dict.subDict("ribbonPatch").lookup("stackAxis")
        );
                
        // Avoid problems if case is decomposed and patches are on
        // different processors. Calculate the ribbons at the start.
        // After that, ribbons are read from the master patch
        if (dict.findEntry("sep"))
        {
            dict.lookup("sep") >> sep_;
        }
      
        else
        {            
            patchDivision_ = divisionNames_.read
            (
                dict.subDict("ribbonPatch").lookup("division")
            );

            switch (patchDivision_)
            {
                case UNIFORM:
                {
                    dict.subDict("ribbonPatch").lookup("planes") >> nPlanes_;
                    if (nPlanes_ < 1)
                    {
                        FatalErrorIn
                        (
                            "Foam::mixingInterfacePolyPatch::"
                            "mixingInterfacePolyPatch"
                        )   << "\n    " << name 
                            << ": There must be at least one mixing plane!" 
                            << exit(FatalError);
                    }
                }
                break;

                case MANUAL:
                {
                    dict.subDict("ribbonPatch").lookup("distances") >> orgsep_;
                }
                break;

                case MESH_DEPENDENT:
                {
                    discretisation_ = discretisationNames_.read
                    (
                        dict.subDict("ribbonPatch").lookup("discretisation")
                    );
                    
                    if (discretisation_ == USER_DEFINED)
                    {
                        dict.subDict("ribbonPatch").lookup("planes")
                            >> nPlanes_;
                        
                        dict.subDict("ribbonPatch").lookup("planesBl")
                            >> nPlanesBl_;
                        
                        dict.subDict("ribbonPatch").lookup("gradingIf")
                            >> gradingIf_;
                    }
                }
                break;
                
                default:
                {
                    FatalErrorIn
                    (
                        "Foam::mixingInterfacePolyPatch::"
                        "mixingInterfacePolyPatch"
                    )   << "\n    " << name 
                        << ": Bad division type: " 
                        << divisionNames_[patchDivision_] 
                        << " Available types: uniform, manual, meshDependent"
                        << exit(FatalError);
                }
                break;
            }
        }
    }
     
    // Build the ribbons if we are on the second patch. We can´t build the  
    // ribbons if we are on the first patch, since we have no access to the
    // second patch at this point
    else
    {
        // Get ribbon patch definition from first patch
        stackAxis_ = shadowPatch().stackAxis_;
        
        patchDivision_ = shadowPatch().patchDivision_;
        
        discretisation_ = shadowPatch().discretisation_;
        
        orgsep_ = shadowPatch().orgsep_;
        
        nPlanes_ = shadowPatch().nPlanes_;
        
        nPlanesBl_ = shadowPatch().nPlanesBl_;
        
        gradingIf_ = shadowPatch().gradingIf_;
        
        // Avoid problems if case is decomposed and patches are on
        // different processors. Calculate the ribbons at the start.
        // After that, ribbons are read from the master patch

        // Ribbons are read from master patch?
        if (shadowPatch().sep_.size() != 0)
        {
            sep_ = shadowPatch().sep_;
        }
        
        // Calculate Ribbons at the start
        else
        {
            const direction stackAxis = getStackDirection();
            
            const pointField profile = calcRibbons();
            
            sep_ = profile.component(stackAxis);
            
            // Set ribbons of first patch identical to second patch
            shadowPatch().sep_ = sep_;
            
            if (debug > 1 && patchDivision_ == MESH_DEPENDENT)
            {   
                pointField patchProfile = profile;

                if (axis().component(0) == 1)
                {
                    patchProfile.replace(0, profile.component(2));
                    patchProfile.replace(1, profile.component(0));
                    patchProfile.replace(2,profile.component(1));
                }
                else if (axis().component(1) == 1)
                {
                    patchProfile.replace(0,profile.component(1));
                    patchProfile.replace(1, profile.component(2));
                    patchProfile.replace(2, profile.component(0));
                }


                // Write control points into VTK
                OFstream os
                (
                    shadowPatch_ + ".vtk"
                );  

                Info<< "Writing ribbon control points to " << os.name() << nl;

                os  << "# vtk DataFile Version 2.0" << nl
                    << "ribbonControlPoints" << nl
                    << "ASCII" << nl
                    << "DATASET POLYDATA" << nl;

                const label nPoints = patchProfile.size() - 2;

                os  << "POINTS " << nPoints << " float" << nl;

                for (label i=0; i<patchProfile.size(); i++)
                {
                    if (i > 0 && i < patchProfile.size() - 1)
                    {
                        os  << float(patchProfile[i].component(0)) << ' '
                            << float(patchProfile[i].component(1)) << ' '
                            << float(patchProfile[i].component(2)) << nl;
                    }
                }

                os  << "VERTICES " << nPoints << ' ' << 2*nPoints << nl;
                for (label id = 0; id < nPoints; ++id)
                {   
                    os  << 1 << ' ' << id << nl;
                } 
            }
        }       
    }
    
    // Check valid stackAxis type
    if (stackAxis_ == STACK_UNKNOWN)
    {
        FatalErrorIn
        (
            "Foam::mixingInterfacePolyPatch::"
            "mixingInterfacePolyPatch"
        )   << "\n    " << name 
            << ": Bad stackAxis type: "  << stackAxisNames_[stackAxis_] 
            << " Available types: R, Z"
            << exit(FatalError);
    }
    
    // mixingInterface is not constraint type
    // so add mixingInterface group explicitly
    if (inGroups().find( typeName) == -1)
    {
        inGroups().append(typeName);
    }
    
    if (debug)
    {
        Info << "Created mixingInterface patch " << name 
        << " with the following attributes:" << endl;
        Info << "  origin:         " << origin_ << endl;
        Info << "  axis:           " << axis_ << endl;
        Info << "  division:       " << divisionNames_[patchDivision_] << endl;
        Info << "  stackAxis:      " << stackAxisNames_[stackAxis_] << endl;  
        
        if (patchDivision_ == MESH_DEPENDENT)
        {
            Info << "  discretisation: " 
                << discretisationNames_[discretisation_] << endl;
            if (discretisation_ == USER_DEFINED)
            {
                Info << "  planes:         " << nPlanes_ << endl;
                Info << "  planesBl:       " << nPlanesBl_ << endl;
                Info << "  gradingIf:      " << gradingIf_ << endl;
            }
        }
        if (patchDivision_ == UNIFORM)
            Info << "  planes:         " << nPlanes_ << endl;
        if (sep_.size() != 0)
            Info << "  separation distances list: " << sep_ << endl;
        Info << endl;
    }
}


Foam::mixingInterfacePolyPatch::mixingInterfacePolyPatch
(
    const mixingInterfacePolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    polyPatch(pp, bm),
    shadowPatch_(pp.shadowPatch_),
    axis_(pp.axis_),
    origin_(pp.origin_),
    orgsep_(pp.orgsep_),
    patchDivision_(pp.patchDivision_),
    stackAxis_(pp.stackAxis_),
    discretisation_(pp.discretisation_),
    nPlanes_(pp.nPlanes_),
    nPlanesBl_(pp.nPlanesBl_),
    gradingIf_(pp.gradingIf_),
    sep_(pp.sep_),
    mxpFaces_(pp.mxpFaces_),
    mxpWeights_(pp.mxpWeights_),
    divided_(pp.divided_),
    cartesianToCylindrical_(pp.cartesianToCylindrical_),
    cylindricalToCartesian_(pp.cylindricalToCartesian_),
    rotationTensorsUpdated_(pp.rotationTensorsUpdated_),
    master_(pp.master_)
{}


Foam::mixingInterfacePolyPatch::mixingInterfacePolyPatch
(
    const mixingInterfacePolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    polyPatch(pp, bm, index, newSize, newStart),
    shadowPatch_(pp.shadowPatch_),
    axis_(pp.axis_),
    origin_(pp.origin_),
    orgsep_(pp.orgsep_),
    patchDivision_(pp.patchDivision_),
    stackAxis_(pp.stackAxis_),
    discretisation_(pp.discretisation_),
    nPlanes_(pp.nPlanes_),
    nPlanesBl_(pp.nPlanesBl_),
    gradingIf_(pp.gradingIf_),
    sep_(pp.sep_),
    mxpFaces_(pp.mxpFaces_),
    mxpWeights_(pp.mxpWeights_),
    divided_(pp.divided_),
    cartesianToCylindrical_(pp.cartesianToCylindrical_),
    cylindricalToCartesian_(pp.cylindricalToCartesian_),
    rotationTensorsUpdated_(pp.rotationTensorsUpdated_),
    master_(pp.master_)
{}


Foam::mixingInterfacePolyPatch::mixingInterfacePolyPatch
(
    const mixingInterfacePolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    polyPatch(pp, bm, index, mapAddressing, newStart),
    shadowPatch_(pp.shadowPatch_),
    axis_(pp.axis_),
    origin_(pp.origin_),
    orgsep_(pp.orgsep_),
    patchDivision_(pp.patchDivision_),
    stackAxis_(pp.stackAxis_),
    discretisation_(pp.discretisation_),
    nPlanes_(pp.nPlanes_),
    nPlanesBl_(pp.nPlanesBl_),
    gradingIf_(pp.gradingIf_),
    sep_(pp.sep_),
    mxpFaces_(pp.mxpFaces_),
    mxpWeights_(pp.mxpWeights_),
    divided_(pp.divided_),
    cartesianToCylindrical_(pp.cartesianToCylindrical_),
    cylindricalToCartesian_(pp.cylindricalToCartesian_),
    rotationTensorsUpdated_(pp.rotationTensorsUpdated_),
    master_(pp.master_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mixingInterfacePolyPatch::~mixingInterfacePolyPatch()
{
    
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mixingInterfacePolyPatch::initGeometry(PstreamBuffers& pBufs)
{
    //polyPatch::initCalcGeometry(pBufs);
    polyPatch::initGeometry(pBufs);
}


void Foam::mixingInterfacePolyPatch::calcGeometry(PstreamBuffers& pBufs)
{
    polyPatch::calcGeometry(pBufs);
}


void Foam::mixingInterfacePolyPatch::initMovePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    polyPatch::initMovePoints(pBufs, p);
}


void Foam::mixingInterfacePolyPatch::movePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    polyPatch::movePoints(pBufs, p);
}


void Foam::mixingInterfacePolyPatch::initUpdateMesh(PstreamBuffers& pBufs)
{
    polyPatch::initUpdateMesh(pBufs);
}


void Foam::mixingInterfacePolyPatch::updateMesh(PstreamBuffers& pBufs)
{
    polyPatch::updateMesh(pBufs);
}


Foam::label Foam::mixingInterfacePolyPatch::shadowPatchID() const
{
    return boundaryMesh().findPatchID(shadowPatchName());
}


const Foam::mixingInterfacePolyPatch&
Foam::mixingInterfacePolyPatch::shadowPatch() const
{
    return refCast<const mixingInterfacePolyPatch>
        (boundaryMesh()[shadowPatchID()]);
}


Foam::scalar Foam::mixingInterfacePolyPatch::capArea
(
    vector const & A, 
    vector const & B, 
    Foam::scalar R
) const
{
    // get real distance between A and B
    scalar AB = mag(B - A);
    
    // skip empty surfaces
    if (AB == 0 or R == 0)
        return 0;
    
    // angle between the line segment AB and the axis
    scalar cosAngle = (axis_ & (B - A)) / AB;
    
    // length of the axis-orthogonal projection of the line segment AB
    scalar c = mag(axis_ ^ (B - A));
    
    // calculate surface of the axis-orthogonal projection of the cap
    //      see http://en.wikipedia.org/wiki/Circular_segment
    scalar theta = 2*asin(0.5*c/R);
    scalar area = 0.5*R*R*(theta - Foam::sin(theta));
    
    // scale area by un-projection
    return mag(cosAngle) < 1. ? area / Foam::sqrt(1. - cosAngle*cosAngle) : 0;
}

int Foam::mixingInterfacePolyPatch::intersectLinePlane
(
    vector const & A, vector const & B, // line points
    scalar h,                           // cutting plane position
    vector * intsct                     // intersections
) const
{
    // line equation : X = A + (B - A)t
    // plane equation : n . (X - O) = 0
    // => solve for 't' by elimination of X
    
    // skip tiny lines
    if (mag(B - A) < 1e-30 or mag(axis_ & (B - A)) < 1e-30)
        return 0;
    
    // plane-axis intersection (= plane origin)
    vector O = origin_ + h * axis_;
    
    // get intersection parameter
    scalar t = ((axis_ & (O - A))) / (axis_ & (B - A));
    
    // store intersections
    int n = 0;
    if (0 < t and t < 1) intsct[n++] = A + (B - A) * t;
    
    // return number of intersections
    return n;
}

int Foam::mixingInterfacePolyPatch::intersectLineCylinder
(
    vector const & A, vector const & B, // line points
    scalar r,                           // cutting cylinder radius
    vector * intsct                     // intersections
) const
{
    // line equation : X = A + (B - A)t
    // circle equation : |(X - O) x Axis| = r
    // => solve for 't' by elimination of X
    
    // prepare some variables for the quadratic equation
    vector u = axis_ ^ (B - A);
    vector v = axis_ ^ (A - origin_);
    scalar a = magSqr(u);
    scalar b = 2. * (u & v);
    scalar c = magSqr(v) - r*r;
    scalar D = b*b - 4*a*c;
    
    // no intersection, or touch only
    if (D <= 0)
        return 0;
    
    // two intersections
    scalar t1 = (-b - sqrt(D))/(2*a);
    scalar t2 = (-b + sqrt(D))/(2*a);
    
    // store intersections
    int n = 0;
    if (0 < t1 and t1 < 1) intsct[n++] = A + (B - A) * t1;
    if (0 < t2 and t2 < 1) intsct[n++] = A + (B - A) * t2;
    
    // return number of intersections
    return n;
}

Foam::scalar Foam::mixingInterfacePolyPatch::calcOverlap
(
    vector const & a, 
    vector const & b,
    vector const & c, scalar r
) const
{
    // skip null cylinders
    if 
    ( (r == 0 and stackAxis_ == STACK_R) )
        return 0;
    
    // calculate area of the triangle
    scalar area = 0.5 * mag((b-a) ^ (c-a));
    
    // distance of triangle points from axis (in radial case) 
    // or along the axis (in axial case)
    scalar ra = 0, rb = 0, rc = 0;
    
    // get axial/radial distances
    //MXP_RADIAL
    if (stackAxis_ == STACK_R)
    {
        // distances from the axis
        ra = mag((a - origin_) ^ axis_);
        rb = mag((b - origin_) ^ axis_);
        rc = mag((c - origin_) ^ axis_);
    }
    //MXP_AXIAL
    else
    {
        // distances along the axis
        ra = (a - origin_) & axis_;
        rb = (b - origin_) & axis_;
        rc = (c - origin_) & axis_;
    }
    
    // determine position of points w.r.t. the cutting plane
    bool a_inside = (ra <= r), a_outside = (ra >= r);
    bool b_inside = (rb <= r), b_outside = (rb >= r);
    bool c_inside = (rc <= r), c_outside = (rc >= r);
    
    // if all vertices are inside, use the whole area as overlap
    if (a_inside and b_inside and c_inside)
        return area;
    
    // if all vertices are outside AND this is a planar cut, return clean zero
    if (a_outside and b_outside and c_outside and stackAxis_ == STACK_Z)
        return 0;
    
    // intersection info
    #define AB 0x01
    #define BC 0x02
    #define CA 0x04
    
    // overlap polygon points and lines they lay on
    vectorList vertices;
    labelList lines;
    
    // intersect all triangle's lines by the separation surface:
    // plane (for axial) or cylinder (for radial)
    if (stackAxis_ == STACK_Z)
    {
        vector intsct;
        
        if (a_inside)
        {
            vertices.append(a);
            lines.append(CA | AB);
        }
        
        if (intersectLinePlane(a,b,r,&intsct))
        {
            vertices.append(intsct);
            lines.append(AB);
        }
        
        if (b_inside)
        {
            vertices.append(b);
            lines.append(AB | BC);
        }
        
        if (intersectLinePlane(b,c,r,&intsct))
        {
            vertices.append(intsct);
            lines.append(BC);
        }
        
        if (c_inside)
        {
            vertices.append(c);
            lines.append(BC | CA);
        };
        
        if (intersectLinePlane(c,a,r,&intsct))
        {
            vertices.append(intsct);
            lines.append(CA);
        }
    }
    //MXP_RADIAL
    if (stackAxis_ == STACK_R)
    {
        vector intsct[2];
        int N;
        
        if (a_inside)
        {
            vertices.append(a);
            lines.append(CA | AB);
        }
        
        N = intersectLineCylinder(a,b,r,intsct);
        if (N > 0)
        {
            vertices.append(intsct[0]);
            lines.append(AB);
        }
        if (N > 1)
        {
            vertices.append(intsct[1]);
            lines.append(AB);
        }
        
        if (b_inside)
        {
            vertices.append(b);
            lines.append(AB | BC);
        }
        
        N = intersectLineCylinder(b,c,r,intsct);
        if (N > 0)
        {
            vertices.append(intsct[0]);
            lines.append(BC);
        }
        if (N > 1)
        {
            vertices.append(intsct[1]);
            lines.append(BC);
        }
        
        if (c_inside)
        {
            vertices.append(c);
            lines.append(BC | CA);
        }
        
        N = intersectLineCylinder(c,a,r,intsct);
        if (N > 0)
        {
            vertices.append(intsct[0]);
            lines.append(CA);
        }
        if (N > 1)
        {
            vertices.append(intsct[1]);
            lines.append(CA);
        }
    }
    
    // trivial cases that should not happen in accurate arithmetic
    if (vertices.size() < 2)
        return 0;
    
    // two-point quadratic polygon
    if (vertices.size() == 2)
    {
        if ( (stackAxis_ == STACK_R and lines[0] == lines[1]) )
        {
            // the area is equal to the area of the circular cap
            return capArea(vertices[0], vertices[1], r);
        }
        else
        {
            // should not happen in accurate arithmetic
            return 0;
        }
    }
    
    // the result is a convex quadratic polygon; calculate its centre
    vector centre = vector::zero;
    forAll (vertices, I)
        centre += vertices[I];
    centre /= vertices.size();
    
    // for all quadratic polygon edges
    scalar weight = 0;
    for (int i = 0; i < vertices.size(); i++)
    {
        // get edge points
        vector const & v1 = vertices[i];
        vector const & v2 = vertices[(i + 1) % vertices.size()];
        
        // get edge IDs
        int edge1 = lines[i];
        int edge2 = lines[(i + 1) % lines.size()];
        
        // calculate area of the sub-triangle
        weight += 0.5 * mag((v1 - centre) ^ (v2 - centre));
        
        // add optional cap surface (for points on distinct edges)
        if ( ((edge1 & edge2) == 0 and stackAxis_ == STACK_R) )
        {
            weight += capArea(v1, v2, r);
        }
    }
    
    return weight;
}

Foam::scalar Foam::mixingInterfacePolyPatch::calcWeight
(
    scalar r1,
    scalar r2,
    label iface
) const
{
    // get face centre
    vector Cf = this->faceCentres()[iface];
    
    /// DEBUG : sharp interface (no interpolation)
//       if 
//        (
//            mode_ == MXP_AXIAL
//            and r1 <=    ((Cf - origin_) & axis_)
//            and    ((Cf - origin_) & axis_) < r2
//        ) return 1;
//       if 
//        (
//            mode_ == MXP_RADIAL
//            and r1 <= mag((Cf - origin_) ^ axis_)
//            and mag((Cf - origin_) ^ axis_) < r2
//        ) return 1;
//       return 0;
    ///
    
    // get all mesh points
    pointField const & pts = this->boundaryMesh().mesh().points();;
    
    // get point labels of the i-th face
    face const & f = this->boundaryMesh().mesh().faces()[this->start() + iface];
    
    // weights
    scalar area_in_r1 = 0, area_in_r2 = 0, area_full = 0;
    
    // calculate full area of the triangle
    forAll (f, i)
    {
        // get vertices
        vector a = Cf;
        vector b = pts[f[i]];
        vector c = pts[f[(i+1) % f.size()]];
        
        area_full += 0.5 * mag((b-a) ^ (c-a));
    }
    
    // split face into triangle fan
    forAll (f, i)
    {
        // get vertices
        vector a = Cf;
        vector b = pts[f[i]];
        vector c = pts[f[(i+1) % f.size()]];
        
        // calculate intersections
        area_in_r1 += calcOverlap(a, b, c, r1);
        area_in_r2 += calcOverlap(a, b, c, r2);
    }
    
    return (area_in_r2 - area_in_r1) / area_full;
}

void Foam::mixingInterfacePolyPatch::planeFaces
(
    label iplane,
    List<label> & faces,
    List<scalar> & weights
) const
{
    // Check that mixing plane has been correctly initialized
    // (= contains at least single strip and mode).
    if (sep_.size() == 0)
    {
        FatalErrorIn("Foam::mixingInterfacePolyPatch::planeFaces")
            << "\n    Mixing plane has not been initialized." 
            << exit(FatalError);
    }
  
    // Check that the requested strip is available.
    if (iplane >= sep_.size() - 1)
    {
        FatalErrorIn("Foam::mixingInterfacePolyPatch::planeFaces")
            << "\n    The requested mixing plane number " << iplane 
            << " is not available (max " << sep_.size() - 2 << ")." 
            << exit(FatalError);
    }
    
    // Initialize lists.
    faces.clear();
    weights.clear();
    
    // Loop over all faces.
    forAll (this->faceCentres(), i)
    {
        // Calculate area fraction on this mixing plane.
        scalar weight = calcWeight(sep_[iplane], sep_[iplane + 1], i);
        
        // If the area is non-zero, add the face to the list.
        if (weight > 0)
        {
            faces.append(i);
            weights.append(weight);
        }
    }
}

Foam::tensorField Foam::mixingInterfacePolyPatch::calcRotationTensor
(
    const vector& rotationAxis,
    const vectorField& vi,
    const vectorField& vf
) const
{
    
    scalar magRotAxis = mag(rotationAxis);

    const vector k = rotationAxis/magRotAxis;
    const scalar k1 = k[0];
    const scalar k2 = k[1];
    const scalar k3 = k[2];

    // [k]_x
    const tensor rotationT (
          0, -k3,  k2,
         k3,   0, -k1,
        -k2,  k1,   0
    );

    // kk' - I
    const tensor projectionT (
        k1*k1 - 1.0, k1*k2,       k1*k3,
        k2*k1,       k2*k2 - 1.0, k2*k3,
        k3*k1,       k3*k2,       k3*k3 - 1.0);

    // Project both vectors onto the plane defined by the rotation axis
    vectorField nvi = -(projectionT & vi);
    vectorField nvf = -(projectionT & vf);
    nvi = nvi/mag(nvi);
    nvf = nvf/mag(nvf);

    const vectorField crossNviNvf = nvi ^ nvf;
    const scalarField cosTheta = nvi & nvf;
    const scalarField sinTheta = mag(crossNviNvf) * sign(crossNviNvf & k);

    const tensorField I_F(vi.size(), I);
    const tensorField rotationT_F(vi.size(), rotationT);
    const tensorField projectionT_F(vi.size(), projectionT);

    return I_F + sinTheta*rotationT_F + (1 - cosTheta)*projectionT_F;
}

void Foam::mixingInterfacePolyPatch::write(Ostream& os) const
{
    polyPatch::write(os);
    
    os.writeKeyword("shadowPatch") << shadowPatch_ 
        << token::END_STATEMENT << nl;
    
    os.writeKeyword("axis") << axis_ << token::END_STATEMENT << nl;
    
    os.writeKeyword("origin") << origin_ << token::END_STATEMENT << nl;
      
    if (master())
    {
        // Avoid problems if case is decomposed and patches are on
        // different processors. Calculate the ribbons at the start.
        // After that, ribbons are read from the master patch
        if (sep_.size() != 0)
        {
            os.writeKeyword("sep") << "( ";
            for (label i=0; i<sep_.size(); i++)
                os << sep_[i] << " ";
            os << ")" << token::END_STATEMENT << nl;
        }
  
        dictionary ribbonPatchDict("ribbonPatch");
        
        ribbonPatchDict.add("division", divisionNames_[patchDivision_]);
        
        ribbonPatchDict.add("stackAxis", stackAxisNames_[stackAxis_]);
        
        if (patchDivision_ == UNIFORM)
        {
            ribbonPatchDict.add("planes", nPlanes_);
        }
        
        if (patchDivision_ == MANUAL)
        {
            ribbonPatchDict.add("distances", orgsep_);
        }
        
        if (patchDivision_ == MESH_DEPENDENT)
        {
            ribbonPatchDict.add
            (
                "discretisation",
                discretisationNames_[discretisation_]
            );
            
            if (discretisation_ == USER_DEFINED)
            {
                ribbonPatchDict.add("planes", nPlanes_);
                
                ribbonPatchDict.add("planesBl", nPlanesBl_);
                
                ribbonPatchDict.add("gradingIf", gradingIf_);
            }
        }
               
        os.writeKeyword("ribbonPatch") << ribbonPatchDict << nl;
    }
}


// ************************************************************************* //
