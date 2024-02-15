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
    Foam::mixingInterfacePolyPatch

Description
    Patch type to use with "mixingInterface" boundary condition.
    
    Subdict "ribbonPatch" must be defined only for the first patch in file 
    constant/polyMesh/boundary. For the second patch, only the entries "type", 
    "shadowPatch", "axis" and "origin" are required. 
   
    There are three possibilities for the patch division:
    manual, uniform, meshDependent

    In case of the manual division, the distances between the ribbons must be
    provided. 
   
    If the automatic uniform division is requested, the boundary condition
    will determine extremal point coordinates (from or along the axis) of the 
    patch and divide the patch in the correct direction (radial / axial) between
    the extremal positions.
  
    In case of the meshDependet division, the ribbon coordinates are derived 
    from the mesh. There are four discretisation possibilities. In case of 
    discretisation thisPatch, shadowPatch or bothPatches, the ribbons are  
    derived from the first patch the second patch or both patches, respectively.  
    In case of discretisation userDefined, the ribbons are first calculated with 
    discretisation bothPatch and than selected according to the user 
    specification.   
  
    Example 1 (manual division of patch):
    
        stator_rotor_mxp
        {
            type            mixingInterfacePatch;
            shadowPatch     rotor_stator_mxp;
            axis            (0 0 1);    // default: (0 0 1)
            origin          (0 0 0);    // default: (0 0 0)
            ribbonPatch
            {
                division    manual;
                stackAxis   R;  
                distances   ( 0.1  0.2  0.3  0.4  0.5 );          
            }
        }
    
    
    Example 2 (automatic uniform division of the patch):
    
        stator_rotor_mxp
        {
            type            mixingInterfacePatch;
            shadowPatch     rotor_stator_mxp;
            axis            (0 0 1);    // default: (0 0 1)
            origin          (0 0 0);    // default: (0 0 0)
            ribbonPatch
            {
                division    uniform;
                stackAxis   R;  
                planes      20;
            }
        }
  
    Example 3 (mesh dependent division of the patch):
    
        stator_rotor_mxp
        {
            type            mixingInterfacePatch;
            shadowPatch     rotor_stator_mxp;
            axis            (0 0 1);    // default: (0 0 1)
            origin          (0 0 0);    // default: (0 0 0)
            ribbonPatch
            {
                division        meshDependent;
                stackAxis       R;  
                discretisation  thisPatch;  
            }
        }
   
    Example 4 (mesh dependent user-defined division of the patch):
    
        stator_rotor_mxp
        {
            type            mixingInterfacePatch;
            shadowPatch     rotor_stator_mxp;
            axis            (0 0 1);    // default: (0 0 1)
            origin          (0 0 0);    // default: (0 0 0)
            ribbonPatch
            {
                division        meshDependent;
                stackAxis       R;  
                discretisation  userDefined;
                planes          20;     // Total number of planes
                planesBl        4;      // Number of planes in boundary layer 
                gradingIf       true;   // grading between boundary layer and 
                                        // internal field? default false
            }
        }
  
Note
    Developed by CFD support in 2015.
    Extended by IHS, University of Stuttgart, in 2021, 2023.

SourceFiles
    mixingInterfacePolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef mixingInterfacePolyPatch_H
#define mixingInterfacePolyPatch_H

#include "polyPatch.H"
#include "fvPatch.H"
#include "coordinateSystem.H"
#include "cylindricalCS.H"
#include "Enum.H"

#include <list>
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;

/*---------------------------------------------------------------------------*\
                      Class mixingInterfacePolyPatch Declaration
\*---------------------------------------------------------------------------*/

class mixingInterfacePolyPatch
:
    public polyPatch
{
    // Private data types
    
    typedef std::map<Foam::scalar, std::list<Foam::point> > profileHistogram;
    
    // Private Member Functions
    
    void updateRotationTensors() const;
    
    void patchDivision() const;
    
    direction getStackDirection() const;
    
    scalar calcGradingDist
    (
        const scalar& varG,
        const scalar& varX
    ) const;
    
    pointField calcUserDefinedRibbons(const pointField tempProfile) const;
    
    pointField computeProfileFromHistograms
    (
        const profileHistogram& patchHisto,
        const profileHistogram& shadowHisto,
        const scalar halfSizeBin
    ) const;
    
    void updateProfileHistogram
    (
        profileHistogram& histo,
        const point& profileCoord,  // 3D point reference
        const direction dir,        // Sorting dimension 0: x, 1: y, 2: z
        const scalar halfSizeBin    // half size of min width for histogram bins
    ) const;
    
    pointField calcRibbons() const;
            
public:
    
    // Type enumerations
        
        //- Patch division mode
        enum division
        {
            UNIFORM,
            MANUAL,
            MESH_DEPENDENT,
            DIVISION_UNKNOWN
        };
        
        //- Stack axis
        enum stackAxis
        {   
            STACK_R,
            STACK_Z,
            STACK_UNKNOWN
        };
        
        //- Discretisation for patch division mode MESH_DEPENDENT 
        enum discretisation
        {
            BOTH_PATCHES,
            THIS_PATCH,
            SHADOW_PATCH,
            USER_DEFINED,
            DISCRETISATION_UNKNOWN
        };
        
        //- Patch division names
        static const Enum<division> divisionNames_;
        
        //- Stack axis names
        static const Enum<stackAxis> stackAxisNames_;
        
        //- Discretisation names
        static const Enum<discretisation> discretisationNames_;
        
        
protected:

    // User settings
    
  
        //- Shadow patch
        word shadowPatch_;

        //- Rotation axis
        vector axis_;

        //- Rotation origin
        vector origin_;
        
        //- User defined distance
        List<scalar> orgsep_;
        
        //- Patch division Mode
        division patchDivision_;
        
        //- stack axis for mesh dependent ribbons
        stackAxis stackAxis_;
        
        //- discretisation type of mesh dependent ribbons
        discretisation discretisation_;
        
        //- Number of planes.
        label nPlanes_;
        
        //- Number of planes in the boundary layer
        label nPlanesBl_;
        
        //- Grading of the ribbons between boundary layer?
        bool gradingIf_;
        
    // Internal data
        
        //- List of mixing plane separation positions with extreme boundaries appended on ends.
        mutable List<scalar> sep_;
        
        //- List of faces belonging to each mixing plane.
        mutable List< List<label> > mxpFaces_;
        
        //- List of face weights for smooth transitions between the mixing planes.
        mutable List< List<scalar> > mxpWeights_;
        
        mutable bool divided_;
        
        //- Transformation tensors
        mutable tensorField cartesianToCylindrical_;
        
        mutable tensorField cylindricalToCartesian_;
        
        mutable bool rotationTensorsUpdated_;
        
        //- Master side?
        mutable bool master_;
        
        //- Initialise the calculation of the patch geometry
        virtual void initGeometry(PstreamBuffers&);

        //- Calculate the patch geometry
        virtual void calcGeometry(PstreamBuffers&);

        //- Initialise the patches for moving points
        virtual void initMovePoints(PstreamBuffers&, const pointField&);

        //- Correct patches after moving points
        virtual void movePoints(PstreamBuffers&, const pointField&);

        //- Initialise the update of the patch topology
        virtual void initUpdateMesh(PstreamBuffers&);

        //- Update of the patch topology
        virtual void updateMesh(PstreamBuffers&);


public:

    //- Runtime type information
    TypeName("mixingInterfacePatch");


    // Constructors

        //- Construct from components
        mixingInterfacePolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType
        );

        //- Construct from dictionary
        mixingInterfacePolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType
        );

        //- Construct as copy, resetting the boundary mesh
        mixingInterfacePolyPatch
        (
            const mixingInterfacePolyPatch&,
            const polyBoundaryMesh&
        );

        //- Construct given the original patch and resetting the
        //  face list and boundary mesh information
        mixingInterfacePolyPatch
        (
            const mixingInterfacePolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        );

        //- Construct given the original patch and a map
        mixingInterfacePolyPatch
        (
            const mixingInterfacePolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const labelUList& mapAddressing,
            const label newStart
        );

        //- Construct and return a clone, resetting the boundary mesh
        virtual autoPtr<polyPatch> clone(const polyBoundaryMesh& bm) const
        {
            return autoPtr<polyPatch>(new mixingInterfacePolyPatch(*this, bm));
        }

        //- Construct and return a clone, resetting the face list
        //  and boundary mesh
        virtual autoPtr<polyPatch> clone
        (
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        ) const
        {
            return autoPtr<polyPatch>
            (
                new mixingInterfacePolyPatch(*this, bm, index, newSize, newStart)
            );
        }

        //- Construct and return a clone, resetting the face list
        //  and boundary mesh
        virtual autoPtr<polyPatch> clone
        (
            const polyBoundaryMesh& bm,
            const label index,
            const labelUList& mapAddressing,
            const label newStart
        ) const
        {
            return autoPtr<polyPatch>
            (
                new mixingInterfacePolyPatch
                (
                    *this,
                    bm,
                    index,
                    mapAddressing,
                    newStart
                )
            );
        }


    //- Destructor
    virtual ~mixingInterfacePolyPatch();


    // Member functions
    
        //- Get roration axis.
        vector axis() const { return axis_; }
        
        //- Get position of origin.
        vector origin() const { return origin_; }
        
        //- Get shadow patch name.
        word shadowPatchName() const { return shadowPatch_; }
        
        //- Get shadow patch ID
        label shadowPatchID() const;
        
        //- Get shadow patch
        const mixingInterfacePolyPatch& shadowPatch() const;
        
        //- Get number of mixing plane strips / rings.
        label numPlanes() const { return sep_.size() - 1; }
        
        //_ master patch?
        bool master() const {return master_;}
        
        //- Get list of labels of faces in a given mixing plane strip / ring.
        void planeFaces
        (
            label iplane,
            List<label> & faces,
            List<scalar> & weights
        ) const;
        
        //- Calculate weight of face in mixing plane bounded by 'a' and 'b'.        
        scalar calcWeight
        (
            scalar a, scalar b,
            label i
        ) const;
        
        //- Calculate triangle sub-area bounded by 'r'.
        scalar calcOverlap
        (
            vector const & a,
            vector const & b,
            vector const & c,
            scalar r
        ) const;
        
        //- Calculate zero or two line-cylinder intersections and store in 'intsct'.
        int intersectLineCylinder
        (
            vector const & A, vector const & B, // line points
            scalar r,                           // cylinder radius
            vector * intsct                     // intersections
        ) const;
        
        //- Calculate zero or one line-plane intersection and store in 'intsct'.
        int intersectLinePlane
        (
            vector const & A, vector const & B, // line points
            scalar r,                           // cylinder radius
            vector * intsct                     // intersections
        ) const;
        
        //- Calculate circular segment (cap) area.
        scalar capArea
        (
            vector const & A, vector const & B, // circular cap vertices (the circle's centre is on the axis).
            scalar r                            // circle radius (centre is in axis)
        ) const;
        
        //- Calculate rotation tensor Field with Rodrigues Rotation Formula
        tensorField calcRotationTensor
        (
            const vector& rotationAxis,
            const vectorField& vi,
            const vectorField& vf
        ) const;
        
        //- Interpolate from patch to ribbon
        template<class Type>
        tmp<Field<Type> > toProfile(const Field<Type>& pf) const;
        
        //- Interpolate from patch to ribbon to patch
        template<class Type>
        tmp<Field<Type> > circumferentialAverage(const Field<Type>& tpf) const;
        
        //- Interpolate from segment patch to ribbon to patch
        template<class Type>
        tmp<Field<Type> > interpolate(const Field<Type>& spf) const;
        
        template<class Type>
        tmp<Field<Type> > transformToCylindrical(const Field<Type>& tpf) const;
        
        template< class Type>
        tmp<Field<Type> > transformToCartesian(const Field<Type>& tpf) const;
        
        //- Write the polyPatch data as a dictionary
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "mixingInterfacePolyPatchTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
