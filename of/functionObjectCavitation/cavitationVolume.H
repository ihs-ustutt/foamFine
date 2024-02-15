/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
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
    Foam::functionObjects::fieldValues::cavitationVolume

Group
    grpFieldFunctionObjects

Description
    Provides the cavitation volume

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    cavitationVolume1
    {
        // Mandatory entries (unmodifiable)
        type            cavitationVolume;
        libs            (cavitationVolume);

        // Mandatory entries (runtime modifiable)
        field          (p);
        operation       <operationType>;
        regionType      <volRegion>;

        // Optional entries 
        patchNameAverage   DT_OUTLET;
        flowDirection      "-z";
        heightAxis         x;
        installDepth       8;
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Req'd | Dflt
      type         | Type name: cavitationVolume        | word |  yes  | -
      libs         | Library name: cavitationVolume     | word |  yes  | -
      field        | Names of operand field (p)         | wordList | yes | -
      regionType   | Face regionType: see below         | word |  yes  | -
      name         | Name for regionType                | word |  yes  | -
      operation    | Operation type: see below          | word |  yes  | -
      patchAverageName | see below                      | word |  no   | -
      flowDirection | see below                         | word |  no   | -
      heightAxis   | see below                          | word |  no   | -
      installDepth | see below                          |scalar|  no   | 0
    \endtable

    The inherited entries are elaborated in:
     - \link fieldValue.H \endlink
     - \link volRegion.H \endlink

    Options for the \c regionType entry:
    \plaintable
      cellZone     | requires a 'name' entry to specify the cellZone
      all          | all cells
    \endplaintable
    
    Options for the \c operation entry:
    \plaintable
      none                 | No operation
      sum                  | Sum
    \endplaintable
    
    patchAverageName:
    Name of the patch, where the area averaged value of the pressure is 
    used for pressure correction. Necessary for example if a tank is added
    to the machine. Needs entry flowDirection.
    
    flowDirection:
    Used to limit the area in which the cavitation volume is calculated. For
    example if a tank is added to the machine, the cavitation volume is not
    calculated in the tank.
    Options for the \c flowDirection entry:
    \plaintable
      x         
      -x           
      y
      -y
      z
      -z          
    \endplaintable
    
    heightAxis:
    Defines the direction in which the water depth decreases. For example
    if the entry is set to "x", the pressure decreases with increasing
    x-Coordinate.
    Options for the \c heightAxis entry:
    \plaintable
      x         
      -x           
      y
      -y
      z
      -z          
    \endplaintable 
 
    installDepth:
    Installation Depth of the machine (machine axis) used to correct the 
    pressure.
 

    Usage by the \c postProcess utility is not available.

See also
    - Foam::functionObject
    - Foam::functionObjects::fieldValues::fieldValue
    - Foam::functionObjects::volRegion

SourceFiles
    cavitationVolume.C
    cavitationVolumeTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_cavitationVolume_H
#define functionObjects_cavitationVolume_H

#include "fieldValue.H"
#include "volRegion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace fieldValues
{

/*---------------------------------------------------------------------------*\
                        Class cavitationVolume Declaration
\*---------------------------------------------------------------------------*/

class cavitationVolume
:
    public fieldValue,
    public volRegion
{
public:

    // Public Data Types

        //- Bitmask values for operation variants
        enum operationVariant
        {
            typeScalar = 0x100,   //!< Operation returns a scalar
        };

        //- Operation type enumeration
        enum operationType
        {
            // Normal operations

            opNone = 0,             //!< No operation
            opSum,                  //!< Sum of values
        };

        //- Operation type names
        static const Enum<operationType> operationTypeNames_;
        
        enum heightAxis
        {
            hXPos,             
            hXNeg,
            hYPos,
            hYNeg,
            hZPos,
            hZNeg,
        };
        
        //- Operation type names
        static const Enum<heightAxis> heightAxisNames_;
        
        enum flowDirection
        {
            diXPos,             
            diXNeg,
            diYPos,
            diYNeg,
            diZPos,
            diZNeg,
        };
        
        //- Operation type names
        static const Enum<flowDirection> flowDirectionNames_;



protected:

    // Protected Data

        //- Operation to apply to values
        operationType operation_;
        
        heightAxis heightAxis_;
        
        flowDirection flowDirection_;

        word patchAverageName_;
        
        bool useHeightAxis_;
        
        scalar installDepth_;
        
        bool patchAverage_;
        

    // Protected Member Functions

        //- True if the operation needs the cell-volume
        //bool usesVol() const noexcept;

        //- True if the operation variant uses mag
        //inline bool is_magOp() const noexcept;

        //- True if the operation variant uses a weight-field
        //inline bool is_weightedOp() const noexcept;

        //- True if field is non-empty on any processor.
        inline bool canWeight(const scalarField& fld) const;

        //- Return true if the field name is valid
        template<class Type>
        bool validField(const word& fieldName) const;

        //- Insert field values into values list
        template<class Type>
        tmp<Field<Type>> getFieldValues
        (
            const word& fieldName,
            const bool mandatory = false
        ) const;

        //- Apply the 'operation' to the values
        template<class Type>
        Type processValues
        (
            const Field<Type>& values,
            const scalarField& V
        ) const;

        //- Helper function to output field values
        label writeAll
        (
            const scalarField& V
        );

        //- Templated helper function to output field values
        template<class Type>
        bool writeValues
        (
            const word& fieldName,
            const scalarField& V
        );

        //- Filter a field according to cellIds
        template<class Type>
        tmp<Field<Type>> filterField(const Field<Type>& field) const;

        //- Output file header information
        virtual void writeFileHeader(Ostream& os) const;


public:

    //- Declare type-name, virtual type (with debug switch)
    TypeName("cavitationVolume");


    // Constructors

        //- Construct from name, Time and dictionary
        cavitationVolume
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Construct from name, objectRegistry and dictionary
        cavitationVolume
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict
        );

        //- No copy construct
        cavitationVolume(const cavitationVolume&) = delete;

        //- No copy assignment
        void operator=(const cavitationVolume&) = delete;


    //- Destructor
    virtual ~cavitationVolume() = default;


    // Member Functions

        //- Read from dictionary
        virtual bool read(const dictionary& dict);

        //- Calculate and write
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fieldValues
} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "cavitationVolumeTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
