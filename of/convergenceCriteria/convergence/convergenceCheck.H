
#ifndef convergenceCheck_H
#define convergenceCheck_H

#include "functionObject.H"
#include "Time.H"
#include "dictionary.H"
#include "DynamicList.H"
#include "OFstream.H"
#include "autoPtr.H"

namespace Foam
{

class convergenceCheck
:
    public functionObject
{
    // Configuration
    const Time& time_;
    
    struct criterionData
    {
        word name;
        scalar toleranceRSD;
        scalar toleranceRel;
        label momentComponent;
    
        DynamicList<scalar> bufferAvg;
        DynamicList<scalar> bufferConv;
        autoPtr<OFstream> filePtr;
        bool converged;
    
        criterionData()
        :
            name(""),
            toleranceRSD(0),
            toleranceRel(0),
            momentComponent(2),
            bufferAvg(),
            bufferConv(),
	    filePtr(nullptr),
            converged(false)
        {}
    };
    
    DynamicList<criterionData> criteria_;

    label bufferSizeAvg_;
    label bufferSizeConv_;
    label minTimeSteps_;
    label stepCounter_;

    // Retrieve monitored properties
    const dictionary& getProperties
    (
        const word& functionObjectName
    ) const;

    // Evaluate convergence for one quantity
    bool checkConvergence
    (
        const scalar currentValue,
        const scalar toleranceRSD,
        const scalar toleranceRel,
        DynamicList<scalar>& bufferAvg,
        DynamicList<scalar>& bufferConv,
        const word& quantityName,
        autoPtr<OFstream>& filePtr
    );

public:

    // Runtime type information
    TypeName("convergenceCheck");

    // Constructor
    convergenceCheck
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict
    );

    // Destructor
    virtual ~convergenceCheck() = default;

    //// Read configuration
    //virtual bool read(const dictionary& dict);
    virtual bool write() override;

    // Execute at the configured simulation interval
    virtual bool execute();
};

} // End namespace Foam

#endif

