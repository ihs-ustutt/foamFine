
#include "convergenceCheck.H"

#include "addToRunTimeSelectionTable.H"
#include "IOdictionary.H"
#include "Pstream.H"
#include "OSspecific.H"
#include "IOobject.H"
#include "mathematicalConstants.H"

namespace Foam
{

defineTypeNameAndDebug(convergenceCheck, 0);

addToRunTimeSelectionTable
(
    functionObject,
    convergenceCheck,
    dictionary
);


// * * * * * * * * * Constructor * * * * * * * * * //

convergenceCheck::convergenceCheck
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObject(name, false),
    time_(runTime),

    bufferSizeAvg_(dict.lookupOrDefault<label>("bufferSizeAvg", 100)),
    bufferSizeConv_(dict.lookupOrDefault<label>("bufferSizeConv", 1)),
    minTimeSteps_(dict.lookupOrDefault<label>("minTimeSteps", 200)),
    stepCounter_(0)
{
    const dictionary& criteriaDict = dict.subDict("criteria");
    
    forAllConstIter(dictionary, criteriaDict, iter)
    {
        const word criterionName = iter().keyword(); 
        const dictionary& criterionDict = iter().dict();
        
	criteria_.setSize(criteria_.size() + 1);
        criterionData& c = criteria_.last();
        //criterionData c;
    
        c.name = criterionName;
    
        c.toleranceRSD =
            criterionDict.get<scalar>("toleranceRSD");
    
        c.toleranceRel =
            criterionDict.get<scalar>("toleranceRel");
    
        // Optional parameter for the moment criterion
        if (criterionName.substr(0,2) == "F_")
        {
            c.momentComponent =
                criterionDict.get<label>("moment");
        }
     
        //criteria_.append(c);
    }
    if (not(criteria_.size() > 0))
    {
        FatalErrorInFunction
            << "no convergence criteria set"
            << exit(FatalError);
    }
    if (bufferSizeAvg_ < 1 or bufferSizeConv_ < 1)
    {
        FatalErrorInFunction
            << "bufferSize must be positive"
            << exit(FatalError);
    }
}

// * * * * * * * * * Retrieve properties * * * * * * * //

const dictionary& convergenceCheck::getProperties
(
    const word& functionObjectName
) const
{
    const IOdictionary& globalProps =
        time_.lookupObject<IOdictionary>
        (
            "functionObjectProperties"
        );

    forAllConstIter(dictionary, globalProps, iter)
    {
        if (!iter().isDict())
        {
            continue;
        }

        const dictionary& group = iter().dict();

        if (group.found(functionObjectName))
        {
            const dictionary& props =
                group.subDict(functionObjectName);

            if (props.found("vector"))
            {
                return props.subDict("vector");
            }

            if (props.found("scalar"))
            {
                return props.subDict("scalar");
            }

            FatalErrorInFunction
                << "Neither 'vector' nor 'scalar' found for "
                << functionObjectName
                << exit(FatalError);
        }
    }

    FatalErrorInFunction
        << "Cannot find " << functionObjectName
        << " in functionObjectProperties"
        << exit(FatalError);

    return globalProps;
}


// * * * * * * * * * Convergence calculation * * * * * //

bool convergenceCheck::checkConvergence
(
    const scalar currentValue,
    const scalar toleranceRSD,
    const scalar toleranceRel,
    DynamicList<scalar>& buffer,
    DynamicList<scalar>& bufferConv,
    const word& quantityName,
    autoPtr<OFstream>& filePtr
)
{
    buffer.append(currentValue);

    if (buffer.size() > bufferSizeAvg_)
    {
        buffer.remove(0);
    }

    scalar meanValue = 0.0;
    scalar meanAbsValue = 0.0;

    scalar minValue = buffer[0];
    scalar maxValue = buffer[0];

    forAll(buffer, i)
    {
        meanValue += buffer[i];
        meanAbsValue += mag(buffer[i]);

        minValue = min(minValue, buffer[i]);
        maxValue = max(maxValue, buffer[i]);
    }

    meanValue /= buffer.size();
    meanAbsValue /= buffer.size();

    const scalar relativeRange =
        (maxValue - minValue)/max(meanAbsValue, SMALL);

    scalar variance = 0.0;

    forAll(buffer, i)
    {
        const scalar difference = buffer[i] - meanValue;
        variance += sqr(difference);
    }

    variance /= buffer.size();

    const scalar standardDeviation = sqrt(variance);

    const scalar relativeStdDev =
        standardDeviation/max(mag(meanValue), SMALL);

    const bool convergedRSD =
        relativeStdDev < toleranceRSD;

    const bool convergedRel =
        relativeRange < toleranceRel;

    const bool converged =
        convergedRSD && convergedRel;
    
    if (converged)
    {
        bufferConv.append(converged);
    }
    else
    {
	bufferConv.clear();
    }

    Info<< quantityName << ':' << nl
	<< "  Current value = " << currentValue << nl
        << "  Relative standard deviation = "
        << relativeStdDev*100 << '%' << nl
        << "  Relative range = "
        << relativeRange*100 << '%' << nl
        << "  Converged = " << converged << nl
        << "  Number of converged time steps: "
        << bufferConv.size() << '/'
        << bufferSizeConv_ << nl
        << endl;

    if (Pstream::master())
    {
        const fileName outputDir =
            time_.globalPath()
           /"postProcessing"
           /"convergence";

        mkDir(outputDir);

        if (!filePtr.valid())
        {
            const word baseName =
                quantityName + "_convergence";

            fileName outputFile =
                outputDir/(baseName + ".dat");

            label fileIndex = 1;

            while (isFile(outputFile))
            {
                outputFile =
                    outputDir
                   /(baseName + "." + Foam::name(fileIndex) + ".dat");

                ++fileIndex;
            }

            filePtr.reset(new OFstream(outputFile));

            filePtr()
                << "# Time RSD relRange" << nl;

            Info<< "Writing convergence data to "
                << outputFile << endl;
        }

        filePtr()
            << time_.value() << ' '
            << relativeStdDev << ' '
            << relativeRange << nl;

        filePtr().flush();
    }

    if (stepCounter_ < minTimeSteps_)
    {
        Info<< "Minimum time steps not yet reached: "
            << stepCounter_ << '/'
            << minTimeSteps_ << endl;

        return false;
    }
    if (buffer.size() < bufferSizeAvg_)
    {
        Info<< quantityName
            << ": collecting values for averaging ("
            << buffer.size() << '/'
            << bufferSizeAvg_ << ')'
            << endl;

        return false;
    }
    if (bufferConv.size() < bufferSizeConv_)
    {
	Info<< "Number of converged time steps not yet reached" << endl;
        return false;
    }

    return converged;
}


// * * * * * * * * * Execute * * * * * * * * * * * * //

bool convergenceCheck::execute()
{
    
    Info<< "\n======================================================================" << endl;
    Info<< "STAUSSPUT - Convergence Criteria" << endl;
    Info<< "======================================================================\n" << endl;

    ++stepCounter_;
     
    forAll(criteria_, c)
    {
        scalar value = 0;
        word label;
	
	if (criteria_[c].name.substr(0,2) == "F_")
	{
    	    // Moment convergence
    	    const dictionary& forceProps =
    	        getProperties(criteria_[c].name);
	    
	    word comp;

    	    vector pressureMoment;
    	    vector viscousMoment;

    	    forceProps.lookup("pressureMoment") >> pressureMoment;
    	    forceProps.lookup("viscousMoment") >> viscousMoment;
	    
	    if (criteria_[c].momentComponent == 0)
	    {
    	    	value = pressureMoment.x() + viscousMoment.x();
		comp = "X";
	    }
	    else if (criteria_[c].momentComponent == 1)
	    {
    	    	value = pressureMoment.y() + viscousMoment.y();
		comp = "Y";
	    }
	    else if (criteria_[c].momentComponent == 2)
	    {
    	    	value = pressureMoment.z() + viscousMoment.z();
		comp = "Z";
	    }
	    else
	    {		
                FatalErrorInFunction
                    << "Invalid moment component option for convergence criteria " 
		    << criteria_[c].name << nl
                    << "Valid options for moment are 0, 1 or 2"
                    << exit(FatalError);
	    }

    	    label = criteria_[c].name + "_Moment" + comp;
	}
	else if (criteria_[c].name.substr(0,3) == "PT_")
	{
    	    // Moment convergence
    	    const dictionary& pressureProps =
    	        getProperties(criteria_[c].name);

            value =
                pressureProps.get<scalar>
                (
                    "weightedAverage("+criteria_[c].name.substr(3)+",ptot)"
                );
    	    
	    label = criteria_[c].name;
	     
	}
	else if (criteria_[c].name.substr(0,2) == "Q_")
	{
    	    // Moment convergence
    	    const dictionary& volProps =
    	        getProperties(criteria_[c].name);

            value =
                volProps.get<scalar>
                (
            	    "sum("+criteria_[c].name.substr(2)+",phi)"
                );
    	    
	    label = criteria_[c].name;
	     
	}
	else
        {		
            FatalErrorInFunction
                << "Invalid pyhsical value for convergence criteria " 
                << criteria_[c].name << nl
                << "Valid value denominators are: " << nl
		<< "F  -> Moment" << nl 
		<< "PT -> Total Pressure" << nl
		<< "Q  -> Discharge"
                << exit(FatalError);
        }

    	criteria_[c].converged =
	    checkConvergence
	    (
	        value,
	        criteria_[c].toleranceRSD,
	        criteria_[c].toleranceRel,
	        criteria_[c].bufferAvg,
	        criteria_[c].bufferConv,
	        label,
	        criteria_[c].filePtr
	    );	
        Info<< "----------------------------------------------------------------------" << endl;
    }

    bool allConverged = true;

    forAll(criteria_, c)
    {
        if (!criteria_[c].converged)
        {
            allConverged = false;
            break;
        }
    }
    
    if (allConverged)
    {
        Info<< "All convergence criteria satisfied. "
            << "Stopping simulation." << endl;
    
        const_cast<Time&>(time_).writeAndEnd();
    }
     
    return true;
}
bool Foam::convergenceCheck::write()
{
    return true;
}
} // End namespace Foam

