#include "GoAT.h"


GoAT::GoAT() :
    useParticleReconstruction(0),
    nEventsWritten(0)
{ 
}

GoAT::~GoAT()
{
}

Bool_t	GoAT::Init(const char* configFile)
{
    cout << endl << "Initialising GoAT analysis..." << endl << endl;

    if(configFile)
        SetConfigFile(configFile);
    std::string config = ReadConfig("Period-Macro");
    if( sscanf(config.c_str(),"%d\n", &period) == 1 ) usePeriodMacro = 1;

	cout << "==========================================================" << endl;	
	cout << "Setting up Data Checks:" << endl;	
	cout << "==========================================================" << endl;	
    if(!GDataChecks::Init())
	{
		cout << "GDataChecks Init failed!" << endl; 
		return kFALSE;
	}

	cout << "==========================================================" << endl;	
	cout << "Setting up sorting criteria:" << endl;	
	cout << "==========================================================" << endl;	
    if(!GSort::Init())
	{
		cout << "GSort Init failed!" << endl; 
		return kFALSE;
	}
	
	cout << "==========================================================" << endl;		
	cout << "Setting up analysis classes:" << endl;	
	cout << "==========================================================" << endl;	
    config = ReadConfig("DO-PARTICLE-RECONSTRUCTION");
    if (strcmp(config.c_str(), "nokey") != 0)
    {
        int buffer=0;
        sscanf( config.c_str(), "%d\n", &buffer);
        useParticleReconstruction = (buffer==1);
    }

    if(useParticleReconstruction)
	{
        if(!GParticleReconstruction::Init())
		{
			cout << "GParticleReconstruction Init failed!" << endl; 
			return kFALSE;
		}
    }

    config = ReadConfig("DO-MESON-RECONSTRUCTION");
    if (strcmp(config.c_str(), "nokey") != 0)
    {
        int buffer=0;
        sscanf( config.c_str(), "%d\n", &buffer);
        useMesonReconstruction = (buffer==1);
    }

    if(useMesonReconstruction)
    {
        if(!GMesonReconstruction::Init())
        {
            cout << "GParticleReconstruction Init failed!" << endl;
            return kFALSE;
        }
    }

	cout << endl;	

	cout << "Initialisation complete." << endl;
	cout << "==========================================================" << endl << endl;
  
	return kTRUE;
}

void	GoAT::ProcessEvent()
{
    if(usePeriodMacro == 1)
    {
        if(GetEventNumber() % period == 0)
            cout << "Event: " << GetEventNumber() << "  Events Accepted: " << nEventsWritten << endl;
    }

    if(SortAnalyseEvent())
    {
        if(useParticleReconstruction)
        {
            if(useMesonReconstruction)
            {
                if(!GParticleReconstruction::ProcessEventWithoutFilling())  return;
                if(!GMesonReconstruction::ProcessEventWithoutFilling())  return;
                if(!SortFillEvent())    return;
                electrons->Fill();
                protons->Fill();
                neutrons->Fill();
                pi0->Fill();
                eta->Fill();
                etap->Fill();
            }
            else
            {
                if(!GParticleReconstruction::ProcessEventWithoutFilling())  return;
                if(!SortFillEvent())    return;
                electrons->Fill();
                protons->Fill();
                neutrons->Fill();
            }
        }
        else if(useMesonReconstruction)
        {
            GMesonReconstruction::ProcessEventWithoutFilling();
            if(!SortFillEvent())    return;
            pi0->Fill();
            eta->Fill();
            etap->Fill();
        }
        eventParameters->SetNReconstructed(GetNReconstructed());
        eventParameters->Fill();
		rootinos->Fill();
        photons->Fill();
        chargedPions->Fill();
        FillReadList();
        nEventsWritten++;
    }
}

Bool_t	GoAT::Start()
{
    if(!IsAcquFile())
    {
        cout << "ERROR: Input File is not a Acqu file." << endl;
        return kFALSE;
    }
    SetAsGoATFile();

    if(useParticleReconstruction)
    {
        if(useMesonReconstruction)
        {
			rootinos->CloseForInput();
            photons->CloseForInput();
            electrons->CloseForInput();
            chargedPions->CloseForInput();
            protons->CloseForInput();
            neutrons->CloseForInput();
            pi0->CloseForInput();
            eta->CloseForInput();
            etap->CloseForInput();
        }
        else
        {
			rootinos->CloseForInput();
            photons->CloseForInput();
            electrons->CloseForInput();
            chargedPions->CloseForInput();
            protons->CloseForInput();
            neutrons->CloseForInput();
        }
    }
    else if(useMesonReconstruction)
    {
        pi0->CloseForInput();
        eta->CloseForInput();
        etap->CloseForInput();
    }

    if(!TraverseValidEvents())		return kFALSE;

    return kTRUE;
}
