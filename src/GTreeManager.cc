#include "GTreeManager.h"

#include <TROOT.h>
#include <TSystemFile.h>
#include <TSystemDirectory.h>
#include <TFileCacheWrite.h>

using namespace std;


GTreeManager::GTreeManager()    :
    GHistManager(),
    GConfigFile(),
    inputFile(0),
    outputFile(0),
    treeList(),
    treeCorreleatedToScalerReadList(),
    readList(),
    readCorreleatedToScalerReadList(),
    writeList(),
    countReconstructed(0),
    tracks(0),
    tagger(0),
    trigger(0),
    scalers(0),
    detectorHits(0),
    setupParameters(0),
    rootinos(0),
    photons(0),
    electrons(0),
    chargedPions(0),
    protons(0),
    neutrons(0),
    neutralPions(0),
    etas(0),
    etaPrimes(0),
#ifdef hasPluto
    linpol(0),
    pluto(NULL),
#else
    linpol(0),
#endif
    geant(NULL)

{
    pdgDB = TDatabasePDG::Instance();

    rootinos = new GTreeParticle(this, "rootinos");
    photons = new GTreeParticle(this, "photons");
    electrons = new GTreeParticle(this, "electrons");
    chargedPions = new GTreeParticle(this, "chargedPions");
    protons = new GTreeParticle(this, "protons");
    neutrons = new GTreeParticle(this, "neutrons");

    neutralPions = new GTreeMeson(this, "neutralPions");
    etas = new GTreeMeson(this, "etas");
    etaPrimes = new GTreeMeson(this, "etaPrimes");

    setupParameters = new GTreeSetupParameters(this);
    eventParameters = new GTreeEventParameters(this);
    detectorHits = new GTreeDetectorHits(this);
    tracks = new GTreeTrack(this, "tracks");
    tagger = new GTreeTagger(this);
    trigger = new GTreeTrigger(this);
    scalers = new GTreeScaler(this);
    linpol = new GTreeLinPol(this);
#ifdef hasPluto
    pluto = new GTreePluto(this);
#endif
    geant = new GTreeA2Geant(this);
}

GTreeManager::~GTreeManager()
{
    while(treeList.GetEntries()>0)
    {
        if((GTree*)treeList[0])
            delete (GTree*)treeList[0];
    }
}

Bool_t  GTreeManager::TraverseEntries(const UInt_t min, const UInt_t max)
{
    if(!inputFile)
        return kFALSE;

    for(UInt_t i=min; i<max; i++)
    {

        GetAndUnpack(i);

        eventParameters->SetEventNumber(i);
        countReconstructed = 0;
        ProcessEvent();
    }

    return kTRUE;
}

Bool_t  GTreeManager::TraverseScalerEntries(const UInt_t min, const UInt_t max)
{
    if(!inputFile)
        return kFALSE;

    if(!scalers)
        return kFALSE;
    if(!scalers->IsOpenForInput())
    {
        if(!scalers->OpenForInput())
        {
            std::cout << "Can not open Scaler tree in input file." << std::endl;
            return kFALSE;
        }
    }

    for(UInt_t i=min; i<=max; i++)
    {
        scalers->GetEntryFast(i);
        countReconstructed = 0;
        ProcessEvent();
    }

    return kTRUE;
}


Bool_t  GTreeManager::TraverseFiles()
{
    Int_t nFiles = GetNFiles();
    for(Int_t i=0; i<nFiles; i++)
    {
        std::string inputFileName = GetInputFile(i);
        std::string outputFileName = GetOutputFile(i);
        if(!StartFile(inputFileName.c_str(), outputFileName.c_str())) cout << "ERROR: Failed on file " << inputFileName << "!" << endl;
    }

    return kTRUE;
}

Bool_t  GTreeManager::StartFile(const char* inputFileName, const char* outputFileName)
{
    if(inputFile)    inputFile->Close();
    inputFile = TFile::Open(inputFileName);
    if(!inputFile)
    {
        cout << "#ERROR: Can not open input file " << inputFileName << "!" << endl;
        return kFALSE;
    }
    cout << "Opened input file " << inputFile->GetName() << "!" << inputFile->GetTitle() << endl;

    for(Int_t l=0; l<treeList.GetEntries(); l++)
    {
        if(inputFile->Get(((GTree*)treeList[l])->GetName()))
            ((GTree*)treeList[l])->OpenForInput();
    }
    for(Int_t l=0; l<treeCorreleatedToScalerReadList.GetEntries(); l++)
    {
        if(inputFile->Get(((GTree*)treeCorreleatedToScalerReadList[l])->GetName()))
            ((GTree*)treeCorreleatedToScalerReadList[l])->OpenForInput();
    }
    for(Int_t l=0; l<treeSingleReadList.GetEntries(); l++)
    {
        if(inputFile->Get(((GTree*)treeSingleReadList[l])->GetName()))
            ((GTree*)treeSingleReadList[l])->OpenForInput();
    }

    if(outputFile)
        outputFile->Close();
    outputFile = TFile::Open(outputFileName, "RECREATE");
    if(!outputFile)
    {
        cout << "#ERROR: Can not create output file " << outputFileName << "!" << endl;
        return kFALSE;
    }
    cout << "Created output file " << outputFile->GetName() << "!" << outputFile->GetTitle() << endl;
    TFileCacheWrite*    cache   = new TFileCacheWrite(outputFile, 104857600);

    isWritten   = kFALSE;
    ClearLinkedHistograms();

    if(!Start())
        return kFALSE;

    if(!isWritten)
        Write();
    cache->Flush();

    for(Int_t l=0; l<treeList.GetEntries(); l++)
        ((GTree*)treeList[l])->Close();
    for(Int_t l=0; l<treeCorreleatedToScalerReadList.GetEntries(); l++)
        ((GTree*)treeCorreleatedToScalerReadList[l])->Close();
    for(Int_t l=0; l<treeSingleReadList.GetEntries(); l++)
        ((GTree*)treeSingleReadList[l])->Close();


    if(inputFile)     inputFile->Close();
    if(outputFile)    outputFile->Close();
    //delete  cache;

    return kTRUE;
}

Bool_t  GTreeManager::Write()
{
    if(!outputFile)   return kFALSE;
    outputFile->cd();

    for(Int_t l=0; l<writeList.GetEntries(); l++)
        ((GTree*)writeList[l])->Write();

    TIter objectIterator(gROOT->GetList());
    TObject *object;
    while((object=(TObject*)objectIterator()))
    {
        object->Write();
        std::cout << "object " << object->GetName() << " has been written to disk." << std::endl;
    }

    WriteLinkedHistograms(outputFile);

    isWritten   = kTRUE;

    return kTRUE;
}

Bool_t  GTreeManager::Write(const TNamed* object)
{
    if(!outputFile)   return kFALSE;
    outputFile->cd();
    object->Write();
    std::cout << "object " << object->GetName() << " has been written to disk." << std::endl;
    return kTRUE;
}

void GTreeManager::GetAndUnpack(const UInt_t event)
{
    // read from all trees
    for(Int_t l=0; l<readList.GetEntriesFast(); l++)
        ((GTree*)readList[l])->GetEntryFast(event);

    // unpack event in all trees
    for(Int_t l=0; l<readList.GetEntriesFast(); l++)
        ((GTree*)readList[l])->Unpack();
}


Bool_t  GTreeManager::TraverseValidEvents_AcquTreeFile()
{
    for(Int_t l=0; l<readSingleReadList.GetEntriesFast(); l++)
    {
        ((GTree*)readSingleReadList[l])->GetEntryFast(0);
        ((GTree*)readSingleReadList[l])->Fill();
        if(!tagger->HasEnergy()) tagger->SetCalibration(setupParameters->GetNTagger(),setupParameters->GetTaggerPhotonEnergy());
    }

    if(!scalers->IsOpenForInput())
    {
        cout << "No Scaler tree available. Expect MC data. Loop over all events" << endl;
        cout << "\tProcess events from " << 0 << " to " << GetNEntries() << "."<< endl;
        TraverseEntries(0, GetNEntries());
        return true;
    }
    if(scalers->GetNEntries()==0)
    {
        cout << "No Scaler tree available. Expect MC data. Loop over all events" << endl;
        cout << "\tProcess events from " << 0 << " to " << GetNEntries() << "."<< endl;
        TraverseEntries(0, GetNEntries());
        return true;
    }

    if(scalers->GetNEntries()<2)
    {
        std::cout << "Less than 2 scaler reads. Can not find events with correct scalers" <<std::endl;
        return kFALSE;
    }

    // find correct shift
    Int_t shift;
    {
        Double_t shiftMean = 0;
        for(Int_t l=1; l<scalers->GetNEntries(); l++)
        {
            scalers->GetEntryFast(l);
            shiftMean    += scalers->GetEventNumber() - scalers->GetEventID();
        }
        shiftMean   /= scalers->GetNEntries()-1;
        Int_t bestIndex = 0;
        scalers->GetEntryFast(0);
        Double_t smallestDifference = shiftMean - (scalers->GetEventNumber() - scalers->GetEventID());
        for(Int_t l=1; l<scalers->GetNEntries(); l++)
        {
            scalers->GetEntryFast(l);
            if((shiftMean - (scalers->GetEventNumber() - scalers->GetEventID())) < smallestDifference)
            {
                bestIndex = l;
                smallestDifference  = shiftMean - (scalers->GetEventNumber() - scalers->GetEventID());
            }
        }
        scalers->GetEntryFast(bestIndex);
        shift   = scalers->GetEventNumber() - scalers->GetEventID();
    }

    outputFile->cd();
    TH1I*   accepted    = new TH1I("CountScalerValid", "Events with correct scalers (all=0,accepted=1,rejected=2)", 3, 0, 3);
    accepted->SetBinContent(1, GetNEntries());

    scalers->GetEntry(scalers->GetNEntries()-1);
    Int_t start = scalers->GetEventNumber();
    scalers->GetEntry(0);
    cout << "Checking scaler reads! Valid events from " << scalers->GetEventNumber() << " to " << start << endl;
    start = scalers->GetEventNumber();

    for(Int_t i=1; i<GetNScalerEntries(); i++)
    {
        for(Int_t l=0; l<readCorreleatedToScalerReadList.GetEntriesFast(); l++)
            ((GTree*)readCorreleatedToScalerReadList[l])->GetEntry(i);
        if(scalers->GetEventNumber() - scalers->GetEventID() == shift)
        {
            currentScalerEntry = i;
            accepted->SetBinContent(2, accepted->GetBinContent(2) + (scalers->GetEventNumber()-start));
            TraverseEntries(start, scalers->GetEventNumber());
            ProcessScalerRead();
            for(Int_t l=0; l<readCorreleatedToScalerReadList.GetEntriesFast(); l++)
                ((GTree*)readCorreleatedToScalerReadList[l])->Fill();
            start = scalers->GetEventNumber();
        }
    }

    accepted->SetBinContent(3, accepted->GetBinContent(1) - accepted->GetBinContent(2));

    if(!Write(accepted))  return kFALSE;

    if(accepted)    delete accepted;
    return kTRUE;
}

Bool_t  GTreeManager::TraverseValidEvents_GoATTreeFile()
{
    for(Int_t l=0; l<readSingleReadList.GetEntriesFast(); l++)
    {
        ((GTree*)readSingleReadList[l])->GetEntryFast(0);
        if(!tagger->HasEnergy()) tagger->SetCalibration(setupParameters->GetNTagger(),setupParameters->GetTaggerPhotonEnergy());
    }

    if(!scalers->IsOpenForInput())
    {
        cout << "No Scaler tree available. Expect MC data. Loop over all events" << endl;
        cout << "\tProcess events from " << 0 << " to " << GetNEntries() << "."<< endl;
        TraverseEntries(0, GetNEntries());
        return true;
    }
    if(scalers->GetNEntries()==0)
    {
        cout << "No Scaler tree available. Expect MC data. Loop over all events" << endl;
        cout << "\tProcess events from " << 0 << " to " << GetNEntries() << "."<< endl;
        TraverseEntries(0, GetNEntries());
        return true;
    }

    Int_t   event       = 0;
    Int_t   start       = 0;
    Int_t   maxEvent    = GetNEntries();
    for(Int_t l=0; l<readList.GetEntriesFast(); l++)
        ((GTree*)readList[l])->GetEntryFast(event);

    cout << GetNScalerEntries() << " scaler reads. " << maxEvent << " events." << endl;

    for(Int_t i=0; i<GetNScalerEntries(); i++)
    {
        for(Int_t l=0; l<readCorreleatedToScalerReadList.GetEntriesFast(); l++)
            ((GTree*)readCorreleatedToScalerReadList[l])->GetEntry(i);
        while(eventParameters->GetEventNumber()<scalers->GetEventNumber())
        {
            event++;
            if(event>=maxEvent)
                break;
            GetAndUnpack(event);
            ProcessEvent();
        }
        if(i!=0)
        {
            if(i%500 == 0)
            {
                cout << "\t" << i << " Scaler reads processed. Events from " << start << " to " << event << "." << endl;
                start   = event;
            }
        }
        ProcessScalerRead();
    }
    cout << "\t" << GetNScalerEntries() << " Scaler reads processed. Events from " << start << " to " << event << "." << endl;
}

UInt_t  GTreeManager::GetNEntries()       const
{
    for(Int_t l=1; l<readList.GetEntriesFast(); l++)
    {
        if(((GTree*)readList[l])->GetNEntries() != ((GTree*)readList[l-1])->GetNEntries())
        {
            cout << "ERROR: input trees have different number of entries!"<< endl;
            return 0;
        }
    }

    if(readList.GetEntriesFast() == 0)
        return 0;

    return  ((GTree*)readList[0])->GetNEntries();
}

UInt_t  GTreeManager::GetNScalerEntries()       const
{
    for(Int_t l=1; l<readCorreleatedToScalerReadList.GetEntriesFast(); l++)
    {
        if(((GTree*)readCorreleatedToScalerReadList[l])->GetNEntries() != ((GTree*)readCorreleatedToScalerReadList[l-1])->GetNEntries())
        {
            cout << "ERROR: input trees have different number of entries!"<< endl;
            return 0;
        }
    }

    if(readCorreleatedToScalerReadList.GetEntriesFast() == 0)
        return 0;

    return  ((GTree*)readCorreleatedToScalerReadList[0])->GetNEntries();
}

void    GTreeManager::SetAsGoATFile()
{
    if(!outputFile)
        return;
    TNamed flag("GoAT_File", "GoAT_File");
    Write(&flag);
}

void    GTreeManager::SetAsPhysicsFile()
{
    if(!outputFile)
        return;
    TNamed flag("Physics_File", "Physics_File");
    Write(&flag);
}

Bool_t  GTreeManager::IsAcquFile()    const
{
    if(IsGoATFile())
        return kFALSE;
    if(IsPhysicsFile())
        return kFALSE;
    return kTRUE;
}

Bool_t  GTreeManager::IsGoATFile()    const
{
    if(!inputFile)
        return kFALSE;
    TNamed* flag    = (TNamed*)inputFile->Get("GoAT_File");
    if(flag)
        return kTRUE;
    return kFALSE;
}

Bool_t  GTreeManager::IsPhysicsFile()    const
{
    if(!inputFile)
        return kFALSE;
    TNamed* flag    = (TNamed*)inputFile->Get("Physics_File");
    if(flag)
        return kTRUE;
    return kFALSE;
}
