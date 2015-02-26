#ifndef __GTree_h__
#define __GTree_h__


#include <iostream>
#include <stdlib.h>

#include <TObject.h>
#include <TFile.h>
#include <TTree.h>



class   GTreeManager;

class  GTree    : protected TObject
{
public:
    enum
    {
        FLAG_CLOSED           = 0,
        FLAG_OPENFORINPUT     = 1,
        FLAG_OPENFOROUTPUT    = 2
    };

    friend  class   GTreeManager;

private:
    TString         name;
    Bool_t          correlatedToScalerRead;
    Bool_t          singleRead;
    Int_t           status;

    void    GetEntryFast(const UInt_t index)    {inputTree->GetEntry(index);}

protected:
    TTree*          inputTree;
    TTree*          outputTree;
    GTreeManager*   manager;
    bool            saveToFile;

    virtual void    SetBranchAdresses() = 0;
    virtual void    SetBranches() = 0;

public:
    GTree(GTreeManager* Manager, const TString& _Name, const Bool_t CorrelatedToScalerRead = kFALSE, const Bool_t SingleRead = kFALSE);
    virtual ~GTree();

    virtual void        Clear() = 0;
            void        CloseForInput();
            void        CloseForOutput();
            void        Close();
            void        Fill();
    inline  Bool_t      GetEntry(const UInt_t index);
    virtual void        Unpack() {}
    const   char*       GetName() const {return name.Data();}
            UInt_t      GetNEntries()   { if(IsOpenForInput()) return inputTree->GetEntries(); return 0;}
            Bool_t      IsClosed()          {return !status;}
            Bool_t      IsOpenForInput()    {return status & FLAG_OPENFORINPUT;}
            Bool_t      IsOpenForOutput()   {return status & FLAG_OPENFOROUTPUT;}
            Bool_t      OpenForInput();
            Bool_t      OpenForOutput();
    virtual void        Print() const;
    virtual Bool_t      Write();

};

Bool_t  GTree::GetEntry(const UInt_t index)
{
    if(index > inputTree->GetEntries())
        return kFALSE;
    inputTree->GetEntry(index);
}


#endif
