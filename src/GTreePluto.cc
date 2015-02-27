#include "GTreePluto.h"
#include "PParticle.h"
#include <iostream>

using namespace std;

GTreePluto::GTreePluto(GTreeManager* Manager):
    GTree(Manager, TString("data")),
    PlutoMCTrue(NULL),
    plutoID(-1),
    plutoRandomID(-1)
{}


GTreePluto::~GTreePluto()
{}

bool IsFinalState(const PParticle* particle) { return particle->GetDaughterIndex() == -1; }
bool IsBeamParticle(const PParticle* particle) { return particle->ID() > 1000; }

void GTreePluto::Unpack()
{
    finalstate.clear();
    allparticles.clear();
    beamparticles.clear();

    for( Long64_t i=0;i<PlutoMCTrue->GetEntries(); ++i) {

        const PParticle* const particle = dynamic_cast<const PParticle*>((*PlutoMCTrue)[i]);

        if(IsBeamParticle(particle))
            beamparticles.emplace_back(particle);
        else {
            allparticles.emplace_back(particle);

            if(IsFinalState(particle))
                finalstate.emplace_back(particle);
        }
    }
}

PParticle* GTreePluto::GetMCTrue(const int idx) const throw(std::exception)
{
    if (!PlutoMCTrue)
        //return nullptr;
        throw noTree;

    return dynamic_cast<PParticle*>(PlutoMCTrue->At(idx));
}

TLorentzVector GTreePluto::GetTrueP4(const int idx) const throw(std::exception)
{
    if (!PlutoMCTrue)
        //return TLorentzVector(0., 0., 0., 0.);
        throw noTree;

    return dynamic_cast<PParticle*>(PlutoMCTrue->At(idx))->Vect4();
}

// Assume proton target if no argument is given
TLorentzVector GTreePluto::GetTrueBeam() const
{
    return GetTrueBeam(TLorentzVector(0., 0., 0., .938272), 14001);
}

TLorentzVector GTreePluto::GetTrueBeam(const TLorentzVector target, const int beamID) const
{

    auto it = std::find_if(beamparticles.begin(), beamparticles.end(), [beamID](PParticle const& p){ return p.ID() == beamID; });
    if (it == beamparticles.end())
        return TLorentzVector(0., 0., 0., 0.);

    return (*it)->Vect4() - target;
}

void    GTreePluto::SetBranchAdresses()
{
    inputTree->SetBranchAddress("Particles", 	&PlutoMCTrue);
    inputTree->SetBranchAddress("plutoID", 	&plutoID);
    inputTree->SetBranchAddress("plutoRandomID", 	&plutoRandomID);
}

void    GTreePluto::SetBranches()
{
    outputTree->Branch("Particles", PlutoMCTrue);
    outputTree->Branch("plutoID", &plutoID);
    outputTree->Branch("plutoRandomID", &plutoRandomID);
}
