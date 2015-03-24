#include "EventManager.h"
#include <stdexcept>
#include "TMath.h"
#include "Rtypes.h"
#ifdef hasPluto
#include "PParticle.h"
#include "PStaticData.h"
#endif

#include "GTreeTrack.h"
#include "Detector.h"
#include "GTreeTagger.h"
#include "GTreePluto.h"

#include <string>

using namespace std;
using namespace ant;

EventManager::EventManager(): maxevents(0), pluto_database(nullptr)
{
#ifdef hasPluto
    pluto_database = makeStaticData();
#endif
}

EventManager::~EventManager()
{
}

Bool_t EventManager::Init(const char *configfile)
{
    // nothing to do...
    return true;
}

void EventManager::Finish()
{
    for( auto& p : physics ) {
        p->Finish();
    }
}

Bool_t EventManager::Start()
{
    SetAsPhysicsFile();

    if(maxevents==0)
        TraverseValidEvents();
    else
        TraverseEntries(0,maxevents);

    return kTRUE;
}

void EventManager::ProcessEvent()
{
    checkMCIDs();

    Event event;

    CopyTracks(GetTracks(), event.TrackStorage());

    CopyParticles(GetPhotons(),         ParticleTypeDatabase::Photon,       event);
    CopyParticles(GetProtons(),         ParticleTypeDatabase::Proton,       event);
    CopyParticles(GetChargedPions(),    ParticleTypeDatabase::PiCharged,    event);
    CopyParticles(GetElectrons(),       ParticleTypeDatabase::eCharged,     event);

    CopyTaggerHits(event.TaggerHitStorage());

#ifdef hasPluto
    CopyPlutoParticles(GetPluto(), event);
#endif

    CopyTriggerInfo(GetTrigger(),    event);

    event.Finalize();

    RunPhysics(event);
}

void EventManager::ProcessScalerRead()
{
}

Bool_t EventManager::Write()
{
    return true;
}

void EventManager::RunPhysics(const Event &event)
{
    for( auto& p : physics ) {
        p->ProcessEvent(event);
    }
}

void EventManager::CopyParticles(GTreeParticle *tree, const ParticleTypeDatabase::Type &type, Event &target)
{
    for(Int_t i=0; i<tree->GetNParticles(); ++i) {

        const TLorentzVector& lv = tree->Particle(i);
        const Int_t trackIndex = tree->GetTrackIndex(i);
        const ant::Track* track = &target.TrackStorage().at(trackIndex);

        target.ParticleStorage().emplace_back(
                    RecParticle(
                        Particle(type,lv),
                        track)
                    );
    }
}


/**
 * @brief map goat apparatus numbers to apparatus_t enum values
 * in case unknown values show up: -> exception and do not sliently ignore
 */
detector_t IntToDetector_t( const int& a ) {
    detector_t d = detector_t::None;
    if(a & GTreeTrack::DETECTOR_NaI) {
        d |= detector_t::NaI;
    }
    if(a & GTreeTrack::DETECTOR_PID) {
        d |= detector_t::PID;
    }
    if(a & GTreeTrack::DETECTOR_MWPC) {
        d |= detector_t::MWPC;
    }
    if(a & GTreeTrack::DETECTOR_BaF2) {
        d |= detector_t::BaF2;
    }
    if(a & GTreeTrack::DETECTOR_PbWO4) {
        d |= detector_t::PbWO4;
    }
    if(a & GTreeTrack::DETECTOR_Veto) {
        d |= detector_t::Veto;
    }
    return d;
}

/**
 * @brief map the cluster sizes from goat to unisgend ints
 * negative values mean no hit in the calorimeter
 * map those to 0
 */
clustersize_t MapClusterSize(const int& size) {
    return size < 0 ? 0 : size;
}

void EventManager::CopyTracks(GTreeTrack *tree, Event::TrackList_t &container)
{
    for(Int_t i=0; i<tree->GetNTracks(); ++i) {

        container.emplace_back(
                    tree->GetClusterEnergy(i),
                    tree->GetTheta(i) * TMath::DegToRad(),
                    tree->GetPhi(i) * TMath::DegToRad(),
                    tree->GetTime(i),
                    MapClusterSize(tree->GetClusterSize(i)),
                    IntToDetector_t(tree->GetDetectors(i)),
                    tree->GetVetoEnergy(i),
                    tree->GetMWPC0Energy(i),
                    tree->GetMWPC1Energy(i)
                    );
    }
}

#ifdef hasPluto
void EventManager::CopyPlutoParticles(GTreePluto *tree, Event& event)
{
    const GTreePluto::ParticleList particles = tree->GetAllParticles();
    event.MCTrueStorage().reserve(particles.size());

    const ParticleTypeDatabase::Type* type=nullptr;
    for( auto& p : particles ) {
        type = ParticleTypeDatabase::GetTypeOfPlutoID(p->ID());
        if(!type) {
            type = ParticleTypeDatabase::AddTempPlutoType(
                        p->ID(),
                        "Pluto_"+to_string(p->ID()),
                        "Pluto_"+ string( pluto_database->GetParticleName(p->ID()) ),
                        pluto_database->GetParticleMass(p->ID())*1000.0,
                        pluto_database->GetParticleCharge(p->ID()) != 0
                    );
            if(!type)
                throw std::out_of_range("Could not create dynamic mapping for Pluto Particle ID "+to_string(p->ID()));
        }

        TLorentzVector lv = *p;
        lv *= 1000.0;   // convert to MeV

        event.MCTrueStorage().emplace_back( MCParticle(
                                                *type,
                                                lv,
                                                (p->GetDaughterIndex()==-1) )
                                            );
    }
}
#endif

void EventManager::CopyTaggerHits(Event::TaggerHitList_t &container)
{
    const GTreeTagger& tagger = *GetTagger();

    for( Int_t i=0; i<tagger.GetNTagged(); ++i) {
        container.emplace_back(
                    TaggerHit(
                        tagger.GetTaggedChannel(i),
                        tagger.GetTaggedEnergy(i),
                        tagger.GetTaggedTime(i))
                    );
    }
}

void EventManager::CopyTriggerInfo(GTreeTrigger *tree, Event &event)
{
    TriggerInfo& ti = event.Trigger();

    ti.CBEenergySum() = tree->GetEnergySum();
    ti.Multiplicity() = tree->GetMultiplicity();

    for( int err=0; err < tree->GetNErrors(); ++err) {
        ti.Errors().emplace_back(
                    tree->GetErrorModuleID()[err],
                    tree->GetErrorModuleIndex()[err],
                    tree->GetErrorCode()[err]);
    }
}

void EventManager::checkMCIDs() {
#ifdef hasPluto
    const bool eventIDmatch = ( GetTrigger()->GetMCTrueEventID() == GetPluto()->GetPlutoID()) || (GetPluto()->GetPlutoID() == -1);
    const bool randomIDmatch = ( GetTrigger()->GetMCRandomID() == GetPluto()->GetPlutoRandomID()) || (GetPluto()->GetPlutoRandomID() == -1);

    if(! (eventIDmatch && randomIDmatch) )
        throw data_check_exception(
                "MC Event ID does not match! Acqu:" +
                to_string(GetTrigger()->GetMCTrueEventID()) + "/" + to_string(GetTrigger()->GetMCRandomID())
                + "vs Pluto:"
                + to_string(GetPluto()->GetPlutoID()) + "/" + to_string(GetPluto()->GetPlutoRandomID()));
#else
    return;
#endif
}
