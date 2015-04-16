#ifndef PION_H
#define PION_H

#include "AntPhysics.h"
#include "base/interval.h"
#include "plot/SmartHist.h"

class TH1D;
class TH2D;

namespace ant {
namespace analysis {

class Pion: public Physics {
protected:

    SmartHist<const TLorentzVector&> p_MM;

    SmartHist<int> pion_rec_multi;

    SmartHist<int> nr_ngamma;
    SmartHist<const TLorentzVector&> nr_4gim;
    SmartHist<const TLorentzVector&> nr_3gim;
    SmartHist<const TLorentzVector&> nr_2gim;

    int n=0;
    SmartHist<const ParticlePtr&> test;

    IntervalD pi0_im_cut;
    IntervalD tagger_energy_cut;

    TLorentzVector target;

    SmartHist<const std::string&> step_levels;

    SmartHist< std::pair<const TLorentzVector&,const TLorentzVector&> > pion_mc_rec_angle;

    SmartHist<const TLorentzVector&> makeInvMassPlot(const std::string& title, const std::string& xlabel, const std::string& ylabel, ant::BinSettings bins, const std::string& name="");
    SmartHist< std::pair<const TLorentzVector&, const TLorentzVector&> > makeAngleDiffPlot(const std::string& title, const std::string& xlabel, const std::string& ylabel, BinSettings bins, const std::string& name);

public:
    Pion(const std::string& name="Pion", const mev_t energy_scale=1000.0);
    virtual ~Pion() {}
    void ProcessEvent(const Event &event);
    void Finish();
    void ShowResult();
};

}
}
#endif
