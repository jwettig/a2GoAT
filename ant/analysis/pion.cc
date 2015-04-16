#include "pion.h"
#include "Particle.h"
#include "Event.h"
#include "TLorentzVector.h"
#include "TH1D.h"
#include "plot/root_draw.h"
#include "plot/Histogram.h"
#include "utils/combinatorics.h"
#include "TaggerHit.h"
#include <string>
#include <iostream>
#include "plot/SmartHist.h"

using namespace std;

ant::SmartHist<const TLorentzVector&> ant::analysis::Pion::makeInvMassPlot(const std::string& title, const std::string& xlabel, const std::string& ylabel, BinSettings bins, const std::string& name) {
    return HistFac.makeHist<const TLorentzVector&>(
                [] (const TLorentzVector& p) { return p.M();},
                title,
                xlabel, ylabel, bins, name);
}


ant::analysis::Pion::Pion(const string &name, const ant::mev_t energy_scale):
    Physics(name),
    tagger_energy_cut(80, 450),
    target(0.0, 0.0, 0.0, ParticleTypeDatabase::Proton.Mass())
{
    const BinSettings energy_bins(1000, 0.0, energy_scale);
    const BinSettings angle_diff_bins(200,0.0,20.0);

    pion_rec_multi = HistFac.makeHist<int>(
                ParticleTypeDatabase::Pi0.PrintName() + " Reconstruction Multiplicity",
                "n",
                "",
                BinSettings(5));

    nr_ngamma = HistFac.makeHist<int>(
                "Not reconstructed: number of photons",
                "number of photons/event",
                "",
                BinSettings(16));

    nr_2gim = HistFac.makeHist<const TLorentzVector&>([] (const TLorentzVector& v) { return v.M(); },
                "Not reconstructed: 2#gamma IM",
                "M_{2#gamma} [MeV]",
                "",
                energy_bins);

    nr_3gim = HistFac.makeHist<const TLorentzVector&>([] (const TLorentzVector& v) { return v.M(); },
                "Not reconstructed: 3#gamma IM",
                "M_{3#gamma} [MeV]",
                "",
                energy_bins);

//    step_levels = HistFac.makeHist<const std::string&>(
//                "Check pass count",
//                "Check",
//                "# passed",
//                BinSettings(10));

//    pion_mc_rec_angle = makeAngleDiffPlot(
//                ParticleTypeDatabase::Pi0.PrintName()+" MC/Rec angle",
//                "angle [#circ]",
//                "# / " + to_string(angle_diff_bins.BinWidth())+" #circ",
//                angle_diff_bins,
//                ParticleTypeDatabase::Eta.Name()+"_mc_rec_angle"
//                );
    n=0;
}

template <class InputIterator, class T>
T sum (InputIterator first, InputIterator last, T init) {
    while (first!=last) {
        init += **first;
        ++first;
    }
    return std::move(init);
}

template <class C, class T>
T sum (const C& data, T init) {
    return std::move(sum(data.begin(), data.end(), init));
}

void ant::analysis::Pion::ProcessEvent(const ant::Event &event)
{

    step_levels.Fill("0 Events Seen");

    if(event.Reconstructed().TriggerInfos().CBEenergySum()<550.0)
        return;

    step_levels.Fill("1 ESum Cut passed");

    const ParticleList& photons = event.Reconstructed().Particles().Get(ParticleTypeDatabase::Photon);

    if(photons.size()<3)
        return;

    step_levels.Fill("2 NPhotons 3+");

    //TODO const ptr?
    ParticlePtr mc_pion;
    for( const ParticlePtr& mcp : event.MCTrue().Intermediates().GetAll() ) {
        if(mcp->Type() == ParticleTypeDatabase::Pi0) {
            if(mc_pion!=nullptr)
                throw string("Multiple pions found in MC True");
            mc_pion = mcp;

         }
    }

    unsigned int n_pion_found = 0;
/*
    for( auto comb = makeCombination(photons,3); !comb.Done(); ++comb) {

        ParticleList ggg;
        ggg.assign(comb.begin(),comb.end());

        TLorentzVector pion = *comb.at(0)+*comb.at(1)+*comb.at(2);

        if( pion_im_cut.Contains(pion.M())) {
            step_levels.Fill("3 #pion IM cut passed");

            for( auto gcomb = makeCombination(ggg,2); !gcomb.Done(); ++gcomb) {

                TLorentzVector g1(*gcomb.at(0));
                TLorentzVector g2(*gcomb.at(1));
                TLorentzVector eta = g1 + g2;

                eta_IM.Fill(eta);

                if(eta_im_cut.Contains(eta.M())) {

                    step_levels.Fill("4 #eta IM cut passed");
                    pion_IM.Fill(pion);
                    n_pion_found++;

                    for( auto& taggerhit : event.Reconstructed().TaggerHits() ) {
                        if( tagger_energy_cut.Contains(taggerhit->PhotonEnergy())) {
                            TLorentzVector p = taggerhit->PhotonBeam() + target - pion;
                            p_MM.Fill(p);
                        }
                    }

                    if(mc_pion) {
                        pion_mc_rec_angle.Fill( make_pair(*mc_pion,pion) );
                    }
                }
            }

        }

    }
 */
    pion_rec_multi.Fill(n_pion_found);

    if(n_pion_found == 0) {
        nr_ngamma.Fill(photons.size());

        for( auto comb = makeCombination(photons,3); !comb.Done(); ++comb) {
            TLorentzVector m = sum(comb, TLorentzVector() );
            nr_3gim.Fill(m);
        }

        for( auto comb = makeCombination(photons,2); !comb.Done(); ++comb) {
            TLorentzVector m = sum(comb, TLorentzVector() );
            nr_2gim.Fill(m);
        }
    }

}


void ant::analysis::Pion::Finish()
{

}


void ant::analysis::Pion::ShowResult()
{
//    canvas("Pion (Reconstructed)") << pion_IM << eta_IM << p_MM << step_levels << pion_rec_multi << pion_mc_rec_angle << endc;
//    canvas("Pion (Not Reconstructed)") << nr_ngamma << nr_2gim << nr_3gim << endc;

}
