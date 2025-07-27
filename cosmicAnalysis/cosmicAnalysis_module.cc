#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft data products
#include "nusimdata/SimulationBase/MCParticle.h"
#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

// ROOT includes
#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"

#include <vector>
#include <string>


namespace duneana{
    class cosmicAnalysis: public art::EDAnalyzer{
        public:
            explicit cosmicAnalysis(fhicl::ParameterSet const& p);
            cosmicAnalysis(cosmicAnalysis const&)                   = delete;
            cosmicAnalysis(cosmicAnalysis&&)                        = delete;
            cosmicAnalysis& operator = (cosmicAnalysis const&)    = delete;
            cosmicAnalysis& operator = (cosmicAnalysis&&)           = delete;

            void beginJob() override;
            void analyze(art::Event const& e) override;

        private:
            art::InputTag fMCParticleTag;
            TTree* fTree;



            //recording data
            int fRun;
            int fSubRun;
            int fEvent;

            //just to make sure the detector geometry is fine
            std::string fDetectorName;
            double fDetHalfWidth;
            double fDetHalfHeight;
            double fDetLength;

            int fPrimPdg;
            double fPrimE;
            double fPrimVx, fPrimVy, fPrimVz;
            double fPrimPx, fPrimPy, fPrimPz;

            double fTrackLenghtInTPC;
            double fTPCEntryX, fTPCEntryY, fTPCEntryZ;
            double fTPCExitX, fTPCExitY, fTPCExitZ;

            int fNSecondaries;
            std::vector<int>    fSecondaryPdg;
            std::vector<double> fSecondaryE;
            std::vector<double> fSecondaryVx, fSecondaryVy, fSecondaryVz;
    };
}
duneana::cosmicAnalysis::cosmicAnalysis(fhicl::ParameterSet const& p)
    :EDAnalyzer(p),
    fMCParticleTag(p.get<art::InputTag>("MCParticleTag"))
{}



void duneana::cosmicAnalysis::beginJob(){
    art::ServiceHandle<art::TFileService> tfs;
    fTree = tfs->make<TTree>("cosmicTree", "Analysis");
    fTree->Branch("run",        &fRun,      "run/I");
    fTree->Branch("subrun",     &fSubRun,   "subrun/I");
    fTree->Branch("event",      &fEvent,     "event/I");

    fTree->Branch("detectorName",   &fDetectorName);
    fTree->Branch("detHalfWidth",   &fDetHalfWidth,     "detHalfWidth/D");
    fTree->Branch("detHalfHeight",  &fDetHalfHeight,    "detHalfHeight/D");
    fTree->Branch("detLength",      &fDetLength,        "detLength/D");

    fTree->Branch("primPdg",    &fPrimPdg,          "primPdg/I");
    fTree->Branch("primE",      &fPrimE,            "primE/D");
    fTree->Branch("primVx",     &fPrimVx,           "primVx/D");
    fTree->Branch("primVy",     &fPrimVy,           "primVy/D");
    fTree->Branch("primVz",     &fPrimVz,           "primVz/D");

    fTree->Branch("primPx",     &fPrimPx,           "primPx/D");
    fTree->Branch("primPy",     &fPrimPy,           "primPy/D");
    fTree->Branch("primPz",     &fPrimPz,           "primPz/D");
    
    fTree->Branch("tracklength",&fTrackLenghtInTPC, "tracklength/D");
    fTree->Branch("tpcEntryX",  &fTPCEntryX,        "tpcEntryX/D");
    fTree->Branch("tpcEntryY",  &fTPCEntryY,        "tpcEntryY/D");
    fTree->Branch("tpcEntryZ",  &fTPCEntryZ,        "tpcEntryZ/D");
    fTree->Branch("tpcExitX",   &fTPCExitX,         "tpcExitX/D");
    fTree->Branch("tpcExitY",   &fTPCExitY,         "tpcExitY/D");
    fTree->Branch("tpcExitZ",   &fTPCExitZ,         "tpcExitZ/D");

    fTree->Branch("nSecondaries",   &fNSecondaries, "nSecondaries/I");
    fTree->Branch("secondaryPdg",   &fSecondaryPdg);
    fTree->Branch("secondaryE",     &fSecondaryE);
    fTree->Branch("secondaryVx",    &fSecondaryVx);
    fTree->Branch("secondaryVy",    &fSecondaryVy);
    fTree->Branch("secondaryVz",    &fSecondaryVz);
}

void duneana::cosmicAnalysis::analyze(art::Event const& e){
    fRun    = e.run();
    fSubRun = e.subRun();
    fEvent  = e.id().event();
    
    art::ServiceHandle<geo::Geometry> geom;
    fDetectorName = geom->DetectorName();
    
    auto mcParticleHandle = e.getValidHandle<std::vector<simb::MCParticle>>(fMCParticleTag);
    auto const& mcParticleList = *mcParticleHandle;

    std::cout<<"\n--------------------------------------Cosmic Track Analysis--------------------------------------\n"<<std::endl;
    std::cout<<"Detector Name: "<<fDetectorName<<" .\n"<<std::endl;

    for (auto const& particle: mcParticleList){
        if(particle.Mother() !=0) continue;
        fSecondaryPdg.clear();
        fSecondaryE.clear();
        fSecondaryVx.clear();
        fSecondaryVy.clear();
        fSecondaryVz.clear();

        fPrimPdg    = particle.PdgCode();
        fPrimE      = particle.E();
        fPrimVx     = particle.Vx();
        fPrimVy     = particle.Vy();
        fPrimVz     = particle.Vz();
        fPrimPx     = particle.Px();
        fPrimPy     = particle.Py();
        fPrimPz     = particle.Pz();
        
        fTrackLengthInTPC       = 0.0;
        TVector3 tpcEntryPoint(-999,-999,-999), tpcExitPoint(-999,-999,-999);

        bool enterTPC           = false;
        const geo::TPCGeo& tpc  = geom->TPC(geom->GetBeginTPC().ID());
        fDetHalfWidth           = tpc.HalfWidth();
        fDetHalfHeight          = tpc.HalfHeight();
        fDetLength              = tpc.Length();

        for (size_t i=0; i<particle.NumberTrajectoryPoints(); ++i){
            TVector3 currentPoint = particle.Position(i).Vect();
            if(tpc.ContainsPosition(currentPoint)){
                tpcEntryPoint   = currentPoint;
                enterTPC      = true;
            }
            tpcExitPoint        = currentPoint;
        }
        if(enterTPC){
            fTPCEntryX          = tpcEntryPoint.X();
            fTPCEntryY          = tpcEntryPoint.Y();
            fTPCEntryZ          = tpcEntryPoint.Z();
            fTPCExitX           = tpcExitPoint.X();
            fTPCExitY           = tpcExitPoint.Y();
            fTPCExitZ           = tpcExitPoint.Z();
            fTrackLengthInTPC     =(tpcExitPoint-tpcEntryPoint).Mag();

        }
        else{
            fTPCEntryX  = -1; fTPCEntryY    = -1; fTPCEntryZ    = -1;
            fTPCExitX   = -1; fTPCExitY     = -1; fTPCExitZ     = -1;
        }
        int primTrackId = particle.TrackId();
        for(auto const& secondary: mcParticleList){
            if(secondary.Mother()==primaryTrackId){
                fSecondaryPdg.push_back(secondary.PdgCode());
                fSecondaryE.push_back(secondary.E());
                fSecondaryVx.push_back(secondary.Vx());
                fSecondaryVy.push_back(secondary.Vy());
                fSecondaryVz.push_back(secondary.Vz());
            }
        }
        fNSecondaries   = fSecondaryPdg.size();
        fTree->Fill();
    
    }
}
DEFINE_ART_MODULE(duneana::cosmicAnalysis)