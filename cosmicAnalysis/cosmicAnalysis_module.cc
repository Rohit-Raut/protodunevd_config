#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "detdataformats/trigger/TriggerPrimitive.hpp"
#include "detdataformats/trigger/TriggerActivityData.hpp"
//#include "triggeralgs/TriggerActivity.hpp"
//#include "trigger/TriggerPrimitive.hpp"
#include "lardataobj/RawData/RawDigit.h"
// LArSoft data products
#include "nusimdata/SimulationBase/MCParticle.h"
#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcore/Geometry/WireReadout.h"
// ROOT includes
#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include <array>
#include <vector>
#include <string>
#include <cmath>

namespace duneana{
    class cosmicAnalysis: public art::EDAnalyzer{
        public:
            explicit cosmicAnalysis(fhicl::ParameterSet const& p);
            cosmicAnalysis(cosmicAnalysis const&)                   = delete;
            cosmicAnalysis(cosmicAnalysis&&)                        = delete;
            cosmicAnalysis& operator = (cosmicAnalysis const&)	    = delete;
            cosmicAnalysis& operator = (cosmicAnalysis&&)           = delete;
            void beginJob() override;
            void analyze(art::Event const& e) override;

        private:
            art::InputTag fMCParticleTag;
	        art::InputTag fTPLabel;
	        art::InputTag fTALabel;
            TTree* fTree;
	        TTree* fTreeTP;
	        TTree* fTreeTA;
            geo::WireReadoutGeom const& fWireReadoutGeom=art::ServiceHandle<geo::WireReadout>()->Get();

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

            double fTrackLengthInTPC;
            double fTPCEntryX, fTPCEntryY, fTPCEntryZ;
            double fTPCExitX, fTPCExitY, fTPCExitZ;

            int fNSecondaries;
            std::vector<int>    fSecondaryPdg;
            std::vector<double> fSecondaryE;
            std::vector<double> fSecondaryVx, fSecondaryVy, fSecondaryVz;
		
	    //TP data storage
	    double fTPTimeStart;
	    double fTPTimePeak;
	    double fTPTimeOverThreshold;
	    int fTPChannel;
	    double fTPADCPeak;
	    double fTPADCSum;
	    double fTPDetId;

	    //TA information
	    int fTANTPs;
	    double fTATimeStart;
	    double fTATimeEnd;
	    double fTATimePeak;
	    double fTAADCPeak;
	    double fTAADCSum;
	    int fTAChannelStart;
	    int fTAChannelEnd;
	    int fTAChannelPeak;
    };
}
duneana::cosmicAnalysis::cosmicAnalysis(fhicl::ParameterSet const& p)
    :EDAnalyzer(p),
    fMCParticleTag(p.get<art::InputTag>("MCParticleTag")),
    fTPLabel(p.get<art::InputTag>("TPLabel")),
    fTALabel(p.get<art::InputTag>("TALabel"))
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
    
    fTree->Branch("tracklength",&fTrackLengthInTPC, "tracklength/D");
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

    fTreeTP = tfs->make<TTree>("TP", "analysis");

    fTreeTP->Branch("run", 		&fRun, 		"run/I");
    fTreeTP->Branch("subrun", 		&fSubRun, 	"subrun/I");
    fTreeTP->Branch("event", 		&fEvent, 	"event/I");
    fTreeTP->Branch("TPStart", 		&fTPTimeStart, 	"timestart/D");
    fTreeTP->Branch("TPPeak", 		&fTPTimePeak, 	"TPPeak/D");
    fTreeTP->Branch("TPSum", 		&fTPADCSum, 	"TPSum/D");
    fTreeTP->Branch("TPTimeOverThreshold", &fTPTimeOverThreshold, "TPTimeOverThreshold/D");
    fTreeTP->Branch("TPChannel", 	&fTPChannel, 	"TPChannel/I");
    fTreeTP->Branch("TPADCPeak", 	&fTPADCPeak, 	"TPADCPeak/D");
    fTreeTP->Branch("TPDetId", 		&fTPDetId,	"TPDetId/D");
    
    //TA information is also recorded in TPTree
    fTreeTA = tfs->make<TTree>("TA", "analysis");
    fTreeTA->Branch("event", 		&fEvent, 	"event/I");
    fTreeTA->Branch("TAnum", 		&fTANTPs, 	"TAnum/I");
    fTreeTA->Branch("TAStart", 		&fTATimeStart, 	"TAStart/D");
    fTreeTA->Branch("TAEnd", 		&fTATimeEnd, 	"TAEnd/D");
    fTreeTA->Branch("TAPeak", 		&fTATimePeak, 	"TAPeak/D");	
    fTreeTA->Branch("TAADCPeak", 	&fTAADCPeak, 	"TAADCPeak/D");
    fTreeTA->Branch("TASum", 		&fTAADCSum, 	"TASum/D");
    fTreeTA->Branch("TAChannelStart", 	&fTAChannelStart,"TAChannelStart/D");
    fTreeTA->Branch("TAChannelEnd",	&fTAChannelEnd,  "TAChannelEnd/D");
    fTreeTA->Branch("TAChannelPeak", 	&fTAChannelPeak ,"TAChannelPeak/D");



}


void duneana::cosmicAnalysis::analyze(art::Event const& e){
    fRun    = e.run();
    fSubRun = e.subRun();
    fEvent  = e.id().event();
    
    art::ServiceHandle<geo::Geometry> geom;
    // geo::BoxBoundedGeo active_volume = geom->ActiveBoundedBox();
    fDetectorName = geom->DetectorName();
    
    auto mcParticleHandle = e.getValidHandle<std::vector<simb::MCParticle>>(fMCParticleTag);
    auto const& mcParticleList = *mcParticleHandle;

    std::cout<<"\n--------------------------------------Cosmic Track Analysis--------------------------------------\n"<<std::endl;
    std::cout<<"Detector Name: "<<fDetectorName<<" .\n"<<std::endl;

    //Testing with getting detector geometry and bound
    double det_top      = -9999.0;
    double det_bottom   =  9999.0;
    double det_front    =  9999.0;
    double det_back     =  -9999.0;
    double det_width    =  0.0;
    for(geo::TPCGeo const& TPC: geom->Iterate<geo::TPCGeo>()){
         auto const center   = TPC.GetCenter();
         double tpc_top      = center.Y()+TPC.HalfHeight();
         double tpc_bottom   = center.Y()-TPC.HalfHeight();
         double tpc_front    = center.Z()-TPC.HalfLength();
         double tpc_back     = center.Z()+TPC.HalfLength();
         if(tpc_top      >   det_top)    det_top     = tpc_top;
         if(tpc_bottom   <   det_bottom) det_bottom  = tpc_bottom;
         if(tpc_front    <   det_front)  det_front   = tpc_front;
         if(tpc_back     >   det_back)   det_back    = tpc_back;
         det_width    = TPC.DriftDistance();
    }
    fDetHalfHeight 	= std::abs(det_top)+std::abs(det_bottom);
    fDetLength  	= std::abs(det_front)+std::abs(det_back);
    fDetHalfWidth 	= det_width;
    std::cout<<"The detector geometry, height: ("<<det_top<<" , "<<det_bottom<<"), Detector Length: ("<<det_front<<" , "<<det_back<<" ). Drift Distance: "<<det_width<<". \n"<<std::endl;
    int primaryCount = 0;
    
    //Debugginf some of the geometry information
    std::array<size_t, 3> tpCountPerPlane = {0, 0, 0};

    //
    //TP information extraction
    art::Handle<std::vector<dunedaq::trgdataformats::TriggerPrimitive>> tpHandle = e.getHandle<std::vector<dunedaq::trgdataformats::TriggerPrimitive>>(fTPLabel);
    if (tpHandle.isValid()){
    	for (const dunedaq::trgdataformats::TriggerPrimitive &tp: *tpHandle){
            fTPTimeStart 	= tp.time_start;
            fTPTimePeak  	= tp.time_peak;
            fTPTimeOverThreshold 	= tp.time_over_threshold;
            fTPChannel		= tp.channel;

            int plane = fWireReadoutGeom.View(tp.channel);            
            if(plane<0 || plane >2){
                mf::LogWarning("CosmicAnalysis")<<"There should not be more than given plane";
                continue;
            }
            ++tpCountPerPlane[plane];
            fTPADCPeak 		= tp.adc_peak;
            fTPDetId		= tp.detid;
            fTreeTP->Fill();
      }
    }
    else{
      mf::LogWarning("TP Analysis")<< "No Trigger Primitive Found:  "<<fTPLabel<<"..........\n";
    }
    mf::LogInfo("CosmicAnalysis")<< "TP counts per plane → "<< "U=" << tpCountPerPlane[0] << ", V=" << tpCountPerPlane[1] << ", "<< "Z=" << tpCountPerPlane[2];

    //flushing out the Trigger Activity information
    art::Handle<std::vector<dunedaq::trgdataformats::TriggerActivityData>> taHandle = e.getHandle<std::vector<dunedaq::trgdataformats::TriggerActivityData>>(fTALabel);
    if(taHandle.isValid()){
      for (const auto& ta: *taHandle){
        //fTANTPs 	= ta.num_tps;
        fTATimeStart	= ta.time_start;
        fTATimeEnd	= ta.time_end;
        fTATimePeak	= ta.time_peak;
        fTAADCPeak	= ta.adc_peak;
        fTAADCSum	= ta.adc_integral;
        fTAChannelStart = ta.channel_start;
        fTAChannelEnd	= ta.channel_end;
        fTAChannelPeak	= ta.channel_peak;
        fTreeTA->Fill();


      }
    }
    else{
      mf::LogWarning("TA analysis")<<"No Trigger Activity Detected: "<<fTALabel<<" .........\n";
    }


    for (auto const& particle: mcParticleList){

      if(particle.Mother() ==0)
      {		
        primaryCount++;
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
        for (size_t i=0; i<particle.NumberTrajectoryPoints(); ++i){
          TVector3 currentPoint = particle.Position(i).Vect();
          bool isInside = false;

          //checking if the point are inside any of the TPCs
          for(geo::TPCGeo const& tpc: geom->Iterate<geo::TPCGeo>()){
            if(tpc.ContainsPosition(currentPoint)){
              isInside      = true;
              break;
            }
          }
          if(isInside){
            if(!enterTPC){
              tpcEntryPoint 	= currentPoint;
              enterTPC 	=true;
            }
            tpcExitPoint = currentPoint;
          }
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
          if(secondary.Mother()==primTrackId){
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
    std::cout<<"\n Number of Primary Recorded: "<<primaryCount<<".\n"<<std::endl;

}
DEFINE_ART_MODULE(duneana::cosmicAnalysis)
