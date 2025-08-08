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
#include "lardata/Utilities/LArFFT.h"
#include "lardataobj/RawData/raw.h"
// ROOT includes
#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include <array>
#include <vector>
#include <string>
#include <cmath>
#include <limits>
constexpr int kMaxRawChannels   = 13000;
constexpr int kMaxRawTicks      = 9600;

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
            void reset();
        private:
            art::InputTag fMCParticleTag;
	        art::InputTag fTPLabel;
	        art::InputTag fTALabel;
            art::InputTag fRawDigitLabel;
            TTree* fTree;
	        TTree* fTreeTP;
	        TTree* fTreeTA;
            TTree* fRawDigitTree;
            geo::WireReadoutGeom const& fWireReadoutGeom=art::ServiceHandle<geo::WireReadout>()->Get();

            //recording data
            int fRun;
            int fSubRun;
            int fEvent;
            
            //double rawMax = std::numeric_limits<double>::lowest();
            //just to make sure the detector geometry is fine
            std::string fDetectorName;
            double fDetHalfWidth;
            double fDetHalfHeight;
            double fDetLength;

            int fNPrimaries;
            std::vector<int> fPrimPdg;
            std::vector<double> fPrimE;
            std::vector<double> fPrimVx, fPrimVy, fPrimVz;
            std::vector<double> fPrimPx, fPrimPy, fPrimPz;
            std::vector<double> fTrackLengthInTPC;
            std::vector<double> fTPCEntryX, fTPCEntryY, fTPCEntryZ;
            std::vector<double> fTPCExitX, fTPCExitY, fTPCExitZ;
            int fNSecondaries;
            std::vector<int>    fSecondaryPdg;
            std::vector<double> fSecondaryE;
            std::vector<double> fSecondaryVx, fSecondaryVy, fSecondaryVz;
		    std::vector<int>    fTPCParticlePdg;
            std::vector<double> fTPCParticleE;
            //TP data storage
    	    uint64_t fTPTimeStart;
    	    uint64_t fTPTimePeak;
    	    uint64_t fTPTimeOverThreshold;
    	    uint32_t    fTPChannel;
    	    uint64_t fTPADCPeak;
    	    uint64_t fTPADCSum;
            uint64_t fTPDetId;
            double fADCIntegralDAQ; 
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

            
            int fRawPlane, fRawChannel;
            int fRawTimeStart, fRawTimeEnd, fRawTimePeak;
            double fRawAdcPeak, fRawAdcIntegral;
            double fRawAdc;
            //Storing the raw didigt data
           // int fRaw_nChan;
           // int fRaw_channel[kMaxRawChannels];
           // int fRaw_plane[kMaxRawChannels];
           // int fRaw_ADCs[kMaxRawChannels][kMaxRawTicks];
    };
}
duneana::cosmicAnalysis::cosmicAnalysis(fhicl::ParameterSet const& p)
    :EDAnalyzer(p),
    fMCParticleTag(p.get<art::InputTag>("MCParticleTag")),
    fTPLabel(p.get<art::InputTag>("TPLabel")),
    fTALabel(p.get<art::InputTag>("TALabel")),
    fRawDigitLabel(p.get<art::InputTag>("RawDigitLabel"))
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

    
    // Store ALL primaries as vectors
    fTree->Branch("nPrimaries", &fNPrimaries, "nPrimaries/I");
    fTree->Branch("primPdg",    &fPrimPdg);
    fTree->Branch("primE",      &fPrimE);
    fTree->Branch("primVx",     &fPrimVx);
    fTree->Branch("primVy",     &fPrimVy);
    fTree->Branch("primVz",     &fPrimVz);
    fTree->Branch("primPx",     &fPrimPx);
    fTree->Branch("primPy",     &fPrimPy);
    fTree->Branch("primPz",     &fPrimPz);
    fTree->Branch("tracklength",&fTrackLengthInTPC);
    fTree->Branch("tpcEntryX",  &fTPCEntryX);
    fTree->Branch("tpcEntryY",  &fTPCEntryY);
    fTree->Branch("tpcEntryZ",  &fTPCEntryZ);
    fTree->Branch("tpcExitX",   &fTPCExitX);
    fTree->Branch("tpcExitY",   &fTPCExitY);
    fTree->Branch("tpcExitZ",   &fTPCExitZ);
    fTree->Branch("tpcE",       &fTPCParticleE);


    fTree->Branch("nSecondaries",   &fNSecondaries, "nSecondaries/I");
    fTree->Branch("secondaryPdg",   &fSecondaryPdg);
    fTree->Branch("secondaryE",     &fSecondaryE);
    fTree->Branch("secondaryVx",    &fSecondaryVx);
    fTree->Branch("secondaryVy",    &fSecondaryVy);
    fTree->Branch("secondaryVz",    &fSecondaryVz);
    fTree->Branch("TPCPdg",         &fTPCParticlePdg);
    fTreeTP = tfs->make<TTree>("TP", "analysis");

    fTreeTP->Branch("run", 		&fRun, 		"run/I");
    fTreeTP->Branch("subrun", 		&fSubRun, 	"subrun/I");
    fTreeTP->Branch("event", 		&fEvent, 	"event/I");
    fTreeTP->Branch("TPStart", 		&fTPTimeStart);
    fTreeTP->Branch("TPPeak", 		&fTPTimePeak);
    fTreeTP->Branch("TPSum", 		&fTPADCSum);
    fTreeTP->Branch("TPTimeOverThreshold", &fTPTimeOverThreshold);
    fTreeTP->Branch("TPChannel", 	&fTPChannel);
    fTreeTP->Branch("TPADCPeak", 	&fTPADCPeak);
    fTreeTP->Branch("TPDetId", 		&fTPDetId);
    fTreeTP->Branch("tp_adcIntegral", &fADCIntegralDAQ);
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

    fRawDigitTree = tfs->make<TTree>("rawDigitTree", "Raw Digit Waveform");
    fRawDigitTree->Branch("event",  &fEvent,    "event/I");
    fRawDigitTree->Branch("run",    &fRun,      "run/I");
    fRawDigitTree->Branch("subrun", &fSubRun,   "subrun/I");
    fRawDigitTree->Branch("plane",  &fRawPlane, "plane/I");
    fRawDigitTree->Branch("chan",   &fRawChannel, "chan/I");
    fRawDigitTree->Branch("timestart", &fRawTimeStart,  "timestart/I");
    fRawDigitTree->Branch("timeend", &fRawTimeEnd,      "timeend/I");
    fRawDigitTree->Branch("timepeak", &fRawTimePeak,    "timepeak/I");
    fRawDigitTree->Branch("ADCIntegral", &fRawAdcIntegral, "ADCIntegral/D");
    fRawDigitTree->Branch("ADCPeak",    &fRawAdcPeak,       "ADCPeak/D");
    fRawDigitTree->Branch("adcs",       &fRawAdc,           "adcs");


    //fRawDigitTree->Branch("nChan",  &fRaw_nChan,"nChan/I");
    //fRawDigitTree->Branch("channel",&fRaw_channel,     "channel[nChan]/I");
    //fRawDigitTree->Branch("plane",  &fRaw_plane,    "plane[nChan]/I");
    //fRawDigitTree->Branch("adcs",   &fRaw_ADCs,     Form("adcs[nChan][%d]/S", kMaxRawTicks));

}


void duneana::cosmicAnalysis::analyze(art::Event const& e){
    fRun    = e.run();
    fSubRun = e.subRun();
    fEvent  = e.id().event();
    reset(); 
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
    std::map<unsigned int, unsigned int> apa_to_crp_map;
    for(geo::TPCGeo const& TPC: geom->Iterate<geo::TPCGeo>()){
         auto const center   = TPC.GetCenter();
         double tpc_top      = center.Y()+TPC.HalfHeight();
         double tpc_bottom   = center.Y()-TPC.HalfHeight();
         double tpc_front    = center.Z()-TPC.HalfLength();
         double tpc_back     = center.Z()+TPC.HalfLength();
         unsigned int apa_number = TPC.ID().TPC;
         unsigned int final_crp_number = 0;
         if (center.X()>0){//top of the detector
            if(center.Y()<0)    final_crp_number = 2;
            else                final_crp_number = 3;
         }
         else{ //bottom crp
            if(center.Y()<0)    final_crp_number    = 5;
            else                final_crp_number    = 4;
         }
         apa_to_crp_map[apa_number] = final_crp_number;
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

    size_t nChannels    = fWireReadoutGeom.Nchannels();
    std::map<unsigned int, std::vector<raw::ChannelID_t>> crp_channel_counts;
    mf::LogInfo("CosmicAnalysis")<<"Mapping all "<<nChannels<<" channels to crps.....\n";
    for (raw::ChannelID_t channel =0; channel<nChannels; ++channel){
        std::vector<geo::WireID> wireIDs    = fWireReadoutGeom.ChannelToWire(channel);
        if(wireIDs.empty()){
            continue;        
        }
        unsigned int apa_number             = wireIDs[0].TPC;
        if(apa_to_crp_map.count(apa_number)){
            unsigned int final_crp_number   = apa_to_crp_map.at(apa_number);
            crp_channel_counts[final_crp_number].push_back(channel);
        }
    }
    int primaryCount = 0;
    
    auto rawDigitHandle = *e.getValidHandle<std::vector<raw::RawDigit>>(fRawDigitLabel);
    std::unordered_map<int, const raw::RawDigit*> rdMap;
    rdMap.reserve(rawDigitHandle.size());
    for (auto const& rd:rawDigitHandle){
        rdMap[rd.Channel()] = &rd;
    }


    //Debugginf some of the geometry information
    std::array<size_t, 3> tpCountPerPlane = {0, 0, 0};
    //TP information extraction
    art::Handle<std::vector<dunedaq::trgdataformats::TriggerPrimitive>> tpHandle = e.getHandle<std::vector<dunedaq::trgdataformats::TriggerPrimitive>>(fTPLabel);
    std::array<int, 3> rawCountPerPlane = {0, 0, 0};
    for (auto const &rd: rawDigitHandle){
        int p = fWireReadoutGeom.View(rd.Channel());
        if(p>=0 && p<3) ++rawCountPerPlane[p];
    }
    mf::LogInfo("CosmicAnalysis")<<"Raw Digit data found U-> "<<rawCountPerPlane[0]<< ", V-> "<<rawCountPerPlane[1]<<", Z-> "<<rawCountPerPlane[2];

    if(rawCountPerPlane[0]==0 || rawCountPerPlane[1] ==0){
        mf::LogWarning("CosmicAnalysis")<<"U/V raw digit missing. Debug the tpcrawdecoder fcl.";
    }
   art::Handle<std::vector<dunedaq::trgdataformats::TriggerActivityData>> taHandle = e.getHandle<std::vector<dunedaq::trgdataformats::TriggerActivityData>>(fTALabel); 
    //first try going through all the raw data and then look for the TP and get the raw data at the moment of TP
   //if(taHandle.isValid() && !taHandle->empty()){
   //    auto const& tpVec     = *taHandle; 
   //    for(auto const& rd: rawDigitHandle){
   //         int chan = rd.Channel();
   //         int plane = fWireReadoutGeom.View(chan);
   //         std::vector<short>wf(rd.Samples());
   //         raw::Uncompress(rd.ADCs(), wf, rd.GetPedestal(), rd.Compression());
   //         int tpTick  = tpVec.front().time_peak;
   //         int halfW   = 100;
   //         int start   = std::max(0, tpTick-halfW);
   //         int end     = std::min(int(wf.size())-1, tpTick+halfW);

   //         double ped  = rd.GetPedestal();
   //         double sum  = 0;
   //         double m    = std::numeric_limits<double>::lowest();
   //         for (int t=start; t<=end; ++t){
   //             double v= wf[t]-ped;
   //             sum    +=v;
   //             m       = std::max(m, v);
   //         }
   //         fRawPlane       = plane;
   //         fRawChannel     = chan;
   //         fRawTimeStart   = start;
   //         fRawTimeEnd     = end;
   //         fRawTimePeak    = tpTick;
   //         fRawAdcIntegral = sum;
   //         fRawAdcPeak     = m;
   //         fRawDigitTree->Fill();

   //     }
   // }


    if (tpHandle.isValid()){
    	for (const dunedaq::trgdataformats::TriggerPrimitive &tp: *tpHandle){
            fTPChannel		= tp.channel;
            int plane = fWireReadoutGeom.View(tp.channel);            
            if(plane<0 || plane>2)continue;
            fTPTimeStart 	= tp.time_start;
            fTPTimePeak  	= tp.time_peak;
            fTPTimeOverThreshold 	= tp.time_over_threshold;
            ++tpCountPerPlane[plane];
            fTPADCPeak 		= tp.adc_peak;
            fADCIntegralDAQ = tp.adc_integral;
            fTPDetId		= tp.detid;
            fTreeTP->Fill();
            
            auto it = rdMap.find(tp.channel);
            if(it==rdMap.end()){
                mf::LogWarning("CosmicAnalysis")<<"No RawDigit found for channel: "<<tp.channel;
                continue;
            }
            const int tstart = static_cast<int>(tp.time_start);
            const int tpeak  = static_cast<int>(tp.time_peak);
            const int tot    = static_cast<int>(tp.time_over_threshold);

            const raw::RawDigit& rd = *it->second;
            std::vector<short> wf(rd.Samples());
            raw::Uncompress(rd.ADCs(), wf, rd.GetPedestal(), rd.Compression());
            const double pedestal = rd.GetPedestal();
            // Use TP's own window; make wend exclusive; clamp to [0, wf.size()]
            const int nTicks = static_cast<int>(wf.size());
            int wstart = std::max<int>(0, tstart);
            int wend       = std::min<int>(tstart+tot, nTicks);

            // Fallback if TOT is 0/invalid
            if (wstart >= wend){
              constexpr int W = 100;
              wstart = std::max(0, tpeak - W);
              wend   = std::min(nTicks, tpeak + W + 1);
              if (wstart >= wend) continue;
            }

            double rawSum = 0.0;
            double rawMax = std::numeric_limits<double>::lowest();
            for (int t = wstart; t < wend; ++t){
              const double v = static_cast<double>(wf[t]) - pedestal;
              rawSum += v;
              if (v > rawMax) rawMax = v;
            }
            fRawPlane       = plane;
            fRawChannel     = tp.channel;
            fRawTimeStart   = wstart;
            fRawTimeEnd     = wend-1;
            fRawTimePeak    = tpeak;
            fRawAdcIntegral = rawSum;
            fRawAdcPeak     = rawMax;
            fRawDigitTree->Fill();
            // auto it = rdMap.find(tp.channel);
           // if(it==rdMap.end()){
           //     mf::LogWarning("CosmicAnalysis")<<"No Raw Digit found"<< tp.channel;
           // }
           // auto const& rd = *it->second;
           // std::vector<short> wf(rd.Samples());
           // raw::Uncompress(rd.ADCs(), wf, rd.GetPedestal(), rd.Compression());
           // int wstart  = std::max(0, tp.time_start);
           // int wend    = std::min((int)wf.size()-1, tp.time_start+tp.time_over_threshold);
           // double rawSum = 0;
           // double rawMax = std::numeric_limit<double>::lowest();
           // double pedestal = rd.GetPedestal();
           // for(int t=wstart; t<wend; ++t){
           //     double v = wf[t] - pedestal;
           //     rawSum +=v;
           //     rawMax  = std::max(rawMax, v);

           // }
           // fRawTimeStart   = wstart;
           // fRawTimeEnd     = wend;
           // fRawAdcIntegral = rawSum;
           // fRawAdcPeak     = rawMax;
           // fRawChannel     = tp.channel;
           // fRawPlane       = plane;
           // 

        }
    }
    else{
      mf::LogWarning("TP Analysis")<< "No Trigger Primitive Found:  "<<fTPLabel<<"..........\n";
    }
    mf::LogInfo("CosmicAnalysis")<< "TP counts per plane → "<< "U=" << tpCountPerPlane[0] << ", V=" << tpCountPerPlane[1] << ", "<< "Z=" << tpCountPerPlane[2];

    //flushing out the Trigger Activity information
    std::array<size_t, 3> taCountPerPlane = {};
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
        int planeno = fWireReadoutGeom.View(ta.channel_start);
        if(planeno<0 || planeno>2){
            mf::LogWarning("CosmicAnalysis")<<"Thres is issue with plane reading";
            continue;
        }
        ++taCountPerPlane[planeno];
        
    

        fTreeTA->Fill();


      }
    }
    else{
      mf::LogWarning("TA analysis")<<"No Trigger Activity Detected: "<<fTALabel<<" .........\n";
    }
    mf::LogInfo("CosmicAnalysis")<<"TA Count per plane -> U: "<<taCountPerPlane[0]<<", V: "<<taCountPerPlane[1]<<" , Z: "<<taCountPerPlane[2]; 

    //auto rawDigitHandle = e.getValidHandle<std::vector<raw::RawDigit>>(fRawDigitLabel);
    //if(rawDigitHandle.isValid()){
    //    mf::LogInfo("CosmicAnalysis")<<"Found "<<rawDigitHandle->size()<<" RawDigit in this event.";
    //    fRaw_nChan      = rawDigitHandle->size();
    //    int chan_idx    = 0;
    //    for (raw::RawDigit const& digit: *rawDigitHandle){
    //        if(chan_idx>=kMaxRawChannels){
    //            mf::LogWarning("CosmicAnalysis")<<"More channels than max channels currently set at 12k";
    //            break;
    //        }
    //        uint32_t chan   = digit.Channel();
    //        fRaw_channel[chan_idx]  = chan;
    //        fRaw_plane[chan_idx]    = fWireReadoutGeom.View(chan);
    //        std::vector<short> uncompressed(digit.Samples());
    //        raw::Uncompress(digit.ADCs(), uncompressed,digit.GetPedestal(), digit.Compression());
    //        for(size_t tick=0; tick<uncompressed.size(); ++tick){
    //            if(tick>=kMaxRawTicks) break;
    //            fRaw_ADCs[chan_idx][tick]   = uncompressed[tick];
    //        }
    //        chan_idx++;
    //    }

    //}
    //fRawDigitTree->Fill();
    


    fPrimPdg.clear();
    fPrimE.clear();
    fTPCParticlePdg.clear();
    fPrimVx.clear(); fPrimVy.clear(); fPrimVz.clear();
    fPrimPx.clear(); fPrimPy.clear(); fPrimPz.clear();
    fTrackLengthInTPC.clear();
    fTPCEntryX.clear(); fTPCEntryY.clear(); fTPCEntryZ.clear();
    fTPCExitX.clear(); fTPCExitY.clear(); fTPCExitZ.clear();
    fSecondaryPdg.clear();
    fSecondaryE.clear();
    fSecondaryVx.clear(); fSecondaryVy.clear(); fSecondaryVz.clear();
    fTPCParticleE.clear(); 
    for (auto const& particle: mcParticleList){

      if(particle.Mother() ==0)
      {		
        primaryCount++;
        fPrimPdg.push_back(particle.PdgCode());
        fPrimE.push_back(particle.E());
        fPrimVx.push_back(particle.Vx());
        fPrimVy.push_back(particle.Vy());
        fPrimVz.push_back(particle.Vz());
        fPrimPx.push_back(particle.Px());
        fPrimPy.push_back(particle.Py());
        fPrimPz.push_back(particle.Pz());



        double tracklength       = 0.0;
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
            fTPCParticlePdg.push_back(particle.PdgCode());
            fTPCEntryX.push_back(tpcEntryPoint.X());
            fTPCEntryY.push_back(tpcEntryPoint.Y());
            fTPCEntryZ.push_back(tpcEntryPoint.Z());
            fTPCExitX.push_back(tpcExitPoint.X());
            fTPCExitY.push_back(tpcExitPoint.Y());
            fTPCExitZ.push_back(tpcExitPoint.Z());
            fTPCParticleE.push_back(particle.E());
            tracklength     =(tpcExitPoint-tpcEntryPoint).Mag();
            
            fTrackLengthInTPC.push_back(tracklength);
        }
        //else{
        //    fTPCEntryX.push_back(-999);
        //    fTPCEntryY.push_back(-999);
        //    fTPCEntryZ.push_back(-999);
        //    fTPCExitX.push_back(-999);
        //    fTPCExitY.push_back(-999);
        //    fTPCExitZ.push_back(-999);
        //    fTrackLengthInTPC.push_back(-999);
        //}
    
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
        }  
    }
    std::cout<<"\n Number of Primary Recorded: "<<primaryCount<<".\n"<<std::endl;
    fNPrimaries = primaryCount;
    fNSecondaries = fSecondaryPdg.size();
    fTree->Fill();



}






void duneana::cosmicAnalysis::reset(){
    fRawChannel = 0;
    fRawChannel = -999;
    fRawPlane   = -999;
    fRawAdcIntegral = 0.0;
    fRawAdcPeak  = 0.0;
   // for(int i =0; i<kMaxRawChannels;++i){
   //    fRawChannel[i] = -999;
   //    fRawPlane  [i] = -999;
   //    fRawAdcIntegral[i] = 0.0;
   //    fRawAdcPeak[i]  = 0.0;
   //    //for(int j=0; j<kMaxRawTicks; ++j){
   //     //    fRawAdcPeak[i][j] = 0;
   //     //    fRawAdcIntegral[i][j]=0
   //     //}
   // }
}



DEFINE_ART_MODULE(duneana::cosmicAnalysis)
