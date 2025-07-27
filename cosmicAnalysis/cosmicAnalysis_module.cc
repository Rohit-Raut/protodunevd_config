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
            cosmicAnalysis& operator = (cosmicAnalysis const&())    = delete;
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

            int
    }
}