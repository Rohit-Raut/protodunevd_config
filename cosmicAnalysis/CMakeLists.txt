include_directories("${dunedetdataformats_DIR}/../../../include")
cet_build_plugin(       cosmicAnalysis art::module 
                        LIBRARIES
			lardataobj::Simulation
                        lardataobj::AnalysisBase
                        nusimdata::SimulationBase
                        lardata::ArtDataHelper
                        art::Framework_Core
                        art::Framework_Principal
                        art::Framework_Services_Registry
                        art_root_io::tfile_support
                        art_root_io::TFileService_service
                        art::Persistency_Common
                        canvas::canvas
                        art::Persistency_Provenance
                        art::Utilities
                        messagefacility::MF_MessageLogger
                        cetlib::cetlib
                        cetlib_except::cetlib_except
                        ROOT::Core ROOT::Tree ROOT::Hist
                        duneprototypes::ProtoDUNEDataUtils
                        ProtoDUNEUtilities
                        BASENAME_ONLY
)
install_fhicl()
install_source()
