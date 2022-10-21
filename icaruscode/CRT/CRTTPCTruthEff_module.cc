/////////////////////////////////////////////////////////////////////////////
/// Class:       CRTTPCTruthEff
/// Module Type: analyzer
/// File:        CRTTPCTruthEff_module.cc
///
/// Author:         Tyler Boone
/// E-mail address: tboone@FNAL.gov
///
/// Modified from CRTT0Matching_module by Tyler Boone
/////////////////////////////////////////////////////////////////////////////

// sbndcode includes
#include "sbnobj/Common/CRT/CRTHit.hh"
#include "icaruscode/IcarusObj/CRTTPCMatchingInfo.h"
#include "icaruscode/CRT/CRTUtils/CRTT0MatchAlg.h"
#include "icaruscode/CRT/CRTUtils/CRTBackTracker.h"
#include "icaruscode/CRT/CRTUtils/RecoUtils.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Event.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Persistency/Common/PtrMaker.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"

#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>
#include <iostream>
#include <map>
#include <iterator>
#include <algorithm>

// LArSoft
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RawData/ExternalTrigger.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "sbnobj/Common/Trigger/ExtraTriggerInfo.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
// ROOT
#include "TVector3.h"
#include "TH1.h"
#include "TH2.h"
#include "TVector3.h"
#include "TTree.h"

namespace icarus {
  
  class CRTTPCTruthEff : public art::EDAnalyzer {
  public:

    explicit CRTTPCTruthEff(fhicl::ParameterSet const & p);

    // The destructor generated by the compiler is fine for classes
    // without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    CRTTPCTruthEff(CRTTPCTruthEff const &) = delete;
    CRTTPCTruthEff(CRTTPCTruthEff &&) = delete;
    CRTTPCTruthEff & operator = (CRTTPCTruthEff const &) = delete; 
    CRTTPCTruthEff & operator = (CRTTPCTruthEff &&) = delete;

    // Required functions.
    void analyze(const art::Event & e) override;

    // Selected optional functions.
    void beginJob() override;

    void endJob() override;

    void reconfigure(fhicl::ParameterSet const & p);
 
    void GetAncestorID(int trueid, int &motherid, int &ancestorid, int &layers, std::map<int, simb::MCParticle> all_particles);

    void getCatCrossXYZ(recob::Track trk, double &my_x, double &my_y, double &my_z);

  private:
    art::ServiceHandle<art::TFileService> tfs;
    icarus::crt::CRTBackTracker bt;

    // Params got from fcl file.......
    //    art::InputTag fTpcTrackModuleLabel; 	     ///< name of track producer
    std::vector<art::InputTag> fTpcTrackModuleLabel; ///< name of track producer
    std::vector<art::InputTag> fPFParticleLabel;     ///< labels for source of PFParticle
    art::InputTag              fCrtHitModuleLabel;   ///< name of crt producer
    art::InputTag              fTriggerLabel;        ///< labels for trigger
    art::InputTag 	       fSimModuleLabel;      ///< name of detsim producer
    CRTT0MatchAlg              t0Alg;		     ///< used to call matching functions
    bool                       fVerbose;       	     ///< print information about what's going on
    bool		       fIsData;              ///< switch for if this is data or MC

    geo::GeometryCore const* fGeometryService;   ///< pointer to Geometry provider
    icarus::crt::CRTCommonUtils* fCrtutils;
    //  CRTCommonUtils* fCrtutils;

    TTree* tr_crttpc;
    int fEvent;        			///< number of the event being processed
    int fRun;          			///< number of the run being processed
    int fSubRun;       			///< number of the sub-run being processed
    int crt_region;    			///< CRT hit region code
    int ttl_tpctrks, ttl_crthits; 	///< Total number of TPC tracks and CRT Hits for a given event
    int driftdir;			///< Drift direction of the track
    int num_crt_candidates;		///< Number of CRT Hit candidates for a given track, including the matched hti
    double crt_tpc_dca;			///< Distance of Closest Approach between CRT Hit and projected TPC track
    double t0min, t0max; 		///< Minimum and maximum possible T0 given the TPC trackorientation
    double crttime;			///< Matched CRT Hit timestamp
    double track_t0;			///< Reconstructed T0 of the TPC track if it has one, otherwise set to -99999
    double crt_pes; 			///< Total PEs of the matched CRT Hit
    double crt_extraplen;		///< Extrapolated length of the TPC track to the CRT

    bool has_crtmatch; 			///< Tells if the TPC track has a valid CRT Hit match or not by looking for a valid CRT Hit timestamp in the returned algorithm struct
    bool is_catcross; 			///< True if the track has a reconstructed T0, false otherwise
    bool simple_catcross;		///< True if my own "simple" method to determine if a track is CC says it is CC, false otherwise

    double tpc_trk_start_x, tpc_trk_start_y, tpc_trk_start_z;	///< XYZ coordinates for the TPC track start point
    double tpc_trk_end_x, tpc_trk_end_y, tpc_trk_end_z;	///< XYZ coordinates for the TPC track end point
    double trk_startdir_x, trk_startdir_y, trk_startdir_z;	///< XYZ components of directional vector pointing from TPC track start to the midpoint
    double trk_enddir_x, trk_enddir_y, trk_enddir_z;		///< XYZ components of directional vector pointing from TPC track end to the midpoint
    double catcross_x, catcross_y, catcross_z;			///< From my "simple"method of determining if a track is a cathode-crosser, my calculation of the crossing point for a given track
    double crt_x, crt_y, crt_z;					///< Matched CRT Hit XYZ

    //Simulation-only variables
    int crt_pdg, trk_pdg;  			///< PDG codes of the CRT Hit andd TPC track
    int track_trueID, crt_trueID;		///< Truth IDs of the TPC track and matched CRT Hit
    int track_mother, crt_mother; 		///< First mother of the CRT Hit/TPC track simulated particle truth ID
    int track_ancestor, crt_ancestor;		///< Truth ID of the CRT Hit/TPC Track ancestor, which I define as the ID that returns itself when querying the Mother ID
    int crt_motherlayers, track_motherlayers;   ///< Number of times the track/hit mother was queried before returning itself or 0
    //End simulation-only variables


    //These vectors contain the information for all CRT Hit candidates that were considered valid to be matched with the TPC track, including the one that was matched
    std::vector<int> all_crt_candidate_trueIDs, all_crt_candidate_motherIDs, all_crt_candidate_ancestorIDs;
    std::vector<int> all_crt_regions, all_crt_pdg, all_crt_candidate_motherlayers;
    std::vector<double> all_crt_candidate_x, all_crt_candidate_y, all_crt_candidate_z;
    std::vector<double> crt_start_dca, crt_end_dca, crt_timestamp;
    int best_dca_pos;			///<Contains position of "best match" for the vectors above

    //These may need to be commented out until I can edit sbnobj to expand what's contained in CRTHits
    std::vector<int> triggered_FEBs_mac5s; std::vector<double> triggered_FEBs_pes, triggered_FEBs_timestamps; ///<for by-FEB analysis of the matched CRT Hit

    //add trigger data product vars
    unsigned int m_gate_type;
    std::string  m_gate_name;
    uint64_t     m_trigger_timestamp;
    uint64_t     m_gate_start_timestamp;
    uint64_t     m_trigger_gate_diff;
    uint64_t     m_gate_crt_diff;

  }; // class CRTTPCTruthEff


  CRTTPCTruthEff::CRTTPCTruthEff(fhicl::ParameterSet const & p)
    : EDAnalyzer(p), t0Alg(p.get<fhicl::ParameterSet>("T0Alg"))
    , fCrtutils(new icarus::crt::CRTCommonUtils())
      // Initialize member data here, if know don't want to reconfigure on the fly
  {

    // Call appropriate produces<>() functions here.
//    produces< std::vector<anab::T0>                                              >();
//    produces< art::Assns<recob::Track , anab::T0, icarus::CRTTPCMatchingInfo>    >();
//    produces< art::Assns<sbn::crt::CRTHit, anab::T0, icarus::CRTTPCMatchingInfo> >();

    fGeometryService = lar::providerFrom<geo::Geometry>();
    reconfigure(p);

  } // CRTTPCTruthEff()


  void CRTTPCTruthEff::reconfigure(fhicl::ParameterSet const & p)
  {
    fTpcTrackModuleLabel = p.get< std::vector<art::InputTag>>("TpcTrackModuleLabel", {"pandoraTrackGausCryoE"});
    fCrtHitModuleLabel   = p.get<art::InputTag> ("CrtHitModuleLabel", "crthit"); 
    fPFParticleLabel    =  p.get< std::vector<art::InputTag> >("PFParticleLabel",             {""});
    fTriggerLabel        = p.get<art::InputTag>("TriggerLabel","daqTrigger");
  } // CRTTPCTruthEff::reconfigure()


  void CRTTPCTruthEff::beginJob()
  {
    // Access tfileservice to handle creating and writing histograms
    art::ServiceHandle<art::TFileService> tfs;

    tr_crttpc = tfs->make<TTree>("matchTree","CRTHit - TPC track matching analysis");

    tr_crttpc->Branch("fEvent",&fEvent,"fEvent/I");
    tr_crttpc->Branch("crt_region",&crt_region,"crt_region/I");
    tr_crttpc->Branch("fRun",&fRun,"fRun/I");
    tr_crttpc->Branch("driftdir",&driftdir,"driftdir/I");
    tr_crttpc->Branch("fSubRun",&fSubRun,"fSubRun/I");
    tr_crttpc->Branch("track_trueID",&track_trueID,"track_trueID/I");
    tr_crttpc->Branch("track_mother",&track_mother,"track_mother/I");
    tr_crttpc->Branch("crt_mother",&crt_mother,"crt_mother/I");
    tr_crttpc->Branch("track_ancestor",&track_ancestor,"track_ancestor/I");
    tr_crttpc->Branch("crt_ancestor",&crt_ancestor,"crt_ancestor/I");
    tr_crttpc->Branch("track_motherlayers",&track_motherlayers,"track_motherlayers/I");
    tr_crttpc->Branch("crt_motherlayers",&crt_motherlayers,"crt_motherlayers/I");
    tr_crttpc->Branch("crt_trueID",&crt_trueID,"crt_trueID/I");
    tr_crttpc->Branch("crt_pdg",&crt_pdg,"crt_pdg/I");
    tr_crttpc->Branch("trk_pdg",&trk_pdg,"trk_pdg/I");
    tr_crttpc->Branch("ttl_tpctrks",&ttl_tpctrks,"ttl_tpctrks/I");
    tr_crttpc->Branch("ttl_crthits",&ttl_crthits,"ttl_crthits/I");
    tr_crttpc->Branch("num_crt_candidates",&num_crt_candidates,"num_crt_candidates/I");
    tr_crttpc->Branch("best_dca_pos",&best_dca_pos,"best_dca_pos/I");
    tr_crttpc->Branch("crt_tpc_dca",&crt_tpc_dca,"crt_tpc_dca/D");
    tr_crttpc->Branch("t0min",&t0min,"t0min/D");
    tr_crttpc->Branch("t0max",&t0max,"t0max/D");
    tr_crttpc->Branch("crt_pes",&crt_pes,"crt_pes/D");
    tr_crttpc->Branch("crttime",&crttime,"crttime/D");
    tr_crttpc->Branch("crt_extraplen",&crt_extraplen,"crt_extraplen/D");
    tr_crttpc->Branch("is_catcross",&is_catcross,"is_catcross/O");
    tr_crttpc->Branch("simple_catcross",&simple_catcross,"simple_catcross/O");
    tr_crttpc->Branch("has_crtmatch",&has_crtmatch,"has_crtmatch/O");
    tr_crttpc->Branch("tpc_trk_start_x",&tpc_trk_start_x,"tpc_trk_start_x/D");
    tr_crttpc->Branch("tpc_trk_start_y",&tpc_trk_start_y,"tpc_trk_start_y/D");
    tr_crttpc->Branch("tpc_trk_start_z",&tpc_trk_start_z,"tpc_trk_start_z/D");
    tr_crttpc->Branch("tpc_trk_end_x",&tpc_trk_end_x,"tpc_trk_end_x/D");
    tr_crttpc->Branch("tpc_trk_end_y",&tpc_trk_end_y,"tpc_trk_end_y/D");
    tr_crttpc->Branch("tpc_trk_end_z",&tpc_trk_end_z,"tpc_trk_end_z/D");
    tr_crttpc->Branch("crt_x",&crt_x,"crt_x/D");
    tr_crttpc->Branch("crt_y",&crt_y,"crt_y/D");
    tr_crttpc->Branch("crt_z",&crt_z,"crt_z/D");
    tr_crttpc->Branch("catcross_x",&catcross_x,"catcross_x/D");
    tr_crttpc->Branch("catcross_y",&catcross_y,"catcross_y/D");
    tr_crttpc->Branch("catcross_z",&catcross_z,"catcross_z/D");
    tr_crttpc->Branch("trk_startdir_x",&trk_startdir_x,"trk_startdir_x/D");
    tr_crttpc->Branch("trk_startdir_y",&trk_startdir_y,"trk_startdir_y/D");
    tr_crttpc->Branch("trk_startdir_z",&trk_startdir_z,"trk_startdir_z/D");
    tr_crttpc->Branch("trk_enddir_x",&trk_enddir_x,"trk_enddir_x/D");
    tr_crttpc->Branch("trk_enddir_y",&trk_enddir_y,"trk_enddir_y/D");
    tr_crttpc->Branch("trk_enddir_z",&trk_enddir_z,"trk_enddir_z/D");
    tr_crttpc->Branch("all_crt_candidate_trueIDs",&all_crt_candidate_trueIDs);
    tr_crttpc->Branch("all_crt_candidate_motherIDs",&all_crt_candidate_motherIDs);
    tr_crttpc->Branch("all_crt_candidate_ancestorIDs",&all_crt_candidate_ancestorIDs);
    tr_crttpc->Branch("all_crt_candidate_motherlayers",&all_crt_candidate_motherlayers);
    tr_crttpc->Branch("all_crt_regions",&all_crt_regions);
    tr_crttpc->Branch("all_crt_pdg",&all_crt_pdg);
    tr_crttpc->Branch("all_crt_candidate_x",&all_crt_candidate_x);
    tr_crttpc->Branch("all_crt_candidate_y",&all_crt_candidate_y);
    tr_crttpc->Branch("all_crt_candidate_z",&all_crt_candidate_z);
    tr_crttpc->Branch("crt_start_dca",&crt_start_dca);
    tr_crttpc->Branch("crt_end_dca",&crt_end_dca);
    tr_crttpc->Branch("crt_timestamp",&crt_timestamp);
    tr_crttpc->Branch("track_t0",&track_t0,"track_t0/D");
    tr_crttpc->Branch("triggered_FEBs_mac5s",&triggered_FEBs_mac5s);
    tr_crttpc->Branch("triggered_FEBs_pes",&triggered_FEBs_pes);
    tr_crttpc->Branch("triggered_FEBs_timestamps",&triggered_FEBs_timestamps);


  } // CRTTPCTruthEff::beginJob()

  void CRTTPCTruthEff::analyze(const art::Event & event)
  {

    fEvent  = event.id().event();
    fRun    = event.run();
    fSubRun = event.subRun();

    std::map<int, simb::MCParticle> particles;

    if(!fIsData){
    // Loop over the true particles
    // mcpart_ids.clear();
   	 bt.Initialize(event);
    	    auto particleHandle = event.getValidHandle<std::vector<simb::MCParticle>>(fSimModuleLabel);
    	    for (auto const& particle: (*particleHandle)){
      
              // Make map with ID
              int partID = particle.TrackId();
	      particles[partID] = particle;
	    }//end loop over particles in particleHandle
     }//end(!fIsData)

    //add trigger info
    if( !fTriggerLabel.empty() ) {

      art::Handle<sbn::ExtraTriggerInfo> trigger_handle;
      event.getByLabel( fTriggerLabel, trigger_handle );
      if( trigger_handle.isValid() ) {
	sbn::triggerSource bit = trigger_handle->sourceType;
	m_gate_type            = (unsigned int)bit;
	m_gate_name            = bitName(bit);
	m_trigger_timestamp    = trigger_handle->triggerTimestamp;
	m_gate_start_timestamp = trigger_handle->beamGateTimestamp;
	m_trigger_gate_diff    = trigger_handle->triggerTimestamp - trigger_handle->beamGateTimestamp;
      }//end if( trigger_handle.isValid() )
      else{
//	mf::LogError("CRTTPCTruthEff:") << "No raw::Trigger associated to label: " << fTriggerLabel.label() << "\n" ;
      }//end else
    }//end if( !fTriggerLabel.empty() )
    else {
//      mf::LogError("CRTTPCTruthEff:") << "Trigger Data product " << fTriggerLabel.label() << " not found!\n" ;
    }//end else

    // Retrieve CRT hit list
    art::Handle<std::vector<sbn::crt::CRTHit>> crtListHandle;
    std::vector<art::Ptr<sbn::crt::CRTHit>> crtList;
    if(event.getByLabel(fCrtHitModuleLabel, crtListHandle))
      art::fill_ptr_vector(crtList, crtListHandle);

    std::vector<sbn::crt::CRTHit> crtHits;
    for (auto const& crtHit : crtList){
      crtHits.push_back(*crtHit);
    }//end for (auto const& crtHit : crtList)
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event);

    // Retrieve track list BEGIN LOOP OVER TRACKS IN EVENT
    for(const auto& trackLabel : fTpcTrackModuleLabel){

      auto it = &trackLabel - fTpcTrackModuleLabel.data();

      art::Handle< std::vector<recob::Track> > trackListHandle;
      std::vector<art::Ptr<recob::Track> > trackList;
      if (event.getByLabel(trackLabel,trackListHandle))
      	art::fill_ptr_vector(trackList, trackListHandle);   

//      mf::LogInfo("CRTTPCTruthEff")
//	<<"Number of reconstructed tracks = "<<trackList.size()<<"\n"
//	<<"Number of CRT hits = "<<crtList.size();

      ttl_tpctrks=(int)trackList.size();
      ttl_crthits=(int)crtList.size();

      //Get PFParticles
      auto pfpListHandle = event.getValidHandle<std::vector<recob::PFParticle> >(fPFParticleLabel[it]);
      if (!pfpListHandle.isValid()) continue;

      //Get PFParticle-Track association
      art::FindManyP<recob::PFParticle> fmpfp(trackListHandle, event, trackLabel);

      //Get T0-PFParticle association
      art::FindManyP<anab::T0> fmt0pandora(pfpListHandle, event, fPFParticleLabel[it]);

      
      if (trackListHandle.isValid() && crtListHandle.isValid() ){
	
	auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(event);
	art::FindManyP<recob::Hit> findManyHits(trackListHandle, event, trackLabel);
	
	// will create art pointers to the new T0 objects:
//	art::PtrMaker<anab::T0> const makeT0Ptr{ event };

	// Loop over all the reconstructed tracks 
	for(size_t track_i = 0; track_i < trackList.size(); track_i++) {

	  //clear out vectors for new entries into the TTree
	  all_crt_candidate_trueIDs.clear();
	  all_crt_candidate_motherIDs.clear();
	  all_crt_candidate_ancestorIDs.clear();
	  all_crt_regions.clear();
	  all_crt_pdg.clear();
	  all_crt_candidate_motherlayers.clear();
	  all_crt_candidate_x.clear();
	  all_crt_candidate_y.clear();
	  all_crt_candidate_z.clear();
	  crt_start_dca.clear();
	  crt_end_dca.clear();
	  crt_timestamp.clear();

	  //Set certain output variables to defaults
	  crt_tpc_dca = DBL_MAX; crt_trueID = INT_MAX; 
	  track_trueID = INT_MAX; track_mother = INT_MAX; track_ancestor = INT_MAX; track_motherlayers = INT_MAX;
	  crt_trueID = INT_MAX; crt_mother = INT_MAX; crt_ancestor = INT_MAX; crt_motherlayers = INT_MAX;
	  crttime = DBL_MAX; track_t0=DBL_MAX;
	  best_dca_pos=-99999;

	  //Find PFParticle for track i
	  //art::Ptr::key() gives the index in the vector
	  auto pfps = fmpfp.at(trackList[track_i]->ID());

	  is_catcross = false;
	  if (!pfps.empty()){
	    //Find T0 for PFParticle
	    auto t0s = fmt0pandora.at(pfps[0].key());
	    if (!t0s.empty()){
	      track_t0 = t0s[0]->Time();
	      is_catcross = true;
	      catcross_x = DBL_MAX;	catcross_y = DBL_MAX;	catcross_z = DBL_MAX;
	      getCatCrossXYZ(*trackList[track_i], catcross_x, catcross_y, catcross_z);
//	      t0 = track_t0;   //Get T0
	    }//end if (!t0s.empty())
	  }//end if (!pfps.empty())

	  crt_tpc_dca = DBL_MAX; crt_trueID = INT_MAX; track_trueID = INT_MAX; track_mother = INT_MAX; crttime = DBL_MAX;

	  std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(trackList[track_i]->ID());

	  if (hits.size() == 0) continue;
	  //Get Truth info if applicable:
	  if(!fIsData){
		track_trueID = RecoUtils::TrueParticleIDFromTotalRecoHits(clockData,hits,false);
		GetAncestorID(track_trueID, track_mother, track_ancestor, track_motherlayers, particles);
		trk_pdg = particles[track_trueID].PdgCode();
	  }//end if(!fIsData)

//	  int const cryoNumber = hits[0]->WireID().Cryostat;
	  matchCand closest = t0Alg.GetClosestCRTHit(detProp, *trackList[track_i], hits, crtHits,  m_gate_start_timestamp, fIsData);

	  if(closest.dca >=0 ){
	    
/*	    icarus::CRTTPCMatchingInfo matchInfo {
	        closest.dca       // DCA
	      , closest.extrapLen // extrapLength
	    }//end definition of icarus::CRTTPCMatchingInfo matchInfo;
*/
//	    double sin_angle = -99999;
	    if(closest.dca != -99999){//this if gives green light to start filling TTree variables for matched CRT Hit variables
	      	auto start = trackList[track_i]->Vertex<TVector3>();
		crt_region = fCrtutils->AuxDetRegionNameToNum(closest.thishit.tagger);	    
		crt_tpc_dca = closest.dca;

		driftdir =  closest.driftdir; //TPCGeoUtil::DriftDirectionFromHits(fGeometryService, hits);
		simple_catcross = closest.simple_cathodecrosser;

		best_dca_pos = closest.best_DCA_pos;
		crttime = closest.t0;
		crt_pes= closest.thishit.peshit;
		crt_extraplen = closest.extrapLen;

		tpc_trk_start_x = closest.tpc_track_start.X();
		tpc_trk_start_y = closest.tpc_track_start.Y();
		tpc_trk_start_z = closest.tpc_track_start.Z();

		tpc_trk_end_x = closest.tpc_track_end.X();
		tpc_trk_end_y = closest.tpc_track_end.Y();
		tpc_trk_end_z = closest.tpc_track_end.Z();
	
		trk_startdir_x = closest.startDir.X();
		trk_startdir_y = closest.startDir.Y();
		trk_startdir_z = closest.startDir.Z();
	
		trk_enddir_x = closest.endDir.X();
		trk_enddir_y = closest.endDir.Y();
		trk_enddir_z = closest.endDir.Z();

		crt_x = closest.thishit.x_pos;
		crt_y = closest.thishit.y_pos;
		crt_z = closest.thishit.z_pos;

		t0min = closest.t0min;
		t0max = closest.t0max;

 		std::vector<icarus::match_geometry> all_crt_candidates = t0Alg.GetClosestCRTHit_geo(detProp, *trackList[track_i], hits, crtHits, m_trigger_timestamp, fIsData);

		if(!fIsData){
			auto checkhit = closest.thishit;
			crt_trueID =  bt.TrueIdFromTotalEnergy(event, checkhit);// std::cout << temp_crt_trueID << std::endl;
			GetAncestorID(crt_trueID, crt_mother, crt_ancestor, crt_motherlayers, particles);
			crt_pdg = particles[crt_trueID].PdgCode();
	    	}//end if(!fIsData)

		has_crtmatch=false;
		num_crt_candidates = (int)all_crt_candidates.size();
		if(num_crt_candidates>0){
			has_crtmatch=true;
			for(int i=0; i<num_crt_candidates; i++){

				auto thiscand = all_crt_candidates[i];
				sbn::crt::CRTHit crtcandidate = thiscand.thishit;
		
				if(!fIsData){ 
					int temp_crt_trueid=INT_MAX, temp_crt_motherid=INT_MAX, temp_crt_ancestorid=INT_MAX, temp_crt_motherlayers=INT_MAX;
					temp_crt_trueid = bt.TrueIdFromTotalEnergy(event, crtcandidate);
				 	GetAncestorID(temp_crt_trueid, temp_crt_motherid, temp_crt_ancestorid, temp_crt_motherlayers, particles);
					all_crt_pdg.push_back(particles[temp_crt_trueid].PdgCode()); 
	
					all_crt_candidate_trueIDs.push_back(temp_crt_trueid);
					all_crt_candidate_motherIDs.push_back(temp_crt_motherid);
					all_crt_candidate_ancestorIDs.push_back(temp_crt_ancestorid);
					all_crt_candidate_motherlayers.push_back(temp_crt_motherlayers);

				}//end if(!fIsData)

				all_crt_candidate_x.push_back(thiscand.crt_hit_pos.X());
				all_crt_candidate_y.push_back(thiscand.crt_hit_pos.Y());
				all_crt_candidate_z.push_back(thiscand.crt_hit_pos.Z());

				all_crt_regions.push_back((int)fCrtutils->AuxDetRegionNameToNum(thiscand.thishit.tagger));
				crt_timestamp.push_back(thiscand.crtTime);

				crt_start_dca.push_back(thiscand.simpleDCA_startDir);	
				crt_end_dca.push_back(thiscand.simpleDCA_endDir);	

			}//end loop over CRT candidates
		}//end if(num_crt_candidates>0)

	    }//end if(closest.dca != -99999)

	  } // DCA check
	  
   	    tr_crttpc->Fill();
	} // Loop over tracks  
	
      } // Validity check

    } // all track labels in a vector 
  } // CRTTPCTruthEff::produce()


  void CRTTPCTruthEff::endJob()
  {

  } // CRTTPCTruthEff::endJob()

void icarus::CRTTPCTruthEff::GetAncestorID(int trueid, int &motherid, int &ancestorid, int &layers, std::map<int, simb::MCParticle> all_particles){

	int size_part = (int)all_particles.size();
	if(trueid<size_part) ancestorid = all_particles[trueid].Mother();
	else ancestorid = trueid;
	motherid = ancestorid;
	int last_ancestorid;
	bool hitbottom = false;
	int counter = 0;
	while(ancestorid%10000000!=0 && counter<1000 && !hitbottom){
		last_ancestorid = ancestorid;
		if(ancestorid<size_part) ancestorid = all_particles[ancestorid].Mother();
		hitbottom = last_ancestorid == ancestorid;
		counter++;
	}//end while loop 
	layers = counter;


}//CRTTPCTruthEff::GetAncestorID

void icarus::CRTTPCTruthEff::getCatCrossXYZ(recob::Track trk, double &my_x, double &my_y, double &my_z){

	size_t ntrk = trk.NPoints();
	std::vector<double> x, y, z, x_left_diff, x_right_diff;
	std::pair<double,double> cathode_yz;
	double left_dist_min = DBL_MAX; int left_dist_min_pos = INT_MAX;
	double right_dist_min = DBL_MAX; int right_dist_min_pos = INT_MAX;
	double xsum=0, xavg;

	int numelmnts = 0;

	//Begin by extracting coordinates into arrays
	for(size_t i=0; i<ntrk; i++){
		geo::Point_t thispt = trk.LocationAtPoint((int)i);
		x.push_back(thispt.X()); y.push_back(thispt.Y()); z.push_back(thispt.Z()); numelmnts++;
		if(std::abs(x.back())<210.215) {
			x_left_diff.push_back(std::abs(std::abs(x[i]) - 210.25));
			if(x_left_diff.back()<=left_dist_min) { left_dist_min=x_left_diff.back(); left_dist_min_pos = i; }
		}//end if(std::abs(x[i])<210.215)
		else if(std::abs(x.back())>=210.215) {
			x_right_diff.push_back(std::abs(std::abs(x[i]) - 210.25));
			if(x_right_diff.back()<=right_dist_min) { right_dist_min=x_right_diff.back(); right_dist_min_pos = i; }
		}//end if(std::abs(x[i])<210.215)
		xsum+=x.back();
	}//end loop over track points
	xavg = xsum/ntrk;

	if(left_dist_min_pos<= (numelmnts-1) && right_dist_min_pos<= (numelmnts-1)) {

		TVector3 leftpt(x[left_dist_min_pos],y[left_dist_min_pos],z[left_dist_min_pos]);
		TVector3 rightpt(x[right_dist_min_pos],y[right_dist_min_pos],z[right_dist_min_pos]);

		//calculate parameters for a line that passes through the two points:
		double mxy = (leftpt.Y() - rightpt.Y())/(leftpt.X()-rightpt.X());
		double mxz = (leftpt.Z() - rightpt.Z())/(leftpt.X()-rightpt.X());

		double bxy = leftpt.Y() - mxy*leftpt.X();
		double bxz = leftpt.Z() - mxz*leftpt.X();

		if(xavg>0) my_x = 210.215; 
		else if(xavg<=0) my_x = -210.215;
		my_y = mxy*my_x + bxy;
		my_z = mxz*my_x + bxz;
	}

}//end definition of std::pair<int,int> getCatCrossYZ(std::vector<art::Ptr<recob::Hit>> trk_hits)

  DEFINE_ART_MODULE(CRTTPCTruthEff)

} // sbnd namespace

namespace {

}
