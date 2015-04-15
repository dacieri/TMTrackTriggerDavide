
// system include files
#include <memory>
#include <bitset>

#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefProd.h"
#include "DataFormats/Common/interface/RefVector.h"

// user include files
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "TF1.h"
#include "TFitResult.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"

#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "SimDataFormats/TrackerDigiSimLink/interface/PixelDigiSimLink.h"
#include "DataFormats/Common/interface/DetSet.h"
#include "DataFormats/Common/interface/DetSetVector.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/TrackingGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"

#include "Geometry/Records/interface/StackedTrackerGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/StackedTrackerGeometry.h"
#include <DataFormats/TrackReco/interface/Track.h>
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/Ptr.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTClusterAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTStubAssociationMap.h"

#include "DataFormats/L1TrackTrigger/interface/TTCluster.h"
#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TCanvas.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 

#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>

#include "TMTrackTrigger/TMTrackFinder/interface/TMTrackFinderAlgorithm.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"




class TMTrackFinder : public edm::EDAnalyzer {
public:

  explicit TMTrackFinder(const edm::ParameterSet&);
  ~TMTrackFinder();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  
private:


  typedef edm::Ptr<TrackingParticle> TPptr;

  typedef std::vector<TPptr> TPs;


  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;


  void find_tracking_particles(TMTrackFinderAlgorithm::SegmentCells&, TPs&);
  void efficiency(TPs&);
  // Plot kinematics of offline reconstructed muons.
  void reconstruction0(const edm::Handle <std::vector<reco::Muon> >& );
  // Plot kinematics of generator particles.
  void generation(const edm::Handle <std::vector<reco::GenParticle> >&);
  // Plot kinematics of offline reconstructed electrons.
  void reconstruction1(const edm::Handle <std::vector<reco::GsfElectron> >& );


  // ---------- member data ---------------------------

  double pt_cut_;
  double eta_cut_;
  bool hough_print_;
  bool debug_;
  std::vector<double> eta_regions_;
  double zError_;
  double beta_range_;
  double hough_pt_;
  int hough_bins_x_;
  int hough_bins_y_;
  int sample_type_;
  double pt_low_;
  double pt_high_;
  bool reco_;
  bool eta_register_;
  edm::Service<TFileService> fs;
  int stop, lessthanten;
  edm::ESHandle<StackedTrackerGeometry> StackedGeometryHandle;
  edm::Handle< std::vector<TrackingParticle> > TrackingParticleHandle_;
  edm::Handle< edmNew::DetSetVector< TTCluster<Ref_PixelDigi_> > > TTClusterInclusiveHandle_;
  edm::Handle< TTClusterAssociationMap<Ref_PixelDigi_> > MCTruthTTClusterHandle_;
  edm::Handle< edmNew::DetSetVector< TTStub<Ref_PixelDigi_> > > TTStubHandle_;
  edm::Handle< TTStubAssociationMap<Ref_PixelDigi_> > MCTruthTTStubHandle_;
  edm::Handle< edmNew::DetSetVector< TTCluster<Ref_PixelDigi_> > > TTClusterHandle_;
  edm::Handle< edmNew::DetSetVector< TTStub<Ref_PixelDigi_> > > TTStubRejectedHandle_;
  edm::Handle<std::vector<reco::Muon> > RecoMuons_;
  edm::Handle<std::vector<reco::GenParticle> > GenParticles_;
  edm::Handle<std::vector<reco::GsfElectron> > RecoElectrons_;
  edm::Handle<std::vector<PileupSummaryInfo> > PileUpInfo_;
  const StackedTrackerGeometry* theStackedGeometry_;
  
  edm::ParameterSet conf_;




  // ---------------- histograms ---------------------
 
  std::vector<TH1D*> h_acc_cells;
  std::vector<TH1D*> h_stubs_acc_cells;
  std::vector<TH1D*> h_filt_cells;
  std::vector<TH1D*> h_stubs_filt_cells;
  std::vector<TH1D*> h_fake_cells;

  TH1D *h_acc_cells0,*h_acc_cells1,*h_acc_cells2,*h_acc_cells3,*h_acc_cells4;
  TH1D *h_fake_cells0,*h_fake_cells1,*h_fake_cells2,*h_fake_cells3,*h_fake_cells4;

  TH1D *h_stubs_acc_cells0,*h_stubs_acc_cells1,*h_stubs_acc_cells2,*h_stubs_acc_cells3,*h_stubs_acc_cells4;
  TH1D *h_filt_cells0,*h_filt_cells1,*h_filt_cells2,*h_filt_cells3,*h_filt_cells4;
  TH1D *h_stubs_filt_cells0,*h_stubs_filt_cells1,*h_stubs_filt_cells2,*h_stubs_filt_cells3,*h_stubs_filt_cells4;
  TH1D *h_candidate;
  TH1D *h_hough_pt, *h_cells_pt, *h_tracks_pt,*h_reco_pt,*h_gen_pt,*h_eff0_pt,*h_eff1_pt;
  TH1D *h_stubs_good, *h_stubs_bad;
  TH1D *h_cells_fake, *h_cells_duplicate;
  TH1D *h_cells_total;
  TH2D *h_cells_etaphi, *h_tracks_etaphi,*h_tracks_reg_etaphi ;
  TH2D *h_eff0,*h_eff1,*h_eff2, *h_reco_etaphi;
  TH1D *h_cells_eta, *h_tracks_eta,*h_ghost_eta,*h_ghost_phi;
  TH1D *h_cells_phi, *h_tracks_phi, *h_cells_id;
  TH1D *h_tracks_reg_eta, *h_tracks_reg_phi, *h_tracks_reg_pt;
  TH1D *h_reco_eta, *h_reco_phi;
  TH1D *h_eff0_eta, *h_eff0_phi;
  TH1D *h_eff1_eta, *h_eff1_phi;
  TH1D *h_gen_eta, *h_gen_phi;
  TH2D *h_gen_etaphi;
  TH2D *h_ghost_etaphi;
  TH1D *h_ghost_stubs, *h_ghost_register, *h_ghost;
  TH2D *h_ghost_regpt;
  TH1D *h_ghost_id, *h_tracks_id, *h_tracks_reg_id;
  TH1D *h_fake_segment;

  TH3D *h_cells_EtaPhiPt, *h_gen_EtaPhiPt, *h_ghost_EtaPhiPt;
  TH2D *stub_zr_10,  *stub_zr_11, *stub_zr_12, *stub_zr_13, *stub_zr_14;
  TH1D *h_barcode, *h_ghost_pt;
  TH1D *h_eff2_eta, *h_eff2_phi;
  TH1D *h_goodstubs;
  TH1D *h_overlap, *h_dupstubs;
  TH1D *h_pileup;
};
