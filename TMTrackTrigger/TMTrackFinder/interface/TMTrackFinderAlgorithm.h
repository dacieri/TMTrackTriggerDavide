

// system include files
#include <memory>
#include <bitset>


// user include files
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

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

#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/Ptr.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"

#include "SimTracker/TrackTriggerAssociation/interface/TTClusterAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTStubAssociationMap.h"

#include "DataFormats/L1TrackTrigger/interface/TTCluster.h"
#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TCanvas.h"

#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>


class TMTrackFinderAlgorithm  {

 public:

  TMTrackFinderAlgorithm(const edm::Handle< edmNew::DetSetVector< TTStub<Ref_PixelDigi_> > >&,
			 const edm::Handle< TTStubAssociationMap<Ref_PixelDigi_> >&,
			 const StackedTrackerGeometry*,
			 const edm::ParameterSet&,
			 const edm::Service<TFileService>);

  ~TMTrackFinderAlgorithm();


  
  typedef edm::Ref< edmNew::DetSetVector< TTStub<Ref_PixelDigi_> >, TTStub<Ref_PixelDigi_> > RefStub;

  struct extended_stub {
    int barcode; // unique counter associated with stub?
    TTStub<Ref_PixelDigi_> stub;
    RefStub ref;
    GlobalPoint position;
    double eta;
    double phi;
    double r;
    double z;
    int isReal;
    double dphi;
    int id; // tracker layer number.
    double pt;
    double measured_pt;
  };

  struct StubSort {
    template<typename T>
    bool operator() ( const T t1, const T t2) const {
      return (t1.barcode < t2.barcode);
    }
  } ;
  
  struct StubUnique {
    template<typename T>
    bool operator() ( const T t1, const T t2) const {
      return (t1.barcode == t2.barcode);
    }
  } ;

  typedef std::vector<extended_stub> Stubs;

  typedef std::map<std::pair<int,int>, Stubs> HoughArray;

  typedef std::vector<Stubs> Cells;

  typedef std::pair<int,double> Segment;

  typedef std::map<Segment, Stubs> SegmentStubs;

  typedef std::map<Segment, Cells> SegmentCells;

  typedef std::vector<Cells> Listcells;

  struct Newsegment {
    int segmentsize;
    int etabin;
    int betabin;
    Cells segmentcells;
  };

  typedef std::vector<Newsegment> Vsegment;

  
  const SegmentStubs& getSortedStubs() const;
  // Return stubs stored in each HT cells that has enough stubs to pass a cut. 
  const SegmentCells& getAcceptedCells() const;
  // Ditto, but only returning cells for which all the stubs also have consistent rapidity.
  const SegmentCells& getFilteredCells() const;

  const Listcells& getListCells() const;

  const Vsegment& getVsegment() const;

 private:

  // Main routine. Get stubs, fill HT array and identify interesting HT cells.
  void run();
  // Get info about stubs.
  void collect_stubs();
  // Fill r-phi HT array.
  void rphi_transform(Stubs&, HoughArray&, double);
  // Find HT cells with a good number of stubs. i.e. track candidates.
  void cell_readout(HoughArray&, Cells&);
  // Check that most of the stubs in these cells have consistent eta.
  void eta_filter(Cells&, Cells&);
  void bend_filter(Cells&, Cells&);
  void r_filter(Cells&, Cells&);
 
  
 
  const edm::Handle< edmNew::DetSetVector< TTStub<Ref_PixelDigi_> > > TTStubHandle_;
  const edm::Handle< TTStubAssociationMap<Ref_PixelDigi_> > MCTruthTTStubHandle_;
  const StackedTrackerGeometry* theStackedGeometry_;
  //  TH2D *h_hough;
       //*h_hough_noreg;
  const bool hough_print_;
  const bool debug_;
  const double pt_cut_;
  const double eta_cut_;
  const std::vector<double> eta_regions_;
  const double zError_;
  const double beta_range_;
  const double hough_pt_;
  const int hough_bins_x_;
  const int hough_bins_y_;
  const int sample_type_;
  const double pt_low_;
  const double pt_high_;
  const bool eta_register_;
  const bool broken_layer_;
  const int broken_layer_id_;
  const int nradii_;
  
  SegmentStubs sortedStubs_;
  SegmentCells accepted_cells_perSeg_;
  SegmentCells filtered_cells_perSeg_;
  Listcells AllCells;
  Vsegment v_segment;


};
