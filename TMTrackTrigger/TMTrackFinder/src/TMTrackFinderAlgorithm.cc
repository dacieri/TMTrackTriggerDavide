

#include "TMTrackTrigger/TMTrackFinder/interface/TMTrackFinderAlgorithm.h"


TMTrackFinderAlgorithm::TMTrackFinderAlgorithm(const edm::Handle< edmNew::DetSetVector< TTStub<Ref_PixelDigi_> > >& TTStubHandle,
					       const edm::Handle< TTStubAssociationMap<Ref_PixelDigi_> >& MCTruthTTStubHandle,
					       const StackedTrackerGeometry* theStackedGeometry,
					       const edm::ParameterSet& conf,
					       const edm::Service<TFileService> fs):

  TTStubHandle_(TTStubHandle),
  MCTruthTTStubHandle_(MCTruthTTStubHandle),
  theStackedGeometry_(theStackedGeometry),
  //h_hough(fs->make<TH2D>("h_hough","",32,-0.003,0.003,32,0.0,0.2)),
  //h_hough_noreg(fs->make<TH2D>("h_hough","",32,-0.003,0.003,32,0.0,0.2)),
  hough_print_(conf.getParameter<bool>("HoughPrint")),
  debug_(conf.getParameter<bool>("Debug")),
  pt_cut_(conf.getParameter<double>("StubMinPt")),
  eta_cut_(conf.getParameter<double>("StubEtaRange")),
  eta_regions_(conf.getParameter<std::vector<double> >("EtaRegions")),
  zError_(conf.getParameter<double>("BeamErrorZ")),
  beta_range_(conf.getParameter<double>("BetaRange")),
  hough_pt_(conf.getParameter<double>("MinHoughPt")),
  hough_bins_x_(conf.getParameter<int>("HoughBinsX")),
  hough_bins_y_(conf.getParameter<int>("HoughBinsY")),
  sample_type_(conf.getParameter<int>("SampleType")),
  pt_low_(conf.getParameter<double>("PtLow")),
  pt_high_(conf.getParameter<double>("PtHigh")),
  eta_register_(conf.getParameter<bool>("EtaRegister")),
  broken_layer_(conf.getParameter<bool>("BrokenLayer")),
  broken_layer_id_(conf.getParameter<int>("BrokenLayerId")),
  nradii_(conf.getParameter<int>("NRadii"))
{

  TMTrackFinderAlgorithm::run();

}



TMTrackFinderAlgorithm::~TMTrackFinderAlgorithm()
{

}



void TMTrackFinderAlgorithm::run()
{

  // Main routine. Get stubs and fill HT array and identify interesting HT cells.

  collect_stubs();

  

  // Loop over eta/phi sectors, looking at the stubs in each.
  SegmentStubs::iterator segment = sortedStubs_.begin();
  
  for ( ; segment!=sortedStubs_.end(); segment++) {
    int eta_region = segment->first.first;
    double beta_seg = segment->first.second;
    //std::cout<< "beta_seg = " << beta_seg << std::endl;
    Stubs segment_stubs = segment->second; 


    // Create an empty Hough-Transform array and fill it with the stubs.
    HoughArray rphi_array;
    rphi_transform(segment_stubs,rphi_array,beta_seg);
    
    // Store, according to their eta/phi segment, all stubs found in each HT cell which passes a cut on the 
    // min. required no. of stubs inside them.
    Cells accepted_cells;
    cell_readout(rphi_array,accepted_cells);
    accepted_cells_perSeg_[segment->first] = accepted_cells;


    // Also require the stubs in each accepted HT cell to have consistent eta.
    Cells filtered_cells;
    eta_filter(accepted_cells,filtered_cells);
    filtered_cells_perSeg_[segment->first] = filtered_cells;

    AllCells.push_back(filtered_cells);
    Newsegment newseg;
    newseg.segmentsize=filtered_cells.size();
    newseg.etabin=eta_region;
    newseg.betabin=(beta_seg+3.2)*10;
    newseg.segmentcells=filtered_cells;
    v_segment.push_back(newseg);
  }
}



const TMTrackFinderAlgorithm::SegmentStubs &TMTrackFinderAlgorithm::getSortedStubs() const
{
  return sortedStubs_;
}



const TMTrackFinderAlgorithm::SegmentCells &TMTrackFinderAlgorithm::getAcceptedCells() const
{
  return accepted_cells_perSeg_;
}



const TMTrackFinderAlgorithm::SegmentCells &TMTrackFinderAlgorithm::getFilteredCells() const
{
  
  return filtered_cells_perSeg_;
}

const TMTrackFinderAlgorithm::Listcells &TMTrackFinderAlgorithm::getListCells() const
{
  return AllCells;
}

const TMTrackFinderAlgorithm::Vsegment &TMTrackFinderAlgorithm::getVsegment() const
{
  return v_segment;
}



void TMTrackFinderAlgorithm::collect_stubs()
{
  // Collect useful info about the stubs, storing it in the "extended_stub" object, and file these according to the eta/phi sector that the stub belongs to.

  int counter = 0; 
  
  edmNew::DetSetVector< TTStub<Ref_PixelDigi_> >::const_iterator module =  TTStubHandle_->begin();
  for ( ; module!=TTStubHandle_->end(); module++) {

    edmNew::DetSet< TTStub<Ref_PixelDigi_> >::const_iterator stub =  module->begin();
    for ( ; stub!=module->end(); stub++) {
      GlobalPoint pos = theStackedGeometry_->findGlobalPosition(stub);
      StackedTrackerDetId stDetId = stub->getDetId();
      RefStub theRef = edmNew::makeRefTo(TTStubHandle_, stub );
      int genuine = (MCTruthTTStubHandle_->isGenuine(theRef)==1) ? 1 : 0;
      // pt of track generating stub only applies to genuine correlations
      // fake correlations are assigned pt = 0
      double pt = 0;
      //int pdg = 0;
      if (genuine) {
      	edm::Ptr< TrackingParticle > tpPtr = MCTruthTTStubHandle_->findTrackingParticlePtr(theRef);
      	pt = tpPtr->pt();
        //std::vector<SimTrack> g4tracks = tpPtr->g4Tracks();
      	//pdg = g4tracks.begin()->type();
      }


      // bend direction (=distance in no. of strips between the two hits in the stub) is reversed for the +ve endcap

      double dir = stub->getTriggerBend() * ((pos.eta()>0)&&(stDetId.isEndcap()) ? -1: 1);


      // extract dphi from bend information and sensor separation in rphi (delta)
      // extract id of layer/disk; barrel [0->6], endcap [7->16]

      double dphi;
      int id;

      // Coordinates of two sensors.
      const GeomDetUnit* det0 = theStackedGeometry_->idToDetUnit( stDetId, 0 );
      const GeomDetUnit* det1 = theStackedGeometry_->idToDetUnit( stDetId, 1 );
      
      double R0 = det0->position().perp();
      double R1 = det1->position().perp();
      double Z0 = det0->position().z();
      double Z1 = det1->position().z();

      double pitch;

      if (theStackedGeometry_->isPSModule(stDetId)) {
      	pitch = 0.010;
      } else {
      	pitch = 0.009;
      }

      
      // Calculate the angle in phi between the stub direction and the normal to the sensor.
      if (stDetId.isBarrel()) {
      	double delta = R1-R0;
      	dphi = ( dir * pitch ) / delta ;
      	id = stDetId.iLayer();
      } 
      else {
      	double delta = (Z1-Z0) * R0 / Z0;
      	dphi = ( dir * pitch ) / delta ;
      	id = (theStackedGeometry_->isPSModule(stDetId)) ? (stDetId.iSide()-1) * 16 + (stDetId.iDisk() + 11) : (stDetId.iSide()-1) * 16 + (stDetId.iDisk() + 6);
      }

      // Estimate phi coordinate of track at beam
      double beta = pos.phi() + dphi;

      // Determine azimuthal segment this belongs to.
      double beta_bin = (beta + 3.1) / beta_range_ ;

      int bin_segment1 = 2 * int(beta_bin) + 1;
      int bin_segment2 = ( ( beta_bin - int(beta_bin) ) < 0.5 ) ? bin_segment1 - 1 : bin_segment1 + 1;

      if ( (bin_segment1 > 63) || (bin_segment2 > 63) || (bin_segment1 < 0) || (bin_segment2 < 0) ) {
	     	std::cout << "Warning: Beta segment out of range!" << std::endl;
      	continue;
      }
      double beta_segment1 = 0.5 * beta_range_ * (bin_segment1 - 1) - 3.1;
      double beta_segment2 = 0.5 * beta_range_ * (bin_segment2 - 1) - 3.1;

      
      //std::cout << "Beta segment 1 " << beta_segment1 << std::endl;
      //std::cout << "Beta segment 2 " << beta_segment2 << std::endl;


      if ( ( pt > pt_cut_ ) && ( fabs(pos.eta()) < eta_cut_ ) ) {
	
      	counter++;
      	extended_stub theStub;

      	theStub.barcode = counter;
      	theStub.stub = *stub;
      	theStub.ref = edmNew::makeRefTo(TTStubHandle_, stub );
      	theStub.position = pos;
      	theStub.eta = pos.eta();
      	theStub.phi = pos.phi();
      	theStub.r = pos.perp();
      	theStub.z = pos.z();
      	theStub.isReal = genuine;
      	theStub.dphi = dphi;
      	theStub.id = id;
      	theStub.pt = pt;
      	theStub.measured_pt = 0.;
      	
        unsigned int region = 0;
	
      	for ( ; region<(eta_regions_.size()-1); region++) {
      	  if(debug_) std::cout << "eta region n. "<<region << std::endl;
      	  
      	  double etaBound = eta_regions_[region];
      	  double etaSlice = eta_regions_[region+1] - etaBound;

      	  double maxR = 110;
      	  
      	  double modMinZ = std::min(Z0,Z1);
      	  double modMaxZ = std::max(Z0,Z1);
      	  double modMinR = std::min(R0,R1);
      	  double modMaxR = std::max(R0,R1);
      	  
      	  // z coordinates at outer radius of tracker of the two ends of this rapidity slice, assuming track comes from centre of CMS.
      	  double etaSliceZ1 = maxR / tan( 2 * atan(exp(-etaBound)) );
      	  double etaSliceZ2 = maxR / tan( 2 * atan(exp(-(etaBound + etaSlice))) );
      	  
      	  // If a track coming from anywhere inside the beam-spot passing through this module would cross the tracker outer radius
      	  // within the just calculated z boundaries, then module is in this rapidity segment.
      	  double etaDist1 =  modMaxZ - ((etaSliceZ1 >= 0 ? modMinR : modMaxR)*(etaSliceZ1 + zError_)/maxR - zError_);
      	  double etaDist2 = -modMinZ + ((etaSliceZ2 >= 0 ? modMaxR : modMinR)*(etaSliceZ2 - zError_)/maxR + zError_);
      	 
                // sensor belongs to this rapidity region.
      	  if( (etaDist1 > 0) && (etaDist2 > 0) ) {

      	    // Note angular region of stub in both rapidity & azimuth, and store all stubs in each specific region.
      	    Segment s1(region,beta_segment1);
      	    Segment s2(region,beta_segment2);
      	    
      	    sortedStubs_[s1].push_back(theStub); 
      	    sortedStubs_[s2].push_back(theStub);
            
      	  }
        }
      }
    }
  }
  std::cout<< "counter " << counter <<std::endl;
}



void TMTrackFinderAlgorithm::rphi_transform(Stubs& stubs, HoughArray& array, double beta_seg)
{

  // Fill the r-phi Hough-Transform array with the stubs.
 
  Stubs::iterator st = stubs.begin();
    
  for ( ; st!=stubs.end(); st++) {
   
    extended_stub stub = *st;
   
    // 0.006 is an approximation to B*c/2e11.
    const double Bc = 0.006;
    double hough_x_min = -Bc / hough_pt_;
    double hough_x_max = Bc / hough_pt_;
    int hough_x_bins = hough_bins_x_;

    double hough_y_min = 0.0;
    double hough_y_max = beta_range_;
    int hough_y_bins = hough_bins_y_;


    // Loop over bins related to 1/Pt.
    for (int j=0; j<hough_x_bins; j++) {

      // size of one bin, and range of bin j.
      double m_scale = (hough_x_max - hough_x_min) / (1.0 * hough_x_bins);

      double m_min = (m_scale * j) + hough_x_min - 0.5*m_scale;
      double m_max = (m_scale * j) + hough_x_min + 0.5*m_scale;

      double c_scale = (hough_y_max - hough_y_min) / (1.0 * hough_y_bins);

      // Calculate the range in track phi0 (=c) that would be consistent with this stub
      // and varying the Pt in a range within this bin (in m).
      double c_min = stub.phi - beta_seg - stub.r * m_max;
      double c_max = stub.phi - beta_seg - stub.r * m_min;

      if ( (c_max>hough_y_min) && (c_min<hough_y_max) ) {
       	int nbins = (int) std::ceil( (c_max - c_min) / c_scale ) ; 

      	for (int bin=0; bin<nbins; bin++) {
          double k = (c_min - hough_y_min) / c_scale + 1.0 * bin;
	  
	        if ( (k>=0) && (k<hough_y_bins) ) {
      	    stub.measured_pt = -Bc/(j*m_scale + hough_x_min)  ;
      	    stub.pt = -Bc/(j*m_scale + hough_x_min)  ;

      	    std::pair<int,int> index(j,abs(k));
      	    array[index].push_back(stub);
      	  }
	    	}
      }
    }
  }
}

void TMTrackFinderAlgorithm::cell_readout(HoughArray& array, Cells& accepted_cells)
{

  // Identify HT cells with a number of stubs exceeding a cut and store the stubs in each of these
  // accepted cells.

  HoughArray::iterator cell = array.begin();

  for ( ; cell!=array.end(); cell++) {
    
    int j = cell->first.first;
    int k = cell->first.second;

    Stubs stubs = cell->second;

    std::vector<int> r_register(40,0);

    Stubs::iterator stub = stubs.begin();
   
    int radii = nradii_;  // Define cut on Min. no. of required stubs per HT cell.
    for ( ; stub!=stubs.end(); stub++) {
      //int r_bin = int((stub->r - 20.0) / 5.0);
      //std::cout << "layer id "<< stub->id << " stub eta" << stub->eta<<  std::endl;
      if(broken_layer_==false || (broken_layer_ == true && stub->id != broken_layer_id_)){
        r_register[stub->id] = 1; // Note which r bins this HT cell has stubs in.
        //std::cout << "layer id = " << stub->id << std::endl;
      }
      if(eta_register_ && nradii_==5 && ((/*stub->eta>-1.1 &&*/ stub->eta < -0.9)||(/*stub->eta<1.1 &&*/ stub->eta > 0.9) ))
        radii = 4; // Looser cut in barrel/endcap transititon region.
    }
    if(debug_) std::cout << "radii = " << radii<< std::endl;

    // Check if number of stubs in HT cell exceeds cut.
    bool accept;
    int r_bins = std::accumulate(r_register.begin(),r_register.end(),0);
    accept = ( r_bins >= radii ) ? true : false ;
    if (accept) {
      accepted_cells.push_back(stubs); // store the stubs in this HT accepted cell.
      if (hough_print_) {
      	std::cout << "*********************" << std::endl;
      	std::cout << "m " << j << " : c " << k << "n. Stubs"<< stubs.size()<< std::endl;
      	if(debug_){
      	  stub = stubs.begin();
      	  for ( ; stub!=stubs.end(); stub++) {
      	    std::cout<<std::setprecision(3)
      		     <<stub->eta<<"\t"
      		     <<stub->phi<<"\t"
      		     <<stub->r<<"\t"
      		     <<stub->z<<"\t"
      		     <<stub->isReal<<"\t"
      		     <<stub->dphi<<"\t"
      		     <<stub->id<<"\t"
      		     <<stub->pt<<"\t"
      		     <<stub->measured_pt<<"\t"
      		     <<(stub->z)/(stub->r)<<std::endl;
      	  }
      	}
      }
    }
  }
}



void TMTrackFinderAlgorithm::eta_filter(Cells& cells, Cells& filtered_cells)
{

  // Like function cell_readout(), but also requires that the stubs should have consistent rapidity
  // to be stored and counted.
  // Loop over stubs in accepted HT cells.
  Cells::iterator stubs = cells.begin();
  for ( ; stubs!=cells.end(); stubs++) {
    // Histograms the eta of the stubs in each HT cell.
    std::vector<int> eta_hist(64,0);
    Stubs::iterator stub = stubs->begin();
    for ( ; stub!=stubs->end(); stub++) {
    	int eta_bin = int((stub->eta + 3.2) / 0.1);
    	if ( ( eta_bin < 0 ) || ( eta_bin > 63) ) {
        std::cout << "Warning: Filter bin " << eta_bin << " out of range!" << std::endl; 
      }
      else {
        eta_hist[eta_bin]++; 
      }
    }
   
    int modebin = -1;
    int maxbin = 0;

    for (int i=0; i<64; i++) {
  	
    	if (eta_hist[i] > maxbin) {
    	  modebin=i;
    	  maxbin=eta_hist[i];
    	}

    }

    if ( modebin==-1 ) std::cout << "Warning: Could not find filter mode bin!" << std::endl; 


    Stubs tempstubs;
    double meaneta = modebin * 0.1 - 3.2;
    double deltaeta = (-0.0775 * fabs(meaneta)) + 0.35; // assumed resolution in eta

    stub = stubs->begin();

    // Count stubs with consistent eta 

    for ( ; stub!=stubs->end(); stub++) {
    	if ( (stub->eta > (meaneta - deltaeta)) && (stub->eta < (meaneta + deltaeta)) ) 
    	  tempstubs.push_back(*stub);
    }
    if (tempstubs.size()>=(unsigned int)nradii_) {
    	filtered_cells.push_back(tempstubs);

    	if (debug_) {
    	  
    	  std::cout << "==============================" << std::endl;
    	  std::cout << "mode eta " << meaneta << " : delta eta " << deltaeta << std::endl;
    	  
    	  stub = tempstubs.begin();
    	  for ( ; stub!=tempstubs.end(); stub++) {
    	    
    	    std::cout<<std::setprecision(3)
    		     <<stub->eta<<"\t"
    		     <<stub->phi<<"\t"
    		     <<stub->r<<"\t"
    		     <<stub->z<<"\t"
    		     <<stub->isReal<<"\t"
    		     <<stub->dphi<<"\t"
    		     <<stub->id<<"\t"
    		     <<stub->pt<<"\t"
    		     <<(stub->z)/(stub->r)<<std::endl;
    	    
    	  }
    	}
    }
  }
}



void TMTrackFinderAlgorithm::bend_filter(Cells& cells, Cells& filtered_cells)
{  
  
  Cells::iterator cell = cells.begin();

  for ( ; cell!=cells.end(); cell++) {

    Stubs stubsPlus;
    Stubs stubsMinus;
    Stubs stubsZero;
         
    Stubs::iterator stub = cell->begin();

    for ( ; stub!=cell->end(); stub++) {
      stub->dphi>0 ? stubsPlus.push_back(*stub) : (stub->dphi<0 ? stubsMinus.push_back(*stub) : stubsZero.push_back(*stub) ) ;
    }

    if (stubsPlus.size()==stubsMinus.size()) continue;

    Stubs stubsSameSign = stubsPlus.size() > stubsMinus.size() ? stubsPlus : stubsMinus ;
    stubsSameSign.insert(stubsSameSign.end(),stubsZero.begin(),stubsZero.end());

    if (stubsSameSign.size()>=5) {
      filtered_cells.push_back(stubsSameSign);
    }
  }
}



void TMTrackFinderAlgorithm::r_filter(Cells& cells, Cells& filtered_cells)
{  

  Cells::iterator cell = cells.begin();
  for ( ; cell!=cells.end(); cell++) {
    std::vector<int> r_register(40,0);
    Stubs::iterator stub = cell->begin();

    for ( ; stub!=cell->end(); stub++) {
      r_register[stub->id] = 1;
    }

    if (std::accumulate(r_register.begin(),r_register.end(),0) >= nradii_) {
      filtered_cells.push_back((*cell));
    }
  }
}
