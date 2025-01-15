#include "include/MarlinMuonID.h"

#include <DD4hep/Detector.h>

#include <EVENT/LCCollection.h>
#include <EVENT/Track.h>
#include <EVENT/CalorimeterHit.h>

#include <UTIL/CellIDDecoder.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/ReconstructedParticleImpl.h>

#include <marlin/AIDAProcessor.h>

#include <math.h>


MarlinMuonID aMarlinMuonID ;

MarlinMuonID::MarlinMuonID()
  : Processor("MarlinMuonID") {

  _description = "MarlinMuonID performs muon ID by matching tracks to muon detector hits." ;

  // --- Register the steering parameters
  registerInputCollection( LCIO::TRACK,
		  	   "InputTrackCollection",
			   "Track input collection",
			   _inputTrackCollection,
			   _inputTrackCollection
			   );

  registerInputCollection( LCIO::CALORIMETERHIT,
		  	   "InputMuonHitCollection",
			   "Muon hit input collection",
			   _inputMuonHitCollection,
			   _inputMuonHitCollection
			   );

  registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			    "OutputMuonCollection", 
			    "Muon output collection", 
			    _outputMuonCollection,
			    _outputMuonCollection 
			    );

  registerProcessorParameter( "DeltaRMatch",
			      "DeltaR for matching tracks with muon detector hits",
			      _deltaRMatch,
			      _deltaRMatch
			      );

  registerProcessorParameter( "NHitsMatch",
			      "Minumum number of matching hits in the muon detector",
			      _nHitsMatch,
			      _nHitsMatch
			      );

  registerProcessorParameter( "FillHistograms",
			      "Fill the diagnostic histograms",
			      _fillHistos,
			      _fillHistos );

}

void MarlinMuonID::init() {

  dd4hep::Detector& theDetector = dd4hep::Detector::getInstance();

  // --- Get the muon detector IDs
  _muonDetBarrel = theDetector.constant<unsigned int>("DetID_Yoke_Barrel");
  _muonDetEndcap = theDetector.constant<unsigned int>("DetID_Yoke_Endcap");
  
  // --- Get the ECAL barrel inner radius and the ECAL endcap minimum z
  _ecalB_inner_r = theDetector.constant<float>("ECalBarrel_inner_radius")/dd4hep::mm;
  _ecalE_min_z = theDetector.constant<float>("ECalEndcap_min_z")/dd4hep::mm;
  
  // --- Get the magnetic field value
  const double pos[3] = {0., 0., 0.}; 
  double bFieldVec[3] = {0., 0., 0.}; 
  theDetector.field().magneticField(pos,bFieldVec);
  _bField = bFieldVec[2]/dd4hep::tesla;

  // --- Book the diagnostic histograms
  if ( _fillHistos ) {
    marlin::AIDAProcessor::histogramFactory(this);
    _hDeltaR_Barrel = new TH1F("hDeltaR_Barrel", "#DeltaR between tracks and muon detector hits (barrel);#DeltaR [rad]",
			       100, 0., 0.5);
    _hDeltaR_Endcap = new TH1F("hDeltaR_Endcap", "#DeltaR between tracks and muon detector hits (endcap);#DeltaR [rad]",
			       100, 0., 0.5);

    for (unsigned int ih=0; ih<_muonDetBarrelLayers; ++ih){
      std::string hname = "hDeltaR_B" + std::to_string(ih);
      std::string htitle = "#DeltaR between tracks and muon detector hits (layer B" + std::to_string(ih) + ");#DeltaR [rad]";
      _hDeltaR_B.push_back(new TH1F(hname.c_str(), htitle.c_str(), 100, 0., 0.5));
    }

    for (unsigned int ih=0; ih<_muonDetEndcapLayers; ++ih){
      std::string hname = "hDeltaR_E" + std::to_string(ih);
      std::string htitle = "#DeltaR between tracks and muon detector hits (layer E" + std::to_string(ih) + ");#DeltaR [rad]";
      _hDeltaR_E.push_back(new TH1F(hname.c_str(), htitle.c_str(), 100, 0., 0.5));
    }

    _hDeltaR_vs_Pt = new TH2F("hDeltaR_vs_Pt", "#DeltaR vs p_{T}^{trk};p_{T}^{trk} [GeV];#DeltaR [rad]",
				 100, 0., 100., 100, 0., 0.5);
    _hDeltaR_vs_Theta = new TH2F("hDeltaR_vs_Theta", "#DeltaR vs #theta_{trk};#theta_{trk} [rad];#DeltaR [rad]",
				 100, 0., M_PI, 100, 0., 0.5);

    _hNhits = new TH1F("hNhits", "Number of muon detector hits matched to tracks;N_{hits}", 20, 0., 20.);
  }
    
  // --- Print the initial parameters
  printParameters() ;

}

void MarlinMuonID::processRunHeader( LCRunHeader* /*run*/) {}

void MarlinMuonID::processEvent( LCEvent * evt ) {

  // --- Get the track collection
  LCCollection* InputTrackCollection = evt->getCollection(_inputTrackCollection);
  if( InputTrackCollection->getTypeName() != lcio::LCIO::TRACK )
    { throw EVENT::Exception( "Invalid collection type: " + InputTrackCollection->getTypeName() ) ; }

  // --- Get the muon detector hit collection
  LCCollection* InputMuonHitCollection = evt->getCollection(_inputMuonHitCollection);
  if( InputMuonHitCollection->getTypeName() != lcio::LCIO::CALORIMETERHIT )
    { throw EVENT::Exception( "Invalid collection type: " + InputMuonHitCollection->getTypeName() ) ; }

  // --- Make the output muon collection
  LCCollectionVec* OutputMuonCollection = new IMPL::LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);

  // --- Define the cell ID decoder
  UTIL::CellIDDecoder<lcio::CalorimeterHit> decoder(InputMuonHitCollection);


  // --- Loop over the reconstructed tracks
  for(int itrk=0; itrk<InputTrackCollection->getNumberOfElements(); ++itrk) {

    EVENT::Track* trk = static_cast<EVENT::Track*>(InputTrackCollection->getElementAt(itrk));

    // --- Find the track intersection point on the calorimeter inner surface
    const lcio::TrackState* ts_atCAL = trk->getTrackState( lcio::TrackState::AtCalorimeter );
    if ( !ts_atCAL ){
      streamlog_out(WARNING) << " MarlinMuonID::processEvent(): skipped track #" << itrk
			     << ", no track state at calorimeter" << std::endl;
      continue;
    }
    
    const float trk_cotTheta = ts_atCAL->getTanLambda();
    const float trk_pt = 0.000299792458 * _bField / fabs(ts_atCAL->getOmega());
    const float trk_phi = ts_atCAL->getPhi();
    const float trk_theta = M_PI_2 - std::atan(trk_cotTheta);

    float x_ecal, y_ecal, z_ecal;
    z_ecal = _ecalB_inner_r * trk_cotTheta;
    if ( fabs(z_ecal) > _ecalE_min_z ){
      z_ecal = ( z_ecal>0. ? _ecalE_min_z : -_ecalE_min_z );
      float r = z_ecal / trk_cotTheta;
      x_ecal = r * std::cos(trk_phi);
      y_ecal = r * std::sin(trk_phi);
    }
    else {
      x_ecal = _ecalB_inner_r * std::cos(trk_phi);
      y_ecal = _ecalB_inner_r * std::sin(trk_phi);
    }

    // --- Loop over the muon detector hits
    int n_matchedHits = 0;
    for(int ihit=0; ihit<InputMuonHitCollection->getNumberOfElements(); ++ihit) {

      EVENT::CalorimeterHit* hit = static_cast<EVENT::CalorimeterHit*>(InputMuonHitCollection->getElementAt(ihit));

      unsigned int system = decoder(hit)["system"];
      unsigned int layer = decoder(hit)["layer"];
      
      float x_centered = hit->getPosition()[0] - x_ecal;
      float y_centered = hit->getPosition()[1] - y_ecal;
      float z_centered = hit->getPosition()[2] - z_ecal;

      float hit_r = std::sqrt(x_centered*x_centered + y_centered*y_centered);
      float hit_phi = std::atan2(y_centered, x_centered); 
      float hit_theta = std::atan2(hit_r, z_centered);

      float deltaR = std::sqrt( (hit_phi-trk_phi)*(hit_phi-trk_phi) + (hit_theta-trk_theta)*(hit_theta-trk_theta) );

      if ( _fillHistos ){
	if ( system == _muonDetBarrel){
	  _hDeltaR_Barrel->Fill(deltaR);
	  _hDeltaR_B[layer]->Fill(deltaR);
	}
	else if ( system == _muonDetEndcap){
	  _hDeltaR_Endcap->Fill(deltaR);
	  _hDeltaR_E[layer]->Fill(deltaR);
	}
	else {
	  continue;
	}
	_hDeltaR_vs_Pt->Fill(trk_pt,deltaR);
	_hDeltaR_vs_Theta->Fill(trk_theta,deltaR);
      }

      // --- Check if the hit is matched to the track
      if ( deltaR < _deltaRMatch ) 
	n_matchedHits++;
	
    } // ihit loop

    if ( _fillHistos ){
      _hNhits->Fill(n_matchedHits);
    }
      
    // --- Declare the track a muon if more than _nHitsMatch hits of muon detectors are associated to it
    if ( n_matchedHits > _nHitsMatch ){

      // --- Create a ReconstructedParticle to store the identified muon
      IMPL::ReconstructedParticleImpl *const muon(new ReconstructedParticleImpl());

      const float charge = ( ts_atCAL->getOmega()>0. ? 1. : -1 );
      const int pdg = ( charge < 0. ? _muonPDG : -_muonPDG ); 
      const float momentum[3] = { trk_pt * std::cos(trk_phi),
	                          trk_pt * std::sin(trk_phi),
				  trk_pt * std::cos(trk_theta) };
      muon->setMomentum(momentum);
      muon->setEnergy(std::sqrt(momentum[0]*momentum[0] + momentum[1]*momentum[1] +
				momentum[2]*momentum[2] + _muonMass*_muonMass));
      muon->setMass(_muonMass);
      muon->setCharge(charge);
      muon->setType(pdg);

      muon->addTrack(trk);

      OutputMuonCollection->addElement(muon);

    }

  } // itrk loop
  
  // --- Save the output muon collection
  evt->addCollection(OutputMuonCollection, _outputMuonCollection);  

}

void MarlinMuonID::check( LCEvent * /*evt*/ ){}

void MarlinMuonID::end(){}
