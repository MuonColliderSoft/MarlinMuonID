#pragma once

#include <marlin/Processor.h>

#include <TH1F.h>
#include <TH2F.h>

/**
 * Performs muon identification by matching tracks to muon detector hits.
 *
 * @parameter InputTrackCollection Name of the track input collection
 * @parameter InputMuonHitCollection Name of the muon hit input collection 
 * @parameter OutputMuonCollection Name of the reconstructed muon output collection
 *
 * @parameter DeltaRMatch DeltaR for matching tracks with muon detector hits in barrel and endcaps
 * @parameter NHitsMatch Minumum number of matching hits in the muon detector
 *
 * @parameter FillHistograms Flag to enable diagnostic histograms
 *
 * @author M. Casarsa
 * @date  11 December 2024
 * @version $Id: $
 */
class MarlinMuonID : public marlin::Processor
{
public:
  virtual Processor*  newProcessor() { return new MarlinMuonID ; }

  MarlinMuonID(const MarlinMuonID &) = delete ;
  MarlinMuonID& operator =(const MarlinMuonID &) = delete ;
  MarlinMuonID() ;

  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;

  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;

  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent * evt ) ; 

  virtual void check( LCEvent * evt ) ; 

  /** Called after data processing for clean up.
   */
  virtual void end() ;  

private:

  static constexpr float _muonMass = 0.1056583745; // in GeV
  static constexpr int _muonPDG = 13;

  static constexpr unsigned int _muonDetBarrelLayers = 7;
  static constexpr unsigned int _muonDetEndcapLayers = 6;

  unsigned int _muonDetBarrel;
  unsigned int _muonDetEndcap;

  float _ecalB_inner_r;
  float _ecalE_min_z;

  float _bField = 5.;

  std::vector<float> _deltaRMatch = {0.2, 0.3}; // 0 --> barrel, 1 --> endcaps
  int _nHitsMatch = 4;

  bool _fillHistos = false;

  TH1F* _hDeltaR_Barrel;
  TH1F* _hDeltaR_Endcap;
  std::vector<TH1F*> _hDeltaR_B;
  std::vector<TH1F*> _hDeltaR_E;
  TH2F* _hDeltaR_vs_Pt;
  TH2F* _hDeltaR_vs_Theta;
  TH1F* _hNhits;
  
  //! Input track collection
  std::string _inputTrackCollection {};

  //! Input muon detector hit collection
  std::string _inputMuonHitCollection {};

  //! Output muon collection
  std::string _outputMuonCollection {};
  
};
