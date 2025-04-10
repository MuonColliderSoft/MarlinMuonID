#pragma once

#include "marlin/Processor.h"

#include "gsl/gsl_rng.h"

#include "TFormula.h"
#include "TH1F.h"
#include "TH2F.h"

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

  static constexpr float _lightSpeed = 299.792458; // [mm/ns]
  static constexpr float _muonMass = 0.1056583745; // [GeV]
  static constexpr int _muonPDG = 13;

  static constexpr unsigned int _muonDetBarrelLayers = 7;
  static constexpr unsigned int _muonDetEndcapLayers = 6;

  gsl_rng* _rng {nullptr};

  unsigned int _muonDetBarrel;
  unsigned int _muonDetEndcap;

  float _ecalB_inner_r;
  float _ecalE_min_z;

  float _bField = 5.;

  float _pt_min = 1.;  // [GeV]
  float _d0_max = 0.1; // [mm]
  float _z0_max = 0.1; // [mm]

  std::vector<float> _xyzResolution = {10., 10., 10.}; // [mm]
  float _timeResolution = 0.1; // [ns]

  TFormula* _tof_correction = nullptr;
  std::string _formulaStr = "(x < 0.5735 || x > 2.57) ? 0.015058 : "
    "(x >= 0.5735 && x < 0.626) ? (-9.63278 + 16.8561*x) : "
    "(x >= 0.626 && x < 2.516) ? (2.63371 - 4.14427*x + 3.42931*x*x - 1.3498*x*x*x + 0.216549*x*x*x*x) : "
    "(x >= 2.516 && x < 2.57) ? (43.5031 - 16.9276*x) : 0";

  std::vector<float> _deltaRMatch = { 0.2, 0.3}; // 0 --> barrel, 1 --> endcaps
  std::vector<float> _timeWindowB = {-0.3, 0.3}; // [ns]
  std::vector<float> _timeWindowE = {-0.3, 0.3}; // [ns]

  int _nHitsMatch = 4;

  bool _fillHistos = false;

  TH1D* _hDeltaR_Barrel;
  TH1D* _hDeltaR_Endcap;
  TH1D* _hDeltaT_Barrel;
  TH1D* _hDeltaT_Endcap;
  std::vector<TH1D*> _hDeltaR_B;
  std::vector<TH1D*> _hDeltaT_B;
  std::vector<TH1D*> _hDeltaR_E;
  std::vector<TH1D*> _hDeltaT_E;
  TH2D* _hDeltaR_vs_Pt;
  TH2D* _hDeltaT_vs_Pt;
  TH2D* _hDeltaR_vs_Theta;
  TH2D* _hDeltaT_vs_Theta;
  TH1D* _hNhits;
  
  //! Input track collection
  std::string _inputTrackCollection {};

  //! Input muon detector hit collection
  std::string _inputMuonHitCollection {};

  //! Output muon collection
  std::string _outputMuonCollection {};
  
};
