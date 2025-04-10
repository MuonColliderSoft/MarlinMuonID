# MarlinMuonID
### Marlin Muon ID Processor

This Marlin processor implements a simple muon identification algorithm, matching tracks to hits in the outer muon detectors:

- tracks are extrapolated to the electromagnetic calorimeter inner surface;
- a cone is opened along the particle's flight direction;
- hits within the cone are matched to the track, if their time is compatible with the particle time of flight.

The identified muons are saved in a ReconstructedParticle collection.

The processor is configured as follows:
```
muonID = MarlinProcessorWrapper("MarlinMuonID")
muonID.OutputLevel = WARNING 
muonID.ProcessorType = "MarlinMuonID" 
muonID.Parameters = {
    "InputTrackCollection": ["SiTracks"],
    "InputMuonHitCollection": ["Muon_digi"],
    "OutputMuonCollection": ["RecoMuons"],
    "TrackPtMin": ["1"], # [GeV]
    "TrackD0Max": ["0.1"], # [mm]
    "TrackZ0Max": ["0.1"], # [mm]
    "MuonHitsSpatialResolution": ["10", "10", "10"], # [mm]
    "MuonHitsTimeResolution": ["0.1"], # [ns]
    "TrackTofCorrFunction": [], # string in ROOT's TFormula format (with no spaces!)
    "DeltaRMatch": ["0.2", "0.3"], # for barrel and endcaps
    "BarrelHitsTimeWindow": ["-0.3", "0.3"], # [ns]
    "EndcapHitsTimeWindow": ["-0.3", "0.3"], # [ns]
    "NHitsMatch": ["4"],
    "FillHistograms":["false"]
}
```
