# MarlinMuonID
Marlin Muon ID Processor

This Marlin processor implements a simple muon identification algorithm, matching tracks to hits in the outer muon detectors:

- tracks are extrapolated to the electromagnetic calorimeter inner surface;
- a cone is opened along the particle's flight direction;
- hits within the cone are matched to the track.

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
    "DeltaRMatch": ["0.1"],
    "NHitsMatch": ["5"],
    "FillHistograms":["false"]
}
```
