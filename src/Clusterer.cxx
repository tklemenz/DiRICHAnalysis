#include "Clusterer.h"
#include "Utility.h"

#include <vector> 
#include <algorithm>

ClassImp(Clusterer);

static Double_t cluster_fib_range  = 2;   // The allowed range in fiber distance for cluster building
static Double_t cluster_time_range = 2;   // The allowed time window for cluster building [ns], since the calibrated data is already in ns
static Double_t tot_uppercut       = 23;  // Make sure there is no total BS
static Double_t tot_lowercut       = 8;   // Get rid of noise
static Int_t completeBSThreshold   = 500;

void Clusterer::findClusters(CTSEvent& event, const ParticleType& particleType)
{
  Int_t fullTrackCounter = 0;
  // Data available within the CTS event
  std::vector<Module> modules{};
  for(Int_t module=0; module < event.getModules().size(); module++) {
    modules.emplace_back(event.getModules().at(module));
  }

  std::vector<Double_t> distances{};
  std::vector<Double_t> time_distances{};
  Double_t min_dist_index = -1;
  Double_t min_time_index = -1; 
  Int_t x = -1;
  Int_t y = -1;
  Int_t layer = -1;
  Int_t fiberNr = -1;

  ULong_t signalCounter     = 0;
  ULong_t goodSignalCounter = 0;


  if (particleType == ParticleType::Cosmic) {
    std::vector<std::vector<Signal>> signalsVec{};        /// intermediate storage for signals per layer
    std::vector<std::vector<Cluster>> clustersVec{};
    for (Int_t layer=0; layer<16; layer++) {
      signalsVec.emplace_back(std::vector<Signal>{});
      clustersVec.emplace_back(std::vector<Cluster>{});
    }

    /// Fill all signals in vectors
    /// One vector per layer
    /// Then one can loop layer by layer to find clusters
    for (auto& module : modules) {
      for (auto& fiber : module.getFibers()) {
        if (fiber.getNSignals() > 0) {
          layer = fiber.getLayer();
          for (auto& signal : fiber.getSignals()) {
            if (signal.getSignalNr() != 1) { continue; }
            signalCounter++;
            if (signal.getToT() == 0 || signal.getToT() < 0 || signal.getToT() > completeBSThreshold) { continue; }
            if (signal.getToT() >= tot_lowercut && signal.getToT() <= tot_uppercut) { goodSignalCounter++; }
            if (mapping::getModule(signal.getConfiguration(), signal.getTDCID()) == 1) { layer += 8; }
            signalsVec.at(layer-1).emplace_back(signal);
          }
        }
      }
    }

    Int_t layerCounter = 0;
    std::vector<Double_t> layerToT{};
    Int_t maxToTIndex = -1;

    /// The assumption for cosmics data is that per event and layer there can only one proper cluster and not more
    for (auto& layer : signalsVec) {                                                      /// loop over signals layer by layer
      if (layer.size() == 0) { continue; }
      for (auto& signal : layer) {
        if (signal.getToT() < tot_uppercut) {
          layerToT.emplace_back(signal.getToT());                                         /// fill all ToT values in the layer to the vector
        }
      }
      maxToTIndex = std::max_element(layerToT.begin(),layerToT.end())-layerToT.begin();   /// find the largest ToT in the layer
      //printf("max tot index: %d\t max ToT: %g\n",maxToTIndex, layer.at(maxToTIndex).getToT());
      clustersVec.at(layerCounter).emplace_back(Cluster());                               /// create a cluster for it
      clustersVec.at(layerCounter).back().addSignal(layer.at(maxToTIndex));               /// and add the signal to the cluster
      layer.at(maxToTIndex).setIsUsedInCluster();

      Int_t maxFiber = mapping::getFiberNr(layer.at(maxToTIndex).getConfiguration(), layer.at(maxToTIndex).getChannelID(), layer.at(maxToTIndex).getTDCID());
      Double_t maxTime  = layer.at(maxToTIndex).getTimeStamp();

      for (auto& signal : layer) {
        if (signal.isUsedInCluster()) { continue; }
        if (std::abs(mapping::getFiberNr(signal.getConfiguration(), signal.getChannelID(), signal.getTDCID()) - maxFiber) < cluster_fib_range) {
          if (std::abs(signal.getTimeStamp() - maxTime) < cluster_time_range) {
            clustersVec.at(layerCounter).back().addSignal(signal);
          }
        }
      }
      layerToT.clear();
      maxToTIndex = -1;
      layerCounter++;
    }

    for (auto& layerClus : clustersVec) {
      for (auto& cluster : layerClus) {
        mClusterVec.emplace_back(std::move(cluster));
      }
    }

    hNSignalsEvent->Fill(signalCounter);
    hNGoodSignalsEvent->Fill(goodSignalCounter);
  }

  else if (particleType == ParticleType::Unknown) { return; } /// Unknown particle type 

  else {
    printf("Please provide proper ParticleType!\n");
    return;
  }

  return;
}
