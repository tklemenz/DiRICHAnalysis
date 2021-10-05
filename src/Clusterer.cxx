#include "Clusterer.h"
#include "Utility.h"

#include <vector> 
#include <algorithm>

ClassImp(Clusterer);

static const std::string process = "[CLUSTERER]";

static Double_t cluster_fib_range  = 1;   // The allowed range in fiber distance for cluster building
static Double_t cluster_time_range = 5;  // The allowed time window for cluster building [ns], since the calibrated data is already in ns
static Double_t tot_uppercut       = 23;  // Make sure there is no total BS
static Double_t tot_lowercut       = 6;   // Get rid of noise
static Int_t completeBSThreshold   = 500;

void Clusterer::findClusters(CTSEvent& event, const ParticleType& particleType, const bool& debug)
{
  if (debug) { printf("%s%s %sStart of Event %lu.%s\n",text::LBLU, process.c_str(), text::LCYN, mEventCounter, text::RESET); }

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
    Int_t maxToTIndex = -1;
    Int_t goodToTInLayer = 0;
    std::vector<Double_t> goodToTs{};

    /// The assumption for cosmics data is that per event and layer there can only one proper cluster and not more
    /// This assumption need to be changed. Lots of background radiation in the lab so more than one cluster can be within one layer
    for (auto& layer : signalsVec) {                                                        /// loop over signals layer by layer
      if (layer.size() == 0) { 
        layerCounter++;
        continue;
      }
      for (auto& signal : layer) {
        if (signal.getToT() >= tot_lowercut && signal.getToT() <= tot_uppercut) { goodToTInLayer++; goodToTs.emplace_back(signal.getToT()); }
      }
      if (debug) {
        printf("%s%s %sNumber of good signals in layer %d: %d. Good means %g <= ToT <= %g.%s\n",text::LBLU, process.c_str(), text::YEL, layerCounter+1, goodToTInLayer, tot_lowercut, tot_uppercut, text::RESET);
        printf("%s%s %sNumber of      signals in layer %d: %lu.%s\n",text::LBLU, process.c_str(), text::YEL, layerCounter+1, layer.size(), text::RESET);
        printf("%s%s %sThe ToT values are:%s\n",text::LBLU, process.c_str(), text::GRN, text::RESET);
        for (auto& signal : layer) {
          printf("%s%s %s%g%s\n",text::LBLU, process.c_str(), text::GRN, signal.getToT(), text::RESET);
        }
      }
      while (goodToTInLayer > 0) {                                                          /// look for clusters with ToT max > tot_lowercut
        if (debug) {
          printf("%s%s %sLooking for Cluster... %d good signal(s) to be found...%s\n",text::LBLU, process.c_str(), text::CYN, goodToTInLayer, text::RESET);
         }
        Double_t maxToT = 0;
        ULong_t index = -1;
        ULong_t indexCounter = 0;
        for (auto& signal : layer) {                                                        /// Look for largest ToT in layer
          if (signal.isUsedInCluster()) {
            if (debug) { printf("%s%s %sSkipping used ToT %g%s\n",text::LBLU, process.c_str(), text::CYN, signal.getToT(), text::RESET); }
            indexCounter++;
            continue;
          }
          if (signal.getToT() > tot_uppercut) {
            indexCounter++;
            continue; }
          if (signal.getToT() > maxToT) {
            maxToT = signal.getToT();
            index = indexCounter;
          }
          indexCounter++;
        }

        /// remove this line after bug is found!!
        if (index == -1) {
          if (debug) { printf("%s%s %s%sSKIPPED A CLUSTER THAT SHOULD NOT BE SKIPPED!!%s\n",text::LBLU, process.c_str(), text::BLINK, text::LRED, text::RESET); }
        }

        clustersVec.at(layerCounter).emplace_back(Cluster());                               /// create a cluster for it
        clustersVec.at(layerCounter).back().addSignal(layer.at(index));                     /// and add the signal to the cluster
        layer.at(index).setIsUsedInCluster();
        goodToTInLayer--;

        Int_t maxFiber = mapping::getFiberNr(layer.at(index).getConfiguration(), layer.at(index).getChannelID(), layer.at(index).getTDCID());
        Double_t maxTime  = layer.at(index).getTimeStamp();

        if (debug) { printf("%s%s %sMax ToT Cluster found. ToT: %g%s\n",text::LBLU, process.c_str(), text::CYN, layer.at(index).getToT(), text::RESET); }

        for (auto& signal : layer) {                                                       /// add fitting signals to the cluster
          if (signal.isUsedInCluster()) { continue; }

          if (std::abs(mapping::getFiberNr(signal.getConfiguration(), signal.getChannelID(), signal.getTDCID()) - maxFiber) <= cluster_fib_range) {
            if (std::abs(signal.getTimeStamp() - maxTime) < cluster_time_range) {
              clustersVec.at(layerCounter).back().addSignal(signal);
              signal.setIsUsedInCluster();
              if (signal.getToT() >= tot_lowercut) { goodToTInLayer--; }
            }
          }
        }
        index = -1;
        if (debug) { printf("%s%s %s%d good signal(s) left...%s\n",text::LBLU, process.c_str(), text::CYN, goodToTInLayer, text::RESET); }
      } // while good signals can be found
      layerCounter++;
      goodToTInLayer = 0;
      goodToTs.clear();
    } // loop over layers

    for (auto& layerClus : clustersVec) {
      for (auto& cluster : layerClus) {
        mClusterVec.emplace_back(std::move(cluster));
      }
    }

    hNSignalsEvent->Fill(signalCounter);
    hNGoodSignalsEvent->Fill(goodSignalCounter);

    if (debug) { printf("%s%s %sEnd of Event %lu.%s\n\n\n",text::LBLU, process.c_str(), text::LCYN, mEventCounter, text::RESET); }
    mEventCounter++;
  } // ParticleType::Cosmic

  else if (particleType == ParticleType::Unknown) { return; } /// Unknown particle type 

  else {
    printf("Please provide proper ParticleType!\n");
    return;
  }

  return;
}
