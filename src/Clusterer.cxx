#include "Clusterer.h"
#include "Utility.h"

#include <vector> 
#include <algorithm>

ClassImp(Clusterer);

static Double_t cluster_fib_range  = 2;   // The allowed range in fiber distance for cluster building
static Double_t cluster_time_range = 2;  // The allowed time window for cluster building [ns], since the calibrated data is already in ns
static Double_t tot_uppercut       = 23;  // Make sure there is no total BS
static Double_t tot_lowercut       = 10;  // Get rid of noise

Int_t Clusterer::findClusters(CTSEvent& event, const ParticleType& particleType)
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
            if (signal.getToT() == 0 || signal.getToT() < 0) { continue; }
            layer = signal.getLayer();
            if (mapping::getModule(signal.getConfiguration(), signal.getTDCID()) == 2) { layer += 8; }
            signalsVec.at(layer-1).emplace_back(signal);
          }
        }
      }
    }

    if (signalsVec.at(0).size() > 0 && signalsVec.at(2).size() > 0 && signalsVec.at(4).size() > 0) { fullTrackCounter++; }  /// Count how many events have signals in layer 1, 3, 5
                                                                                                                            /// to see if the module is actually hit by cosmics

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
        //printf("Cluster size: %d\n", cluster.getNSignals());
        mHnSignal->Fill(cluster.getNSignals());
        mClusterVec.emplace_back(std::move(cluster));
      }
    }
  }


/*  else if (particleType == ParticleType::Unknown) {
  // Clusters can only contain signals from the same layer
  // Depending on the padiwa configuration a few of the vectors will always be empty
  std::vector<Cluster> layer1clusters{};
  std::vector<Cluster> layer2clusters{};
  std::vector<Cluster> layer3clusters{};
  std::vector<Cluster> layer4clusters{};
  std::vector<Cluster> layer5clusters{};
  std::vector<Cluster> layer6clusters{};
  std::vector<Cluster> layer7clusters{};
  std::vector<Cluster> layer8clusters{};

  for(auto& fiber : fibers) {
    if(fiber.getNSignals() > 0) {
      layer = fiber.getLayer();
      x     = fiber.getX();
      y     = fiber.getY();
      for(auto& signal : fiber.getSignals()) {
        if(signal.getSignalNr()==1) { 
          if(signal.getToT()<tot_lowercut || signal.getToT()>tot_uppercut) { continue; }  //allow only good signals
          // Since the signals are separeted in fibers but not layer,
          // this switch makes sure only signals from same layer are compared
          switch(layer){  
            case 1: if(layer1clusters.size()==0){
                      layer1clusters.emplace_back(Cluster());
                      layer1clusters.back().addSignal(signal);
                    }
                    else{
                      distances.clear();
                      min_dist_index=-1;
                      time_distances.clear();
                      min_time_index=-1;
                      for(auto& cluster : layer1clusters){ 
                        distances.emplace_back(std::abs(cluster.getMeanFiber()-(x==0 ? y : x)));
                        time_distances.emplace_back(std::abs(cluster.getMeanTimeStamp()-signal.getTimeStamp()));
                      }
                      min_dist_index = std::min_element(distances.begin(),distances.end())-distances.begin();
                      min_time_index = std::min_element(time_distances.begin(),time_distances.end())-time_distances.begin();
                      MhTimeDiff->Fill(time_distances.at(min_time_index));
                      MhSpaceDiff->Fill(distances.at(min_dist_index));
                      if(distances.at(min_dist_index)<cluster_fib_range && time_distances.at(min_time_index)<cluster_time_range/* && min_dist_index==min_time_index*//*){
                        layer1clusters.at(min_dist_index).addSignal(signal);
                      }
                      else{
                        layer1clusters.emplace_back(Cluster());
                        layer1clusters.back().addSignal(signal);
                      }
                    }
                    break;
            case 2: if(layer2clusters.size()==0){
                      layer2clusters.emplace_back(Cluster());
                      layer2clusters.back().addSignal(signal);
                    }
                    else{
                      distances.clear();
                      min_dist_index=-1;
                      time_distances.clear();
                      min_time_index=-1;
                      for(auto& cluster : layer2clusters){ 
                        distances.emplace_back(std::abs(cluster.getMeanFiber()-(x==0 ? y : x)));
                        time_distances.emplace_back(std::abs(cluster.getMeanTimeStamp()-signal.getTimeStamp()));
                      }
                      min_dist_index = std::min_element(distances.begin(),distances.end())-distances.begin();
                      min_time_index = std::min_element(time_distances.begin(),time_distances.end())-time_distances.begin();
                      MhTimeDiff->Fill(time_distances.at(min_time_index));
                      MhSpaceDiff->Fill(distances.at(min_dist_index));
                      if(distances.at(min_dist_index)<cluster_fib_range && time_distances.at(min_time_index)<cluster_time_range/* && min_dist_index==min_time_index*//*){
                        layer2clusters.at(min_dist_index).addSignal(signal);
                      }
                      else{
                        layer2clusters.emplace_back(Cluster());
                        layer2clusters.back().addSignal(signal);
                      }
                    }
                    break;
            case 3: if(layer3clusters.size()==0){
                      layer3clusters.emplace_back(Cluster());
                      layer3clusters.back().addSignal(signal);
                    }
                    else{
                      distances.clear();
                      min_dist_index=-1;
                      time_distances.clear();
                      min_time_index=-1;
                      for(auto& cluster : layer3clusters){ 
                        distances.emplace_back(std::abs(cluster.getMeanFiber()-(x==0 ? y : x)));
                        time_distances.emplace_back(std::abs(cluster.getMeanTimeStamp()-signal.getTimeStamp()));
                      }
                      min_dist_index = std::min_element(distances.begin(),distances.end())-distances.begin();
                      min_time_index = std::min_element(time_distances.begin(),time_distances.end())-time_distances.begin();
                      MhTimeDiff->Fill(time_distances.at(min_time_index));
                      MhSpaceDiff->Fill(distances.at(min_dist_index));
                      if(distances.at(min_dist_index)<cluster_fib_range && time_distances.at(min_time_index)<cluster_time_range/* && min_dist_index==min_time_index*//*){
                        layer3clusters.at(min_dist_index).addSignal(signal);
                      }
                      else{
                        layer3clusters.emplace_back(Cluster());
                        layer3clusters.back().addSignal(signal);
                      }
                    }
                    break;
            case 4: if(layer4clusters.size()==0){
                      layer4clusters.emplace_back(Cluster());
                      layer4clusters.back().addSignal(signal);
                    }
                    else{
                      distances.clear();
                      min_dist_index=-1;
                      time_distances.clear();
                      min_time_index=-1;
                      for(auto& cluster : layer4clusters){ 
                        distances.emplace_back(std::abs(cluster.getMeanFiber()-(x==0 ? y : x)));
                        time_distances.emplace_back(std::abs(cluster.getMeanTimeStamp()-signal.getTimeStamp()));
                      }
                      min_dist_index = std::min_element(distances.begin(),distances.end())-distances.begin();
                      min_time_index = std::min_element(time_distances.begin(),time_distances.end())-time_distances.begin();
                      MhTimeDiff->Fill(time_distances.at(min_time_index));
                      MhSpaceDiff->Fill(distances.at(min_dist_index));
                      if(distances.at(min_dist_index)<cluster_fib_range && time_distances.at(min_time_index)<cluster_time_range/* && min_dist_index==min_time_index*//*){
                        layer4clusters.at(min_dist_index).addSignal(signal);
                      }
                      else{
                        layer4clusters.emplace_back(Cluster());
                        layer4clusters.back().addSignal(signal);
                      }
                    }
                    break;
            case 5: if(layer5clusters.size()==0){
                      layer5clusters.emplace_back(Cluster());
                      layer5clusters.back().addSignal(signal);
                    }
                    else{
                      distances.clear();
                      min_dist_index=-1;
                      time_distances.clear();
                      min_time_index=-1;
                      for(auto& cluster : layer5clusters){ 
                        distances.emplace_back(std::abs(cluster.getMeanFiber()-(x==0 ? y : x)));
                        time_distances.emplace_back(std::abs(cluster.getMeanTimeStamp()-signal.getTimeStamp()));
                      }
                      min_dist_index = std::min_element(distances.begin(),distances.end())-distances.begin();
                      min_time_index = std::min_element(time_distances.begin(),time_distances.end())-time_distances.begin();
                      MhTimeDiff->Fill(time_distances.at(min_time_index));
                      MhSpaceDiff->Fill(distances.at(min_dist_index));
                      if(distances.at(min_dist_index)<cluster_fib_range && time_distances.at(min_time_index)<cluster_time_range/* && min_dist_index==min_time_index*//*){
                        layer5clusters.at(min_dist_index).addSignal(signal);
                      }
                      else{
                        layer5clusters.emplace_back(Cluster());
                        layer5clusters.back().addSignal(signal);
                      }
                    }
                    break;
            case 6: if(layer6clusters.size()==0){
                      layer6clusters.emplace_back(Cluster());
                      layer6clusters.back().addSignal(signal);
                    }
                    else{
                      distances.clear();
                      min_dist_index=-1;
                      time_distances.clear();
                      min_time_index=-1;
                      for(auto& cluster : layer6clusters){ 
                        distances.emplace_back(std::abs(cluster.getMeanFiber()-(x==0 ? y : x)));
                        time_distances.emplace_back(std::abs(cluster.getMeanTimeStamp()-signal.getTimeStamp()));
                      }
                      min_dist_index = std::min_element(distances.begin(),distances.end())-distances.begin();
                      min_time_index = std::min_element(time_distances.begin(),time_distances.end())-time_distances.begin();
                      MhTimeDiff->Fill(time_distances.at(min_time_index));
                      MhSpaceDiff->Fill(distances.at(min_dist_index));
                      if(distances.at(min_dist_index)<cluster_fib_range && time_distances.at(min_time_index)<cluster_time_range/* && min_dist_index==min_time_index*//*){
                        layer6clusters.at(min_dist_index).addSignal(signal);
                      }
                      else{
                        layer6clusters.emplace_back(Cluster());
                        layer6clusters.back().addSignal(signal);
                      }
                    }
                    break;
            case 7: if(layer7clusters.size()==0){
                      layer7clusters.emplace_back(Cluster());
                      layer7clusters.back().addSignal(signal);
                    }
                    else{
                      distances.clear();
                      min_dist_index=-1;
                      time_distances.clear();
                      min_time_index=-1;
                      for(auto& cluster : layer7clusters){ 
                        distances.emplace_back(std::abs(cluster.getMeanFiber()-(x==0 ? y : x)));
                        time_distances.emplace_back(std::abs(cluster.getMeanTimeStamp()-signal.getTimeStamp()));
                      }
                      min_dist_index = std::min_element(distances.begin(),distances.end())-distances.begin();
                      min_time_index = std::min_element(time_distances.begin(),time_distances.end())-time_distances.begin();
                      MhTimeDiff->Fill(time_distances.at(min_time_index));
                      MhSpaceDiff->Fill(distances.at(min_dist_index));
                      if(distances.at(min_dist_index)<cluster_fib_range && time_distances.at(min_time_index)<cluster_time_range/* && min_dist_index==min_time_index*//*){
                        layer7clusters.at(min_dist_index).addSignal(signal);
                      }
                      else{
                        layer7clusters.emplace_back(Cluster());
                        layer7clusters.back().addSignal(signal);
                      }
                    }
                    break;
            case 8: if(layer8clusters.size()==0){
                      layer8clusters.emplace_back(Cluster());
                      layer8clusters.back().addSignal(signal);
                    }
                    else{
                      distances.clear();
                      min_dist_index=-1;
                      time_distances.clear();
                      min_time_index=-1;
                      for(auto& cluster : layer8clusters){ 
                        distances.emplace_back(std::abs(cluster.getMeanFiber()-(x==0 ? y : x)));
                        time_distances.emplace_back(std::abs(cluster.getMeanTimeStamp()-signal.getTimeStamp()));
                      }
                      min_dist_index = std::min_element(distances.begin(),distances.end())-distances.begin();
                      min_time_index = std::min_element(time_distances.begin(),time_distances.end())-time_distances.begin();
                      MhTimeDiff->Fill(time_distances.at(min_time_index));
                      MhSpaceDiff->Fill(distances.at(min_dist_index));
                      if(distances.at(min_dist_index)<cluster_fib_range && time_distances.at(min_time_index)<cluster_time_range/* && min_dist_index==min_time_index*//*){
                        layer8clusters.at(min_dist_index).addSignal(signal);
                      }
                      else{
                        layer8clusters.emplace_back(Cluster());
                        layer8clusters.back().addSignal(signal);
                      }
                    }
                    break;
          }
        }
      } /// loop over signals in fiber
    }
  } /// loop over fibers in module

  // Add clusters to final vector if they are not empty
  if(layer1clusters.size()>0) { for(auto& cluster : layer1clusters) { mClusterVec.emplace_back(std::move(cluster)); } }
  if(layer2clusters.size()>0) { for(auto& cluster : layer2clusters) { mClusterVec.emplace_back(std::move(cluster)); } }
  if(layer3clusters.size()>0) { for(auto& cluster : layer3clusters) { mClusterVec.emplace_back(std::move(cluster)); } }
  if(layer4clusters.size()>0) { for(auto& cluster : layer4clusters) { mClusterVec.emplace_back(std::move(cluster)); } }
  if(layer5clusters.size()>0) { for(auto& cluster : layer5clusters) { mClusterVec.emplace_back(std::move(cluster)); } }
  if(layer6clusters.size()>0) { for(auto& cluster : layer6clusters) { mClusterVec.emplace_back(std::move(cluster)); } }
  if(layer7clusters.size()>0) { for(auto& cluster : layer7clusters) { mClusterVec.emplace_back(std::move(cluster)); } }
  if(layer8clusters.size()>0) { for(auto& cluster : layer8clusters) { mClusterVec.emplace_back(std::move(cluster)); } }
  ///-----------------------------------------

  } /// Unknown particle type */

  else {
    printf("Please provide proper ParticleType!\n");
    return -1;
  }

  return fullTrackCounter;
}
