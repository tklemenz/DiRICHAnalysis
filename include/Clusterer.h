#ifndef CLUSTERER_H
#define CLUSTERER_H

#include "CTSEvent.h"
#include "Cluster.h"
#include <TH1.h>
#include <TH2.h>
#include <vector>

/// Clusterer class to assign signals from CTSEvents to clusters.
/// For ParticleTyple::Cosmic the clusterer assumes that there is not more than one cluster per layer, which is a fair assumption.

/// @todo implement proper clustering for the case of multiple modules

class Clusterer
{
 public:
  Clusterer() = default;
  ~Clusterer() = default;

  /// Use this if you already know which clusters/signals you want to have
  inline void setClusters (std::vector<Cluster> &clusterVec) { mClusterVec = clusterVec; }

  /// If it is not clear which signals should be taken into your clusters, use this
  /// This function loops through the signals in the given CTS Event,
  /// applies some cuts: allowed distance cut in position (fibers) and time (Work In Progress; Not Yet Implemented) 
  /// And gives back the vector of found clusters.
  void findClusters(CTSEvent& event, const ParticleType& particleType = ParticleType::Cosmic);

  std::vector<Cluster>&       getClusters()       { return mClusterVec; }
  const std::vector<Cluster>& getClusters() const { return mClusterVec; }

  /// Useful when calling this function iteratively in macro so that your cluster container doesn't grow indefinetely
  /// Try to guess how I found out that this was needed
  inline void reset() { mClusterVec.clear(); }

  /// number of signals per event
  /// does not count signalNr > 1
  TH1D* hNSignalsEvent     = new TH1D("hNSignalsEvent","n Signals in Event",200,0,200);

  /// number of good signals in the event
  /// nSignals from above with reasonable ToT cuts
  TH1D* hNGoodSignalsEvent = new TH1D("hNGoodSignalsEvent","n good Signals in Event",200,0,200);

 private:
  std::vector<Cluster> mClusterVec{}; //holds clusters

  ClassDef(Clusterer,2);
};

#endif