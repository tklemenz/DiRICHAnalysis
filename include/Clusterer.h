#ifndef CLUSTERER_H
#define CLUSTERER_H

#include "CTSEvent.h"
#include "Cluster.h"
#include <TH1.h>
#include <vector>

/// Clusterer class to assign signals from CTSEvents to clusters.

class Clusterer
{
 public:
  Clusterer() = default;
  ~Clusterer() = default;

  // Use this if you already know which clusters/signals you want to have
  inline void setClusters (std::vector<Cluster> &clusterVec) { mClusterVec = clusterVec; }

  // If it is not clear which signals should be taken into your clusters, use this
  // This function loops through the signals in the given CTS Event,
  // applies some cuts: allowed distance cut in position (fibers) and time (Work In Progress; Not Yet Implemented) 
  // And gives back the vector of found clusters.
  Int_t findClusters(CTSEvent& event, const ParticleType& particleType = ParticleType::Cosmic);

  std::vector<Cluster>&       getClusters()       { return mClusterVec; }
  const std::vector<Cluster>& getClusters() const { return mClusterVec; }

  // Useful when calling this function iteratively in macro so that your cluster container doesn't grow indefinetely
  // Try to guess how I found out that this was needed
  inline void reset() { mClusterVec.clear(); }

  TH1D* MhTimeDiff  = new TH1D("hTimeDiff","time difference of signal to closest cluster;timeDiff [ns]",500,0,50);
  TH1D* MhSpaceDiff = new TH1D("hSpaceDiff","space difference of signal to closest cluster;spaceDiff [fiber]",320,0,32);

 private:
  std::vector<Cluster> mClusterVec{}; //holds clusters 

  ClassDef(Clusterer,1);
};

#endif