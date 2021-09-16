#ifndef CTSEVENTCLUSTERS_H
#define CTSEVENTCLUSTERS_H

#include "EventBase.h"
#include "Clusterer.h"

/// This class represents a CTS event containing clusters.
/// It holds:
///          - a vector containing all clusters in the event

class CTSEventClusters : public EventBase
{
 public:
  CTSEventClusters() = default;
  ~CTSEventClusters() = default;
  CTSEventClusters(const CTSEventClusters &event);

  inline void addCluster (Cluster &cluster)  { mClusterVec.emplace_back(cluster); }
  inline void addCluster (std::vector<Cluster> &cluster)  { mClusterVec = cluster; }

  std::vector<Cluster>&       getClusters()       { return mClusterVec; }
  const std::vector<Cluster>& getClusters() const { return mClusterVec; }

 private:
  std::vector<Cluster> mClusterVec {};  ///< contains all clusters in the event

  ClassDef(CTSEventClusters,1);
};

#endif