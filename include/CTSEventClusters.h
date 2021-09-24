#ifndef CTSEVENTCLUSTERS_H
#define CTSEVENTCLUSTERS_H

#include "EventBase.h"
#include "Cluster.h"

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

  inline void setClusters (std::vector<Cluster> &clusters)  { mClusterVec = clusters; }

  std::vector<Cluster>&       getClusters()       { return mClusterVec; }
  const std::vector<Cluster>& getClusters() const { return mClusterVec; }

 private:
  std::vector<Cluster> mClusterVec {};  ///< contains all clusters in the event

  ClassDef(CTSEventClusters,1);
};

#endif