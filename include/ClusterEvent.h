#ifndef CLUSTEREVENT_H
#define CLUSTEREVENT_H

#include "EventBase.h"
#include "CTSEventClusters.h"

/// This class represents a basic Cluster event containing clusters.
///
/// Before an event is written to file the removeEmpty funtion should be
/// called on the module to reduce data size and increase performance.

class ClusterEvent : public EventBase
{
 public:
  ClusterEvent() = default;
  ~ClusterEvent() = default;
  ClusterEvent(const ClusterEvent &event);

  inline void addCluster (Cluster &cluster)  { mClusterVec.emplace_back(cluster); }
  inline void addCluster (std::vector<Cluster> &cluster)  { mClusterVec = cluster; }
  inline void setEventNrCluster (Int_t number) {eventNrCluster = number;}

  std::vector<Cluster>&       getClusters()       { return mClusterVec; }
  const std::vector<Cluster>& getClusters() const { return mClusterVec; }
  Int_t                       getEventNrCluster() const { return eventNrCluster; }

  inline void reset() { mClusterVec.clear(); }

 private:
  std::vector<Cluster> mClusterVec {};  ///< contains all clusters in the event
  Int_t                eventNrCluster;

  ClassDef(ClusterEvent,1);
};

#endif