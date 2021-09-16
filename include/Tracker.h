#ifndef TRACKER_H
#define TRACKER_H

#include "CTSEventClusters.h"
#include "Track.h"
#include "Utility.h"

/// Tracker class to assign clusters from CTSEventClusters to tracks.
/// Ideally a track consists of 4 clusters, one in each layer.
///
/// Option 1: The Tracker takes a CTSEventClusters and find tracks from the clusters.
/// Option 2: The Tracker gets a vector of Clusters and does the tracking.
///           Then the user needs to make sure that the clusters are all valid (e.g. not 5 seconds apart).
///
/// The idea at the moment is to give a set of clusters to the tracker who then puts the proper
/// clusters together and creates a temporary instance of Track and copies it to the Track vector
/// member of the Tracker class. 
///
/// ParticleType of the track is not so important for us but if we want to have it
/// we need to think about how to set the proper type. In general it is possible to
/// have pions in the proton data! Otherwise it would be easy since we know which
/// particles we expect in every run. For now I would set unknown by default.

class Tracker
{
 public:
  Tracker() = default;
  ~Tracker() = default;

  /// run the tracking
  void run(CTSEventClusters& event);
  void run(std::vector<Cluster>& clusters);

  std::vector<Track>&       getTracks()       { return mTrackVec; }
  const std::vector<Track>& getTracks() const { return mTrackVec; }

 private:

  std::vector<Track> mTrackVec{}; ///< contains all found tracks

  ClassDef(Tracker,1);
};

#endif