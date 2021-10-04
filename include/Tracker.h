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

/// For cosmic particles:
/// vertical fibers
///           /|
///          /a|
///         /  |  4mm (2 fibers, distance between "neighboring" layers; 1->3 4mm, 1->5 8mm, 1->7 12mm)
///        /   |
///       /____|
///         dx

/// dx: distance of coordinates of clusters in "neighboring" same direction fibers e.g. layer 1 and layer 3, dx is negative in this example
/// a = atan(dx/4), a is negative in this example

/// horizontal fibers
///           /|
///          /b|
///         /  |  4mm (2 fibers, distance between "neighboring" layers; 2->4 4mm, 2->6 8mm, 2->8 12mm)
///        /   |
///       /____|
///         dy

/// dy: distance of coordinates of clusters in "neighboring" same direction fibers e.g. layer 2 and layer 4, dy is negative in this example
/// b = atan(dy/4), b is negative in this example

class Tracker
{
 public:
  Tracker() = default;
  ~Tracker() = default;

  /// run the tracking
  /// @todo Properly implement tracking in second module
  void run(CTSEventClusters& event, const ParticleType& particletype, const bool& debug = false);
  void run(std::vector<Cluster>& clusters, const ParticleType& particletype, const bool& debug = false);

  std::vector<Track>&       getTracks()       { return mTrackVec; }
  const std::vector<Track>& getTracks() const { return mTrackVec; }

  inline void reset() { mTrackVec.clear(); }

 private:

  std::vector<Track> mTrackVec{}; ///< contains all found tracks
  ULong_t mEventCounter = 0;

  ClassDef(Tracker,1);
};

#endif