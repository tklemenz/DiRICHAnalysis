#include "Tracker.h"
#include <experimental/random>

ClassImp(Tracker);

//________________________________________________________________________________
void Tracker::run(std::vector<Cluster>& clusters)
{
  /// - find the proper clusters and put them into a vector
  /// - put Track object in Track vector
  /// - somehow take care of the possibility that there could be more than one track
  ///   NOT considered in this example!!
  ///
  /// actual position in mm (with (0/0) coordinate at bottom left if looking directly at the module (like in 0 deg runs))
  /// can be obtained by the mapping::getCoord(Float_t meanFiber), where meanFiber is taken from the Cluster object.
  ///
  /// something like
  /// 
  /// std::vector<Cluster> inTrack{}; //(maybe the brackets are not needed)
  /// for (auto &cluster : clusters) {
  ///   if (condition) {
  ///     inTrack.emplace_back(cluster);
  ///   }
  /// }
  /// mTrackVec.emplace_back(Track(inTrack, ParticleType(Unknown)));


  /// This is only a very basic playing around
  /// Is checked with macro tryAroundTracker
  ///-----------------------------------------
  bool clustersLeft = true;
  Int_t nClusters = clusters.size();
  Int_t used = 0;

  std::vector<Cluster> inTrack{};

  while (clustersLeft) {
    for (auto &cluster : clusters) {
      if (!(cluster.isUsed()) && (std::experimental::randint(1,10) <= 5)) {
        inTrack.emplace_back(cluster);
        cluster.setIsUsed();
        used++;
      }
    }

    if (inTrack.size() > 0) {
      mTrackVec.emplace_back(Track(std::move(inTrack), ParticleType(ParticleType::Unknown)));
      inTrack.clear();
    }

    if (used == nClusters) { clustersLeft = false; }
  }

  ///-----------------------------------------
}

//________________________________________________________________________________
void Tracker::run(CTSEventClusters& event)
{
  Tracker::run(event.getClusters());
}
