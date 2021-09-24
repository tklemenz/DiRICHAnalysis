#include "Tracker.h"
#include <experimental/random>

ClassImp(Tracker);

const Float_t  allowedAngle = 60;  //[deg]
const Double_t qMaxThresh   = 8;   //[ToT]
const Double_t timeThresh   = 5;   //[ns]

//________________________________________________________________________________
void Tracker::run(std::vector<Cluster>& clusters, const ParticleType& particletype)
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
  /*bool clustersLeft = true;
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
  }*/

  ///-----------------------------------------

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

  const Double_t toDeg = 180.0/std::acos(-1);

  Cluster* seedCl  = nullptr;
  Cluster* seedCl2 = nullptr;
  std::vector<Cluster> inTrack{};

  Int_t seedFiber   = -1;
  Int_t clFiber     = -1;
  Double_t seedTime = -999;
  Double_t clTime   = -990;
  Int_t seedLayer   = -1;
  Int_t clLayer     = -1;
  Float_t seedCoord = -1;
  Float_t clCoord   = -1;

  Int_t seedFiber2   = -1;
  Int_t clFiber2     = -1;
  Double_t seedTime2 = -999;
  Double_t clTime2   = -990;
  Int_t seedLayer2   = -1;
  Int_t clLayer2     = -1;
  Float_t seedCoord2 = -1;
  Float_t clCoord2   = -1;

  Float_t alpha = -999;
  Float_t beta  = -999;
  Float_t angle = -999;

  Int_t oddCounter  = 0;
  Int_t evenCounter = 0;

  Float_t clusterDistXY = -999;
  Int_t clusterDistZ    = -999;
  Double_t timeDist     = -999;

  bool foundOddSeed  = false; // for odd layers
  bool foundEvenSeed = false; // for even layers

  if (particletype == ParticleType::Cosmic) {
    if (clusters.size() <= 4) { return; }
    //clusters.at(0).setIsUsed();
    if (clusters.at(0).getQMax() >= qMaxThresh) {                                              // take the first cluster in the event as seed if qMax exceeds the threshold
      clusters.at(0).setIsUsed();
      seedCl = &clusters.at(0);
      seedFiber = seedCl->getMeanFiber();
      seedTime  = seedCl->getFirstTimeStamp();
      seedLayer = seedCl->getLayer();
      seedCoord = mapping::getCoord(seedFiber, seedLayer);
      isEven(seedLayer) ? foundEvenSeed = true : foundOddSeed = true;
      if (foundEvenSeed) { printf("%s[TRACKER]%s %s(First seed)%s Set first seed, %sEVEN%s\n",text::LBLU,text::RESET, text::LCYN, text::RESET, text::LGRN, text::RESET); }
      else if (foundOddSeed) { printf("%s[TRACKER]%s %s(First seed)%s Set first seed, %sODD%s\n",text::LBLU,text::RESET, text::LCYN, text::RESET, text::LGRN, text::RESET); }
      else { printf("%s[TRACKER]%s %s(First seed)%s %sNO FIRST SEED FOUND%s\n",text::LBLU,text::RESET, text::LCYN, text::RESET, text::LGRN, text::RESET); }
    }

    for ( auto& cluster : clusters ) {                                                       // loop over the rest of the clusters
      if (cluster.isUsed()) { continue; }
      //if (cluster.getQMax() < qMaxThresh) { continue; }
      clFiber = cluster.getMeanFiber();
      clTime  = cluster.getFirstTimeStamp();
      clLayer = cluster.getLayer();
      clCoord = mapping::getCoord(clFiber, clLayer);

      clusterDistXY = std::abs(clCoord - seedCoord);                                          // distance of clusters in module plane
      clusterDistZ  = 2*(clLayer-seedLayer);                                                  // distance of clusters in z-direction
      timeDist      = clTime - seedTime;                                                      // time distance between the clusters

      angle = std::atan(clusterDistXY/clusterDistZ)*toDeg;

      if (foundOddSeed && !isEven(clLayer) && clLayer < 9) {                                  // fibers are in same direction as seed cluster, seed cluster is in layer 1
        if (angle <= allowedAngle && timeDist <= timeThresh) {
          angle = angle;
          oddCounter++;
          cluster.setIsUsed();
          printf("%s[TRACKER]%s Found odd layer cluster\n",text::LBLU,text::RESET);
        }
      }
      else if (foundEvenSeed && isEven(clLayer) && clLayer < 9) {
        if (angle <= allowedAngle && timeDist <= timeThresh) {
          beta = angle;
          evenCounter++;
          cluster.setIsUsed();
          printf("%s[TRACKER]%s Found even layer cluster\n",text::LBLU,text::RESET);
        }
      }
      else if (clLayer > 8) {
        printf("%s[TRACKER]%s %s(First seed)%s The cluster is in the second module!\n",text::LBLU,text::RESET, text::LCYN, text::RESET);
      }
      else {
        printf("%s[TRACKER]%s %s(First seed)%s No fitting cluster found\n",text::LBLU,text::RESET, text::LCYN, text::RESET);
      }
    }

    printf("%s[TRACKER]%s Finished tracking with first seed!\n",text::LBLU,text::RESET);

    for ( auto& cluster : clusters ) {                                                       // loop over the rest of the clusters again to find the seed for the other direction
      if (cluster.isUsed()) { continue; }
      if (cluster.getQMax() < qMaxThresh) { continue; }
      if (!foundOddSeed && !isEven(cluster.getLayer())) {
        cluster.setIsUsed();
        seedCl2 = &cluster;
        seedFiber2 = seedCl->getMeanFiber();
        seedTime2  = seedCl->getFirstTimeStamp();
        seedLayer2 = seedCl->getLayer();
        seedCoord2 = mapping::getCoord(seedFiber2, seedLayer2);
        foundOddSeed = true;
      }
      if (!foundEvenSeed && isEven(cluster.getLayer())) {
        cluster.setIsUsed();
        seedCl2    = &cluster;
        seedFiber2 = seedCl->getMeanFiber();
        seedTime2  = seedCl->getFirstTimeStamp();
        seedLayer2 = seedCl->getLayer();
        seedCoord2 = mapping::getCoord(seedFiber2, seedLayer2);
        foundEvenSeed = true;
      }
    }

    for ( auto& cluster : clusters ) {                                                       // loop over the rest of the clusters
      if (cluster.isUsed()) { continue; }
      //if (cluster.getQMax() < qMaxThresh) { continue; }
      clFiber = cluster.getMeanFiber();
      clTime  = cluster.getFirstTimeStamp();
      clLayer = cluster.getLayer();
      clCoord = mapping::getCoord(clFiber, clLayer);

      clusterDistXY = std::abs(clCoord - seedCoord2);                                        // distance of clusters in module plane
      clusterDistZ  = 2*(clLayer-seedLayer2);                                                // distance of clusters in z-direction
      timeDist      = clTime - seedTime2;                                                    // time distance between the clusters

      angle = std::atan(clusterDistXY/clusterDistZ)*toDeg;

      if (foundOddSeed && !isEven(clLayer) && clLayer < 9) {                                 // fibers are in same direction as seed cluster, seed cluster is in layer 1
        if (alpha <= allowedAngle && timeDist <= timeThresh) {
          alpha = angle;
          oddCounter++;
          cluster.setIsUsed();
        }
      }
      else if (foundEvenSeed && isEven(clLayer) && clLayer < 9) {
        if (alpha <= allowedAngle && timeDist <= timeThresh) {
          beta = angle;
          evenCounter++;
          cluster.setIsUsed();
        }
      }
      else if (clLayer > 8) {
        printf("%s[TRACKER]%s %s(Second seed)%s The cluster is in the second module!\n",text::LBLU,text::RESET, text::LCYN, text::RESET);
      }
      else {
        printf("%s[TRACKER]%s %s(Second seed)%s Something weird happened!\n",text::LBLU,text::RESET, text::LCYN, text::RESET);
      }
    }

    printf("%s[TRACKER]%s Finished tracking with second seed!\n",text::LBLU,text::RESET);
  }

  for ( auto& cluster : clusters ) {
    if (cluster.isUsed()) { inTrack.emplace_back(cluster); }                                 // Now write the clusters marked as used into the vector
  }

  Float_t xVertex = -1;
  Float_t yVertex = -1;

  isEven(seedLayer)  ? yVertex = seedFiber  : xVertex = seedFiber;
  isEven(seedLayer2) ? yVertex = seedFiber2 : xVertex = seedFiber2;

  if (xVertex == -1 || yVertex == -1) {
    printf("%s[TRACKER]%s Both seeds are in the same direction!!\n",text::LBLU,text::RESET);
    return;
  }

  /// Check for alpha and beta before creating tracks

  else {
    printf("%s[TRACKER]%s This particle type is not implemented for the tracker!\n",text::LBLU,text::RESET);
  }
}

//________________________________________________________________________________
void Tracker::run(CTSEventClusters& event, const ParticleType& particletype)
{
  Tracker::run(event.getClusters(), particletype);
}
