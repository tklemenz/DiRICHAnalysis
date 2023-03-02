#include "Tracker.h"

#include <cmath>

ClassImp(Tracker);

static const std::string process = "[TRACKER]";

const Float_t  allowedAngle = 60;  //[deg]
const Double_t qMaxThresh   = 8;   //[ToT]
const Double_t timeThresh   = 5;   //[ns]
const Int_t skipEventThr    = 5;   // event has to have at least this amount of clusters to try tracking

//________________________________________________________________________________
void Tracker::run(std::vector<Cluster>& clusters, const ParticleType& particletype, const bool& debug)
{
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
  Cluster* tmpSeed = nullptr;
  std::vector<Cluster> inTrack{};

  Float_t seedFiber   = -1;
  Float_t clFiber     = -1;
  Double_t seedTime = -999;
  Double_t clTime   = -990;
  Int_t seedLayer   = -1;
  Int_t clLayer     = -1;
  Float_t seedCoord = -1;
  Float_t clCoord   = -1;

  Float_t seedFiber2   = -1;
  Double_t seedTime2 = -999;
  Int_t seedLayer2   = -1;
  Float_t seedCoord2 = -1;

  Float_t tmpSeedFiber   = -1;
  Int_t tmpSeedLayer   = -1;
  Float_t tmpSeedCoord = -1;

  Float_t alpha = -999;
  Float_t beta  = -999;
  Float_t angle = -999;

  Int_t oddCounter  = 0;
  Int_t evenCounter = 0;

  Float_t clusterDistXY = -999;
  Float_t clusterDistZ  = -999;
  Double_t timeDist     = -999;
  Double_t seedTimeDist = -999;

  Int_t nOddSeedCandidatesL1  = 0;
  Int_t nOddSeedCandidatesL3  = 0;
  Int_t nEvenSeedCandidatesL2 = 0;
  Int_t nEvenSeedCandidatesL4 = 0;

  bool foundOddSeed  = false; // for odd layers
  bool foundEvenSeed = false; // for even layers

  if (particletype == ParticleType::Cosmic) {
    if (debug) { printf("%s%s %sStart of Event %lu.%s\n",text::LBLU, process.c_str(), text::LCYN, mEventCounter, text::RESET); }
    if (clusters.size() < skipEventThr) {
      if (debug) { 
        printf("%s%s %sLess than %i Clusters in Event.%s\n",text::LBLU, process.c_str(), text::LCYN, skipEventThr, text::RESET);
        printf("%s%s %sSkipped Event %lu.%s\n\n\n",text::LBLU, process.c_str(), text::LCYN, mEventCounter, text::RESET);
      }
      mEventCounter++;
      return;
    }

  //=============================================================================================================================================================================== START FIRST SEED
    for (auto& cluster : clusters) {
      if (cluster.getQMax() < qMaxThresh) { continue; }
      switch(cluster.getLayer()) {
        case 1:
          nOddSeedCandidatesL1++;
          break;
        case 2:
          nEvenSeedCandidatesL2++;
          break;
        case 3:
          nOddSeedCandidatesL3++;
          break;
        case 4:
          nEvenSeedCandidatesL4++;
          break;
        default:
          break;
      }
    }
    if (debug) {
      printf("%s%s %sNumber of odd  seed candidates: %i (Layer 1).%s\n",text::LBLU, process.c_str(), text::GRN, nOddSeedCandidatesL1, text::RESET);
      printf("%s%s %sNumber of even seed candidates: %i (Layer 2).%s\n",text::LBLU, process.c_str(), text::GRN, nEvenSeedCandidatesL2, text::RESET);
      printf("%s%s %sNumber of odd  seed candidates: %i (Layer 3).%s\n",text::LBLU, process.c_str(), text::GRN, nOddSeedCandidatesL3, text::RESET);
      printf("%s%s %sNumber of even seed candidates: %i (Layer 4).%s\n",text::LBLU, process.c_str(), text::GRN, nEvenSeedCandidatesL4, text::RESET);
    }

    if (clusters.front().getQMax() >= qMaxThresh) {                                              // take the first cluster in the event as seed if qMax exceeds the threshold
      clusters.front().setIsCandidate();
      seedCl = &clusters.front();
      seedFiber = seedCl->getMeanFiber();
      seedTime  = seedCl->getFirstTimeStamp();
      seedLayer = seedCl->getLayer();
      seedCoord = mapping::getCoord(seedFiber, seedLayer);
      isEven(seedLayer) ? foundEvenSeed = true : foundOddSeed = true;
      clusters.front().setIsSeed();
      clusters.front().setIsUsed();

      if(debug) {
        if (foundEvenSeed && !foundOddSeed) { printf("%s%s%s %s(First seed)%s Set first seed, %sEVEN%s in layer %d\n",text::LBLU, process.c_str(), text::RESET, text::LCYN, text::RESET, text::LGRN, text::RESET, seedLayer); }
        else if (foundOddSeed && !foundEvenSeed) { printf("%s%s%s %s(First seed)%s Set first seed, %sODD%s in layer %d\n",text::LBLU, process.c_str(), text::RESET, text::LCYN, text::RESET, text::LGRN, text::RESET, seedLayer); }
        else { printf("%s%s%s %s(First seed)%s %sNO FIRST SEED FOUND%s\n",text::LBLU, process.c_str(), text::RESET, text::LCYN, text::RESET, text::LGRN, text::RESET); }
        printf("%s%s %s(First seed)%s %sSeed info:%s\n",text::LBLU, process.c_str(), text::LCYN, text::RESET, text::LBLU, text::RESET);
        printf("%s%s %s(First seed)%s\tmean fiber  : %g\n",text::LBLU, process.c_str(), text::LCYN, text::RESET, seedFiber);
        printf("%s%s %s(First seed)%s\ttime stamp  : %g ns\n",text::LBLU, process.c_str(), text::LCYN, text::RESET, seedTime);
        printf("%s%s %s(First seed)%s\tlayer       : %d\n",text::LBLU, process.c_str(), text::LCYN, text::RESET, seedLayer);
        printf("%s%s %s(First seed)%s\tcoordinate  : %g mm\n",text::LBLU, process.c_str(), text::LCYN, text::RESET, seedCoord);
      }
    }

    for ( auto& cluster : clusters ) {                                                       // loop over the rest of the clusters
      if (cluster.isCandidate()) { continue; }
      //if (cluster.getQMax() < qMaxThresh) { continue; }
      clFiber = cluster.getMeanFiber();
      clTime  = cluster.getFirstTimeStamp();
      clLayer = cluster.getLayer();
      clCoord = mapping::getCoord(clFiber, clLayer);

      clusterDistXY = std::abs(clCoord - seedCoord);                                          // distance of clusters in module plane
      clusterDistZ  = 2*(clLayer-seedLayer);                                                  // distance of clusters in z-direction
      timeDist      = std::abs(clTime - seedTime);                                            // time distance between the clusters

      angle = std::atan(clusterDistXY/clusterDistZ)*toDeg;

      if(debug) {
        printf("%s%s %s(First seed)%s %sProcessing Cluster...%s\n",text::LBLU, process.c_str(), text::LCYN, text::RESET, text::LBLU, text::RESET);
        printf("%s%s %s(First seed)%s\tmean fiber  : %g\n",text::LBLU, process.c_str(), text::LCYN, text::RESET, clFiber);
        printf("%s%s %s(First seed)%s\ttime stamp  : %g ns\n",text::LBLU, process.c_str(), text::LCYN, text::RESET, clTime);
        printf("%s%s %s(First seed)%s\tlayer       : %d\n",text::LBLU, process.c_str(), text::LCYN, text::RESET, clLayer);
        printf("%s%s %s(First seed)%s\tcoordinate  : %g mm\n",text::LBLU, process.c_str(), text::LCYN, text::RESET, clCoord);
        printf("%s%s %s(First seed)%s\tseed dist xy: %g mm\n",text::LBLU, process.c_str(), text::LCYN, text::RESET, clusterDistXY);
        printf("%s%s %s(First seed)%s\tseed dist z : %g mm\n",text::LBLU, process.c_str(), text::LCYN, text::RESET, clusterDistZ);
        printf("%s%s %s(First seed)%s\tseed dist t : %g ns\n",text::LBLU, process.c_str(), text::LCYN, text::RESET, timeDist);
        printf("%s%s %s(First seed)%s\tangle       : %g deg\n",text::LBLU, process.c_str(), text::LCYN, text::RESET, angle);
      }

      if (foundOddSeed && !isEven(clLayer) && clLayer < 9) {                                  // fibers are in same direction as seed cluster, seed cluster is in layer 1
        if (angle <= allowedAngle && timeDist <= timeThresh) {
          angle = angle;
          oddCounter++;
          cluster.setIsCandidate();
          if(debug) { printf("%s%s%s %s(First seed)%s Found odd layer cluster\n",text::LBLU,process.c_str(), text::RESET, text::LCYN, text::RESET); }
        }
      }
      else if (foundEvenSeed && isEven(clLayer) && clLayer < 9) {
        if (angle <= allowedAngle && timeDist <= timeThresh) {
          beta = angle;
          evenCounter++;
          cluster.setIsCandidate();
          if(debug) { printf("%s%s%s %s(First seed)%s Found even layer cluster\n",text::LBLU, process.c_str(), text::RESET, text::LCYN, text::RESET); }
        }
      }
      else if (clLayer > 8) {
        if(debug) { printf("%s%s%s %s(First seed)%s The cluster is in the second module!\n",text::LBLU, process.c_str(), text::RESET, text::LCYN, text::RESET); }
      }
      else {
        if(debug) { printf("%s%s%s %s(First seed)%s No fitting cluster found\n",text::LBLU, process.c_str(), text::RESET, text::LCYN, text::RESET); }
      }
    }

    if(debug) { printf("%s%s%s %s%sFinished tracking with first seed!%s\n",text::LBLU, process.c_str(), text::RESET, text::BOLD, text::GRN, text::RESET); }

  //=============================================================================================================================================================================== FIRST SEED DONE

    if(debug) {
      for ( auto& cluster : clusters ) {
        if (cluster.isCandidate()) {
          printf("%s%s %s[INFO] Candidate in layer %d.%s\n",text::LBLU, process.c_str(), text::GRN, cluster.getLayer(), text::RESET);
        }
      }
    }

  //=============================================================================================================================================================================== START SECOND SEED

    for ( auto& cluster : clusters ) {                                                       // loop over the rest of the clusters again to find the seed for the other direction
      if (cluster.isCandidate()) {
        if(debug) { printf("%s%s %s(Second seed) %sSkipping used Cluster in layer %d.%s\n",text::LBLU, process.c_str(), text::LCYN, text::LGRN, cluster.getLayer(), text::RESET); }
        continue;
      }
      if (cluster.getQMax() < qMaxThresh) { continue; }
      if (std::abs(seedTime - cluster.getFirstTimeStamp()) > timeThresh) { continue; }
      if (!foundOddSeed && !isEven(cluster.getLayer())) {
        cluster.setIsCandidate();
        seedCl2 = &cluster;
        seedFiber2 = seedCl2->getMeanFiber();
        seedTime2  = seedCl2->getFirstTimeStamp();
        seedLayer2 = seedCl2->getLayer();
        seedCoord2 = mapping::getCoord(seedFiber2, seedLayer2);
        seedTimeDist = std::abs(seedTime - seedTime2);
        foundOddSeed = true;
        cluster.setIsSeed();
        cluster.setIsUsed();

        if(debug) {
          printf("%s%s%s %s(Second seed)%s Set first seed, %sODD%s in layer %d\n",text::LBLU, process.c_str(), text::RESET, text::LCYN, text::RESET, text::LGRN, text::RESET, seedLayer2);
          printf("%s%s %s(Second seed)%s %sSeed info:%s\n",text::LBLU, process.c_str(), text::LCYN, text::RESET, text::LBLU, text::RESET);
          printf("%s%s %s(Second seed)%s\tseed dist t : %g ns\n",text::LBLU, process.c_str(), text::LCYN, text::RESET, seedTimeDist);
          printf("%s%s %s(Second seed)%s\tmean fiber  : %g\n",text::LBLU, process.c_str(), text::LCYN, text::RESET, seedFiber2);
          printf("%s%s %s(Second seed)%s\ttime stamp  : %g ns\n",text::LBLU, process.c_str(), text::LCYN, text::RESET, seedTime2);
          printf("%s%s %s(Second seed)%s\tlayer       : %d\n",text::LBLU, process.c_str(), text::LCYN, text::RESET, seedLayer2);
          printf("%s%s %s(Second seed)%s\tcoordinate  : %g mm\n",text::LBLU, process.c_str(), text::LCYN, text::RESET, seedCoord2);
        }
      }
      if (!foundEvenSeed && isEven(cluster.getLayer())) {
        cluster.setIsCandidate();
        seedCl2    = &cluster;
        seedFiber2 = seedCl2->getMeanFiber();
        seedTime2  = seedCl2->getFirstTimeStamp();
        seedLayer2 = seedCl2->getLayer();
        seedCoord2 = mapping::getCoord(seedFiber2, seedLayer2);
        seedTimeDist = std::abs(seedTime - seedTime2);
        foundEvenSeed = true;
        cluster.setIsSeed();
        cluster.setIsUsed();

        if(debug) {
          printf("%s%s%s %s(Second seed)%s Set second seed, %sEVEN%s in layer %d\n",text::LBLU, process.c_str(), text::RESET, text::LCYN, text::RESET, text::LGRN, text::RESET, seedLayer2);
          printf("%s%s %s(Second seed)%s %sSeed info:%s\n",text::LBLU, process.c_str(), text::LCYN, text::RESET, text::LBLU, text::RESET);
          printf("%s%s %s(Second seed)%s\tseed dist t : %g ns\n",text::LBLU, process.c_str(), text::LCYN, text::RESET, seedTimeDist);
          printf("%s%s %s(Second seed)%s\tmean fiber  : %g\n",text::LBLU, process.c_str(), text::LCYN, text::RESET, seedFiber2);
          printf("%s%s %s(Second seed)%s\ttime stamp  : %g ns\n",text::LBLU, process.c_str(), text::LCYN, text::RESET, seedTime2);
          printf("%s%s %s(Second seed)%s\tlayer       : %d\n",text::LBLU, process.c_str(), text::LCYN, text::RESET, seedLayer2);
          printf("%s%s %s(Second seed)%s\tcoordinate  : %g mm\n",text::LBLU, process.c_str(), text::LCYN, text::RESET, seedCoord2);
        }
      }
    }

    for ( auto& cluster : clusters ) {                                                       // loop over the rest of the clusters
      if (cluster.isCandidate()) { continue; }
      //if (cluster.getQMax() < qMaxThresh) { continue; }
      clFiber = cluster.getMeanFiber();
      clTime  = cluster.getFirstTimeStamp();
      clLayer = cluster.getLayer();
      clCoord = mapping::getCoord(clFiber, clLayer);

      clusterDistXY = std::abs(clCoord - seedCoord2);                                        // distance of clusters in module plane
      clusterDistZ  = 2*(clLayer-seedLayer2);                                                // distance of clusters in z-direction
      timeDist      = std::abs(clTime - seedTime2);                                          // time distance between the clusters

      angle = std::atan(clusterDistXY/clusterDistZ)*toDeg;

      if(debug) {
        printf("%s%s %s(Second seed)%s %sProcessing Cluster...%s\n",text::LBLU, process.c_str(), text::LCYN, text::RESET, text::LBLU, text::RESET);
        printf("%s%s %s(Second seed)%s\tmean fiber  : %g\n",text::LBLU, process.c_str(), text::LCYN, text::RESET, clFiber);
        printf("%s%s %s(Second seed)%s\ttime stamp  : %g ns\n",text::LBLU, process.c_str(), text::LCYN, text::RESET, clTime);
        printf("%s%s %s(Second seed)%s\tlayer       : %d\n",text::LBLU, process.c_str(), text::LCYN, text::RESET, clLayer);
        printf("%s%s %s(Second seed)%s\tcoordinate  : %g mm\n",text::LBLU, process.c_str(), text::LCYN, text::RESET, clCoord);
        printf("%s%s %s(Second seed)%s\tseed dist xy: %g mm\n",text::LBLU, process.c_str(), text::LCYN, text::RESET, clusterDistXY);
        printf("%s%s %s(Second seed)%s\tseed dist z : %g mm\n",text::LBLU, process.c_str(), text::LCYN, text::RESET, clusterDistZ);
        printf("%s%s %s(Second seed)%s\tseed dist t : %g ns\n",text::LBLU, process.c_str(), text::LCYN, text::RESET, timeDist);
        printf("%s%s %s(Second seed)%s\tangle       : %g deg\n",text::LBLU, process.c_str(), text::LCYN, text::RESET, angle);
      }

      if (foundOddSeed && !isEven(clLayer) && clLayer < 9) {                                 // fibers are in same direction as seed cluster, seed cluster is in layer 1
        if (alpha <= allowedAngle && timeDist <= timeThresh) {
          alpha = angle;
          oddCounter++;
          cluster.setIsCandidate();
          if(debug) { printf("%s%s%s %s(Second seed)%s Found odd layer cluster\n",text::LBLU, process.c_str(), text::RESET, text::LCYN, text::RESET); }
        }
      }
      else if (foundEvenSeed && isEven(clLayer) && clLayer < 9) {
        if (alpha <= allowedAngle && timeDist <= timeThresh) {
          beta = angle;
          evenCounter++;
          cluster.setIsCandidate();
          if(debug) { printf("%s%s%s %s(Second seed)%s Found even layer cluster\n",text::LBLU, process.c_str(), text::RESET, text::LCYN, text::RESET); }
        }
      }
      else if (clLayer > 8) {
        if(debug) { printf("%s%s%s %s(Second seed)%s The cluster is in the second module!\n",text::LBLU, process.c_str(), text::RESET, text::LCYN, text::RESET); }
      }
      else {
        if(debug) { printf("%s%s%s %s(Second seed)%s No fitting cluster found\n",text::LBLU, process.c_str(), text::RESET, text::LCYN, text::RESET); }
      }
    }

    if(debug) { printf("%s%s%s %s%sFinished tracking with second seed!%s\n",text::LBLU, process.c_str(), text::RESET, text::BOLD, text::GRN, text::RESET); }

  //=============================================================================================================================================================================== SECOND SEED DONE

  //=============================================================================================================================================================================== CHECK IF CLUSTERS ACTUALLY FIT

    std::vector<Cluster> candidates{};

    for ( auto& cluster : clusters ) {
      if (cluster.isCandidate()) { candidates.emplace_back(cluster); }                                 // Now write the clusters marked as used into the candidates vector
    }

    if (oddCounter == 0) {
      if (debug) {
        printf("%s%s %s[ERROR] No tracklet in odd Layers found!%s\n",text::LBLU, process.c_str(), text::LRED, text::RESET);
        printf("%s%s %sEnd of Event %lu.%s\n\n\n",text::LBLU, process.c_str(), text::LCYN, mEventCounter, text::RESET);
      }
      mEventCounter++;
      return;
    }
    if (evenCounter == 0) {
      if (debug) {
        printf("%s%s %s[ERROR] No tracklet in even Layers found!%s\n",text::LBLU, process.c_str(), text::LRED, text::RESET);
        printf("%s%s %sEnd of Event %lu.%s\n\n\n",text::LBLU, process.c_str(), text::LCYN, mEventCounter, text::RESET);
      }
      mEventCounter++;
      return;
    }

    Float_t xVertex = -1;  // odd layer seed
    Float_t yVertex = -1;  // even layer seed

    isEven(seedLayer)  ? yVertex = seedFiber  : xVertex = seedFiber;
    isEven(seedLayer2) ? yVertex = seedFiber2 : xVertex = seedFiber2;

    if (xVertex == -1 || yVertex == -1) {
      if(debug) {
        printf("%s%s %s[ERROR] Both seeds are in the same direction!!%s\n",text::LBLU, process.c_str(), text::LRED, text::RESET);
        printf("%s%s %sEnd of Event %lu.%s\n\n\n",text::LBLU, process.c_str(), text::LCYN, mEventCounter, text::RESET);
        mEventCounter++;
      }
      return;
    }

    // Check for alpha and beta before creating tracks
    std::vector<Float_t> angles{};
    Float_t meanAngle = 0;
    Float_t sumAngles = 0;
    Float_t varAnglesAlpha = 0;
    Float_t tmpSum    = 0;

    //-----------------------------------------------------------------------------------------------------------------------------------
    !isEven(seedLayer) ? tmpSeed = seedCl : tmpSeed = seedCl2;                                // look at odd layers first --> alpha
    tmpSeedFiber = tmpSeed->getMeanFiber();
    tmpSeedLayer = tmpSeed->getLayer();
    tmpSeedCoord = mapping::getCoord(tmpSeedFiber, tmpSeedLayer);

    for (auto& cluster : candidates) {
      if (isEven(cluster.getLayer()) || cluster.isSeed()) { continue; }
      if (debug) { printf("%s%s %s[INFO] Cluster in layer %d. Should be odd!%s\n",text::LBLU, process.c_str(), text::GRN, cluster.getLayer(), text::RESET); }
      clFiber = cluster.getMeanFiber();
      clLayer = cluster.getLayer();
      clCoord = mapping::getCoord(clFiber, clLayer);

      clusterDistXY = std::abs(clCoord - tmpSeedCoord);                                          // distance of clusters in module plane
      clusterDistZ  = 2*(clLayer-tmpSeedLayer);                                                  // distance of clusters in z-direction

      if (debug) {
        printf("%s%s%s cl dist xy : %g mm\n",text::LBLU, process.c_str(), text::RESET, clusterDistXY);
        printf("%s%s%s cl dist z  : %g mm\n",text::LBLU, process.c_str(), text::RESET, clusterDistZ);
      }

      angle = std::atan(clusterDistXY/clusterDistZ)*toDeg;
      angles.emplace_back(angle);
      sumAngles += angle;
    }

    meanAngle = sumAngles/angles.size();
    for (auto& angle : angles) {
      tmpSum += ((angle-meanAngle)*(angle-meanAngle));
    }
    varAnglesAlpha = tmpSum/(angles.size()-1);
    if(debug) {
      printf("%s%s%s VAR alpha: %g\n",text::LBLU, process.c_str(), text::RESET, varAnglesAlpha);
      printf("%s%s%s Sum angles   : %g\n",text::LBLU, process.c_str(), text::RESET, sumAngles);
      printf("%s%s%s nAngles      : %lu\n",text::LBLU, process.c_str(), text::RESET, angles.size());
      for (auto& angle : angles) {
        printf("%s%s%s angle        : %g deg\n",text::LBLU, process.c_str(), text::RESET, angle);
      }
      printf("%s%s%s mean diff sum: %g\n",text::LBLU, process.c_str(), text::RESET, tmpSum);
    }

    //-----------------------------------------------------------------------------------------------------------------------------------
    angles.clear();
    meanAngle = 0;
    sumAngles = 0;
     Float_t varAnglesBeta = 0;
    tmpSum    = 0;
    tmpSeed = nullptr;

    isEven(seedLayer) ? tmpSeed = seedCl : tmpSeed = seedCl2;                                // look at even layers --> beta
    tmpSeedFiber = tmpSeed->getMeanFiber();
    tmpSeedLayer = tmpSeed->getLayer();
    tmpSeedCoord = mapping::getCoord(tmpSeedFiber, tmpSeedLayer);

    for (auto& cluster : candidates) {
      if (!isEven(cluster.getLayer()) || cluster.isSeed()) { continue; }
      if (debug) { printf("%s%s %s[INFO] Cluster in layer %d. Should be even!%s\n",text::LBLU, process.c_str(), text::GRN, cluster.getLayer(), text::RESET); }
      clFiber = cluster.getMeanFiber();
      clLayer = cluster.getLayer();
      clCoord = mapping::getCoord(clFiber, clLayer);

      clusterDistXY = std::abs(clCoord - tmpSeedCoord);                                          // distance of clusters in module plane
      clusterDistZ  = 2*(clLayer-tmpSeedLayer);                                                  // distance of clusters in z-direction

      if (debug) {
        printf("%s%s%s cl dist xy : %g mm\n",text::LBLU, process.c_str(), text::RESET, clusterDistXY);
        printf("%s%s%s cl dist z  : %g mm\n",text::LBLU, process.c_str(), text::RESET, clusterDistZ);
      }

      angle = std::atan(clusterDistXY/clusterDistZ)*toDeg;
      angles.emplace_back(angle);
      sumAngles += angle;
    }

    meanAngle = sumAngles/angles.size();
    for (auto& angle : angles) {
      tmpSum += ((angle-meanAngle)*(angle-meanAngle));
    }
    varAnglesBeta = tmpSum/(angles.size()-1);
    if(debug) {
      printf("%s%s%s VAR beta     : %g\n",text::LBLU, process.c_str(), text::RESET, varAnglesBeta);
      printf("%s%s%s Sum angles   : %g\n",text::LBLU, process.c_str(), text::RESET, sumAngles);
      printf("%s%s%s nAngles      : %lu\n",text::LBLU, process.c_str(), text::RESET, angles.size());
      for (auto& angle : angles) {
        printf("%s%s%s angle        : %g deg\n",text::LBLU, process.c_str(), text::RESET, angle);
      }
      printf("%s%s%s mean diff sum: %g\n",text::LBLU, process.c_str(), text::RESET, tmpSum);
    }

    //-----------------------------------------------------------------------------------------------------------------------------------
    //=============================================================================================================================================================================== CHECK DONE

    for ( auto& cluster : candidates ) {
      if (cluster.isUsed()) { inTrack.emplace_back(cluster); }                                 // Now write the clusters marked as used into the cluster vector which becomes the track
    }
  } // loop over cosmics

  else {
    printf("%s%s %sThis particle type is not implemented for the tracker!%s\n",text::LBLU, process.c_str(), text::LRED, text::RESET);
  }

  if (debug) { printf("%s%s %sEnd of Event %lu.%s\n\n\n",text::LBLU, process.c_str(), text::LCYN, mEventCounter, text::RESET); }
  mEventCounter++;

  if (inTrack.size() > 0) {                                                                   // Put the tracks in the member variable and clear the tmp container
    mTrackVec.emplace_back(Track(std::move(inTrack), 0, 0, std::pair<Float_t, Float_t>(0,0), ParticleType(ParticleType::Cosmic)));
    inTrack.clear();
  }
}

//________________________________________________________________________________
void Tracker::run(CTSEventClusters& event, const ParticleType& particletype, const bool& debug)
{
  Tracker::run(event.getClusters(), particletype, debug);
}
