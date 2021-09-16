#include "Track.h"

ClassImp(Track);

//________________________________________________________________________________
Track::Track(const Track &track)
: mClusterVec(track.mClusterVec),
  mParticleType(track.mParticleType)
{

}

//________________________________________________________________________________
Track::Track(const std::vector<Cluster> &clusterVec, const ParticleType &type)
: mClusterVec(clusterVec),
  mParticleType(type)
{

}