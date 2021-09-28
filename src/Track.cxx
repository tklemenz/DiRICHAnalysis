#include "Track.h"

ClassImp(Track);

//________________________________________________________________________________
Track::Track(const Track &track)
: mClusterVec(track.mClusterVec),
  mAlpha(track.mAlpha),
  mBeta(track.mBeta),
  mVertex(track.mVertex),
  mParticleType(track.mParticleType)
{

}

//________________________________________________________________________________
Track::Track(const std::vector<Cluster> &clusterVec, const Float_t &alpha, const Float_t &beta, const std::pair<Float_t, Float_t> &vertex, const ParticleType &type)
: mClusterVec(clusterVec),
  mAlpha(alpha),
  mBeta(beta),
  mVertex(vertex),
  mParticleType(type)
{

}