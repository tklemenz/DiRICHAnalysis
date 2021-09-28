#ifndef TRACK_H
#define TRACK_H

#include "Cluster.h"

/// This class represents a track in the module.
/// Ideally a track consists of 4 clusters, one in each layer.
///
/// The idea at the moment is to give a set of clusters to the tracker who then puts the proper
/// clusters together and creates a temporary instance of Track and copies it to the Track vector
/// member of the Tracker class. 
/// --> only the Track constructor which gets the clusterVec and ParticleType is for now
///     intended to be used. Therefore, addCluster and setType are still there if needed
///     but made private to prevent usage outside of the class.
///
/// ParticleType of the track is not so important for us but if we want to have it
/// we need to think about how to set the proper type. In general it is possible to
/// have pions in the proton data! Otherwise it would be easy since we know which
/// particles we expect in every run. For now I would set unknown by default.
/// @todo Think about more info to be added to the track. E.g. angles, position (could be x-y in layer 1+2 and layer 3+4 separately)

class Track
{
 public:

  /// default constructor
  Track() = default;

  /// default destructor
  ~Track() = default;

  /// copy constructor
  Track(const Track &track);

  /// This one would most likely be used. Gets all the track information directly.
  Track(const std::vector<Cluster> &clusterVec, const Float_t &alpha, const Float_t &beta, const std::pair<Float_t, Float_t> &vertex, const ParticleType &type);

  /// Get the clusters in a std::vector
  inline std::vector<Cluster>& getClusters() { return mClusterVec; }

  /// Get the particle type
  inline ParticleType getType() { return mParticleType; }

  /// Get alpha
  inline Float_t getAlpha() { return mAlpha; }

  /// Get beta
  inline Float_t getBeta() { return mBeta; }

  /// Get vertex
  inline std::pair<Float_t, Float_t> getVertex() { return mVertex; }

 private:

  std::vector<Cluster>        mClusterVec{}; ///< holds all clusters that are assigned to the track
  ParticleType                mParticleType; ///< particle type: can be Pion, Proton or Unknown
  Float_t                     mAlpha;        ///< angle between track and vertical line while looking in odd layer fiber direction. Positive alpha: particle is traveling in direction of higher fiber numbers
  Float_t                     mBeta;         ///< angle between track and horizontal line while looking in even layer fiber direction. Positive beta: particle is traveling in direction of higher fiber numbers
  std::pair<Float_t, Float_t> mVertex;       ///< mean fiber numbers of the seed clusters; first: x-direction (odd Layer), second: y-direction (even layer)

  /// Associate a Cluster with the track.
  /// @param Cluster
  inline void addCluster(const Cluster &cluster) { mClusterVec.emplace_back(cluster); }

  /// Set the ParticleType
  /// @param ParticleType
  inline void setType(ParticleType &type) { mParticleType = type; }

  /// Set alpha
  inline void setAlpha(Float_t &alpha) { mAlpha = alpha; }

  /// Set beta
  inline void setBeta(Float_t &beta) { mAlpha = beta; }

  /// Set vertex
  inline void setVertex(std::pair<Float_t, Float_t> &vertex) { mVertex = vertex; }

  ClassDef(Track,1);
};

#endif