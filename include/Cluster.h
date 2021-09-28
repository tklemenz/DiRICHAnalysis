#ifndef CLUSTER_H
#define CLUSTER_H

#include "Signal.h"
#include "Utility.h"
#include <vector>

///!
///! This class represents a cluster of signals.
///! Clusters consist of signals in neighboring fibers that stem from the same particle.
///! The idea is to initialize a cluster object and call the function addSignal
///! for every signal that shall be included in the cluster.
///!

class Cluster
{
 public:
  enum Flags : unsigned short {
    usedInTrack         = 0x1 << 0,         ///< Check if the Cluster is already associated to a Track.
    isCandidateForTrack = 0x1 << 1,         ///< Check if the Cluster is a candidate for a Track
    isSeedForDirection  = 0x1 << 2          ///< Cluster is used as seed
  };

  /// default constructor
  Cluster() = default;

  /// default destructor
  ~Cluster() = default;

  /// Copy constructor
  Cluster(const Cluster &cluster);

  /// This constructor should be used in most cases. It takes all Cluster information.
  Cluster(const Double_t qTot, const Double_t qMax, const Float_t meanFiber, const Float_t sigmaFiber, const Double_t meanTimeStamp,
          const Double_t sigmaTimeStamp, const Double_t firstTimeStamp, const Int_t layer, const Int_t TDCID, const Short_t flags);

  /// Check if the Cluster was already associated to a Track.
  inline bool isUsed() const { return mFlags & usedInTrack; }

  /// Check if the Cluster is a candidate for a Track.
  inline bool isCandidate() const { return mFlags & isCandidateForTrack; }

  /// Check if the Cluster is used as seed.
  inline bool isSeed() const { return mFlags & isSeedForDirection; }

  /// Set the usedInTrack flag for the Tracker.
  inline void setIsUsed() { mFlags |= usedInTrack; }

  /// Set the isCandidate flag for the Tracker.
  inline void setIsCandidate() { mFlags |= isCandidateForTrack; }

  /// Set the isSeedForDirection flag for the Tracker.
  inline void setIsSeed() { mFlags |= isCandidateForTrack; }

  inline void setQTot           (Double_t &qTot)           { mQTot = qTot; }
  inline void setQMax           (Double_t &qMax)           { mQMax = qMax; }
  inline void setMeanFiber      (Float_t &meanFiber)       { mMeanFiber = meanFiber; }
  inline void setSigmaFiber     (Float_t &sigmaFiber)      { mSigmaFiber = sigmaFiber; }
  inline void setMeanTimeStamp  (Double_t &meanTimeStamp)  { mMeanTimeStamp = meanTimeStamp; }
  inline void setSigmaTimeStamp (Double_t &sigmaTimeStamp) { mSigmaTimeStamp = sigmaTimeStamp; }
  inline void setFirstTimeStamp (Double_t &firstTimeStamp) { mFirstTimeStamp = firstTimeStamp; }
  inline void setLayer          (Int_t &layer)             { mLayer = layer; }
  inline void setTDCID          (Int_t &tdcID)             { mTDCID = tdcID; }

  inline Double_t getQTot()           { return mQTot; }
  inline Double_t getQMax()           { return mQMax; }
  inline Float_t  getMeanFiber()      { return mMeanFiber; }
  inline Float_t  getSigmaFiber()     { return mSigmaFiber; }
  inline Double_t getMeanTimeStamp()  { return mMeanTimeStamp; }
  inline Double_t getSigmaTimeStamp() { return mSigmaTimeStamp; }
  inline Double_t getFirstTimeStamp() { return mFirstTimeStamp; }
  inline Int_t    getLayer()          { return mLayer; }
  inline Int_t    getTDCID()          { return mTDCID; }

  /// add a signal to the cluster
  /// calculates all cluster properties and updates members
  /// @param Signal
  /// @todo Figure out how to properly calculate the sigmas of fiber and timeStamp
  void addSignal(const Signal &signal);

  /// get the overall number of signals in the cluster
  inline Int_t getNSignals() { return mSignals.size(); }

  /// returns the vector containing the signals in the cluster
  inline std::vector<Signal>& getSignals() { return mSignals; }

 private:

  std::vector<Signal>  mSignals{};           ///< holds all signals in the cluster
  Double_t             mQTot = 0;            ///< sum of all ToT values (ToT corresponds to charge)
  Double_t             mQMax = 0;            ///< largest ToT among signals in the cluster
  Float_t              mMeanFiber = 0;       ///< arithmetic mean of fiber numbers weighted with the ToT of the corresponding signal
  Float_t              mSigmaFiber = 0;      ///< standard deviation
  Double_t             mMeanTimeStamp = 0;   ///< arithmetic mean of time stamps weighted with the ToT of the corresponding signal
  Double_t             mSigmaTimeStamp = 0;  ///< standard deviation
  Double_t             mFirstTimeStamp = 0;  ///< earliest time stamp among signals
  Int_t                mLayer = 0;           ///< layer in which the cluster is located
  Int_t                mTDCID = 0;           ///< TDC that is connected to the layer
  Short_t              mFlags = 0;           ///< Associated to track or not. This might be useful for tracking.

  ClassDef(Cluster,2);
};

#endif