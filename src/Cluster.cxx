#include "Cluster.h"
#include <cmath>

ClassImp(Cluster);

//________________________________________________________________________________
Cluster::Cluster(const Cluster &cluster)
: mSignals(cluster.mSignals),
  mQTot(cluster.mQTot),
  mQMax(cluster.mQMax),
  mMeanFiber(cluster.mMeanFiber),
  mSigmaFiber(cluster.mSigmaFiber),
  mMeanTimeStamp(cluster.mMeanTimeStamp),
  mSigmaTimeStamp(cluster.mSigmaTimeStamp),
  mFirstTimeStamp(cluster.mFirstTimeStamp),
  mLayer(cluster.mLayer),
  mTDCID(cluster.mTDCID),
  mFlags(cluster.mFlags)
{
  
}

//________________________________________________________________________________
Cluster::Cluster(const Double_t qTot, const Double_t qMax, const Float_t meanFiber, const Float_t sigmaFiber, const Double_t meanTimeStamp,
                 const Double_t sigmaTimeStamp, const Double_t firstTimeStamp, const Int_t layer, const Int_t TDCID, const Short_t flags)
: mQTot(qTot),
  mQMax(qMax),
  mMeanFiber(meanFiber),
  mSigmaFiber(sigmaFiber),
  mMeanTimeStamp(meanTimeStamp),
  mSigmaTimeStamp(sigmaTimeStamp),
  mFirstTimeStamp(firstTimeStamp),
  mLayer(layer),
  mTDCID(TDCID),
  mFlags(flags)
{
  
}

//________________________________________________________________________________
void Cluster::addSignal(const Signal &signal)
{
  Double_t  meanFiber = 0;
  Double_t  meanTimeStamp = 0;

  Double_t weights = 0;

  mSignals.emplace_back(signal);

  for (auto &sig : mSignals) {
    if (sig.getLayer() != signal.getLayer()) {
      printf("Signal is from another layer than the previous ones added to the cluster!\n Signal not added to the cluster.\n");
      mSignals.pop_back();
      return;
    }

    else {
      weights += sig.getToT();
    }
  }

  for (auto &sig : mSignals) {
    Int_t fiberNr = mapping::getFiberNr(sig.getConfiguration(), sig.getChannelID(), sig.getTDCID());
    meanFiber     += (fiberNr*sig.getToT())/weights;
    meanTimeStamp += (sig.getTimeStamp()*sig.getToT())/weights;
  }

  mLayer = signal.getLayer();
  mTDCID = signal.getTDCID();
  mQTot += signal.getToT();
  mMeanFiber = meanFiber;
  mMeanTimeStamp = meanTimeStamp;
  mSigmaFiber = 1/std::sqrt(weights);
  mSigmaTimeStamp = 1/std::sqrt(weights);

  if (signal.getToT() > mQMax) { mQMax = signal.getToT(); }

  if (((signal.getTimeStamp() < mFirstTimeStamp) && (mFirstTimeStamp != 0)) || (mFirstTimeStamp == 0))  { mFirstTimeStamp = signal.getTimeStamp(); }
}
