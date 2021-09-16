#include "CTSEventClusters.h"

ClassImp(CTSEventClusters);

//________________________________________________________________________________
CTSEventClusters::CTSEventClusters(const CTSEventClusters &event)
: mClusterVec(event.mClusterVec)
{
  this->setEventNr(event.getEventNr());
  this->setPadiwaConfig(event.getPadiwaConfig());
}