#include "ClusterEvent.h"

ClassImp(ClusterEvent);

//________________________________________________________________________________
ClusterEvent::ClusterEvent(const ClusterEvent &event)
: mClusterVec(event.mClusterVec)
{
  this->setEventNr(event.getEventNr());
  this->setPadiwaConfig(event.getPadiwaConfig());
}