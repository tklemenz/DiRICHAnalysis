#include "CTSEventTracks.h"

ClassImp(CTSEventTracks);

//________________________________________________________________________________
CTSEventTracks::CTSEventTracks(const CTSEventTracks &event)
: mTrackVec(event.mTrackVec)
{
  this->setEventNr(event.getEventNr());
  this->setPadiwaConfig(event.getPadiwaConfig());
}