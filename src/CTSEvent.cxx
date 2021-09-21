#include "CTSEvent.h"

ClassImp(CTSEvent);

//________________________________________________________________________________
CTSEvent::CTSEvent(const CTSEvent &event)
: mModuleVec(event.mModuleVec)
{
  this->setEventNr(event.getEventNr());
  this->setPadiwaConfig(event.getPadiwaConfig());
}