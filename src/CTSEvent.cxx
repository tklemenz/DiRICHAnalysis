#include "CTSEvent.h"

ClassImp(CTSEvent);

//________________________________________________________________________________
CTSEvent::CTSEvent(const CTSEvent &event)
: mModule(event.mModule)
{
  this->setEventNr(event.getEventNr());
  this->setPadiwaConfig(event.getPadiwaConfig());
}