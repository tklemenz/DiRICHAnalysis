#ifndef EVENTBASE_H
#define EVENTBASE_H

#include "Fiber.h"

/// This class represents a basic event from which specified event classes 
/// (e.g. CTSEvent or physical event classes) can be derived.
/// 
/// The philosophy still needs a bit of thought...
/// The basic information the event holds is:
///                                          - EventNr
///                                          - PadiwaConfig

class EventBase
{
 public:
  EventBase() = default;
  ~EventBase() = default;
  EventBase(const EventBase &event);

  inline void setEventNr      (Float_t eventNr) { mEventNr = eventNr; }
  inline void setPadiwaConfig (Int_t config)    { mPadiwaConfig = config; }

  inline Float_t getEventNr()      const { return mEventNr; }
  inline Int_t   getPadiwaConfig() const { return mPadiwaConfig; }

 private:

  Float_t mEventNr;       ///< event number
  Int_t   mPadiwaConfig;  ///< padiwa configuration, 0-4

  ClassDef(EventBase,1);
};

#endif