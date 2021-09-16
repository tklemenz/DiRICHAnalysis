#ifndef CTSEVENT_H
#define CTSEVENT_H

#include "EventBase.h"
#include "Module.h"

/// This class represents a basic CTS event containing signals.
///
/// The module object contains all signal information.
///
/// Before an event is written to file the removeEmpty funtion should be
/// called on the module to reduce data size and increase performance.

class CTSEvent : public EventBase
{
 public:
  CTSEvent() = default;
  ~CTSEvent() = default;
  CTSEvent(const CTSEvent &event);

  inline void setModule (Module &module)  { mModule = module; }

  Module&       getModule()       { return mModule; }
  const Module& getModule() const { return mModule; }

 private:
  Module mModule;  ///< contains all fibers and signals

  ClassDef(CTSEvent,1);
};

#endif