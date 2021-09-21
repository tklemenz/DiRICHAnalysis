#ifndef CTSEVENT_H
#define CTSEVENT_H

#include "EventBase.h"
#include "Module.h"

/// This class represents a basic CTS event containing signals.
///
/// The module objects contain all signal information.
///
/// Before an event is written to file the removeEmpty funtion should be
/// called on the module to reduce data size and increase performance.

class CTSEvent : public EventBase
{
 public:
  CTSEvent() = default;
  ~CTSEvent() = default;
  CTSEvent(const CTSEvent &event);

  inline void setModule (std::vector<Module> &moduleVec)  { mModuleVec = moduleVec; }

  std::vector<Module>&       getModules()       { return mModuleVec; }
  const std::vector<Module>& getModules() const { return mModuleVec; }

 private:
  std::vector<Module> mModuleVec;  ///< contains all fibers and signals

  ClassDef(CTSEvent,2);
};

#endif