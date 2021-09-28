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

  inline void setModules(std::vector<Module> &moduleVec) { mModuleVec = moduleVec; }

  inline void setNSignals(ULong_t nSignals) { mNSignals = nSignals; }

  std::vector<Module>&       getModules()       { return mModuleVec; }
  const std::vector<Module>& getModules() const { return mModuleVec; }

  ULong_t&       getNSignals()       { return mNSignals; }
  const ULong_t& getNSignals() const { return mNSignals; }

 private:
  std::vector<Module> mModuleVec;  ///< contains all fibers and signals
  ULong_t mNSignals;               ///< number of signals in event; only counts signals with signalNr == 1

  ClassDef(CTSEvent,3);
};

#endif