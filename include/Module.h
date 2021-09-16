#ifndef MODULE_H
#define MODULE_H

#include "Fiber.h"
#include "Utility.h"

/// The Module class represents the whole module. 
/// It has a vector member holding 256 Fiber Objects (whole module).
///
/// Before an event is written to file the removeEmpty funtion should be
/// called on the module to reduce data size and increase performance.

class Module
{
 public:
  Module();
  ~Module() = default;
  Module(const Module &module);

  /// add a signal to the module
  void addSignal(Signal &signal);

  /// get the overall number of signals in the module
  Float_t getNSignals();

  /// get the number of fibers that have a signal
  Int_t   getNFibers();

  /// get the fiber vector
  inline std::vector<Fiber>& getFibers() { return mFibers; }

  /// remove all empty fibers from the module
  void removeEmpty();

  /// remove all fibers from the module and make a fresh init()
  void reset();

  /// remove all signals from the module
  void refresh();

 private:

  void init();

  std::vector<Fiber>  mFibers{}; ///< holds all fibers of the module

  ClassDef(Module,1);
};

#endif