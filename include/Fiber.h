#ifndef FIBER_H
#define FIBER_H

#include <vector>

#include "Signal.h"

/// The Fiber class represents a fiber in the module.
/// It holds:
///           - a vector containing all signals in the fiber
///           - layer, x-, y-coordinate

class Fiber
{
 public:
  Fiber() = default;
  ~Fiber() = default;
  Fiber(const Fiber &fiber);
  Fiber(const Int_t layer, const Int_t x, const Int_t y);

  inline void setLayer   (Int_t layer)                    { mLayer = layer; }
  inline void setX       (Int_t x)                        { mX = x; }
  inline void setY       (Int_t y)                        { mY = y; }
  inline void setSignals (std::vector<Signal> &signalVec) { mSignalVec = signalVec; }              /// set the whole signal vector

  inline void addSignal  (Signal &signal)                 { mSignalVec.emplace_back(signal); }    /// add a single signal to the fiber

  inline Int_t getLayer()    const { return mLayer; }
  inline Int_t getX()        const { return mX; }
  inline Int_t getY()        const { return mY; }
  inline Int_t getNSignals() const { return mSignalVec.size(); }

  std::vector<Signal>&       getSignals()       { return mSignalVec; }
  const std::vector<Signal>& getSignals() const { return mSignalVec; }

  /// clears the signal vector
  void reset() { mSignalVec.clear(); }

 private:

  std::vector<Signal> mSignalVec{}; ///< holds all signals recorded in the fiber
  Int_t               mLayer;       ///< layer the fiber is located in
  Int_t               mX;           ///< xth fiber in vertical (odd) layers, counts from left to right
  Int_t               mY;           ///< yth fiber in horizontal (even) layers, counts from bottom to top

  ClassDef(Fiber,1);
};

#endif