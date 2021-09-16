#include "Fiber.h"

ClassImp(Fiber);

//________________________________________________________________________________
Fiber::Fiber(const Fiber &fiber)
: mSignalVec(fiber.mSignalVec),
  mLayer(fiber.mLayer),
  mX(fiber.mX),
  mY(fiber.mY)
{
  
}

//________________________________________________________________________________
Fiber::Fiber(const Int_t layer, const Int_t x, const Int_t y)
: mLayer(layer),
  mX(x),
  mY(y)
{
  
}