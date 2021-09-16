#include <algorithm>

#include "Module.h"

ClassImp(Module);

Module::Module()
{
  Module::init();
}

//________________________________________________________________________________
Module::Module(const Module &module)
: mFibers(module.mFibers)
{
  
}

//________________________________________________________________________________
void Module::addSignal(Signal &signal)
{
  mFibers.at(mapping::getModuleSpot(signal.getLayer(), mapping::getFiberNr(signal.getConfiguration(),signal.getChannelID(),signal.getTDCID()))).addSignal(signal);
}

//________________________________________________________________________________
Float_t Module::getNSignals()
{
  Float_t nSignals = 0;

  for (auto &fiber : mFibers) {
    nSignals += fiber.getNSignals();
  }

  return nSignals;
}

//________________________________________________________________________________
Int_t Module::getNFibers()
{
  Int_t nFibers = 0;

  for (auto &fiber : mFibers) {
    if (fiber.getNSignals() > 0) {
      nFibers += 1;
    }
  }

  return nFibers;
}

//________________________________________________________________________________
void Module::removeEmpty()
{
  auto end = std::remove_if(mFibers.begin(), mFibers.end(), [](Fiber &fiber) {return fiber.getNSignals() == 0;});

  mFibers.erase(end, mFibers.end());
}

//________________________________________________________________________________
void Module::reset()
{
  mFibers.clear();
  Module::init();
}

//________________________________________________________________________________
void Module::refresh()
{
  for(auto& fiber : mFibers) {
    fiber.reset();
  }
}

//________________________________________________________________________________
void Module::init()
{
  for(Int_t i=0; i<256; i++) {

    std::pair<Int_t, Int_t> layFib = mapping::getFiberInfoFromModSpot(i);

    mFibers.emplace_back(Fiber(layFib.first, mapping::getX(layFib.first, layFib.second), mapping::getY(layFib.first, layFib.second)));
  }
}