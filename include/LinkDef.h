#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class Dummy +;

#pragma link C++ function mapping::getModuleSpot +;
#pragma link C++ function mapping::getFiberInfoFromModSpot +;
#pragma link C++ function mapping::invodd0 +;
#pragma link C++ function mapping::invodd1 +;
#pragma link C++ function mapping::even0 +;
#pragma link C++ function mapping::even1 +;
#pragma link C++ function mapping::getFiberNr +;
#pragma link C++ function mapping::getLayerNr +;
#pragma link C++ function mapping::getCoord +;
#pragma link C++ function mapping::getModule +;

#pragma link C++ function beautify::setStyle +;
#pragma link C++ function beautify::setStyleHisto +;

#endif
