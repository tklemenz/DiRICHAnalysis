#ifndef UTILITY_H
#define UTILITY_H

#include <Rtypes.h>
#include <utility>

class TObjArray;
class TString;
class TChain;

/// ParticleType enum class
enum class ParticleType : char {
  Pion,       ///< particle is a pion
  Proton,     ///< particle is a proton
  Cosmic,     ///< particle is a cosmic (lab)
  Unknown     ///< unknown particle
};

/// DiRICH enum class
enum class DiRich : char {
  d1200,  ///< TDCID 0
  d1201,  ///< TDCID 1
  d1202,  ///< TDCID 2
  d1203,  ///< TDCID 3
  d1204,  ///< TDCID 4
  d1205,  ///< TDCID 5
  d1206,  ///< TDCID 6
  d1207,  ///< TDCID 7
  d1208,  ///< TDCID 8
  d1209,  ///< TDCID 9
  d120a,  ///< TDCID 10
  d120b,  ///< TDCID 11
};

/// Padiwa enum class
enum class Padiwa : char {
  p1500_0,    ///< TDC 0, ch 1-16
  p1500_1,    ///< TDC 0, ch 17-32
  p1510_0,    ///< TDC 1, ch 1-16
  p1510_1,    ///< TDC 1, ch 17-32
  p1520_0,    ///< TDC 2, ch 1-16
  p1520_1,    ///< TDC 2, ch 17-32
  p1530_0,    ///< TDC 3, ch 1-16
  p1530_1     ///< TDC 3, ch 17-32
};

enum class PadiwaSocket : char {
  L1even,
  L2even,
  L3even,
  L4even,
  L5even,
  L6even,
  L7even,
  L8even,
  L1odd,
  L2odd,
  L3odd,
  L4odd,
  L5odd,
  L6odd,
  L7odd,
  L8odd,
};

/// This namespace holds the mapping functions.
namespace mapping
{

  /// get the fiber number (0 - 255) that the layer+fiber combination takes in the module 
  inline Int_t getModuleSpot(Int_t layer, Int_t fiber) { return ((fiber + (layer-1)*32)-1); }

  /// get the layer and fiber number from the module spot
  /// returns std::pair<layer,fiber>
  std::pair<Int_t, Int_t> getFiberInfoFromModSpot(Int_t modSpot);

  /// mapping functions from channel number to fiber
  /// invodd maps channel number to odd numbers from 31-1
  /// even maps channel number to even numbers from 2-32
  /// 0 is for channel number 1-16
  /// 1 is for channel number 17-32
  inline UInt_t invodd0(UInt_t chID)    { return 32-(2*chID-1); }
  inline UInt_t invodd1(UInt_t chID)    { return 32-(2*(chID-16)-1); }
  inline UInt_t even0(UInt_t chID)      { return 2*chID; }
  inline UInt_t even1(UInt_t chID)      { return 2*(chID-16); }

  /// get the fiber number within the layer from the DiRich configuration, channel ID, and TDC ID
  Int_t getFiberNr(UInt_t configuration, UInt_t chID, UInt_t tdcID);

  /// get the layer number from the DiRich configuration, channel ID, and TDC ID
  Int_t getLayerNr(UInt_t configuration, UInt_t chID, UInt_t tdcID);

  /// get the fiber number in x-direction from the fiber number and layer number. returns 0 for even layers.
  Int_t getX(Int_t layer, Int_t fiber);

  /// get the fiber number in y-direction from the fiber number and layer number. returns 0 for odd layers.
  Int_t getY(Int_t layer, Int_t fiber);
  
  /// get the module from the DiRich configuration and the TDC ID, 1 corresponds to top module, 2 corresponds to bottom module
  Int_t getModule(UInt_t configuration, UInt_t tdcID);
   
  /// get the actual spatial coordinate in mm
  inline Float_t getCoord(Float_t meanFiber) { return 2*meanFiber-1; }



} // namespace mapping

namespace beautify
{
  void setStyle();

  template <class T>
  void setStyleHisto(T* histo);

} // namespace beautify

/// This namespace holds some useless text modifications for terminal output.
namespace text
{

  const char* const BLK   = "\e[30m";
  const char* const RED   = "\e[31m";
  const char* const GRN   = "\e[32m";
  const char* const YEL   = "\e[33m";
  const char* const BLU   = "\e[34m";
  const char* const MAG   = "\e[35m";
  const char* const CYN   = "\e[36m";
  const char* const GRY   = "\e[37m";
  const char* const LBLK  = "\e[90m";
  const char* const LRED  = "\e[91m";
  const char* const LGRN  = "\e[92m";
  const char* const LYEL  = "\e[93m";
  const char* const LBLU  = "\e[94m";
  const char* const LMAG  = "\e[95m";
  const char* const LCYN  = "\e[96m";
  const char* const WHT   = "\e[97m";
  const char* const RESET = "\e[0m";

  const char* const BOLD  = "\e[1m";
  const char* const ULINE = "\e[4m";
  const char* const BLINK = "\e[5m";

} // namespace text

namespace fileHandling
{

  void makeChain(TChain& chain, const TString& input);

  std::vector<std::string> splitString(std::string inString, const char* delimiter = "/");

} // namespace fileHandling

#endif
