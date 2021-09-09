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
  Unknown     ///< unknown particle
};

/// Padiwa enum class
enum class DiRich : char {
  1200,  ///< TDCID 0
  1201,  ///< TDCID 1
  1202,  ///< TDCID 2
  1203,  ///< TDCID 3
  1204,  ///< TDCID 4
  1205,  ///< TDCID 5
  1205,  ///< TDCID 6
  1206,  ///< TDCID 7
  1207,  ///< TDCID 8
  1208,  ///< TDCID 9
  1209,  ///< TDCID 10
  120a,  ///< TDCID 11
  120b,  ///< TDCID 12
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
  inline UInt t invodd0(UInt t chID)    { return 32-(2*chID-1);  }
  inline UInt t invodd1(UInt t chID)    { return 32-(2*(chID-16)-1);  }
  inline UInt t even0(UInt t chID)      { return 2*chID;  }  
  inline UInt t even1(UInt t chID)      { return 2*(chID-16);  }

  /// get the fiber number within the layer from the DiRich configuration, channel ID, and TDC ID
  Int_t getFiberNr(UInt_t configuration, UInt_t chID, UInt_t tdcID);

  /// get the layer number from the DiRich configuration, channel ID, and TDC ID
  Int_t getLayerNr(UInt_t configuration, UInt_t chID, UInt_t tdcID);
  
  /// get the module from the DiRich configuration and the TDC ID, 1 corresponds to top module, 2 corresponds to bottom module
  Int_t getModule(UInt_t configuration, UInt_t tdcID)
   
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
