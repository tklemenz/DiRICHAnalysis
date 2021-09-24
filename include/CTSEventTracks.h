#ifndef CTSEVENTTRACKS_H
#define CTSEVENTTRACKS_H

#include "EventBase.h"
#include "Track.h"

/// This class represents a CTS event containing tracks.
/// It holds:
///          - a vector containing all tracks in the event

class CTSEventTracks : public EventBase
{
 public:
  CTSEventTracks() = default;
  ~CTSEventTracks() = default;
  CTSEventTracks(const CTSEventTracks &event);

  inline void addTrack  (Track &track)  { mTrackVec.emplace_back(track); }

  inline void setTracks (std::vector<Track> &tracks)  { mTrackVec = tracks; }

  std::vector<Track>&       getTracks()       { return mTrackVec; }
  const std::vector<Track>& getTracks() const { return mTrackVec; }

 private:
  std::vector<Track> mTrackVec {};  ///< contains all tracks in the event

  ClassDef(CTSEventTracks,1);
};

#endif