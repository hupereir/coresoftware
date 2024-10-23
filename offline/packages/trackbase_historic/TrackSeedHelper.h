#ifndef TRACKBASEHISTORIC_TRACKSEEDHELPER_H
#define TRACKBASEHISTORIC_TRACKSEEDHELPER_H

/*!
 * \file TrackSeedHelper.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@lanl.gov>
 */

#include <trackbase/TrkrDefs.h>

#include <Acts/Definitions/Algebra.hpp>

#include <cstdint>
#include <map>

class TrackSeed;

class TrackSeedHelper
{

  public:

  using position_map_t = std::map<TrkrDefs::cluskey, Acts::Vector3>;

  static float get_phi(TrackSeed const*, const position_map_t&);
  static float get_phi_fastsim(TrackSeed const*);

  static void circleFitByTaubin(
    TrackSeed*, const position_map_t& positions,
    uint8_t startLayer = 0,
    uint8_t endLayer = 58);

  static void lineFit(
    TrackSeed*, const position_map_t& positions,
    uint8_t startLayer = 0,
    uint8_t endLayer = 58);

  static float get_x(TrackSeed const*);
  static float get_y(TrackSeed const*);
  static float get_z(TrackSeed const*);
  static Acts::Vector3 get_xyz(TrackSeed const*);

  protected:
  static std::pair<float, float> findRoot(TrackSeed const*);

};



#endif
