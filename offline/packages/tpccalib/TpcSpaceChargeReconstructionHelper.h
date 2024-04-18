#ifndef TPCCALIB_TPCSPACECHARGERECONSTRUCTIONHELPER_H
#define TPCCALIB_TPCSPACECHARGERECONSTRUCTIONHELPER_H
/**
 * \file TpcSpaceChargeReconstructionHelper.h
 * \brief performs simple histogram manipulations for generating space charge distortion map suitable for correcting TPC clusters
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include <array>

// forward declarations
class TH2;
class TH3;
class TString;

class TpcSpaceChargeReconstructionHelper
{

  public:

  /// create TPOT acceptance mask, using provided binning
  /**
   * this histogram contains 1 in cells that match the TPOT acceptance
   * only Micromegas in the central sector (TPC sectors 9 and 21) are considered
   */
  static void create_tpot_mask( TH3* /*source*/ );

  /// z extrapolation
  /**
   * interpolate between micromegas in the fully equiped sector
   */
  static void extrapolate_z( TH3* /*source*/, TH3* /*mask*/ );

  /// first phi extrapolation
  /**
   * copy the full z dependence of reference sector to all other sectors, separately for positive and negative z,
   * normalized by the measurement from provided micromegas, at the appropriate z
   */
  static void extrapolate_phi1( TH3* /*source*/, TH2* /*source_cm*/, TH3* /*mask*/ );

  /// second phi extrapolation
  /**
   * for each r, z and phi bin, linearly extrapolate between neighbor phi sector measurements
   */
  static void extrapolate_phi2(  TH3* /*source*/, TH3* /*mask*/ );

  /// separate positive and negative z histograms
  /**
   * split histograms in two, the first with negative z values only, the second with positive z values
   * this must be done before adding guarding bins around each axis, in order to prevent artifacts during calls to Interpolate
   * at the central membrane (z = 0)
   */
  static std::tuple<TH3*, TH3*> split( TH3* hin );

  /**
   * copy input histogram into output, with new name, while adding two "guarding bins" on
   * each axis, with identical content and error as the first and last bin of the original histogram
   * this is necessary for being able to call TH3->Interpolate() when using these histograms
   * to correct for the space charge distortions.
   */
  static TH3* add_guarding_bins( TH3* hin, const TString& name );

};

#endif
