#ifndef TRACKRECO_PHACTSTRKPROP_H
#define TRACKRECO_PHACTSTRKPROP_H

#include "PHTrackPropagating.h"

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>

#include <Acts/Utilities/BinnedArray.hpp>
#include <Acts/Utilities/Definitions.hpp>
#include <Acts/Utilities/Logger.hpp>

//#include <Acts/Propagator/detail/DebugOutputActor.hpp>
//#include <Acts/Propagator/detail/StandardAborters.hpp>
#include <Acts/Propagator/detail/SteppingLogger.hpp>
#include <Acts/Propagator/MaterialInteractor.hpp>

#include <ACTFW/EventData/TrkrClusterSourceLink.hpp>
#include <ACTFW/EventData/Track.hpp>

#include <memory>
#include <string>
#include <map>

struct ActsTrack;

class MakeActsGeometry;
class SvtxTrack;
class SvtxTrackMap;

namespace Acts
{
  class Surface;
  class PerigeeSurface;
}




using RecordedMaterial = Acts::MaterialInteractor::result_type;
using PropagationOutput
    = std::pair<std::vector<Acts::detail::Step>, RecordedMaterial>;

using PerigeeSurface = std::shared_ptr<const Acts::PerigeeSurface>;
using Surface = std::shared_ptr<const Acts::Surface>;
using SourceLink = FW::Data::TrkrClusterSourceLink;

class PHActsTrkProp : public PHTrackPropagating
{
 public:
  /// Default constructor
  PHActsTrkProp(const std::string& name = "PHActsTrkProp");

  /// Destructor
  ~PHActsTrkProp();

  /// End, write and close files
  int End();

  /// Get and create nodes
  int Setup(PHCompositeNode* topNode);

  /// Process each event by calling the fitter
  int Process();

 private:
  /// Event counter
  int m_event;

  /// Get all the nodes
  int getNodes(PHCompositeNode *topNode);

  /// Create new nodes
  void createNodes(PHCompositeNode *topNode);

  Acts::BoundSymMatrix getActsCovMatrix(const SvtxTrack *track);

  /// Run the propagation algorithm in acts. Returns result of propagation
  PropagationOutput propagate(FW::TrackParameters parameters);

  /// The acts geometry constructed in MakeActsGeometry
  MakeActsGeometry *m_actsGeometry;

  /// Minimum track pT to propagate, for acts propagator, units of GeV
  double m_minTrackPt;

  /// Propagation step size, units of mm
  double m_maxStepSize; 

  /// Track map with Svtx objects
  SvtxTrackMap *m_trackMap;

  /// Cluster-surface association maps created in MakeActsGeometry
  std::map<TrkrDefs::hitsetkey, Surface> *m_clusterSurfaceMap;
  std::map<TrkrDefs::cluskey, Surface> *m_clusterSurfaceTpcMap;

  /// Acts proto tracks to be put on the node tree by this module
  std::vector<ActsTrack> *m_actsProtoTracks;

  /// Acts source links created by PHActsSourceLinks
  /// SourceLink is defined as TrkrClusterSourceLink elsewhere
  std::map<unsigned int, SourceLink> *m_sourceLinks;

 
};

#endif
