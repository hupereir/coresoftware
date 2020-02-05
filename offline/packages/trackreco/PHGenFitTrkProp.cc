/*!
 *  \file PHGenFitTrkProp.C
 *  \brief Progressive pattern recgnition based on GenFit Kalman filter
 *  \detail using Alan Dion's HelixHough for seeding, GenFit Kalman filter to do track propagation
 *  \author Christof Roland & Haiwang Yu
 */

#include "PHGenFitTrkProp.h"

#include "AssocInfoContainer.h"

#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxVertex.h>
#include <trackbase_historic/SvtxVertexMap.h>

// sPHENIX Geant4 includes
#include <g4detectors/PHG4CylinderCellGeom.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>
#include <g4detectors/PHG4CylinderGeom.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>
//

#include <intt/CylinderGeomIntt.h>
#include <intt/InttDefs.h>

#include <mvtx/CylinderGeom_Mvtx.h>
#include <mvtx/MvtxDefs.h>

#include <g4bbc/BbcVertexMap.h>

#include <phfield/PHFieldUtility.h>

#include <phgeom/PHGeomUtility.h>

// sPHENIX includes
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/PHTimer.h>
#include <phool/phool.h>

// GenFit
#include <GenFit/EventDisplay.h>
#include <GenFit/MeasuredStateOnPlane.h>
#include <GenFit/RKTrackRep.h>
#include <GenFit/Track.h>

#include <phgenfit/Fitter.h>
#include <phgenfit/Measurement.h>
#include <phgenfit/PlanarMeasurement.h>
#include <phgenfit/Track.h>

//ROOT includes for debugging
#include <TFile.h>
#include <TMatrixDSymfwd.h>
#include <TMatrixTSym.h>
#include <TMatrixTUtils.h>
#include <TNtuple.h>
#include <TVector3.h>
#include <TVectorDfwd.h>
#include <TVectorT.h>

// standard includes
#include <cassert>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>

#define LogDebug(exp) std::cout << "DEBUG: " << __FILE__ << ": " << __LINE__ << ": " << exp
#define LogError(exp) std::cout << "ERROR: " << __FILE__ << ": " << __LINE__ << ": " << exp
#define LogWarning(exp) std::cout << "WARNING: " << __FILE__ << ": " << __LINE__ << ": " << exp

//#define _DEBUG_
//#define _USE_CONSTANT_SEARCH_WIN_
//#define _DO_FULL_FITTING_

#ifdef _DEBUG_
ofstream fout_kalman_pull("kalman_pull.txt");
ofstream fout_chi2("chi2.txt");
#endif

//___________________________________________
template< class T > T square( T x ) { return x*x; }

//___________________________________________
PHGenFitTrkProp::PHGenFitTrkProp(
    const std::string& name,
    unsigned int nlayers_maps,
    unsigned int nlayers_intt,
    unsigned int nlayers_tpc,
    unsigned int nlayers_outer)
  : PHTrackPropagating(name)
  , _nlayers_maps(nlayers_maps)
  , _nlayers_intt(nlayers_intt)
  , _nlayers_tpc(nlayers_tpc)
  , _nlayers_outer(nlayers_outer)
  , _nlayers_all(_nlayers_maps + _nlayers_intt + _nlayers_tpc + _nlayers_outer)
  , _firstlayer_maps(0)
  , _firstlayer_intt(_firstlayer_maps + _nlayers_maps)
  , _firstlayer_tpc(_firstlayer_intt + _nlayers_intt)
  , _firstlayer_outer(_firstlayer_tpc + _nlayers_tpc)
{}

//___________________________________________
int PHGenFitTrkProp::Setup(PHCompositeNode* topNode)
{

  int ret = PHTrackPropagating::Setup(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  ret = InitializeGeometry(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK)
    return ret;

  ret = InitializePHGenFit(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  if (_analyzing_mode)
  {
    std::cout << "Ana Mode, creating ntuples! " << std::endl;
    _analyzing_file = new TFile("./PHGenFitTrkProp.root", "RECREATE");
    _analyzing_ntuple = new TNtuple("ana_nt", "ana_nt", "pt:kappa:d:phi:dzdl:z0:nhit:ml:rec:dt");
    std::cout << "Done" << std::endl;
  }

  /*!
   * Initilize parameters
   */
  for (int layer = 0; layer < _nlayers_all; ++layer)
  {
    _search_wins_phi.insert(std::make_pair(layer, _search_win_phi));
    _search_wins_theta.insert(std::make_pair(layer, _search_win_theta));
    _max_incr_chi2s.insert(std::make_pair(layer, _max_incr_chi2));
  }

  #ifdef _DEBUG_
  for (int layer = 0; layer < _nlayers_all; ++layer)
  {
    std::cout
      << __LINE__
      << ": layer: " << layer
      << ": search_wins_rphi: " << _search_wins_phi[layer]
      << ": search_wins_z: " << _search_wins_theta[layer]
      << ": max_incr_chi2: " << _max_incr_chi2s[layer]
      << std::endl;
  }
  #endif

  _t_seeds_cleanup = new PHTimer("_t_seeds_cleanup");
  _t_seeds_cleanup->stop();

  _t_translate_to_PHGenFitTrack = new PHTimer("_t_translate_to_PHGenFitTrack");
  _t_translate_to_PHGenFitTrack->stop();

  _t_translate1 = new PHTimer("_t_translate1");
  _t_translate1->stop();
  _t_translate2 = new PHTimer("_t_translate2");
  _t_translate2->stop();
  _t_translate3 = new PHTimer("_t_translate3");
  _t_translate3->stop();

  _t_kalman_pat_rec = new PHTimer("_t_kalman_pat_rec");
  _t_kalman_pat_rec->stop();

  _t_search_clusters = new PHTimer("_t_search_clusters");
  _t_search_clusters->stop();

  _t_search_clusters_encoding = new PHTimer("_t_search_clusters_encoding");
  _t_search_clusters_encoding->stop();

  _t_search_clusters_map_iter = new PHTimer("_t_search_clusters_map_iter");
  _t_search_clusters_map_iter->stop();

  _t_track_propagation = new PHTimer("_t_track_propergation");
  _t_track_propagation->stop();

  _t_full_fitting = new PHTimer("_t_full_fitting");
  _t_full_fitting->stop();

  _t_output_io = new PHTimer("_t_output_io");
  _t_output_io->stop();

  return ret;
}

//_________________________________________________________________________
int PHGenFitTrkProp::Process()
{
  if (Verbosity() > 10)
  {
    std::cout << "PHGenFitTrkProp::process_event -- entered" << std::endl;
    std::cout << "nMapsLayers = " << _nlayers_maps << std::endl;
    std::cout << "nInttLayers = " << _nlayers_intt << std::endl;
    std::cout << "nTPCLayers = " << _nlayers_tpc << std::endl;
  }
  // start fresh
  _gftrk_hitkey_map.clear();

  _vertex.clear();
  _vertex_error.clear();

  if (_vertex_map)
  {
    for(unsigned int ivert=0; ivert<_vertex_map->size(); ++ivert)
    {
      auto vertex = _vertex_map->get(ivert);
      _vertex.push_back( { vertex->get_x(), vertex->get_y(), vertex->get_z() } );
      _vertex_error.push_back( {
        std::sqrt(vertex->get_error(0,0)),
        std::sqrt(vertex->get_error(1,1)),
        std::sqrt(vertex->get_error(2,2)) } );
    }
  }

  {
    //-----------------------------------
    // Kalman track propagating
    //-----------------------------------
    if (Verbosity() >= 1) _t_kalman_pat_rec->restart();
    int ret = KalmanTrkProp();
    if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;
    if (Verbosity() >= 1) _t_kalman_pat_rec->stop();
  }
  if (Verbosity() > 1) print_timers();

  ++_event;

  return Fun4AllReturnCodes::EVENT_OK;
}

//_________________________________________________________________________
void PHGenFitTrkProp::print_timers()
{
  std::cout << "=============== PHGenFitTrkProp::print_timers: ===============" << std::endl;
  std::cout << "\t - Seeds Cleanup:          " << _t_seeds_cleanup->get_accumulated_time() / 1000. << " sec" << std::endl;
  std::cout << "CPUSCALE Pattern recognition time:    " << _t_kalman_pat_rec->get_accumulated_time() / 1000. << " sec" << std::endl;
  std::cout << "\t - Track Translation time: " << _t_translate_to_PHGenFitTrack->get_accumulated_time() / 1000. << " sec" << std::endl;
  std::cout << "\t -    - Translation1 time: " << _t_translate1->get_accumulated_time() / 1000. << " sec" << std::endl;
  std::cout << "\t -    - Translation2 time: " << _t_translate2->get_accumulated_time() / 1000. << " sec" << std::endl;
  std::cout << "\t -    - Translation3 time: " << _t_translate3->get_accumulated_time() / 1000. << " sec" << std::endl;
  std::cout << "\t - Cluster searching time: " << _t_search_clusters->get_accumulated_time() / 1000. << " sec" << std::endl;
  std::cout << "\t\t - Encoding time:        " << _t_search_clusters_encoding->get_accumulated_time() / 1000. << " sec" << std::endl;
  std::cout << "\t\t - Map iteration:        " << _t_search_clusters_map_iter->get_accumulated_time() / 1000. << " sec" << std::endl;
  std::cout << "\t - Kalman updater time:    " << _t_track_propagation->get_accumulated_time() / 1000. << " sec" << std::endl;
  std::cout << "Full fitting time:           " << _t_full_fitting->get_accumulated_time() / 1000. << " sec" << std::endl;
  std::cout << "Output IO time:              " << _t_output_io->get_accumulated_time() / 1000. << " sec" << std::endl;
  std::cout << "=======================================" << std::endl;
}

//______________________________________________________________
int PHGenFitTrkProp::End()
{

  if (_do_evt_display) _fitter->displayEvent();
  delete _t_seeds_cleanup;
  delete _t_translate_to_PHGenFitTrack;
  delete _t_translate1;
  delete _t_translate2;
  delete _t_translate3;
  delete _t_kalman_pat_rec;
  delete _t_full_fitting;
  delete _t_search_clusters;
  delete _t_track_propagation;
  delete _t_output_io;

#ifdef _DEBUG_
  LogDebug("Leaving End \n");
#endif

#ifdef _DEBUG_
  fout_kalman_pull.close();
  fout_chi2.close();
#endif

  if (_analyzing_mode)
  {
    std::cout << " cleaning up " << std::endl;
    _analyzing_file->cd();
    _analyzing_ntuple->Write();
    _analyzing_file->Close();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//_______________________________________________________________________________-
int PHGenFitTrkProp::InitializePHGenFit(PHCompositeNode* topNode)
{

  // create fitter
  auto tgeo_manager = PHGeomUtility::GetTGeoManager(topNode);
  auto field = PHFieldUtility::GetFieldMapNode(nullptr, topNode);
  _fitter.reset( PHGenFit::Fitter::getInstance(tgeo_manager, field, _track_fitting_alg_name, "RKTrackRep", _do_evt_display) );

  if (!_fitter)
  {
    std::cerr << PHWHERE << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  #ifdef _DEBUG_
  _fitter->set_verbosity(10);
  #endif

  return Fun4AllReturnCodes::EVENT_OK;
}

//______________________________________________________________
int PHGenFitTrkProp::GetNodes(PHCompositeNode* topNode)
{
  _bbc_vertexes = findNode::getClass<BbcVertexMap>(topNode, "BbcVertexMap");
  _geom_container_intt = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_INTT");
  _geom_container_maps = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MVTX");
  return Fun4AllReturnCodes::EVENT_OK;
}

//______________________________________________________________
int PHGenFitTrkProp::check_track_exists(MapPHGenFitTrack::iterator iter, SvtxTrackMap::Iter phtrk_iter)
{
  //Loop over hitIDs on current track and check if they have been used
  unsigned int n_clu = iter->second->get_cluster_keys().size();
  unsigned int n_clu_used = 0;
  int n = 0;
  const auto& clusterkeys = iter->second->get_cluster_keys();
  for (TrkrDefs::cluskey iCluId = 0; iCluId < clusterkeys.size(); ++iCluId)
  {
    auto cluster_ID = clusterkeys[iCluId];

    if(_gftrk_hitkey_map.count(iCluId)>0) n_clu_used++;
    if (Verbosity() >= 10)
    {
      std::cout << " trk map size: " << _gftrk_hitkey_map.count(iCluId) << std::endl;
      std::cout << "n: " << n << "#Clu_g = " << iCluId << " ntrack match: "  << _assoc_container->GetTracksFromCluster(cluster_ID).size()
        << " layer: " << (float)TrkrDefs::getLayer(cluster_ID)
        << " r: " << std::sqrt( square( _cluster_map->findCluster(cluster_ID)->getX()) + square(_cluster_map->findCluster(cluster_ID)->getY()) )
        << " used: " << n_clu_used
        << std::endl;
      ++n;
    }
  }

  if (((float) n_clu_used / n_clu) > 0.3)
  {
    if (Verbosity() >= 1) std::cout << "Found duplicate track. n_clu: " << n_clu << " c_clu_used: " << n_clu_used << std::endl;
    return 0;
  }
  return 1;
}

//_________________________________________________________________________
int PHGenFitTrkProp::KalmanTrkProp()
{
  // _track_map is created in PHHoughSeeding::ExportOutput at line 1714 and written to SvtxTrack node.
  // PHTrackPropagating (base class of this one) gets it from the node tree and calls this "Process" method.
  // The process method then calls this method.

  #ifdef _DEBUG_
  std::cout << "=========================" << std::endl;
  std::cout << "PHGenFitTrkProp::KalmanTrkProp: Start: Event: " << _event << std::endl;
  std::cout << "Total Raw Tracks: " << _track_map->size() << std::endl;
  std::cout << "=========================" << std::endl;
  #endif

  // A map is needed for each vertex location. This is a vector of maps
  _layer_thetaID_phiID_clusterID.clear();
  /*  for(unsigned int ivert=0; ivert<_vertex.size(); ++ivert)
    {
      BuildLayerZPhiHitMap(ivert);
    }
  */
  BuildLayerZPhiHitMap(0);

  std::vector<genfit::Track*> evt_disp_copy;

  _PHGenFitTracks.clear();
  if (Verbosity() > 1)
  { std::cout << " found " << _track_map->size() << " track seeds " << std::endl; }

  for (auto phtrk_iter = _track_map->begin(); phtrk_iter != _track_map->end();)
  {
    auto tracklet = phtrk_iter->second;
    if (Verbosity() >= 10)
    {
      std::cout
        << __LINE__
        << ": Processing seed itrack: " << phtrk_iter->first
        << ": Total tracks: " << _track_map->size()
        << std::endl;
    }

    // Translate sPHENIX track To PHGenFitTracks
    if (Verbosity() > 1) _t_translate_to_PHGenFitTrack->restart();
    SvtxTrackToPHGenFitTracks(tracklet);

    // Handle track propagation, termination, output and evt disp.
    bool is_splitting_track = false;
    #ifdef _DEBUG_
    int i = 0;
    #endif
    if( _PHGenFitTracks.empty() )
    {
      std::cout << "Warning: Conversion of SvtxTrack tracklet " <<  phtrk_iter->first << " to PHGenFitTrack failed, moving to next tracklet " << std::endl;
      ++phtrk_iter;
      continue;
    }

    for (auto gftrk_iter = _PHGenFitTracks.begin(); gftrk_iter != _PHGenFitTracks.end(); ++gftrk_iter)
    {
      if(Verbosity() > 10)
      {
        std::cout
          << __LINE__
          << ": propagating Genfit track: " << phtrk_iter->first << std::endl;
      }

      // associate this track with the same vertex as the seed track
      unsigned int ivert = tracklet->get_vertex_id();
      gftrk_iter->second->set_vertex_id(ivert);
      if(ivert > _vertex.size())
      {
        std::cout << PHWHERE << " Track vertex is screwed up, have to quit! " << std::endl;
        return Fun4AllReturnCodes::ABORTRUN;
      }

      std::vector<TrkrDefs::cluskey> clusterkeys = gftrk_iter->second->get_cluster_keys();

      unsigned int init_layer = UINT_MAX;

      if (!is_splitting_track)
      {
        if( _init_direction == 1 )
        {
          init_layer = TrkrDefs::getLayer(clusterkeys.front());
          TrackPropPatRec(ivert, gftrk_iter, init_layer, _nlayers_all, true);
          TrackPropPatRec(ivert, gftrk_iter, init_layer, 0, false);
        } else {
          init_layer = TrkrDefs::getLayer(clusterkeys.back());
          TrackPropPatRec(ivert, gftrk_iter, init_layer, 0, true);
          TrackPropPatRec(ivert, gftrk_iter, init_layer, _nlayers_all, false);
        }
        is_splitting_track = true;

      } else {

        if( _init_direction == 1 )
        {
          init_layer = TrkrDefs::getLayer(clusterkeys.front());
          TrackPropPatRec(ivert, gftrk_iter, init_layer, _nlayers_all, false);
        } else {
          init_layer = TrkrDefs::getLayer(clusterkeys.back());
          TrackPropPatRec(ivert, gftrk_iter, init_layer, 0, false);
        }
      }

      #ifdef _DEBUG_
      std::cout
        << __LINE__
        << ": tracki: " << i
        << ": clusterkeys size:  " << gftrk_iter->second->get_cluster_keys().size()
        << ": quality: " << gftrk_iter->first
        << std::endl;
      ++i;
      #endif

    }

    #ifdef _DEBUG_
    i = 0;
    for (auto iter = _PHGenFitTracks.begin(); iter != _PHGenFitTracks.end(); ++iter)
    {
      std::cout
        << __LINE__
        << ": track: " << i++
        << ": clusterkeys size:  " << iter->second->get_cluster_keys().size()
        << ": quality: " << iter->first
        << std::endl;
    }
    #endif

    _PHGenFitTracks.sort();

    #ifdef _DEBUG_
    for (auto iter = _PHGenFitTracks.begin(); iter != _PHGenFitTracks.end(); ++iter)
    {
      std::cout
        << __LINE__
        << ": clusterkeys size:  " << iter->second->get_cluster_keys().size()
        << ": quality: " << iter->first
        << std::endl;
    }
    #endif

    auto gftrk_iter_best = _PHGenFitTracks.begin();
    int track_exists = check_track_exists(gftrk_iter_best,phtrk_iter);

    if (gftrk_iter_best->second->get_cluster_keys().size() >= _min_good_track_hits && track_exists)
    {
      OutputPHGenFitTrack(gftrk_iter_best, phtrk_iter);
      #ifdef _DEBUG_
      std::cout << __LINE__ << std::endl;
      #endif
      if (_do_evt_display)
      {
        evt_disp_copy.push_back( new genfit::Track(*gftrk_iter_best->second->getGenFitTrack()));
      }

    } else {

      auto key = phtrk_iter->first;
      ++phtrk_iter;
      _track_map->erase(key);
      continue;
    }

    _PHGenFitTracks.clear();
    ++phtrk_iter;
  }

  if(Verbosity() > 1)
  {
    std::cout << "=========================" << std::endl;
    std::cout << "PHGenFitTrkProp::KalmanTrkProp: End: Event: " << _event << std::endl;
    std::cout << "PHGenFitTrkProp::KalmanTrkProp: End: Event: " << _event << std::endl;
    std::cout << "Total Final Tracks: " << _track_map->size() << std::endl;
  }

  if (_do_evt_display) _fitter->getEventDisplay()->addEvent(evt_disp_copy);
  else evt_disp_copy.clear();

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int PHGenFitTrkProp::OutputPHGenFitTrack( MapPHGenFitTrack::iterator gftrk_iter, SvtxTrackMap::Iter phtrk_iter)
{
  if(Verbosity() > 10)
  {
    std::cout << "=========================" << std::endl;
    std::cout << __LINE__ << ": iPHGenFitTrack: " << phtrk_iter->first << std::endl;
    std::cout << __LINE__ << ": _track_map->size(): " << _track_map->size() << std::endl;
    std::cout << "Contains: " << gftrk_iter->second->get_cluster_keys().size() << " clusters." << std::endl;
    std::cout << "=========================" << std::endl;
  }

  // get the track from the node tree and extract the track and vertex id
  auto track = phtrk_iter->second;
  auto track_id = track->get_id();
  auto vertex_id = track->get_vertex_id();

  // now reset the track info and rewrite it using info from the Genfit track
  track->Reset();
  track->set_id(track_id);
  track->set_vertex_id(vertex_id);

  #ifdef _DO_FULL_FITTING_
  if (Verbosity() >= 1) _t_full_fitting->restart();
  if (_fitter->processTrack(gftrk_iter->second.get(), false) != 0)
  {
    if (Verbosity() >= 1) LogWarning("Track fitting failed\n");
    return -1;
  }
  if (Verbosity() >= 1) _t_full_fitting->stop();
  if (Verbosity() >= 1) _t_output_io->restart();

  track.set_chisq(gftrk_iter->second->get_chi2());
  track.set_ndf(gftrk_iter->second->get_ndf());

  // Use fitted vertex
  TVector3 vertex_position(_vertex[vertex_id][0], _vertex[vertex_id][1], _vertex[vertex_id][2]);
  std::unique_ptr<genfit::MeasuredStateOnPlane> gf_state_vertex_ca = nullptr;
  try
  {

    gf_state_vertex_ca = std::unique_ptr<genfit::MeasuredStateOnPlane>(gftrk_iter->second->extrapolateToPoint(vertex_position));

  } catch (...) {

    if (Verbosity() >= 2) LogWarning("extrapolateToPoint failed!");

  }

  if (!gf_state_vertex_ca) return -1;

  auto mom = gf_state_vertex_ca->getMom();
  auto pos = gf_state_vertex_ca->getPos();
  auto cov = gf_state_vertex_ca->get6DCov();
  #else
  auto state = gftrk_iter->second->getGenFitTrack()->getStateSeed();
  TVector3 pos(state(0), state(1), state(2));
  TVector3 mom(state(3), state(4), state(5));
  #endif
  track->set_px(mom.Px());
  track->set_py(mom.Py());
  track->set_pz(mom.Pz());

  track->set_x(pos.X());
  track->set_y(pos.Y());
  track->set_z(pos.Z());

  for (TrkrDefs::cluskey cluster_key : gftrk_iter->second->get_cluster_keys())
  {
    if(Verbosity() > 10) std::cout << " track id: " << phtrk_iter->first <<  " adding clusterkey " << cluster_key << std::endl;
    _gftrk_hitkey_map.insert(std::make_pair(cluster_key, phtrk_iter->first));
    track->insert_cluster_key(cluster_key);
  }

  //Check track quality
  Int_t n_maps = 0;
  Int_t n_intt = 0;
  Int_t n_tpc = 0;
  Int_t n_outer = 0;

  for( auto iter = track->begin_cluster_keys(); iter != track->end_cluster_keys(); ++iter)
  {
    auto cluster_key = *iter;
    auto layer = TrkrDefs::getLayer(cluster_key);
    if( is_maps_layer( layer ) ) {

      ++n_maps;

    } else if( is_intt_layer( layer ) ) {

      ++n_intt;
      if (n_intt > 4)
      {
        std::cout << PHWHERE << " Can not have more than 4 INTT layers, quit!" << std::endl;
        exit(1);
      }

    } else if( is_tpc_layer( layer ) ) {

      ++n_tpc;

    } else if( is_outertracker_layer( layer ) ) {

      ++n_outer;

    }
  }

  // Add the cluster-track association to the association table for later use
  for( const auto& cluster_key : gftrk_iter->second->get_cluster_keys())
  { _assoc_container->SetClusterTrackAssoc(cluster_key, track->get_id()); }

  if (Verbosity() > 10)
  {
    std::cout
      << "Output propagated track " << track->get_id() << " vertex " << track->get_vertex_id()<< " quality = " << track->get_quality()
      << " clusters mvtx: " << n_maps << " intt: " << n_intt << " tpc:  " << n_tpc << " outer tracker: " << n_outer
      << std::endl;
    std::cout << "px = " << track->get_px() << " py = " << track->get_py()  << " pz = " << track->get_pz() << std::endl;
  }

  if (Verbosity() >= 1) _t_output_io->stop();

  return 0;
}

//______________________________________________________________
int PHGenFitTrkProp::SvtxTrackToPHGenFitTracks(const SvtxTrack* svtxtrack)
{
  // clean up working array for each event
  _PHGenFitTracks.clear();

  double time1 = 0;
  double time2 = 0;

  //FIXME used in analysis ntuple
  float kappa = 0;
  float d = 0;
  float phi = 0;
  float dzdl = 0;
  float z0 = 0;
  float nhit = 0;
  float ml = 0;
  float rec = 0;
  float dt = 0;

  if (Verbosity() > 1)
  {
    time1 = _t_translate1->get_accumulated_time();
    time2 = _t_translate1->get_accumulated_time();
    _t_translate1->restart();
  }

  TVector3 seed_pos(
      svtxtrack->get_x(),
      svtxtrack->get_y(),
      svtxtrack->get_z());

  TVector3 seed_mom(
      svtxtrack->get_px(),
      svtxtrack->get_py(),
      svtxtrack->get_pz());

  const float blowup_factor = 1.;

  TMatrixDSym seed_cov(6);
  for (int i = 0; i < 6; i++)
  {
    for (int j = 0; j < 6; j++)
    {
      seed_cov[i][j] = blowup_factor * svtxtrack->get_error(i, j);
    }
  }

  if(Verbosity() > 10) std::cout << PHWHERE << "Converting SvtxTrack to PHGenFit track: track ID " << svtxtrack->get_id()
       << " track z " << svtxtrack->get_z()
       << " vertex ID " << svtxtrack->get_vertex_id()
       << " vertex z " << _vertex[svtxtrack->get_vertex_id()][2]
       << std::endl;

  auto rep = new genfit::RKTrackRep(_primary_pid_guess);
  std::shared_ptr<PHGenFit::Track> track( new PHGenFit::Track(rep, seed_pos, seed_mom, seed_cov));

  std::multimap<float, TrkrDefs::cluskey> m_r_clusterID;
  for (auto hit_iter = svtxtrack->begin_cluster_keys(); hit_iter != svtxtrack->end_cluster_keys(); ++hit_iter)
  {

    auto clusterkey = *hit_iter;
    auto cluster = _cluster_map->findCluster(clusterkey);
    float r = std::sqrt(
      square( cluster->getPosition(0) ) +
      square( cluster->getPosition(1) ) );

    m_r_clusterID.insert(std::make_pair(r, clusterkey));
    if (Verbosity() >= 10)
    {
      std::cout << PHWHERE << " inserted r " << r << " clusterkey " << clusterkey
        << " layer: " << (float)TrkrDefs::getLayer(clusterkey)
        <<  std::endl;
    }
  }

  std::vector<PHGenFit::Measurement*> measurements;
  for (auto iter = m_r_clusterID.begin();  iter != m_r_clusterID.end();  ++iter)
  {

    auto cluster_key = iter->second;
    auto cluster = _cluster_map->findCluster(cluster_key);
    ml += TrkrDefs::getLayer(cluster_key);
    if (!cluster)
    {
      LogError("No cluster Found!\n");
      continue;
    }

    auto meas = TrkrClusterToPHGenFitMeasurement(cluster);
    if (meas) measurements.push_back(meas);
  }
  track->addMeasurements(measurements);

  if (Verbosity() > 1) _t_translate2->stop();
  if (Verbosity() > 1) _t_translate3->restart();

  if (_fitter->processTrack(track.get(), false) != 0)
  {
    if (Verbosity() >= 1) LogWarning("Seed fitting failed") << std::endl;
    if (Verbosity() > 1) _t_translate3->stop();
    if (Verbosity() > 1)
    {
      _t_translate1->stop();
      time2 = _t_translate1->get_accumulated_time();
    }
    dt = time2 - time1;
    if (_analyzing_mode == true)
    { _analyzing_ntuple->Fill(svtxtrack->get_pt(), kappa, d, phi, dzdl, z0, nhit, ml / nhit, rec, dt); }
    return -1;
  }

  int nhits = track->get_cluster_keys().size();
  float chi2 = track->get_chi2();
  float ndf = track->get_ndf();

  if (nhits > 0 && chi2 > 0 && ndf > 0)
  { _PHGenFitTracks.push_back( std::make_pair( PHGenFitTrkProp::TrackQuality(nhits, chi2, ndf, 0, nhits, 0, 0), track)); }

  if (Verbosity() > 1) _t_translate3->stop();
  if (Verbosity() > 1)
  {
    _t_translate1->stop();
    time2 = _t_translate1->get_accumulated_time();
  }
  dt = time2 - time1;
  rec = 1;
  if (_analyzing_mode == true)
    _analyzing_ntuple->Fill(svtxtrack->get_pt(), kappa, d, phi, dzdl, z0, nhit, rec, dt);

  return Fun4AllReturnCodes::EVENT_OK;
}

//______________________________________________________________
int PHGenFitTrkProp::TrackPropPatRec(
  const unsigned int ivert,
  MapPHGenFitTrack::iterator& track_iter,
  unsigned int init_layer, unsigned int end_layer,
  const bool use_fitted_state_once)
{
  #ifdef _DEBUG_
  std::cout
      << __LINE__
      << " TrackPropPatRec"
      << " : init_layer: " << init_layer
      << " : end_layer: " << end_layer
      << " : use_fitted_state_once: " << use_fitted_state_once
      << std::endl;
  #endif

  auto& track = track_iter->second;
  if(Verbosity() > 10) std::cout << std::endl << PHWHERE << " Entering TrackPropPatRec for track with vertex ID " << ivert << std::endl;
  if(ivert >= _vertex.size())
  {
    std::cout << PHWHERE << " WARNING: vertex ID out of range, something wrong, quit! " << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  int direction = end_layer >= init_layer ? 1 : -1;
  int first_extrapolate_base_TP_id = -1;

  bool use_fitted_state = use_fitted_state_once;
  float blowup_factor = use_fitted_state ? _blowup_factor : 1.;
  unsigned int layer_occupied[_nlayers_all];

  for(int i = 0;i<_nlayers_all;i++) layer_occupied[i] = 0;

  /*
  * Find the last layer of with TrackPoint (TP)
  * Asumming measuremnts are sorted by radius
  * and cluster keys are syncronized with the TP IDs
  */
  const auto clusterkeys = track->get_cluster_keys();
  for (unsigned int i = 0; i < clusterkeys.size(); ++i)
  {
    unsigned int layer = TrkrDefs::getLayer(clusterkeys[i]);
    layer_occupied[layer] = 1;
    if(layer == init_layer)
    {
      first_extrapolate_base_TP_id = i;
      break;
    }
  }

  if (first_extrapolate_base_TP_id < 0)
  {
    if (Verbosity() > 0) LogError("first_extrapolate_base_TP_id < 0");
    return -1;
  }

  int extrapolate_base_TP_id = first_extrapolate_base_TP_id;

  unsigned int consecutive_missing_layer = 0;
  for( unsigned int layer = init_layer + direction; layer != end_layer + direction; layer += direction)
  {
    // layer is unsigned int, check for >=0 is meaningless
    if (layer >= (unsigned int) _nlayers_all) break;
    if (layer_occupied[layer]) continue;

    /*!
     * if miss too many layers terminate track propagating
     */
    if (consecutive_missing_layer > _max_consecutive_missing_layer)
    {
      if (Verbosity() > 1)
      {
        LogWarning("consecutive_missing_layer > ") << _max_consecutive_missing_layer << std::endl;
      }
      if (track->get_cluster_keys().size() >= _min_good_track_hits)
        return 0;
      else
        return -1;
    }

    bool layer_updated = false;
    float layer_r = _radii_all[_layer_ilayer_map_all[layer]];

    if(Verbosity() > 10)
    {
      std::cout << "=========================" << std::endl;
      std::cout << __LINE__ << ": Event: " << _event << ": _PHGenFitTracks.size(): " << _PHGenFitTracks.size() << ": layer: " << layer << std::endl;
      std::cout << "=========================" << std::endl;
    }



    #ifdef _DEBUG_
    {
      unsigned int tempIdx = extrapolate_base_TP_id >= 0 ? extrapolate_base_TP_id : extrapolate_base_TP_id + track->get_cluster_keys().size();
      std::cout
        << __LINE__
        << " tempIdx: " << tempIdx
        << std::endl;
      // tempIdx is unsigned int, checking for >=0 is meaningless
      if (tempIdx < track->get_cluster_keys().size())
      {
        auto extrapolate_base_cluster_id = track->get_cluster_keys()[tempIdx];
        auto extrapolate_base_cluster = _cluster_map->findCluster(extrapolate_base_cluster_id);
        int from_layer = TrkrDefs::getLayer(extrapolate_base_cluster_id);
        std::cout
          << __LINE__
          << ": Target layer: { " << layer
          << ", " << layer_r
          << "} : From layer: { " << from_layer
          << ", " << std::sqrt( square(extrapolate_base_cluster->getX()) + square(extrapolate_base_cluster->getY()))
          << "} : ID: " << extrapolate_base_cluster_id
          << std::endl;
      }
    }
    #endif

    std::unique_ptr<genfit::MeasuredStateOnPlane> state;
    try
    {
      state = std::unique_ptr<genfit::MeasuredStateOnPlane>(track->extrapolateToCylinder(
        layer_r, TVector3(0, 0, 0),
        TVector3(0, 0, 1), extrapolate_base_TP_id, direction));
    } catch (...) {

      // if (Verbosity() > 1)
      LogWarning("Can not extrapolate to Cylinder") << " from " << extrapolate_base_TP_id << " to layer: " << layer << " (r= " << layer_r << ")" << std::endl;
      continue;

    }

    if (!state)
    {
      // if (Verbosity() > 1)
      LogWarning("Can not extrapolate to Cylinder") << " from " << extrapolate_base_TP_id << " to layer: " << layer << " (r= " << layer_r << ")" << std::endl;
      continue;
    }

    // The vertex position is used here to caculate theta and phi for the
    // extrapolated track position. These values of theta and phi are compared
    // with the cluster theta and phi positions in the map  "_layer_thetaID_phiID_clusterID"

    // Need to use the correct vertex position for this track as well as for the theta-phi map.

    auto pos = state->getPos();
    if(Verbosity() > 10) std::cout << PHWHERE << " pos z " << pos.Z() << " _vertex z " << _vertex[ivert][2] << std::endl;
    pos.SetXYZ(
      pos.X() - _vertex[ivert][0],
      pos.Y() - _vertex[ivert][1],
      pos.Z() - _vertex[ivert][2]);

    float phi_center = pos.Phi();
    float theta_center = pos.Theta();

    #ifdef _USE_CONSTANT_SEARCH_WIN_
    float phi_window = 25e-4;
    float theta_window = 25e-4;

    if (layer >= 3 && layer <= 6)
    {
      phi_window = 300e-4;
      theta_window = 0.2;
    }

    if (layer <= 2)
    {
      phi_window = 3000e-4;
      theta_window = 3000e-4;
    }
    #else
    TMatrixDSym cov = state->get6DCov();
    float phi_window = _search_wins_phi[layer] * std::sqrt(cov[0][0] + cov[1][1] + cov[0][1] + cov[1][0]) / pos.Perp();
    float theta_window = _search_wins_theta[layer] * std::sqrt(cov[2][2]) / pos.Perp();

    // compare to built-in windows
    if( is_maps_layer( layer ) )
    {

      if (phi_window > _max_search_win_phi_maps) phi_window = _max_search_win_phi_maps;
      if (phi_window < _min_search_win_phi_maps) phi_window = _min_search_win_phi_maps;
      if (theta_window > _max_search_win_theta_maps) theta_window = _max_search_win_theta_maps;
      if (theta_window < _min_search_win_theta_maps) theta_window = _min_search_win_theta_maps;

    } else if( is_intt_layer( layer ) ) {

      const auto layer_intt = layer - _firstlayer_intt;
      if (phi_window > _max_search_win_phi_intt[layer_intt]) phi_window = _max_search_win_phi_intt[layer_intt];
      if (phi_window < _min_search_win_phi_intt[layer_intt]) phi_window = _min_search_win_phi_intt[layer_intt];
      if (theta_window > _max_search_win_theta_intt[layer_intt]) theta_window = _max_search_win_theta_intt[layer_intt];
      if (theta_window < _min_search_win_theta_intt[layer_intt]) theta_window = _min_search_win_theta_intt[layer_intt];

    } else if( is_tpc_layer( layer ) ) {

      if (phi_window > _max_search_win_phi_tpc) phi_window = _max_search_win_phi_tpc;
      if (phi_window < _min_search_win_phi_tpc) phi_window = _min_search_win_phi_tpc;
      if (theta_window > _max_search_win_theta_tpc) theta_window = _max_search_win_theta_tpc;
      if (theta_window < _min_search_win_theta_tpc) theta_window = _min_search_win_theta_tpc;

    } else if( is_outertracker_layer( layer ) ) {

      if (phi_window > _max_search_win_phi_outer) phi_window = _max_search_win_phi_outer;
      if (phi_window < _min_search_win_phi_outer) phi_window = _min_search_win_phi_outer;
      if (theta_window > _max_search_win_theta_outer) theta_window = _max_search_win_theta_outer;
      if (theta_window < _min_search_win_theta_outer) theta_window = _min_search_win_theta_outer;
    }
    #endif

    #ifdef _DEBUG_
    std::cout << __LINE__ << ": ";
    printf("layer: %u: r: %f: phi: %f +- %f; theta: %f +- %f\n",
      layer, pos.Perp(),
      phi_center, phi_window,
      theta_center, theta_window);
    #endif

    if (Verbosity() >= 1) _t_search_clusters->restart();
    std::vector<TrkrDefs::cluskey> new_cluster_keys = SearchHitsNearBy(
      ivert, layer,
      theta_center, phi_center, theta_window, phi_window);
    if (Verbosity() >= 1) _t_search_clusters->stop();

    #ifdef _DEBUG_
    std::cout << __LINE__ << ": new_cluster_keys size: " << new_cluster_keys.size() << std::endl;
    #endif

    std::vector<PHGenFit::Measurement*> measurements;
    for (TrkrDefs::cluskey cluster_key : new_cluster_keys)
    {
      auto cluster = _cluster_map->findCluster(cluster_key);
      if (!cluster)
      {
        LogError("No cluster Found!\n");
        continue;
      }

      auto meas = TrkrClusterToPHGenFitMeasurement(cluster);
      if (meas) measurements.push_back(meas);
    }

    std::map<double, std::shared_ptr<PHGenFit::Track> > incr_chi2s_new_tracks;

    #ifdef _DEBUG_
    std::cout << __LINE__ << ": measurements.size(): " << measurements.size() << std::endl;
    #endif

    if (Verbosity() >= 1) _t_track_propagation->restart();
    track->updateOneMeasurementKalman(measurements, incr_chi2s_new_tracks, extrapolate_base_TP_id, direction, blowup_factor, use_fitted_state);
    use_fitted_state = false;
    blowup_factor = 1.;
    if (Verbosity() >= 1) _t_track_propagation->stop();

    #ifdef _DEBUG_
    std::cout << __LINE__ << ": incr_chi2s_new_tracks.size(): " << incr_chi2s_new_tracks.size() << std::endl;
    #endif

    PHGenFitTrkProp::TrackQuality tq(track_iter->first);

    // Update first track candidate
    if (incr_chi2s_new_tracks.size() > 0)
    {
      auto iter = incr_chi2s_new_tracks.begin();

      if (iter->first < _max_incr_chi2s[layer] && iter->first > 0)
      {
        #ifdef _DEBUG_
        std::cout
          << __LINE__
          << ": iPHGenFitTrack: " << iPHGenFitTrack << std::endl
          << ": First accepted IncrChi2: " << iter->first << std::endl
          << "; before update: " << track->get_cluster_keys().back()
          << std::endl;
        #endif

        track_iter->first.nhits = tq.nhits + 1;
        track_iter->first.chi2 = tq.chi2 + iter->first;
        track_iter->first.ndf = tq.ndf + 2;
        track_iter->first.nouter = tq.nouter + (is_outertracker_layer( layer ) ? 1 : 0);
        track_iter->first.ntpc = tq.ntpc + (is_tpc_layer( layer ) ? 1 : 0);
        track_iter->first.nintt = tq.nintt + (is_intt_layer( layer ) ? 1 : 0);
        track_iter->first.nmaps = tq.nmaps + (is_maps_layer( layer ) ? 1 : 0);

        track_iter->second = iter->second;

        consecutive_missing_layer = 0;
        layer_updated = true;
        extrapolate_base_TP_id = -1;

        #ifdef _DEBUG_
        std::cout
          << __LINE__
          << ": after update: " << track->get_cluster_keys().back()
          << std::endl;

        fout_chi2
          << _event << "\t"
          << iPHGenFitTrack << "\t"
          << layer << "\t "
          << iter->first
          << std::endl;
        #endif
      }
    }

    // Update other candidates
    if( incr_chi2s_new_tracks.size() > 1 && is_intt_layer( layer ) )
    {
      for (auto iter = (++incr_chi2s_new_tracks.begin()); iter != incr_chi2s_new_tracks.end(); ++iter)
      {

        if (!(iter->first < _max_splitting_chi2 && iter->first > 0))
        { break; }

        #ifdef _DEBUG_
        std::cout << __LINE__ << ": "
          << "Track Spliting with "
          << "IncrChi2: " << iter->first << std::endl;
        #endif

        // this is a new track, have to associate it with the vertex
        iter->second->set_vertex_id(ivert);
        _PHGenFitTracks.push_back(
          std::make_pair(
          PHGenFitTrkProp::TrackQuality(
          tq.nhits + 1,
          tq.chi2 + iter->first,
          tq.ndf + 2,
          tq.nouter + (is_outertracker_layer( layer ) ? 1 : 0),
          tq.ntpc + (is_tpc_layer( layer ) ? 1 : 0),
          tq.nintt + (is_intt_layer( layer ) ? 1 : 0),
          tq.nmaps + (is_maps_layer( layer ) ? 1 : 0) ),
          iter->second));
      }

      #ifdef _DEBUG_
      std::cout
        << __LINE__ << ": "
        << "_PHGenFitTracksSize: " << _PHGenFitTracks.size() << std::endl;
      std::cout << __LINE__ << ": " << track_iter->second->get_cluster_keys().back() << std::endl;
      #endif
    }

    #ifdef _DEBUG_
    std::cout<<__LINE__<<": updateOneMeasurementKalman:"<<std::endl;
    std::cout<<"iPHGenFitTrack: "<<iPHGenFitTrack
      <<", layer: "<<layer
      <<", #meas: "<<measurements.size()
      <<", #tracks: "<<incr_chi2s_new_tracks.size()
      <<std::endl;

    for (auto iter = incr_chi2s_new_tracks.begin(); iter != incr_chi2s_new_tracks.end(); iter++)
    { std::cout << __LINE__ << ": IncrChi2: " << iter->first << std::endl; }
    #endif

    if (!layer_updated) ++consecutive_missing_layer;
  }

  #ifdef _DEBUG_
  std::cout
    << __LINE__
    << ": cluster keys size:  " << track->get_cluster_keys().size()
    << std::endl;
  #endif

  //! Track succesfully propagated and return 0
  return 0;
}

//________________________________________________________________________________
PHGenFit::Measurement* PHGenFitTrkProp::TrkrClusterToPHGenFitMeasurement( const TrkrCluster* cluster)
{

  if (!cluster) return nullptr;

  TVector3 pos(cluster->getPosition(0), cluster->getPosition(1), cluster->getPosition(2));
  TVector3 n(cluster->getPosition(0), cluster->getPosition(1), 0);

  // get the trkrid
  const auto cluster_id = cluster->getClusKey();
  unsigned int trkrid = TrkrDefs::getTrkrId(cluster_id);
  int layer = TrkrDefs::getLayer(cluster_id);

  // todo: check if this is also applicable to outer detector
  if(trkrid == TrkrDefs::mvtxId)
  {

    int stave_index = MvtxDefs::getStaveId(cluster_id);
    int chip_index = MvtxDefs::getChipId(cluster_id);

    auto geom = dynamic_cast<CylinderGeom_Mvtx*>(_geom_container_maps->GetLayerGeom(layer));
    std::array<double,3> ladder_location;
    geom->find_sensor_center(stave_index, 0, 0, chip_index, &ladder_location[0]);
    n.SetXYZ(ladder_location[0], ladder_location[1], 0);
    n.RotateZ(geom->get_stave_phi_tilt());

  } else if(trkrid == TrkrDefs::inttId) {

    auto geom = dynamic_cast<CylinderGeomIntt*>(_geom_container_intt->GetLayerGeom(layer));
    std::array<double,3> hit_location;
    geom->find_segment_center(InttDefs::getLadderZId(cluster_id), InttDefs::getLadderPhiId(cluster_id), &hit_location[0]);
    n.SetXYZ(hit_location[0], hit_location[1], 0);
    n.RotateZ(geom->get_strip_phi_tilt());

  }

  auto meas = new PHGenFit::PlanarMeasurement(pos, n, cluster->getRPhiError(), cluster->getZError());
  meas->set_cluster_key(cluster->getClusKey());

  #ifdef _DEBUG_
  int layer_out = TrkrDefs::getLayer(cluster->getClusKey());
  std::cout
    << __LINE__
    << ": ID: " << cluster->getClusKey()
    << ": layer: " << layer_out
    << ": pos: {" << pos.X() << ", " << pos.Y() << ", " << pos.Z() << "}"
    << ": n: {" << n.X() << ", " << n.Y() << ", " << n.Z() << "}"
    << ": r*phi_error: " << cluster->getRPhiError()
    << ": z error: " << cluster->getZError()
    << std::endl;
  #endif

  return meas;
}

//_____________________________________________________________________________
int PHGenFitTrkProp::BuildLayerZPhiHitMap(unsigned int ivert)
{
  // make this map for each collision vertex
  std::multimap<unsigned int, TrkrDefs::cluskey> this_layer_thetaID_phiID_clusterID;

  auto clusrange = _cluster_map->getClusters();
  for(TrkrClusterContainer::ConstIterator clusiter = clusrange.first; clusiter != clusrange.second; ++clusiter)
  {
    auto cluster = clusiter->second;
    auto cluskey = clusiter->first;

    // This z-phi map relative to the collision vertex includes all clusters,
    // we do not know whch clusters belong to which vertices yet
    // make a map for every vertex?

    unsigned int layer = TrkrDefs::getLayer(cluskey);
    float x = cluster->getPosition(0) - _vertex[ivert][0];
    float y = cluster->getPosition(1) - _vertex[ivert][1];
    float z = cluster->getPosition(2) - _vertex[ivert][2];

    float phi = atan2(y, x);
    float r = std::sqrt( square(x) + square(y) );
    float theta = std::atan2(r, z);

    unsigned int idx = encode_cluster_index(layer, theta, phi);
    this_layer_thetaID_phiID_clusterID.insert(std::make_pair(idx, cluster->getClusKey()));
  }

  _layer_thetaID_phiID_clusterID.push_back(this_layer_thetaID_phiID_clusterID);

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________________
std::vector<TrkrDefs::cluskey> PHGenFitTrkProp::SearchHitsNearBy(
  const unsigned int ivert, const unsigned int layer,
  const float theta_center, const float phi_center, const float theta_window,
  const float phi_window)
{
  if(Verbosity() > 10) std::cout << "SearchHitsNearBy for ivert = " << ivert << std::endl;

  std::vector<TrkrDefs::cluskey> cluster_keys;

  const unsigned int max_phi_bin = 16383;  //2^14 - 1
  const unsigned int max_z_bin = 2047;     // 2^11 - 1

  float lower_phi = phi_center - phi_window + _half_max_phi;
  float upper_phi = phi_center + phi_window + _half_max_phi;
  float lower_z = theta_center - theta_window + _half_max_theta;
  float upper_z = theta_center + theta_window + _half_max_theta;

  unsigned int lower_phi_bin = (unsigned int) ((lower_phi) / _layer_thetaID_phiID_clusterID_phiSize);
  unsigned int upper_phi_bin = (unsigned int) ((upper_phi) / _layer_thetaID_phiID_clusterID_phiSize);
  unsigned int lower_z_bin = (unsigned int) ((lower_z) / _layer_thetaID_phiID_clusterID_zSize);
  unsigned int upper_z_bin = (unsigned int) ((upper_z) / _layer_thetaID_phiID_clusterID_zSize);

  if (lower_phi < 0) lower_phi_bin = 0;
  if (upper_phi_bin > max_phi_bin) upper_phi_bin = max_phi_bin;

  if (lower_z < 0) lower_z_bin = 0;
  if (upper_z_bin > max_z_bin) upper_z_bin = max_z_bin;

  for (unsigned int iz = lower_z_bin; iz <= upper_z_bin; ++iz)
  {
    for (unsigned int irphi = lower_phi_bin; irphi <= upper_phi_bin; ++irphi)
    {
      if (Verbosity() >= 2) _t_search_clusters_encoding->restart();
      unsigned int idx = encode_cluster_index(layer, iz, irphi);
      if (Verbosity() >= 2) _t_search_clusters_encoding->stop();

      if(Verbosity() > 10)
      {
        if(_layer_thetaID_phiID_clusterID[ivert].count(idx) > 0)
        {
          std::cout
            <<__LINE__<<": "
            <<"{ "
            <<layer <<", "
            <<iz <<", "
            <<irphi << "} =>"
            <<idx << ": size: "
            <<_layer_thetaID_phiID_clusterID[ivert].count(idx)
            << std::endl;
        }
      }

      if (Verbosity() >= 2) _t_search_clusters_map_iter->restart();
      for (auto iter = _layer_thetaID_phiID_clusterID[ivert].lower_bound(idx); iter != _layer_thetaID_phiID_clusterID[ivert].upper_bound(idx); ++iter)
      {
        if(Verbosity() > 10) std::cout << "      adding cluster with key " << iter->second << std::endl;
        cluster_keys.push_back(iter->second);
      }

      if (Verbosity() >= 2) _t_search_clusters_map_iter->stop();
    }
  }

  if(Verbosity() > 10)
  {
    std::cout
      << __LINE__ << ": "
      << "layer: " << layer
      << ", rphi: {" << lower_phi_bin
      << ", " << upper_phi_bin
      << "}, z: {" << lower_z_bin
      << ", " << upper_z_bin
      << "}, found #clusters: " << cluster_keys.size()
      << std::endl;
  }

  return cluster_keys;
}

//_______________________________________________________________________
unsigned int PHGenFitTrkProp::encode_cluster_index(const unsigned int layer, const float theta, const float phi)
{
  unsigned int idx = UINT_MAX;

  if (layer >= 128)
  {
    LogError("layer >= 128\n");
    return idx;
  }

  if ((theta + _half_max_theta) < 0)
  {
    LogError("(theta + _half_max_theta) < 0 \n");
    return idx;
  }
  unsigned int itheta = (theta + _half_max_theta) / _layer_thetaID_phiID_clusterID_zSize;
  if ((phi + _half_max_phi) < 0)
  {
    LogError("(rphi + _half_max_rphi) < 0 \n");
    return idx;
  }
  unsigned int irphi = (phi + _half_max_phi) / _layer_thetaID_phiID_clusterID_phiSize;
  return encode_cluster_index(layer, itheta, irphi);
}

//_______________________________________________________________________________
unsigned int PHGenFitTrkProp::encode_cluster_index(const unsigned int layer, const unsigned int iz, const unsigned int irphi)
{
  if (layer >= 128)
  {
    LogError("layer >= 128: ") << layer << std::endl;
    return UINT_MAX;
  }

  if (iz >= 2048)
  {
    LogError("iz >= 2048: ") << iz << std::endl;
    return UINT_MAX;
  }

  if (irphi >= 16384)
  {
    LogError("irphi >= 16384: ") << irphi << std::endl;
    return UINT_MAX;
  }

  unsigned int index = 0;

  index |= (layer << 25);
  index |= (iz << 14);
  index |= (irphi);

  return index;
}

//_______________________________________________________________________________
int PHGenFitTrkProp::InitializeGeometry(PHCompositeNode* topNode)
{
  auto outergeos = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_OuterTracker");
  auto cellgeos = findNode::getClass<PHG4CylinderCellGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  auto laddergeos = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_INTT");
  auto mapsladdergeos = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MVTX");

  std::map<float, int> radius_layer_map;

  _radii_all.assign(_nlayers_all, 0.0);
  _layer_ilayer_map_all.clear();

  // outer tracker
  if (outergeos)
  {
    const auto range = outergeos->get_begin_end();
    for (auto layeriter = range.first; layeriter != range.second; ++layeriter)
    {
      std::cout
        << "PHGenFitTrkProp::InitializeGeometry -"
        << " outer-tracker adding layer: " << layeriter->second->get_layer()
        << " radius: " << layeriter->second->get_radius() << std::endl;
      radius_layer_map.insert( std::make_pair(layeriter->second->get_radius(), layeriter->second->get_layer()));
    }
  }

  // tpc
  if (cellgeos)
  {
    const auto range = cellgeos->get_begin_end();
    for (auto layeriter = range.first; layeriter != range.second; ++layeriter)
    { radius_layer_map.insert( std::make_pair(layeriter->second->get_radius(), layeriter->second->get_layer())); }
  }

  // intt
  if (laddergeos)
  {
    const auto range = laddergeos->get_begin_end();
    for (auto layeriter = range.first; layeriter != range.second; ++layeriter)
    { radius_layer_map.insert( std::make_pair(layeriter->second->get_radius(), layeriter->second->get_layer())); }
  }

  // maps
  if (mapsladdergeos)
  {
    const auto range = mapsladdergeos->get_begin_end();
    for (auto layeriter = range.first; layeriter != range.second; ++layeriter)
    {
      std::cout
        << "PHGenFitTrkProp::InitializeGeometry -"
        << " mvtx adding layer: " << layeriter->second->get_layer()
        << " radius: " << layeriter->second->get_radius() << std::endl;
      radius_layer_map.insert( std::make_pair(layeriter->second->get_radius(), layeriter->second->get_layer()));

    }
  }

  // map layer id to index
  for( auto iter = radius_layer_map.begin(); iter != radius_layer_map.end(); ++iter)
  { _layer_ilayer_map_all.insert( std::make_pair(iter->second, _layer_ilayer_map_all.size())); }

  // now we extract the information from all geometries

  // outer tracker
  if (outergeos)
  {
    const auto range = outergeos->get_begin_end();
    for( auto iter = range.first; iter != range.second; iter++)
    {
      auto geo = iter->second;
      _radii_all[_layer_ilayer_map_all[geo->get_layer()]] = geo->get_radius() + 0.5 * geo->get_thickness();
    }
  }

  // tpc
  if (cellgeos)
  {
    const auto range = cellgeos->get_begin_end();
    for( auto iter = range.first; iter != range.second; iter++)
    {
      auto geo = iter->second;
      _radii_all[_layer_ilayer_map_all[geo->get_layer()]] = geo->get_radius() + 0.5 * geo->get_thickness();
    }
  }

  // intt
  if (laddergeos)
  {
    const auto range = laddergeos->get_begin_end();
    for( auto iter = range.first; iter != range.second; iter++)
    {
      auto geo = iter->second;
      _radii_all[_layer_ilayer_map_all[geo->get_layer()]] = geo->get_radius() + 0.5 * geo->get_thickness();
    }
  }

  // maps
  if (mapsladdergeos)
  {
    const auto range = mapsladdergeos->get_begin_end();
    for( auto iter = range.first; iter != range.second; iter++)
    {
      auto geo = iter->second;
      _radii_all[_layer_ilayer_map_all[geo->get_layer()]] = geo->get_radius();
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
