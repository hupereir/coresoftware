#include "TrackingEvaluator_hp.h"

#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/PHTimer.h>
#include <iostream>

ClassImp(ClusterStruct)

//_____________________________________________________________________
namespace
{

  // square
  template<class T> T square( T x ) { return x*x; }

  // radius
  template<class T> T get_r( T x, T y ) { return std::sqrt( square(x) + square(y) ); }

  // phi
  template<class T> T get_phi( T x, T y ) { return std::atan2( y, x ); }

  // create cluster struct from svx cluster
  ClusterStruct create_cluster( TrkrDefs::cluskey key, TrkrCluster* cluster )
  {
    ClusterStruct clusterStruct;
    clusterStruct._layer = TrkrDefs::getLayer(key);
    clusterStruct._x = cluster->getX();
    clusterStruct._y = cluster->getY();
    clusterStruct._z = cluster->getZ();
    clusterStruct._r = get_r( cluster->getX(), cluster->getY() );
    clusterStruct._phi = get_phi( cluster->getX(), cluster->getY() );
    return clusterStruct;
  }

  // add track information
  void add_trk_information( ClusterStruct& cluster, SvtxTrackState* state )
  {
    cluster._trk_x = state->get_x();
    cluster._trk_y = state->get_y();
    cluster._trk_z = state->get_z();
    cluster._trk_r = get_r( state->get_x(), state->get_y() );
    cluster._trk_phi = get_phi( state->get_x(), state->get_y() );
  }

}

//_____________________________________________________________________
TrackingEvaluator_hp::TrackingEvaluator_hp( const std::string& name, const std::string& filename ):
  SubsysReco( name),
  _filename( filename )
{
  std::cout << "TrackingEvaluator_hp::TrackingEvaluator_hp." << std::endl;
}

//_____________________________________________________________________
int TrackingEvaluator_hp::Init(PHCompositeNode* )
{
  std::cout << "TrackingEvaluator_hp::Init." << std::endl;

  _tfile.reset( new TFile( _filename.c_str(), "RECREATE" ) );

  const Int_t kSplitlevel = 98;
  const Int_t kBufsize = 32000;
  _clusterArray = new TClonesArray( "ClusterStruct" );
  _clusterArray->SetName( "ClusterArray" );
  _clusterArray->SetOwner( kTRUE );

  _tree = new TTree( "Tree", "Clusters" );
  _tree->Branch( "clusterArray", "TClonesArray", &_clusterArray, kBufsize, kSplitlevel );

  // initialize timer
  _timer.reset( new PHTimer("_tracking_evaluator_hp_timer") );
  _timer->restart();

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int TrackingEvaluator_hp::InitRun(PHCompositeNode* )
{
  std::cout << "TrackingEvaluator_hp::InitRun." << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int TrackingEvaluator_hp::process_event(PHCompositeNode* topNode)
{
  // print event number
  if( _ievent % 100 == 0 )
  { std::cout << "TrackingEvaluator_hp::process_event - Event = " << _ievent << std::endl; }
  ++_ievent;

  // load nodes
  auto res =  load_nodes(topNode);
  if( res != Fun4AllReturnCodes::EVENT_OK ) return res;

  // evaluate_clusters();
  evaluate_tracks();

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int TrackingEvaluator_hp::End(PHCompositeNode* )
{
  std::cout << "TrackingEvaluator_hp::End." << std::endl;

  // write tree to tfile and close
  if( _tfile )
  {
    _tfile->cd();
    if( _tree ) _tree->Write();
    _tfile->Close();
  }

  // print timer information
  _timer->stop();
  std::cout
    << "TrackingEvaluator_hp::End -"
    << " time per event:" << _timer->get_accumulated_time()/(1000.*_ievent) << " sec"
    << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int TrackingEvaluator_hp::load_nodes( PHCompositeNode* topNode )
{

  // get necessary nodes
  _trackMap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if( !_trackMap ) return Fun4AllReturnCodes::ABORTEVENT;

  // cluster map
  _clusterMap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if( !_clusterMap ) return Fun4AllReturnCodes::ABORTEVENT;

  return Fun4AllReturnCodes::EVENT_OK;

}

//_____________________________________________________________________
void TrackingEvaluator_hp::evaluate_clusters()
{

  if( !( _clusterMap && _clusterArray ) ) return;

  // clear array
  _clusterArray->Clear();
  _clusterCount = 0;

  auto range = _clusterMap->getClusters();
  for( auto clusterIter = range.first; clusterIter != range.second; ++clusterIter )
  {

    const auto& key = clusterIter->first;
    const auto& cluster = clusterIter->second;

    // create cluster structure and add in array
    auto clusterStruct = create_cluster( key, cluster );
    new((*_clusterArray)[_clusterCount++]) ClusterStruct( clusterStruct );

  }

  // fill tree
  if( _tree ) _tree->Fill();

}

//_____________________________________________________________________
void TrackingEvaluator_hp::evaluate_tracks()
{
  if( !( _trackMap && _clusterMap && _clusterArray ) ) return;

  // clear array
  _clusterArray->Clear();
  _clusterCount = 0;

  for( auto trackIter = _trackMap->begin(); trackIter != _trackMap->end(); ++trackIter )
  {

    auto track = trackIter->second;
    // print_track( track );

    // loop over clusters
    auto stateIter = track->begin_states();
    for( auto keyIter = track->begin_cluster_keys(); keyIter != track->end_cluster_keys(); ++keyIter )
    {
      auto cluster = _clusterMap->findCluster( *keyIter );
      if( !cluster )
      {
        std::cout << "TrackingEvaluator_hp::print_track - unable to find cluster for key " << *keyIter << std::endl;
        continue;
      }

      // create new cluster struct
      // and add in array
      auto clusterStruct = create_cluster( *keyIter, cluster );
      const auto radius( clusterStruct._r );

      // find first track state element whose radius is larger than current
      using StatePair = SvtxTrack::StateMap::value_type;
      stateIter = std::find_if( stateIter, track->end_states(), [radius](const StatePair& pair )
         { return get_r( pair.second->get_x(), pair.second->get_y() ) > radius; } );

      // copy and compare distance also to previous state to get the closest
      auto iter = stateIter;
      if( iter == track->end_states() ) --iter;
      else if( iter !=  track->begin_states() )
      {

        // get previous
        auto previous = iter; previous--;

        // compare radius difference
        const auto delta_r =  std::abs( radius - get_r(iter->second->get_x(), iter->second->get_y() ) );
        const auto delta_r_prev =  std::abs( radius - get_r(previous->second->get_x(), previous->second->get_y() ) );
        if( delta_r_prev < delta_r ) iter = previous;
      }

      // store track state in cluster struct
      add_trk_information( clusterStruct, iter->second );

      // add to array
      new((*_clusterArray)[_clusterCount++]) ClusterStruct( clusterStruct );

    }

  }

  // fill tree
  if( _tree ) _tree->Fill();

}

//_____________________________________________________________________
void TrackingEvaluator_hp::print_track(SvtxTrack* track)
{

  std::cout << "TrackingEvaluator_hp::print_track - track: " << track << std::endl;
  if( !track ) return;

  // loop over cluster keys
  if( _clusterMap )
  {
    for( auto keyIter = track->begin_cluster_keys(); keyIter != track->end_cluster_keys(); ++keyIter )
    {
      auto cluster = _clusterMap->findCluster( *keyIter );
      if( !cluster )
      {
        std::cout << "TrackingEvaluator_hp::print_track - unable to find cluster for key " << *keyIter << std::endl;
        continue;
      }

      std::cout
        << "TrackingEvaluator_hp::print_track -"
        << " cluster layer: "  << (int)TrkrDefs::getLayer(*keyIter)
        << " position: (" << cluster->getX() << ", " << cluster->getY() << ", " << cluster->getZ() << ")"
        << " polar: (" << get_r( cluster->getX(), cluster->getY() ) << ", " << get_phi( cluster->getX(), cluster->getY() ) << "," << cluster->getZ() << ")"
        << std::endl;

    }

  }

  // loop over track states
  for( auto stateIter = track->begin_states(); stateIter != track->end_states(); ++ stateIter )
  {
    auto state = stateIter->second;
    if( !state ) return;

    std::cout
        << "TrackingEvaluator_hp::print_track -"
        << " state pathLength: " << stateIter->first
        << " position: (" << state->get_x() << ", " << state->get_y() << ", " << state->get_z() << ")"
        << " polar: (" << get_r( state->get_x(), state->get_y() ) << ", " << get_phi( state->get_x(), state->get_y() ) << "," << state->get_z() << ")"
        << std::endl;
  }

  std::cout << "TrackingEvaluator_hp::print_track - done." << std::endl;

}
