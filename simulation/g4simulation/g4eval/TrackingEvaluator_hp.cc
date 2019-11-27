#include "TrackingEvaluator_hp.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHTimer.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrHitTruthAssoc.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>

#include <iostream>
#include <algorithm>

ClassImp(ClusterStruct)
ClassImp(TrackStruct)
ClassImp(Container)

//_____________________________________________________________________
namespace
{

  /// square
  template<class T> T square( T x ) { return x*x; }

  /// radius
  template<class T> T get_r( T x, T y ) { return std::sqrt( square(x) + square(y) ); }

  /// radius
  float get_r( PHG4Hit* hit, int i )
  {  return get_r( hit->get_x(i), hit->get_y(i) ); }

  /// phi
  template<class T> T get_phi( T x, T y ) { return std::atan2( y, x ); }

  /// phi
  float get_phi( PHG4Hit* hit, int i )
  {  return get_phi( hit->get_x(i), hit->get_y(i) ); }

  /// calculate the average of member function called on all members in collection
  template< float (PHG4Hit::*accessor)(int) const>
  float interpolate( std::set<PHG4Hit*> hits, float rextrap )
  {
    // calculate all terms needed for the interpolation
    // need to use double everywhere here due to numerical divergences
    double s1 = 0;
    double sr = 0;
    double sr2 = 0;
    double sphi = 0;
    double srphi = 0;

    for( const auto& hit:hits )
    {
      const double r0 = get_r( hit, 0 );
      const double r1 = get_r( hit, 1 );

      const double phi0 = (hit->*accessor)(0);
      const double phi1 = (hit->*accessor)(1);

      s1 += 2;
      sr += r0 + r1;
      sr2 += square(r0) + square(r1);
      sphi += phi0 + phi1;
      srphi += r0*phi0 + r1*phi1;
    }

    const auto alpha = (s1*srphi - sr*sphi);
    const auto beta = (sr2*sphi - sr*srphi);
    const auto denom = (s1*sr2 - square(sr));

    return ( alpha*rextrap + beta )/denom;
  }

  /// create cluster struct from svx cluster
  TrackStruct create_track( SvtxTrack* track )
  {
    TrackStruct trackStruct;

    trackStruct._charge = track->get_charge();
    trackStruct._x = track->get_x();
    trackStruct._y = track->get_y();
    trackStruct._z = track->get_z();
    trackStruct._r = get_r( trackStruct._x, trackStruct._y );
    trackStruct._phi = get_phi( trackStruct._x, trackStruct._y );

    trackStruct._px = track->get_px();
    trackStruct._py = track->get_py();
    trackStruct._pz = track->get_pz();
    trackStruct._pt = std::sqrt( square(trackStruct._px) + square( trackStruct._py ) );
    trackStruct._p = std::sqrt( square(trackStruct._px) + square( trackStruct._py ) + square( trackStruct._pz ) );

    const auto p = std::sqrt( square( trackStruct._pt ) + square( trackStruct._pz ) );
    trackStruct._eta = std::log( (p+trackStruct._pz)/(p-trackStruct._pz) )/2;

    return trackStruct;
  }

  /// create cluster struct from svx cluster
  ClusterStruct create_cluster( TrkrDefs::cluskey key, TrkrCluster* cluster )
  {
    ClusterStruct clusterStruct;
    clusterStruct._layer = TrkrDefs::getLayer(key);
    clusterStruct._x = cluster->getX();
    clusterStruct._y = cluster->getY();
    clusterStruct._z = cluster->getZ();
    clusterStruct._r = get_r( clusterStruct._x, clusterStruct._y );
    clusterStruct._phi = get_phi( clusterStruct._x, clusterStruct._y );
    return clusterStruct;
  }

  /// add track information
  void add_trk_information( ClusterStruct& cluster, SvtxTrackState* state )
  {

    // need to extrapolate to the right r
    const auto trk_r = get_r( state->get_x(), state->get_y() );
    const auto dr = cluster._r - trk_r;
    const auto trk_drdt = get_r( state->get_px(), state->get_py() );
    const auto trk_dxdr = state->get_px()/trk_drdt;
    const auto trk_dydr = state->get_py()/trk_drdt;
    const auto trk_dzdr = state->get_pz()/trk_drdt;

    // store
    cluster._trk_x = state->get_x() + dr*trk_dxdr;
    cluster._trk_y = state->get_y() + dr*trk_dydr;
    cluster._trk_z = state->get_z() + dr*trk_dzdr;
    cluster._trk_r = get_r( cluster._trk_x, cluster._trk_y );
    cluster._trk_phi = get_phi( cluster._trk_x, cluster._trk_y );

  }

  /// number of hits associated to cluster
  unsigned int cluster_size( TrkrDefs::cluskey key, TrkrClusterHitAssoc* cluster_hit_map )
  {
    if( cluster_hit_map )
    {
      const auto range = cluster_hit_map->getHits(key);
      return std::distance( range.first, range.second );
    } else {
      return 0;
    }
  }

  // add truth information
  void add_truth_information( ClusterStruct& cluster, std::set<PHG4Hit*> hits )
  {
    const auto rextrap = get_r( cluster._x, cluster._y );

    cluster._truth_size = hits.size();
    cluster._truth_x = interpolate<&PHG4Hit::get_x>( hits, rextrap );
    cluster._truth_y = interpolate<&PHG4Hit::get_y>( hits, rextrap );
    cluster._truth_z = interpolate<&PHG4Hit::get_z>( hits, rextrap );
    cluster._truth_r = get_r( cluster._truth_x, cluster._truth_y );
    cluster._truth_phi = get_phi( cluster._truth_x, cluster._truth_y );
  }

  // print to stream
  std::ostream& operator << (std::ostream& out, const ClusterStruct& cluster )
  {
    out << "ClusterStruct" << std::endl;
    out << "  cluster: (" << cluster._x << "," << cluster._y << "," << cluster._z << ")" << std::endl;
    out << "  track:   (" << cluster._trk_x << "," << cluster._trk_y << "," << cluster._trk_z << ")" << std::endl;
    out << "  truth:   (" << cluster._truth_x << "," << cluster._truth_y << "," << cluster._truth_z << ")" << std::endl;
    return out;
  }

}

//_____________________________________________________________________
Container::Container()
{
  _clusters = new TClonesArray( "ClusterStruct" );
  _clusters->SetName( "ClusterArray" );
  _clusters->SetOwner( kTRUE );

  _tracks = new TClonesArray( "TrackStruct" );
  _tracks->SetName( "TrackArray" );
  _tracks->SetOwner( kTRUE );
}

//_____________________________________________________________________
Container::~Container()
{
  delete _clusters;
  delete _tracks;
}

//_____________________________________________________________________
void Container::Reset()
{
  _clusters->Clear();
  _tracks->Clear();
}

//_____________________________________________________________________
TrackingEvaluator_hp::TrackingEvaluator_hp( const std::string& name ):
  SubsysReco( name)
{
  std::cout << "TrackingEvaluator_hp::TrackingEvaluator_hp." << std::endl;
}

//_____________________________________________________________________
int TrackingEvaluator_hp::Init(PHCompositeNode* topNode )
{

  std::cout << "TrackingEvaluator_hp::Init." << std::endl;

  // find DST node
  PHNodeIterator iter(topNode);
  auto dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << "TrackingEvaluator_hp::Init - DST Node missing" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // get EVAL node
  iter = PHNodeIterator(dstNode);
  auto evalNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "EVAL"));
  if( !evalNode )
  {
    // create
    std::cout << "TrackingEvaluator_hp::Init - EVAL node missing - creating" << std::endl;
    evalNode = new PHCompositeNode( "EVAL" );
    dstNode->addNode(evalNode);
  }

  _container = new Container;
  auto newNode = new PHIODataNode<PHObject>( _container, "Container","PHObject");
  evalNode->addNode(newNode);

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
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int TrackingEvaluator_hp::load_nodes( PHCompositeNode* topNode )
{

  // get necessary nodes
  _track_map = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");

  // cluster map
  _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");

  // cluster hit association map
  _cluster_hit_map = findNode::getClass<TrkrClusterHitAssoc>(topNode, "TRKR_CLUSTERHITASSOC");

  // cluster hit association map
  _hit_truth_map = findNode::getClass<TrkrHitTruthAssoc>(topNode,"TRKR_HITTRUTHASSOC");

  // local container
  _container = findNode::getClass<Container>(topNode, "Container");

  // g4hits
  _g4hits_tpc = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_TPC");
  _g4hits_intt = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_INTT");
  _g4hits_mvtx = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_MVTX");
  _g4hits_outertracker = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_OuterTracker");

  return Fun4AllReturnCodes::EVENT_OK;

}

//_____________________________________________________________________
void TrackingEvaluator_hp::evaluate_clusters()
{

  if( !( _cluster_map && _container ) ) return;

  // clear array
  _container->clusters()->Clear();
  _cluster_count = 0;

  auto range = _cluster_map->getClusters();
  for( auto clusterIter = range.first; clusterIter != range.second; ++clusterIter )
  {

    const auto& key = clusterIter->first;
    const auto& cluster = clusterIter->second;

    // create cluster structure and add in array
    auto clusterStruct = create_cluster( key, cluster );
    clusterStruct._size = cluster_size( key, _cluster_hit_map );
    new((*_container->clusters())[_cluster_count++]) ClusterStruct( clusterStruct );

  }

}

//_____________________________________________________________________
void TrackingEvaluator_hp::evaluate_tracks()
{
  if( !( _track_map && _cluster_map && _container ) ) return;

  // clear array
  _container->clusters()->Clear();
  _cluster_count = 0;

  _container->tracks()->Clear();
  _track_count = 0;

  for( auto trackIter = _track_map->begin(); trackIter != _track_map->end(); ++trackIter )
  {

    auto track = trackIter->second;
    auto trackStruct = create_track( track );

    // add to array
    new((*_container->tracks())[_track_count++]) TrackStruct( trackStruct );

    // loop over clusters
    auto state_iter = track->begin_states();
    for( auto key_iter = track->begin_cluster_keys(); key_iter != track->end_cluster_keys(); ++key_iter )
    {

      const auto& cluster_key = *key_iter;
      auto cluster = _cluster_map->findCluster( cluster_key );
      if( !cluster )
      {
        std::cout << "TrackingEvaluator_hp::evaluate_tracks - unable to find cluster for key " << cluster_key << std::endl;
        continue;
      }

      // create new cluster struct
      auto clusterStruct = create_cluster( cluster_key, cluster );
      clusterStruct._size = cluster_size( cluster_key, _cluster_hit_map );

      const auto radius( clusterStruct._r );

      // find track state that is the closest to cluster
      /* this assumes that both clusters and states are sorted along r */
      float dr_min = -1;
      for( auto iter = state_iter; iter != track->end_states(); ++iter )
      {
        const auto dr = std::abs( radius - get_r( iter->second->get_x(), iter->second->get_y() ) );
        if( dr_min < 0 || dr < dr_min )
        {
          state_iter = iter;
          dr_min = dr;
        } else break;
      }

      // store track state in cluster struct
      add_trk_information( clusterStruct, state_iter->second );

      // truth information
      add_truth_information( clusterStruct, find_g4hits( cluster_key ) );

      // add to array
      new((*_container->clusters())[_cluster_count++]) ClusterStruct( clusterStruct );

    }

  }

}

//_____________________________________________________________________
void TrackingEvaluator_hp::print_track(SvtxTrack* track) const
{

  if( !track ) return;

  // print track position and momentum
  std::cout << "TrackingEvaluator_hp::print_track - id: " << track->get_id() << std::endl;
  std::cout << "TrackingEvaluator_hp::print_track - position: (" << track->get_x() << ", " << track->get_y() << ", " << track->get_z() << ")" << std::endl;
  std::cout << "TrackingEvaluator_hp::print_track - momentum: (" << track->get_px() << ", " << track->get_py() << ", " << track->get_pz() << ")" << std::endl;
  std::cout << "TrackingEvaluator_hp::print_track - clusters: " << track->size_cluster_keys() << ", states: " << track->size_states() << std::endl;

  // loop over cluster keys
  if( false && _cluster_map )
  {
    for( auto key_iter = track->begin_cluster_keys(); key_iter != track->end_cluster_keys(); ++key_iter )
    {

      const auto& cluster_key = *key_iter;
      auto cluster = _cluster_map->findCluster( cluster_key );
      if( !cluster )
      {
        std::cout << "TrackingEvaluator_hp::print_track - unable to find cluster for key " << cluster_key << std::endl;
        continue;
      }

      std::cout
        << "TrackingEvaluator_hp::print_track -"
        << " cluster layer: "  << (int)TrkrDefs::getLayer(cluster_key)
        << " position: (" << cluster->getX() << ", " << cluster->getY() << ", " << cluster->getZ() << ")"
        << " polar: (" << get_r( cluster->getX(), cluster->getY() ) << ", " << get_phi( cluster->getX(), cluster->getY() ) << "," << cluster->getZ() << ")"
        << std::endl;

    }
  }

  // loop over track states
  if( false )
  {
    for( auto state_iter = track->begin_states(); state_iter != track->end_states(); ++ state_iter )
    {
      auto state = state_iter->second;
      if( !state ) return;

      std::cout
        << "TrackingEvaluator_hp::print_track -"
        << " state pathLength: " << state_iter->first
        << " position: (" << state->get_x() << ", " << state->get_y() << ", " << state->get_z() << ")"
        << " polar: (" << get_r( state->get_x(), state->get_y() ) << ", " << get_phi( state->get_x(), state->get_y() ) << "," << state->get_z() << ")"
        << std::endl;
    }
  }

  std::cout << std::endl;

}

//_____________________________________________________________________
void TrackingEvaluator_hp::print_cluster( TrkrDefs::cluskey key, TrkrCluster* cluster ) const
{

     std::cout
        << "TrackingEvaluator_hp::print_cluster -"
        << " layer: " << (int)TrkrDefs::getLayer(key)
        << " position: (" << cluster->getX() << "," << cluster->getY() << "," << cluster->getZ() << ")"
        << " polar: (" << get_r( cluster->getX(), cluster->getY()) << "," << get_phi( cluster->getX(), cluster->getY()) << "," << cluster->getZ() << ")"
        << std::endl;

      // get associated g4 hist
      auto g4hits = find_g4hits( key );
      for( const auto& g4hit:g4hits )
      {

        std::cout
          << "TrackingEvaluator_hp::print_cluster -"
          << " in: (" << g4hit->get_x(0) << "," << g4hit->get_y(0) << "," << g4hit->get_z(0) << ")"
          << " out: (" << g4hit->get_x(1) << "," << g4hit->get_y(1) << "," << g4hit->get_z(1) << ")"
          << " polar in: (" << get_r( g4hit->get_x(0), g4hit->get_y(0) ) << "," << get_phi( g4hit->get_x(0), g4hit->get_y(0) ) << "," << g4hit->get_z(0) << ")"
          << " polar out: (" << get_r( g4hit->get_x(1), g4hit->get_y(1) ) << "," << get_phi( g4hit->get_x(1), g4hit->get_y(1) ) << "," << g4hit->get_z(1) << ")"
          << std::endl;
      }

      // interpolate g4hits positions at the same radius as the cluster to get resolution
      const auto rextrap = get_r( cluster->getX(), cluster->getY());
      const auto xextrap = interpolate<&PHG4Hit::get_x>( g4hits, rextrap );
      const auto yextrap = interpolate<&PHG4Hit::get_y>( g4hits, rextrap );
      const auto zextrap = interpolate<&PHG4Hit::get_z>( g4hits, rextrap );

      // print interpolation
      std::cout
        << "TrackingEvaluator_hp::print_cluster -"
        << " interpolation: (" << xextrap << "," << yextrap << "," << zextrap << ")"
        << " polar: (" << get_r( xextrap, yextrap ) << "," << get_phi( xextrap, yextrap ) << "," << zextrap << ")"
        << std::endl;

}

//_____________________________________________________________________
TrackingEvaluator_hp::G4HitSet TrackingEvaluator_hp::find_g4hits( TrkrDefs::cluskey cluster_key ) const
{

  G4HitSet out;

  // check maps
  if( !( _cluster_hit_map && _hit_truth_map ) ) return out;

  // find hitset associated to cluster
  const auto hitset_key = TrkrDefs::getHitSetKeyFromClusKey(cluster_key);

  // loop over hits associated to clusters
  const auto range = _cluster_hit_map->getHits(cluster_key);
  for( auto iter = range.first; iter != range.second; ++iter )
  {

    // hit key
    const auto& hit_key = iter->second;

    // store hits to g4hit associations
    TrkrHitTruthAssoc::MMap g4hit_map;
    _hit_truth_map->getG4Hits( hitset_key, hit_key, g4hit_map );

    // find corresponding g4 hist
    for( auto truth_iter = g4hit_map.begin(); truth_iter != g4hit_map.end(); ++truth_iter )
    {

      const auto g4hit_key = truth_iter->second.second;
      PHG4Hit* g4hit = nullptr;

      switch( TrkrDefs::getTrkrId( hitset_key ) )
      {
        case TrkrDefs::tpcId:
        if( _g4hits_tpc ) g4hit = _g4hits_tpc->findHit( g4hit_key );
        break;

        case TrkrDefs::inttId:
        if( _g4hits_intt ) g4hit = _g4hits_intt->findHit( g4hit_key );
        break;

        case TrkrDefs::outertrackerId:
        if( _g4hits_outertracker ) g4hit = _g4hits_outertracker->findHit( g4hit_key );
        break;

        case TrkrDefs::mvtxId:
        if( _g4hits_mvtx ) g4hit = _g4hits_mvtx->findHit( g4hit_key );
        break;

        default: break;
      }

      if( g4hit ) out.insert( g4hit );
      else std::cout << "TrackingEvaluator_hp::find_g4hits - g4hit not found " << g4hit_key << std::endl;

    }
  }

  return out;

}
