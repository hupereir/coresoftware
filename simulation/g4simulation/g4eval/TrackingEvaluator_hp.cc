#include "TrackingEvaluator_hp.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <intt/InttDefs.h>
#include <micromegas/CylinderGeomMicromegas.h>
#include <micromegas/MicromegasDefs.h>
#include <mvtx/MvtxDefs.h>
#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>
#include <tpc/TpcDefs.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrHitTruthAssoc.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>

#include <TVector3.h>

#include <algorithm>
#include <bitset>
#include <cassert>
#include <iostream>
#include <numeric>

//_____________________________________________________________________
namespace
{

  //! range adaptor to be able to use range-based for loop
  template<class T> class range_adaptor
  {
    public:
    range_adaptor( const T& range ):m_range(range){}
    const typename T::first_type& begin() {return m_range.first;}
    const typename T::second_type& end() {return m_range.second;}
    private:
    T m_range;
  };
  
  /// square
  template<class T> inline constexpr T square( const T& x ) { return x*x; }

  //! radius
  template<class T> inline constexpr T get_r( T x, T y ) { return std::sqrt( square(x) + square(y) ); }

  //! pt
  template<class T> T get_pt( T px, T py ) { return std::sqrt( square(px) + square(py) ); }

  //! p
  template<class T> T get_p( T px, T py, T pz ) { return std::sqrt( square(px) + square(py) + square(pz) ); }

  //! eta
  template<class T> T get_eta( T p, T pz ) { return std::log( (p+pz)/(p-pz) )/2; }

  //! radius
  float get_r( PHG4Hit* hit, int i )
  {  return get_r( hit->get_x(i), hit->get_y(i) ); }

  //! calculate the interpolation of member function called on all members in collection to the provided r_extrap
  template< float (PHG4Hit::*accessor)(int) const>
  float interpolate( const std::set<PHG4Hit*>& hits, float rextrap )
  {
    
    // calculate all terms needed for the interpolation
    // need to use double everywhere here due to numerical divergences
    double sw = 0;
    double swr = 0;
    double swr2 = 0;
    double swx = 0;
    double swrx = 0;

    bool valid( false );
    for( const auto& hit:hits )
    {

      const double x0 = (hit->*accessor)(0);
      const double x1 = (hit->*accessor)(1);
      if( std::isnan( x0 ) || std::isnan( x1 ) ) continue;

      const double w = hit->get_edep();
      if( w < 0 ) continue;

      valid = true;
      const double r0 = get_r( hit, 0 );
      const double r1 = get_r( hit, 1 );

      sw += w*2;
      swr += w*(r0 + r1);
      swr2 += w*(square(r0) + square(r1));
      swx += w*(x0 + x1);
      swrx += w*(r0*x0 + r1*x1);
    }

    if( !valid ) return NAN;

    const auto alpha = (sw*swrx - swr*swx);
    const auto beta = (swr2*swx - swr*swrx);
    const auto denom = (sw*swr2 - square(swr));

    return ( alpha*rextrap + beta )/denom;
  }

  //! needed for weighted linear interpolation
  struct interpolation_data_t 
  {
    using list = std::vector<interpolation_data_t>;
    double x() const { return position.x(); }
    double y() const { return position.y(); }
    double z() const { return position.z(); }
    
    double px() const { return momentum.x(); }
    double py() const { return momentum.y(); }
    double pz() const { return momentum.z(); }
    
    TVector3 position;
    TVector3 momentum;
    double weight = 1;
  };
  
  //! calculate the interpolation of member function called on all members in collection to the provided y_extrap
  template<double (interpolation_data_t::*accessor)() const>
  double interpolate( const interpolation_data_t::list& hits, double y_extrap )
  {
    
    // calculate all terms needed for the interpolation
    // need to use double everywhere here due to numerical divergences
    double sw = 0;
    double swy = 0;
    double swy2 = 0;
    double swx = 0;
    double swyx = 0;

    bool valid( false );
    for( const auto& hit:hits )
    {

      const double x = (hit.*accessor)();
      const double w = hit.weight;
      if( w <= 0 ) continue;

      valid = true;
      const double y = hit.y();
        
      sw += w;
      swy += w*y;
      swy2 += w*square(y);
      swx += w*x;
      swyx += w*x*y;
    }

    if( !valid ) return NAN;

    const auto alpha = (sw*swyx - swy*swx);
    const auto beta = (swy2*swx - swy*swyx);
    const auto denom = (sw*swy2 - square(swy));

    return ( alpha*y_extrap + beta )/denom;
  }
  
  //! true if a track is a primary
  inline int is_primary( PHG4Particle* particle )
  { return particle->get_parent_id() == 0; }

  //! get mask from track clusters
  int64_t get_mask( SvtxTrack* track )
  { return std::accumulate( track->begin_cluster_keys(), track->end_cluster_keys(), int64_t(0),
      []( int64_t value, const TrkrDefs::cluskey& key ) {
        return TrkrDefs::getLayer(key)<64 ? value|(1LL<<TrkrDefs::getLayer(key)) : value;
      } );
  }

  //! return number of clusters of a given type
  template<int type>
    int get_clusters( SvtxTrack* track )
  {
    return std::count_if( track->begin_cluster_keys(), track->end_cluster_keys(),
      []( const TrkrDefs::cluskey& key ) { return TrkrDefs::getTrkrId(key) == type; } );
  }


  //! true if given layer is in mask
  [[maybe_unused]] bool has_layer( int64_t mask, int layer )
  { return mask & (1LL<<layer); }

  //! create track struct from struct from svx track
  TrackingEvaluator_hp::TrackStruct create_track( SvtxTrack* track )
  {
    TrackingEvaluator_hp::TrackStruct trackStruct;

    trackStruct._charge = track->get_charge();
    trackStruct._nclusters = track->size_cluster_keys();
    trackStruct._mask = get_mask( track );
    trackStruct._nclusters_mvtx = get_clusters<TrkrDefs::mvtxId>( track );
    trackStruct._nclusters_intt = get_clusters<TrkrDefs::inttId>( track );
    trackStruct._nclusters_tpc = get_clusters<TrkrDefs::tpcId>( track );
    trackStruct._nclusters_micromegas = get_clusters<TrkrDefs::micromegasId>( track );

    trackStruct._chisquare = track->get_chisq();
    trackStruct._ndf = track->get_ndf();

    trackStruct._x = track->get_x();
    trackStruct._y = track->get_y();
    trackStruct._z = track->get_z();
    trackStruct._r = get_r( trackStruct._x, trackStruct._y );
    trackStruct._phi = std::atan2( trackStruct._y, trackStruct._x );

    trackStruct._px = track->get_px();
    trackStruct._py = track->get_py();
    trackStruct._pz = track->get_pz();
    trackStruct._pt = get_pt( trackStruct._px, trackStruct._py );
    trackStruct._p = get_p( trackStruct._px, trackStruct._py, trackStruct._pz );
    trackStruct._eta = get_eta( trackStruct._p, trackStruct._pz );

    return trackStruct;
  }

  //! create track struct from struct from svx track
  TrackingEvaluator_hp::TrackPairStruct create_track_pair( SvtxTrack* first, SvtxTrack* second )
  {
    TrackingEvaluator_hp::TrackPairStruct trackpair_struct;

    trackpair_struct._charge = first->get_charge() + second->get_charge();
    trackpair_struct._px = first->get_px() + second->get_px();
    trackpair_struct._py = first->get_py() + second->get_py();
    trackpair_struct._pz = first->get_pz() + second->get_pz();
    trackpair_struct._pt = get_pt( trackpair_struct._px, trackpair_struct._py );
    trackpair_struct._p = get_p( trackpair_struct._px, trackpair_struct._py, trackpair_struct._pz );
    trackpair_struct._eta = get_eta( trackpair_struct._p, trackpair_struct._pz );

    // electron mass
    static constexpr double electronMass = 0.511e-3;
    auto firstE = std::sqrt( square( electronMass ) + square( get_p( first->get_px(), first->get_py(), first->get_pz() ) ) );
    auto secondE = std::sqrt( square( electronMass ) + square( get_p( second->get_px(), second->get_py(), second->get_pz() ) ) );
    trackpair_struct._e = firstE + secondE;
    trackpair_struct._m = std::sqrt( square( trackpair_struct._e ) - square( trackpair_struct._p ) );

    trackpair_struct._trk_pt[0] = get_pt( first->get_px(), first->get_py()  );
    trackpair_struct._trk_pt[1] = get_pt( second->get_px(), second->get_py()  );

    return trackpair_struct;
  }

  //! create cluster struct from svx cluster
  TrackingEvaluator_hp::ClusterStruct create_cluster( TrkrDefs::cluskey key, TrkrCluster* cluster )
  {
    TrackingEvaluator_hp::ClusterStruct cluster_struct;
    cluster_struct._layer = TrkrDefs::getLayer(key);
    cluster_struct._x = cluster->getX();
    cluster_struct._y = cluster->getY();
    cluster_struct._z = cluster->getZ();
    cluster_struct._r = get_r( cluster_struct._x, cluster_struct._y );
    cluster_struct._phi = std::atan2( cluster_struct._y, cluster_struct._x );
    cluster_struct._phi_error = cluster->getPhiError();
    cluster_struct._z_error = cluster->getZError();

    if(TrkrDefs::getTrkrId(key) == TrkrDefs::micromegasId)
    { cluster_struct._tileid = MicromegasDefs::getTileId(key); }
    
    return cluster_struct;
  }

  //! number of hits associated to cluster
  void add_cluster_size( TrackingEvaluator_hp::ClusterStruct& cluster, TrkrDefs::cluskey clus_key, TrkrClusterHitAssoc* cluster_hit_map )
  {
    if( !cluster_hit_map ) return;
    const auto range = cluster_hit_map->getHits(clus_key);

    // store full size
    cluster._size =  std::distance( range.first, range.second );

    const auto detId = TrkrDefs::getTrkrId(clus_key);
    if(detId == TrkrDefs::micromegasId)
    {

      // for micromegas the directional cluster size depends on segmentation type
      auto segmentation_type = MicromegasDefs::getSegmentationType(clus_key);
      if( segmentation_type == MicromegasDefs::SegmentationType::SEGMENTATION_Z ) cluster._z_size = cluster._size;
      else cluster._phi_size = cluster._size;

    } else {

      // for other detectors, one must loop over the constituting hits
      std::set<int> phibins;
      std::set<int> zbins;
      for(const auto& [first, hit_key]:range_adaptor(range))
      {
        switch( detId )
        {
          default: break;
          case TrkrDefs::mvtxId:
          {
            phibins.insert( MvtxDefs::getRow( hit_key ) );
            zbins.insert( MvtxDefs::getCol( hit_key ) );
            break;
          }
          case TrkrDefs::inttId:
          {
            phibins.insert( InttDefs::getRow( hit_key ) );
            zbins.insert( InttDefs::getCol( hit_key ) );
            break;
          }
          case TrkrDefs::tpcId:
          {
            phibins.insert( TpcDefs::getPad( hit_key ) );
            zbins.insert( TpcDefs::getTBin( hit_key ) );
            break;
          }
        }
      }
      cluster._phi_size = phibins.size();
      cluster._z_size = zbins.size();
    }
  }


  //! hit energy for a given cluster
  void add_cluster_energy( TrackingEvaluator_hp::ClusterStruct& cluster, TrkrDefs::cluskey clus_key,
    TrkrClusterHitAssoc* cluster_hit_map,
    TrkrHitSetContainer* hitsetcontainer )
  {

    // check container
    if(!(cluster_hit_map && hitsetcontainer)) return;

    // for now this is only filled for micromegas
    const auto detId = TrkrDefs::getTrkrId(clus_key);
    if(detId != TrkrDefs::micromegasId) return;

    const auto hitset_key = TrkrDefs::getHitSetKeyFromClusKey(clus_key);
    const auto hitset = hitsetcontainer->findHitSet( hitset_key );
    if( !hitset ) return;

    const auto range = cluster_hit_map->getHits(clus_key);
    cluster._energy_max = 0;
    cluster._energy_sum = 0;

    for( const auto& pair:range_adaptor(range))
    {
      const auto hit = hitset->getHit( pair.second );
      if( hit )
      {
        const auto energy = hit->getEnergy();
        cluster._energy_sum += energy;
        if( energy > cluster._energy_max ) cluster._energy_max = energy;
      }
    }

  }

  // add truth information
  void add_truth_information( TrackingEvaluator_hp::TrackStruct& track, PHG4Particle* particle )
  {
    if( particle )
    {
      track._is_primary = is_primary( particle );
      track._pid = particle->get_pid();
      track._truth_px = particle->get_px();
      track._truth_py = particle->get_py();
      track._truth_pz = particle->get_pz();
      track._truth_pt = get_pt( track._truth_px, track._truth_py );
      track._truth_p = get_p( track._truth_px, track._truth_py, track._truth_pz );
      track._truth_eta = get_eta( track._truth_p, track._truth_pz );
    }
  }

  // print to stream
  [[maybe_unused]] std::ostream& operator << (std::ostream& out, const TrackingEvaluator_hp::ClusterStruct& cluster )
  {
    out << "ClusterStruct" << std::endl;
    out << "  cluster: (" << cluster._x << "," << cluster._y << "," << cluster._z << ")" << std::endl;
    out << "  track:   (" << cluster._trk_x << "," << cluster._trk_y << "," << cluster._trk_z << ")" << std::endl;
    out << "  truth:   (" << cluster._truth_x << "," << cluster._truth_y << "," << cluster._truth_z << ")" << std::endl;
    return out;
  }

  [[maybe_unused]] std::ostream& operator << (std::ostream& out, const TVector3& position)
  {
    out << "(" << position.x() << ", " << position.y() << ", " << position.z() << ")";
    return out;
  }
  
}

//_____________________________________________________________________
void TrackingEvaluator_hp::Container::Reset()
{
  _events.clear();
  _clusters.clear();
  _tracks.clear();
  _track_pairs.clear();
}

//_____________________________________________________________________
TrackingEvaluator_hp::TrackingEvaluator_hp( const std::string& name ):
  SubsysReco( name)
{}

//_____________________________________________________________________
int TrackingEvaluator_hp::Init(PHCompositeNode* topNode )
{
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

  auto newNode = new PHIODataNode<PHObject>( new Container, "TrackingEvaluator_hp::Container","PHObject");
  evalNode->addNode(newNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int TrackingEvaluator_hp::InitRun(PHCompositeNode* )
{ return Fun4AllReturnCodes::EVENT_OK; }

//_____________________________________________________________________
int TrackingEvaluator_hp::process_event(PHCompositeNode* topNode)
{
  // load nodes
  auto res =  load_nodes(topNode);
  if( res != Fun4AllReturnCodes::EVENT_OK ) return res;

  // cleanup output
  if( m_container ) m_container->Reset();

  if(m_flags&PrintClusters) print_clusters();
  if(m_flags&PrintTracks) print_tracks();

  if(m_flags&EvalEvent) evaluate_event();
  if(m_flags&EvalClusters) evaluate_clusters();
  if(m_flags&EvalTracks) evaluate_tracks();
  if(m_flags&EvalTrackPairs) evaluate_track_pairs();

  // clear maps
  m_g4hit_map.clear();
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int TrackingEvaluator_hp::End(PHCompositeNode* )
{ return Fun4AllReturnCodes::EVENT_OK; }

//_____________________________________________________________________
int TrackingEvaluator_hp::load_nodes( PHCompositeNode* topNode )
{

  // get necessary nodes
  m_track_map = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");

  // cluster map
  m_cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");

  // cluster hit association map
  m_cluster_hit_map = findNode::getClass<TrkrClusterHitAssoc>(topNode, "TRKR_CLUSTERHITASSOC");

  // cluster hit association map
  m_hit_truth_map = findNode::getClass<TrkrHitTruthAssoc>(topNode,"TRKR_HITTRUTHASSOC");

  // local container
  m_container = findNode::getClass<Container>(topNode, "TrackingEvaluator_hp::Container");

  // hitset container
  m_hitsetcontainer = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");

  // g4hits
  m_g4hits_tpc = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_TPC");
  m_g4hits_intt = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_INTT");
  m_g4hits_mvtx = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_MVTX");
  m_g4hits_micromegas = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_MICROMEGAS");

  // g4 truth info
  m_g4truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

  // micromegas geometry
  m_micromegas_geonode = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MICROMEGAS_FULL" );
  assert( m_micromegas_geonode );
  
  return Fun4AllReturnCodes::EVENT_OK;

}

//_____________________________________________________________________
void TrackingEvaluator_hp::evaluate_event()
{
  if(!m_container) return;

  // TPC threshold (adc)
  // from PHG4TpcDigitizer::InitRun - ADCThreshold * ADCNoiseConversionGain*1024/2200
  static constexpr float tpc_threshold = 78.4645;

  // create event struct
  EventStruct event;

  if( m_hitsetcontainer )
  {
    // fill hit related information
    for(const auto& [hitsetkey,hitset]:range_adaptor(m_hitsetcontainer->getHitSets()))
    {
      const auto trkrId = TrkrDefs::getTrkrId(hitsetkey);
      const auto layer = TrkrDefs::getLayer(hitsetkey);
      assert(layer<EventStruct::max_layer);

      const auto hit_range = hitset->getHits();

      // nhits per layer
      if( trkrId == TrkrDefs::tpcId )
      {
        const auto accepted = std::count_if( hit_range.first, hit_range.second,
          []( const TrkrHitSet::Map::value_type pair )
          { return pair.second->getAdc() >= tpc_threshold; });
        event._nhits[layer] += accepted;
      } else {
        event._nhits[layer] += std::distance( hit_range.first, hit_range.second );
      }

      // nhits per layer, with no threshold cut for the TPC
      event._nhits_raw[layer] += std::distance( hit_range.first, hit_range.second );

      // also count corresponding clusters
      if(m_cluster_map)
      {

        // fill cluster related information
        for( const auto[key,cluster]:range_adaptor(m_cluster_map->getClusters(hitsetkey)) )
        {
          const auto trkrid = TrkrDefs::getTrkrId(key);
          switch( trkrid )
          {
            case TrkrDefs::mvtxId: ++event._nclusters_mvtx; break;
            case TrkrDefs::inttId: ++event._nclusters_intt; break;
            case TrkrDefs::tpcId: ++event._nclusters_tpc; break;
            case TrkrDefs::micromegasId: ++event._nclusters_micromegas; break;
          }

          const auto layer = static_cast<size_t>(TrkrDefs::getLayer(key));
          assert(layer<EventStruct::max_layer);
          ++event._nclusters[layer];
        }
      }
    }
  }

  // store
  m_container->addEvent(event);
}

//_____________________________________________________________________
void TrackingEvaluator_hp::evaluate_clusters()
{

  if(!(m_cluster_map&&m_hitsetcontainer&&m_container)) return;

  // clear array
  m_container->clearClusters();

  // first loop over hitsets
  for( const auto& [hitsetkey,hitset]:range_adaptor(m_hitsetcontainer->getHitSets()))
  {
    for( const auto& [key,cluster]:range_adaptor(m_cluster_map->getClusters(hitsetkey)))
    {
      // create cluster structure
      auto cluster_struct = create_cluster( key, cluster );
      add_cluster_size( cluster_struct, key, m_cluster_hit_map );
      add_cluster_energy( cluster_struct, key, m_cluster_hit_map, m_hitsetcontainer );

      // truth information
      const auto g4hits = find_g4hits( key );
      const bool is_micromegas( TrkrDefs::getTrkrId(key) == TrkrDefs::micromegasId );
      if( is_micromegas ) add_truth_information_micromegas( cluster_struct, g4hits );
      else add_truth_information( cluster_struct, g4hits );

      // add in array
      m_container->addCluster( cluster_struct );
    }
  }
}

//_____________________________________________________________________
void TrackingEvaluator_hp::evaluate_tracks()
{
  if( !( m_track_map && m_cluster_map && m_container ) ) return;

  // clear array
  m_container->clearTracks();

  for( const auto& trackpair:*m_track_map )
  {

    const auto track = trackpair.second;
    auto track_struct = create_track( track );

    // truth information
    const auto pair = get_max_contributor( track );
    const auto id = pair.first;
    track_struct._contributors = pair.second;

    auto particle = m_g4truthinfo->GetParticle(id);
    track_struct._embed = get_embed(particle);
    ::add_truth_information(track_struct, particle);

    // loop over clusters
    for( auto key_iter = track->begin_cluster_keys(); key_iter != track->end_cluster_keys(); ++key_iter )
    {

      const auto& cluster_key = *key_iter;
      auto cluster = m_cluster_map->findCluster( cluster_key );
      if( !cluster )
      {
        std::cout << "TrackingEvaluator_hp::evaluate_tracks - unable to find cluster for key " << cluster_key << std::endl;
        continue;
      }

      // create new cluster struct
      auto cluster_struct = create_cluster( cluster_key, cluster );
      add_cluster_size( cluster_struct, cluster_key, m_cluster_hit_map );
      add_cluster_energy( cluster_struct, cluster_key, m_cluster_hit_map, m_hitsetcontainer );

      // truth information
      const auto g4hits = find_g4hits( cluster_key );

      const bool is_micromegas( TrkrDefs::getTrkrId(cluster_key) == TrkrDefs::micromegasId );
      if( is_micromegas ) add_truth_information_micromegas( cluster_struct, g4hits );
      else add_truth_information( cluster_struct, g4hits );

      // find track state that is the closest to cluster
      /* this assumes that both clusters and states are sorted along r */
      const auto radius( cluster_struct._r );

      float dr_min = -1;
      auto state_iter = track->begin_states();
      for( auto iter = state_iter; iter != track->end_states(); ++iter )
      {
        const auto dr = std::abs( radius - get_r( iter->second->get_x(), iter->second->get_y() ) );
        if( dr_min < 0 || dr < dr_min )
        {
          state_iter = iter;
          dr_min = dr;
        }
      }

      // store track state in cluster struct
      if( is_micromegas ) add_trk_information_micromegas( cluster_struct, state_iter->second );
      else add_trk_information( cluster_struct, state_iter->second );

      // add to track
      track_struct._clusters.push_back( cluster_struct );
    }
    m_container->addTrack( track_struct );
  }
}

//_____________________________________________________________________
void TrackingEvaluator_hp::evaluate_track_pairs()
{
  if( !( m_track_map && m_container ) ) return;

  // clear array
  m_container->clearTrackPairs();
  for( auto firstIter = m_track_map->begin(); firstIter != m_track_map->end(); ++firstIter )
  {
    const auto first = firstIter->second;
    for( auto  secondIter = m_track_map->begin(); secondIter != firstIter; ++secondIter )
    {
      const auto second = secondIter->second;
      auto trackpair_struct = create_track_pair( first, second );
      trackpair_struct._contributors[0] = get_max_contributor( first ).second;
      trackpair_struct._contributors[1] = get_max_contributor( second ).second;
      m_container->addTrackPair( trackpair_struct );
    }
  }
}

//_____________________________________________________________________
void TrackingEvaluator_hp::print_clusters() const
{

  if(!(m_cluster_map && m_hitsetcontainer)) return;

  for(const auto& [hitsetkey,hitset]:range_adaptor(m_hitsetcontainer->getHitSets()))
  {
    // get corresponding clusters
    for(const auto& [clusterkey,cluster]:range_adaptor(m_cluster_map->getClusters(hitsetkey)))
    {
      // only print for TPC ids
      if(TrkrDefs::getTrkrId(clusterkey) == TrkrDefs::tpcId)
      { print_cluster( clusterkey, cluster ); }
    }
  }
}

//_____________________________________________________________________
void TrackingEvaluator_hp::print_tracks() const
{
  if( !m_track_map ) return;
  for(const auto& trackpair:*m_track_map)
  { print_track( trackpair.second ); }
}

//_____________________________________________________________________
void TrackingEvaluator_hp::print_cluster( TrkrDefs::cluskey cluster_key, TrkrCluster* cluster ) const
{
  // get detector type
  const auto trkrId = TrkrDefs::getTrkrId( cluster_key );
  std::cout
    << "TrackingEvaluator_hp::print_cluster -"
    << " layer: " << (int)TrkrDefs::getLayer(cluster_key)
    << " type: " << (int) trkrId
    // << " position: (" << cluster->getX() << "," << cluster->getY() << "," << cluster->getZ() << ")"
    << " polar: (" << get_r( cluster->getX(), cluster->getY()) << "," << std::atan2( cluster->getY(), cluster->getX()) << "," << cluster->getZ() << ")"
    << " errors: (" << cluster->getPhiError() << ", " << cluster->getZError() << ")"
    << std::endl;

  // get associated hits
  if( false )
  {

    // loop over hits associated to clusters
    const auto range = m_cluster_hit_map->getHits(cluster_key);
    std::cout
      << "TrackingEvaluator_hp::print_cluster -"
      << " hit_count: " << std::distance( range.first, range.second )
      << std::endl;

    for(const auto& [hitsetkey,hitkey]:range_adaptor(range))
    {

      switch( trkrId )
      {
        case TrkrDefs::tpcId:
        std::cout
          << "TrackingEvaluator_hp::print_cluster -"
          << " hit: " << hitkey
          << " sector: " << TpcDefs::getSectorId( cluster_key )
          << " pad: " << TpcDefs::getPad( hitkey )
          << " time bin: " << TpcDefs::getTBin( hitkey )
          << std::endl;
        break;

        default: break;

      }
    }
  }

  // get associated g4 hits
  if( false )
  {

    auto g4hits = find_g4hits( cluster_key );
    for( const auto& g4hit:g4hits )
    {

      std::cout
        << "TrackingEvaluator_hp::print_cluster -"
        << " layer: " << g4hit->get_layer()
        << " track: " << g4hit->get_trkid()
        << " in: (" << g4hit->get_x(0) << "," << g4hit->get_y(0) << "," << g4hit->get_z(0) << ")"
        << " out: (" << g4hit->get_x(1) << "," << g4hit->get_y(1) << "," << g4hit->get_z(1) << ")"
        << " polar in: (" << get_r( g4hit->get_x(0), g4hit->get_y(0) ) << "," << std::atan2( g4hit->get_y(0), g4hit->get_x(0) ) << "," << g4hit->get_z(0) << ")"
        << " polar out: (" << get_r( g4hit->get_x(1), g4hit->get_y(1) ) << "," << std::atan2( g4hit->get_y(1), g4hit->get_x(1) ) << "," << g4hit->get_z(1) << ")"
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
      << " polar: (" << get_r( xextrap, yextrap ) << "," << std::atan2( yextrap, xextrap ) << "," << zextrap << ")"
      << std::endl;

  }

  std::cout << std::endl;

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
  if( false && m_cluster_map )
  {
    for( auto key_iter = track->begin_cluster_keys(); key_iter != track->end_cluster_keys(); ++key_iter )
    {

      const auto& cluster_key = *key_iter;
      auto cluster = m_cluster_map->findCluster( cluster_key );
      if( !cluster )
      {
        std::cout << "TrackingEvaluator_hp::print_track - unable to find cluster for key " << cluster_key << std::endl;
        continue;
      }

      std::cout
        << "TrackingEvaluator_hp::print_track -"
        << " cluster layer: "  << (int)TrkrDefs::getLayer(cluster_key)
        << " position: (" << cluster->getX() << ", " << cluster->getY() << ", " << cluster->getZ() << ")"
        << " polar: (" << get_r( cluster->getX(), cluster->getY() ) << ", " << std::atan2( cluster->getY(), cluster->getX() ) << "," << cluster->getZ() << ")"
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
        << " polar: (" << get_r( state->get_x(), state->get_y() ) << ", " << std::atan2( state->get_y(), state->get_x() ) << "," << state->get_z() << ")"
        << std::endl;
    }
  }

  std::cout << std::endl;

}

//_____________________________________________________________________
TrackingEvaluator_hp::G4HitSet TrackingEvaluator_hp::find_g4hits( TrkrDefs::cluskey cluster_key ) const
{

  // check maps
  if( !( m_cluster_hit_map && m_hit_truth_map ) ) return G4HitSet();

  // check if in map
  auto map_iter = m_g4hit_map.lower_bound( cluster_key );
  if( map_iter != m_g4hit_map.end() && cluster_key == map_iter->first )
  { return map_iter->second; }

  // find hitset associated to cluster
  G4HitSet out;
  const auto hitset_key = TrkrDefs::getHitSetKeyFromClusKey(cluster_key);
  for(const auto& [first,hit_key]:range_adaptor( m_cluster_hit_map->getHits(cluster_key)))
  {

    // store hits to g4hit associations
    TrkrHitTruthAssoc::MMap g4hit_map;
    m_hit_truth_map->getG4Hits( hitset_key, hit_key, g4hit_map );

    // find corresponding g4 hist
    for( auto truth_iter = g4hit_map.begin(); truth_iter != g4hit_map.end(); ++truth_iter )
    {

      const auto g4hit_key = truth_iter->second.second;
      PHG4Hit* g4hit = nullptr;

      switch( TrkrDefs::getTrkrId( hitset_key ) )
      {
        case TrkrDefs::mvtxId:
        if( m_g4hits_mvtx ) g4hit = m_g4hits_mvtx->findHit( g4hit_key );
        break;

        case TrkrDefs::inttId:
        if( m_g4hits_intt ) g4hit = m_g4hits_intt->findHit( g4hit_key );
        break;

        case TrkrDefs::tpcId:
        if( m_g4hits_tpc ) g4hit = m_g4hits_tpc->findHit( g4hit_key );
        break;

        case TrkrDefs::micromegasId:
        if( m_g4hits_micromegas ) g4hit = m_g4hits_micromegas->findHit( g4hit_key );
        break;

        default: break;
      }

      if( g4hit ) out.insert( g4hit );
      else std::cout << "TrackingEvaluator_hp::find_g4hits - g4hit not found " << g4hit_key << std::endl;

    }
  }

  // insert in map and return
  return m_g4hit_map.insert( map_iter, std::make_pair( cluster_key, std::move( out ) ) )->second;

}

//_____________________________________________________________________
std::pair<int,int> TrackingEvaluator_hp::get_max_contributor( SvtxTrack* track ) const
{
  if(!(m_track_map && m_cluster_map && m_g4truthinfo)) return {0,0};

  // maps MC track id and number of matching g4hits
  using IdMap = std::map<int,int>;
  IdMap contributor_map;

  // loop over clusters
  for( auto key_iter = track->begin_cluster_keys(); key_iter != track->end_cluster_keys(); ++key_iter )
  {
    const auto& cluster_key = *key_iter;
    for( const auto& hit:find_g4hits( cluster_key ) )
    {
      const int trkid = hit->get_trkid();
      auto iter = contributor_map.lower_bound( trkid );
      if( iter == contributor_map.end() || iter->first != trkid )
      {
        contributor_map.insert(iter, std::make_pair(trkid,1));
      } else ++iter->second;
    }
  }

  if( contributor_map.empty() ) return {0,0};
  else return *std::max_element(
    contributor_map.cbegin(), contributor_map.cend(),
    []( const IdMap::value_type& first, const IdMap::value_type& second )
    { return first.second < second.second; } );

}

//_____________________________________________________________________
int TrackingEvaluator_hp::get_embed( PHG4Particle* particle ) const
{ return (m_g4truthinfo && particle) ? m_g4truthinfo->isEmbeded( particle->get_primary_id() ):0; }

//_____________________________________________________________________
void TrackingEvaluator_hp::add_trk_information( TrackingEvaluator_hp::ClusterStruct& cluster, SvtxTrackState* state ) const
{
  // need to extrapolate to the right r
  const auto trk_r = get_r( state->get_x(), state->get_y() );
  const auto dr = cluster._r - trk_r;
  const auto trk_drdt = (state->get_x()*state->get_px() + state->get_y()*state->get_py())/trk_r;
  const auto trk_dxdr = state->get_px()/trk_drdt;
  const auto trk_dydr = state->get_py()/trk_drdt;
  const auto trk_dzdr = state->get_pz()/trk_drdt;
  
  // store state position
  cluster._trk_x = state->get_x() + dr*trk_dxdr;
  cluster._trk_y = state->get_y() + dr*trk_dydr;
  cluster._trk_z = state->get_z() + dr*trk_dzdr;
  cluster._trk_r = get_r( cluster._trk_x, cluster._trk_y );
  cluster._trk_phi = std::atan2( cluster._trk_y, cluster._trk_x );
  
  /* store local momentum information */
  cluster._trk_px = state->get_px();
  cluster._trk_py = state->get_py();
  cluster._trk_pz = state->get_pz();
  
  /*
  store state angles in (r,phi) and (r,z) plans
  they are needed to study space charge distortions
  */
  const auto cosphi( std::cos( cluster._trk_phi ) );
  const auto sinphi( std::sin( cluster._trk_phi ) );
  const auto trk_pphi = -state->get_px()*sinphi + state->get_py()*cosphi;
  const auto trk_pr = state->get_px()*cosphi + state->get_py()*sinphi;
  const auto trk_pz = state->get_pz();
  cluster._trk_alpha = std::atan2( trk_pphi, trk_pr );
  cluster._trk_beta = std::atan2( trk_pz, trk_pr );
  cluster._trk_phi_error = state->get_phi_error();
  cluster._trk_z_error = state->get_z_error();
  
}

//_____________________________________________________________________
void TrackingEvaluator_hp::add_trk_information_micromegas( TrackingEvaluator_hp::ClusterStruct& cluster, SvtxTrackState* state ) const
{

  // get geometry cylinder from layer
  const auto layer = cluster._layer;
  const auto tileid = cluster._tileid;
  const auto layergeom = dynamic_cast<CylinderGeomMicromegas*>(m_micromegas_geonode->GetLayerGeom(layer));
  assert( layergeom );
  
  // convert cluster position to local tile coordinates
  const TVector3 cluster_world( cluster._x, cluster._y, cluster._z );
  const TVector3 cluster_local = layergeom->get_local_from_world_coords( tileid, cluster_world );

  // convert track position to local tile coordinates
  TVector3 track_world( state->get_x(), state->get_y(), state->get_z() );
  TVector3 track_local = layergeom->get_local_from_world_coords( tileid, track_world );

  // convert direction to local tile coordinates
  const TVector3 direction_world( state->get_px(), state->get_py(), state->get_pz() );
  const TVector3 direction_local = layergeom->get_local_from_world_vect( tileid, direction_world );
  
  // extrapolate to same local y (should be zero) as cluster
  const auto delta_y = cluster_local.y() - track_local.y();
  track_local += TVector3( 
    delta_y*direction_local.x()/direction_local.y(),
    delta_y, 
    delta_y*direction_local.z()/direction_local.y() );

  // convert back to global coordinates
  track_world = layergeom->get_world_from_local_coords( tileid, track_local );
  
  // store state position
  cluster._trk_x = track_world.x();
  cluster._trk_y = track_world.y();
  cluster._trk_z = track_world.z();
  cluster._trk_r = get_r( cluster._trk_x, cluster._trk_y );
  cluster._trk_phi = std::atan2( cluster._trk_y, cluster._trk_x );
  
  /* store local momentum information */
  cluster._trk_px = state->get_px();
  cluster._trk_py = state->get_py();
  cluster._trk_pz = state->get_pz();
  
  /*
  store state angles in (r,phi) and (r,z) plans
  they are needed to study space charge distortions
  */
  const auto cosphi( std::cos( cluster._trk_phi ) );
  const auto sinphi( std::sin( cluster._trk_phi ) );
  const auto trk_pphi = -state->get_px()*sinphi + state->get_py()*cosphi;
  const auto trk_pr = state->get_px()*cosphi + state->get_py()*sinphi;
  const auto trk_pz = state->get_pz();
  cluster._trk_alpha = std::atan2( trk_pphi, trk_pr );
  cluster._trk_beta = std::atan2( trk_pz, trk_pr );
  cluster._trk_phi_error = state->get_phi_error();
  cluster._trk_z_error = state->get_z_error();
  
}

//_____________________________________________________________________
void TrackingEvaluator_hp::add_truth_information( TrackingEvaluator_hp::ClusterStruct& cluster, std::set<PHG4Hit*> hits ) const
{
  const auto rextrap = cluster._r;
  cluster._truth_size = hits.size();
  cluster._truth_x = interpolate<&PHG4Hit::get_x>( hits, rextrap );
  cluster._truth_y = interpolate<&PHG4Hit::get_y>( hits, rextrap );
  cluster._truth_z = interpolate<&PHG4Hit::get_z>( hits, rextrap );
  cluster._truth_r = get_r( cluster._truth_x, cluster._truth_y );
  cluster._truth_phi = std::atan2( cluster._truth_y, cluster._truth_x );
  
  /* add truth momentum information */
  cluster._truth_px = interpolate<&PHG4Hit::get_px>( hits, rextrap );
  cluster._truth_py = interpolate<&PHG4Hit::get_py>( hits, rextrap );
  cluster._truth_pz = interpolate<&PHG4Hit::get_pz>( hits, rextrap );
  
  /*
  store state angles in (r,phi) and (r,z) plans
  they are needed to study space charge distortions
  */
  const auto cosphi( std::cos( cluster._truth_phi ) );
  const auto sinphi( std::sin( cluster._truth_phi ) );
  const auto truth_pphi = -cluster._truth_px*sinphi + cluster._truth_py*cosphi;
  const auto truth_pr = cluster._truth_px*cosphi + cluster._truth_py*sinphi;
  
  cluster._truth_alpha = std::atan2( truth_pphi, truth_pr );
  cluster._truth_beta = std::atan2( cluster._truth_pz, truth_pr );
  
}

//_____________________________________________________________________
void TrackingEvaluator_hp::add_truth_information_micromegas( TrackingEvaluator_hp::ClusterStruct& cluster, std::set<PHG4Hit*> g4hits ) const
{
  const auto layer = cluster._layer;
  const auto tileid = cluster._tileid;
  const auto layergeom = dynamic_cast<CylinderGeomMicromegas*>(m_micromegas_geonode->GetLayerGeom(layer));
  assert( layergeom );

  // convert cluster position to local tile coordinates
  const TVector3 cluster_world( cluster._x, cluster._y, cluster._z );
  const TVector3 cluster_local = layergeom->get_local_from_world_coords( tileid, cluster_world );
  
  // std::cout << "TrackingEvaluator_hp::add_truth_information_micromegas - cluster_local: " << cluster_local << std::endl;
  
  // convert hits to list of interpolation_data_t
  interpolation_data_t::list hits;
  for( const auto& g4hit:g4hits )
  {
    const auto weight = g4hit->get_edep();
    for( int i = 0; i < 2; ++i )
    {
      
      // convert position to local
      TVector3 g4hit_world(g4hit->get_x(i), g4hit->get_y(i), g4hit->get_z(i));
      TVector3 g4hit_local = layergeom->get_local_from_world_coords( tileid, g4hit_world );
      
      // convert momentum to local
      TVector3 momentum_world(g4hit->get_px(i), g4hit->get_py(i), g4hit->get_pz(i));
      TVector3 momentum_local = layergeom->get_local_from_world_vect( tileid, momentum_world );
      
      hits.push_back( {.position = g4hit_local, .momentum = momentum_local, .weight = weight } );
    }
  }
  
  // do position interpolation
  const auto y_extrap = cluster_local.y();
  const TVector3 interpolation_local( 
    interpolate<&interpolation_data_t::x>( hits, y_extrap ),
    interpolate<&interpolation_data_t::y>( hits, y_extrap ),
    interpolate<&interpolation_data_t::z>( hits, y_extrap ) );
    
  const TVector3 interpolation_world = layergeom->get_world_from_local_coords( tileid, interpolation_local );

  // do momentum interpolation
  const TVector3 momentum_local( 
    interpolate<&interpolation_data_t::px>( hits, y_extrap ),
    interpolate<&interpolation_data_t::py>( hits, y_extrap ),
    interpolate<&interpolation_data_t::pz>( hits, y_extrap ) );
  
  const TVector3 momentum_world = layergeom->get_world_from_local_vect( tileid, momentum_local );

  // update cluster structure
  cluster._truth_size = hits.size();
  cluster._truth_x = interpolation_world.x();
  cluster._truth_y = interpolation_world.y();
  cluster._truth_z = interpolation_world.z();
  cluster._truth_r = get_r( cluster._truth_x, cluster._truth_y );
  cluster._truth_phi = std::atan2( cluster._truth_y, cluster._truth_x );
  
  /* add truth momentum information */
  cluster._truth_px = momentum_world.x();
  cluster._truth_py = momentum_world.y();
  cluster._truth_pz = momentum_world.z();
  
  /*
  store state angles in (r,phi) and (r,z) plans
  they are needed to study space charge distortions
  */
  const auto cosphi( std::cos( cluster._truth_phi ) );
  const auto sinphi( std::sin( cluster._truth_phi ) );
  const auto truth_pphi = -cluster._truth_px*sinphi + cluster._truth_py*cosphi;
  const auto truth_pr = cluster._truth_px*cosphi + cluster._truth_py*sinphi;
  
  cluster._truth_alpha = std::atan2( truth_pphi, truth_pr );
  cluster._truth_beta = std::atan2( cluster._truth_pz, truth_pr );
  
}
