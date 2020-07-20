#include "SimEvaluator_hp.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>
#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>

#include <Geant4/G4SystemOfUnits.hh>

#include <algorithm>
#include <bitset>
#include <iostream>
#include <numeric>

//_____________________________________________________________________
namespace
{

  /// square
  template<class T> T square( T x ) { return x*x; }

  /// pt
  template<class T> T get_pt( T px, T py ) { return std::sqrt( square(px) + square(py) ); }

  /// p
  template<class T> T get_p( T px, T py, T pz ) { return std::sqrt( square(px) + square(py) + square(pz) ); }

  /// eta
  template<class T> T get_eta( T p, T pz ) { return std::log( (p+pz)/(p-pz) )/2; }

  /// true if particle is primary
  inline bool is_primary( PHG4Particle* particle )
  { return particle->get_parent_id() == 0; }

  //_____________________________________________________________________
  SimEvaluator_hp::EventStruct create_event(PHG4TruthInfoContainer* container)
  {
    SimEvaluator_hp::EventStruct eventStruct;
    if( container ) 
    {
      std::set<int> flags;
      const auto range = container->GetEmbeddedTrkIds();
      for( auto iter = range.first; iter != range.second; ++iter )
      { flags.insert( iter->second ); }

      // print
      std::cout << "::create_event - flags: ";
      for( const auto& flag:flags ) std::cout << flag << " ";
      std::cout << std::endl;

      for( const auto& flag:flags )
      {
        ++eventStruct._nevt;
        if( flag < 0 ) ++eventStruct._nevt_bg;
        else ++eventStruct._nevt_active;
      }
    }
    
    return eventStruct;
    
  }
  
  //_____________________________________________________________________
  /// create track struct from struct from svx track
  SimEvaluator_hp::VertexStruct create_vertex( PHG4VtxPoint* vertex )
  {
    SimEvaluator_hp::VertexStruct vertexStruct;
    vertexStruct._x = vertex->get_x();
    vertexStruct._y = vertex->get_y();
    vertexStruct._z = vertex->get_z();
    vertexStruct._t = vertex->get_t();
    return vertexStruct;
  }

  //_____________________________________________________________________
  /// create track struct from struct from svx track
  SimEvaluator_hp::ParticleStruct create_particle( PHG4Particle* particle )
  {
    SimEvaluator_hp::ParticleStruct particleStruct;
    particleStruct._pid = particle->get_pid();
    particleStruct._charge = particle->get_IonCharge()/eplus;
    particleStruct._is_primary = is_primary( particle );
    particleStruct._px = particle->get_px();
    particleStruct._py = particle->get_py();
    particleStruct._pz = particle->get_pz();
    particleStruct._pt = get_pt( particle->get_px(), particle->get_py() );
    particleStruct._p = get_p( particle->get_px(), particle->get_py(), particle->get_pz() );
    particleStruct._eta = get_eta( particleStruct._p, particleStruct._pz );
    return particleStruct;
  }

  //_____________________________________________________________________
  std::ostream& operator << (std::ostream& out, const PHG4VtxPoint& vertex )
  {
    out << "( " << vertex.get_x() << ", " << vertex.get_y() << ", " << vertex.get_z() << ", " << vertex.get_t() << ")";
    return out;
  }

}

//_____________________________________________________________________
void SimEvaluator_hp::Container::Reset()
{
  _vertex_list.clear();
  _particle_list.clear();
}

//_____________________________________________________________________
SimEvaluator_hp::SimEvaluator_hp( const std::string& name ):
  SubsysReco( name)
{}

//_____________________________________________________________________
int SimEvaluator_hp::Init(PHCompositeNode* topNode )
{

  // find DST node
  PHNodeIterator iter(topNode);
  auto dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << "SimEvaluator_hp::Init - DST Node missing" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // get EVAL node
  iter = PHNodeIterator(dstNode);
  auto evalNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "EVAL"));
  if( !evalNode )
  {
    // create
    std::cout << "SimEvaluator_hp::Init - EVAL node missing - creating" << std::endl;
    evalNode = new PHCompositeNode( "EVAL" );
    dstNode->addNode(evalNode);
  }

  auto newNode = new PHIODataNode<PHObject>( new Container, "SimEvaluator_hp::Container", "PHObject" );
  evalNode->addNode(newNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int SimEvaluator_hp::InitRun(PHCompositeNode* )
{ return Fun4AllReturnCodes::EVENT_OK; }

//_____________________________________________________________________
int SimEvaluator_hp::process_event(PHCompositeNode* topNode)
{
  // load nodes
  auto res =  load_nodes(topNode);
  if( res != Fun4AllReturnCodes::EVENT_OK ) return res;

  fill_event();
  fill_vertices();
  fill_g4particle_map();
  fill_particles();
  // print_vertices();

  m_g4particle_map.clear();

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int SimEvaluator_hp::End(PHCompositeNode* )
{ return Fun4AllReturnCodes::EVENT_OK; }

//_____________________________________________________________________
int SimEvaluator_hp::load_nodes( PHCompositeNode* topNode )
{

  // local container
  m_container = findNode::getClass<Container>(topNode, "SimEvaluator_hp::Container");

  // g4 truth info
  m_g4truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

  // g4hits
  m_g4hits_tpc = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_TPC");
  m_g4hits_intt = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_INTT");
  m_g4hits_mvtx = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_MVTX");
  m_g4hits_micromegas = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_MICROMEGAS");

  return Fun4AllReturnCodes::EVENT_OK;

}


//_____________________________________________________________________
void SimEvaluator_hp::fill_g4particle_map()
{
  m_g4particle_map.clear();
  for( const auto& container: {m_g4hits_tpc, m_g4hits_intt, m_g4hits_mvtx, m_g4hits_micromegas} )
  {
    if( !container ) continue;

    // loop over hits
    const auto range = container->getHits();
    std::cout << "SimEvaluator_hp::fill_g4particle_map - counts: " << std::distance( range.first, range.second ) << std::endl;
    for( auto iter = range.first; iter != range.second; ++iter )
    {
      const auto map_iter = m_g4particle_map.lower_bound( iter->second->get_trkid() );
      if( map_iter != m_g4particle_map.end() && map_iter->first == iter->second->get_trkid() )
      {

        map_iter->second |= (1LL<<iter->second->get_layer());

      } else {

        m_g4particle_map.insert( map_iter, std::make_pair( iter->second->get_trkid(), 1LL<<iter->second->get_layer() ) );

      }

    }

  }

}

//_____________________________________________________________________
void SimEvaluator_hp::fill_event()
{

  if( !( m_container && m_g4truthinfo ) )
  {
    std::cerr << "SimEvaluator_hp::fill_event - nodes not found." << std::endl;
    return;
  }

  // clear vertices from previous event
  m_container->clearEventList();
  m_container->addEvent( create_event( m_g4truthinfo ) );

}

//_____________________________________________________________________
void SimEvaluator_hp::fill_vertices()
{

  if( !( m_container && m_g4truthinfo ) )
  {
    std::cerr << "SimEvaluator_hp::fill_vertices - nodes not found." << std::endl;
    return;
  }

  // clear vertices from previous event
  m_container->clearVertexList();

  // get main primary vertex id
  const auto main_vertex_id = m_g4truthinfo->GetPrimaryVertexIndex();

  auto range = m_g4truthinfo->GetPrimaryVtxRange();
  for( auto iter = range.first; iter != range.second; ++iter )
  {
    auto vertex = iter->second;
    if( vertex )
    {
      auto vertexStruct = create_vertex( vertex );
      vertexStruct._is_main_vertex = ( vertex->get_id() == main_vertex_id );
      m_container->addVertex( vertexStruct );
    }

  }

}

//_____________________________________________________________________
void SimEvaluator_hp::fill_particles()
{

  if( !( m_container && m_g4truthinfo ) )
  {
    std::cerr << "SimEvaluator_hp::fill_vertices - nodes not found." << std::endl;
    return;
  }

  // clear vertices from previous event
  m_container->clearParticleList();

  auto range = m_g4truthinfo->GetPrimaryParticleRange();
  for( auto iter = range.first; iter != range.second; ++iter )
  {
    auto particle = iter->second;
    if( particle )
    {
      auto particleStruct = create_particle( particle );

      // embed index
      particleStruct._embed = get_embed( particle );

      // hit mask
      const auto iter( m_g4particle_map.find( particle->get_track_id() ) );
      if( iter !=  m_g4particle_map.cend() )
      { particleStruct._mask = iter->second; }

      m_container->addParticle( particleStruct );
    }

  }

}

//_____________________________________________________________________
void SimEvaluator_hp::print_vertices()
{
  if( !m_g4truthinfo )
  {
    std::cerr << "SimEvaluator_hp::print_vertices - nodes not found." << std::endl;
    return;
  }

  // get main primary vertex id
  const auto main_vertex_id = m_g4truthinfo->GetPrimaryVertexIndex();
  const auto vertex = m_g4truthinfo->GetPrimaryVtx( main_vertex_id );
  if( vertex )
  {

    std::cout << "SimEvaluator_hp::print_vertices - main primary vertex: " << *vertex << std::endl;

  } else {

    std::cerr << "SimEvaluator_hp::print_vertices - no main primary vertex found." << std::endl;

  }

  auto range = m_g4truthinfo->GetPrimaryVtxRange();
  std::cout << "SimEvaluator_hp::print_vertices - primary vertex count: " << std::distance( range.first, range.second ) << std::endl;
  for( auto iter = range.first; iter != range.second; ++iter )
  {
    auto vertex = iter->second;
    if( vertex && vertex->get_id() != main_vertex_id )
    { std::cout << "SimEvaluator_hp::print_vertices - primary vertex: " << *vertex << std::endl; }
  }

}

//_____________________________________________________________________
int SimEvaluator_hp::get_embed( PHG4Particle* particle ) const
{ return (m_g4truthinfo && particle) ? m_g4truthinfo->isEmbeded( particle->get_primary_id() ):0; }
