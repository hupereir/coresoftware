#include "SimEvaluator_hp.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>
#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>

#include <Geant4/G4SystemOfUnits.hh>

#include <iostream>
#include <algorithm>
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

  //_____________________________________________________________________
  /// create track struct from struct from svx track
  VertexStruct create_vertex( PHG4VtxPoint* vertex )
  {
    VertexStruct vertexStruct;
    vertexStruct._x = vertex->get_x();
    vertexStruct._y = vertex->get_y();
    vertexStruct._z = vertex->get_z();
    vertexStruct._t = vertex->get_t();
    return vertexStruct;
  }

  //_____________________________________________________________________
  /// create track struct from struct from svx track
  ParticleStruct create_particle( PHG4Particle* particle )
  {
    ParticleStruct particleStruct;
    particleStruct._pid = particle->get_pid();
    particleStruct._charge = particle->get_IonCharge()/eplus;
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
{
  std::cout << "SimEvaluator_hp::SimEvaluator_hp." << std::endl;
}

//_____________________________________________________________________
int SimEvaluator_hp::Init(PHCompositeNode* topNode )
{

  std::cout << "SimEvaluator_hp::Init." << std::endl;

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

  _container = new Container;
  auto newNode = new PHIODataNode<PHObject>( _container, "SimEvaluator_hp::Container", "PHObject" );
  evalNode->addNode(newNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int SimEvaluator_hp::InitRun(PHCompositeNode* )
{
  std::cout << "SimEvaluator_hp::InitRun." << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int SimEvaluator_hp::process_event(PHCompositeNode* topNode)
{
  // load nodes
  auto res =  load_nodes(topNode);
  if( res != Fun4AllReturnCodes::EVENT_OK ) return res;

  fill_vertices();
  fill_particles();
  // print_vertices();

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int SimEvaluator_hp::End(PHCompositeNode* )
{
  std::cout << "SimEvaluator_hp::End." << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int SimEvaluator_hp::load_nodes( PHCompositeNode* topNode )
{

  // local container
  _container = findNode::getClass<Container>(topNode, "SimEvaluator_hp::Container");

  // g4 truth info
  _g4truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

  return Fun4AllReturnCodes::EVENT_OK;

}

//_____________________________________________________________________
void SimEvaluator_hp::fill_vertices()
{

  if( !( _container && _g4truthinfo ) )
  {
    std::cerr << "SimEvaluator_hp::fill_vertices - nodes not found." << std::endl;
    return;
  }

  // clear vertices from previous event
  _container->clearVertexList();

  // get main primary vertex id
  const auto main_vertex_id = _g4truthinfo->GetPrimaryVertexIndex();

  auto range = _g4truthinfo->GetPrimaryVtxRange();
  for( auto iter = range.first; iter != range.second; ++iter )
  {
    auto vertex = iter->second;
    if( vertex )
    {
      auto vertexStruct = create_vertex( vertex );
      vertexStruct._is_main_vertex = ( vertex->get_id() == main_vertex_id );
      _container->addVertex( vertexStruct );
    }

  }

}

//_____________________________________________________________________
void SimEvaluator_hp::fill_particles()
{

  if( !( _container && _g4truthinfo ) )
  {
    std::cerr << "SimEvaluator_hp::fill_vertices - nodes not found." << std::endl;
    return;
  }

  // clear vertices from previous event
  _container->clearParticleList();

  auto range = _g4truthinfo->GetPrimaryParticleRange();
  for( auto iter = range.first; iter != range.second; ++iter )
  {
    auto particle = iter->second;
    if( particle )
    {
      auto particleStruct = create_particle( particle );
      particleStruct._embed = get_embed( particle );
      _container->addParticle( particleStruct );
    }

  }

}

//_____________________________________________________________________
void SimEvaluator_hp::print_vertices()
{
  if( !_g4truthinfo )
  {
    std::cerr << "SimEvaluator_hp::print_vertices - nodes not found." << std::endl;
    return;
  }

  // get main primary vertex id
  const auto main_vertex_id = _g4truthinfo->GetPrimaryVertexIndex();
  const auto vertex = _g4truthinfo->GetPrimaryVtx( main_vertex_id );
  if( vertex )
  {

    std::cout << "SimEvaluator_hp::print_vertices - main primary vertex: " << *vertex << std::endl;

  } else {

    std::cerr << "SimEvaluator_hp::print_vertices - no main primary vertex found." << std::endl;

  }

  auto range = _g4truthinfo->GetPrimaryVtxRange();
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
{ return (_g4truthinfo && particle) ? _g4truthinfo->isEmbeded( particle->get_primary_id() ):0; }
