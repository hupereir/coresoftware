#include "SpaceChargeEvaluator_hp.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>

#include <TFile.h>
#include <THnSparse.h>
#include <TString.h>

#include <iostream>
#include <algorithm>
#include <numeric>

//_____________________________________________________________________
namespace
{

  /// square
  template<class T> T square( T x ) { return x*x; }

  /// get average between entrance and exit point
  template< float (PHG4Hit::*accessor)(int) const>
    float average( PHG4Hit* hit )
  { return 0.5*((hit->*accessor)(0) + (hit->*accessor)(1)); }

}

//_____________________________________________________________________
SpaceChargeEvaluator_hp::SpaceChargeEvaluator_hp( const std::string& name ):
  SubsysReco( name)
{
  std::cout << "SpaceChargeEvaluator_hp::SpaceChargeEvaluator_hp." << std::endl;
}

//_____________________________________________________________________
int SpaceChargeEvaluator_hp::Init(PHCompositeNode* topNode )
{
  std::cout << "SpaceChargeEvaluator_hp::Init." << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int SpaceChargeEvaluator_hp::InitRun(PHCompositeNode* )
{
  std::cout << "SpaceChargeEvaluator_hp::InitRun." << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int SpaceChargeEvaluator_hp::process_event(PHCompositeNode* topNode)
{
  // load nodes
  auto res =  load_nodes(topNode);
  if( res != Fun4AllReturnCodes::EVENT_OK ) return res;
  fill_space_charge_map();

  // increment event
  ++_ievent;

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int SpaceChargeEvaluator_hp::End(PHCompositeNode* )
{
  std::cout << "SpaceChargeEvaluator_hp::End." << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int SpaceChargeEvaluator_hp::load_nodes( PHCompositeNode* topNode )
{

  // g4hits
  _g4hits_tpc = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_TPC");
  if( !_g4hits_tpc )
  {
    std::cout << "SpaceChargeEvaluator_hp::load_nodes - could not find node G4HIT_TPC" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  return Fun4AllReturnCodes::EVENT_OK;

}

//_____________________________________________________________________
void SpaceChargeEvaluator_hp::fill_space_charge_map()
{

  if( !_g4hits_tpc ) return;

  // create histogram
  std::array<int,3> bins = {{ _z_bins, _r_bins, _phi_bins }};
  std::array<double,3> xmin = {{ 0, 20, 0 }};
  std::array<double,3> xmax = {{ 105, 78, 2*M_PI }};
  THnSparseF h( "h_sc", "h_sc", 3, &bins[0], &xmin[0], &xmax[0] );

  // define z max as the center of the last bin
  const double zmax = xmax[0] - (xmax[0]-xmin[0])/(2*_z_bins);

  // loop over g4 hits
  const auto range = _g4hits_tpc->getHits();
  for( auto iter = range.first; iter != range.second; ++iter )
  {
    // get g4 hit
    const auto hit = iter->second;

    // get r, phi and z
    const double x = square(average<&PHG4Hit::get_x>(hit));
    const double y = square(average<&PHG4Hit::get_y>(hit));
    const double z = square(average<&PHG4Hit::get_z>(hit));
    const double r = std::sqrt( square(x) + square(y) );
    const double phi = std::atan2( y, x );

    // get number of primary electrons/ions
    const double nprimary =  hit->get_eion()*1e9/_eion;

    // get number of secondary ions
    const double nsecondary = nprimary*_gain*_ibf;

    // fill histogram
    { std::array<double,3> x = {{z, r, phi}}; h.Fill( &x[0], nprimary ); }
    { std::array<double,3> x = {{zmax, r, phi}}; h.Fill( &x[0], nsecondary ); }
  }

  const TString filename = Form( "%s_%02i.root", _basefilename.c_str(), _ievent );
  std::cout << "SpaceChargeEvaluator_hp::fill_space_charge_map - writing map to " << filename << std::endl;
  TFile f( filename, "RECREATE" );
  f.cd();
  h.Write();
  f.Close();

}

