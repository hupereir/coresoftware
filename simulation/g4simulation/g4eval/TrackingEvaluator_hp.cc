#include "TrackingEvaluator_hp.h"

#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase_historic/SvtxCluster.h>
#include <trackbase_historic/SvtxClusterMap.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/PHTimer.h>
#include <iostream>

ClassImp(ClusterStruct)

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
  _tree = new TTree( "Tree", "Clusters" );
  _tree->Branch( "clusterArray", "TClonesArray", &_clusterArray, kBufsize, kSplitlevel );

  // initialize timer
  _timer.reset( new PHTimer("_tracking_evaluator_hp_timer") );
  _timer->stop();

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
  // if( Verbosity>0 && (_ievent % 100 == 0))
  // if( _ievent % 100 == 0 )
  { std::cout << "TrackingEvaluator_hp::process_event - Event = " << _ievent << std::endl; }
  ++_ievent;

  // get necessary nodes
  auto trackMap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if( !trackMap ) return Fun4AllReturnCodes::ABORTEVENT;

//   // cluster map
//   auto clusterMap = findNode::getClass<SvtxClusterMap>( topNode, "SvtxClusterMap" );
//   if( !clusterMap ) return Fun4AllReturnCodes::ABORTEVENT;

  // cluster map
  auto clusterMap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if( !clusterMap ) return Fun4AllReturnCodes::ABORTEVENT;

  // clear array
  _clusterArray->Clear();
  _clusterCount = 0;

  auto range = clusterMap->getClusters();
  for( auto clusterIter = range.first; clusterIter != range.second; ++clusterIter )
  {

    const auto cluster = clusterIter->second;

    // create new object in array
    ClusterStruct clusterStruct;
    clusterStruct._x = cluster->getX();
    clusterStruct._y = cluster->getY();
    clusterStruct._z = cluster->getZ();

    // add in array
    if( _clusterArray )
    { new((*_clusterArray)[_clusterCount++]) ClusterStruct( clusterStruct ); }

  }

  // fill tree
  if( _tree )
  {
    std::cout << "TrackingEvaluator_hp::process_event - filling tree." << std::endl;
    _tree->Fill();
  }

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

  return Fun4AllReturnCodes::EVENT_OK;
}

