#include "TrackingEvaluator_hp.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHTimer.h>

#include <iostream>

TrackingEvaluator_hp::TrackingEvaluator_hp( const std::string& name )
  : SubsysReco( name)
{
  std::cout << "TrackingEvaluator_hp::TrackingEvaluator_hp." << std::endl;
}

int TrackingEvaluator_hp::Init(PHCompositeNode* topNode)
{
  std::cout << "TrackingEvaluator_hp::Init." << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

int TrackingEvaluator_hp::InitRun(PHCompositeNode* topNode)
{
  std::cout << "TrackingEvaluator_hp::InitRun." << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

int TrackingEvaluator_hp::process_event(PHCompositeNode* topNode)
{
  if ((Verbosity() > 0) && (_ievent % 100 == 0))
  { std::cout << "TrackingEvaluator_hp::process_event - Event = " << _ievent << std::endl; }
  ++_ievent;
  return Fun4AllReturnCodes::EVENT_OK;
}

int TrackingEvaluator_hp::End(PHCompositeNode* topNode)
{
  std::cout << "TrackingEvaluator_hp::End." << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

