/*!
 *  \file		PHActsTrkFitter.C
 *  \brief		Refit SvtxTracks with PHActs.
 *  \details	Refit SvtxTracks with PHActs.
 *  \author	        Tony Frawley <afrawley@fsu.edu>
 */

#include "PHActsTrkFitter.h"
#include "MakeActsGeometry.h"
#include "ActsTrack.h"
#include "ActsCovarianceRotater.h"

/// Tracking includes
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <Acts/EventData/TrackParameters.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Acts/Surfaces/PlaneSurface.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/EventData/MultiTrajectory.hpp>

#include <ACTFW/EventData/Track.hpp>
#include <ACTFW/Framework/AlgorithmContext.hpp>

#include <cmath>
#include <iostream>
#include <vector>
#include <chrono>
using namespace std::chrono;

PHActsTrkFitter::PHActsTrkFitter(const std::string& name)
  : PHTrackFitting(name)
  , m_event(0)
  , m_actsProtoTracks(nullptr)
  , m_tGeometry(nullptr)
  , m_trackMap(nullptr)
  , m_timeAnalysis(false)
  , m_timeFile(nullptr)
  , h_eventTime(nullptr)
{
  Verbosity(0);
}

PHActsTrkFitter::~PHActsTrkFitter()
{
}

int PHActsTrkFitter::Setup(PHCompositeNode* topNode)
{
  if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
    return Fun4AllReturnCodes::ABORTEVENT;
  
  fitCfg.fit = FW::TrkrClusterFittingAlgorithm::makeFitterFunction(
               m_tGeometry->tGeometry,
	       m_tGeometry->magField,
	       Acts::Logging::INFO);

  if(m_timeAnalysis)
    {
      m_timeFile = new TFile("ActsTimeFile.root","RECREATE");
      h_eventTime = new TH1F("h_eventTime",";time [ms]",100,0,100);
    }		 
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsTrkFitter::Process()
{
  auto startTime = high_resolution_clock::now();
  m_event++;

  if (Verbosity() > 1)
  {
    std::cout << PHWHERE << "Events processed: " << m_event << std::endl;
    std::cout << "Start PHActsTrkFitter::process_event" << std::endl;
  }


  /// Construct a perigee surface as the target surface
  /// This surface is what Acts fits with respect to. So we put it
  /// at 0 so that the fitter is fitting with respect to the global 
  /// position. Presumably we could put this as the zvertex
  auto pSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(
	          Acts::Vector3D{0., 0., 0.});

  std::map<unsigned int, ActsTrack>::iterator trackIter;

  for (trackIter = m_actsProtoTracks->begin();
       trackIter != m_actsProtoTracks->end();
       ++trackIter)
  {
    ActsTrack track = trackIter->second;
    /// Can correlate with the SvtxTrackMap with the key
    const unsigned int trackKey = trackIter->first;

    std::vector<SourceLink> sourceLinks = track.getSourceLinks();
    FW::TrackParameters trackSeed = track.getTrackParams();
      
    if(Verbosity() > 10)
      {
	std::cout << " Processing proto track with position:" 
		  << trackSeed.position() << std::endl 
		  << "momentum: " << trackSeed.momentum() << std::endl;

      }
    /// Call KF now. Have a vector of sourceLinks corresponding to clusters
    /// associated to this track and the corresponding track seed which
    /// corresponds to the PHGenFitTrkProp track seeds
    Acts::KalmanFitterOptions<Acts::VoidOutlierFinder> kfOptions(
      m_tGeometry->geoContext,
      m_tGeometry->magFieldContext,
      m_tGeometry->calibContext,
      Acts::VoidOutlierFinder(),
      &(*pSurface));
  
    auto result = fitCfg.fit(sourceLinks, trackSeed, kfOptions);

    /// Check that the track fit result did not return an error
    if (result.ok())
    {  
      const Acts::KalmanFitterResult<SourceLink>& fitOutput = result.value();
      if (fitOutput.fittedParameters)
      {
	/// Get position, momentum from the Acts output. Update the values of
	/// the proto track
        updateSvtxTrack(fitOutput, trackKey);

        if (Verbosity() > 10)
        {
	  const auto& params = fitOutput.fittedParameters.value();
          std::cout << "Fitted parameters for track" << std::endl;
          std::cout << " position : " << params.position().transpose()
                    << std::endl;
          std::cout << " momentum : " << params.momentum().transpose()
                    << std::endl;
        }
      }
    }
    

  }

  auto stopTime = high_resolution_clock::now();
  auto eventTime = duration_cast<milliseconds>(stopTime - startTime);

  if(m_timeAnalysis)
    h_eventTime->Fill(eventTime.count());

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsTrkFitter::End(PHCompositeNode* topNode)
{
  if(m_timeAnalysis)
    {
      m_timeFile->cd();
      h_eventTime->Write();
      m_timeFile->Write();
      m_timeFile->Close();
    } 

  if (Verbosity() > 10)
  {
    std::cout << "Finished PHActsTrkFitter" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void PHActsTrkFitter::updateSvtxTrack(const Acts::KalmanFitterResult<SourceLink>& fitOutput, const unsigned int trackKey)
{
  const auto& params = fitOutput.fittedParameters.value();
 
  SvtxTrackMap::Iter trackIter = m_trackMap->find(trackKey);
  SvtxTrack *track = trackIter->second;

  /// Acts default unit is mm. So convert to cm
  track->set_x(params.position()(0) / Acts::UnitConstants::cm);
  track->set_y(params.position()(1) / Acts::UnitConstants::cm);
  track->set_z(params.position()(2) / Acts::UnitConstants::cm);
  track->set_px(params.momentum()(0));
  track->set_py(params.momentum()(1));
  track->set_pz(params.momentum()(2));

  /// This will contain the chi2/ndf
  Acts::MultiTrajectory<SourceLink> states = fitOutput.fittedStates;

  if(params.covariance())
    {
      ActsCovarianceRotater *rotater = new ActsCovarianceRotater();
      Acts::BoundSymMatrix rotatedCov = 
	rotater->rotateActsCovToSvtxTrack(fitOutput);
      
      for(int i = 0; i < 6; i++)
	{
	  for(int j = 0; j < 6; j++)
	    {
	      track->set_error(i,j, rotatedCov(i,j));
	    }
	}
    }
 
  return;

}

int PHActsTrkFitter::createNodes(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsTrkFitter::getNodes(PHCompositeNode* topNode)
{
  m_actsProtoTracks = findNode::getClass<std::map<unsigned int, ActsTrack>>(topNode, "ActsTrackMap");

  if (!m_actsProtoTracks)
  {
    std::cout << "Acts proto tracks not on node tree. Exiting."
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_tGeometry = findNode::getClass<ActsTrackingGeometry>(topNode, "ActsTrackingGeometry");
  if(!m_tGeometry)
    {
      std::cout << "ActsTrackingGeometry not on node tree. Exiting."
		<< std::endl;
      
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  m_trackMap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  
  if(!m_trackMap)
    {
      std::cout << PHWHERE << "SvtxTrackMap not found on node tree. Exiting."
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  return Fun4AllReturnCodes::EVENT_OK;
}
