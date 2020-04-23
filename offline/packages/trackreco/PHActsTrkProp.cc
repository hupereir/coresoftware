
#include "PHActsTrkProp.h"
#include "MakeActsGeometry.h"
#include "ActsTrack.h"

/// Fun4All includes
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/phool.h>

/// Tracking includes
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>

#include <Acts/MagneticField/InterpolatedBFieldMap.hpp>
#include <Acts/MagneticField/SharedBField.hpp>
#include <Acts/MagneticField/ConstantBField.hpp>
#include <Acts/MagneticField/MagneticFieldContext.hpp>

#include <Acts/EventData/ChargePolicy.hpp>
#include <Acts/EventData/SingleCurvilinearTrackParameters.hpp>
#include <Acts/EventData/TrackParameters.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>

#include <Acts/Propagator/EigenStepper.hpp> // this include causes seg fault when put in header file for some reason
#include <Acts/Propagator/Propagator.hpp>
#include <Acts/Propagator/Navigator.hpp>
#include <Acts/Propagator/AbortList.hpp>
#include <Acts/Propagator/ActionList.hpp>
#include <Acts/Utilities/Helpers.hpp>
#include <Acts/Utilities/Units.hpp>

#include <ACTFW/Plugins/BField/BFieldOptions.hpp>
#include <ACTFW/Plugins/BField/ScalableBField.hpp>
#include <ACTFW/Framework/ProcessCode.hpp>
#include <ACTFW/Framework/WhiteBoard.hpp>
#include <ACTFW/EventData/Track.hpp>
#include <ACTFW/Framework/AlgorithmContext.hpp>

/// Setup aliases for creating propagator
/// For some reason putting these in the header file, with appropriate headers
/// causes seg fault. Propagator also must be instantiated immediately
using ConstantBField = Acts::ConstantBField;
using Stepper = Acts::EigenStepper<ConstantBField>;
using Propagator = Acts::Propagator<Stepper, Acts::Navigator>;
Propagator *propagator;

#include <TMatrixDSym.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <utility>

PHActsTrkProp::PHActsTrkProp(const std::string& name)
  : PHTrackPropagating(name)
  , m_event(0)
  , m_actsGeometry(nullptr)
  , m_minTrackPt(0.15)
  , m_maxStepSize(3.)
  , m_trackMap(nullptr)
  , m_clusterSurfaceMap(nullptr)
  , m_clusterSurfaceTpcMap(nullptr)
  , m_actsProtoTracks(nullptr)
  , m_sourceLinks(nullptr)
{
  Verbosity(0);
}

 PHActsTrkProp::~PHActsTrkProp()
{
}

int PHActsTrkProp::Setup(PHCompositeNode* topNode)
{
  createNodes(topNode);
  
  if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
    return Fun4AllReturnCodes::ABORTEVENT;

  /// Get the magnetic field and tracking geometry to setup the Acts::Stepper
  //FW::Options::BFieldVariant bFieldVar = m_actsGeometry->magField;
  Acts::Navigator navigator(m_actsGeometry->getTGeometry());

  /// Just use the default magnetic field for now. Can access BField
  /// from m_actsGeometry using std::visit, but can't figure out how to 
  /// get necessary information out of lambda function
  ConstantBField bField (0, 0, 1.4 * Acts::UnitConstants::T);
  Stepper stepper(bField);
  propagator = new Propagator(std::move(stepper), std::move(navigator));

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsTrkProp::Process()
{

  m_event++;

  if (Verbosity() > 1)
  {
    std::cout << PHWHERE << "Events processed: " << m_event << std::endl;
    std::cout << "Start PHActsTrkProp::process_event" << std::endl;
  }


  PerigeeSurface surface = Acts::Surface::makeShared<Acts::PerigeeSurface>
    (Acts::Vector3D(0., 0., 0.));


  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsTrkProp::End()
{
  if (Verbosity() > 10)
  {
    std::cout << "Finished PHActsTrkProp" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}



void PHActsTrkProp::createNodes(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);

  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));

  if (!dstNode)
  {
    std::cerr << "DST node is missing, quitting" << std::endl;
    throw std::runtime_error("Failed to find DST node in PHActsTracks::createNodes");
  }

  PHCompositeNode *svtxNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "SVTX"));

  if (!svtxNode)
  {
    svtxNode = new PHCompositeNode("SVTX");
    dstNode->addNode(svtxNode);
  }


  m_actsProtoTracks = findNode::getClass<std::vector<ActsTrack>>(topNode, "ActsProtoTracks");

  if(!m_actsProtoTracks)
    {
      m_actsProtoTracks = new std::vector<ActsTrack>;

      PHDataNode<std::vector<ActsTrack>> *protoTrackNode =
        new PHDataNode<std::vector<ActsTrack>>(m_actsProtoTracks, "ActsProtoTracks");
      
      svtxNode->addNode(protoTrackNode);

    }


  return;
}

int PHActsTrkProp::getNodes(PHCompositeNode* topNode)
{

  m_sourceLinks = findNode::getClass<std::map<unsigned int, SourceLink>>(topNode, "TrkrClusterSourceLinks");

  if (!m_sourceLinks)
    {
      std::cout << PHWHERE << "TrkrClusterSourceLinks node not found on node tree. Exiting."
		<< std::endl;
      
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  
  m_actsGeometry = findNode::getClass<MakeActsGeometry>(topNode, "MakeActsGeometry");

  if (!m_actsGeometry)
  {
    std::cout << "ActsGeometry not on node tree. Exiting."
              << std::endl;

    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_trackMap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if (!m_trackMap)
    {
      std::cout << PHWHERE << " ERROR: Can't find SvtxTrackMap. Exiting " 
		<< std::endl;

      return Fun4AllReturnCodes::ABORTEVENT;
    }

  m_clusterSurfaceMap = findNode::getClass<std::map<TrkrDefs::hitsetkey, Surface>>(topNode, "HitSetKeySurfaceActsMap");

  if(!m_clusterSurfaceMap)
    {
      std::cout << PHWHERE <<"Can't find cluster-acts surface map Exiting."
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  m_clusterSurfaceTpcMap = findNode::getClass<std::map<TrkrDefs::cluskey, Surface>>(topNode, "ClusterTpcSurfaceActsMap");

  if(m_clusterSurfaceTpcMap)
    {
      std::cout << PHWHERE <<"Can't find TPC cluster-acts surface map Exiting."
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }


  return Fun4AllReturnCodes::EVENT_OK;
}


Acts::BoundSymMatrix PHActsTrkProp::getActsCovMatrix(const SvtxTrack *track)
{
  Acts::BoundSymMatrix matrix = Acts::BoundSymMatrix::Zero();
  const double px = track->get_px();
  const double py = track->get_py();
  const double pz = track->get_pz();
  const double p = sqrt(px * px + py * py + pz * pz);

  // Get the track seed covariance matrix
  // These are the variances, so the std devs are sqrt(seedCov[i][j])
  TMatrixDSym seedCov(6);
  for (int i = 0; i < 6; i++)
  {
    for (int j = 0; j < 6; j++)
    {
      seedCov[i][j] = track->get_error(i, j);
    }
  }

  const double sigmap = sqrt(px * px * seedCov[3][3] + py * py * seedCov[4][4] + pz * pz * seedCov[5][5]) / p;

  // Need to convert seedCov from x,y,z,px,py,pz basis to Acts basis of
  // x,y,phi/theta of p, qoverp, time
  double phi = track->get_phi();

  const double pxfracerr = seedCov[3][3] / (px * px);
  const double pyfracerr = seedCov[4][4] / (py * py);
  const double phiPrefactor = fabs(py) / (fabs(px) * (1 + (py / px) * (py / px)));
  const double sigmaPhi = phi * phiPrefactor * sqrt(pxfracerr + pyfracerr);
  const double theta = acos(pz / p);
  const double thetaPrefactor = ((fabs(pz)) / (p * sqrt(1 - (pz / p) * (pz / p))));
  const double sigmaTheta = thetaPrefactor * sqrt(sigmap * sigmap / (p * p) + seedCov[5][5] / (pz * pz));
  const double sigmaQOverP = sigmap / (p * p);

  // Just set to 0 for now?
  const double sigmaTime = 0;

  if (Verbosity() > 10)
  {
    std::cout << "Track (px,py,pz,p) = (" << px << "," << py
              << "," << pz << "," << p << ")" << std::endl;
    std::cout << "Track covariance matrix: " << std::endl;

    for (int i = 0; i < 6; i++)
    {
      for (int j = 0; j < 6; j++)
      {
        std::cout << seedCov[i][j] << ", ";
      }
      std::cout << std::endl;
    }
    std::cout << "Corresponding uncertainty calculations: " << std::endl;
    std::cout << "perr: " << sigmap << std::endl;
    std::cout << "phi: " << phi << std::endl;
    std::cout << "pxfracerr: " << pxfracerr << std::endl;
    std::cout << "pyfracerr: " << pyfracerr << std::endl;
    std::cout << "phiPrefactor: " << phiPrefactor << std::endl;
    std::cout << "sigmaPhi: " << sigmaPhi << std::endl;
    std::cout << "theta: " << theta << std::endl;
    std::cout << "thetaPrefactor: " << thetaPrefactor << std::endl;
    std::cout << "sigmaTheta: " << sigmaTheta << std::endl;
    std::cout << "sigmaQOverP: " << sigmaQOverP << std::endl;
  }

  /// Seed covariances are already variances, so don't need to square them
  matrix(Acts::eLOC_0, Acts::eLOC_0) = seedCov[0][0];
  matrix(Acts::eLOC_1, Acts::eLOC_1) = seedCov[1][1];
  matrix(Acts::ePHI, Acts::ePHI) = sigmaPhi * sigmaPhi;
  matrix(Acts::eTHETA, Acts::eTHETA) = sigmaTheta * sigmaTheta;
  matrix(Acts::eQOP, Acts::eQOP) = sigmaQOverP * sigmaQOverP;
  matrix(Acts::eT, Acts::eT) = sigmaTime;

  return matrix;
}
