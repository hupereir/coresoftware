
#include "PHSiliconSeedMerger.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/PHObject.h>
#include <phool/PHTimer.h>

#include <trackbase/TrkrDefs.h>

#include <trackbase_historic/TrackSeedContainer.h>
#include <trackbase_historic/TrackSeed.h>


#include <phool/PHCompositeNode.h>

//____________________________________________________________________________..
PHSiliconSeedMerger::PHSiliconSeedMerger(const std::string &name):
 SubsysReco(name)
{}

//____________________________________________________________________________..
PHSiliconSeedMerger::~PHSiliconSeedMerger()
{
}

//____________________________________________________________________________..
int PHSiliconSeedMerger::Init(PHCompositeNode*)
{

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHSiliconSeedMerger::InitRun(PHCompositeNode *topNode)
{
  int ret = getNodes(topNode);
  return ret;
}

//____________________________________________________________________________..
int PHSiliconSeedMerger::process_event(PHCompositeNode *)
{

  std::multimap<unsigned int, std::set<TrkrDefs::cluskey>> matches;
  std::set<unsigned int> seedsToDelete;

  if(Verbosity() > 0)
  {
    std::cout << "Silicon seed track container has " << m_siliconTracks->size() << std::endl;
  }

  for(unsigned int track1ID = 0; track1ID < m_siliconTracks->size(); ++track1ID)
  {
    auto track1 = m_siliconTracks->get(track1ID);

    if(!track1) { continue; }
    if(seedsToDelete.find(track1ID) != seedsToDelete.end()) { continue; }

    std::set<TrkrDefs::cluskey> mvtx1Keys;
    std::copy_if( track1->begin_cluster_keys(), track1->end_cluster_keys(), std::inserter(mvtx1Keys,mvtx1Keys.end()),
      [this]( const TrkrDefs::cluskey& ckey ) { return (!m_mvtxOnly) || (TrkrDefs::getTrkrId(ckey) == TrkrDefs::TrkrId::mvtxId); } );

    for(unsigned int track2ID = track1ID+1; track2ID < m_siliconTracks->size(); ++track2ID)
    {
      auto track2 = m_siliconTracks->get(track2ID);
      if(!track2) { continue; }
      if(seedsToDelete.find(track2ID) != seedsToDelete.end()) { continue; }

      std::set<TrkrDefs::cluskey> mvtx2Keys;
      std::copy_if( track2->begin_cluster_keys(), track2->end_cluster_keys(), std::inserter(mvtx2Keys,mvtx2Keys.end()),
        [this]( const TrkrDefs::cluskey& ckey ) { return (!m_mvtxOnly) || (TrkrDefs::getTrkrId(ckey) == TrkrDefs::TrkrId::mvtxId); } );

      std::vector<TrkrDefs::cluskey> intersection;
      std::set_intersection(mvtx1Keys.begin(),
        mvtx1Keys.end(),
        mvtx2Keys.begin(),
        mvtx2Keys.end(),
        std::back_inserter(intersection));

      /// If we have two clusters in common in the triplet, it is likely
      /// from the same track
      if(intersection.size() > m_clusterOverlap)
      {
        if(Verbosity() > 2)
        {
          std::cout << "Track " << track1ID << " keys " << std::endl;
          for(auto& key : mvtx1Keys)
          { std::cout << "   ckey: " << key << std::endl; }

          std::cout << "Track " << track2ID << " keys " << std::endl;
          for(auto& key : mvtx2Keys)
          { std::cout << "   ckey: " << key << std::endl; }

          std::cout << "Intersection keys " << std::endl;

          for(auto& key : intersection)
          { std::cout << "   ckey: " << key << std::endl; }
        }

	      for(auto& key : mvtx2Keys)
        {
          mvtx1Keys.insert(key);
        }

        if(Verbosity() > 2)
        {
          std::cout << "Match IDed"<<std::endl;
          for(auto& key : mvtx1Keys)
          { std::cout << "  total track keys " << key << std::endl; }
        }

        matches.insert(std::make_pair(track1ID, mvtx1Keys));
        seedsToDelete.insert(track2ID);
        // break;
      }
    }
  }

  // include missing clusters from the matched seed to the kept seed
  for(const auto& [trackKey, mvtxKeys] : matches)
  {
    auto track = m_siliconTracks->get(trackKey);
    if(Verbosity() > 2)
    { std::cout << "original track: " << std::endl; track->identify(); }

    // this is very strange:
    // with this algorithm you can end up with a seed that has two MVTX clusters from the same layer
    for(auto& key : mvtxKeys)
    {
      if(track->find_cluster_key(key) == track->end_cluster_keys())
      {
        track->insert_cluster_key(key);
        if(Verbosity() > 2)
        { std::cout << "adding " << key << std::endl; }
      }
    }
  }

  // delete all seeds marked as duplicates
  for(const auto& key : seedsToDelete)
  {
    if(Verbosity() > 2 )
    { std::cout << "Erasing track " << key << std::endl; }
    m_siliconTracks->erase(key);
  }

  if(Verbosity() > 2)
  {
    for(const auto& seed : *m_siliconTracks)
    {
      if (!seed) continue;
      seed->identify();
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHSiliconSeedMerger::ResetEvent(PHCompositeNode *)
{

  return Fun4AllReturnCodes::EVENT_OK;
}



//____________________________________________________________________________..
int PHSiliconSeedMerger::End(PHCompositeNode *)
{

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHSiliconSeedMerger::getNodes(PHCompositeNode *topNode)
{
  m_siliconTracks = findNode::getClass<TrackSeedContainer>(topNode, m_trackMapName.c_str());
  if(!m_siliconTracks)
    {
      std::cout << PHWHERE << "No silicon track container, can't merge seeds"
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

return Fun4AllReturnCodes::EVENT_OK;
}

