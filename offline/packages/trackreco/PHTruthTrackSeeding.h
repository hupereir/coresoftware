/*!
 *  \file		PHTruthTrackSeeding.h
 *  \brief		Vertexing using truth info
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#ifndef TRACKRECO_PHTRUTHTRACKSEEDING_H
#define TRACKRECO_PHTRUTHTRACKSEEDING_H

#include "PHTrackSeeding.h"

#include <set>
#include <string>  // for string

// forward declarations
class PHCompositeNode;
class PHG4TruthInfoContainer;
class PHG4HitContainer;
class TrkrClusterHitAssoc;
class TrkrHitTruthAssoc;

//class SvtxHitMap;
//class PHG4CellContainer;

/// \class PHTruthTrackSeeding
///
/// \brief Vertexing using truth info
///

class PHTruthTrackSeeding : public PHTrackSeeding
{
 public:

  PHTruthTrackSeeding(const std::string& name = "PHTruthTrackSeeding");

  unsigned int get_min_clusters_per_track() const
  { return _min_clusters_per_track; }

  void set_min_clusters_per_track(unsigned int minClustersPerTrack)
  { _min_clusters_per_track = minClustersPerTrack; }

  const std::set<unsigned int>& get_seeding_layers() const
  { return _seeding_layers; }

  void set_seeding_layers(const unsigned int a[], const unsigned int n)
  {
    _seeding_layers.clear();
    for (unsigned int i = 0; i < n; ++i) _seeding_layers.insert(a[i]);
  }

 protected:
  int Setup(PHCompositeNode* topNode) override;

  int Process(PHCompositeNode* topNode) override;

  int End() override;

 private:
  /// fetch node pointers
  int GetNodes(PHCompositeNode* topNode);

  PHG4TruthInfoContainer* _g4truth_container = nullptr;

  PHG4HitContainer* phg4hits_tpc = nullptr;
  PHG4HitContainer* phg4hits_intt = nullptr;
  PHG4HitContainer* phg4hits_mvtx = nullptr;

  TrkrHitTruthAssoc* hittruthassoc = nullptr;
  TrkrClusterHitAssoc* clusterhitassoc = nullptr;

  //SvtxHitMap* hitsmap;
  //PHG4CellContainer* cells_svtx = nullptr;
  //PHG4CellContainer* cells_intt = nullptr;
  //PHG4CellContainer* cells_maps = nullptr;

  /// seeding layers
  std::set<unsigned int> _seeding_layers;

  unsigned int _min_clusters_per_track = 0;

};

#endif
