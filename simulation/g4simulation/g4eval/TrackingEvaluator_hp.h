#ifndef G4EVAL_TRACKINGEVALUATOR_HP_H
#define G4EVAL_TRACKINGEVALUATOR_HP_H

#include <fun4all/SubsysReco.h>
#include <phool/PHObject.h>
#include <phool/PHTimer.h>
#include <trackbase/TrkrDefs.h>

#include <TClonesArray.h>

#include <set>
#include <string>

class PHG4Hit;
class PHG4HitContainer;
class SvtxTrack;
class SvtxTrackMap;
class TrkrCluster;
class TrkrClusterContainer;
class TrkrClusterHitAssoc;
class TrkrHitTruthAssoc;

// cluster information to be stored in tree
class ClusterStruct: public TObject
{
  public:

  const char* GetName() const override
  { return "ClusterStruct"; }

  unsigned int _layer = 0;
  unsigned int _size = 0;

  ///@name cluster position
  //@{
  float _x = 0;
  float _y = 0;
  float _z = 0;
  float _r = 0;
  float _phi = 0;
  //@}

  ///@name track position
  //@{
  float _trk_x = 0;
  float _trk_y = 0;
  float _trk_z = 0;
  float _trk_r = 0;
  float _trk_phi = 0;
  //@}

  ///@name truth position
  //@{
  float _truth_x = 0;
  float _truth_y = 0;
  float _truth_z = 0;
  float _truth_r = 0;
  float _truth_phi = 0;
  //@}

  ClassDef(ClusterStruct,1)

};

/// cluster container
class ClusterContainer: public PHObject
{

  public:

  /// constructor
  ClusterContainer();

  /// copy constructor
  explicit ClusterContainer(const ClusterContainer &) = delete;

  /// assignment operator
  ClusterContainer& operator = ( const ClusterContainer& ) = delete;

  /// destructor
  ~ClusterContainer() override;

  /// reset
  void Reset() override;

  /// accessor
  TClonesArray* get() const
  { return _array; }

  private:

  /// cluster array
  TClonesArray* _array = nullptr;

  ClassDef(ClusterContainer,1)

};

class TrackingEvaluator_hp : public SubsysReco
{
  public:

  /// constructor
  TrackingEvaluator_hp( const std::string& = "TRACKINGEVALUATOR_HP" );

  /// global initialization
  int Init(PHCompositeNode*) override;

  /// run initialization
  int InitRun(PHCompositeNode*)  override;

  /// event processing
  int process_event(PHCompositeNode*) override;

  /// end of processing
  int End(PHCompositeNode*) override;

  private:

  // load nodes
  int load_nodes( PHCompositeNode* );

  // evaluate clusters
  void evaluate_clusters();

  // evaluate clusters
  void evaluate_tracks();

  // print track content
  void print_track( SvtxTrack* ) const;

  // print cluster and association
  void print_cluster( TrkrDefs::cluskey, TrkrCluster* ) const;

  // get geant hits associated to a cluster
  using G4HitSet = std::set<PHG4Hit*>;
  G4HitSet find_g4hits( TrkrDefs::cluskey ) const;

  // cluster array
  ClusterContainer* _cluster_container = nullptr;
  int _clusterCount = 0;

  // nodes
  SvtxTrackMap* _track_map = nullptr;
  TrkrClusterContainer* _cluster_map = nullptr;
  TrkrClusterHitAssoc* _cluster_hit_map = nullptr;
  TrkrHitTruthAssoc* _hit_truth_map = nullptr;

  PHG4HitContainer* _g4hits_tpc = nullptr;
  PHG4HitContainer* _g4hits_intt = nullptr;
  PHG4HitContainer* _g4hits_mvtx = nullptr;
  PHG4HitContainer* _g4hits_outertracker = nullptr;

};

#endif  // G4EVAL_TRACKINGEVALUATOR_HP_H
