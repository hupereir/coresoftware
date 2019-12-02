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
class PHG4TruthInfoContainer;
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

  ///@name track momentum
  float _trk_px = 0;
  float _trk_py = 0;
  float _trk_pz = 0;
  float _trk_pt = 0;
  float _trk_p = 0;
  float _trk_eta = 0;
  //@}

  /// number of G4Hits contributing to this cluster
  unsigned int _truth_size = 0;

  ///@name truth position
  //@{
  float _truth_x = 0;
  float _truth_y = 0;
  float _truth_z = 0;
  float _truth_r = 0;
  float _truth_phi = 0;
  //@}

  ///@name truth momentum
  //@{
  float _truth_px = 0;
  float _truth_py = 0;
  float _truth_pz = 0;
  float _truth_pt = 0;
  float _truth_p = 0;
  float _truth_eta = 0;
  //@}

  ClassDef(ClusterStruct,1)

};

// cluster information to be stored in tree
class TrackStruct: public TObject
{
  public:

  const char* GetName() const override
  { return "TrackStruct"; }

  int _charge = 0;

  ///@name position
  //@{
  float _x = 0;
  float _y = 0;
  float _z = 0;
  float _r = 0;
  float _phi = 0;
  //@}

  ///@name momentum
  //@{
  float _px = 0;
  float _py = 0;
  float _pz = 0;
  float _pt = 0;
  float _p = 0;
  float _eta = 0;
  //@}

  ClassDef(TrackStruct,1)

};

/// track container
class Container: public PHObject
{

  public:

  /// constructor
  Container();

  /// copy constructor
  explicit Container(const Container &) = delete;

  /// assignment operator
  Container& operator = ( const Container& ) = delete;

  /// destructor
  ~Container() override;

  /// reset
  void Reset() override;

  ///@name accessors
  //@{
  TClonesArray* clusters() const
  { return _clusters; }

  TClonesArray* tracks() const
  { return _tracks; }

  TClonesArray* mc_tracks() const
  { return _mc_tracks; }

  //@}

  private:

  /// clusters array
  TClonesArray* _clusters = nullptr;

  /// tracks array
  TClonesArray* _tracks = nullptr;

  /// tracks array
  TClonesArray* _mc_tracks = nullptr;

  ClassDef(Container,1)

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

  // evaluate clusters
  void evaluate_mc_tracks();

  // print track content
  void print_track( SvtxTrack* ) const;

  // print cluster and association
  void print_cluster( TrkrDefs::cluskey, TrkrCluster* ) const;

  // get geant hits associated to a cluster
  using G4HitSet = std::set<PHG4Hit*>;
  G4HitSet find_g4hits( TrkrDefs::cluskey ) const;

  // cluster array
  Container* _container = nullptr;
  int _cluster_count = 0;
  int _track_count = 0;
  int _mc_track_count = 0;

  // nodes
  SvtxTrackMap* _track_map = nullptr;
  TrkrClusterContainer* _cluster_map = nullptr;
  TrkrClusterHitAssoc* _cluster_hit_map = nullptr;
  TrkrHitTruthAssoc* _hit_truth_map = nullptr;

  PHG4HitContainer* _g4hits_tpc = nullptr;
  PHG4HitContainer* _g4hits_intt = nullptr;
  PHG4HitContainer* _g4hits_mvtx = nullptr;
  PHG4HitContainer* _g4hits_outertracker = nullptr;

  PHG4TruthInfoContainer* _g4truthinfo = nullptr;

};

#endif  // G4EVAL_TRACKINGEVALUATOR_HP_H
