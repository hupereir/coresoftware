#ifndef G4EVAL_TRACKINGEVALUATOR_HP_H
#define G4EVAL_TRACKINGEVALUATOR_HP_H

#include <fun4all/SubsysReco.h>
// #include <g4main/PHG4HitDefs.h>
#include <phool/PHObject.h>
#include <phool/PHTimer.h>

#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>

#include <set>
#include <string>
#include <memory>

class PHG4Hit;
class SvtxTrack;
class SvtxTrackMap;
class TrkrClusterContainer;
class TrkrCluster;
class TrkrHitTruthAssoc;

// cluster information to be stored in tree
class ClusterStruct: public TObject
{
  public:

  const char* GetName() const override
  { return "ClusterStruct"; }

  unsigned int _layer = 0;

  // cluster position
  float _x = 0;
  float _y = 0;
  float _z = 0;
  float _r = 0;
  float _phi = 0;

  // track position
  float _trk_x = 0;
  float _trk_y = 0;
  float _trk_z = 0;
  float _trk_r = 0;
  float _trk_phi = 0;

  ClassDef(ClusterStruct,1)

};

class TrackingEvaluator_hp : public SubsysReco
{
  public:

  /// constructor
  TrackingEvaluator_hp(
    const std::string &name = "TRACKINGEVALUATOR_HP",
    const std::string &filename = "trackevaluator_hp.root" );

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

  // get geant hits associated to a cluster
  using G4HitSet = std::set<PHG4Hit*>;
  G4HitSet find_g4_hits( TrkrCluster* ) const;

  // event counter
  unsigned int _ievent = 0;
  std::unique_ptr<PHTimer> _timer;

  // filename
  std::string _filename;
  std::unique_ptr<TFile> _tfile;
  TTree* _tree = nullptr;

  // cluster array
  TClonesArray* _clusterArray = nullptr;
  int _clusterCount = 0;

  // nodes
  SvtxTrackMap* _trackMap = nullptr;
  TrkrClusterContainer* _clusterMap = nullptr;
  TrkrHitTruthAssoc* _hitTruthAssociationMap = nullptr;

};

#endif  // G4EVAL_TRACKINGEVALUATOR_HP_H
