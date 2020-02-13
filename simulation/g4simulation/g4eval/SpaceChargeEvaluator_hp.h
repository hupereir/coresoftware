#ifndef G4EVAL_SpaceChargeEvaluator_hp_H
#define G4EVAL_SpaceChargeEvaluator_hp_H

#include <fun4all/SubsysReco.h>
#include <phool/PHObject.h>
#include <phool/PHTimer.h>
#include <trackbase/TrkrDefs.h>

#include <TClonesArray.h>

#include <map>
#include <set>
#include <string>

class PHG4Hit;
class PHG4HitContainer;
class PHG4Particle;
class PHG4TruthInfoContainer;
class SvtxTrack;
class SvtxTrackMap;
class TrkrCluster;
class TrkrClusterContainer;
class TrkrClusterHitAssoc;
class TrkrHitTruthAssoc;


class SpaceChargeEvaluator_hp : public SubsysReco
{
  public:

  /// constructor
  SpaceChargeEvaluator_hp( const std::string& = "SPACECHARGEEVALUATOR_HP" );

  /// global initialization
  virtual int Init(PHCompositeNode*);

  /// run initialization
  virtual int InitRun(PHCompositeNode*);

  /// event processing
  virtual int process_event(PHCompositeNode*);

  /// end of processing
  virtual int End(PHCompositeNode*);

  ///@name modifiers
  //@{

  void set_offset( int value )
  { _offset = value; }

  void set_basefilename( const std::string& value )
  { _basefilename = value; }

  //@}

  private:

  /// load nodes
  int load_nodes( PHCompositeNode* );

  /// fill MC track map
  void fill_space_charge_map();

  // g4 hits in the TPC
  PHG4HitContainer* _g4hits_tpc = nullptr;

  // base file name
  std::string _basefilename = "spacechargemap_%05i.root";

  // map size
  unsigned int _ievent = 0;
  unsigned int _offset = 0;
  int _phi_bins = 72;
  int _r_bins = 48;
  int _z_bins = 50;
  float _eion = 38; // ev

};

#endif  // G4EVAL_SpaceChargeEvaluator_hp_H
