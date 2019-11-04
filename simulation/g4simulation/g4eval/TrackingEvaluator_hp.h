#ifndef G4EVAL_TRACKINGEVALUATOR_HP_H
#define G4EVAL_TRACKINGEVALUATOR_HP_H

#include <fun4all/SubsysReco.h>
#include <phool/PHObject.h>
#include <phool/PHTimer.h>

#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>

#include <string>
#include <memory>

class PHCompositeNode;
class SvtxEvalStack;
class TFile;
class TNtuple;

// cluster information to be stored in tree
class ClusterStruct: public TObject
{
  public:

  const char* GetName() const override
  { return "ClusterStruct"; }

  unsigned int _layer = 0;
  size_t _size = 0;
  float _x = 0;
  float _y = 0;
  float _z = 0;

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

};

#endif  // G4EVAL_TRACKINGEVALUATOR_HP_H
