#ifndef G4EVAL_TRACKINGEVALUATOR_HP_H
#define G4EVAL_TRACKINGEVALUATOR_HP_H

#include <fun4all/SubsysReco.h>

#include <string>
#include <memory>

class PHCompositeNode;
class PHTimer;
class SvtxEvalStack;
class TFile;
class TNtuple;

class TrackingEvaluator_hp : public SubsysReco
{
  public:

  /// constructor
  TrackingEvaluator_hp( const std::string &name = "TRACKINGEVALUATOR_HP" );

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

};

#endif  // G4EVAL_TRACKINGEVALUATOR_HP_H
