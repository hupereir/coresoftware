// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PHSILICONSEEDMERGER_H
#define PHSILICONSEEDMERGER_H

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>
#include <algorithm>

class PHCompositeNode;
class TrackSeedContainer;

class PHSiliconSeedMerger : public SubsysReco
{
 public:

  PHSiliconSeedMerger(const std::string &name = "PHSiliconSeedMerger");

  ~PHSiliconSeedMerger() override = default;

  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int ResetEvent(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  void trackMapName(const std::string& name) { m_trackMapName = name; }
  void clusterOverlap(const unsigned int nclusters) { m_clusterOverlap = nclusters; }
  void searchIntt() { m_mvtxOnly = false; }

 private:

  //! merge seeds based on strong matching between cluster keys
  /**
   * all cluster keys are expected to match.
   * For MVTX clusters, the strobe part of the cluster key is ignored
   */
  int strong_merger();

  //! merge seeds based on loose matching between cluster keys
  int weak_merger();

  int getNodes(PHCompositeNode *topNode);

  TrackSeedContainer *m_siliconTracks {nullptr};
  std::string m_trackMapName {"SiliconTrackSeedContainer"};
  unsigned int m_clusterOverlap {1};
  bool m_mvtxOnly {true};

  //!@name counters
  //@{
  uint64_t m_seed_counter_total = 0;
  uint64_t m_seed_counter_deleted = 0;
  //@}
};

#endif // PHSILICONSEEDMERGER_H
