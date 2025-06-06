// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef TRIGGER_LL1PACKETGETTER_H
#define TRIGGER_LL1PACKETGETTER_H

#include "TriggerDefs.h"

#include <fun4all/SubsysReco.h>

#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>

class PHCompositeNode;
class LL1Out;
class TriggerPrimitive;
class TriggerPrimitiveContainer;

class LL1PacketGetter : public SubsysReco
{
 public:
  explicit LL1PacketGetter(const std::string &name = "LL1PacketGetter", const std::string &trigger = "NONE", const std::string &ll1 = "NONE");
  ~LL1PacketGetter() override = default;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int CreateNodeTree(PHCompositeNode *topNode);

  void setTriggerType(const std::string &name);

  void set_dataflag(bool flag)
  {
    m_isdata = flag;
    return;
  }

 private:
  std::string m_trigger{"NONE"};
  std::string m_ll1{"NONE"};

  std::string m_triggerprimitive_nodename{"NONE"};
  std::string m_triggerprimitive_ll1_nodename{"NONE"};
  std::string m_ll1_nodename{"NONE"};

  LL1Out *m_ll1out{nullptr};
  TriggerPrimitiveContainer *m_trigger_primitives{nullptr};
  TriggerPrimitiveContainer *m_trigger_primitives_ll1{nullptr};
  TriggerPrimitive *_trigger_primitive{nullptr};
  std::map<unsigned int, std::vector<unsigned int>> *_trigger_words{nullptr};

  std::map<unsigned int, std::pair<int, int>> m_prim_sum_map{};
  std::map<unsigned int, std::pair<int, int>> m_packet_map{};
  std::map<unsigned int, int> m_word_map{};
  TriggerDefs::TriggerId m_triggerid{TriggerDefs::TriggerId::noneTId};
  unsigned int m_triggerkey{0};
  TriggerDefs::DetectorId m_detectorid{TriggerDefs::DetectorId::noneDId};
  TriggerDefs::PrimitiveId m_primitiveid{TriggerDefs::PrimitiveId::nonePId};
  std::vector<unsigned int> *_sum{nullptr};

  int m_packet_low{std::numeric_limits<int>::min()};
  int m_packet_high{std::numeric_limits<int>::min()};

  int m_nchannels{256};
  int m_nchannels_per_primitive{0};
  int m_nprimitives{0};
  int m_ntriggerwords{0};
  bool m_isdata{true};
  bool m_no_ll1out{false};
};

#endif  // TRIGGER_LL1PACKETGETTER_H
