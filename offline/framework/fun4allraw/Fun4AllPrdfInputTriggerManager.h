// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FUN4ALLRAW_FUN4ALLPRDFINPUTTRIGGERMANAGER_H
#define FUN4ALLRAW_FUN4ALLPRDFINPUTTRIGGERMANAGER_H

#include "InputManagerType.h"

#include <fun4all/Fun4AllInputManager.h>

#include <Event/phenixTypes.h>

#include <map>
#include <string>
#include <utility>
#include <vector>

class Event;
class SinglePrdfInput;
class CaloPacket;
class Gl1Packet;
class LL1Packet;
class PHCompositeNode;
class SingleTriggerInput;
class SyncObject;
class OfflinePacket;

class Fun4AllPrdfInputTriggerManager : public Fun4AllInputManager
{
 public:
  Fun4AllPrdfInputTriggerManager(const std::string &name = "DUMMY", const std::string &prdfnodename = "PRDF", const std::string &topnodename = "TOP");
  ~Fun4AllPrdfInputTriggerManager() override;

  int fileopen(const std::string & /* filenam */) override { return 0; }
  // cppcheck-suppress virtualCallInConstructor
  int fileclose() override;
  int run(const int nevents = 0) override;

  void Print(const std::string &what = "ALL") const override;
  int ResetEvent() override;
  int PushBackEvents(const int i) override;
  int GetSyncObject(SyncObject **mastersync) override;
  int SyncIt(const SyncObject *mastersync) override;
  int HasSyncObject() const override { return 1; }
  std::string GetString(const std::string &what) const override;
  void registerTriggerInput(SingleTriggerInput *prdfin, InputManagerType::enu_subsystem system);
  void UpdateEventFoundCounter(const int evtno);
  void UpdateDroppedPacket(const int packetid);
  void AddBeamClock(const int evtno, const int bclk, SinglePrdfInput *prdfin);
  void SetReferenceClock(const int evtno, const int bclk);
  void SetReferenceInputMgr(SinglePrdfInput *inp) { m_RefPrdfInput = inp; }
  void CreateBclkOffsets();
  uint64_t CalcDiffBclk(const uint64_t bclk1, const uint64_t bclk2);
  void DitchEvent(const int eventno);
  void Resynchronize();
  void ClearAllEvents();
  void SetPoolDepth(unsigned int d) { m_PoolDepth = d; }
  int FillCemc(const unsigned int nEvents = 2);
  int MoveCemcToNodeTree();
  void AddCemcPacket(int eventno, CaloPacket *pkt);
  int FillGl1(const unsigned int nEvents = 2);
  int MoveGl1ToNodeTree();
  void AddGl1Packet(int eventno, Gl1Packet *gl1pkt);
  int FillLL1(const unsigned int nEvents = 2);
  int MoveLL1ToNodeTree();
  void AddLL1Packet(int eventno, LL1Packet *pkt);
  int FillMbd(const unsigned int nEvents = 2);
  int MoveMbdToNodeTree();
  void AddMbdPacket(int eventno, CaloPacket *mbdpkt);
  int FillHcal(const unsigned int nEvents = 2);
  int MoveHcalToNodeTree();
  void AddHcalPacket(int eventno, CaloPacket *pkt);
  int FillZdc(const unsigned int nEvents = 2);
  int MoveZdcToNodeTree();
  void AddZdcPacket(int eventno, CaloPacket *pkt);
  // the sepd is read together with the zdc in the FillZdc method
  int MoveSEpdToNodeTree();
  void AddSEpdPacket(int eventno, CaloPacket *pkt);
  void InitialPoolDepth(unsigned int n) {m_InitialPoolDepth = n; m_PoolDepth = n;}
  void DetermineReferenceEventNumber();
  void ClockSyncCheck();

 private:
  struct SinglePrdfInputInfo
  {
    uint64_t bclkoffset{0};
  };

  struct Gl1PacketInfo
  {
    std::vector<Gl1Packet *> Gl1PacketVector;
    std::map<int, uint64_t> BcoMap;
    unsigned int EventFoundCounter{0};
  };

  struct MbdPacketInfo
  {
    std::vector<CaloPacket *> MbdPacketVector;
    unsigned int EventFoundCounter{0};
  };

  struct CemcPacketInfo
  {
    std::vector<CaloPacket *> CemcPacketVector;
    unsigned int EventFoundCounter{0};
  };

  struct HcalPacketInfo
  {
    std::vector<CaloPacket *> HcalPacketVector;
    unsigned int EventFoundCounter{0};
  };

  struct LL1PacketInfo
  {
    std::vector<LL1Packet *> LL1PacketVector;
    unsigned int EventFoundCounter{0};
  };

  struct SEpdPacketInfo
  {
    std::map<int, CaloPacket *> SEpdSinglePacketMap;
    std::map<int, uint64_t> BcoMap;
    unsigned int EventFoundCounter{0};
  };
  struct ZdcPacketInfo
  {
    std::vector<CaloPacket *> ZdcPacketVector;
    std::map<int, uint64_t> BcoMap;
    unsigned int EventFoundCounter{0};
  };

  int m_RunNumber{0};
  int m_RefEventNo{std::numeric_limits<int>::min()};
  bool m_gl1_registered_flag{false};
  bool m_mbd_registered_flag{false};
  bool m_cemc_registered_flag{false};
  bool m_hcal_registered_flag{false};
  bool m_ll1_registered_flag{false};
  bool m_zdc_registered_flag{false};
  unsigned int m_InitialPoolDepth = 10;
  unsigned int m_DefaultPoolDepth = 2;
  unsigned int m_PoolDepth = m_InitialPoolDepth;
  std::vector<SinglePrdfInput *> m_PrdfInputVector;
  std::vector<SingleTriggerInput *> m_TriggerInputVector;
  std::vector<SingleTriggerInput *> m_Gl1InputVector;
  std::vector<SingleTriggerInput *> m_CemcInputVector;
  std::vector<SingleTriggerInput *> m_HcalInputVector;
  std::vector<SingleTriggerInput *> m_LL1InputVector;
  std::vector<SingleTriggerInput *> m_MbdInputVector;
  std::vector<SingleTriggerInput *> m_SEpdInputVector;
  std::vector<SingleTriggerInput *> m_ZdcInputVector;
  SyncObject *m_SyncObject {nullptr};
  PHCompositeNode *m_topNode {nullptr};
  SinglePrdfInput *m_RefPrdfInput {nullptr};
  std::map<int, Gl1PacketInfo> m_Gl1PacketMap;
  std::map<int, MbdPacketInfo> m_MbdPacketMap;
  std::map<int, CemcPacketInfo> m_CemcPacketMap;
  std::map<int, HcalPacketInfo> m_HcalPacketMap;
  std::map<int, LL1PacketInfo> m_LL1PacketMap;
  std::map<int, SEpdPacketInfo> m_SEpdPacketMap;
  std::map<int, ZdcPacketInfo> m_ZdcPacketMap;
  std::map<int, int> m_DroppedPacketMap;
  std::map<int, std::vector<std::pair<int, SinglePrdfInput *>>> m_ClockCounters;
  std::map<int, int> m_RefClockCounters;
  std::map<SinglePrdfInput *, SinglePrdfInputInfo> m_SinglePrdfInputInfo;
};

#endif /* FUN4ALL_FUN4ALLPRDFINPUTPOOLMANAGER_H */
