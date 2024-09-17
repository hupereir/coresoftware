#include "GlobalQA.h"

// Calo includes
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoDefs.h>

#include <mbd/MbdPmtContainer.h>
#include <mbd/MbdPmtHit.h>

#include <zdcinfo/ZdcReco.h>
#include <zdcinfo/Zdcinfo.h>

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/MbdVertex.h>
#include <globalvertex/MbdVertexMap.h>

#include <qautils/QAHistManagerDef.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/phool.h> // for PHWHERE

#include <ffaobjects/EventHeader.h>
#include <ffarawobjects/Gl1Packet.h>

#include <cdbobjects/CDBTTree.h> // for CDBTTree
#include <ffamodules/CDBInterface.h>

#include <TH1.h>
#include <TH2.h>
#include <TLorentzVector.h>
#include <TProfile.h>
#include <TSystem.h>

#include <boost/format.hpp>

#include <cassert>
#include <cmath> // for log10, pow, sqrt, abs, M_PI
#include <cstdint>
#include <iostream> // for operator<<, endl, basic_...
#include <limits>
#include <map> // for operator!=, _Rb_tree_con...
#include <string>
#include <utility> // for pair

GlobalQA::GlobalQA(const std::string &name)
    : SubsysReco(name), detector("HCALIN") {
  evtcount = 0;
}

GlobalQA::~GlobalQA() = default;

int GlobalQA::Init(PHCompositeNode * /*unused*/) {
  if (m_debug) {
    std::cout << "In GlobalQA::Init" << std::endl;
  }

  if (!m_overrideSEPDMapName) {
    m_sEPDMapName = "SEPD_CHANNELMAP";
  }
  if (!m_overrideSEPDFieldName) {
    m_sEPDfieldname = "epd_channel_map";
  }
  std::string calibdir = CDBInterface::instance()->getUrl(m_sEPDMapName);
  if (!calibdir.empty()) {
    cdbttree = new CDBTTree(calibdir);
  } else {
    std::cout << "GlobalQA::::InitRun No SEPD mapping file for domain "
              << m_sEPDMapName << " found" << std::endl;
    exit(1);
  }
  v.clear();
  for (int i = 0; i < 768; i++) {

    int keymap = cdbttree->GetIntValue(i, m_sEPDfieldname);
    if (keymap == 999) {
      continue;
    }

    key = TowerInfoDefs::encode_epd(keymap);
    v.push_back(key);
  }

  createHistos();

  if (m_debug) {
    std::cout << "Leaving GlobalQA::Init" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int GlobalQA::process_event(PHCompositeNode *topNode) {
  _eventcounter++;
  process_towers(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

int GlobalQA::process_towers(PHCompositeNode *topNode) {
  if (m_debug) {
    std::cout << _eventcounter << std::endl;
  }

  //--------------------------- trigger and GL1-------------------------------//
  Gl1Packet *gl1PacketInfo =
      findNode::getClass<Gl1Packet>(topNode, "GL1Packet");
  if (!gl1PacketInfo) {
    std::cout << PHWHERE << "GlobalQA::process_event: GL1Packet node is missing"
              << std::endl;
  }
  _eventcounter++;

  uint64_t triggervec = 0;
  if (gl1PacketInfo) {
    triggervec = gl1PacketInfo->getScaledVector();
    for (int i = 0; i < 64; i++) {
      bool trig_decision = ((triggervec & 0x1U) == 0x1U);

      if (trig_decision) {
        h_GlobalQA_triggerVec->Fill(i);
      }
      triggervec = (triggervec >> 1U) & 0xffffffffU;
    }
    triggervec = gl1PacketInfo->getScaledVector();
  }

  if ((triggervec >> 0xAU) & 0x1U) {
    //--------------------------- MBD vertex------------------------------//
    MbdVertexMap *mbdmap =
        findNode::getClass<MbdVertexMap>(topNode, "MbdVertexMap");
    MbdVertex *bvertex = nullptr;
    float mbd_zvtx = std::numeric_limits<float>::quiet_NaN();
    if (mbdmap) {
      for (MbdVertexMap::ConstIter mbditer = mbdmap->begin();
           mbditer != mbdmap->end(); ++mbditer) {
        bvertex = mbditer->second;
      }
      if (bvertex) {
        mbd_zvtx = bvertex->get_z();
      }
    }
    h_GlobalQA_mbd_zvtx->Fill(mbd_zvtx);
    h_GlobalQA_mbd_zvtx_wide->Fill(mbd_zvtx);
    if (!std::isfinite(mbd_zvtx)) {
      h_GlobalQA_mbd_zvtxq->SetBinContent(
          1, h_GlobalQA_mbd_zvtxq->GetBinContent(1) + 1);
    } else {
      h_GlobalQA_mbd_zvtxq->SetBinContent(
          2, h_GlobalQA_mbd_zvtxq->GetBinContent(2) + 1);
    }
  }

  if ((triggervec >> 0xAU) & 0x1U) {
    //--------------------------- sEPD ------------------------------//
    TowerInfoContainer *_sepd_towerinfo =
        findNode::getClass<TowerInfoContainer>(topNode, "TOWERS_SEPD");
    unsigned int ntowers = 0;
    if (_sepd_towerinfo) {
      ntowers = _sepd_towerinfo->size();
    }
    if (ntowers != 744) {
      std::cout << "sEPD container has unexpected size - exiting now!"
                << std::endl;
      exit(1);
    }
    float sepdsouthadcsum = 0.;
    float sepdnorthadcsum = 0.;
    if (_sepd_towerinfo) {
      for (unsigned int i = 0; i < ntowers; i++) {
        float _time =
            _sepd_towerinfo->get_tower_at_channel(i)->get_time_float();
        float _e = _sepd_towerinfo->get_tower_at_channel(i)->get_energy();
        int arm = TowerInfoDefs::get_epd_arm(v[i]);
        if (_time > 0.) {
          h_GlobalQA_sEPD_tile[i]->Fill(_e);

          if (arm == 0) {
            sepdsouthadcsum += _e;
          } else if (arm == 1) {
            if (i != 29) // problematic channel
            {
              sepdnorthadcsum += _e;
            }
          }
        }
      }

      h_GlobalQA_sEPD_adcsum_s->Fill(sepdsouthadcsum);
      h_GlobalQA_sEPD_adcsum_n->Fill(sepdnorthadcsum);
      h2_GlobalQA_sEPD_adcsum_ns->Fill(sepdsouthadcsum, sepdnorthadcsum);
    }
  }

  if ((triggervec >> 0x3U) & 0x1U) {
    // ------------------------------------- ZDC
    // -----------------------------------------//

    Zdcinfo *_zdcinfo = findNode::getClass<Zdcinfo>(topNode, "Zdcinfo");
    float totalzdcsouthcalib = 0.;
    float totalzdcnorthcalib = 0.;
    if (_zdcinfo) {
      totalzdcsouthcalib = _zdcinfo->get_zdc_energy(0);
      totalzdcnorthcalib = _zdcinfo->get_zdc_energy(1);
      zdc_zvtx = _zdcinfo->get_zvertex();
      h_GlobalQA_zdc_zvtx->Fill(zdc_zvtx);
      h_GlobalQA_zdc_zvtx_wide->Fill(zdc_zvtx);
      h_GlobalQA_zdc_energy_s->Fill(totalzdcsouthcalib);
      h_GlobalQA_zdc_energy_n->Fill(totalzdcnorthcalib);
    }
  }

  if ((triggervec >> 0xAU) & 0x1U) {
    //--------------------------- MBD ----------------------------------------//
    MbdPmtContainer *bbcpmts =
        findNode::getClass<MbdPmtContainer>(topNode, "MbdPmtContainer");
    if (!bbcpmts) {
      std::cout << "GlobalQA::process_event: Could not find MbdPmtContainer,"
                << std::endl;
      // return Fun4AllReturnCodes::ABORTEVENT;
    }

    int hits = 0;
    int hits_n = 0;
    int hits_s = 0;
    int hits_n_t = 0;
    int hits_s_t = 0;
    std::vector<float> time_sum_s = {};
    std::vector<float> time_sum_n = {};
    float sum_s = 0.;
    float sum_n = 0.;
    float sum_s2 = 0.;
    float sum_n2 = 0.;
    float tot_charge_s = 0.;
    float tot_charge_n = 0.;

    float charge_thresh = 0.4;
    if (bbcpmts) {
      int nPMTs = bbcpmts->get_npmt();
      for (int i = 0; i < nPMTs; i++) {
        MbdPmtHit *mbdpmt = bbcpmts->get_pmt(i);
        float q = mbdpmt->get_q();
        float t = mbdpmt->get_time();
        if (i < 64) {
          tot_charge_s += q;
          if (q > charge_thresh) {
            hits_s++;
          }
          if (i == 56 || isnan(t)) {
            continue;
          }
          hits_s_t++;
          time_sum_s.push_back(t);
          sum_s += t;
          sum_s2 += t * t;
        } else if (i >= 64) {
          tot_charge_n += q;
          if (q > charge_thresh) {
            hits_n++;
          }
          if (i == 120 || isnan(t)) {
            continue;
          }
          hits_n_t++;
          time_sum_n.push_back(t);
          sum_n += t;
          sum_n2 += t * t;
        }
        if (q > charge_thresh) {
          hits++;
        }
      }
    }

    // Calculating the zvtx
    std::sort(time_sum_n.begin(), time_sum_n.end());
    std::sort(time_sum_s.begin(), time_sum_s.end());
    unsigned length_s = time_sum_s.size();
    unsigned length_n = time_sum_n.size();
    float mean_north = 999;
    float mean_south = 999;
    int central_cut = 4;
    float sigma_cut = 1.5;

    if (hits_s_t >= central_cut) {
      mean_south = sum_s / static_cast<float>(hits_s_t);
      float rms_s = sqrt(sum_s2 / static_cast<float>(hits_s_t) -
                         TMath::Power(mean_south, 2));
      int nhit_s_center = 0;
      float sum_s_center = 0.;

      for (unsigned int is = 0; is < length_s; is++) {
        if (fabs(time_sum_s.at(is) - mean_south) < sigma_cut * rms_s) {
          sum_s_center += time_sum_s.at(is);
          nhit_s_center++;
        }
      }

      if (nhit_s_center > 0) {
        float mean_south_center =
            sum_s_center / static_cast<float>(nhit_s_center);
        mean_south = mean_south_center;
      }
    } else if (hits_s >= 2 && (hits_s_t >= 1)) {
      mean_south = sum_s / static_cast<float>(hits_s_t);
    }

    if (hits_n_t >= central_cut) {
      mean_north = sum_n / static_cast<float>(hits_n_t);
      float rms_n = sqrt(sum_n2 / static_cast<float>(hits_n_t) -
                         TMath::Power(mean_north, 2));
      int nhit_n_center = 0;
      float sum_n_center = 0.;

      for (unsigned int ino = 0; ino < length_n; ino++) {
        if (fabs(time_sum_n.at(ino) - mean_north) < sigma_cut * rms_n) {
          sum_n_center += time_sum_n.at(ino);
          nhit_n_center++;
        }
      }

      if (nhit_n_center > 0) {
        float mean_north_center =
            sum_n_center / static_cast<float>(nhit_n_center);
        mean_north = mean_north_center;
      }
    } else if (hits_n >= 2 && hits_n_t >= 1) {
      mean_north = sum_n / static_cast<float>(hits_n_t);
    }
    float calc_zvtx;
    if (mean_north != 999 && mean_south != 999) {
      calc_zvtx = 15 * (mean_south - mean_north);
    } else {
      calc_zvtx = 999;
    }

    h_GlobalQA_calc_zvtx->Fill(calc_zvtx);
    h_GlobalQA_calc_zvtx_wide->Fill(calc_zvtx);
    h_GlobalQA_mbd_charge_s->Fill(tot_charge_s);
    h_GlobalQA_mbd_charge_n->Fill(tot_charge_n);
    h_GlobalQA_mbd_nhit_s->Fill(hits_s);
    h_GlobalQA_mbd_nhit_n->Fill(hits_n);
    h_GlobalQA_mbd_charge_sum->Fill(tot_charge_s + tot_charge_n);
    h2_GlobalQA_mbd_charge_NS_correlation->Fill(tot_charge_s, tot_charge_n);
    h2_GlobalQA_mbd_nhits_NS_correlation->Fill(hits_s, hits_n);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int GlobalQA::End(PHCompositeNode * /*topNode*/) {

  return Fun4AllReturnCodes::EVENT_OK;
}

void GlobalQA::createHistos() {
  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  // MBD QA
  h_GlobalQA_mbd_zvtxq =
      new TH1D("h_GlobalQA_mbd_zvtxq",
               ";Scaled Trigger 10: MBD Coincidence    Has zvtx?;percentage", 2,
               -0.5, 1.5);
  h_GlobalQA_mbd_zvtx = new TH1D(
      "h_GlobalQA_mbd_zvtx", ";Scaled Trigger 10: MBD Coincidence    zvtx [cm]",
      100, -50, 50);
  h_GlobalQA_mbd_zvtx_wide = new TH1D(
      "h_GlobalQA_mbd_zvtx_wide",
      ";Scaled Trigger 10: MBD Coincidence    zvtx [cm]", 100, -300, 300);
  h_GlobalQA_calc_zvtx = new TH1D(
      "h_GlobalQA_calc_zvtx",
      ";Scaled Trigger 10: MBD Coincidence    zvtx [cm]", 100, -50, 50);
  h_GlobalQA_calc_zvtx_wide = new TH1D(
      "h_GlobalQA_calc_zvtx_wide",
      ";Scaled Trigger 10: MBD Coincidence    zvtx [cm]", 100, -300, 300);
  h_GlobalQA_mbd_charge_s =
      new TH1D("h_GlobalQA_mbd_charge_s",
               ";Scaled Trigger 10: MBD Coincidence    charge", 100, 0, 10);
  h_GlobalQA_mbd_charge_n =
      new TH1D("h_GlobalQA_mbd_charge_n",
               ";Scaled Trigger 10: MBD Coincidence    charge", 100, 0, 10);
  h_GlobalQA_mbd_nhit_s =
      new TH1D("h_GlobalQA_mbd_nhit_s",
               ";Scaled Trigger 10: MBD Coincidence    nhit", 30, -0.5, 29.5);
  h_GlobalQA_mbd_nhit_n =
      new TH1D("h_GlobalQA_mbd_nhit_n",
               ";Scaled Trigger 10: MBD Coincidence    nhit", 30, -0.5, 29.5);

  h_GlobalQA_mbd_charge_sum =
      new TH1F("h_GlobalQA_mbd_charge_sum ", " ; MBD total charge ; counts",
               100, 0., 20);

  h2_GlobalQA_mbd_charge_NS_correlation = new TH2F(
      "h2_GlobalQA_mbd_charge_NS_correlation",
      "MBD NS charge correlation ; total charge (south); total charge (north)",
      100, 0, 10, 100, 0, 10);

  h2_GlobalQA_mbd_nhits_NS_correlation =
      new TH2F("h2_GlobalQA_mbd_nhits_NS_correlation",
               "MBD NS number of hits correlation ; number of hits (south); "
               "number of hits (north)",
               70, 0., 70, 70, 0., 70);

  hm->registerHisto(h_GlobalQA_mbd_zvtx);
  hm->registerHisto(h_GlobalQA_mbd_zvtxq);
  hm->registerHisto(h_GlobalQA_mbd_zvtx_wide);
  hm->registerHisto(h_GlobalQA_calc_zvtx);
  hm->registerHisto(h_GlobalQA_calc_zvtx_wide);
  hm->registerHisto(h_GlobalQA_mbd_charge_s);
  hm->registerHisto(h_GlobalQA_mbd_charge_n);
  hm->registerHisto(h_GlobalQA_mbd_nhit_s);
  hm->registerHisto(h_GlobalQA_mbd_nhit_n);
  hm->registerHisto(h_GlobalQA_mbd_charge_sum);
  hm->registerHisto(h2_GlobalQA_mbd_charge_NS_correlation);
  hm->registerHisto(h2_GlobalQA_mbd_nhits_NS_correlation);

  // sEPD QA
  h_GlobalQA_sEPD_adcsum_s = new TH1D(
      "h_GlobalQA_sEPD_adcsum_s", " ; sEPD ADC sum ; counts", 100, -10, 50000);

  h_GlobalQA_sEPD_adcsum_n = new TH1D(
      "h_GlobalQA_sEPD_adcsum_n", " ; sEPD ADC sum ; counts", 100, -10, 50000);

  h2_GlobalQA_sEPD_adcsum_ns =
      new TH2D("h2_GlobalQA_sEPD_adcsum_ns",
               "sEPD NS ADC sum correlation ; ADC sum (south); ADC sum (north)",
               100, -10, 50000, 100, -10, 50000);

  for (int tile = 0; tile < 744; tile++) {
    h_GlobalQA_sEPD_tile[tile] = new TH1D(
        boost::str(boost::format("h_GlobalQA_sEPD_tile%d") % tile).c_str(), "",
        1000, -20000, 20000);

    hm->registerHisto(h_GlobalQA_sEPD_tile[tile]);
  }

  hm->registerHisto(h_GlobalQA_sEPD_adcsum_s);
  hm->registerHisto(h_GlobalQA_sEPD_adcsum_n);
  hm->registerHisto(h2_GlobalQA_sEPD_adcsum_ns);

  // ZDC QA
  h_GlobalQA_zdc_zvtx = new TH1D(
      "h_GlobalQA_zdc_zvtx", ";Scaled Trigger 3: ZDC Coincidence    zvtx [cm]",
      100, -1000, 1000);
  h_GlobalQA_zdc_zvtx_wide = new TH1D(
      "h_GlobalQA_zdc_zvtx_wide",
      ";Scaled Trigger 3: ZDC Coincidence    zvtx [cm]", 100, -2000, 2000);
  h_GlobalQA_zdc_energy_s = new TH1D(
      "h_GlobalQA_zdc_energy_s",
      ";Scaled Trigger 3: ZDC Coincidence    Energy [Gev]", 100, 10, 510);
  h_GlobalQA_zdc_energy_n = new TH1D(
      "h_GlobalQA_zdc_energy_n",
      ";Scaled Trigger 3: ZDC Coincidence    Energy [Gev]", 100, 10, 510);
  hm->registerHisto(h_GlobalQA_zdc_zvtx);
  hm->registerHisto(h_GlobalQA_zdc_zvtx_wide);
  hm->registerHisto(h_GlobalQA_zdc_energy_s);
  hm->registerHisto(h_GlobalQA_zdc_energy_n);

  h_GlobalQA_triggerVec = new TH1F("h_GlobalQA_triggerVec", "", 65, -0.5, 64.5);
  h_GlobalQA_triggerVec->SetDirectory(nullptr);
  hm->registerHisto(h_GlobalQA_triggerVec);
}
