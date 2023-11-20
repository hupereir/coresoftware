#ifndef MBD_MBDCALIB_H
#define MBD_MBDCALIB_H

#include "MbdDefs.h"

#include <fun4all/Fun4AllBase.h>

#include <phool/recoConsts.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>

class TTree;
class CDBInterface;

class MbdCalib : public Fun4AllBase
{
 public:
  MbdCalib();

  // MbdCalib(MbdCalib &other) = delete;
  // void operator=(const MbdCalib &) = delete;

  virtual ~MbdCalib() {}

  float get_qgain(const int ipmt) const { return _qfit_mpv[ipmt]; }
  float get_tq0(const int ipmt) const { return _tqfit_t0mean[ipmt]; }
  int get_sampmax(const int ifeech) const { return _sampmax[ifeech]; }

  int Download_Gains(const std::string& dbfile);
  int Download_TQT0(const std::string& dbfile);
  int Download_SampMax(const std::string& dbfile);

  int Download_All();

  int StoreInDatabase();

  // void Dump_to_file(const std::string& what = "ALL");

  void Reset();
  // void Print(Option_t* option) const;

 private:
  CDBInterface* _cdb{nullptr};
  recoConsts* _rc{nullptr};

  int _status{0};
  // int          _run_number {0};
  // uint64_t     _timestamp {0};
  std::string _dbfilename;

  // Assumes Landau fit
  std::array<float, MbdDefs::BBC_N_PMT> _qfit_integ{};
  std::array<float, MbdDefs::BBC_N_PMT> _qfit_mpv{};
  std::array<float, MbdDefs::BBC_N_PMT> _qfit_sigma{};
  std::array<float, MbdDefs::BBC_N_PMT> _qfit_integerr{};
  std::array<float, MbdDefs::BBC_N_PMT> _qfit_mpverr{};
  std::array<float, MbdDefs::BBC_N_PMT> _qfit_sigmaerr{};
  std::array<float, MbdDefs::BBC_N_PMT> _qfit_chi2ndf{};

  // T0 offsets, charge channels
  std::array<float, MbdDefs::BBC_N_PMT> _tqfit_t0mean{};
  std::array<float, MbdDefs::BBC_N_PMT> _tqfit_t0meanerr{};
  std::array<float, MbdDefs::BBC_N_PMT> _tqfit_t0sigma{};
  std::array<float, MbdDefs::BBC_N_PMT> _tqfit_t0sigmaerr{};

  // Peak of waveform
  std::array<int, MbdDefs::BBC_N_FEECH> _sampmax{};
};

#endif  // MBD_MBDCALIB_H
