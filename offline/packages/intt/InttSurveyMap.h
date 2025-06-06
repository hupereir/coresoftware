#ifndef INTT_SURVEY_MAP_H
#define INTT_SURVEY_MAP_H

#include "InttMapping.h"

#ifdef __clang__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wundefined-internal"
#include <Eigen/Geometry>
#pragma GCC diagnostic pop
#else
#include <Eigen/Geometry>
#endif

#include <cstddef>  // for size_t
#include <iostream>
#include <map>
#include <string>

class CDBTTree;

class InttSurveyMap
{
 public:
  typedef std::map<InttNameSpace::Offline_s, Eigen::Affine3d> map_t;
  typedef InttNameSpace::Offline_s key_t;
  typedef Eigen::Affine3d val_t;

  InttSurveyMap() = default;
  virtual ~InttSurveyMap();

  int LoadFromFile(std::string const& = "InttSurveyMap.root");
  int LoadFromCDB(std::string const& = "InttSurveyMap");

  int GetStripTransform(key_t const&, val_t&) const;
  int GetSensorTransform(key_t const&, val_t&) const;
  int GetLadderTransform(key_t const&, val_t&) const;

  virtual void identify(std::ostream& = std::cout) const;
  virtual std::size_t size() const;

  virtual val_t const* GetAbsoluteTransform(key_t const&) const;
  int const static Wildcard = 0xffff;

 protected:
  virtual int v_LoadFromCDBTTree(CDBTTree&);

 private:
  map_t* m_absolute_transforms = nullptr;

};

#endif  // INTT_SURVEY_MAP_H
