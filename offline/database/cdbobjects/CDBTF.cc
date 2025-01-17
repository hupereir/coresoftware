#include "CDBTF.h"

#include <TClass.h>       // for TClass
#include <TCollection.h>  // for TIter
#include <TDirectory.h>   // for TDirectoryAtomicAdapter, TDirectory, gDirec...
#include <TFile.h>
#include <TF1.h>
#include <TKey.h>
#include <TList.h>    // for TList
#include <TObject.h>  // for TObject
#include <TROOT.h>
#include <TSystem.h>

#include <iostream>
#include <utility>  // for pair, make_pair

CDBTF::CDBTF(const std::string &fname)
  : m_Filename(fname)
{
}

CDBTF::~CDBTF()
{
  for (auto &iter : m_TFMap)
  {
    delete iter.second;
  }
  m_TFMap.clear();
}

void CDBTF::WriteCDBTF()
{
  if (m_TFMap.empty())
  {
    std::cout << __PRETTY_FUNCTION__ << " no TFs to be saved " << std::endl;
    return;
  }
  std::string currdir = gDirectory->GetPath();
  TFile *f = TFile::Open(m_Filename.c_str(), "RECREATE");
  for (auto &iter : m_TFMap)
  {
    iter.second->Write();
  }
  f->Close();
  gROOT->cd(currdir.c_str());  // restore previous directory
}

void CDBTF::LoadCalibrations()
{
  std::string currdir = gDirectory->GetPath();
  TFile *fin = TFile::Open(m_Filename.c_str());
  if (fin == nullptr)
  {
    std::cout << __PRETTY_FUNCTION__ << " Could not open " << m_Filename << std::endl;
    return;
  }
  TList *list = fin->GetListOfKeys();
  if (!list)
  {
    std::cout << __PRETTY_FUNCTION__ << " No keys found in " << m_Filename << std::endl;
    fin->Close();
    gROOT->cd(currdir.c_str());  // restore previous directory
    return;
  }
  TIter next(list);
  TKey *key;
  TObject *obj;
  TF1*t1 = nullptr;
  while ((key = (TKey *) next()))
  {
    obj = key->ReadObj();
    if ((obj->InheritsFrom("TF1")))
    {
      fin->GetObject(obj->GetName(), t1);
      //t1->SetDirectory(nullptr);
      m_TFMap.insert(std::make_pair(obj->GetName(), t1));
    }
  }
  fin->Close();
  gROOT->cd(currdir.c_str());  // restore previous directory
}

void CDBTF::Print() const
{
  for (auto &iter : m_TFMap)
  {
    std::cout << "TF " << iter.first << ", type "
              << iter.second->IsA()->GetName() << std::endl;
  }
  return;
}

void CDBTF::registerTF(TF1 *t1)
{
  const auto iter = m_TFMap.find(t1->GetName());
  if (iter != m_TFMap.end())
  {
    std::cout << __PRETTY_FUNCTION__ << " TF " << t1->GetName() << " already registered, use a different name and try again" << std::endl;
    gSystem->Exit(1);
  }
  m_TFMap.insert(std::make_pair(t1->GetName(), t1));
  return;
}

TF1 *CDBTF::getTF(const std::string &name, bool printerror)
{
  const auto iter = m_TFMap.find(name);
  if (iter == m_TFMap.end())
  {
    if (printerror)
    {
      std::cout << __PRETTY_FUNCTION__ << ": TF " << name << " not found" << std::endl;
    }
    return nullptr;
  }
  return iter->second;
}
