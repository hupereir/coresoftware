// Tell emacs that this is a C++ source
// -*- C++ -*-.
#ifndef G4MVTX_PHG4OUTERTRACKERDETECTOR_H
#define G4MVTX_PHG4OUTERTRACKERDETECTOR_H

#include <g4main/PHG4Detector.h>

#include <array>
#include <cmath>                 // for M_PI
#include <map>
#include <set>
#include <string>
#include <tuple>                  // for tuple

class G4AssemblyVolume;
class G4LogicalVolume;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4OuterTrackerDisplayAction;
class PHG4OuterTrackerSubsystem;
class PHParametersContainer;
class PHParameters;

class PHG4OuterTrackerDetector : public PHG4Detector
{
 public:
  //! constructor
  PHG4OuterTrackerDetector(PHG4OuterTrackerSubsystem* subsys, const int layer, PHCompositeNode* Node, const PHParametersContainer* _paramsContainer, const std::string& dnam = "OTRACK");

  //! construct
  void ConstructMe(G4LogicalVolume* world) override;

  //!@name volume accessors
  //@{
  int IsInOuterTracker(G4VPhysicalVolume*) const;
  int IsBlackHole(G4VPhysicalVolume* volume) const;
  //@}

  void SuperDetector(const std::string& name) { m_SuperDetector = name; }
  const std::string SuperDetector() const { return m_SuperDetector; }
  void Detector(const std::string& name) { m_Detector = name; }
  const std::string Detector() const { return m_Detector; }

 private:
  void AddGeometryNode();
  int ConstructOuterTracker(G4LogicalVolume* sandwich);
  void SetDisplayProperty(G4AssemblyVolume* av);
  void SetDisplayProperty(G4LogicalVolume* lv);

  PHG4OuterTrackerDisplayAction* m_DisplayAction = nullptr;
  const PHParameters *params = nullptr;

  int layer = 0;
  int nseg_phi = 0;
  int nseg_z = 0;
  double inner_radius = 0;
  double outer_radius = 0;
  double length = 0;

  std::string m_Detector;
  std::string m_SuperDetector;
  std::string m_StaveGeometryFile;

  std::set<G4VPhysicalVolume *> activevols;

};

#endif
