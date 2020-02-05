#include "PHG4EICForwardEcalDetector.h"

#include "PHG4ForwardEcalDetector.h"  // for PHG4ForwardEcalDetector
#include "PHG4ForwardEcalDisplayAction.h"

#include <phparameter/PHParameters.h>

#include <Geant4/G4Box.hh>
#include <Geant4/G4Cons.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4RotationMatrix.hh>  // for G4RotationMatrix
#include <Geant4/G4String.hh>          // for G4String
#include <Geant4/G4SystemOfUnits.hh>   // for cm, mm
#include <Geant4/G4ThreeVector.hh>     // for G4ThreeVector
#include <Geant4/G4Transform3D.hh>     // for G4Transform3D
#include <Geant4/G4Types.hh>           // for G4double, G4int

#include <TSystem.h>

#include <cmath>                          // for M_PI
#include <fstream>
#include <iostream>
#include <sstream>
#include <utility>  // for pair, make_pair

class G4VSolid;
class PHCompositeNode;

using namespace std;

//_______________________________________________________________________
PHG4EICForwardEcalDetector::PHG4EICForwardEcalDetector(PHG4Subsystem* subsys, PHCompositeNode* Node, PHParameters* parameters, const std::string& dnam)
  : PHG4ForwardEcalDetector(subsys, Node, parameters, dnam)
  , _tower_dx(30 * mm)
  , _tower_dy(30 * mm)
  , _tower_dz(170.0 * mm)
  , _materialScintillator("G4_POLYSTYRENE")
  , _materialAbsorber("G4_Pb")
{
}

//_______________________________________________________________________
void PHG4EICForwardEcalDetector::ConstructMe(G4LogicalVolume* logicWorld)
{
  if (Verbosity() > 0)
  {
    cout << "PHG4EICForwardEcalDetector: Begin Construction" << endl;
  }

  /* Read parameters for detector construction and mappign from file */
  ParseParametersFromTable();

  /* Create the cone envelope = 'world volume' for the crystal calorimeter */
  G4Material* Air = G4Material::GetMaterial("G4_AIR");

  G4VSolid* ecal_envelope_solid = new G4Cons("hEcal_envelope_solid",
                                             GetRMin(0), GetRMax(0),
                                             GetRMin(1), GetRMax(1),
                                             GetdZ() / 2.,
                                             0, 2 * M_PI);

  G4LogicalVolume* ecal_envelope_log = new G4LogicalVolume(ecal_envelope_solid, Air, G4String("hEcal_envelope"), 0, 0, 0);

  /* Define visualization attributes for envelope cone */
  GetDisplayAction()->AddVolume(ecal_envelope_log, "Envelope");

  /* Define rotation attributes for envelope cone */
  G4RotationMatrix ecal_rotm;
  ecal_rotm.rotateX(GetXRot());
  ecal_rotm.rotateY(GetYRot());
  ecal_rotm.rotateZ(GetZRot());

  /* Place envelope cone in simulation */
  ostringstream name_envelope;
  name_envelope.str("");
  name_envelope << TowerLogicNamePrefix() << "_envelope";

  new G4PVPlacement(G4Transform3D(ecal_rotm, G4ThreeVector(GetPlaceX(), GetPlaceY(), GetPlaceZ())),
                    ecal_envelope_log, name_envelope.str().c_str(), logicWorld, 0, false, OverlapCheck());

  /* Construct single calorimeter tower */
  G4LogicalVolume* singletower = ConstructTower();

  /* Place calorimeter tower within envelope */
  PlaceTower(ecal_envelope_log, singletower);

  return;
}

//_______________________________________________________________________
G4LogicalVolume*
PHG4EICForwardEcalDetector::ConstructTower()
{
  if (Verbosity() > 0)
  {
    cout << "PHG4EICForwardEcalDetector: Build logical volume for single tower..." << endl;
  }

  /* create logical volume for single tower */
  G4Material* material_air = G4Material::GetMaterial("G4_AIR");

  G4VSolid* single_tower_solid = new G4Box(G4String("single_tower_solid"),
                                           _tower_dx / 2.0,
                                           _tower_dy / 2.0,
                                           _tower_dz / 2.0);

  G4LogicalVolume* single_tower_logic = new G4LogicalVolume(single_tower_solid,
                                                            material_air,
                                                            "single_tower_logic",
                                                            0, 0, 0);

  /* create geometry volumes for scintillator and absorber plates to place inside single_tower */
  G4int nlayers = 60;
  G4double thickness_layer = _tower_dz / (float) nlayers;
  G4double thickness_absorber = thickness_layer / 3.0 * 2.0;      // 2/3rd absorber
  G4double thickness_scintillator = thickness_layer / 3.0 * 1.0;  // 1/3rd scintillator
  // PHENIX EMCal JGL 12/27/2015
  // G4int nlayers = 66;
  // G4double thickness_layer = _tower_dz/(float)nlayers;
  // G4double thickness_absorber = thickness_layer*0.23; // 27% absorber by length
  // G4double thickness_scintillator = thickness_layer*0.73; // 73% scintillator by length

  G4VSolid* solid_absorber = new G4Box(G4String("single_plate_absorber_solid"),
                                       _tower_dx / 2.0,
                                       _tower_dy / 2.0,
                                       thickness_absorber / 2.0);

  G4VSolid* solid_scintillator = new G4Box(G4String("single_plate_scintillator"),
                                           _tower_dx / 2.0,
                                           _tower_dy / 2.0,
                                           thickness_scintillator / 2.0);

  /* create logical volumes for scintillator and absorber plates to place inside single_tower */
  G4Material* material_scintillator = G4Material::GetMaterial(_materialScintillator.c_str());
  G4Material* material_absorber = G4Material::GetMaterial(_materialAbsorber.c_str());

  G4LogicalVolume* logic_absorber = new G4LogicalVolume(solid_absorber,
                                                        material_absorber,
                                                        "single_plate_absorber_logic",
                                                        0, 0, 0);

  G4LogicalVolume* logic_scint = new G4LogicalVolume(solid_scintillator,
                                                     material_scintillator,
                                                     "hEcal_scintillator_plate_logic",
                                                     0, 0, 0);
  AbsorberLogicalVolSetInsert(logic_absorber);
  ScintiLogicalVolSetInsert(logic_scint);
  GetDisplayAction()->AddVolume(logic_absorber, "Absorber");
  GetDisplayAction()->AddVolume(logic_scint, "Scintillator");

  /* place physical volumes for absorber and scintillator plates */
  G4double xpos_i = 0;
  G4double ypos_i = 0;
  G4double zpos_i = (-1 * _tower_dz / 2.0) + thickness_absorber / 2.0;

  ostringstream name_absorber;
  name_absorber.str("");
  name_absorber << TowerLogicNamePrefix() << "_single_plate_absorber";

  ostringstream name_scintillator;
  name_scintillator.str("");
  name_scintillator << TowerLogicNamePrefix() << "_single_plate_scintillator";

  for (int i = 1; i <= nlayers; i++)
  {
    new G4PVPlacement(0, G4ThreeVector(xpos_i, ypos_i, zpos_i),
                      logic_absorber,
                      name_absorber.str().c_str(),
                      single_tower_logic,
                      0, 0, OverlapCheck());

    zpos_i += (thickness_absorber / 2. + thickness_scintillator / 2.);

    new G4PVPlacement(0, G4ThreeVector(xpos_i, ypos_i, zpos_i),
                      logic_scint,
                      name_scintillator.str().c_str(),
                      single_tower_logic,
                      0, 0, OverlapCheck());

    zpos_i += (thickness_absorber / 2. + thickness_scintillator / 2.);
  }

  GetDisplayAction()->AddVolume(single_tower_logic, "ScintillatorSingleTower");

  if (Verbosity() > 0)
  {
    cout << "PHG4EICForwardEcalDetector: Building logical volume for single tower done." << endl;
  }

  return single_tower_logic;
}

int PHG4EICForwardEcalDetector::PlaceTower(G4LogicalVolume* ecalenvelope, G4LogicalVolume* singletower)
{
  /* Loop over all tower positions in vector and place tower */
  typedef std::map<std::string, towerposition>::iterator it_type;

  for (it_type iterator = _map_tower.begin(); iterator != _map_tower.end(); ++iterator)
  {
    if (Verbosity() > 0)
    {
      cout << "PHG4EICForwardEcalDetector: Place tower " << iterator->first
           << " at x = " << iterator->second.x << " , y = " << iterator->second.y << " , z = " << iterator->second.z << endl;
    }

    new G4PVPlacement(0, G4ThreeVector(iterator->second.x, iterator->second.y, iterator->second.z),
                      singletower,
                      iterator->first.c_str(),
                      ecalenvelope,
                      0, 0, OverlapCheck());
  }

  return 0;
}

int PHG4EICForwardEcalDetector::ParseParametersFromTable()
{
  /* Open the datafile, if it won't open return an error */
  ifstream istream_mapping;
  istream_mapping.open(GetParams()->get_string_param("mapping_file"));
  if (!istream_mapping.is_open())
  {
    cout << "ERROR in PHG4EICForwardEcalDetector: Failed to open mapping file " << GetParams()->get_string_param("mapping_file") << endl;
    gSystem->Exit(1);
  }

  /* loop over lines in file */
  string line_mapping;
  while (getline(istream_mapping, line_mapping))
  {
    /* Skip lines starting with / including a '#' */
    if (line_mapping.find("#") != string::npos)
    {
      if (Verbosity() > 0)
      {
        cout << "PHG4EICForwardEcalDetector: SKIPPING line in mapping file: " << line_mapping << endl;
      }
      continue;
    }

    istringstream iss(line_mapping);

    /* If line starts with keyword Tower, add to tower positions */
    if (line_mapping.find("Tower ") != string::npos)
    {
      unsigned idx_j, idx_k, idx_l;
      G4double pos_x, pos_y, pos_z;
      G4double size_x, size_y, size_z;
      G4double rot_x, rot_y, rot_z;
      G4double dummy;
      string dummys;

      /* read string- break if error */
      if (!(iss >> dummys >> dummy >> idx_j >> idx_k >> idx_l >> pos_x >> pos_y >> pos_z >> size_x >> size_y >> size_z >> rot_x >> rot_y >> rot_z))
      {
        cerr << "ERROR in PHG4EICForwardEcalDetector: Failed to read line in mapping file " << GetParams()->get_string_param("mapping_file") << endl;
        gSystem->Exit(1);
      }

      /* Construct unique name for tower */
      /* Mapping file uses cm, this class uses mm for length */
      ostringstream towername;
      towername.str("");
      towername << TowerLogicNamePrefix() << "_j_" << idx_j << "_k_" << idx_k;

      /* Add Geant4 units */
      pos_x = pos_x * cm;
      pos_y = pos_y * cm;
      pos_z = pos_z * cm;

      /* insert tower into tower map */
      towerposition tower_new;
      tower_new.x = pos_x;
      tower_new.y = pos_y;
      tower_new.z = pos_z;
      _map_tower.insert(make_pair(towername.str(), tower_new));
    }
    else
    {
      /* If this line is not a comment and not a tower, save parameter as string / value. */
      string parname;
      G4double parval;

      /* read string- break if error */
      if (!(iss >> parname >> parval))
      {
        cerr << "ERROR in PHG4EICForwardEcalDetector: Failed to read line in mapping file " << GetParams()->get_string_param("mapping_file") << endl;
        gSystem->Exit(1);
      }
      InsertParam(parname, parval);
    }
  }

  /* Update member variables for global parameters based on parsed parameter file */
  std::map<string, double>::const_iterator parit;

  parit = FindIter("Gtower_dx");
  if (parit != EndIter())
    _tower_dx = parit->second * cm;

  parit = FindIter("Gtower_dy");
  if (parit != EndIter())
    _tower_dy = parit->second * cm;

  parit = FindIter("Gtower_dz");
  if (parit != EndIter())
    _tower_dz = parit->second * cm;

  ostringstream rad;
  for (int i = 0; i < 2; i++)
  {
    int index = i + 1;
    rad.str("");
    rad << "Gr" << index << "_inner";
    parit = FindIter(rad.str());
    if (parit != EndIter())
    {
      SetRMin(i, parit->second * cm);
    }
    rad.str("");
    rad << "Gr" << index << "_outer";
    parit = FindIter(rad.str());
    if (parit != EndIter())
    {
      SetRMax(i, parit->second * cm);
    }
  }
  parit = FindIter("Gdz");
  if (parit != EndIter())
  {
    SetdZ(parit->second * cm);
  }
  parit = FindIter("Gx0");
  if (parit != EndIter())
  {
    SetPlaceX(parit->second * cm);
  }

  parit = FindIter("Gy0");
  if (parit != EndIter())
  {
    SetPlaceY(parit->second * cm);
  }
  parit = FindIter("Gz0");
  if (parit != EndIter())
  {
    SetPlaceZ(parit->second * cm);
  }
  parit = FindIter("Grot_x");
  if (parit != EndIter())
    SetXRot(parit->second);

  parit = FindIter("Grot_y");
  if (parit != EndIter())
    SetYRot(parit->second);

  parit = FindIter("Grot_z");
  if (parit != EndIter())
    SetZRot(parit->second);

  return 0;
}
