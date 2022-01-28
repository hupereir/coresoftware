/*!
 * \file PHG4MicromegasDetector.cc
 * \brief strongly inspired by code from Qinhua Huang <qinhua.huang@cea.fr>
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "PHG4MicromegasDetector.h"

#include <phparameter/PHParameters.h>

#include <g4detectors/PHG4CylinderGeomContainer.h>

#include <g4main/PHG4Detector.h>
#include <g4main/PHG4Subsystem.h>
#include <micromegas/CylinderGeomMicromegas.h>

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>                           // for PHNode
#include <phool/PHNodeIterator.h>                   // for PHNodeIterator
#include <phool/PHObject.h>                         // for PHObject
#include <phool/recoConsts.h>

#include <Geant4/G4Tubs.hh>
#include <Geant4/G4Box.hh>
#include <Geant4/G4Color.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4VisAttributes.hh>
#include <Geant4/G4String.hh>                       // for G4String
#include <Geant4/G4ThreeVector.hh>                  // for G4ThreeVector
#include <Geant4/G4Types.hh>                        // for G4double
#include <Geant4/G4VPhysicalVolume.hh>              // for G4VPhysicalVolume
#include <Geant4/G4VSolid.hh>                       // for G4VSolid

#include <cmath>
#include <iostream>
#include <numeric>
#include <tuple>                                    // for make_tuple, tuple
#include <utility>                                  // for pair, make_pair
#include <vector>                                   // for vector

//____________________________________________________________________________..
PHG4MicromegasDetector::PHG4MicromegasDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam)
  : PHG4Detector(subsys, Node, dnam)
  , m_Params(parameters)
{}

//_______________________________________________________________
bool PHG4MicromegasDetector::IsInDetector(G4VPhysicalVolume *volume) const
{ return m_activeVolumes.find( volume ) != m_activeVolumes.end(); }

//_______________________________________________________________
int PHG4MicromegasDetector::get_layer(G4VPhysicalVolume *volume) const
{
  const auto iter = m_activeVolumes.find( volume );
  return iter == m_activeVolumes.end() ? -1:iter->second;
}

//_______________________________________________________________
void PHG4MicromegasDetector::ConstructMe(G4LogicalVolume* logicWorld)
{
  create_materials();
  construct_micromegas(logicWorld);
  add_geometry_node();
}

//_______________________________________________________________
void PHG4MicromegasDetector::Print(const std::string &what) const
{
  std::cout << "PHG4Micromegas Detector:" << std::endl;
  if (what == "ALL" || what == "VOLUME")
  {
    std::cout << "Version 0.1" << std::endl;
    std::cout << "Parameters:" << std::endl;
    m_Params->Print();
  }
  return;
}

//_______________________________________________________________
void PHG4MicromegasDetector::create_materials() const
{
  // get the list of NIST materials
  // ---------------------------------
  auto G4_N = GetDetectorMaterial("G4_N");
  auto G4_O = GetDetectorMaterial("G4_O");
  auto G4_C = GetDetectorMaterial("G4_C");
  auto G4_H = GetDetectorMaterial("G4_H");
  auto G4_Si = GetDetectorMaterial("G4_Si");
  auto G4_Ar = GetDetectorMaterial("G4_Ar");
  auto G4_Cr = GetDetectorMaterial("G4_Cr");
  auto G4_Fe = GetDetectorMaterial("G4_Fe");
  auto G4_Mn = GetDetectorMaterial("G4_Mn");
  auto G4_Ni = GetDetectorMaterial("G4_Ni");
  auto G4_Cu = GetDetectorMaterial("G4_Cu");

  // combine elements
  // ----------------
  static const G4double temperature = 298.15*kelvin;
  static const G4double pressure = 1.*atmosphere;

  // FR4
  if (!GetDetectorMaterial("mmg_FR4", false))
  {
    auto mmg_FR4 = new G4Material( "mmg_FR4", 1.860*g/cm3, 4, kStateSolid);
    mmg_FR4->AddMaterial( G4_C,  0.43550 );
    mmg_FR4->AddMaterial( G4_H,  0.03650 );
    mmg_FR4->AddMaterial( G4_O,  0.28120 );
    mmg_FR4->AddMaterial( G4_Si, 0.24680 );
  }

  // Kapton
  if (!GetDetectorMaterial("mmg_Kapton", false))
  {
    auto mmg_Kapton = new G4Material( "mmg_Kapton", 1.420*g/cm3, 4, kStateSolid);
    mmg_Kapton->AddMaterial( G4_C, 0.6911330 );
    mmg_Kapton->AddMaterial( G4_H, 0.0263620 );
    mmg_Kapton->AddMaterial( G4_N, 0.0732700 );
    mmg_Kapton->AddMaterial( G4_O, 0.2092350);
  }

  // MMgas
  if (!GetDetectorMaterial("mmg_Gas", false))
  {
    auto mmg_Gas = new G4Material( "mmg_Gas", 0.00170335*g/cm3, 3, kStateGas, temperature, pressure);
    mmg_Gas->AddMaterial( G4_Ar, 0.900 );
    mmg_Gas->AddMaterial( G4_C,  0.0826586 );
    mmg_Gas->AddMaterial( G4_H,  0.0173414 );
  }

  // MMMesh
  if (!GetDetectorMaterial("mmg_Mesh", false))
  {
    auto mmg_Mesh = new G4Material( "mmg_Mesh", 2.8548*g/cm3, 5, kStateSolid);
    mmg_Mesh->AddMaterial( G4_Cr, 0.1900 );
    mmg_Mesh->AddMaterial( G4_Fe, 0.6800 );
    mmg_Mesh->AddMaterial( G4_Mn, 0.0200 );
    mmg_Mesh->AddMaterial( G4_Ni, 0.1000 );
    mmg_Mesh->AddMaterial( G4_Si, 0.0100 );
  }

  // MMStrips
  if (!GetDetectorMaterial("mmg_Strips", false))
  { new G4Material( "mmg_Strips", 5.248414*g/cm3, G4_Cu, kStateSolid); }

  // MMResistivePaste
  if (!GetDetectorMaterial("mmg_ResistPaste", false))
  { new G4Material( "mmg_ResistPaste", 0.77906*g/cm3, G4_C, kStateSolid); }

}

//_______________________________________________________________
void PHG4MicromegasDetector::construct_micromegas(G4LogicalVolume* logicWorld)
{
  // components enumeration
  /*
  this describes all the detector onion layers for a single side
  note that the detector is two sided
  */
  enum class Component
  {
    PCB,
    CuStrips,
    KaptonStrips,
    ResistiveStrips,
    Gas1,
    Mesh,
    Gas2,
    DriftCuElectrode,
    DriftKapton,
    DriftCarbon
  };
  
  
  // layer definition
  struct LayerStruct
  {
    // constructor
    LayerStruct( float thickness, G4Material* material, G4Colour color ):
      m_thickness( thickness ),
      m_material( material ),
      m_color( color )
    {}
    
    // thickness
    float m_thickness = 0;
    
    // material
    G4Material* m_material = nullptr;
    
    // color
    G4Colour m_color = 0;    
  };

  // define all layers
  const std::map<Component, LayerStruct> layer_map = 
  {
    { Component::PCB, {1.*mm, GetDetectorMaterial("mmg_FR4"), G4Colour::Green() }},
    { Component::CuStrips, { 12.*micrometer, GetDetectorMaterial("mmg_Strips"), G4Colour::Brown() }},
    { Component::KaptonStrips, { 50.*micrometer, GetDetectorMaterial("mmg_Kapton"), G4Colour::Brown() }},
		{ Component::ResistiveStrips, { 20.*micrometer, GetDetectorMaterial("mmg_ResistPaste" ), G4Colour::Black() }},
		{ Component::Gas1, { 120.*micrometer, GetDetectorMaterial( "mmg_Gas" ), G4Colour::Grey() }},
		{ Component::Mesh, { 18.*0.8*micrometer,  GetDetectorMaterial("mmg_Mesh"), G4Colour::White()} }, // 0.8 correction factor to thickness is to account for the mesh denstity@18/45
		{ Component::Gas2, { 3.*mm, GetDetectorMaterial( "mmg_Gas" ), G4Colour::Grey()}},
		{ Component::DriftCuElectrode, { 15.*micrometer, GetDetectorMaterial("G4_Cu"), G4Colour::Brown() }},
		{ Component::DriftKapton, { 50.*micrometer, GetDetectorMaterial("mmg_Kapton"), G4Colour::Brown() }},
		{ Component::DriftCarbon, { 1.*mm, GetDetectorMaterial("G4_C"), G4Colour(150/255., 75/255., 0) }}
  };
  
  // setup layers in the correct order, going outwards from beam axis
  /* same compoment can appear multiple times. Layer names must be unique */
  using LayerDefinition = std::tuple<Component,std::string>;
  const std::vector<LayerDefinition> layer_stack =
  {
    // inner side
    std::make_tuple( Component::DriftCarbon, "DriftCarbon_inner" ),
    std::make_tuple( Component::DriftKapton, "DriftKapton_inner" ),
    std::make_tuple( Component::DriftCuElectrode, "DriftCuElectrode_inner" ),
    std::make_tuple( Component::Gas2, "Gas2_inner" ),
    std::make_tuple( Component::Mesh, "Mesh_inner" ),
    std::make_tuple( Component::Gas1, "Gas1_inner" ),
    std::make_tuple( Component::ResistiveStrips, "ResistiveStrips_inner" ),
    std::make_tuple( Component::KaptonStrips, "KaptonStrips_inner" ),
    std::make_tuple( Component::CuStrips, "CuStrips_inner"  ),
    
    // PCB
    std::make_tuple( Component::PCB, "PCB" ),

    // outer side (= inner side, mirrored)
    std::make_tuple( Component::CuStrips, "CuStrips_outer" ),
    std::make_tuple( Component::KaptonStrips, "KaptonStrips_outer" ),
    std::make_tuple( Component::ResistiveStrips, "ResistiveStrips_outer" ),
    std::make_tuple( Component::Gas1, "Gas1_outer" ),
    std::make_tuple( Component::Mesh, "Mesh_outer" ),
    std::make_tuple( Component::Gas2, "Gas2_outer" ),
    std::make_tuple( Component::DriftCuElectrode, "DriftCuElectrode_outer" ),
    std::make_tuple( Component::DriftKapton, "DriftKapton_outer" ),
    std::make_tuple( Component::DriftCarbon, "DriftCarbon_outer" )
  };

  // start seting up volumes
  // get initial radius
  const double radius = m_Params->get_double_param("mm_radius")*cm;
  const double length = m_Params->get_double_param("mm_tilelength")*cm;
  const double width =  m_Params->get_double_param("mm_width")*cm;
  const double cyllength =  m_Params->get_double_param("mm_cyllength")*cm;
  
  //tiles coordinates x,y,z
  const int ntiles = 8;
  //first tiles along z (bottom) then each side
  const double centerPhi[ntiles] = {0.,0.,0.,0.,-M_PI/12.,-M_PI/12.,M_PI/12.,M_PI/12.};
  // 4 tiles equally distributed along tpc length then 2 on the middle on each sid << std::endl;
  const double centerZ[ntiles] = {cyllength*0.5/4., cyllength*1.5/4., cyllength*2.5/4., cyllength*3.5/4., 
                                  cyllength*1.5/4., cyllength*2.5/4., cyllength*1.5/4., cyllength*2.5/4.};  

  // get total thickness
  const double thickness = std::accumulate(
    layer_stack.begin(), layer_stack.end(), 0.,
    [&layer_map](double value, LayerDefinition layer )
    { return value + layer_map.at(std::get<0>(layer)).m_thickness; } );

  std::cout << "PHG4MicromegasDetector::ConstructMe - detector thickness is " << thickness/cm << " cm" << std::endl;

  // create mother volume
  auto cylinder_solid = new G4Tubs( G4String(GetName()+ "_m"), radius - 0.001*mm, sqrt( radius*radius + thickness*thickness/2. ) + 0.001*mm, cyllength / 2., 0, M_PI*2);

  recoConsts* rc = recoConsts::instance();
  auto cylinder_logic = new G4LogicalVolume( cylinder_solid, GetDetectorMaterial(rc->get_StringFlag("WorldMaterial")), G4String(GetName())+"_m" );
  auto vis = new G4VisAttributes(G4Color(G4Colour::Grey()));
  vis->SetForceSolid(true);
  vis->SetVisibility(false);
  //cylinder_logic->SetVisAttributes(vis);

  // add placement
  PHG4Subsystem *mysys = GetMySubsystem();
  mysys->SetLogicalVolume(cylinder_logic);
  new G4PVPlacement( nullptr, G4ThreeVector(0,0,0), cylinder_logic, G4String(GetName()), logicWorld, false, 0, OverlapCheck() );
 
 for( int idtile=0; idtile<ntiles; idtile++ )
{ 
    // make the tiles
    // ----------
    //auto tile_o = new G4Tubs(G4String(GetName())+"_tile0", radius - 0.001*mm, radius + thickness + 0.001*mm, length, -width/radius *rad, width/radius *rad);
    auto tile_o = new G4Box(G4String(GetName())+"_tile_"+ G4String(idtile), thickness / 2., width / 2., length / 2.);  
    auto tile_o_logic = new G4LogicalVolume(tile_o, GetDetectorMaterial(rc->get_StringFlag("WorldMaterial")), G4String(GetName())+"_tile0");
    tile_o_logic->SetVisAttributes(vis);
    auto zRot = new G4RotationMatrix();
    zRot->rotateY(centerPhi[idtile]);
    new G4PVPlacement( zRot, G4ThreeVector(0,0,centerZ[idtile]), tile_o_logic, G4String(GetName()) +"_tile_"+ G4String(idtile), logicWorld, false, 0, OverlapCheck() );// no rotation for tile0

    // keep track of current layer
  int layer_index = m_FirstLayer;

  // create detector
  /* we loop over registered layers and create volumes for each */
  auto current_radius = radius;
  for( const auto& layer:layer_stack )
  {
    const Component& type = std::get<0>(layer);
    const std::string& name = std::get<1>(layer);

    // layer name
    G4String cname = G4String(GetName()) + "_" + name;

    // get thickness, material and name
    const auto& thickness = layer_map.at(type).m_thickness;
    const auto& material = layer_map.at(type).m_material;
    const auto& color = layer_map.at(type).m_color;

    auto component_solid = new G4Tubs(cname+"_solid", current_radius, current_radius+thickness, length/2, 0, M_PI*2);
    auto component_logic = new G4LogicalVolume( component_solid, material, cname+"_logic");
    auto vis = new G4VisAttributes( color );
    vis->SetForceSolid(true);
    vis->SetVisibility(true);
    component_logic->SetVisAttributes(vis);

    auto component_phys = new G4PVPlacement( nullptr, G4ThreeVector(0,0,0), component_logic, cname+"_phys", cylinder_logic, false, 0, OverlapCheck() );

    // store active volume
    if( type == Component::Gas2 ) m_activeVolumes.insert( std::make_pair( component_phys, layer_index++ ) );
    else m_passiveVolumes.insert( component_phys );

      // update radius
      current_radius += thickness;
    }
  }

  // print physical layers
  std::cout << "PHG4MicromegasDetector::ConstructMe - first layer: " << m_FirstLayer << std::endl;
  for( const auto& pair:m_activeVolumes )
  {  std::cout << "PHG4MicromegasDetector::ConstructMe - layer: " << pair.second << " volume: " << pair.first->GetName() << std::endl; }

  return;
}

//_______________________________________________________________
void PHG4MicromegasDetector::add_geometry_node()
{
  // do nothing if detector is inactive
  if( !m_Params->get_int_param("active")) return;

  // find or create geometry node
  std::string geonode_name = std::string( "CYLINDERGEOM_" ) + m_SuperDetector;
  auto geonode = findNode::getClass<PHG4CylinderGeomContainer>(topNode(), geonode_name);
  if (!geonode)
  {
    geonode = new PHG4CylinderGeomContainer();
    PHNodeIterator iter(topNode());
    auto runNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "RUN"));
    auto newNode = new PHIODataNode<PHObject>(geonode, geonode_name, "PHObject");
    runNode->addNode(newNode);
  }

  // add cylinder objects
  /* one single cylinder is defined and contains all the tiles. The dimensions cover the tile corners. */
  for( const auto& pair:m_activeVolumes )
  {
    // store layer and volume
    const int layer = pair.second;
    const G4VPhysicalVolume* volume_phys = pair.first;

    // get solid volume, cast to a box
    const auto box = dynamic_cast<const G4Box*>( volume_phys->GetLogicalVolume()->GetSolid() );

    // create cylinder and match geometry
    /* note: cylinder segmentation type and pitch is set in PHG4MicromegasHitReco */
    auto cylinder = new CylinderGeomMicromegas(layer);
    const double innerRadius = m_Params->get_double_param("mm_radius")*cm;
    const double outerRadius = sqrt( innerRadius/cm*innerRadius/cm + box->GetYHalfLength()/cm*box->GetYHalfLength()/cm );
    cylinder->set_radius( ( innerRadius/cm + outerRadius/cm )/2. );
    cylinder->set_thickness( box->GetXHalfLength()*2./cm );
    cylinder->set_zmin( -box->GetZHalfLength()/cm );
    cylinder->set_zmax( box->GetZHalfLength()/cm );
    std::cout << "PHG4MicromegasDetector:: added geom node ok" <<std::endl;
    geonode->AddLayerGeom(layer, cylinder);
  }

}
