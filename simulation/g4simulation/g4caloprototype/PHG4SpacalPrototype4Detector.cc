#include "PHG4SpacalPrototype4Detector.h"

#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <g4detectors/PHG4CylinderGeom_Spacalv1.h>  // for PHG4CylinderGeom_Spacalv1:...

#include <phparameter/PHParameters.h>

#include <g4main/PHG4Utils.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>

#include <Geant4/G4Box.hh>
#include <Geant4/G4DisplacedSolid.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4PhysicalConstants.hh>
#include <Geant4/G4String.hh>  // for G4String
#include <Geant4/G4SubtractionSolid.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>  // for G4ThreeVector
#include <Geant4/G4Transform3D.hh>  // for G4Transform3D, G4Translate3D
#include <Geant4/G4Trap.hh>
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4Types.hh>  // for G4double
#include <Geant4/G4UserLimits.hh>
#include <Geant4/G4Vector3D.hh>
#include <Geant4/G4VisAttributes.hh>

#include <boost/foreach.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>  // for operator<<, basic_ostream
#include <numeric>   // std::accumulate
#include <sstream>
#include <string>  // std::string, std::to_string

class PHG4CylinderGeom;

using namespace std;

//_______________________________________________________________
//note this inactive thickness is ~1.5% of a radiation length
PHG4SpacalPrototype4Detector::PHG4SpacalPrototype4Detector(PHCompositeNode* Node, PHParameters* parameters, const std::string& dnam)
  : PHG4Detector(Node, dnam)
  , construction_params(parameters)
  , cylinder_solid(nullptr)
  , cylinder_logic(nullptr)
  , cylinder_physi(nullptr)
  ,  //
  active(0)
  , absorberactive(0)
  ,  //
  step_limits(nullptr)
  , clading_step_limits(nullptr)
  , fiber_core_step_limits(nullptr)
  ,  //
  _geom(nullptr)
{
  SetActive(construction_params->get_int_param("active"));
  SetAbsorberActive(construction_params->get_int_param("absorberactive"));
  SuperDetector(dnam);
}

PHG4SpacalPrototype4Detector::~PHG4SpacalPrototype4Detector(void)
{
  // deleting nullptr pointers is allowed (results in NOOP)
  // so checking for not null before deleting is not needed
  delete step_limits;
  delete clading_step_limits;
  delete fiber_core_step_limits;
  delete _geom;
}

//_______________________________________________________________
int PHG4SpacalPrototype4Detector::IsInCylinderActive(
    const G4VPhysicalVolume* volume)
{
  //  cout << "checking detector" << endl;
  if (active && fiber_core_vol.find(volume) != fiber_core_vol.end())
  {
    //      return fiber_core_vol.find(volume)->second;
    return FIBER_CORE;
  }
  if (absorberactive)
  {
    if (fiber_vol.find(volume) != fiber_vol.end())
      return FIBER_CLADING;

    if (block_vol.find(volume) != block_vol.end())
      return ABSORBER;

    if (calo_vol.find(volume) != calo_vol.end())
      return SUPPORT;
  }
  return INACTIVE;
}

//_______________________________________________________________
void PHG4SpacalPrototype4Detector::Construct(G4LogicalVolume* logicWorld)
{
  PHNodeIterator iter(topNode());
  PHCompositeNode* parNode = dynamic_cast<PHCompositeNode*>(iter.findFirst(
      "PHCompositeNode", "RUN"));
  assert(parNode);

  if (!_geom)
  {
    _geom = new SpacalGeom_t();
  }
  assert(_geom);

  _geom->set_config(
      SpacalGeom_t::kFullProjective_2DTaper_Tilted_SameLengthFiberPerTower);
  //  _geom->load_demo_sector_tower_map4();
  //    _geom->set_construction_verbose(2);
  _geom->ImportParameters(*construction_params);

  _geom->subtower_consistency_check();

  //  step_limits = new G4UserLimits(_geom->get_calo_step_size() * cm);
  //
  //  clading_step_limits = new G4UserLimits(
  //      _geom->get_fiber_clading_step_size() * cm);
  step_limits = nullptr;
  clading_step_limits = nullptr;

  fiber_core_step_limits = new G4UserLimits(
      _geom->get_fiber_core_step_size() * cm);

  const G4double z_rotation = construction_params->get_double_param(
                                  "z_rotation_degree") *
                              degree;
  const G4double enclosure_x = construction_params->get_double_param(
                                   "enclosure_x") *
                               cm;
  const G4double enclosure_x_shift = construction_params->get_double_param(
                                         "enclosure_x_shift") *
                                     cm;
  const G4double enclosure_y = construction_params->get_double_param(
                                   "enclosure_y") *
                               cm;
  const G4double enclosure_z = construction_params->get_double_param(
                                   "enclosure_z") *
                               cm;

  const G4double box_x_shift = (_geom->get_radius() + 0.5 * _geom->get_thickness()) * cm + enclosure_x_shift;
  const G4double box_z_shift = 0.5 * (_geom->get_zmin() + _geom->get_zmax()) * cm;

  //  G4Tubs* _cylinder_solid = new G4Tubs(G4String(GetName()),
  //      _geom->get_radius() * cm, _geom->get_max_radius() * cm,
  //      _geom->get_length() * cm / 2.0, 0, twopi);
  //  G4VSolid* cylinder_solid = new G4Box(G4String(GetName()),
  //      _geom->get_thickness() * 1.1 * 0.5 * cm, //
  //      twopi / _geom->get_azimuthal_n_sec() * _geom->get_sector_map().size()
  //          * 0.5 * _geom->get_max_radius() * cm, //
  //      _geom->get_length() * cm / 2.0);
  G4VSolid* cylinder_solid = new G4Box(G4String(GetName()),
                                       enclosure_x * 0.5,  //
                                       enclosure_y * 0.5,  //
                                       enclosure_z * 0.5);

  cylinder_solid = new G4DisplacedSolid(G4String(GetName() + "_displaced"),
                                        cylinder_solid, 0,  //
                                        G4ThreeVector(box_x_shift, 0, box_z_shift));

  G4Material* cylinder_mat = G4Material::GetMaterial("G4_AIR");
  assert(cylinder_mat);

  cylinder_logic = new G4LogicalVolume(cylinder_solid, cylinder_mat,
                                       G4String(GetName()), 0, 0, 0);
  G4VisAttributes* VisAtt = new G4VisAttributes();
  VisAtt->SetColour(G4Colour::White());
  VisAtt->SetVisibility(true);
  VisAtt->SetForceSolid(false);
  VisAtt->SetForceWireframe(true);
  cylinder_logic->SetVisAttributes(VisAtt);

  G4Transform3D cylinder_place(
      G4Translate3D(_geom->get_xpos() * cm, _geom->get_ypos() * cm,
                    _geom->get_zpos() * cm)  //
      * G4RotateZ3D(z_rotation)              //
      * G4Translate3D(
            -(_geom->get_radius() + 0.5 * _geom->get_thickness()) * cm, 0,
            0));

  cylinder_physi = new G4PVPlacement(cylinder_place, cylinder_logic,
                                     G4String(GetName()), logicWorld, false, 0, OverlapCheck());

  // install sectors
  std::pair<G4LogicalVolume*, G4Transform3D> psec = Construct_AzimuthalSeg();
  G4LogicalVolume* sec_logic = psec.first;
  const G4Transform3D& sec_trans = psec.second;
  BOOST_FOREACH (const SpacalGeom_t::sector_map_t::value_type& val, _geom->get_sector_map())
  {
    const int sec = val.first;
    const double rot = val.second;

    G4Transform3D sec_place = G4RotateZ3D(rot) * sec_trans;

    ostringstream name;
    name << GetName() << "_sec" << sec;

    G4PVPlacement* calo_phys = new G4PVPlacement(sec_place, sec_logic,
                                                 G4String(name.str()), cylinder_logic, false, sec,
                                                 OverlapCheck());
    calo_vol[calo_phys] = sec;
  }
  _geom->set_nscint(_geom->get_nscint() * _geom->get_sector_map().size());

  // install electronics
  const G4double electronics_thickness = construction_params->get_double_param(
                                             "electronics_thickness") *
                                         cm;
  const string electronics_material = construction_params->get_string_param(
      "electronics_material");

  G4VSolid* electronics_solid = new G4Box(G4String(GetName()),
                                          electronics_thickness / 2.0,  //
                                          sin(
                                              twopi / _geom->get_azimuthal_n_sec() * _geom->get_sector_map().size() / 2) *
                                              (_geom->get_radius() - _geom->get_max_lightguide_height()) * cm,  //
                                          _geom->get_length() * cm / 2.0);

  electronics_solid = new G4DisplacedSolid(G4String(GetName() + "_displaced"),
                                           electronics_solid,
                                           0,  //
                                           G4ThreeVector(
                                               cos(
                                                   twopi / _geom->get_azimuthal_n_sec() * _geom->get_sector_map().size() / 2) *
                                                       (_geom->get_radius() - _geom->get_max_lightguide_height()) * cm -
                                                   electronics_thickness,
                                               0, box_z_shift));

  G4Material* electronics_mat = G4Material::GetMaterial(electronics_material);
  assert(electronics_mat);

  G4LogicalVolume* electronics_logic = new G4LogicalVolume(electronics_solid,
                                                           electronics_mat, G4String(GetName()) + "_Electronics", 0, 0, 0);
  VisAtt = new G4VisAttributes();
  PHG4Utils::SetColour(VisAtt, electronics_material);
  VisAtt->SetVisibility(true);
  VisAtt->SetForceSolid(not _geom->is_virualize_fiber());
  electronics_logic->SetVisAttributes(VisAtt);

  new G4PVPlacement(G4Transform3D::Identity, electronics_logic,
                    G4String(GetName()) + "_Electronics", cylinder_logic, false, 0,
                    OverlapCheck());

  // install the enclosure box
  const G4double enclosure_thickness = construction_params->get_double_param(
                                           "enclosure_thickness") *
                                       cm;
  const string enclosure_material = construction_params->get_string_param(
      "enclosure_material");

  G4VSolid* outer_enclosur_solid = new G4Box(
      G4String(GetName()) + "_outer_enclosur_solid", enclosure_x * 0.5,  //
      enclosure_y * 0.5,                                                 //
      enclosure_z * 0.5);
  G4VSolid* inner_enclosur_solid = new G4Box(
      G4String(GetName()) + "_inner_enclosur_solid",
      enclosure_x * 0.5 - enclosure_thickness,  //
      enclosure_y * 0.5 - enclosure_thickness,  //
      enclosure_z * 0.5 - enclosure_thickness);
  G4VSolid* enclosure_solid = new G4SubtractionSolid(
      G4String(GetName()) + "_enclosure_solid",  //
      outer_enclosur_solid, inner_enclosur_solid);
  enclosure_solid = new G4DisplacedSolid(
      G4String(GetName() + "_enclosure_solid_displaced"), enclosure_solid, 0,  //
      G4ThreeVector(box_x_shift, 0, box_z_shift));

  G4Material* enclosure_mat = G4Material::GetMaterial(enclosure_material);
  assert(enclosure_mat);

  G4LogicalVolume* enclosure_logic = new G4LogicalVolume(enclosure_solid,
                                                         enclosure_mat, G4String(GetName()) + "_enclosure", 0, 0, 0);
  VisAtt = new G4VisAttributes();
  PHG4Utils::SetColour(VisAtt, enclosure_material);
  VisAtt->SetVisibility(true);
  VisAtt->SetForceSolid(not _geom->is_virualize_fiber());
  enclosure_logic->SetVisAttributes(VisAtt);

  new G4PVPlacement(G4Transform3D::Identity, enclosure_logic,
                    G4String(GetName()) + "_enclosure", cylinder_logic, false, 0,
                    OverlapCheck());

  if (active)
  {
    ostringstream geonode;
    if (superdetector != "NONE")
    {
      geonode << "CYLINDERGEOM_" << superdetector;
    }
    else
    {
      geonode << "CYLINDERGEOM_" << detector_type;
    }
    PHG4CylinderGeomContainer* geo = findNode::getClass<
        PHG4CylinderGeomContainer>(topNode(), geonode.str());
    if (!geo)
    {
      geo = new PHG4CylinderGeomContainer();
      PHNodeIterator iter(topNode());
      PHCompositeNode* runNode =
          dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode",
                                                        "RUN"));
      PHIODataNode<PHObject>* newNode = new PHIODataNode<PHObject>(geo,
                                                                   geonode.str(), "PHObject");
      runNode->addNode(newNode);
    }
    // here in the detector class we have internal units, convert to cm
    // before putting into the geom object
    PHG4CylinderGeom* mygeom = new SpacalGeom_t(*_geom);
    geo->AddLayerGeom(0, mygeom);
    if (_geom->get_construction_verbose() >= 1)
    {
      cout << "PHG4SpacalPrototype4Detector::Construct::" << GetName()
           << " - Print Layer Geometry:" << endl;
      geo->identify();
    }
  }

  if (absorberactive)
  {
    ostringstream geonode;
    if (superdetector != "NONE")
    {
      geonode << "CYLINDERGEOM_ABSORBER_" << superdetector;
    }
    else
    {
      geonode << "CYLINDERGEOM_ABSORBER_" << detector_type << "_" << 0;
    }
    PHG4CylinderGeomContainer* geo = findNode::getClass<PHG4CylinderGeomContainer>(topNode(), geonode.str());
    if (!geo)
    {
      geo = new PHG4CylinderGeomContainer();
      PHNodeIterator iter(topNode());
      PHCompositeNode* runNode =
          dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode",
                                                        "RUN"));
      PHIODataNode<PHObject>* newNode = new PHIODataNode<PHObject>(geo, geonode.str(), "PHObject");
      runNode->addNode(newNode);
    }
    // here in the detector class we have internal units, convert to cm
    // before putting into the geom object
    PHG4CylinderGeom* mygeom = new SpacalGeom_t(*_geom);
    geo->AddLayerGeom(0, mygeom);
    if (_geom->get_construction_verbose() >= 1)
    {
      cout << "PHG4SpacalPrototype4Detector::Construct::" << GetName()
           << " - Print Absorber Layer Geometry:" << endl;
      geo->identify();
    }
  }

  if (_geom->get_construction_verbose() >= 1)
  {
    cout << "PHG4SpacalPrototype4Detector::Construct::" << GetName()
         << " - Completed. Print G4 Geometry:" << endl;
    Print();
    cout << "box_x_shift = " << box_x_shift << endl;
    cout << "box_z_shift = " << box_z_shift << endl;
  }
}

std::pair<G4LogicalVolume*, G4Transform3D>
PHG4SpacalPrototype4Detector::Construct_AzimuthalSeg()
{
  assert(_geom);
  assert(_geom->get_azimuthal_n_sec() > 4);

  // basic tilt geometry
  const G4double half_chord_backend =
      _geom->get_max_radius() * cm * tan(pi / _geom->get_azimuthal_n_sec())  //
      + fabs(_geom->get_thickness() * cm * 0.5 * tan(_geom->get_azimuthal_tilt()));
  const G4double reduced_outer_radius = sqrt(pow(_geom->get_max_radius() * cm, 2) - half_chord_backend * half_chord_backend);
  const G4double enclosure_depth = reduced_outer_radius - _geom->get_radius() * cm;
  const G4double enclosure_center = 0.5 * (reduced_outer_radius + _geom->get_radius() * cm);
  const G4double enclosure_half_height_half_width = enclosure_center * tan(pi / _geom->get_azimuthal_n_sec());

  const G4double width_adj1 = tan(_geom->get_azimuthal_tilt() - pi / _geom->get_azimuthal_n_sec()) * enclosure_depth * 0.5;
  const G4double width_adj2 = tan(_geom->get_azimuthal_tilt() + pi / _geom->get_azimuthal_n_sec()) * enclosure_depth * 0.5;

  const G4double center_adj = (width_adj1 + width_adj2) * 0.5;
  const G4double center_tilt_angle = atan2(center_adj, enclosure_depth * 0.5);
  const G4double inner_half_width = enclosure_half_height_half_width + 0.5 * (width_adj1 - width_adj2);
  const G4double outter_half_width = enclosure_half_height_half_width + 0.5 * (-width_adj1 + width_adj2);

  // enclosure walls
  const G4double edge1_tilt_angle = atan2(width_adj1, enclosure_depth * 0.5);
  const G4double edge2_tilt_angle = atan2(width_adj2, enclosure_depth * 0.5);
//  const G4double edge1_half_depth = sqrt(width_adj1 * width_adj1 + enclosure_depth * enclosure_depth * .25);
//  const G4double edge2_half_depth = sqrt(width_adj2 * width_adj2 + enclosure_depth * enclosure_depth * .25);

  // projective center
  const G4double half_projection_ratio = 0.5 * (-width_adj1 + width_adj2) / enclosure_half_height_half_width;
  const G4double projection_center_y = enclosure_center - ((enclosure_depth * 0.5) / half_projection_ratio);
  const G4double projection_center_x = center_adj / half_projection_ratio;

  // blocks azimuthal segmentation
  const int phi_bin_in_sec = _geom->get_max_phi_bin_in_sec();
  assert(phi_bin_in_sec >= 1);
  const G4double block_azimuth_angle = (edge2_tilt_angle - edge1_tilt_angle) / phi_bin_in_sec;
  assert(block_azimuth_angle > 0);
  if (!(fabs(block_azimuth_angle - M_PI * 2 / _geom->get_azimuthal_n_sec() / phi_bin_in_sec) < M_PI * numeric_limits<G4double>::epsilon()))
  {
    cout << "angle/nsec out of range: " << M_PI * numeric_limits<G4double>::epsilon() << endl;
    exit(1);
  }
  const G4double block_edge1_half_width = enclosure_half_height_half_width - (_geom->get_sidewall_thickness() * cm + _geom->get_sidewall_outer_torr() * cm + 2.0 * _geom->get_assembly_spacing() * cm) / cos(edge1_tilt_angle);
  const G4double block_edge2_half_width = enclosure_half_height_half_width - (_geom->get_sidewall_thickness() * cm + _geom->get_sidewall_outer_torr() * cm + 2.0 * _geom->get_assembly_spacing() * cm) / cos(edge2_tilt_angle);
  G4double block_width_ratio = 0;
  for (int s = 0; s < phi_bin_in_sec; ++s)
  {
    block_width_ratio += 1 / cos(block_azimuth_angle * (0.5 + s) + edge1_tilt_angle);
  }
  const G4double block_half_height_width = (block_edge1_half_width + block_edge2_half_width) / block_width_ratio;
  assert(block_half_height_width > 0);

  // write out the azimuthal block geometry
  // block azimuth geometry records
  struct block_azimuth_geom
  {
    G4double angle;
    G4double projection_center_y;
    G4double projection_center_x;
    G4double projection_length;
  };
  vector<block_azimuth_geom> block_azimuth_geoms(phi_bin_in_sec,
                                                 block_azimuth_geom{
                                                     numeric_limits<double>::signaling_NaN(),
                                                     numeric_limits<double>::signaling_NaN(),
                                                     numeric_limits<double>::signaling_NaN(),
                                                     numeric_limits<double>::signaling_NaN()});  // [phi-bin in sector] -> azimuth geometry
  G4double block_x_edge1 = block_edge1_half_width;
  for (int s = 0; s < phi_bin_in_sec; ++s)
  {
    block_azimuth_geom& geom = block_azimuth_geoms[s];

    geom.angle = block_azimuth_angle * (0.5 + s) + edge1_tilt_angle;
    const G4double block_x_size = block_half_height_width / cos(geom.angle);
    assert(block_x_size > 0);
    const G4double x_center = block_x_edge1 - 0.5 * block_x_size;

    // projection center per block
    geom.projection_length = block_half_height_width / 2. / tan(block_azimuth_angle / 2.);
    assert(geom.projection_length > 0);
    geom.projection_center_y = enclosure_center - geom.projection_length * cos(geom.angle);
    geom.projection_center_x = x_center + geom.projection_length * sin(geom.angle);

    // next step
    block_x_edge1 -= block_x_size;
  }

  //write out the azimuthal block divider's geometry
  struct block_divider_azimuth_geom
  {
    G4double angle;  //! rotation angle
    G4double projection_center_y;
    G4double projection_center_x;
    G4double thickness;            // thickness in the approximate azimuth direction
    G4double radial_displacement;  //! displacement along the width direction, which is the radial direction if tilt = 0
    G4double width;                //! wdith along the approximate radial direction
  };
  assert(phi_bin_in_sec >= 1);
  vector<block_divider_azimuth_geom> divider_azimuth_geoms(phi_bin_in_sec - 1,
                                                           block_divider_azimuth_geom{
                                                               numeric_limits<double>::signaling_NaN(),
                                                               numeric_limits<double>::signaling_NaN(),
                                                               numeric_limits<double>::signaling_NaN(),
                                                               numeric_limits<double>::signaling_NaN(),
                                                               numeric_limits<double>::signaling_NaN(),
                                                               numeric_limits<double>::signaling_NaN()});

  if (_geom->get_sidewall_thickness() > 0)
  {
    for (int s = 0; s < phi_bin_in_sec - 1; ++s)
    {
      block_divider_azimuth_geom& geom = divider_azimuth_geoms[s];

      geom.angle = 0.5 * (block_azimuth_geoms[s].angle + block_azimuth_geoms[s + 1].angle);
      geom.projection_center_y = 0.5 * (block_azimuth_geoms[s].projection_center_y + block_azimuth_geoms[s + 1].projection_center_y);
      geom.projection_center_x = 0.5 * (block_azimuth_geoms[s].projection_center_x + block_azimuth_geoms[s + 1].projection_center_x);
      geom.radial_displacement = 0.5 * (block_azimuth_geoms[s].projection_length + block_azimuth_geoms[s + 1].projection_length);

      geom.thickness = 2.0 * _geom->get_assembly_spacing() * cm * cos(block_azimuth_angle / 2.) - 2 * um;
      geom.width = _geom->get_divider_width() * cm;
    }
  }

  if (fabs(block_x_edge1 - (-block_edge2_half_width)) > _geom->get_assembly_spacing() * cm)
  {
    cout << "PHG4SpacalPrototype4Detector::Construct_AzimuthalSeg - ERROR - " << endl
         << "\t block_x_edge1 = " << block_x_edge1 << endl
         << "\t block_edge2_half_width = " << block_edge2_half_width << endl
         << "\t fabs(block_x_edge1 - (-block_edge2_half_width)) = " << fabs(block_x_edge1 - (-block_edge2_half_width)) << endl
         << "\t _geom->get_assembly_spacing() * cm = " << _geom->get_assembly_spacing() * cm << endl;
  }
  if (!(fabs(block_x_edge1 - (-block_edge2_half_width)) < _geom->get_assembly_spacing() * cm))  // closure check
  {
    cout << "PHG4SpacalPrototype4Detector::Construct_AzimuthalSeg - ERROR - closure check failed: " << fabs(block_x_edge1 - (-block_edge2_half_width)) << endl;
    exit(1);
  }

  if (Verbosity())
  {
    cout << "PHG4SpacalPrototype4Detector::Construct_AzimuthalSeg - " << endl
         << "\t edge1_tilt_angle = " << edge1_tilt_angle << endl
         << "\t edge2_tilt_angle = " << edge2_tilt_angle << endl
         << "\t projection_center_y = " << projection_center_y << endl
         << "\t projection_center_x = " << projection_center_x << endl
         << "\t block_azimuth_angle = " << block_azimuth_angle << endl
         << "\t block_edge1_half_width = " << block_edge1_half_width << endl
         << "\t block_edge2_half_width = " << block_edge2_half_width << endl
         << "\t block_width_ratio = " << block_width_ratio << endl
         << "\t block_half_height_width = " << block_half_height_width << endl;

    for (int s = 0; s < phi_bin_in_sec; ++s)
    {
      cout << "\t block[" << s << "].angle = " << block_azimuth_geoms[s].angle << endl;
      cout << "\t block[" << s << "].projection_center_y = " << block_azimuth_geoms[s].projection_center_y << endl;
      cout << "\t block[" << s << "].projection_center_x = " << block_azimuth_geoms[s].projection_center_x << endl;
    }
    for (int s = 0; s < phi_bin_in_sec - 1; ++s)
    {
      cout << "\t divider[" << s << "].angle = " << divider_azimuth_geoms[s].angle << endl;
      cout << "\t divider[" << s << "].projection_center_x = " << divider_azimuth_geoms[s].projection_center_x << endl;
      cout << "\t divider[" << s << "].projection_center_y = " << divider_azimuth_geoms[s].projection_center_y << endl;
      cout << "\t divider[" << s << "].radial_displacement = " << divider_azimuth_geoms[s].radial_displacement << endl;
      cout << "\t divider[" << s << "].thickness = " << divider_azimuth_geoms[s].thickness << endl;
      cout << "\t divider[" << s << "].width = " << divider_azimuth_geoms[s].width << endl;
    }
  }

  assert(enclosure_depth > 10 * cm);

  G4VSolid* sec_solid = new G4Trap(
      /*const G4String& pName*/ G4String(GetName() + string("_sec_trap")),
      enclosure_depth * 0.5,                                                              // G4double pDz,
      center_tilt_angle, halfpi,                                                          // G4double pTheta, G4double pPhi,
      inner_half_width, _geom->get_length() * cm / 2.0, _geom->get_length() * cm / 2.0,   // G4double pDy1, G4double pDx1, G4double pDx2,
      0,                                                                                  // G4double pAlp1,
      outter_half_width, _geom->get_length() * cm / 2.0, _geom->get_length() * cm / 2.0,  // G4double pDy2, G4double pDx3, G4double pDx4,
      0                                                                                   // G4double pAlp2 //
  );
  const G4double box_z_shift = 0.5 * (_geom->get_zmin() + _geom->get_zmax()) * cm;
  G4Transform3D sec_solid_transform =
      G4TranslateZ3D(box_z_shift) * G4TranslateY3D(enclosure_center) * G4RotateY3D(halfpi) * G4RotateX3D(-halfpi);
  G4VSolid* sec_solid_place = new G4DisplacedSolid(
      G4String(GetName() + string("_sec")), sec_solid, sec_solid_transform);

  G4Material* cylinder_mat = G4Material::GetMaterial("G4_AIR");
  assert(cylinder_mat);

  G4LogicalVolume* sec_logic = new G4LogicalVolume(sec_solid_place, cylinder_mat,
                                                   G4String(G4String(GetName() + string("_sec"))), 0, 0, nullptr);


  G4VisAttributes* VisAtt = new G4VisAttributes();
  VisAtt->SetColor(.5, .9, .5, .1);
  VisAtt->SetVisibility(true);
  VisAtt->SetForceSolid(false);
  VisAtt->SetForceWireframe(true);
  sec_logic->SetVisAttributes(VisAtt);

//  // construct walls
//
//  G4Material* wall_mat = G4Material::GetMaterial(_geom->get_sidewall_mat());
//  assert(wall_mat);
//
//  if (_geom->get_sidewall_thickness() > 0)
//  {
//    // end walls
//    if (_geom->get_construction_verbose() >= 1)
//    {
//      cout << "PHG4SpacalPrototype4Detector::Construct_AzimuthalSeg::" << GetName()
//           << " - construct end walls." << endl;
//    }
//    //        G4Tubs* wall_solid = new G4Tubs(G4String(GetName() + string("_EndWall")),
//    //            _geom->get_radius() * cm + _geom->get_sidewall_outer_torr() * cm,
//    //            _geom->get_max_radius() * cm - _geom->get_sidewall_outer_torr() * cm,
//    //            _geom->get_sidewall_thickness() * cm / 2.0,
//    //            halfpi - pi / _geom->get_azimuthal_n_sec(),
//    //            twopi / _geom->get_azimuthal_n_sec());
//    const G4double side_wall_half_thickness = _geom->get_sidewall_thickness() * cm / 2.0;
//    G4VSolid* wall_solid = new G4Trap(G4String(GetName() + string("_EndWall_trap")),
//                                      enclosure_depth * 0.5,                                                  // G4double pDz,
//                                      center_tilt_angle, halfpi,                                              // G4double pTheta, G4double pPhi,
//                                      inner_half_width, side_wall_half_thickness, side_wall_half_thickness,   // G4double pDy1, G4double pDx1, G4double pDx2,
//                                      0,                                                                      // G4double pAlp1,
//                                      outter_half_width, side_wall_half_thickness, side_wall_half_thickness,  // G4double pDy2, G4double pDx3, G4double pDx4,
//                                      0                                                                       // G4double pAlp2 //
//    );
//    G4VSolid* wall_solid_place = new G4DisplacedSolid(
//        G4String(GetName() + string("_EndWall")), wall_solid, sec_solid_transform);
//
//    G4LogicalVolume* wall_logic = new G4LogicalVolume(wall_solid_place, wall_mat,
//                                                      G4String(G4String(GetName() + string("_EndWall"))), 0, 0,
//                                                      nullptr);
//    GetDisplayAction()->AddVolume(wall_logic, "Wall");
//
//    typedef map<int, double> z_locations_t;
//    z_locations_t z_locations;
//    z_locations[000] = _geom->get_sidewall_thickness() * cm / 2.0 + _geom->get_assembly_spacing() * cm;
//    z_locations[001] = _geom->get_length() * cm / 2.0 - (_geom->get_sidewall_thickness() * cm / 2.0 + _geom->get_assembly_spacing() * cm);
//    z_locations[100] = -(_geom->get_sidewall_thickness() * cm / 2.0 + _geom->get_assembly_spacing() * cm);
//    z_locations[101] = -(_geom->get_length() * cm / 2.0 - (_geom->get_sidewall_thickness() * cm / 2.0 + _geom->get_assembly_spacing() * cm));
//
//    BOOST_FOREACH (z_locations_t::value_type& val, z_locations)
//    {
//      if (_geom->get_construction_verbose() >= 2)
//        cout << "PHG4SpacalPrototype4Detector::Construct_AzimuthalSeg::"
//             << GetName() << " - constructed End Wall ID " << val.first
//             << " @ Z = " << val.second << endl;
//
//      G4Transform3D wall_trans = G4TranslateZ3D(val.second);
//
//      G4PVPlacement* wall_phys = new G4PVPlacement(wall_trans, wall_logic,
//                                                   G4String(GetName()) + G4String("_EndWall_") + to_string(val.first), sec_logic,
//                                                   false, val.first, OverlapCheck());
//
//      calo_vol[wall_phys] = val.first;
//      assert(gdml_config);
//      gdml_config->exclude_physical_vol(wall_phys);
//    }
//  }
//  //
//  if (_geom->get_sidewall_thickness() > 0)
//  {
//    // side walls
//    if (_geom->get_construction_verbose() >= 1)
//    {
//      cout << "PHG4SpacalPrototype4Detector::Construct_AzimuthalSeg::" << GetName()
//           << " - construct side walls." << endl;
//    }
//
//    typedef map<int, pair<int, int> > sign_t;
//    sign_t signs;
//    signs[100] = make_pair(+1, +1);
//    signs[101] = make_pair(+1, -1);
//    signs[200] = make_pair(-1, +1);
//    signs[201] = make_pair(-1, -1);
//
//    BOOST_FOREACH (sign_t::value_type& val, signs)
//    {
//      const int sign_z = val.second.first;
//      const int sign_azimuth = val.second.second;
//
//      if (_geom->get_construction_verbose() >= 2)
//        cout << "PHG4SpacalPrototype4Detector::Construct_AzimuthalSeg::"
//             << GetName() << " - constructed Side Wall ID " << val.first
//             << " with"
//             << " Shift X = "
//             << sign_azimuth * (_geom->get_sidewall_thickness() * cm / 2.0 + _geom->get_sidewall_outer_torr() * cm)
//             << " Rotation Z = "
//             << sign_azimuth * pi / _geom->get_azimuthal_n_sec()
//             << " Shift Z = " << sign_z * (_geom->get_length() * cm / 4)
//             << endl;
//
//      const G4double azimuth_roate =
//          sign_azimuth > 0 ? edge1_tilt_angle : edge2_tilt_angle;
//      const G4double edge_half_depth = -_geom->get_sidewall_thickness() * cm - _geom->get_sidewall_outer_torr() * cm +
//                                       (sign_azimuth > 0 ? edge1_half_depth : edge2_half_depth);
//
//      G4Box* wall_solid = new G4Box(G4String(GetName() + G4String("_SideWall_") + to_string(val.first)),
//                                    _geom->get_sidewall_thickness() * cm / 2.0,
//                                    edge_half_depth,
//                                    (_geom->get_length() / 2. - 2 * (_geom->get_sidewall_thickness() + 2. * _geom->get_assembly_spacing())) * cm * .5);
//
//      G4LogicalVolume* wall_logic = new G4LogicalVolume(wall_solid, wall_mat,
//                                                        G4String(G4String(GetName() + G4String("_SideWall_") + to_string(val.first))), 0, 0,
//                                                        nullptr);
//      GetDisplayAction()->AddVolume(wall_logic, "Wall");
//
//      const G4Transform3D wall_trans =
//          G4TranslateZ3D(sign_z * (_geom->get_length() * cm / 4)) *
//          G4TranslateY3D(enclosure_center) *
//          G4TranslateX3D(sign_azimuth * enclosure_half_height_half_width) *
//          G4RotateZ3D(azimuth_roate) *
//          G4TranslateX3D(-sign_azimuth * (_geom->get_sidewall_thickness() * cm / 2.0 + _geom->get_sidewall_outer_torr() * cm));
//
//      G4PVPlacement* wall_phys = new G4PVPlacement(wall_trans, wall_logic,
//                                                   G4String(GetName()) + G4String("_SideWall_") + to_string(val.first), sec_logic,
//                                                   false, val.first, OverlapCheck());
//
//      calo_vol[wall_phys] = val.first;
//      assert(gdml_config);
//      gdml_config->exclude_physical_vol(wall_phys);
//    }
//  }
//
//  // construct dividers
//  if (_geom->get_divider_width() > 0)
//  {
//    if (_geom->get_construction_verbose() >= 1)
//    {
//      cout << "PHG4SpacalPrototype4Detector::Construct_AzimuthalSeg::" << GetName()
//           << " - construct dividers" << endl;
//    }
//
//    G4Material* divider_mat = G4Material::GetMaterial(_geom->get_divider_mat());
//    assert(divider_mat);
//
//    int ID = 300;
//    for (const auto& geom : divider_azimuth_geoms)
//    {
//      G4Box* divider_solid = new G4Box(G4String(GetName() + G4String("_Divider_") + to_string(ID)),
//                                       geom.thickness / 2.0,
//                                       geom.width / 2.,
//                                       (_geom->get_length() / 2. - 2 * (_geom->get_sidewall_thickness() + 2. * _geom->get_assembly_spacing())) * cm * .5);
//
//      G4LogicalVolume* wall_logic = new G4LogicalVolume(divider_solid, divider_mat,
//                                                        G4String(G4String(GetName() + G4String("_Divider_") + to_string(ID))), 0, 0,
//                                                        nullptr);
//      GetDisplayAction()->AddVolume(wall_logic, "Divider");
//
//      for (int sign_z = -1; sign_z <= 1; sign_z += 2)
//      {
//        G4Transform3D wall_trans =
//            G4TranslateX3D(geom.projection_center_x) *
//            G4TranslateY3D(geom.projection_center_y) *
//            G4RotateZ3D(geom.angle) *
//            G4TranslateY3D(geom.radial_displacement) *
//            G4TranslateZ3D(sign_z * (_geom->get_length() * cm / 4));
//
//        G4PVPlacement* wall_phys = new G4PVPlacement(wall_trans, wall_logic,
//                                                     G4String(GetName()) + G4String("_Divider_") + to_string(ID), sec_logic,
//                                                     false, ID, OverlapCheck());
//
//        calo_vol[wall_phys] = ID;
//        //        assert(gdml_config);
//        //        gdml_config->exclude_physical_vol(wall_phys);
//
//        if (Verbosity())
//        {
//          cout << "PHG4SpacalPrototype4Detector::Construct_AzimuthalSeg - placing divider " << wall_phys->GetName() << " copy ID " << ID << endl;
//        }
//
//        ++ID;
//      }
//    }  //    for (const auto & geom : divider_azimuth_geoms)
//  }

  //  // construct towers
  //
  BOOST_FOREACH (const SpacalGeom_t::tower_map_t::value_type& val, _geom->get_sector_tower_map())
  {
    SpacalGeom_t::geom_tower g_tower = val.second;

    const int tower_id = g_tower.id;
    const int tower_phi_id_in_sec = tower_id % 10;
    assert(tower_phi_id_in_sec >= 0);
    assert(tower_phi_id_in_sec < phi_bin_in_sec);

    const auto& block_azimuth_geom = block_azimuth_geoms.at(tower_phi_id_in_sec);

    G4LogicalVolume* LV_tower = Construct_Tower(g_tower);

    G4Transform3D block_trans =
        G4TranslateX3D(block_azimuth_geom.projection_center_x) *
        G4TranslateY3D(block_azimuth_geom.projection_center_y) *
        G4RotateZ3D(block_azimuth_geom.angle) *
        G4TranslateX3D(g_tower.centralX * cm) *
        G4TranslateY3D(g_tower.centralY * cm) *
        G4TranslateZ3D(g_tower.centralZ * cm) *
        G4RotateX3D(g_tower.pRotationAngleX * rad);

    const bool overlapcheck_block = OverlapCheck() and (_geom->get_construction_verbose() >= 2);

    G4PVPlacement* block_phys = new G4PVPlacement(block_trans, LV_tower,
                                                  G4String(GetName()) + G4String("_Tower_") + to_string(g_tower.id), sec_logic, false,
                                                  g_tower.id, overlapcheck_block);
    block_vol[block_phys] = g_tower.id;
    //    assert(gdml_config);
    //    gdml_config->exclude_physical_vol(block_phys);

    if (g_tower.LightguideHeight > 0)
    {
      // also build a light guide

      for (int ix = 0; ix < g_tower.NSubtowerX; ix++)
      //  int ix = 0;
      {
        for (int iy = 0; iy < g_tower.NSubtowerY; iy++)
        //        int iy = 0;
        {
          G4LogicalVolume* LV_lg = Construct_LightGuide(g_tower, ix,
                                                        iy);

          G4PVPlacement* lg_phys = new G4PVPlacement(block_trans, LV_lg, LV_lg->GetName(),
                                                     sec_logic, false, g_tower.id, overlapcheck_block);

          block_vol[lg_phys] = g_tower.id * 100 + ix * 10 + iy;

          //          assert(gdml_config);
          //          gdml_config->exclude_physical_vol(lg_phys);
        }
      }
    }
  }

  cout << "PHG4SpacalPrototype4Detector::Construct_AzimuthalSeg::" << GetName()
       << " - constructed " << _geom->get_sector_tower_map().size()
       << " unique towers" << endl;

  return make_pair(sec_logic, G4Transform3D::Identity);
}

G4LogicalVolume*
PHG4SpacalPrototype4Detector::Construct_Fiber(const G4double length,
                                              const string& id)
{
  G4Tubs* fiber_solid = new G4Tubs(G4String(GetName() + string("_fiber") + id),
                                   0, _geom->get_fiber_outer_r() * cm, length / 2.0, 0, twopi);

  G4Material* clading_mat = G4Material::GetMaterial(
      _geom->get_fiber_clading_mat());
  assert(clading_mat);

  G4LogicalVolume* fiber_logic = new G4LogicalVolume(fiber_solid, clading_mat,
                                                     G4String(G4String(GetName() + string("_fiber") + id)), 0, 0,
                                                     clading_step_limits);

  {
    G4VisAttributes* VisAtt = new G4VisAttributes();
    PHG4Utils::SetColour(VisAtt, "G4_POLYSTYRENE");
    VisAtt->SetVisibility(_geom->is_virualize_fiber());
    VisAtt->SetForceSolid(_geom->is_virualize_fiber());
    fiber_logic->SetVisAttributes(VisAtt);
  }

  G4Tubs* core_solid = new G4Tubs(
      G4String(GetName() + string("_fiber_core") + id), 0,
      _geom->get_fiber_core_diameter() * cm / 2, length / 2.0, 0, twopi);

  G4Material* core_mat = G4Material::GetMaterial(_geom->get_fiber_core_mat());
  assert(core_mat);

  G4LogicalVolume* core_logic = new G4LogicalVolume(core_solid, core_mat,
                                                    G4String(G4String(GetName() + string("_fiber_core") + id)), 0, 0,
                                                    fiber_core_step_limits);

  {
    G4VisAttributes* VisAtt = new G4VisAttributes();
    PHG4Utils::SetColour(VisAtt, "G4_POLYSTYRENE");
    VisAtt->SetVisibility(false);
    VisAtt->SetForceSolid(false);
    core_logic->SetVisAttributes(VisAtt);
  }

  const bool overlapcheck_fiber = OverlapCheck() and (_geom->get_construction_verbose() >= 3);
  G4PVPlacement* core_physi = new G4PVPlacement(0, G4ThreeVector(), core_logic,
                                                G4String(G4String(GetName() + string("_fiber_core") + id)), fiber_logic,
                                                false, 0, overlapcheck_fiber);
  fiber_core_vol[core_physi] = 0;

  return fiber_logic;
}

void PHG4SpacalPrototype4Detector::Print(const std::string& what) const
{
  assert(_geom);

  cout << "PHG4SpacalPrototype4Detector::Print::" << GetName()
       << " - Print Geometry:" << endl;
  _geom->Print();

  return;
}

//! Fully projective spacal with 2D tapered modules. To speed up construction, same-length fiber is used cross one tower
int PHG4SpacalPrototype4Detector::Construct_Fibers_SameLengthFiberPerTower(
    const PHG4SpacalPrototype4Detector::SpacalGeom_t::geom_tower& g_tower,
    G4LogicalVolume* LV_tower)
{
  assert(_geom);

  // construct fibers

  // first check out the fibers geometry

  typedef map<int, pair<G4Vector3D, G4Vector3D> > fiber_par_map;
  fiber_par_map fiber_par;
  G4double min_fiber_length = g_tower.pDz * cm * 4;

  G4Vector3D v_zshift = G4Vector3D(tan(g_tower.pTheta) * cos(g_tower.pPhi),
                                   tan(g_tower.pTheta) * sin(g_tower.pPhi), 1) *
                        g_tower.pDz;
  //  int fiber_ID = 0;
  for (int ix = 0; ix < g_tower.NFiberX; ix++)
  //  int ix = 0;
  {
    const double weighted_ix = static_cast<double>(ix) / (g_tower.NFiberX - 1.);

    const double weighted_pDx1 = (g_tower.pDx1 - g_tower.ModuleSkinThickness - _geom->get_fiber_outer_r()) * (weighted_ix * 2 - 1);
    const double weighted_pDx2 = (g_tower.pDx2 - g_tower.ModuleSkinThickness - _geom->get_fiber_outer_r()) * (weighted_ix * 2 - 1);

    const double weighted_pDx3 = (g_tower.pDx3 - g_tower.ModuleSkinThickness - _geom->get_fiber_outer_r()) * (weighted_ix * 2 - 1);
    const double weighted_pDx4 = (g_tower.pDx4 - g_tower.ModuleSkinThickness - _geom->get_fiber_outer_r()) * (weighted_ix * 2 - 1);

    for (int iy = 0; iy < g_tower.NFiberY; iy++)
    //        int iy = 0;
    {
      if ((ix + iy) % 2 == 1)
        continue;  // make a triangle pattern

      const double weighted_iy = static_cast<double>(iy) / (g_tower.NFiberY - 1.);

      const double weighted_pDy1 = (g_tower.pDy1 - g_tower.ModuleSkinThickness - _geom->get_fiber_outer_r()) * (weighted_iy * 2 - 1);
      const double weighted_pDy2 = (g_tower.pDy2 - g_tower.ModuleSkinThickness - _geom->get_fiber_outer_r()) * (weighted_iy * 2 - 1);

      const double weighted_pDx12 = weighted_pDx1 * (1 - weighted_iy) + weighted_pDx2 * (weighted_iy) + weighted_pDy1 * tan(g_tower.pAlp1);
      const double weighted_pDx34 = weighted_pDx3 * (1 - weighted_iy) + weighted_pDx4 * (weighted_iy) + weighted_pDy1 * tan(g_tower.pAlp2);

      G4Vector3D v1 = G4Vector3D(weighted_pDx12, weighted_pDy1, 0) - v_zshift;
      G4Vector3D v2 = G4Vector3D(weighted_pDx34, weighted_pDy2, 0) + v_zshift;

      G4Vector3D vector_fiber = (v2 - v1);
      vector_fiber *= (vector_fiber.mag() - _geom->get_fiber_outer_r()) / vector_fiber.mag();  // shrink by fiber boundary protection
      G4Vector3D center_fiber = (v2 + v1) / 2;

      // convert to Geant4 units
      vector_fiber *= cm;
      center_fiber *= cm;

      const int fiber_ID = g_tower.compose_fiber_id(ix, iy);
      fiber_par[fiber_ID] = make_pair(vector_fiber, center_fiber);

      const G4double fiber_length = vector_fiber.mag();

      min_fiber_length = min(fiber_length, min_fiber_length);

      //          ++fiber_ID;
    }
  }

  int fiber_count = 0;

  const G4double fiber_length = min_fiber_length;
  vector<G4double> fiber_cut;

  ostringstream ss;
  ss << string("_Tower") << g_tower.id;
  G4LogicalVolume* fiber_logic = Construct_Fiber(fiber_length, ss.str());

  BOOST_FOREACH (const fiber_par_map::value_type& val, fiber_par)
  {
    const int fiber_ID = val.first;
    G4Vector3D vector_fiber = val.second.first;
    G4Vector3D center_fiber = val.second.second;
    const G4double optimal_fiber_length = vector_fiber.mag();

    const G4Vector3D v1 = center_fiber - 0.5 * vector_fiber;

    // keep a statistics
    assert(optimal_fiber_length - fiber_length >= 0);
    fiber_cut.push_back(optimal_fiber_length - fiber_length);

    center_fiber += (fiber_length / optimal_fiber_length - 1) * 0.5 * vector_fiber;
    vector_fiber *= fiber_length / optimal_fiber_length;

    //      const G4Vector3D v1_new = center_fiber - 0.5 *vector_fiber;

    if (_geom->get_construction_verbose() >= 3)
      cout
          << "PHG4SpacalPrototype4Detector::Construct_Fibers_SameLengthFiberPerTower::"
          << GetName() << " - constructed fiber " << fiber_ID << ss.str()
          //
          << ", Length = " << optimal_fiber_length << "-"
          << (optimal_fiber_length - fiber_length) << "mm, "  //
          << "x = " << center_fiber.x() << "mm, "             //
          << "y = " << center_fiber.y() << "mm, "             //
          << "z = " << center_fiber.z() << "mm, "             //
          << "vx = " << vector_fiber.x() << "mm, "            //
          << "vy = " << vector_fiber.y() << "mm, "            //
          << "vz = " << vector_fiber.z() << "mm, "            //
          << endl;

    const G4double rotation_angle = G4Vector3D(0, 0, 1).angle(vector_fiber);
    const G4Vector3D rotation_axis =
        rotation_angle == 0 ? G4Vector3D(1, 0, 0) : G4Vector3D(0, 0, 1).cross(vector_fiber);

    G4Transform3D fiber_place(
        G4Translate3D(center_fiber.x(), center_fiber.y(), center_fiber.z()) * G4Rotate3D(rotation_angle, rotation_axis));

    ostringstream name;
    name << GetName() + string("_Tower") << g_tower.id << "_fiber"
         << ss.str();

    const bool overlapcheck_fiber = OverlapCheck() and (_geom->get_construction_verbose() >= 3);
    G4PVPlacement* fiber_physi = new G4PVPlacement(fiber_place, fiber_logic,
                                                   G4String(name.str()), LV_tower, false, fiber_ID,
                                                   overlapcheck_fiber);
    fiber_vol[fiber_physi] = fiber_ID;

    fiber_count++;
  }

  if (_geom->get_construction_verbose() >= 2)
    cout
        << "PHG4SpacalPrototype4Detector::Construct_Fibers_SameLengthFiberPerTower::"
        << GetName() << " - constructed tower ID " << g_tower.id << " with "
        << fiber_count << " fibers. Average fiber length cut = "
        << accumulate(fiber_cut.begin(), fiber_cut.end(), 0.0) / fiber_cut.size() << " mm" << endl;

  return fiber_count;
}

//! a block along z axis built with G4Trd that is slightly tapered in x dimension
G4LogicalVolume*
PHG4SpacalPrototype4Detector::Construct_Tower(
    const PHG4SpacalPrototype4Detector::SpacalGeom_t::geom_tower& g_tower)
{
  assert(_geom);

  if (_geom->get_construction_verbose() >= 2)
  {
    cout << "PHG4SpacalPrototype4Detector::Construct_Tower::" << GetName()
         << " - constructed tower ID " << g_tower.id
         << " with geometry parameter: ";
    g_tower.identify(cout);
  }

  std::ostringstream sout;
  sout << "_" << g_tower.id;
  const G4String sTowerID(sout.str());

  //Processed PostionSeeds 1 from 1 1

  G4Trap* block_solid = new G4Trap(
      /*const G4String& pName*/ G4String(GetName()) + sTowerID,
      g_tower.pDz * cm,                                         // G4double pDz,
      g_tower.pTheta * rad, g_tower.pPhi * rad,                 // G4double pTheta, G4double pPhi,
      g_tower.pDy1 * cm, g_tower.pDx1 * cm, g_tower.pDx2 * cm,  // G4double pDy1, G4double pDx1, G4double pDx2,
      g_tower.pAlp1 * rad,                                      // G4double pAlp1,
      g_tower.pDy2 * cm, g_tower.pDx3 * cm, g_tower.pDx4 * cm,  // G4double pDy2, G4double pDx3, G4double pDx4,
      g_tower.pAlp2 * rad                                       // G4double pAlp2 //
  );

  G4Material* cylinder_mat = G4Material::GetMaterial(
      _geom->get_absorber_mat());
  assert(cylinder_mat);

  G4LogicalVolume* block_logic = new G4LogicalVolume(block_solid, cylinder_mat,
                                                     G4String(G4String(GetName()) + string("_Tower") + sTowerID), 0, 0,
                                                     step_limits);

  G4VisAttributes* VisAtt = new G4VisAttributes();
  //  PHG4Utils::SetColour(VisAtt, "W_Epoxy");
  VisAtt->SetColor(.3, .3, .3, .3);
  VisAtt->SetVisibility(
      _geom->is_azimuthal_seg_visible() or _geom->is_virualize_fiber());
  VisAtt->SetForceSolid(not _geom->is_virualize_fiber());
  block_logic->SetVisAttributes(VisAtt);

  // construct fibers

  int fiber_count = 0;

  fiber_count = Construct_Fibers_SameLengthFiberPerTower(g_tower, block_logic);

  if (_geom->get_construction_verbose() >= 2)
    cout << "PHG4SpacalPrototype4Detector::Construct_Tower::" << GetName()
         << " - constructed tower ID " << g_tower.id << " with " << fiber_count
         << " fibers using Construct_Fibers_SameLengthFiberPerTower" << endl;

  return block_logic;
}

G4LogicalVolume*
PHG4SpacalPrototype4Detector::Construct_LightGuide(
    const PHG4SpacalPrototype4Detector::SpacalGeom_t::geom_tower& g_tower,
    const int index_x, const int index_y)
{
  assert(_geom);

  std::ostringstream sout;
  sout << "_Lightguide_" << g_tower.id << "_" << index_x << "_" << index_y;
  const G4String sTowerID(sout.str());

  assert(g_tower.LightguideHeight > 0);

  // light guide parameters in PHENIX units
  const double weight_x1 = 1 - (double) index_y / g_tower.NSubtowerY;
  const double weight_x2 = 1 - (double) (index_y + 1) / g_tower.NSubtowerY;
  const double weight_xcenter = 1 - (double) (index_y + 0.5) / g_tower.NSubtowerY;

  assert(weight_x1 >= 0 and weight_x1 <= 1);
  assert(weight_x2 >= 0 and weight_x2 <= 1);
  assert(weight_xcenter >= 0 and weight_xcenter <= 1);

  const double lg_pDx1 = (g_tower.pDx1 * weight_x1  //
                          + g_tower.pDx2 * (1 - weight_x1)) /
                         g_tower.NSubtowerX;
  const double lg_pDx2 = (g_tower.pDx1 * weight_x2  //
                          + g_tower.pDx2 * (1 - weight_x2)) /
                         g_tower.NSubtowerX;
  const double lg_pDy1 = g_tower.pDy1 / g_tower.NSubtowerY;
  const double lg_Alp1 = atan(
      (g_tower.pDx2 - g_tower.pDx1) * (-g_tower.NSubtowerX + 1. + 2 * index_x) / (double) (g_tower.NSubtowerX) / (2. * g_tower.pDy1) + tan(g_tower.pAlp1));

  const double shift_xcenter = (g_tower.pDx1 * weight_xcenter           //
                                + g_tower.pDx2 * (1 - weight_xcenter))  //
                               *                                        //
                               (-g_tower.NSubtowerX + 1. + 2 * index_x) / (double) (g_tower.NSubtowerX);
  const double shift_ycenter = g_tower.pDy1  //
                               *             //
                               (-g_tower.NSubtowerY + 1. + 2 * index_y) / (double) (g_tower.NSubtowerY);

  G4VSolid* block_solid = new G4Trap(
      /*const G4String& pName*/ G4String(GetName()) + sTowerID,
      0.5 * g_tower.LightguideHeight * cm,  // G4double pDz,
      0 * rad, 0 * rad,                     // G4double pTheta, G4double pPhi,
      g_tower.LightguideTaperRatio * lg_pDy1 * cm,
      g_tower.LightguideTaperRatio * lg_pDx1 * cm,
      g_tower.LightguideTaperRatio * lg_pDx2 * cm,  // G4double pDy1, G4double pDx1, G4double pDx2,
      lg_Alp1 * rad,                                // G4double pAlp1,
      lg_pDy1 * cm, lg_pDx1 * cm, lg_pDx2 * cm,     // G4double pDy2, G4double pDx3, G4double pDx4,
      lg_Alp1 * rad                                 // G4double pAlp2 //
  );

  block_solid = new G4DisplacedSolid(G4String(GetName() + "_displaced"),
                                     block_solid, 0,                                           //
                                     G4ThreeVector(                                            //
                                         tan(g_tower.pTheta * rad) * cos(g_tower.pPhi * rad),  //
                                         tan(g_tower.pTheta * rad) * sin(g_tower.pPhi * rad),  //
                                         1) *                                                  // G4ThreeVector
                                             -(g_tower.pDz) *
                                             cm                                                       //
                                         + G4ThreeVector(shift_xcenter * cm, shift_ycenter * cm, 0)   // shit in subtower direction
                                         + G4ThreeVector(0, 0, -0.5 * g_tower.LightguideHeight * cm)  //shift in the light guide height
  );

  G4Material* cylinder_mat = G4Material::GetMaterial(
      g_tower.LightguideMaterial);
  assert(cylinder_mat);

  G4LogicalVolume* block_logic = new G4LogicalVolume(block_solid, cylinder_mat,
                                                     G4String(G4String(GetName()) + string("_Tower") + sTowerID), 0, 0,
                                                     step_limits);

  G4VisAttributes* VisAtt = new G4VisAttributes();
  PHG4Utils::SetColour(VisAtt, g_tower.LightguideMaterial);
  //  VisAtt->SetColor(.3, .3, .3, .3);
  VisAtt->SetVisibility(
      _geom->is_azimuthal_seg_visible() or _geom->is_virualize_fiber());
  VisAtt->SetForceSolid(not _geom->is_virualize_fiber());
  block_logic->SetVisAttributes(VisAtt);

  return block_logic;
}
