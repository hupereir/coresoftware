#ifndef PHHEPMC_PHHEPMCGENEVENT_H
#define PHHEPMC_PHHEPMCGENEVENT_H

#include <phool/PHObject.h>

#include <phool/phool.h>

#include <HepMC/SimpleVector.h>

#include <CLHEP/Vector/LorentzRotation.h>

#include <iostream>  // for cout, ostream

#include <map>

namespace HepMC
{
  class GenEvent;
}

class PHHepMCGenEvent : public PHObject
{
 public:

  //! default constructor
  explicit PHHepMCGenEvent() = default;

  //! copy constructor
  PHHepMCGenEvent(const PHHepMCGenEvent& event);

  //! assignment operator
  PHHepMCGenEvent& operator=(const PHHepMCGenEvent& event);

  //! destructor
  ~PHHepMCGenEvent() override;


  //! \brief reaction plane angles thrown by hijing flowAfterburner 
  // ! only will be filled if the flowAfterburner enabled and input HepMC is from hijing
  // ! it is a map of n to psi_n[n]
  typedef std::map<unsigned int, float> FlowAfterburner_PsiMap;


  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override
  {
    return (getEvent() != nullptr) ? 1 : 0;
  }
  PHObject* CloneMe() const override { return new PHHepMCGenEvent(*this); }

  HepMC::GenEvent* getEvent();
  const HepMC::GenEvent* getEvent() const;

  //! embedding ID for the event
  //! positive ID is the embedded event of interest, e.g. jetty event from pythia
  //! negative IDs are backgrounds, .e.g out of time pile up collisions
  //! Usually, ID = 0 means the primary Au+Au collision background
  int get_embedding_id() const { return _embedding_id; }

  //! embedding ID for the event
  //! positive ID is the embedded event of interest, e.g. jetty event from pythia
  //! negative IDs are backgrounds, .e.g out of time pile up collisions
  //! Usually, ID = 0 means the primary Au+Au collision background
  void set_embedding_id(int id) { _embedding_id = id; }

  //! whether this event has been processed in Geant4 simulation
  bool is_simulated() const { return _isSimulated; }

  //! whether this event has been processed in Geant4 simulation
  void is_simulated(bool v) { _isSimulated = v; }


  //! collision vertex position in the Hall coordinate system, use PHENIX units of cm, ns
  const HepMC::FourVector& get_collision_vertex() const { return _collisionVertex; }

  //! collision vertex position in the Hall coordinate system, use PHENIX units of cm, ns
  void set_collision_vertex(const HepMC::FourVector& v) { _collisionVertex = v; }

  //! boost beta vector for Lorentz Transform, part of composition of a LorentzRotation to translate from hepmc event frame to lab frame
  virtual const HepMC::ThreeVector& get_boost_beta_vector() const
  {
    PHOOL_VIRTUAL_WARNING;
    static HepMC::ThreeVector dummy_vec(0, 0, 0);
    return dummy_vec;
  }

  //! boost beta vector for Lorentz Transform, part of composition of a LorentzRotation to translate from hepmc event frame to lab frame
  virtual void set_boost_beta_vector(const HepMC::ThreeVector&) { PHOOL_VIRTUAL_WARNING; }

  //! rotation axis vector, part of composition of a LorentzRotation to translate from hepmc event frame to lab frame
  virtual const HepMC::ThreeVector& get_rotation_vector() const
  {
    PHOOL_VIRTUAL_WARNING;
    static HepMC::ThreeVector dummy_vec(0, 0, 1);
    return dummy_vec;
  }

  //! rotation axis vector, part of composition of a LorentzRotation to translate from hepmc event frame to lab frame
  virtual void set_rotation_vector(const HepMC::ThreeVector&) { PHOOL_VIRTUAL_WARNING; }

  //! rotation angle, part of composition of a LorentzRotation to translate from hepmc event frame to lab frame
  virtual double get_rotation_angle() const
  {
    PHOOL_VIRTUAL_WARNING;
    return 0;
  }

  //! rotation angle, part of composition of a LorentzRotation to translate from hepmc event frame to lab frame
  virtual void set_rotation_angle(const double) { PHOOL_VIRTUAL_WARNING; }

  //!LorentzRotation to translate from hepmc event frame to lab frame
  virtual CLHEP::HepLorentzRotation get_LorentzRotation_EvtGen2Lab() const { return CLHEP::HepLorentzRotation::IDENTITY; }

  //!LorentzRotation to translate from lab frame to hepmc event frame
  virtual CLHEP::HepLorentzRotation get_LorentzRotation_Lab2EvtGen() const { return CLHEP::HepLorentzRotation::IDENTITY; }


  //! reaction plane angles thrown by hijing flowAfterburner
  virtual float get_flow_psi(unsigned int /* n */) const
  {
    PHOOL_VIRTUAL_WARNING;
    return 0.0f;
  }

  // ! get the reaction plane angle psi_n[n]
  virtual const FlowAfterburner_PsiMap& get_flow_psi_map() const 
  { 
    PHOOL_VIRTUAL_WARNING;
    static FlowAfterburner_PsiMap dummy_map;
    return dummy_map; 
  }

  //! set the reaction plane angle psi_n[n]
  virtual void set_flow_psi(unsigned int /*n*/, float /*psi*/) { PHOOL_VIRTUAL_WARNING; }

  virtual void set_flow_psi_map(const FlowAfterburner_PsiMap& /*psi_map*/)
  {
    PHOOL_VIRTUAL_WARNING;
  }


  //! host an HepMC event
  bool addEvent(HepMC::GenEvent* evt);
  bool addEvent(HepMC::GenEvent& evt);
  bool swapEvent(HepMC::GenEvent*& evt);
  void clearEvent();

  //! move the collision vertex position in the Hall coordinate system, use PHENIX units of cm, ns
  virtual void moveVertex(double x, double y, double z, double t = 0);

  // the number of entries in the array of particles
  virtual int size(void) const;
  virtual int vertexSize(void) const;

  void print(std::ostream& os = std::cout) const;

  void PrintEvent();

 protected:
  //! \brief Embedding ID for this generated event
  //! positive ID is the embedded event of interest, e.g. jetty event from pythia
  //! negative IDs are backgrounds, .e.g out of time pile up collisions
  //! Usually, ID = 0 means the primary Au+Au collision background
  int _embedding_id{0};

  //! whether this event has been processed in Geant4 simulation
  bool _isSimulated{false};

  //! collision vertex position in the Hall coordinate system, use PHENIX units of cm, ns
  HepMC::FourVector _collisionVertex{0,0,0,0};

  //! The HEP MC record from event generator. Note the units are recorded in GenEvent
  HepMC::GenEvent* _theEvt{};

  ClassDefOverride(PHHepMCGenEvent, 5)
};

#endif  // PHHEPMC_PHHEPMCEVENT_H
