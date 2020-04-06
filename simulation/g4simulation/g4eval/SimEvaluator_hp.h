#ifndef G4EVAL_SIMEVALUATOR_HP_H
#define G4EVAL_SIMEVALUATOR_HP_H

#include <fun4all/SubsysReco.h>
#include <phool/PHObject.h>

#include <map>
#include <set>
#include <vector>

class PHG4Hit;
class PHG4HitContainer;
class PHG4Particle;
class PHG4TruthInfoContainer;

// vertex information
class VertexStruct
{
  public:

  using List = std::vector<VertexStruct>;

  float _x = 0;
  float _y = 0;
  float _z = 0;
  float _t = 0;

  bool _is_main_vertex = false;

};

// vertex information
class ParticleStruct
{
  public:
  using List = std::vector<ParticleStruct>;

  int _charge = 0;
  int _pid = 0;
  int _embed = 0;
  bool _is_primary = false;
  int64_t _mask = 0;


  float _px = 0;
  float _py = 0;
  float _pz = 0;
  float _pt = 0;
  float _p = 0;
  float _eta = 0;

};


class SimEvaluator_hp : public SubsysReco
{
  public:

  /// constructor
  SimEvaluator_hp( const std::string& = "SimEvaluator_hp" );

  /// global initialization
  virtual int Init(PHCompositeNode*);

  /// run initialization
  virtual int InitRun(PHCompositeNode*);

  /// event processing
  virtual int process_event(PHCompositeNode*);

  /// end of processing
  virtual int End(PHCompositeNode*);

  class Container: public PHObject
  {

    public:

    /// constructor
    explicit Container() = default;

    /// copy constructor
    explicit Container(const Container &) = delete;

    /// assignment operator
    Container& operator = ( const Container& ) = delete;

    /// reset
    virtual void Reset();

    ///@name accessors
    //@{
    const VertexStruct::List& vertexList() const
    { return _vertex_list; }

    const ParticleStruct::List& particleList() const
    { return _particle_list; }

    //@}

    ///@name modifiers
    //@{

    void addVertex( const VertexStruct& vertex )
    { _vertex_list.push_back( vertex ); }

    void addParticle( const ParticleStruct& particle )
    { _particle_list.push_back( particle ); }

    void clearVertexList()
    { _vertex_list.clear(); }

    void clearParticleList()
    { _particle_list.clear(); }

    //@}

    private:

    //* vertex list
    VertexStruct::List _vertex_list;

    //* particles
    ParticleStruct::List _particle_list;

    ClassDef(Container,1)

  };

  private:

  /// load nodes
  int load_nodes( PHCompositeNode* );

  /// fill MC track map
  void fill_mc_track_map();

  /// fill vertices
  void fill_vertices();

  /// fill particles
  void fill_particles();

  /// print vertices
  void print_vertices();

  // get embedded id for given g4track
  int get_embed(PHG4Particle*) const;

  //* data container
  Container* _container = nullptr;

  PHG4HitContainer* _g4hits_tpc = nullptr;
  PHG4HitContainer* _g4hits_intt = nullptr;
  PHG4HitContainer* _g4hits_mvtx = nullptr;
  PHG4HitContainer* _g4hits_outertracker = nullptr;

  //* truth information
  PHG4TruthInfoContainer* _g4truthinfo = nullptr;

  // map trk_id to list of g4 hits
  using G4HitSet = std::set<PHG4Hit*>;
  using G4ParticleMap = std::map<int,G4HitSet>;
  G4ParticleMap _g4particle_map;

};

#endif  // G4EVAL_SIMEVALUATOR_HP_H
