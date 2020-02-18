#ifndef G4EVAL_SIMEVALUATOR_HP_H
#define G4EVAL_SIMEVALUATOR_HP_H

#include <fun4all/SubsysReco.h>
#include <phool/PHObject.h>

#include <TClonesArray.h>

class PHG4TruthInfoContainer;

// vertex information
class VertexStruct: public TObject
{
  public:

  virtual const char* GetName() const
  { return "VertexStruct"; }

  void Clear()
  {
    _x = 0;
    _y = 0;
    _z = 0;
    _t = 0;
    _is_main_vertex = false;
  }

  float _x = 0;
  float _y = 0;
  float _z = 0;
  float _t = 0;

  bool _is_main_vertex = false;

  ClassDef(VertexStruct,1)

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
    Container();

    /// copy constructor
    explicit Container(const Container &) = delete;

    /// assignment operator
    Container& operator = ( const Container& ) = delete;

    /// destructor
    virtual ~Container();

    /// reset
    virtual void Reset();

    TClonesArray* primary_vertex_list()
    { return _vertex_list; }
    //@}

    private:

    //* background vertex list
    TClonesArray* _vertex_list = nullptr;

    ClassDef(Container,1)

  };

  private:

  /// load nodes
  int load_nodes( PHCompositeNode* );

  /// fill vertices
  void fill_vertices();

  /// print vertices
  void print_vertices();

  //* data container
  Container* _container = nullptr;
  int _vertex_count = 0;

  //* truth information
  PHG4TruthInfoContainer* _g4truthinfo = nullptr;

};

#endif  // G4EVAL_SIMEVALUATOR_HP_H
