#ifndef TRACKRECO_PHSPACECHARGERECONSTRUCTION_H
#define TRACKRECO_PHSPACECHARGERECONSTRUCTION_H

#include <fun4all/SubsysReco.h>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <array>

// forward declaration
class SvtxTrack;
class SvtxTrackMap;
class TrkrCluster;
class TrkrClusterContainer;

class PHSpaceChargeReconstruction: public SubsysReco
{
  public:

  /// constructor
  PHSpaceChargeReconstruction( const std::string& = "PHSPACECHARGERECONSTRUCTION" );

  /// global initialization
  virtual int Init(PHCompositeNode*);

  /// run initialization
  virtual int InitRun(PHCompositeNode*);

  /// event processing
  virtual int process_event(PHCompositeNode*);

  /// end of processing
  virtual int End(PHCompositeNode*);

  private:

  /// load nodes
  int load_nodes( PHCompositeNode* );

  /// process tracks
  void process_tracks();

  /// process track
  void process_track( SvtxTrack* );

  /// calculate distortions
  void calculate_distortions();

  /// get cell from z, r and phi index
  int get_cell( int iz, int ir, int iphi ) const;

  /// get relevant cell for a given cluster
  int get_cell( TrkrCluster* ) const;

  /// output file
  std::string m_outputfile = "PHSpaceChargeReconstruction.root";

  ///@name grid size
  //@{
  static constexpr int m_z_size = 1;
  static constexpr int m_phi_size = 1;
  static constexpr int m_r_size = 12;
//   static constexpr int m_z_size = 50;
//   static constexpr int m_phi_size = 72;
//   static constexpr int m_r_size = 48;
  static constexpr int m_total_size = m_z_size*m_phi_size*m_r_size;
  //@}

  // shortcut for relevant eigen matrices
  static constexpr int m_ncoord = 3;
  using matrix_t = Eigen::Matrix<float, m_ncoord, m_ncoord >;
  using column_t = Eigen::Matrix<float, m_ncoord, 1 >;

  /// left hand side matrices for distortion inversions
  std::array<matrix_t, m_total_size> m_lhs;

  /// right hand side matrices for distortion inversions
  std::array<column_t, m_total_size> m_rhs;

  /// keep track of how many clusters are used per cell
  std::array<int, m_total_size> m_cluster_count;

  ///@name nodes
  //@{
  SvtxTrackMap* m_track_map = nullptr;
  TrkrClusterContainer* m_cluster_map = nullptr;
  //@}

};

#endif
