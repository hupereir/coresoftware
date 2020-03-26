#include "PHSpaceChargeReconstruction.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>

#include <TFile.h>
#include <TGraphErrors.h>

namespace
{

  /// square
  template<class T> T square( T x ) { return x*x; }

  /// calculate delta_phi between -pi and pi
  template< class T>
    T delta_phi( const T& phi )
  {
    if( phi > M_PI ) return phi - 2*M_PI;
    else if( phi <= -M_PI ) return phi + 2*M_PI;
    else return phi;
  }

  /// radius
  template<class T> T get_r( T x, T y ) { return std::sqrt( square(x) + square(y) ); }

  /// phi
  template<class T> T get_phi( T x, T y ) { return std::atan2( y, x ); }

  /// get rphi error of a given track state
  float get_rphi_error( SvtxTrackState* state )
  {
    using matrix_t = Eigen::Matrix<float, 3, 3>;
    matrix_t covar;
    for (unsigned int i = 0; i < 3; ++i)
      for (unsigned int j = 0; j < 3; ++j)
    { covar(i,j) = state->get_error(i, j); }

    const auto phi = -get_phi( state->get_x(), state->get_y() );
    const auto cosphi = std::cos( phi );
    const auto sinphi = std::sin( phi );
    matrix_t rotation;
    rotation(0,0) = cosphi;
    rotation(0,1) = -sinphi;
    rotation(0,2) = 0;
    rotation(1,0) = sinphi;
    rotation(1,1) = cosphi;
    rotation(1,2) = 0;
    rotation(2,0) = 0;
    rotation(2,1) = 0;
    rotation(2,2) = 1;

    const auto transformed = rotation*covar*rotation.transpose();
    return std::sqrt(transformed(1,1));
  }
}

//_____________________________________________________________________
PHSpaceChargeReconstruction::PHSpaceChargeReconstruction( const std::string& name ):
  SubsysReco( name)
{}

//_____________________________________________________________________
int PHSpaceChargeReconstruction::Init(PHCompositeNode* topNode )
{
  // initialize all arrays
  for( auto&& matrix:m_lhs ) { matrix = matrix_t::Zero(); }
  for( auto&& column:m_rhs ) { column = column_t::Zero(); }
  for( auto&& count:m_cluster_count ) { count = 0; }

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int PHSpaceChargeReconstruction::InitRun(PHCompositeNode* )
{ return Fun4AllReturnCodes::EVENT_OK; }

//_____________________________________________________________________
int PHSpaceChargeReconstruction::process_event(PHCompositeNode* topNode)
{
  // load nodes
  const auto res =  load_nodes(topNode);
  if( res != Fun4AllReturnCodes::EVENT_OK ) return res;

  process_tracks();
  return Fun4AllReturnCodes::EVENT_OK;

}

//_____________________________________________________________________
int PHSpaceChargeReconstruction::End(PHCompositeNode* )
{
  calculate_distortions();
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int PHSpaceChargeReconstruction::load_nodes( PHCompositeNode* topNode )
{
  // get necessary nodes
  m_track_map = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  m_cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
void PHSpaceChargeReconstruction::process_tracks()
{
  if( !( m_track_map && m_cluster_map ) ) return;
  std::cout << PHWHERE << " tracks: " << m_track_map->size() << " clusters: " << m_cluster_map->size() << std::endl;

  for( auto iter = m_track_map->begin(); iter != m_track_map->end(); ++iter )
  { process_track( iter->second ); }
}

//_____________________________________________________________________
void PHSpaceChargeReconstruction::process_track( SvtxTrack* track )
{

  // running track state
  auto state_iter = track->begin_states();

  // loop over clusters
  for( auto key_iter = track->begin_cluster_keys(); key_iter != track->end_cluster_keys(); ++key_iter )
  {

    const auto& cluster_key = *key_iter;
    auto cluster = m_cluster_map->findCluster( cluster_key );
    if( !cluster )
    {
      std::cout << PHWHERE << " unable to find cluster for key " << cluster_key << std::endl;
      continue;
    }

    // should make sure that cluster belongs to TPC

    // cluster r, phi and z
    const auto cluster_r = get_r( cluster->getX(), cluster->getY() );
    const auto cluster_phi = get_phi( cluster->getX(), cluster->getY() );
    const auto cluster_z = cluster->getZ();

    // cluster errors
    const auto cluster_rphi_error = cluster->getRPhiError();
    const auto cluster_z_error = cluster->getZError();

    // find track state that is the closest to cluster
    /* this assumes that both clusters and states are sorted along r */
    float dr_min = -1;
    for( auto iter = state_iter; iter != track->end_states(); ++iter )
    {
      const auto dr = std::abs( cluster_r - get_r( iter->second->get_x(), iter->second->get_y() ) );
      if( dr_min < 0 || dr < dr_min )
      {
        state_iter = iter;
        dr_min = dr;
      } else break;
    }

    // track r, phi and z
    auto state = state_iter->second;
    const auto track_phi = get_phi(  state->get_x(), state->get_y() );
    const auto track_z = state->get_z();

    // track errors
    const auto track_rphi_error = get_rphi_error( state );
    const auto track_z_error = std::sqrt( state->get_error(2,2) );

    //       std::cout << PHWHERE << " cluster errors: " << cluster_rphi_error << ", " << cluster_z_error << std::endl;
    //       std::cout << PHWHERE << " track errors: " << track_rphi_error << ", " << track_z_error << std::endl;

    /*
    remove clusters with too small errors since they are likely pathological
    and have a large contribution to the chisquare
    TODO: make these cuts configurable
    */
    if( cluster_rphi_error < 0.015 ) continue;
    if( cluster_z_error < 0.05 ) continue;

    // also cut on track errors
    // warning: smaller errors are probably needed when including outer tracker
    if( track_rphi_error < 0.015 ) continue;
    if( track_z_error < 0.1 ) continue;

    //       // get residual errors squared
    //       const auto erp = square(track_rphi_error) + square(cluster_rphi_error);
    //       const auto ez = square(track_z_error) + square(cluster_z_error);

    // get residual errors squared
    const auto erp = 1;
    const auto ez = 1;

    // sanity check
    if( std::isnan( erp ) ) continue;
    if( std::isnan( ez ) ) continue;

    // get residuals
    const auto drp = cluster_r*delta_phi( track_phi - cluster_phi );
    const auto dz = track_z - cluster_z;

    // get track angles
    const auto cosphi( std::cos( track_phi ) );
    const auto sinphi( std::sin( track_phi ) );
    const auto track_pphi = -state->get_px()*sinphi + state->get_py()*cosphi;
    const auto track_pr = state->get_px()*cosphi + state->get_py()*sinphi;
    const auto track_pz = state->get_pz();
    const auto talpha = -track_pphi/track_pr;
    const auto tbeta = -track_pz/track_pr;

    // check against limits
    // TODO: make this configurable
    static constexpr float max_talpha = 0.6;
    static constexpr float max_residual = 5;
    if( std::abs( talpha ) > max_talpha ) continue;
    if( std::abs( drp ) > max_residual ) continue;

    // get cell
    const auto i = get_cell( cluster );

    if( i < 0 || i >= m_total_size ) continue;

    m_lhs[i](0,0) += 1./erp;
    m_lhs[i](0,1) += 0;
    m_lhs[i](0,2) += talpha/erp;

    m_lhs[i](1,0) += 0;
    m_lhs[i](1,1) += 1./ez;
    m_lhs[i](1,2) += tbeta/ez;

    m_lhs[i](2,0) += talpha/erp;
    m_lhs[i](2,1) += tbeta/ez;
    m_lhs[i](2,2) += square(talpha)/erp + square(tbeta)/ez;

    m_rhs[i](0,0) += drp/erp;
    m_rhs[i](1,0) += dz/ez;
    m_rhs[i](2,0) += talpha*drp/erp + tbeta*dz/ez;

    ++m_cluster_count[i];

  }

}

//_____________________________________________________________________
void  PHSpaceChargeReconstruction::calculate_distortions()
{
  // calculate distortions in each volume elements
  std::array<column_t, m_total_size> delta;
  std::array<matrix_t, m_total_size> cov;
  for( int i = 0; i < m_total_size; ++i )
  {

    std::cout
      << PHWHERE << " lhs " << i << std::endl
      << m_lhs[i] << std::endl;

    cov[i] = m_lhs[i].inverse();
    delta[i] = m_lhs[i].partialPivLu().solve( m_rhs[i] );
  }

  // create tgraphs
  using TGraphPointer = std::unique_ptr<TGraphErrors>;
  std::array<TGraphPointer, m_ncoord> tg;
  for( int i = 0; i < m_ncoord; ++i )
  {
    tg[i].reset( new TGraphErrors() );

    constexpr int iz = 0;
    constexpr int iphi = 0;
    tg[i]->SetName( Form( "tg_%i_%i_%i", iz, iphi, i ) );

    for( int ir = 0; ir < m_r_size; ++ir )
    {
      // get radius from index (see ::get_cell)
      static constexpr float r_min = 30;
      static constexpr float r_max = 80;
      const float r = r_min + (0.5+ir) *(r_max-r_min)/m_r_size;

      int index = get_cell( iz, ir, iphi );


      tg[i]->SetPoint( ir, r, delta[index](i,0) );
      tg[i]->SetPointError( ir, 0, std::sqrt(cov[index](i,i)) );
    }
  }

  // save to root file
  {
    std::unique_ptr<TFile> outputfile( TFile::Open( m_outputfile.c_str(), "RECREATE" ) );
    outputfile->cd();
    for( int i = 0; i < m_ncoord; ++i )
    { tg[i]->Write(); }
    outputfile->Close();
  }
}

//_____________________________________________________________________
int PHSpaceChargeReconstruction::get_cell( int iz, int ir, int iphi ) const
{
  if( ir < 0 || ir >= m_r_size ) return -1;
  if( iphi < 0 || iphi >= m_phi_size ) return -1;
  if( iz < 0 || iz >= m_z_size ) return -1;
  return iz + m_z_size*( ir + m_r_size*iphi );
}

//_____________________________________________________________________
int PHSpaceChargeReconstruction::get_cell( TrkrCluster* cluster ) const
{
  // radius
  // todo: define rbin as a function of layer index rather than position
  const auto cluster_r = get_r( cluster->getX(), cluster->getY() );
  static constexpr float r_min = 30;
  static constexpr float r_max = 80;
  const int ir = m_r_size*(cluster_r-r_min)/(r_max-r_min);

  // azimuth
  auto cluster_phi = get_phi( cluster->getX(), cluster->getY() );
  if( cluster_phi >= M_PI ) cluster_phi -= 2*M_PI;
  const int iphi = m_phi_size*(cluster_phi + M_PI)/(2.*M_PI);

  // z
  const auto cluster_z = cluster->getZ();
  static constexpr float z_min = -212/2;
  static constexpr float z_max = 212/2;
  const int iz = m_z_size*(cluster_z-z_min)/(z_max-z_min);

  return get_cell( iz, ir, iphi );
}
