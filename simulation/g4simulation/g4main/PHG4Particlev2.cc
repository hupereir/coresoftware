#include "PHG4Particlev2.h"

#include "PHG4Particle.h"  // for PHG4Particle

#include <cmath>

PHG4Particlev2::PHG4Particlev2(const std::string &name, const int pid, const double px, const double py, const double pz)
  : PHG4Particlev1(name, pid, px, py, pz)
{
}

PHG4Particlev2::PHG4Particlev2(const PHG4Particle *in)
  : PHG4Particlev1(in)
  , trkid(in->get_track_id())
  , vtxid(in->get_vtx_id())
  , parentid(in->get_parent_id())
  , primaryid(in->get_primary_id())
  , fe(in->get_e())
{
}

void PHG4Particlev2::identify(std::ostream &os) const
{
  if (!fname.empty())
  {
    os << "PHG4Particlev2 name: " << fname << ", ";
  }
  else
  {
    os << "PHG4Particlev2 name: missing, ";
  }

  os << "track id: " << trkid
     << ", vtxid: " << vtxid
     << ", parent id: " << parentid
     << ", primary id: " << primaryid
     << ", pid: " << fpid
     << ", px: " << fpx
     << ", py: " << fpy
     << ", pz: " << fpz
     << ", phi: " << atan2(fpy, fpx)
     << ", eta: " << -1 * log(tan(0.5 * acos(fpz / sqrt((fpx * fpx) + (fpy * fpy) + (fpz * fpz)))))
     << ", e: " << fe << std::endl;
  return;
}
