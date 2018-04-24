#include "types_imu.h"

namespace g2o {

  VertexImu::VertexImu() : BaseVertex<15, ImuState>()
  {
  }


  bool VertexImu::read(std::istream& is)
  {
    return true;
  }

  bool VertexImu::write(std::ostream& os) const
  {
    return true;
  }



} // end namespace
