#include "common.hpp"
#include "Fragment.hpp"
#include "FragmentInterEnergyGrid.hpp"
#include "Point3d.hpp"

#ifndef OPEN_DX_H_
#define OPEN_DX_H_
namespace opendx {
  class OpenDx {
  private:
    // fragdock::Fragment fragment;
    fragdock::FragmentInterEnergyGrid fgrid;
    int heavy_num;
  public:
    OpenDx(const fragdock::FragmentInterEnergyGrid& fgrid, const int heavy_num) : fgrid(fgrid), heavy_num(heavy_num) {}
    void writeDx(std::ofstream& ofs, const fltype scale) const;
    void write(const std::string& filename, const fltype scale) const;
  };

  template<typename T>
  std::string point3d_to_string(const fragdock::Point3d<T>& p);
}
#endif