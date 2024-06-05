#include "common.hpp"
#include "InterEnergyGrid.hpp"
// #include "EnergyCalculator.hpp"
#include "AtomInterEnergyGrid.hpp"
#include "Fragment.hpp"
#include "Point3d.hpp"
#include "infile_reader.hpp"

#ifndef FRAGMENT_ENERGY_GRID_H_
#define FRAGMENT_ENERGY_GRID_H_
namespace fragdock {
  typedef std::pair<fragdock::Vector3d, fltype> FragPose;

  class FragmentInterEnergyGrid /* : public InterEnergyGrid */ {
    InterEnergyGrid grid;
    // void parse(const std::string& filename, int rot_size);
  public:
    int frag_idx;
    FragmentInterEnergyGrid() { frag_idx = -1; }
    // FragmentInterEnergyGrid(int frag_id, const std::string& filename, int rot_size)
    // : frag_id(frag_id) { parse(filename, rot_size); }
    FragmentInterEnergyGrid(const Fragment& orig_frag,
                       const std::vector<Vector3d>& rot_angles,
                       const std::vector<AtomInterEnergyGrid>& atom_grids,
                       const InterEnergyGrid& distance_grid);
    ~FragmentInterEnergyGrid() {}
    const InterEnergyGrid& getGrid() const { return grid; }
    // void writeFile(const std::string& filename) const;

    /**
     * Get the energy of the grid at the given position.
     * @param[in] x x coordinate
     * @param[in] y y coordinate
     * @param[in] z z coordinate
     * @param[in] scale a scale factor for the energy
    */
    fltype getSearchScore(int x, int y, int z, const fltype scale=1) const { return getGrid().getInterEnergy(x, y, z) / scale; }
    fltype getSearchScore(int x, int y, int z) const { return getGrid().getInterEnergy(x, y, z); }

    /**
     * 
     * @param[in] priority priority configuration
     * @param[in] scale a scale factor for the energy
    */
    std::vector<FragPose> getBestPoses(const format::Priority& priority, const fltype scale) const;
    std::vector<FragPose> getBestPoses(const format::Priority& priority) const { return getBestPoses(priority, 1); }
  };
}
#endif

