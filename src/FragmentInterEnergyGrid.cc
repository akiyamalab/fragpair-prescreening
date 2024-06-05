#include "FragmentInterEnergyGrid.hpp"

namespace fragdock {
  FragmentInterEnergyGrid::FragmentInterEnergyGrid(const Fragment& orig_frag,
                                         const std::vector<Vector3d>& rot_angles,
                                         const std::vector<AtomInterEnergyGrid>& atom_grids,
                                         const InterEnergyGrid& distance_grid) : frag_idx(orig_frag.getIdx()) {
    using namespace std;
    if(atom_grids.empty()) {
      cerr << "atom_grids is empty" << endl;
      return;
    }
    assert(orig_frag.getCenter().abs() < EPS);
    const Point3d<int>& num = atom_grids[0].getNum();
    grid = InterEnergyGrid(atom_grids[0].getCenter(), atom_grids[0].getPitch(), num, LIMIT_ENERGY);
    int rotsz = rot_angles.size();
    if (orig_frag.size() <= 1) rotsz = 1;

    fltype radius = orig_frag.getRadius();

    for(int rot_id = 0; rot_id < rotsz; rot_id++) {
      Fragment frag = orig_frag;
      frag.rotate(rot_angles[rot_id]);
      const vector<Atom>& atoms = frag.getAtoms();
      for(int x = 0; x < num.x; x++) {
        for(int y = 0; y < num.y; y++) {
          for(int z = 0; z < num.z; z++) {
            if (rot_id & 7) {
              fltype collision = 2;
              fltype far = 6;
              fltype dist = distance_grid.getInterEnergy(x, y, z);
              if (dist < collision) {
                // collision
                continue;
              }
              if (dist > radius + far) {
                // too far
                continue;
              }
            }

            fltype sum = 0.0;
            for (const Atom& atom : atoms) {
              if (atom.getXSType() == XS_TYPE_H) continue;

              const AtomInterEnergyGrid& agrid = atom_grids[atom.getXSType()];
              // atom += agrid.convert(x, y, z);
              fltype diff = agrid.getInterEnergy(atom + agrid.convert(x, y, z));
              // atom -= agrid.convert(x, y, z);

              sum += diff;
              if(sum >= LIMIT_ENERGY) {
                break;
              }
            }
            if (sum < grid.getInterEnergy(x, y, z)) {
              grid.setInterEnergy(x, y, z, sum);
            }
          }
        }
      }
    }
  }

  std::vector<FragPose> FragmentInterEnergyGrid::getBestPoses(const format::Priority& priority, const fltype scale) const {
    using namespace std;
    vector<FragPose> coord_scores;
    vector<FragPose> ret(priority.promising_pose);
    Point3d<int> score_num = getGrid().getNum();
    
    coord_scores.reserve(score_num.x * score_num.y * score_num.z);
    for (int x = 0; x < score_num.x; ++x) {
      for (int y = 0; y < score_num.y; ++y) {
        for (int z = 0; z < score_num.z; ++z) {
          // push back FragPose
          coord_scores.emplace_back(getGrid().convert(x, y, z), getSearchScore(x, y, z, scale));
        }
      }
    }

    // sort by score
    sort(coord_scores.begin(), coord_scores.end(), [](const FragPose& a, const FragPose& b) {
        return a.second < b.second;
      });

    int ret_size = 0;
    // get the best poses
    for (int i = 0; i < coord_scores.size(); ++i) {
      fltype min_dist = HUGE_VAL;
      // check if the distance is too close to the existing (ret) poses
      for (int r = 0; r < ret_size; ++r) {
        min_dist = min(min_dist, (ret[r].first - coord_scores[i].first).abs());
      }
      // if the distance is not too close, add the pose to the return value
      if (min_dist > priority.cluster_size) {
        ret[ret_size] = coord_scores[i];
        ret_size++;
      }

      if (ret_size == ret.size()) break;
    }

    return ret;
  }
}
