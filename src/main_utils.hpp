#include "common.hpp"
#include "Point3d.hpp"
#include "Vector3d.hpp"
#include "InterEnergyGrid.hpp"
#include "Molecule.hpp"
#include "Fragment.hpp"
#include "OBMol.hpp"
#include "infile_reader.hpp"

#include <cmath>

#ifndef MAIN_UTILS_H_
#define MAIN_UTILS_H_

namespace main_utils {
  typedef format::DockingConfiguration::FragmentEfficiency FragEffi;

	void makeFolder(std::string folderName);
	std::string getDate();
	fragdock::Point3d<int> to_score_num(int k, 
                                      const fragdock::Point3d<int>& score_num, 
                                      const fragdock::Point3d<int>& search_num, 
                                      const fragdock::Point3d<int>& ratio);
	fragdock::Vector3d operator/(fragdock::Vector3d v, fragdock::Point3d<fltype> p);
	fragdock::Point3d<int> round(const fragdock::Vector3d& v);
	
  /**
   * Convert OpenBabel molecules to Fragment.
   * @param obmols a vector of OpenBabel molecules
   * @return a vector of Fragment
  */
	std::vector<fragdock::Fragment> convert_fragments(std::vector<OpenBabel::OBMol>& obmols);

  /**
   * Calculate the fragment efficiency scale.
   * @param heavy_num the number of heavy atoms in the fragment
   * @param efficiency the efficiency type
   * @return scale of the score
  */
	const fltype calcFragmentEfficiency(const fltype heavy_num, const FragEffi efficiency);

  /**
   * Get the description of the option.
   * @param main the main description
   * @param options a vector of option descriptions (option name, option description)
   * @return the description of the option
  */
	std::string option_desc(const std::string main, const std::vector<std::pair<std::string, std::string> >& options);
	
  /**
   * Struct for displaying progress bar.
  */
  struct progress_bar {
    const int bar_size = 100;
    const int total;
    const std::string desc;
    const bool show_bar = false; // default is not to show progress bar
    progress_bar(int total, bool show_bar) : total(total), desc(""), show_bar(show_bar) { display(0); }
    progress_bar(int total, std::string desc, bool show_bar) : total(total), desc(desc+":"), show_bar(show_bar) { display(0); }
    void display(int progress) const;
    void clear() const;
  };
} // namespace main_utils

namespace fragdock {
	InterEnergyGrid makeDistanceGrid(const Point3d<fltype>& center,
                                   const Point3d<fltype>& pitch,
                                   const Point3d<int>& num,
                                   const Molecule& receptor_mol);
} // namespace fragdock

#endif