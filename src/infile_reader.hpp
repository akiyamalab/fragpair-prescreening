#include "common.hpp"
#include "InterEnergyGrid.hpp"
#include <string>
#include <vector>

#ifndef INFILE_READER_H_
#define INFILE_READER_H_

namespace format {
  struct SearchGrid {
    fragdock::Point3d<fltype> center, outer_width, inner_width, search_pitch, score_pitch;
  };

  struct Priority {
    enum Function {MIN, SIMPLE, FRAC1, FRAC2};

    int64_t promising_pose = 40;
    fltype cluster_size = 1.0;
    fltype distance_width = 0.1;
    Function function = Function::FRAC2;

    const std::string getFunctionString() const {
      switch (function) {
        case Priority::Function::MIN   : return "FUNCTION_MIN";
        case Priority::Function::SIMPLE: return "FUNCTION_SIMPLE";
        case Priority::Function::FRAC1 : return "FUNCTION_FRAC1";
        case Priority::Function::FRAC2 : return "FUNCTION_FRAC2";
      }
      assert(0);
    }
  };

  struct DockingConfiguration {
    enum ReuseStrategy {OFFLINE, ONLINE, NONE};
    enum FragmentEfficiency {OFF, VOLUME, SURFACE};

    SearchGrid grid;
    Priority priority;
    std::vector<std::string> ligand_files;
    std::string receptor_file, output_file;
    std::string log_file, grid_folder, dx_folder;
    std::string rotangs_file;
    std::string fragments_file;
    ReuseStrategy reuse_grid = ReuseStrategy::OFFLINE;
    FragmentEfficiency efficiency = FragmentEfficiency::OFF;
    bool reorder = true;
    bool ob_err_log = true; // for debug
    bool show_bar = false; // for debug in local

    const std::string getReuseGridString() {
      switch (reuse_grid) {
        case ReuseStrategy::OFFLINE: return "REUSE_OFFLINE";
        case ReuseStrategy::ONLINE : return "REUSE_ONLINE";
        case ReuseStrategy::NONE   : return "REUSE_NONE";
      }
      assert(0);
    }
    const std::string getFragmentEfficiencyString() {
      switch (efficiency) {
        case FragmentEfficiency::OFF    : return "EFFICIENCY_OFF";
        case FragmentEfficiency::VOLUME : return "EFFICIENCY_VOLUME";
        case FragmentEfficiency::SURFACE: return "EFFICIENCY_SURFACE";
      }
      assert(0);
    }
    const std::string getPriorityFunctionString() { return priority.getFunctionString(); }
    const FragmentEfficiency parseFragmentEfficiency(const std::string &str) {
      if (str == "volume") {
          return DockingConfiguration::FragmentEfficiency::VOLUME;
        }
        else if (str == "surface") {
          return DockingConfiguration::FragmentEfficiency::SURFACE;
        }
        else {
          return DockingConfiguration::FragmentEfficiency::OFF;
        }
    }
    const Priority::Function parsePriorityFunctin(const std::string &str) {
      if (str == "min") {
        return Priority::Function::MIN;
      }
      else if (str == "simple") {
        return Priority::Function::SIMPLE;
      }
      else if (str == "frac1") {
        return Priority::Function::FRAC1;
      }
      else {
        return Priority::Function::FRAC2;
      }
    }
    void checkConfigValidity() const;
  };

  DockingConfiguration ParseInFile(const char *filename);
} // namespace format

#endif