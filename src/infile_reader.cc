#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include "infile_reader.hpp"


namespace format {
  DockingConfiguration ParseInFile(const char *filename){
    std::ifstream ifs(filename);
    if (ifs.fail()){
      std::cerr << "opening grid file failed:" << filename << std::endl;
      abort();
    }
    DockingConfiguration conf;
    std::string buffer;
    while(!ifs.eof()){
      std::getline(ifs, buffer);
      if (boost::algorithm::starts_with(buffer, "OUTERBOX ")) {
        std::string str = buffer.substr(9);
        std::vector<std::string> vals;
        boost::algorithm::split(vals, str, boost::algorithm::is_any_of(","));
        boost::algorithm::trim(vals[0]);
        boost::algorithm::trim(vals[1]);
        boost::algorithm::trim(vals[2]);
        conf.grid.outer_width = fragdock::Point3d<fltype>(
          boost::lexical_cast<fltype>(vals[0]),
          boost::lexical_cast<fltype>(vals[1]),
          boost::lexical_cast<fltype>(vals[2]));
      }
      else if (boost::algorithm::starts_with(buffer, "BOX_CENTER ")) {
        std::string str = buffer.substr(11);
        std::vector<std::string> vals;
        boost::algorithm::split(vals, str, boost::algorithm::is_any_of(","));
        boost::algorithm::trim(vals[0]);
        boost::algorithm::trim(vals[1]);
        boost::algorithm::trim(vals[2]);
        conf.grid.center = fragdock::Point3d<fltype>(
          boost::lexical_cast<fltype>(vals[0]),
          boost::lexical_cast<fltype>(vals[1]),
          boost::lexical_cast<fltype>(vals[2]));
      }
      else if (boost::algorithm::starts_with(buffer, "SCORING_PITCH ")) {
        std::string str = buffer.substr(14);
        std::vector<std::string> vals;
        boost::algorithm::split(vals, str, boost::algorithm::is_any_of(","));
        boost::algorithm::trim(vals[0]);
        boost::algorithm::trim(vals[1]);
        boost::algorithm::trim(vals[2]);
        conf.grid.score_pitch = fragdock::Point3d<fltype>(
          boost::lexical_cast<fltype>(vals[0]),
          boost::lexical_cast<fltype>(vals[1]),
          boost::lexical_cast<fltype>(vals[2]));
      }
      else if (boost::algorithm::starts_with(buffer, "RECEPTOR ")) {
        conf.receptor_file = buffer.substr(9);
      }
      else if (boost::algorithm::starts_with(buffer, "OUTPUT ")) {
        conf.output_file = buffer.substr(7);
      }
      else if (boost::algorithm::starts_with(buffer, "LOG ")) {
        conf.log_file = buffer.substr(4);
      }
      else if (boost::algorithm::starts_with(buffer, "GRID_FOLDER ")) {
        conf.grid_folder = buffer.substr(12);
      }
      else if (boost::algorithm::starts_with(buffer, "ROTANGS ")) {
        conf.rotangs_file = buffer.substr(8);
      }
      else if (boost::algorithm::starts_with(buffer, "FRAGMENTS ")) {
        // FRAGMENTS [FRAGMENTS]
        conf.fragments_file = buffer.substr(10);
      }
      else if (boost::algorithm::starts_with(buffer, "DX_FOLDER ")) {
        conf.dx_folder = buffer.substr(10);
      }
      else if (boost::algorithm::starts_with(buffer, "OB_ERROR_LOG ")) { // for debug
        std::string str = buffer.substr(13);
        boost::algorithm::trim(str);
        std::transform(str.begin(), str.end(), str.begin(), ::tolower);
        conf.ob_err_log = (str != "false");
      }
      else if (boost::algorithm::starts_with(buffer, "FRAG_EFFICIENCY ")) {
        std::string str = buffer.substr(16);
        boost::algorithm::trim(str);
        std::transform(str.begin(), str.end(), str.begin(), ::tolower);
        conf.efficiency = conf.parseFragmentEfficiency(str);
      }
      else if (boost::algorithm::starts_with(buffer, "PROMISING_POSE ")) {
        std::string str = buffer.substr(10);
        boost::algorithm::trim(str);
        conf.priority.promising_pose = boost::lexical_cast<int64_t>(str);
      }
      else if (boost::algorithm::starts_with(buffer, "CLUSTER_SIZE ")) {
        std::string str = buffer.substr(13);
        boost::algorithm::trim(str);
        conf.priority.cluster_size = boost::lexical_cast<fltype>(str);
      }
      else if (boost::algorithm::starts_with(buffer, "DISTANCE_WIDTH ")) {
        std::string str = buffer.substr(15);
        boost::algorithm::trim(str);
        conf.priority.distance_width = boost::lexical_cast<fltype>(str);
      }
      else if (boost::algorithm::starts_with(buffer, "FUNCTION ")) {
        std::string str = buffer.substr(9);
        boost::algorithm::trim(str);
        std::transform(str.begin(), str.end(), str.begin(), ::tolower);
        conf.priority.function = conf.parsePriorityFunctin(str);
      }
    }

    return conf;
  }
} // namespace format