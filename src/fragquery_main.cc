#include "common.hpp"
#include "main_utils.hpp"
#include "utils.hpp"
#include "OBMol.hpp"
#include "Point3d.hpp"
#include "infile_reader.hpp"
#include "log_writer_stream.hpp"
#include "AtomInterEnergyGrid.hpp"
#include "FragmentInterEnergyGrid.hpp"
#include "EnergyCalculator.hpp"
#include "QueryPriority.hpp"

#include <iostream>
#include <iomanip>
#include <string>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <cstdlib>
#include <chrono>
#include <stdexcept>

namespace {
  format::DockingConfiguration parseArgs(int argc, char **argv){
    using namespace boost::program_options;
    options_description options("Options");
    options_description hidden("Hidden options");
    positional_options_description pos_desc;
    hidden.add_options()
      ("conf-file", value<std::string>(), "configuration file");
    pos_desc.add("conf-file", 1);
    options.add_options()
      ("help,h", "show help")
      ("output,o", value<std::string>(), "output file (.txt file)")
      ("receptor,r", value<std::string>(), "receptor file (.pdb file)")
      ("grid,g", value<std::string>(), "grid folder")
      ("log", value<std::string>(), "log file")
      ("fragments,f", value<std::string>(), "fragments file (.sdf file)")
      ("efficiency,e", value<std::string>(), main_utils::option_desc("efficiency mode", {std::make_pair("\"volume\"", "volume efficiency"), 
                                                                                         std::make_pair("\"surface\"", "surface efficiency"), 
                                                                                         std::make_pair("default/else", "off")}
                                                                    ).c_str())
      ("oberror", value<std::string>(), "ob error log mode")
      ("show_bar", value<bool>(), "show progress bar")
      ("promising_pose", value<int64_t>(), "The number of promising pose each fragment (default: 40)")
      ("cluster_size", value<fltype>(), "cluster size (Å) of promising poses (default: 1.0)")
      ("distance_width", value<fltype>(), "discretization width (Å) of relative distance (default: 0.1)")
      ("function", value<std::string>(), main_utils::option_desc("priority funcion", {std::make_pair("\"min\"", "only min"),
                                                                                      std::make_pair("\"simple\"", "no weight"),
                                                                                      std::make_pair("\"frac1\"", "1/rank weight"),
                                                                                      std::make_pair("\"frac2\"/default", "1/(rank)^2 weight")}
                                                                ).c_str());
    options_description desc;
    desc.add(options).add(hidden);
    variables_map vmap;
    store(command_line_parser(argc, argv).
	  options(desc).positional(pos_desc).run(), vmap);
    notify(vmap);

    if (!vmap.count("conf-file") || vmap.count("help")){
      if (!vmap.count("conf-file") && !vmap.count("help")){
	std::cout << "too few arguments" << std::endl;
      }
      std::cout << "Usage: ligandock conf-file [options]\n"
		<< options << std::endl;
      std::exit(0);
    }
    format::DockingConfiguration conf = format::ParseInFile(vmap["conf-file"].as<std::string>().c_str());
    if (vmap.count("receptor")) conf.receptor_file = vmap["receptor"].as<std::string>();
    if (vmap.count("fragments")) conf.fragments_file = vmap["fragments"].as<std::string>();
    if (vmap.count("output")) conf.output_file = vmap["output"].as<std::string>();
    if (vmap.count("grid")) conf.grid_folder = vmap["grid"].as<std::string>();
    if (vmap.count("log")) conf.log_file = vmap["log"].as<std::string>();
    if (vmap.count("oberror")) conf.ob_err_log = vmap["oberror"].as<std::string>() != "false";
    if (vmap.count("show_bar")) conf.show_bar = vmap["show_bar"].as<bool>();
    if (vmap.count("efficiency")) conf.efficiency = conf.parseFragmentEfficiency(vmap["efficiency"].as<std::string>());
    if (vmap.count("promising_pose")) conf.priority.promising_pose = vmap["promising_pose"].as<int64_t>();
    if (vmap.count("cluster_size")) conf.priority.cluster_size = vmap["cluster_size"].as<fltype>();
    if (vmap.count("distance_width")) conf.priority.distance_width = vmap["distance_width"].as<fltype>();
    if (vmap.count("function")) conf.priority.function = conf.parsePriorityFunctin(vmap["function"].as<std::string>());
    return conf;
  }

  void logConfig(const format::DockingConfiguration config){
    logs::lout << "Receptor file name  : "+config.receptor_file   << std::endl;
    logs::lout << "Fragments file name  : "+config.fragments_file << std::endl;
    logs::lout << "Output file name    : "+config.output_file     << std::endl;
    logs::lout << "Grid folder name    : "+config.grid_folder     << std::endl;
  }

  struct query_param {
    std::string smi1, smi2;
    fltype dist_min;
    fltype dist_max;
    fltype priority;
    query_param(const std::string& smi1, const std::string& smi2, fltype dist_min, fltype dist_max, fltype priority) 
      : smi1(smi1), smi2(smi2), dist_min(dist_min), dist_max(dist_max), priority(priority) {}
    query_param(const std::string& smi1, const std::string& smi2, std::pair<fltype, fltype> dist_range, fltype priority)
      : smi1(smi1), smi2(smi2), dist_min(dist_range.first), dist_max(dist_range.second), priority(priority) {}
    
    // 1st: priority, 2nd: id, 3rd: dist
    bool operator<(const query_param& o) {
      if (std::abs(priority - o.priority) < EPS) {
        if (smi1+smi2 == o.smi1+o.smi2) return dist_min < o.dist_min;
        else return smi1+smi2 < o.smi1+o.smi2;
      }
      return priority < o.priority;
    }
    friend std::ostream& operator<<(std::ostream& os, const query_param& qp);
  };

  std::ostream& operator<<(std::ostream& os, const query_param& qp) {
    os << qp.smi1 << " " << qp.smi2 << " " << (qp.dist_max + qp.dist_min) / 2 << " " << (qp.dist_max - qp.dist_min) / 2 << " " << qp.priority;
    return os;
  }

  void prioDebug(const format::Priority& priority) {
    logs::lout << logs::debug << "config.priority.promising_pose: " << priority.promising_pose      << std::endl;
    logs::lout << logs::debug << "config.priority.cluster_size  : " << priority.cluster_size        << std::endl;
    logs::lout << logs::debug << "config.priority.distance_width: " << priority.distance_width      << std::endl;
    logs::lout << logs::debug << "config.priority.function      : " << priority.getFunctionString() << std::endl;
  }
} // namespace

int main(int argc, char **argv){
  using namespace std;
  using namespace fragdock;
  using namespace main_utils;

  std::chrono::milliseconds whole_time(0);
  auto t0 = std::chrono::system_clock::now();

  format::DockingConfiguration config = parseArgs(argc, argv);

  if(config.log_file == ""){
    config.log_file = config.output_file + "__" + getDate() + ".log";
  }
  logs::log_init(config.log_file, true);
  logConfig(config);

  if (!config.ob_err_log) {
    OpenBabel::obErrorLog.StopLogging();
    logs::lout << "obErrorLog.StopLogging" << endl;
  }

  // parse receptor file
  const OpenBabel::OBMol receptor = format::ParseFileToOBMol(config.receptor_file.c_str())[0];
  const Molecule receptor_mol = format::toFragmentMol(receptor);

  // parse fragments file
  vector<OpenBabel::OBMol> fragments = format::ParseFileToOBMol(config.fragments_file);

  int frags_sz = fragments.size(); /* the number of fragments */
  logs::lout << "number of fragments: " << frags_sz << endl;

  // ================================================================
  // prepare atomgrids and rotations
  // ================================================================
  logs::lout << logs::info << "[start] read energy grids" << endl;
  vector<AtomInterEnergyGrid> atom_grids = AtomInterEnergyGrid::readAtomGrids(config.grid_folder);
  logs::lout << logs::info << "[ end ] read energy grids" << endl;
  // logs::lout << "atom grid size: " << atom_grids.size() << endl;

  const Point3d<int>& score_num = atom_grids[0].getNum(); // # of grid points of score grids (atom grids, fragment grids)


  logs::lout << logs::info << "[TIME STAMP] START MOLECULE OBJECT CONVERSION" << endl;
   vector<Fragment> fragments_frag = convert_fragments(fragments);
  logs::lout << logs::info << "[TIME STAMP] END   MOLECULE OBJECT CONVERSION" << endl;


  InterEnergyGrid distance_grid = makeDistanceGrid(atom_grids[0].getCenter(), atom_grids[0].getPitch(), atom_grids[0].getNum(), receptor_mol);

  // ---------------------------------------

  logs::lout << logs::debug << "config.efficiency : " << config.getFragmentEfficiencyString() << endl;
  prioDebug(config.priority);

  logs::lout << logs::info << "[TIME STAMP] START CALCULATING BY FRAGGRID" << endl;

  std::chrono::milliseconds query_time(0);
  std::chrono::milliseconds pose_time(0);
  int query_cnt = 0;

  auto t1 = std::chrono::system_clock::now();

  // vector of (position, score) of each fragment
  vector<vector<FragPose> > best_poses(frags_sz, vector<FragPose>());

  progress_bar grid_prog(frags_sz, "Searching fragment best poses", config.show_bar);

  for (int f_ind = 0; f_ind < frags_sz; ++f_ind) {
    FragmentInterEnergyGrid frag_grid(fragments_frag[f_ind], makeRotations60(), atom_grids, distance_grid);
    best_poses[f_ind] = frag_grid.getBestPoses(config.priority, calcFragmentEfficiency(fragments_frag[f_ind].heavyNum(), config.efficiency));
    grid_prog.display(f_ind+1);
  }
  grid_prog.clear();

  logs::lout << logs::info << "[TIME STAMP] END   CALCULATING BY FRAGGRID" << endl;


  auto t2 = std::chrono::system_clock::now();

  pose_time += std::chrono::duration_cast< std::chrono::milliseconds >(t2 - t1);
  logs::lout << logs::info << "[TIME STAMP] pose_time  : " << pose_time.count() << endl;

  std::chrono::milliseconds priority_time(0);


  logs::lout << logs::info << "[TIME STAMP] START CALCULATING DISTANCE PRIORITIES" << endl;

  progress_bar output_prog(frags_sz*frags_sz, "Calculating distance priorities", config.show_bar);

  vector<query_param> queries;
  for (int i = 0; i < frags_sz; ++i) {
    for (int j = i; j < frags_sz; ++j) {
      const vector<FragPose>& poses1 = best_poses[i];
      const vector<FragPose>& poses2 = best_poses[j];
      query::QueryPriority query_priority(config.priority);
      for (int k = 0; k < poses1.size(); ++k) {
        for (int l = (i==j)?k+1:0; l < poses2.size(); ++l) {
          // distance between two fragments
          fltype distance = (poses1[k].first - poses2[l].first).abs();
          query_priority.append(distance, poses1[k].second + poses2[l].second);
        }
      }
      vector<query::output_query> priorities = query_priority.getPriorities();
      for (int d = 0; d < priorities.size(); ++d) {
        queries.push_back(query_param(fragments_frag[i].getSmiles(), //fragments_frag[i].getFragmentID(),
                                      fragments_frag[j].getSmiles(), //fragments_frag[j].getFragmentID(), 
                                      priorities[d].dist_min,
                                      priorities[d].dist_max,
                                      priorities[d].priority));
        query_cnt++;
      }
      output_prog.display(i*frags_sz+j+1);
    }
  }
  output_prog.clear();

  std::sort(queries.begin(), queries.end());
  ofstream outputsh(config.output_file);
  int query_size = (NUM_QUERIES == -1) ? queries.size() : min(NUM_QUERIES, (int)queries.size());
  for (int i = 0; i < query_size; ++i) {
    outputsh << queries[i] << endl;
  }
  outputsh.close();

  auto t3 = std::chrono::system_clock::now();

  priority_time += std::chrono::duration_cast< std::chrono::milliseconds >(t3 - t2);
  query_time += std::chrono::duration_cast< std::chrono::milliseconds >(t3 - t1);

  logs::lout << logs::info << "[TIME STAMP] END   CALCULATING DISTANCE PRIORITIES" << endl;
  logs::lout << logs::info << "[TIME STAMP] priority_time  : " << priority_time.count() << endl;
  logs::lout << logs::info << "[TIME STAMP] query_time : " << query_time.count() << endl;
  logs::lout << logs::debug << "query_cnt  : " << query_cnt << endl;

  auto t4 = std::chrono::system_clock::now();
  whole_time += std::chrono::duration_cast< std::chrono::milliseconds >(t4 - t0);

  logs::lout << logs::info << "################ Program end ################" << endl;
  logs::lout << logs::info << "[TIME STAMP] whole_time : " << whole_time.count() << endl;

  logs::close();
  return 0;
}