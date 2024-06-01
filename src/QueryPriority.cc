#include "QueryPriority.hpp"

#include <set>

namespace query {
  int QueryPriority::distToIndex(fltype dist) const {
    int q = dist/width;
    fltype r = dist - q*width;
    if (r >= width/2) q++;
		return (q < 0) ? 0 : q;
	}

  fltype QueryPriority::getPriority(const std::vector<fltype>& scores) const {
    std::vector<fltype> sorted_scores = scores;
		sort(sorted_scores.begin(), sorted_scores.end());
		fltype prio = 0;
    switch (function) {
      case format::Priority::Function::MIN:
        if (scores.size() > 0) prio = sorted_scores[0];
        break;
      case format::Priority::Function::SIMPLE:
        for (int i = 0; i < sorted_scores.size(); i++) {
          prio += sorted_scores[i];
        }
        break;
      case format::Priority::Function::FRAC1:
        for (int i = 0; i < sorted_scores.size(); i++) {
          fltype si = sorted_scores[i] / (i+1);
          if (abs(si) < EPS/10) break;
          prio += si;
        }
        break;
      case format::Priority::Function::FRAC2:
      default:
        for (int i = 0; i < sorted_scores.size(); i++) {
          fltype si = sorted_scores[i] / ((i+1)*(i+1));
          if (abs(si) < EPS/10) break;
          prio += si;
        }
        break;
    }
		return prio;
	}

  void QueryPriority::append(fltype dist, fltype score) {
    eves.push_back(dist_event(dist - radius*2, score, EventType::START));
    eves.push_back(dist_event(dist + radius*2 + width, score, EventType::END)); // dist-radius*2 <= range <"=" dist+radius*2+width
  }

  std::vector<output_query> QueryPriority::getPriorities() {
    std::multiset<fltype> scores;
    std::vector<output_query> ret;

    std::sort(eves.begin(), eves.end());

    bool event = false;
    int e_ind = 0, d_ind = 0;
    while (e_ind < eves.size()) {
      while (distToIndex(eves[e_ind].dist) <= d_ind && e_ind < eves.size()) {
        if (eves[e_ind].type == EventType::START) {
          scores.insert(eves[e_ind].score);
        } else {
          scores.erase(scores.find(eves[e_ind].score));
        }
        e_ind++;
        event = true;
      }
      
      if (!scores.empty()) {
        fltype priority = getPriority(std::vector<fltype>(scores.begin(), scores.end()));
        if (!ret.empty() && !event) {
          ret.back().dist_max = indexToMaxDist(d_ind);
        } else {
          ret.push_back(output_query(indexToMinDist(d_ind), indexToMaxDist(d_ind), priority));
        }
      }
      d_ind++;
      event = false;
    }
    
    for (int i = 1; i < ret.size(); i++) {
      if (ret[i-1].is_next(ret[i])) {
        ret[i-1].dist_max = ret[i].dist_max;
        ret.erase(ret.begin() + i);
        i--;
      }
    }

    return ret;
  }
}