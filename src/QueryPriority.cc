#include "QueryPriority.hpp"

#include <set>

namespace query {
  int QueryPriority::distToIndex(fltype dist) const {
    int q = dist/width;
    // round
    fltype r = dist - q*width;
    if (r >= width/2) q++;

    // return (>0)
		return (q < 0) ? 0 : q;
	}

  fltype QueryPriority::getPriority(const std::vector<fltype>& scores) const {
    // sort scores
    std::vector<fltype> sorted_scores = scores;
		sort(sorted_scores.begin(), sorted_scores.end());

		fltype prio = 0;
    // calculate priority using function
    switch (function) {
      case format::Priority::Function::MIN:
        // min(scores)
        if (scores.size() > 0) prio = sorted_scores[0];
        break;
      case format::Priority::Function::SIMPLE:
        // sum score_i
        for (int i = 0; i < sorted_scores.size(); i++) {
          prio += sorted_scores[i];
        }
        break;
      case format::Priority::Function::FRAC1:
        // sum score_i/(i+1)
        for (int i = 0; i < sorted_scores.size(); i++) {
          fltype si = sorted_scores[i] / (i+1);
          if (abs(si) < EPS/10) break;
          prio += si;
        }
        break;
      case format::Priority::Function::FRAC2:
      default:
        // sum score_i/((i+1)^2)
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
    // dist-radius*2 <= range <"=" dist+radius*2+width
    eves.push_back(dist_event(dist - radius*2, score, EventType::START));
    eves.push_back(dist_event(dist + radius*2 + width, score, EventType::END));
  }

  std::vector<output_query> QueryPriority::getPriorities() {
   // score multiset (set including duplicates) updated in real time
    std::multiset<fltype> scores;
    // return value
    std::vector<output_query> ret;

    // sort events by distance
    std::sort(eves.begin(), eves.end());

    bool event = false;
    int e_ind = 0, d_ind = 0;
    // iterate over events
    while (e_ind < eves.size()) {
      // update score multiset by events
      while (distToIndex(eves[e_ind].dist) <= d_ind && e_ind < eves.size()) {
        if (eves[e_ind].type == EventType::START) {
          scores.insert(eves[e_ind].score);
        } else {
          scores.erase(scores.find(eves[e_ind].score));
        }
        e_ind++;
        event = true;
      }
      
      // calculate output_query
      if (!scores.empty()) {
        // calculate priority (score: multiset -> vector)
        fltype priority = getPriority(std::vector<fltype>(scores.begin(), scores.end()));
        // update ret
        if (!ret.empty() && !event) { // if the scores is the same as the last query
          // merge adjacent queries
          ret.back().dist_max = indexToMaxDist(d_ind);
        } else {
          // append a new query
          ret.push_back(output_query(indexToMinDist(d_ind), indexToMaxDist(d_ind), priority));
        }
      }
      d_ind++;
      event = false;
    }
    
    // merge adjacent queries (if the priority is the same)
    // rounding may equal priority even if scores are different
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