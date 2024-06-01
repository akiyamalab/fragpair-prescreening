#ifndef QPRIO_H_
#define QPRIO_H_

#include "common.hpp"
#include "infile_reader.hpp"
#include "utils.hpp"

#include <vector>

namespace query {
  enum EventType { START, END };

  struct dist_event {
    fltype dist, score;
    EventType type;
    dist_event(fltype dist, fltype score, EventType type) : dist(dist), score(score), type(type) {}
    bool operator<(const dist_event& o) const { return dist < o.dist; }
  };

  struct output_query {
    fltype dist_min, dist_max, priority;
    output_query(fltype dist_min, fltype dist_max, fltype priority) : dist_min(dist_min), dist_max(dist_max), priority(priority) {}
    bool is_next(const output_query& o) const { return abs(dist_max - o.dist_min) < EPS && abs(priority - o.priority) < EPS; }
  };

	class QueryPriority {
    const fltype width, radius;
		std::vector<dist_event> eves;
    format::Priority::Function function;
		int distToIndex(fltype dist) const;
    int indexToDist(int index) const { return index*width; }
    fltype indexToMinDist(int index) const { return index*width - width/2; }
    fltype indexToMaxDist(int index) const { return index*width + width/2; }
    fltype getPriority(const std::vector<fltype>& scores) const;
	public:
		QueryPriority(const format::Priority& priority) : width(priority.distance_width), radius(priority.cluster_size), function(priority.function) {}
		void append(fltype dist, fltype score);
		std::vector<output_query> getPriorities();
	};
}
#endif