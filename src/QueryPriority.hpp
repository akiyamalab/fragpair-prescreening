#ifndef QPRIO_H_
#define QPRIO_H_

#include "common.hpp"
#include "infile_reader.hpp"
#include "utils.hpp"

#include <vector>

namespace query {
  enum EventType { START, END };

  /**
   * Struct for storing distance and score of an event (start or end).
   * It is used for "event sort".
   * @param dist distance of the event
   * @param score score of the event
   * @param type type of the event (START or END)
  */
  struct dist_event {
    fltype dist, score;
    EventType type;
    
    dist_event(fltype dist, fltype score, EventType type) : dist(dist), score(score), type(type) {}
    bool operator<(const dist_event& o) const { return dist < o.dist; }
  };

  /**
   * Struct for storing distance range and priority of an output query.
   * @param dist_min minimum distance of the query
   * @param dist_max maximum distance of the query
   * @param priority priority of the query
  */
  struct output_query {
    fltype dist_min, dist_max, priority;

    output_query(fltype dist_min, fltype dist_max, fltype priority) : dist_min(dist_min), dist_max(dist_max), priority(priority) {}

    // check if the next query is the same as this query
    bool is_next(const output_query& o) const { return abs(dist_max - o.dist_min) < EPS && abs(priority - o.priority) < EPS; }
  };

  /**
   * Class for storing and calculating query priority.
   * It is used for calculating query priority.
   * 
  */
	class QueryPriority {
    const fltype width, radius;
		std::vector<dist_event> eves;
    format::Priority::Function function;

		int distToIndex(fltype dist) const;
    int indexToDist(int index) const { return index*width; }
    fltype indexToMinDist(int index) const { return index*width - width/2; }
    fltype indexToMaxDist(int index) const { return index*width + width/2; }

    /**
     * Get the priority of a distance range using the function.
     * @param[in] scores a vector of scores (not necessarily sorted when input)
     * @return the priority using function
    */
    fltype getPriority(const std::vector<fltype>& scores) const;
	public:
		QueryPriority(const format::Priority& priority) : width(priority.distance_width), radius(priority.cluster_size), function(priority.function) {}

    // append a distance and a score to the priority
		void append(fltype dist, fltype score);

    /**
     * Get the priority of each distance range.
     * @return a vector of output_query
    */
		std::vector<output_query> getPriorities();
	};
}
#endif