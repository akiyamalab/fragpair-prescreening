#ifndef COMMON_H_
#define COMMON_H_
typedef float fltype;

// below codes are used by fragdock and exhdock
#include "log_writer_stream.hpp"
const fltype INF_ENERGY = 1e9;
const fltype LIMIT_ENERGY = 1e2;
const fltype EPS = 1e-4;

const int NUM_QUERIES = -1;
const std::string ID_KEY = "fragmentID";

#endif

