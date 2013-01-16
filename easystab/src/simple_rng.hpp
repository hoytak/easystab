#ifndef _SIMPLE_RNG_H_
#define _SIMPLE_RNG_H_

#include <stdint.h>
#include <vector>

using namespace std;

// Simple lcg rng; super fast and good enough for us here
class CheapRNG { 
public:
  CheapRNG(size_t seed)
    : state( ( (uint64_t) seed) + 0xcb63b83e3798bbfeull)
  {}

  size_t operator()() {
    state = 6364136223846793005ull * state + 1442695040888963407ull;
    return (uint32_t)( (state >> 32) ^ state);
  }

  size_t operator()(size_t top) {
    return (*this)() % top;
  }

  double inRange(double r1, double r2) {
    double r_01 = ((double)((*this)() )) / size_t(-1);
    return r1 + (r2 - r1) * r_01;
  } 

private:
  uint64_t state;
};



class Shuffler { 
public:
  Shuffler(unsigned long seed, size_t max_idx) 
    : values(max_idx)
  {
    for(size_t i = 0; i < values.size(); ++i) 
      values[i] = i;

    CheapRNG rng(seed);

    random_shuffle(values.begin(), values.end(), rng);
  }

  size_t operator[](size_t idx) const {
    return values[idx];
  }

private:
  vector<size_t> values;
};


#endif /* _SIMPLE_RNG_H_ */
