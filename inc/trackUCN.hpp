#ifndef TRACKUCN_H
#define TRACKUCN_H

#include <vector>

typedef struct result {
    double energy;
    double theta;
    double t;
    double ePerp;
    double x;
    double y;
    double z;
    double zOff;
    int nHit;
    int nHitHouseLow;
    int nHitHouseHigh;
    double eStart;
    double deathTime;
} result;

result fixedEffDaggerHitTime(std::vector<double> state, double dt);

#endif /* TRACKUCN_H */