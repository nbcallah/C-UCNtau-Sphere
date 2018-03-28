#ifndef TRACKUCN_H
#define TRACKUCN_H

#include <vector>
#define NRECORDS 50

typedef struct noabsResult {
    double energy;
    double theta;
    float times[NRECORDS];
    float ePerps[NRECORDS];
} noabsResult;

typedef struct fixedResult {
    double energy;
    double theta;
    double t;
    double settlingT;
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
} fixedResult;

fixedResult fixedEffDaggerHitTime(std::vector<double> state, double dt);
noabsResult daggerHitTimes(std::vector<double> state, double dt);

#endif /* TRACKUCN_H */