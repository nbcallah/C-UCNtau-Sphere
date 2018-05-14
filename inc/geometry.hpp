#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <vector>

std::vector<double> cross(std::vector<double> a, std::vector<double> b);

void normalize(std::vector<double> &a);

double zOffDipCalc(double t);

void reflect(std::vector<double> &state, std::vector<double> norm, std::vector<double> tang);

bool checkDagHit(double x, double y, double z, double zOff);

int checkClean(std::vector<double> state, std::vector<double> prevState, double cleanHeight);

double calcDagZeta(double x, double y, double z, double zOff);

bool checkHouseHitLow(double x, double y, double z, double zOff);

bool checkHouseHitHigh(double x, double y, double z, double zOff);

std::vector<double> initializeLyapState(std::vector<double> ref);

void resetStates(std::vector<double> ref, std::vector<double> &pair);

double distance(std::vector<double> ref, std::vector<double> pair);


#endif /* GEOMETRY_H */