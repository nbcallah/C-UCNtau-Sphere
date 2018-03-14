#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <vector>

std::vector<double> cross(std::vector<double> a, std::vector<double> b);

double zOffDipCalc(double t);

void reflect(std::vector<double> &state, std::vector<double> norm, std::vector<double> tang);

bool checkDagHit(double x, double y, double z, double zOff);

bool checkHouseHitLow(double x, double y, double z, double zOff);

bool checkHouseHitHigh(double x, double y, double z, double zOff);

#endif /* GEOMETRY_H */