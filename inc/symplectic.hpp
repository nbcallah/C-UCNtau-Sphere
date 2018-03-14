#ifndef SYMPLECTIC_H
#define SYMPLECTIC_H

#include <vector>

void symplecticStep(std::vector<double> &state, double deltaT, double &energy);

#endif /* SYMPLECTIC_H */