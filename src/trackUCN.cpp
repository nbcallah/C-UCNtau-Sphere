#include "../inc/trackUCN.hpp"
#include "../inc/constants.h"
#include "../inc/symplectic.hpp"
#include "../inc/geometry.hpp"
#include "../inc/quant_refl.hpp"
#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>
#include <assert.h>

extern "C" {
    #include "../inc/fields_nate.h"
    #include "../inc/xorshift.h"
}

result fixedEffDaggerHitTime(std::vector<double> state, double dt) {
    std::vector<double> tang = {0.0, 0.0, 1.0};
    std::vector<double> normPlus = {0.0, 1.0, 0.0};
    std::vector<double> normMinus = {0.0, -1.0, 0.0};
    result res;
    res.theta = acos(state[5]/sqrt(state[3]*state[3] + state[4]*state[4] + state[5]*state[5]));
    double t = 0;
    double freq = 0;
    potential(&state[0], &state[1], &state[2], &(res.eStart), &t, &freq);
    res.eStart = res.eStart - MINU + (state[3]*state[3] + state[4]*state[4] + state[5]*state[5])/(2*MASS_N);
    
    double deathTime = -877.7*log(nextU01());
    
    double settlingTime = 50.0;
    int numSteps = settlingTime/dt;
    double energy;
    for(int i = 0; i < numSteps; i++) {
        symplecticStep(state, dt, energy);
        t = t + dt;
    }
    
    int nHit = 0;
    int nHitHouseLow = 0;
    int nHitHouseHigh = 0;
    std::vector<double> prevState(6);
    while(true) {
        prevState = state;
        symplecticStep(state, dt, energy);
        t = t + dt;
        if(t - settlingTime > deathTime) {
            res.energy = energy;
            res.t = t-settlingTime;
            res.ePerp = state[4]*state[4]/(2*MASS_N);
            res.x = state[0];
            res.y = state[1];
            res.z = state[2];
            res.zOff = -1;
            res.nHit = nHit;
            res.nHitHouseLow = nHitHouseLow;
            res.nHitHouseHigh = nHitHouseHigh;
            res.deathTime = deathTime;
            return res;
        }
        if((prevState[1] < 0 && state[1] > 0) || (prevState[1] > 0 && state[1] < 0)) {
            double fracTravel = fabs(prevState[1])/(fabs(state[1]) + fabs(prevState[1]));
            double predX = prevState[0] + fracTravel * (state[0] - prevState[1]);
            double predZ = prevState[2] + fracTravel * (state[2] - prevState[2]);
            
            double zOff = zOffDipCalc(t - settlingTime);
            
            if(checkDagHit(predX, 0.0, predZ, zOff)) {
                nHit += 1;
                if(absorbMultilayer(state[4]*state[4]/(2*MASS_N), 4.6)) {
                    res.energy = energy;
                    res.t = t-settlingTime;
                    res.ePerp = state[4]*state[4]/(2*MASS_N);
                    res.x = predX;
                    res.y = 0.0;
                    res.z = predZ;
                    res.zOff = zOff;
                    res.nHit = nHit;
                    res.nHitHouseLow = nHitHouseLow;
                    res.nHitHouseHigh = nHitHouseHigh;
                    res.deathTime = deathTime;
                    return res;
                }
                if(prevState[1] > 0 && prevState[4] < 0) {
                    reflect(prevState, normPlus, tang);
                }
                else {
                    reflect(prevState, normMinus, tang);
                }
                state = prevState;
            }
            else if(checkHouseHitLow(predX, 0.0, predZ, zOff)) {
                nHitHouseLow += 1;
                if(prevState[1] > 0 && prevState[4] < 0) {
                    reflect(prevState, normPlus, tang);
                }
                else {
                    reflect(prevState, normMinus, tang);
                }
                state = prevState;
            }
            else if(checkHouseHitHigh(predX, 0.0, predZ, zOff)) {
                nHitHouseHigh += 1;
                if(prevState[1] > 0 && prevState[4] < 0) {
                    reflect(prevState, normPlus, tang);
                }
                else {
                    reflect(prevState, normMinus, tang);
                }
                state = prevState;
            }
        }
    }
}