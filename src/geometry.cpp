#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>

extern "C" {
    #include "../inc/xorshift.h"
}

std::vector<double> cross(std::vector<double> a, std::vector<double> b) {
    std::vector<double> res(3);
    res[0] = a[1]*b[2]-a[2]*b[1];
    res[1] = a[0]*b[2]-a[2]*b[0];
    res[2] = a[0]*b[1]-a[1]*b[0];
    return res;
}

double zOffDipCalc(double t) {
    double z;
    int i;
    
    double holdT = 100.0;
    double speed;
    
//    int nDips = 10;
//    double dipHeights[10] = {0.49, 0.380, 0.250, 0.180, 0.140, 0.110, 0.080, 0.060, 0.040, 0.010}; //9 dip
//    double dipEnds[10] = {0.0,  40.0,  80.0,  100.0, 120.0, 140.0, 160.0, 180.0, 200.0, 300.0}; //9 dip
    
    int nDips = 4;
    double dipHeights[4] = {0.49, 0.380, 0.250, 0.010}; //3 dip
    double dipEnds[4] = {holdT, holdT+20.0, holdT+40.0, holdT+140.0}; //3 dip
    
//    int nDips = 12;
//    double dipHeights[12] = {0.49, 0.250, 0.49, 0.380, 0.250, 0.180, 0.140, 0.110, 0.080, 0.060, 0.040, 0.010}; //9Dip PSE
//    double dipEnds[12] = {0.0,  200.0,  200.0+holdT, 200.0+holdT+20.0, 200.0+holdT+40.0, 200.0+holdT+50.0,
//                200.0+holdT+60.0, 200.0+holdT+70.0, 200.0+holdT+80.0, 200.0+holdT+90.0,
//                200.0+holdT+100.0, 200.0+holdT+120.0}; //9 dip PSE
    
    if(t > dipEnds[nDips-1]) {
        return 0.01;
    }
    
    for(i = 0; i < nDips; i++) {
        if(dipEnds[i] > t) {
            break;
        }
    }
    
    if(i == 0) {
        return dipHeights[0];
    }
    
    if(dipHeights[i-1] > dipHeights[i]) {
        speed = -0.49/13.0;
        z = dipHeights[i-1] + speed*(t-dipEnds[i-1]);
        return z > dipHeights[i] ? z : dipHeights[i];
    }
    else {
        speed = 0.49/13.0;
        z = dipHeights[i-1] + speed*(t-dipEnds[i-1]);
        return z < dipHeights[i] ? z : dipHeights[i];
    }
}

void reflect(std::vector<double> &state, std::vector<double> norm, std::vector<double> tang) {
    double pTarget = sqrt(state[3]*state[3] + state[4]*state[4] + state[5]*state[5]);
    std::vector<double> tangPrime = cross(norm, tang);
    
    double u1 = nextU01();
    double u2 = nextU01();
    
    double theta = asin(sqrt(u1));
    double phi = 2 * M_PI * u2;
    
    double pN = cos(theta);
    double pT = sin(theta)*cos(phi);
    double pTprime = sin(theta)*sin(phi);
    std::vector<double> newPdir(3);
    newPdir[0] = pN*norm[0] + pT*tang[0] + pTprime*tangPrime[0];
    newPdir[1] = pN*norm[1] + pT*tang[1] + pTprime*tangPrime[1];
    newPdir[2] = pN*norm[2] + pT*tang[2] + pTprime*tangPrime[2];
    
    double pLen = sqrt(newPdir[0]*newPdir[0] + newPdir[1]*newPdir[1] + newPdir[2]*newPdir[2]);
    
    state[3] = newPdir[0] * pTarget / pLen;
    state[4] = newPdir[1] * pTarget / pLen;
    state[5] = newPdir[2] * pTarget / pLen;
}

bool checkDagHit(double x, double y, double z, double zOff) {
    double zeta;
    if(x > 0) {
        zeta = 0.5 - sqrt(x*x + pow(fabs(z - zOff) - 1.0, 2));
    }
    else {
        zeta = 1.0 - sqrt(x*x + pow(fabs(z - zOff) - 0.5, 2));
    }
    if((x > -0.3524) && (x < 0.0476) && (zeta > 0.0) && (z < (-1.5 + zOff + 0.2))) {
        return true;
    }
    return false;
}

bool checkHouseHitLow(double x, double y, double z, double zOff) {
    if(z >= (-1.5 + zOff + 0.2) && z < (-1.5 + zOff + 0.2 + 0.14478) && fabs(x + 0.1524) < (0.40 + 2.0179*(z + 1.5 - zOff - 0.2))/2.0) {
        return true;
    }
    return false;
}

bool checkHouseHitHigh(double x, double y, double z, double zOff) {
    if(z >= (-1.5 + zOff + 0.2 + 0.14478) && z < (-1.5 + zOff + 0.2 + 0.2667) && fabs(x + 0.1524) < 0.69215/2.0) {
        return true;
    }
    return false;
}
