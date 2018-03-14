#include <vector>
#include <stdio.h>
#include "inc/track_gen.hpp"
#include "inc/trackUCN.hpp"
#include <cmath>
#include <mpi.h>
//#include "inc/geometry.hpp"

extern "C" {
    #include "inc/xorshift.h"
}

/*typedef struct result {
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
} result;*/

int main(int argc, char** argv) {
    int ierr = MPI_Init(&argc, &argv);
    int nproc;
    int rank;
    
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    
    initxorshift(0);
    
//    printf("%d - %d\n", rank*(4096000/nproc), (rank+1)*(4096000/nproc));
    
    std::vector<std::vector<double>> traj;
    std::vector<int> is;
    traj.reserve(4096000/nproc);
    for(int i = 0; i < 4096000; i++) {
        if(i >= rank*(4096000/nproc) && i < (rank+1)*(4096000/nproc)) {
            traj.push_back(randomPointTrapOptimum());
            is.push_back(i);
        }
        else {
            randomPointTrapOptimum();
        }
    }
    
    for(int i = 0; i < rank; i++) {
        jump();
    }
        
    for(auto it = traj.begin(); it < traj.end(); it++) {
        result res = fixedEffDaggerHitTime(*it, 0.0005);
    }
    
    ierr = MPI_Finalize();
    
    return 0;
}