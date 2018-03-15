#include <vector>
#include <stdio.h>
#include "inc/track_gen.hpp"
#include "inc/trackUCN.hpp"
#include <cmath>
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <getopt.h>
#include <cstring>

#include "setup.h"

//#include "inc/geometry.hpp"

extern "C" {
    #include "inc/xorshift.h"
}

/*typedef struct fixedResult {
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
} fixedResult;*/

void writeFixedRes(std::ofstream &binfile, fixedResult res) {
    const size_t buff_len = sizeof(unsigned int) + 8*sizeof(double) + 3*sizeof(int) + 2*sizeof(double) + sizeof(unsigned int);
    char buf[buff_len];
    if(!binfile.is_open()) {
        fprintf(stderr, "Error! file closed\n");
        return;
    }
    *((unsigned int *)(&buf[0])) = buff_len - 2*sizeof(unsigned int);
    *((double *)(&buf[0] + sizeof(unsigned int))) = res.energy;
    *((double *)(&buf[0] + sizeof(unsigned int) + 1*sizeof(double))) = res.theta;
    *((double *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(double))) = res.t;
    *((double *)(&buf[0] + sizeof(unsigned int) + 3*sizeof(double))) = res.ePerp;
    *((double *)(&buf[0] + sizeof(unsigned int) + 4*sizeof(double))) = res.x;
    *((double *)(&buf[0] + sizeof(unsigned int) + 5*sizeof(double))) = res.y;
    *((double *)(&buf[0] + sizeof(unsigned int) + 6*sizeof(double))) = res.z;
    *((double *)(&buf[0] + sizeof(unsigned int) + 7*sizeof(double))) = res.zOff;
    *((int *)(&buf[0] + sizeof(unsigned int) + 8*sizeof(double))) = res.nHit;
    *((int *)(&buf[0] + sizeof(unsigned int) + 8*sizeof(double) + 1*sizeof(int))) = res.nHitHouseLow;
    *((int *)(&buf[0] + sizeof(unsigned int) + 8*sizeof(double) + 2*sizeof(int))) = res.nHitHouseHigh;
    *((double *)(&buf[0] + sizeof(unsigned int) + 8*sizeof(double) + 3*sizeof(int))) = res.eStart;
    *((double *)(&buf[0] + sizeof(unsigned int) + 8*sizeof(double) + 3*sizeof(int) + 1*sizeof(double))) = res.deathTime;
    *((unsigned int *)(&buf[0] + sizeof(unsigned int) + 8*sizeof(double) + 3*sizeof(int) + 2*sizeof(double))) = buff_len - 2*sizeof(unsigned int);
    binfile.write(buf, buff_len);
}

void writeNoabsRes(std::ofstream &binfile, noabsResult res) {
    const size_t buff_len = sizeof(unsigned int) + 2*sizeof(double) + 2*NRECORDS*sizeof(float) + sizeof(unsigned int);
    char buf[buff_len];
    if(!binfile.is_open()) {
        fprintf(stderr, "Error! file closed\n");
        return;
    }
    *((unsigned int *)(&buf[0])) = buff_len - 2*sizeof(unsigned int);
    *((double *)(&buf[0] + sizeof(unsigned int))) = res.energy;
    *((double *)(&buf[0] + sizeof(unsigned int) + 1*sizeof(double))) = res.theta;
    std::memcpy((void *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(double)), (void *)&(res.times[0]), NRECORDS*sizeof(float));
    std::memcpy((void *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(double) + NRECORDS*sizeof(float)), (void *)&(res.ePerps[0]), NRECORDS*sizeof(float));
    *((unsigned int *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(double) + 2*NRECORDS*sizeof(float))) = buff_len - 2*sizeof(unsigned int);
    binfile.write(buf, buff_len);
}

int main(int argc, char** argv) {
    int ierr = MPI_Init(&argc, &argv);
    int nproc;
    int rank;
    
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    
    
    
    int c;
    
    char fName[1024];
    fName[0] = '\0';
    double dt = 0.0;
    double nTraj = 0;

    while (1) {
        static struct option long_options[] = {
            {"file", required_argument, 0, 'f'},
            {"dt", required_argument, 0, 'd'},
            {"ntraj", required_argument, 0, 'n'},
            {0, 0, 0, 0},
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long(argc, argv, "f:d:n:", long_options, &option_index);

        /* Detect the end of the options. */
        if(c == -1) {
            break;
        }

        switch(c) {
            case 'f':
                strncpy(fName, optarg, 1024-1);
                fName[1024-1] = '\0';
                break;
            case 'd':
                dt = atof(optarg);
                break;
            case 'n':
                nTraj = atoi(optarg);
                break;
            case '?':
                /* getopt_long already printed an error message. */
                break;
            default:
                exit(1);
        }
    }
    
    if(fName[0] == '\0' || dt == 0 || nTraj == 0) {
        fprintf(stderr, "Error! Usage: ./arrival_time_fixed_eff --file=fName --dt=timestep --ntraj=N\n");
        exit(1);
    }
    
    char fNameRank[1024];
    snprintf(fNameRank, 1024, "%s%d", fName, rank);
    std::ofstream binfile(fNameRank, std::ios::out | std::ios::binary);
    
    
    initxorshift(INITBLOCK);
    
//    printf("%d - %d\n", rank*(nTraj/nproc), (rank+1)*(nTraj/nproc));
    
    std::vector<std::vector<double>> traj;
    traj.reserve(nTraj/nproc);
    for(int i = 0; i < nTraj; i++) {
        if(i >= rank*(nTraj/nproc) && i < (rank+1)*(nTraj/nproc)) {
            traj.push_back(TRACKGENERATOR());
        }
        else {
            TRACKGENERATOR();
        }
    }
    
    for(int i = 0; i < rank; i++) {
        jump();
    }
        
    for(auto it = traj.begin(); it < traj.end(); it++) {
        auto res = TRACKER(*it, dt);
//        writeFixedRes(binfile, res);
        WRITER(binfile, res);
    }
    
    binfile.close();
    ierr = MPI_Finalize();
    
    return 0;
}
