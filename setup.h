#define INITBLOCK 7

//#define TRACKGENERATOR randomPointTrapEdE
//#define TRACKGENERATOR randomPointTrapOptimum
#define TRACKGENERATOR randomPointTrapOptimumCleanable
//#define TRACKGENERATOR randomPointTrapEdECleanable

//#define TRACKER daggerHitTimes
#define TRACKER fixedEffDaggerHitTime

//#define WRITER writeNoabsRes
#define WRITER writeFixedRes

#define CLEANINGTIME 200

#define FIRSTDIPTIME 20

#define HOLDTIME 20

//9 Dip
//#define NDIPS 10
//#define HEIGHTS {0.49, 0.380, 0.250, 0.180, 0.140, 0.110, 0.080, 0.060, 0.040, 0.010}
//#define ENDTIMES {holdT,  holdT+40.0,  holdT+80.0,  holdT+100.0, holdT+120.0, holdT+140.0, holdT+160.0, holdT+180.0, holdT+200.0, holdT+300.0}

//#define NDIPS 4
//#define HEIGHTS {0.49, 0.380, 0.250, 0.010}
//#define ENDTIMES {holdT, holdT+20.0, holdT+40.0, holdT+140.0}

//#define NDIPS 12
//#define HEIGHTS {0.49, 0.250, 0.49, 0.380, 0.250, 0.180, 0.140, 0.110, 0.080, 0.060, 0.040, 0.010}
//#define ENDTIMES {0.0,  200.0,  200.0+holdT, 200.0+holdT+20.0, 200.0+holdT+40.0, 200.0+holdT+50.0, 200.0+holdT+60.0, 200.0+holdT+70.0, 200.0+holdT+80.0, 200.0+holdT+90.0, 200.0+holdT+100.0, 200.0+holdT+120.0}

//dipHeights = (/0.49_8, 0.380_8, 0.250_8, 0.01_8/)
//dipEnds =     (/0.0_8,  40.0_8,  400.0_8, 500.0_8/)
#define NDIPS 10
#define HEIGHTS {0.49, 0.38, 0.25, 0.18, 0.14, 0.11, 0.08, 0.06, 0.04, 0.01}
#define ENDTIMES {holdT, holdT+40.0, holdT+80, holdT+100, holdT+120, holdT+140, holdT+160, holdT+180, holdT+200.0, holdT+500.0}

#define XFNAME "/N/u/nbcallah/BigRed2/ChaoticTrap/C-UCNtau-Trap-Sims/xvals.bin"
#define YFNAME "/N/u/nbcallah/BigRed2/ChaoticTrap/C-UCNtau-Trap-Sims/yvals.bin"
#define ZFNAME "/N/u/nbcallah/BigRed2/ChaoticTrap/C-UCNtau-Trap-Sims/zvals.bin"
//#define XFNAME "./xvals.bin"
//#define YFNAME "./yvals.bin"
//#define ZFNAME "./zvals.bin"
//#define AMPLITUDE 0.00002
//#define FREQ 60.0
