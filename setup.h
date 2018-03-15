#define INITBLOCK 0

#define TRACKGENERATOR randomPointTrapEdE
//#define TRACKGENERATOR randomPointTrapOptimum

#define TRACKER daggerHitTimes
//#define TRACKER fixedEffDaggerHitTime

#define WRITER writeNoabsRes
//#define WRITER writeFixedRes

#define SETTLINGTIME 150

#define FIRSTDIPTIME 20

#define HOLDTIME 20

//9 Dip
#define NDIPS 10
#define HEIGHTS {0.49, 0.380, 0.250, 0.180, 0.140, 0.110, 0.080, 0.060, 0.040, 0.010}
#define ENDTIMES {holdT,  holdT+40.0,  holdT+80.0,  holdT+100.0, holdT+120.0, holdT+140.0, holdT+160.0, holdT+180.0, holdT+200.0, holdT+300.0}

//#define NDIPS 4
//#define HEIGHTS {0.49, 0.380, 0.250, 0.010}
//#define ENDTIMES {holdT, holdT+20.0, holdT+40.0, holdT+140.0}

//#define NDIPS 12
//#define HEIGHTS {0.49, 0.250, 0.49, 0.380, 0.250, 0.180, 0.140, 0.110, 0.080, 0.060, 0.040, 0.010}
//#define ENDTIMES {0.0,  200.0,  200.0+holdT, 200.0+holdT+20.0, 200.0+holdT+40.0, 200.0+holdT+50.0, 200.0+holdT+60.0, 200.0+holdT+70.0, 200.0+holdT+80.0, 200.0+holdT+90.0, 200.0+holdT+100.0, 200.0+holdT+120.0}