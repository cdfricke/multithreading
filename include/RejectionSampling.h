#ifndef REJECTIONSAMPLING_H
#define REJECTIONSAMPLING_H

    #include <iostream>
    #include <random>
    #include <cmath>
    #include <vector>
    #include "Calculus.h"
    #include "Stats.h"
    #ifdef WRITE_DATA
        #include <fstream>
        #include <string>
    #endif

    /* SHAPES */
    #define SPHERE 0
    #define CUBE 1
    #define CUBOID 2

    /* DISTRIBUTIONS */
    #define CE3 0
    #define CE4 1

    /* MACROS */
    #define RAND distrib(gen)   // RANDOM double between 0.0 and 1.0 (generator must be seeded beforehand)

    /* CUSTOM TYPES */
    struct Coord {  // used for position of rocks and dimensions of cuboids
        double x;
        double y;
        double z;
    };
    struct ProgramParams {  // stores all params of the main program.
        // CONTROL DEFAULT PARAMS HERE 
        int SHAPE = SPHERE;
        int DIST = CE3;
        double DMIN = 0.05;
        double DMAX = 3.0;
        double DSTEP = 0.001;
        double DEPTH = 10.0;
        double AREA = 50.0*50.0;
    };

    /* TYPEDEFS */
    typedef std::vector<double> Vector;
    typedef std::vector<Coord> CoordVector;

    /* HELPER FUNCTION PROTOTYPES */
    ProgramParams setParams(int argc, char** argv);
    void setGeoConst(int SHAPE);
    void setDist(int DIST);
    double Fk(double D);
    double nDpm2(double D);
    double nDpm3(double D);
    double DnDpm3(double D);
    double expintA(double x);
    double NgtD4(double D);
    double hpD(double x);
    double bpD(double x);
    double F_SVII(double Dmin, double Dmax);
    double n_a_SVII(double D);
    double n_V_SVII(double D);
    Vector N_SVII(double Dmin, double Dmax);
    Vector NgtD4_Di(double D, const double INF, const int N);
    Vector CheckNgtD(Vector& rocks, double Dmin, double Dmax, double Dstep, double A);
    Vector CheckNgtD(CoordVector& rocks, double Dmin, double Dmax, double Dstep, double A);
    Vector CheckFgtD(Vector& rocks, double Dmin, double Dmax, double Dstep, double A, const int SHAPE);
    Vector CheckFgtD(CoordVector& rocks, double Dmin, double Dmax, double Dstep, double A);

#endif
