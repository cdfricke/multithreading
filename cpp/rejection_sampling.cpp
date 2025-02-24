/*
    // Original (Rejection_Sampling_V6.py) by Payton Linton
    Purpose: This code is meant to simulate populating
            a volume according to the CE4 distribution.
            The CE4 distribution is an area density and I
            want to populate a volume SUCH THAT if I take any
            slice in the z-direction, on average I can get back
            the fractional area distribution reported by CE4
            in Wu et al. (2021)

    NOTE: Cuboids are slightly different. Spheres, and Cubes
         are completely defined by D, but cuboids need a
         diameter, D, a width, b, and a height, h

    NOTE: I am not assuming anything fancy like giving the shapes
         some random orientation. So for Cubes and Cuboids, the assumption
         is that their faces point along the x, y, z directions. Therefore
         when determining if a particular rock intersects a plane at point z0
         the only relevant variable is its height, and z position.
*/

// Programmer: Connor Fricke (fricke.59@osu.edu)
// file: rejection_sampling.cpp
// Date: 19-FEB-2025
// Desc: C++ translation and multithreading version of Rejection_Sampling_V6.py by Payton Linton (linton.93@osu.edu)

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#ifdef WRITE_DATA
    #include <fstream>
    #include <string>
#endif
#include "Calculus.h"
#include "Stats.h"

/* DEFINE SHAPES */
#define SPHERE 0
#define CUBE 1
#define CUBOID_A 2
#define CUBOID_B 3

/* DEFINE DISTRIBUTIONS */
#define CE3 0
#define CE4 1

/* MACROS */
#define RAND (double(std::rand())/RAND_MAX)   // RANDOM double between 0.0 and 1.0

/* CUSTOM TYPES */
struct Coord
{
    double x;
    double y;
    double z;
};
typedef std::vector<double> Vector;
typedef std::vector<Coord> CoordVector;

/* DEFINE GLOBAL VARS */
static const double PI = 4.0 * std::atan(1);    // PI
static double GeoConst = 0.0;                   // INIT WITH setGeoConst();
static double k = 0.0;                          // INIT WITH setDist();
static double qk = 0.0;                         // INIT WITH setDist();

/* HELPER FUNCTION PROTOTYPES */
void setGeoConst(int SHAPE);
void setDist(int DIST);
double Fk(double D);
//double Fk_Di(double D);
double nDpm2(double D);
double nDpm3(double D);
double DnDpm3(double D);
double expintA(double x);
double NgtD4(double D, const int SHAPE);
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

/* MAIN FUNCTION PROTOTYPES */
int theBigOne_NonCuboid(const int SHAPE, const int DIST, const double DMIN, const double DMAX, const double DSTEP, const double DEPTH, const double AREA);
int theBigOne_Cuboid(const int SHAPE, const int DIST, const double DMIN, const double DMAX, const double DSTEP, const double DEPTH, const double AREA);

/* MAIN PROGRAM */
int main(int argc, char **argv)
{
    int SHAPE, DIST;
    double DMIN, DMAX, DEPTH, AREA;
    if (argc == 1) { // defaults
        std::cout << "Using DEFAULT Parameters..." << std::endl;
        SHAPE = SPHERE;
        DIST = CE3;
        DMIN = 0.05;
        DMAX = 3.0;
        DEPTH = 10.0;
        AREA = 50.0*50.0;
    }
    else if (argc == 7) { // user specified
        std::cout << "Using GIVEN Parameters..." << std::endl;
        SHAPE = std::atoi(argv[1]);
        DIST = std::atoi(argv[2]);
        DMIN = std::atof(argv[3]);
        DMAX = std::atof(argv[4]);
        DEPTH = std::atof(argv[5]);
        AREA = std::atof(argv[6]);
    }
    else {
        std::cout << "Something went wrong! (Error 05)\n";
        exit(5);
    }

    // CRITICAL SPLIT IN FUNCTIONALITY HERE
    if (SHAPE == CUBOID_A || SHAPE == CUBOID_B)
    {
        // TODO: IMPLEMENT THIS BRANCH
        return theBigOne_Cuboid(SHAPE, DIST, DMIN, DMAX, 0.001, DEPTH, AREA);
    }
    else {
        return theBigOne_NonCuboid(SHAPE, DIST, DMIN, DMAX, 0.001, DEPTH, AREA);
    }
}


/* MAIN FUNCTION DEFINITIONS */
int theBigOne_NonCuboid(const int SHAPE, const int DIST, const double DMIN, const double DMAX, const double DSTEP, const double DEPTH, const double AREA) {

    setGeoConst(SHAPE);                 // geoConst fixed for duration of main() based on SHAPE
    setDist(DIST);                      // k and qk fixed for duration of main() based on DIST

    // calculate average diam in the space
    using calculus::integral::monteCarlo;
    double D_ave = monteCarlo(DnDpm3, DMIN, DMAX, 10000) / monteCarlo(nDpm3, DMIN, DMAX, 10000);
    std::cout << "Average Diameter: " << D_ave << std::endl;

    double OccVol = 0.0, Nrocks = 0.0;

    Nrocks = monteCarlo(nDpm3, DMIN, DMAX, 10000) * DEPTH * AREA;
    std::cout << "NSphere: " << Nrocks << std::endl;

    // main arrays to build
    Vector rocks = {};
    CoordVector rockCoords = {}; // xyz
    std::vector<Vector> Sampled_NgtD = {};
    std::vector<Vector> Sampled_FgtD = {};
    Vector Avg_Sampled_NgtD = {}, Avg_Sampled_FgtD = {};

    // sample the function and populate the volume with rocks
    for (int i = 0; i < int(Nrocks); i++) {
        if (i % 1000 == 0) {std::cout << i << '\n';}

        // THESE VARIABLE NAMES ARE HORRIBLE STOP
        double fx = -99.0, fxmax = nDpm3(DMIN);
        double W = fxmax, x = 0.0;
        while (W > fx) {
            W = fxmax * RAND;
            x = (DMAX - DMIN) * RAND + DMIN;
            fx = nDpm3(x);
        }

        // generate random coords for this iteration
        Coord coords = { std::sqrt(AREA) * RAND, std::sqrt(AREA) * RAND, DEPTH * RAND };
        rockCoords.push_back(coords);
        rocks.push_back(x);

        if (SHAPE == CUBE) {
            OccVol += x*x*x;
        }
        else if (SHAPE == SPHERE) {
            OccVol += (4.0/3)*PI*(x*x*x) / 8;
        }
        else {
            std::cout << "Something went wrong! (Error 06)\n";
            exit(6);
        }
    }

    // Now we need to check that if we take any random slice in z, 
    // on average the original distributions are preserved
    for (int i = 0; i < 1000; i++) {
        if (i % 50 == 0) {std::cout << '\t' << i << '\n';}
            
        double z0 = DEPTH * RAND;
        Vector rocksInSlice = {};
        for (int j = 0; j < int(rocks.size()); j++) {
            if (z0 > (rockCoords[j].z - rocks[j]/2.0) && z0 < (rockCoords[j].z + rocks[j]/2.0)) {
                rocksInSlice.push_back(rocks[j]);
            }
        }
            
        Vector SAMPLED_NGTD = CheckNgtD(rocksInSlice, DMIN, DMAX, 0.001, AREA);
        Vector SAMPLED_FGTD = CheckFgtD(rocksInSlice, DMIN, DMAX, 0.001, AREA, SHAPE);
        Sampled_NgtD.push_back(SAMPLED_FGTD);
        Sampled_FgtD.push_back(SAMPLED_NGTD);
                
    }
    Avg_Sampled_NgtD = stats::mean(Sampled_NgtD, 0);
    Avg_Sampled_FgtD = stats::mean(Sampled_FgtD, 0);


    #ifdef WRITE_DATA
        using std::string, std::to_string;
        std::cout << "Writing arrays to files..." << std::endl;
        string shapeStr;
        if (SHAPE == SPHERE) {shapeStr = "SPHERE";}
        else if (SHAPE == CUBE) {shapeStr = "CUBE";}
        else {
            std::cout << "Something went wrong! (Error 07)\n";
            exit(7);
        }
        string fileName = "./data/Rock_Data_" + to_string(AREA) + "XY_" + to_string(DEPTH) + "Z_" + shapeStr + ".csv";
        std::ofstream fout(fileName);
        while (fout.is_open()) {
            fout << "Diameter (m),x (m),y (m), z (m)\n";
            for (int i = 0; i < int(rocks.size()); i++) {
                fout << rocks[i] << ',' << rockCoords[i].x << ',' << rockCoords[i].y << ',' << rockCoords[i].z << '\n';
            }
            fout.close();
        }
        std::cout << "Data written to: " << fileName << std::endl;
    #endif

    return EXIT_SUCCESS;
}

// TODO: CORRECT FUNCTIONALITY FOR CUBOID_A OR CUBOID_B
int theBigOne_Cuboid(const int SHAPE, const int DIST, const double DMIN, const double DMAX, const double DSTEP, const double DEPTH, const double AREA) {

    setGeoConst(SHAPE);                 // geoConst fixed for duration of main() based on SHAPE
    setDist(DIST);                      // k and qk fixed for duration of main() based on DIST

    // calculate average diam in the space
    using calculus::integral::monteCarlo;
    double D_ave = monteCarlo(DnDpm3, DMIN, DMAX, 10000) / monteCarlo(nDpm3, DMIN, DMAX, 10000);
    std::cout << "Average Diameter: " << D_ave << std::endl;

    double OccVol = 0.0, Nrocks = 0.0;

    Nrocks = monteCarlo(nDpm3, DMIN, DMAX, 10000) * DEPTH * AREA;
    std::cout << "NSphere: " << Nrocks << std::endl;

    // main arrays to build
    Vector rocks = {};
    CoordVector rockCoords = {}; // xyz
    std::vector<Vector> Sampled_NgtD = {};
    std::vector<Vector> Sampled_FgtD = {};
    Vector Avg_Sampled_NgtD = {}, Avg_Sampled_FgtD = {};

    // sample the function and populate the volume with rocks
    for (int i = 0; i < int(Nrocks); i++) {
        if (i % 1000 == 0) {std::cout << i << '\n';}

        // THESE VARIABLE NAMES ARE HORRIBLE STOP
        double fx = -99.0, fxmax = nDpm3(DMIN);
        double W = fxmax, x = 0.0;
        while (W > fx) {
            W = fxmax * RAND;
            x = (DMAX - DMIN) * RAND + DMIN;
            fx = nDpm3(x);
        }

        // generate random coords for this iteration
        Coord coords = { std::sqrt(AREA) * RAND, std::sqrt(AREA) * RAND, DEPTH * RAND };
        rockCoords.push_back(coords);
        rocks.push_back(x);

        if (SHAPE == CUBE) {
            OccVol += x*x*x;
        }
        else if (SHAPE == SPHERE) {
            OccVol += (4.0/3)*PI*(x*x*x) / 8;
        }
        else {
            std::cout << "Something went wrong! (Error 06)\n";
            exit(6);
        }
    }

    // Now we need to check that if we take any random slice in z, 
    // on average the original distributions are preserved
    for (int i = 0; i < 1000; i++) {
        if (i % 50 == 0) {std::cout << '\t' << i << '\n';}
            
        double z0 = DEPTH * RAND;
        Vector rocksInSlice = {};
        for (int j = 0; j < int(rocks.size()); j++) {
            if (z0 > (rockCoords[j].z - rocks[j]/2.0) && z0 < (rockCoords[j].z + rocks[j]/2.0)) {
                rocksInSlice.push_back(rocks[j]);
            }
        }
            
        Vector SAMPLED_NGTD = CheckNgtD(rocksInSlice, DMIN, DMAX, 0.001, AREA);
        Vector SAMPLED_FGTD = CheckFgtD(rocksInSlice, DMIN, DMAX, 0.001, AREA, SHAPE);
        Sampled_NgtD.push_back(SAMPLED_FGTD);
        Sampled_FgtD.push_back(SAMPLED_NGTD);
                
    }
    Avg_Sampled_NgtD = stats::mean(Sampled_NgtD, 0);
    Avg_Sampled_FgtD = stats::mean(Sampled_FgtD, 0);


    #ifdef WRITE_DATA
        using std::string, std::to_string;
        std::cout << "Writing arrays to files..." << std::endl;
        string shapeStr;
        if (SHAPE == SPHERE) {shapeStr = "SPHERE";}
        else if (SHAPE == CUBE) {shapeStr = "CUBE";}
        else {
            std::cout << "Something went wrong! (Error 07)\n";
            exit(7);
        }
        string fileName = "./data/Rock_Data_" + to_string(AREA) + "XY_" + to_string(DEPTH) + "Z_" + shapeStr + ".csv";
        std::ofstream fout(fileName);
        while (fout.is_open()) {
            fout << "Diameter (m),x (m),y (m), z (m)\n";
            for (int i = 0; i < int(rocks.size()); i++) {
                fout << rocks[i] << ',' << rockCoords[i].x << ',' << rockCoords[i].y << ',' << rockCoords[i].z << '\n';
            }
            fout.close();
        }
        std::cout << "Data written to: " << fileName << std::endl;
    #endif

    return EXIT_SUCCESS;
}

/* HELPER FUNCTIONS */

/* setGeoConst() */
/* sets GeoConst global var. */
void setGeoConst(int SHAPE) {
    if (SHAPE == SPHERE) GeoConst = 4.0 / PI;
    else if (SHAPE == CUBE) GeoConst = 1.0;
    else if (SHAPE == CUBOID_A) GeoConst = 1.0 / (0.8);
    else if (SHAPE == CUBOID_B) GeoConst = 1.0 / (0.8 * 0.54);
    else { 
        std::cout << "Something went wrong! (Error 01)\n"; 
        exit(1);
    }
}

/* setDist() */
/* sets Distribution (either CE3 or CE4) */
/* sets global vars k, qk */
void setDist(int DIST)
{
    if (DIST == CE3)
    {
        k = 0.0125;
        qk = 1.743;
    }
    else if (DIST == CE4)
    {
        k = 0.0021;
        qk = 0.5648 + 0.01258 / k;
    }
    else {
        std::cout << "Something went wrong! (Error 04)\n";
        exit(4);
    }
}

/* Fk() */
/* k and qk set externally with global vars. Must be set prior to using this func with setDist() */
double Fk(double D) {
    return k * std::exp(-qk * D);
}

/* nDpm2() */
/* We integrate this function, so GeoConst, qk, k need to be set externally! */
/* These consts are controlled by the Distribution and Shape being passed to setDist() and setGeoConst() */
double nDpm2(double D) {
    double k = 0.0135;
    double qk = 1.734;
    return GeoConst * k * qk * std::exp(-qk * D) / (D * D);
}


/* nDpm3() */
/* We integrate this function, so GeoConst, qk, k need to be set externally! */
double nDpm3(double D) {
    return GeoConst * k * qk * std::exp(-qk * D) / (D*D*D);
}

/* DnDpm3() */
/* We integrate this function, so GeoConst, qk, k need to be set externally! */
double DnDpm3(double D) {
    return D * nDpm3(D);
}


/* expintA() */
double expintA(double x) {
    double A = std::log((0.56146 / x + 0.65) * (1.0 + x));
    double B = (x*x*x*x) * std::exp(7.7 * x) * std::pow((2.0 + x), 3.7);
    return std::pow((std::pow(A, -7.7) + B), -0.13);
}

/* NgtD4() */
double NgtD4(double D, const int SHAPE) {    
    setGeoConst(SHAPE);
    double x = 6.555 * D;
    return GeoConst * 0.0137608 * (std::exp(-x) / D - 6.555 * expintA(x));
}

double hpD(double x) {
    double mean = 0.54;
    double std = 0.03;
    double expo = 0.5 * ((x - mean) / std) * ((x - mean) / std);
    return std::exp(-expo);
}

double bpD(double x) {
    double mean = 0.8;
    double std = 0.16;
    double expo = 0.5 * ((x - mean) / std) * ((x - mean) / std);
    return std::exp(-expo);
}


/* NgtD4_Di() */
/* INF requires a numerical approximation for infinity. np.inf aint 'round here partner */
/* N is the number of samples for the integration. Larger INF requires more samples, so choose wisely */
Vector NgtD4_Di(double D, const double INF, const int N) {

    using calculus::integral::monteCarlo;
    setGeoConst(SPHERE); // TODO: CHECK WHETHER THIS IS THE RIGHT MOVE HERE
    // combined two for loops into one from Rejection_Sampling_V6.py
    Vector results = {};
    for (double d = D; d < 3.0; d += 0.001) {
        double integral = monteCarlo(nDpm2, d, INF, N);
        results.push_back(integral);
    }

    return results;
}

/* CheckNgtD() */
/* NON-CUBOID VERSION */
Vector CheckNgtD(Vector& rocks, double Dmin, double Dmax, double Dstep, double A)
{
    Vector N_rocks = {};
    // first off, just iterate based on params
    for (double d = Dmin; d < Dmax; d += Dstep)
    {
        int cuml = 0;
        // "rock" is just another diameter
        for (double rock : rocks) {
            if (rock >= d) cuml++;
        }
        N_rocks.push_back(cuml / A);
    }
    return N_rocks;
}
/* CUBOID VERSION */
/* Coord is just a general struct for a set of x, y, and z values. */
/* It is valid for defining a cuboid shape as well, then. */
Vector CheckNgtD(CoordVector& rocks, double Dmin, double Dmax, double Dstep, double A)
{
    Vector N_rocks = {};
    // first off, just iterate based on params
    for (double d = Dmin; d < Dmax; d += Dstep) {
        int cuml = 0;
        // "rock" is now a Cuboid (Coord), with x, y, and z values
        for (Coord rock : rocks) {
            if (rock.x >= d) cuml++;
        }
        N_rocks.push_back(cuml / A);
    }
    return N_rocks;
}

/* CheckFgtD() */
/* NON-CUBOID VERSION -> this one requires specifying the shape */
Vector CheckFgtD(Vector& rocks, double Dmin, double Dmax, double Dstep, double A, const int SHAPE)
{
    Vector F_area = {};
    for (double d = Dmin; d < Dmax; d += Dstep) {
        double cuml = 0.0;
        // once again, "rock" is just a diameter
        for (double rock : rocks) {
            if (SHAPE == CUBE) cuml += rock * rock;
            else if (SHAPE == SPHERE) cuml += (PI / 4.0) * rock * rock;
            else {
                std::cout << "Something went wrong! (Error 02)\n";
                exit(2);
            }
        }
        F_area.push_back(cuml / A);
    }
    return F_area;
}
/* CUBOID VERSION */
Vector CheckFgtD(CoordVector& rocks, double Dmin, double Dmax, double Dstep, double A)
{
    Vector F_area = {};
    for (double d = Dmin; d < Dmax; d += Dstep) {
        double cuml = 0.0;
        // "rock" is now a Cuboid (Coord), with x, y, and z values
        for (Coord rock : rocks) {
            cuml += rock.x * rock.y;
        }
        F_area.push_back(cuml / A);
    }
    return F_area;
}

/* UNUSED FUNCS */
double F_SVII(double Dmin, double Dmax) {
    double gamma = 1.8;
    double K = 7.9e3 / std::pow(1000., gamma);
    double coeff = gamma * K / (2 - gamma);
    return coeff * (std::pow(Dmax,(2 - gamma)) - std::pow(Dmin, (2 - gamma)));
}
double n_a_SVII(double D) {
    double gamma = 1.8;
    double K = 7.9e3 / std::pow(1000., gamma);
    return gamma * K * std::pow(D, -(gamma + 1));
}
double n_V_SVII(double D) {
    double gamma = 1.8;
    double K = 7.9e3 / std::pow(1000., gamma);
    return gamma * K * std::pow(D, -(gamma + 2));
}
Vector N_SVII(double Dmin, double Dmax) {
    double gamma = 1.8;
    double K = 7.9e3 / std::pow(1000., gamma);
    Vector N = {};
    for (double d = Dmin; d < Dmax; d += 0.001) {
        N.push_back(K * std::pow(d, -gamma));
    }
    return N;
}