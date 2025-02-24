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
#include "Calculus.h"
#include "Stats.h"
#ifdef WRITE_DATA
    #include <fstream>
    #include <string>
#endif

/* SHAPES */
#define SPHERE 0
#define CUBE 1
#define CUBOID_A 2
#define CUBOID_B 3

/* DISTRIBUTIONS */
#define CE3 0
#define CE4 1

/* MACROS */
#define RAND (double(std::rand())/RAND_MAX)   // RANDOM double between 0.0 and 1.0

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
typedef std::vector<double> Vector;
typedef std::vector<Coord> CoordVector;

/* GLOBAL VARS */
static const double PI = 4.0 * std::atan(1);    // PI
static double GeoConst = 0.0;                   // INIT WITH setGeoConst();
static double k = 0.0;                          // INIT WITH setDist();
static double qk = 0.0;                         // INIT WITH setDist();

/* HELPER FUNCTION PROTOTYPES */
ProgramParams setParams(int argc, char** argv);
void setGeoConst(int SHAPE);
void setDist(int DIST);
double Fk(double D);
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
int theBigOne_NonCuboid(const ProgramParams& P, std::pair<int, int> loopRange);
int theBigOne_Cuboid(const ProgramParams& P, std::pair<int, int> loopRange);

/* MAIN PROGRAM */
int main(int argc, char **argv)
{
    ProgramParams P = setParams(argc, argv);

    // CRITICAL SPLIT IN FUNCTIONALITY HERE
    if (P.SHAPE == CUBOID_A || P.SHAPE == CUBOID_B)
    {
        return theBigOne_Cuboid(P, {0,1000});
    }
    else {
        return theBigOne_NonCuboid(P, {0,1000});
    }
}


/* MAIN FUNCTION DEFINITIONS */
int theBigOne_NonCuboid(const ProgramParams &P, std::pair<int, int> loopRange) {

    setGeoConst(P.SHAPE);                 // geoConst fixed for duration of main() based on SHAPE
    setDist(P.DIST);                      // k and qk fixed for duration of main() based on DIST

    // calculate average diam in the space
    using calculus::integral::monteCarlo;
    double D_ave = monteCarlo(DnDpm3, P.DMIN, P.DMAX, 10000) / monteCarlo(nDpm3, P.DMIN, P.DMAX, 10000);
    std::cout << "Average Diameter: " << D_ave << std::endl;

    double OccVol = 0.0, Nrocks = 0.0;

    Nrocks = monteCarlo(nDpm3, P.DMIN, P.DMAX, 10000) * P.DEPTH * P.AREA;
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
        double fx = -99.0, fxmax = nDpm3(P.DMIN);
        double W = fxmax, x = 0.0;
        while (W > fx) {
            W = fxmax * RAND;
            x = (P.DMAX - P.DMIN) * RAND + P.DMIN;
            fx = nDpm3(x);
        }

        // generate random coords for this iteration
        Coord coords = { std::sqrt(P.AREA) * RAND, std::sqrt(P.AREA) * RAND, P.DEPTH * RAND };
        rockCoords.push_back(coords);
        rocks.push_back(x);

        if (P.SHAPE == CUBE) {
            OccVol += x*x*x;
        }
        else if (P.SHAPE == SPHERE) {
            OccVol += (4.0/3)*PI*(x*x*x) / 8;
        }
        else {
            std::cout << "Something went wrong! (Error 06)\n";
            exit(6);
        }
    }

    // Now we need to check that if we take any random slice in z, 
    // on average the original distributions are preserved
    for (int i = loopRange.first; i < loopRange.second; i++) {
        if (i % 50 == 0) {std::cout << '\t' << i << '\n';}
            
        double z0 = P.DEPTH * RAND;
        Vector rocksInSlice = {};
        for (int j = 0; j < int(rocks.size()); j++) {
            if (z0 > (rockCoords[j].z - rocks[j]/2.0) && z0 < (rockCoords[j].z + rocks[j]/2.0)) {
                rocksInSlice.push_back(rocks[j]);
            }
        }

        Vector SAMPLED_NGTD = CheckNgtD(rocksInSlice, P.DMIN, P.DMAX, P.DSTEP, P.AREA);
        Vector SAMPLED_FGTD = CheckFgtD(rocksInSlice, P.DMIN, P.DMAX, P.DSTEP, P.AREA, P.SHAPE);
        Sampled_NgtD.push_back(SAMPLED_FGTD);
        Sampled_FgtD.push_back(SAMPLED_NGTD);
                
    }
    Avg_Sampled_NgtD = stats::mean(Sampled_NgtD, 0);
    Avg_Sampled_FgtD = stats::mean(Sampled_FgtD, 0);

    std::cout << "Occupied Volume Fraction = " << OccVol / (P.AREA * P.DEPTH) << std::endl;

    #ifdef WRITE_DATA
        using std::string, std::to_string;
        std::cout << "Writing arrays to files..." << std::endl;
        string shapeStr;
        if (P.SHAPE == SPHERE)
        {
            shapeStr = "SPHERE";
        }
        else if (P.SHAPE == CUBE)
        {
            shapeStr = "CUBE";
        }
        else {
            std::cout << "Something went wrong! (Error 07)\n";
            exit(7);
        }
        string fileName = "./data/Rock_Data_" + to_string(int(P.AREA)) + "XY_" + to_string(int(P.DEPTH)) + "Z_" + shapeStr + ".csv";
        std::ofstream fout(fileName);
        while (fout.is_open()) {
            fout << "Diameter (m),x (m),y (m),z (m)\n";
            for (int i = 0; i < int(rocks.size()); i++) {
                fout << rocks[i] << ',' << rockCoords[i].x << ',' << rockCoords[i].y << ',' << rockCoords[i].z << '\n';
            }
            fout.close();
        }
        std::cout << "Data written to: " << fileName << std::endl;

    #endif

    return EXIT_SUCCESS;
}

int theBigOne_Cuboid(const ProgramParams& P, std::pair<int, int> loopRange) {

    setGeoConst(P.SHAPE);               // geoConst fixed for duration of main() based on SHAPE
    setDist(P.DIST);                    // k and qk fixed for duration of main() based on DIST

    // calculate average diam in the space
    using calculus::integral::monteCarlo;
    double D_ave = monteCarlo(DnDpm3, P.DMIN, P.DMAX, 10000) / monteCarlo(nDpm3, P.DMIN, P.DMAX, 10000);
    std::cout << "Average Diameter: " << D_ave << std::endl;

    double OccVol = 0.0, Nrocks = 0.0;
    double INF = 1000.0;
    Nrocks = monteCarlo(nDpm3, P.DMIN, INF, 100000) * P.DEPTH * P.AREA;
    std::cout << "NSphere: " << Nrocks << std::endl;

    // main arrays to build
    CoordVector rocks = {};
    CoordVector rockCoords = {}; // xyz
    std::vector<Vector> Sampled_NgtD = {};
    std::vector<Vector> Sampled_FgtD = {};
    Vector Avg_Sampled_NgtD = {}, Avg_Sampled_FgtD = {};

    // sample the function and populate the volume with rocks
    for (int i = 0; i < int(Nrocks); i++) {
        if (i % 1000 == 0) {std::cout << i << '\n';}

        // This loop chooses x
        double fx = -99.0, fxmax = nDpm3(P.DMIN);
        double W = fxmax, x = 0.0;
        while (W > fx) {
            W = fxmax * RAND;
            x = (P.DMAX - P.DMIN) * RAND + P.DMIN;
            fx = nDpm3(x);
        }

        // These loops choose h and b
        double b = 0.0, bx = 0.0, h = 0.0, hx = 0.0;
        double fb = -99.0, fh = -99.0;
        double bymax = 1.0, hymax = 1.0;
        while (bymax > fb) {
            bx = RAND;
            bymax = RAND;
            fb = bpD(bx);
        }
        while (hymax > fh)
        {
            hx = RAND;
            hymax = RAND;
            fh = hpD(hx);
        }
        h = hx * x;
        b = bx * x;
        // create rock and give it random coords within the volume
        Coord ROCK = {x, b, h};
        Coord coords = {std::sqrt(P.AREA) * RAND, std::sqrt(P.AREA) * RAND, P.DEPTH * RAND};
        rockCoords.push_back(coords);
        rocks.push_back(ROCK);

        OccVol += OccVol + (x*b*h);
    }

    // Now we need to check that if we take any random slice in z, 
    // on average the original distributions are preserved
    for (int i = loopRange.first; i < loopRange.second; i++) {
        if (i % 50 == 0) {std::cout << '\t' << i << '\n';}

        double z0 = P.DEPTH * RAND;
        CoordVector rocksInSlice = {};
        for (int j = 0; j < int(rocks.size()); j++) {
            if (z0 > (rockCoords[j].z - rocks[j].z/2.0) && z0 < (rockCoords[j].z + rocks[j].z/2.0)) {
                rocksInSlice.push_back(rocks[j]);
            }
        }

        Vector SAMPLED_NGTD = CheckNgtD(rocksInSlice, P.DMIN, P.DMAX, P.DSTEP, P.AREA);
        Vector SAMPLED_FGTD = CheckFgtD(rocksInSlice, P.DMIN, P.DMAX, P.DSTEP, P.AREA);
        Sampled_NgtD.push_back(SAMPLED_FGTD);
        Sampled_FgtD.push_back(SAMPLED_NGTD);
    }
    Avg_Sampled_NgtD = stats::mean(Sampled_NgtD, 0);
    Avg_Sampled_FgtD = stats::mean(Sampled_FgtD, 0);

    std::cout << "Occupied Volume Fraction = " << OccVol / (P.AREA * P.DEPTH) << std::endl;

    #ifdef WRITE_DATA
        using std::string, std::to_string;
        std::cout << "Writing arrays to files..." << std::endl;
        string shapeStr;
        if (P.SHAPE == CUBOID_A)
        {
            shapeStr = "CUBOID_A";
        }
        else if (P.SHAPE == CUBOID_B)
        {
            shapeStr = "CUBOID_B";
        }
        else {
            std::cout << "Something went wrong! (Error 07)\n";
            exit(7);
        }
        string fileName = "./data/Rock_Data_" + to_string(int(P.AREA)) + "XY_" + to_string(int(P.DEPTH)) + "Z_" + shapeStr + ".csv";
        std::ofstream fout(fileName);
        while (fout.is_open()) {
            fout << "Diameter (m),b (m),h (m),x (m),y (m),z (m)\n";
            for (int i = 0; i < int(rocks.size()); i++) {
                fout << rocks[i].x << ',' 
                    << rocks[i].y << ',' 
                    << rocks[i].z << ',' 
                    << rockCoords[i].x << ',' 
                    << rockCoords[i].y << ',' 
                    << rockCoords[i].z << '\n';
            }
            fout.close();
        }
        std::cout << "Data written to: " << fileName << std::endl;
    #endif

    return EXIT_SUCCESS;
}

/* HELPER FUNCTIONS */

/* setParams() */
ProgramParams setParams(int argc, char** argv) {

    ProgramParams P;

    if (argc == 1) { // defaults
        std::cout << "Using DEFAULT Parameters..." << std::endl;
    }
    else if (argc == 8) { // user specified
        std::cout << "Using GIVEN Parameters..." << std::endl;
        P.SHAPE = std::atoi(argv[1]);
        P.DIST = std::atoi(argv[2]);
        P.DMIN = std::atof(argv[3]);
        P.DMAX = std::atof(argv[4]);
        P.DSTEP = std::atof(argv[5]);
        P.DEPTH = std::atof(argv[6]);
        P.AREA = std::atof(argv[7]);
        return P;
    }
    else {
        std::cout << "Something went wrong! (Error 05)\n";
        exit(5);
    }
    return P;
}

/* setGeoConst() */
/* sets GeoConst global var. */
void setGeoConst(int SHAPE)
{
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
    if (DIST == CE3) {
        k = 0.0125;
        qk = 1.743;
    }
    else if (DIST == CE4) {
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