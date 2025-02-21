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

/* DEFINE SHAPES */
#define SPHERE 0
#define CUBE 1
#define CUBOID_A 2
#define CUBOID_B 3

/* DEFINE DISTRIBUTIONS */
#define CE3 0
#define CE4 1

struct Cuboid
{
    double x, y, z;
};

typedef std::vector<double> Vector;
typedef std::vector<Cuboid> CuboidVector;

/* DEFINE GLOBAL VARS */
static const double PI = 4.0 * std::atan(1);    // PI
static double GeoConst = 0.0;                   // controlled by individual funcs with setGeoConst()

/* FUNCTION PROTOTYPES */
void setGeoConst(int SHAPE);
double Fk(double D);
double Fk_Di(double D);
double nDpm2(double D);
double nDpm3(double D, const int SHAPE, const int DIST);
double DnDpm3(double D, const int SHAPE, const int DIST);
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
Vector CheckNgtD(CuboidVector& rocks, double Dmin, double Dmax, double Dstep, double A);
Vector CheckFgtD(Vector& rocks, double Dmin, double Dmax, double Dstep, double A, const int SHAPE);
Vector CheckFgtD(CuboidVector& rocks, double Dmin, double Dmax, double Dstep, double A);


/* MAIN PROGRAM */
int main(int argc, char **argv)
{
    // Program Parameter control
    int SHAPE;
    double DMIN, VOL;
    if (argc == 1) {
        // DEFAULT PARAMS
        SHAPE = SPHERE;
        DMIN = 0.01;
        VOL = 100.*100.*100.;
    }
    else if (argc == 4)
    {
        // PASSED PARAMS
        std::cout << "Using Passed Parameters..." << std::endl;
        SHAPE = std::atoi(argv[1]);
        DMIN = std::atof(argv[2]);
        VOL = std::atof(argv[3]);
    }
    

    return EXIT_SUCCESS;
}

/* FUNCTION DEFINITIONS */

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


/* Fk() */
double Fk(double D) {
    double k = 0.0021;
    double qk = 0.5648 + (0.01258 / k);
    return k * std::exp(-qk * D);
}


/* Fk_Di() */
double Fk_Di(double D) {
    double k = 0.0135;
    double qk = 1.734;
    return k * std::exp(-qk * D);
}


/* nDpm2() */
/* We integrate this function, so GeoConst needs to be set externally! */
double nDpm2(double D) {
    double k = 0.0135;
    double qk = 1.734;
    return GeoConst * k * qk * std::exp(-qk * D) / (D * D);
}


/* nDpm3() */
double nDpm3(double D, const int SHAPE, const int DIST) {
    setGeoConst(SHAPE);
    if (DIST == CE3)
    {
        double k = 0.0135;
        double qk = 1.734;
        return GeoConst * k * qk * std::exp(-qk * D) / (D*D*D);
    }
    else if (DIST == CE4)
    {
        double k = 0.0021;
        double qk = 0.5648 + 0.01258 / k;
        return GeoConst * k * qk * std::exp(-qk * D) / (D*D*D);
    }
    else {
        std::cout << "Something went wrong! (Error 03)\n";
        exit(3);
    }
}

/* DnDpm3() */
double DnDpm3(double D, const int SHAPE, const int DIST) {
    return D * nDpm3(D, SHAPE, DIST);
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
Vector CheckNgtD(CuboidVector& rocks, double Dmin, double Dmax, double Dstep, double A)
{
    Vector N_rocks = {};
    // first off, just iterate based on params
    for (double d = Dmin; d < Dmax; d += Dstep) {
        int cuml = 0;
        // "rock" is now a Cuboid, with x, y, and z values
        for (Cuboid rock : rocks) {
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
Vector CheckFgtD(CuboidVector& rocks, double Dmin, double Dmax, double Dstep, double A)
{
    Vector F_area = {};
    for (double d = Dmin; d < Dmax; d += Dstep) {
        double cuml = 0.0;
        // "rock" is now a Cuboid
        for (Cuboid rock : rocks) {
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