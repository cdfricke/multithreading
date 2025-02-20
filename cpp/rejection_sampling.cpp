/*
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
#include <thread>
#include <cmath>
#include <vector>
#include <string>
#include "Calculus.h"

/* define shapes */
#define SPHERE 0
#define TORUS 1
#define CUBE 2
#define PYRAMID 3
#define CUBOID_A 4
#define CUBOID_B 5
/* define distributions */
#define CE3 0
#define CE4 1

static const double PI = 4.0 * std::atan(1);    // PI
static double GeoConst = 0.0;                   // value set by setGeoConst()

/* FUNCTION PROTOTYPES */
void setGeoConst(const int SHAPE);
double Fk(const double D);
double Fk_Di(const double D);
double nDpm2(double D);
double DnDpm2(double D);
std::vector<double> NgtD4_Di(const double D, const double INF, const int N);

double f(double x) {
    return std::exp(-x);
}

int main() {

    setGeoConst(SPHERE);

    std::cout << calculus::integral::monteCarlo(f, 0., 1., 1000) << std::endl;

    return EXIT_SUCCESS;
}



/* FUNCTION DEFINITIONS */


/* setGeoConst() */
/* sets GeoConst global var. */
void setGeoConst(const int SHAPE) {
    if (SHAPE == SPHERE || SHAPE == PYRAMID) GeoConst = 4.0 / PI;
    else if (SHAPE == CUBE || SHAPE == PYRAMID) GeoConst = 1.0;
    else if (SHAPE == CUBOID_A) GeoConst = 1.0 / (0.8);
    else if (SHAPE == CUBOID_B) GeoConst = 1.0 / (0.8 * 0.54);
    else { 
        std::cout << "Something went wrong! (Error 01)\n"; 
        exit(1);
    }
}


/* Fk() */
double Fk(const double D) {
    const double k = 0.0021;
    const double qk = 0.5648 + (0.01258 / k);
    return k * std::exp(-qk * D);
}


/* Fk_Di() */
double Fk_Di(const double D) {
    const double k = 0.0125;
    const double qk = 1.743;
    return k * std::exp(-qk * D);
}


/* nDpm2() */
double nDpm2(double D) {
    const double k = 0.0135;
    const double qk = 1.734;
    return GeoConst * k * qk * std::exp(-qk * D) / (D * D);
}


/* DnDpm2() */
double DnDpm2(double D) {
    return D * nDpm2(D);
}


/* NgtD4_Di() */
/* INF requires a numerical approximation for infinity. np.inf aint 'round here partner */
/* ITER is the number of samples for the monteCarlo integration */
std::vector<double> NgtD4_Di(const double D, const double INF, const int N) {

    using calculus::integral::monteCarlo;

    std::vector<double> Diam_Array = {};
    for (double d = D; d < 3.0; d += 0.001) {
        Diam_Array.push_back(d);
    }

    std::vector<double> results = {};
    for (int i = 0; i < Diam_Array.size(); i++) {
        double integral = monteCarlo(nDpm2, Diam_Array[i], INF, N);
        results.push_back(integral);
    }
    return results;
}