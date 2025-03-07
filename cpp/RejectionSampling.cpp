#include "RejectionSampling.h"

static const double PI = 4.0 * std::atan(1.0);

double GeoConst;
double qk;
double k;

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
double NgtD4(double D) {    
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

    using calculus::integral::simpsons;
    setGeoConst(SPHERE); // TODO: CHECK WHETHER THIS IS THE RIGHT MOVE HERE
    // combined two for loops into one from Rejection_Sampling_V6.py
    Vector results = {};
    for (double d = D; d < 3.0; d += 0.001) {
        double integral = simpsons(nDpm2, d, INF, N);
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
Vector CheckFgtD(Vector& rocks, double Dmin, double Dmax, double Dstep, double A, const int SHAPE) {
    Vector F_area = {};
    for (double d = Dmin; d < Dmax; d += Dstep) {
        double cuml = 0.0;
        // once again, "rock" is just a diameter
        for (double diam : rocks) {
            if (diam >= d) {
                if (SHAPE == CUBE) cuml += diam * diam;
                else if (SHAPE == SPHERE) cuml += (PI / 4.0) * diam * diam;
                else {
                    std::cout << "Something went wrong! (Error 02)\n";
                    exit(2);
                }
            }
        }
        F_area.push_back(cuml / A);
    }
    return F_area;
}
/* CUBOID VERSION */
Vector CheckFgtD(CoordVector& rocks, double Dmin, double Dmax, double Dstep, double A) {
    Vector F_area = {};
    for (double d = Dmin; d < Dmax; d += Dstep) {
        double cuml = 0.0;
        // "rock" is now a Cuboid (Coord), with x, y, and z values
        for (Coord rock : rocks) {
            if (rock.x >= d) cuml += rock.x * rock.y;
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