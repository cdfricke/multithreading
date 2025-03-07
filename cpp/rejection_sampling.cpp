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
// Date: 24-FEB-2025
// Desc: C++ translation and multithreading version of Rejection_Sampling_V6.py by Payton Linton (linton.93@osu.edu)

#include "RejectionSampling.h"

static const double PI = 4.0 * std::atan(1.0);  // PI

/* MAIN FUNCTION PROTOTYPES */
int theBigOne_NonCuboid(const ProgramParams& P, std::pair<int, int> loopRange);
int theBigOne_Cuboid(const ProgramParams& P, std::pair<int, int> loopRange);

/* MAIN PROGRAM */
int main(int argc, char **argv)
{
    ProgramParams P = setParams(argc, argv);

    // CRITICAL SPLIT IN FUNCTIONALITY HERE
    if (P.SHAPE == CUBOID) {
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

    // Seed the RNG
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> distrib(0.0, 1.0);

    // calculate average diam in the space
    using calculus::integral::simpsons;
    double D_ave = simpsons(DnDpm3, P.DMIN, P.DMAX, 10000) / simpsons(nDpm3, P.DMIN, P.DMAX, 10000);
    std::cout << "Average Diameter: " << D_ave << std::endl;

    double OccVol = 0.0, Nrocks = 0.0;

    Nrocks = simpsons(nDpm3, P.DMIN, P.DMAX, 10000) * P.DEPTH * P.AREA;
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
            OccVol += (4.0/3)*(x*x*x) / 8;
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
        for (size_t j = 0; j < rocks.size(); j++) {
            if (z0 > (rockCoords[j].z - rocks[j]/2.0) && z0 < (rockCoords[j].z + rocks[j]/2.0)) {
                rocksInSlice.push_back(rocks[j]);
            }
        }

        Vector SAMPLED_NGTD = CheckNgtD(rocksInSlice, P.DMIN, P.DMAX, P.DSTEP, P.AREA);
        Vector SAMPLED_FGTD = CheckFgtD(rocksInSlice, P.DMIN, P.DMAX, P.DSTEP, P.AREA, P.SHAPE);
        Sampled_NgtD.push_back(SAMPLED_NGTD);
        Sampled_FgtD.push_back(SAMPLED_FGTD);
                
    }
    Avg_Sampled_NgtD = stats::mean(Sampled_NgtD, 0);
    Avg_Sampled_FgtD = stats::mean(Sampled_FgtD, 0);

    std::cout << "Occupied Volume Fraction = " << OccVol / (P.AREA * P.DEPTH) << std::endl;

    #ifdef WRITE_DATA
        using std::string, std::to_string;
        std::cout << "Writing arrays to files..." << std::endl;
        string shapeStr;
        if (P.SHAPE == SPHERE) {
            shapeStr = "SPHERE";
        }
        else if (P.SHAPE == CUBE) {
            shapeStr = "CUBE";
        }
        else {
            std::cout << "Something went wrong! (Error 07)\n";
            exit(7);
        }
        
        /* WRITE ROCK DATA */
        string fileName = "./data/Rock_Data_" + to_string(int(P.AREA)) + "XY_" + to_string(int(P.DEPTH)) + "Z_" + shapeStr + ".csv";
        std::ofstream fout(fileName);
        while (fout.is_open()) {
            fout << "Diameter (m),x (m),y (m),z (m)\n";
            for (size_t i = 0; i < rocks.size(); i++) {
                fout << rocks[i] << ',' << rockCoords[i].x << ',' << rockCoords[i].y << ',' << rockCoords[i].z << '\n';
            }
            fout.close();
        }
        std::cout << "Data written to: " << fileName << std::endl;

        /* WRITE SAMPLED NgtD DATA */
        fout.open("./data/Avg_Sampled_NgtD.csv", std::ios::out);
        if (!fout.is_open()) {
            std::cout << "Something went wrong! Failed to open file. (Error 08)" << std::endl;
            exit(8);
        }
        for (size_t i = 0; i < Avg_Sampled_NgtD.size(); i++) {
            fout << Avg_Sampled_NgtD[i] << '\n';
        }
        fout.close();
        std::cout << "Data written to: ./data/Avg_Sampled_NgtD.csv (" << Avg_Sampled_NgtD.size() << " values)" << std::endl;

        /* WRITE SAMPLED FgtD DATA */
        fout.open("./data/Avg_Sampled_FgtD.csv", std::ios::out);
        if (!fout.is_open()) {
            std::cout << "Something went wrong! Failed to open file. (Error 08)" << std::endl;
            exit(8);
        }
        for (size_t i = 0; i < Avg_Sampled_FgtD.size(); i++) {
            fout << Avg_Sampled_FgtD[i] << '\n';
        }
        fout.close();
        std::cout << "Data written to: ./data/Avg_Sampled_FgtD.csv (" << Avg_Sampled_FgtD.size() << " values)"<< std::endl;

    #endif

    return EXIT_SUCCESS;
}

int theBigOne_Cuboid(const ProgramParams& P, std::pair<int, int> loopRange) {

    setGeoConst(P.SHAPE);               // geoConst fixed for duration of main() based on SHAPE
    setDist(P.DIST);                    // k and qk fixed for duration of main() based on DIST

    // Seed the RNG
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> distrib(0.0, 1.0);

    // calculate average diam in the space
    using calculus::integral::simpsons;
    double D_ave = simpsons(DnDpm3, P.DMIN, P.DMAX, 10000) / simpsons(nDpm3, P.DMIN, P.DMAX, 10000);
    std::cout << "Average Diameter: " << D_ave << std::endl;

    double OccVol = 0.0, Nrocks = 0.0;
    double INF = 1000.0;

    Nrocks = NgtD4(P.DMIN) * P.DEPTH * P.AREA / (0.54 * D_ave);
    std::cout << "NSphere: " << Nrocks << std::endl;
    Nrocks = simpsons(nDpm3, P.DMIN, INF, 1000000) * P.DEPTH * P.AREA;
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

        OccVol += (x*b*h);
    }   

    // Now we need to check that if we take any random slice in z, 
    // on average the original distributions are preserved
    for (int i = loopRange.first; i < loopRange.second; i++) {
        if (i % 50 == 0) {std::cout << '\t' << i << '\n';}

        double z0 = P.DEPTH * RAND;
        CoordVector rocksInSlice = {};
        for (size_t j = 0; j < rocks.size(); j++) {
            if (z0 > (rockCoords[j].z - rocks[j].z/2.0) && z0 < (rockCoords[j].z + rocks[j].z/2.0)) {
                rocksInSlice.push_back(rocks[j]);
            }
        }

        Vector SAMPLED_NGTD = CheckNgtD(rocksInSlice, P.DMIN, P.DMAX, P.DSTEP, P.AREA);
        Vector SAMPLED_FGTD = CheckFgtD(rocksInSlice, P.DMIN, P.DMAX, P.DSTEP, P.AREA);
        Sampled_NgtD.push_back(SAMPLED_NGTD);
        Sampled_FgtD.push_back(SAMPLED_FGTD);
    }
    Avg_Sampled_NgtD = stats::mean(Sampled_NgtD, 0);
    Avg_Sampled_FgtD = stats::mean(Sampled_FgtD, 0);

    std::cout << "Occupied Volume Fraction = " << OccVol / (P.AREA * P.DEPTH) << std::endl;

    #ifdef WRITE_DATA
        using std::string, std::to_string;
        std::cout << "Writing arrays to files..." << std::endl;
        string shapeStr = "CUBOID";

        /* WRITE ROCK DATA */
        string fileName = "./data/Rock_Data_" + to_string(int(P.AREA)) + "XY_" + to_string(int(P.DEPTH)) + "Z_" + shapeStr + ".csv";
        std::ofstream fout(fileName, std::ios::out);
        if (!fout.is_open()) {
            std::cout << "Something went wrong! Failed to open file. (Error 08)" << std::endl;
            exit(8);
        }
        fout << "Diameter (m),b (m),h (m),x (m),y (m),z (m)\n";
        for (size_t i = 0; i < rocks.size(); i++) {
            fout << rocks[i].x << ',' 
                << rocks[i].y << ',' 
                << rocks[i].z << ',' 
                << rockCoords[i].x << ',' 
                << rockCoords[i].y << ',' 
                << rockCoords[i].z << '\n';
        }
        fout.close();
        std::cout << "Data written to: " << fileName << std::endl;

        /* WRITE SAMPLED NgtD DATA */
        fout.open("./data/Avg_Sampled_NgtD.csv", std::ios::out);
        if (!fout.is_open()) {
            std::cout << "Something went wrong! Failed to open file. (Error 08)" << std::endl;
            exit(8);
        }
        for (size_t i = 0; i < Avg_Sampled_NgtD.size(); i++) {
            fout << Avg_Sampled_NgtD[i] << '\n';
        }
        fout.close();
        std::cout << "Data written to: ./data/Avg_Sampled_NgtD.csv (" << Avg_Sampled_NgtD.size() << " values)" << std::endl;

        /* WRITE SAMPLED FgtD DATA */
        fout.open("./data/Avg_Sampled_FgtD.csv", std::ios::out);
        if (!fout.is_open()) {
            std::cout << "Something went wrong! Failed to open file. (Error 08)" << std::endl;
            exit(8);
        }
        for (size_t i = 0; i < Avg_Sampled_FgtD.size(); i++) {
            fout << Avg_Sampled_FgtD[i] << '\n';
        }
        fout << Avg_Sampled_FgtD.back();
        fout.close();
        std::cout << "Data written to: ./data/Avg_Sampled_FgtD.csv (" << Avg_Sampled_FgtD.size() << " values)" << std::endl;
    
    #endif

    return EXIT_SUCCESS;
}

/* HELPER FUNCTIONS */

