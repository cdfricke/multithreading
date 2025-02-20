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

/*
trapezoidal()
    numerical integration technique which integrates a passed function over 
*/
double trapezoidal( double (*func)(double), const double a, const double b, int n)
{
    const double h = (b - a) / n;
    double sum = 0.5 * (func(a) + func(b));
    for (int i = 1; i < n; ++i) {
        sum += func(a + i * h);
    }
    return sum * h;
}

/*
This function is the cumulative fractional area covered by
rocks with diameter > D per area
Basically says look at the ground, every square area Fk(D) of the
area is covered by rocks that are at least D in diameter
This is the distribution reported in Wu et al. (2021)
*/
// by Payton Linton
double Fk(const double D) {
    
    const double k = 0.0021;
    const double qk = 0.5648 + (0.01258 / k);
    return k * std::exp(-qk * D);
}

double Fk_Di(const double D) {
    /*
    This function is the cumulative fractional area covered by
    rocks with diameter > D per area
    Basically says look at the ground, every square area Fk(D) of the
    area is covered by rocks that are at least D in diameter
    This is the distribution reported in Wu et al. (2021)
    */
    const double k = 0.0125;
    const double qk = 1.743;
    return k * std::exp(-qk * D);
}

std::vector<double> NgtD4_Di(const double D) {
    std::vector<double> Diam_Array = {};
    for (double d = D; d < 3.0; d += 0.001) {
        Diam_Array.push_back(d);
    }

    std::vector<double> N = {};
    for (int i = 0; i < Diam_Array.size(); i++) {
        double integral = quad(nDpm2, Diam_Array[i], std::inf )
        N.push_back()
    }

}

int main() {



    return EXIT_SUCCESS;
}