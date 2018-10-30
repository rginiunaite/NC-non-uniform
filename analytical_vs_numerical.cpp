//analytical solution of reaction equation on a growing domain

#include <numeric>
#include <random>
#include <map>
#include <Eigen/Core>
#include <algorithm> // std::unique_copy
#include <iostream>// writing on a text file
#include <fstream>
#include <math.h>
#include <assert.h>

using namespace std;
using namespace Eigen; // objects VectorXd, MatrixXd


int main() {

    std:vector<double> an = {1, 2, 3, 4};

    double sum_of_elems;

    sum_of_elems = std::accumulate(an.begin(), an.end(), 0);

cout << "sum " << sum_of_elems << endl;


}