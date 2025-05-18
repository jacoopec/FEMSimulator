
#include <vector>
#include "imports.h"

using namespace std;

struct Element {
    int node1, node2;
    double length;
};

int main() {
    double elLength, elArea, E;
    int num_elements;
    std::vector<double> elForces;
    std::vector<bool> is_fixed;
    std::vector<double> elDisplacements;
    readInput("input.txt", elLength, elArea, E, num_elements, elForces, is_fixed, elDisplacements);


    int num_nodes = num_elements + 1;
    double element_length = elLength / num_elements;

    // Global stiffness matrix
    vector<vector<double>> K(num_nodes, vector<double>(num_nodes, 0.0));

    // Assemble K matrix
    assembleStiffnessMatrix(K, E, elArea, element_length, num_elements);

    //Printing K
    printStiffnessMatrix(K);

    // Apply boundary conditions
    applyBoundaryConditions(K, elForces, is_fixed, elDisplacements, num_nodes);

    // Solve system
    solveSystem(K, elForces, elDisplacements, num_nodes);

    // Output displacements
    // Write output
    writeOutput("output.txt", elDisplacements, element_length, E, num_elements, num_nodes);


    return 0;
}
