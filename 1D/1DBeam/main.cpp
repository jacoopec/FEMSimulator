// main.cpp
#include <iostream>
#include <fstream>
#include <vector>
#include "fem_utils.h"

int main() {
    std::ifstream infile("input.txt");
    double length, EI;
    int num_elements;
    infile >> length >> num_elements >> EI;

    int num_nodes = num_elements + 1;
    int total_dofs = 2 * num_nodes;
    double element_length = length / num_elements;

    std::vector<double> forces(total_dofs);
    for (int i = 0; i < total_dofs; ++i)
        infile >> forces[i];

    std::vector<bool> is_fixed(total_dofs, false);
    std::vector<double> displacements(total_dofs, 0.0);

    int bc_dof;
    double bc_value;
    while (infile >> bc_dof >> bc_value) {
        is_fixed[bc_dof] = true;
        displacements[bc_dof] = bc_value;
    }
    infile.close();

    std::vector<std::vector<double>> K(total_dofs, std::vector<double>(total_dofs, 0.0));

    for (int i = 0; i < num_elements; ++i)
        assembleBeamElement(K, EI, element_length, i);

    printStiffnessMatrix(K);
    applyBoundaryConditions(K, forces, is_fixed, displacements);
    solveSystem(K, forces, displacements);
    writeOutput("output.txt", displacements, element_length, EI, num_elements);

    return 0;
}
