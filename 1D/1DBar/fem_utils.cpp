#include "imports.h"

#include <iostream>

void printStiffnessMatrix(const std::vector<std::vector<double>>& K) {
    std::cout << "Global Stiffness Matrix [K]:" << std::endl;
    int num_nodes = K.size();
    for (int i = 0; i < num_nodes; ++i) {
        for (int j = 0; j < num_nodes; ++j) {
            std::cout << K[i][j] << "\t";
        }
        std::cout << std::endl;
    }
}

void assembleStiffnessMatrix(std::vector<std::vector<double>>& K, double E, double A, double element_length, int num_elements) {
    for (int i = 0; i < num_elements; ++i) {
        double k = E * A / element_length;
        K[i][i]     += k;
        K[i][i + 1] -= k;
        K[i + 1][i] -= k;
        K[i + 1][i + 1] += k;
    }
}

void readInput(const std::string& filename, double& length, double& area, double& E,
               int& num_elements, std::vector<double>& forces,
               std::vector<bool>& is_fixed, std::vector<double>& displacements){
    std::ifstream infile("input.txt");
    
    infile >> length >> num_elements >> area >> E;

    int num_nodes = num_elements + 1;
    double element_length = length / num_elements;

    for (int i = 0; i < num_nodes; ++i)
        infile >> forces[i];

    int bc_node;
    double bc_value;
    while (infile >> bc_node >> bc_value) {
        is_fixed[bc_node] = true;
        displacements[bc_node] = bc_value;
    }
    infile.close();
}


void applyBoundaryConditions(std::vector<std::vector<double>>& K,
                             std::vector<double>& forces,
                             const std::vector<bool>& is_fixed,
                             const std::vector<double>& displacements,
                            int num_nodes){
    for (int i = 0; i < num_nodes; ++i) {
        if (is_fixed[i]) {
            for (int j = 0; j < num_nodes; ++j) {
                K[i][j] = 0;
                K[j][i] = 0;
            }
            K[i][i] = 1;
            forces[i] = displacements[i];
        }
    }
}


void solveSystem(std::vector<std::vector<double>>& K,
                 std::vector<double>& forces,
                 std::vector<double>& displacements, int num_nodes){
    // Solve using Gaussian elimination
    // Forward elimination
    for (int i = 0; i < num_nodes; ++i) {
        double pivot = K[i][i];
        for (int j = 0; j < num_nodes; ++j)
            K[i][j] /= pivot;
        forces[i] /= pivot;

        for (int k = i + 1; k < num_nodes; ++k) {
            double factor = K[k][i];
            for (int j = 0; j < num_nodes; ++j)
                K[k][j] -= factor * K[i][j];
            forces[k] -= factor * forces[i];
        }
    }

    // Back substitution
    for (int i = num_nodes - 1; i >= 0; --i) {
        displacements[i] = forces[i];
        for (int j = i + 1; j < num_nodes; ++j)
            displacements[i] -= K[i][j] * displacements[j];
    }
}

void writeOutput(const std::string& filename,
                 const std::vector<double>& displacements,
                 double element_length, double E, int num_elements, int num_nodes){
    std::ofstream outfile("output.txt");
    // Output displacements
    outfile << "Nodal Displacements:\n";
    for (int i = 0; i < num_nodes; ++i)
        outfile << "Node " << i << ": " << displacements[i] << "\n";

    // Compute stresses
    outfile << "\nElement Stresses:\n";
    for (int i = 0; i < num_elements; ++i) {
        double strain = (displacements[i + 1] - displacements[i]) / element_length;
        double stress = E * strain;
        outfile << "Element " << i << ": " << stress << "\n";
    }
    
    
    outfile.close();
}
