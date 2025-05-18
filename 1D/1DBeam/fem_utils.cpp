#include "fem_utils.h"
#include <iostream>
#include <fstream>
#include <vector>

void printStiffnessMatrix(const std::vector<std::vector<double>>& K) {
    std::cout << "Global Stiffness Matrix [K]:\n";
    for (const auto& row : K) {
        for (double val : row)
            std::cout << val << "\t";
        std::cout << std::endl;
    }
}

void assembleBeamElement(std::vector<std::vector<double>>& K, double EI, double L, int element_index) {
    double L2 = L * L;
    double L3 = L2 * L;
    double coeff = EI / L3;

    std::vector<std::vector<double>> Ke = {
        { 12,   6*L,   -12,   6*L  },
        { 6*L,  4*L2,  -6*L,  2*L2 },
        {-12,  -6*L,   12,   -6*L  },
        { 6*L,  2*L2, -6*L,  4*L2 }
    };

    int dof_start = element_index * 2;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            K[dof_start + i][dof_start + j] += coeff * Ke[i][j];
}

void applyBoundaryConditions(std::vector<std::vector<double>>& K,
                              std::vector<double>& forces,
                              const std::vector<bool>& is_fixed,
                              const std::vector<double>& displacements) {
    int dofs = K.size();
    for (int i = 0; i < dofs; ++i) {
        if (is_fixed[i]) {
            for (int j = 0; j < dofs; ++j) {
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
                 std::vector<double>& displacements) {
    int n = K.size();

    for (int i = 0; i < n; ++i) {
        double pivot = K[i][i];
        for (int j = 0; j < n; ++j)
            K[i][j] /= pivot;
        forces[i] /= pivot;

        for (int k = i + 1; k < n; ++k) {
            double factor = K[k][i];
            for (int j = 0; j < n; ++j)
                K[k][j] -= factor * K[i][j];
            forces[k] -= factor * forces[i];
        }
    }

    for (int i = n - 1; i >= 0; --i) {
        displacements[i] = forces[i];
        for (int j = i + 1; j < n; ++j)
            displacements[i] -= K[i][j] * displacements[j];
    }
}

void writeOutput(const std::string& filename,
                 const std::vector<double>& displacements,
                 double element_length, double EI, int num_elements) {
    std::ofstream outfile(filename);
    int num_nodes = (int)(displacements.size() / 2);

    outfile << "Nodal Results (displacement, rotation):\n";
    for (int i = 0; i < num_nodes; ++i) {
        outfile << "Node " << i << ": u = " << displacements[2 * i]
                << ", theta = " << displacements[2 * i + 1] << "\n";
    }
    outfile.close();
}
