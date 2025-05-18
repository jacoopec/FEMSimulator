#ifndef FEM_UTILS_H
#define FEM_UTILS_H

#include <vector>
#include <string>

void printStiffnessMatrix(const std::vector<std::vector<double>>& K);

void assembleBeamElement(std::vector<std::vector<double>>& K, double EI, double L, int element_index);

void applyBoundaryConditions(std::vector<std::vector<double>>& K,
                              std::vector<double>& forces,
                              const std::vector<bool>& is_fixed,
                              const std::vector<double>& displacements);

void solveSystem(std::vector<std::vector<double>>& K,
                 std::vector<double>& forces,
                 std::vector<double>& displacements);

void writeOutput(const std::string& filename,
                 const std::vector<double>& displacements,
                 double element_length, double EI, int num_elements);

#endif // FEM_UTILS_H
