#ifndef FEM_UTILS_H
#define FEM_UTILS_H
#include <vector>
#include <iostream>
#include <fstream>

void printStiffnessMatrix(const std::vector<std::vector<double>>& K);
void assembleStiffnessMatrix(std::vector<std::vector<double>>& K, double E, double A, double element_length, int num_elements);
void readInput(const std::string& filename, double& length, double& area, double& E,
               int& num_elements, std::vector<double>& forces,
               std::vector<bool>& is_fixed, std::vector<double>& displacements);
void applyBoundaryConditions(std::vector<std::vector<double>>& K,
                             std::vector<double>& forces,
                             const std::vector<bool>& is_fixed,
                             const std::vector<double>& displacements,
                            int num_nodes);

void solveSystem(std::vector<std::vector<double>>& K,
                 std::vector<double>& forces,
                 std::vector<double>& displacements, int num_nodes);

void writeOutput(const std::string& filename,
                 const std::vector<double>& displacements,
                 double element_length, double E, int num_elements, int num_nodes);


#endif // FEM_UTILS_H