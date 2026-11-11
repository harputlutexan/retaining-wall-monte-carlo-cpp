#include "optimization.hpp"
#include <iostream>
#include <fstream>

/*
 * Retaining wall Monte Carlo simulation and design optimization, C++17.
 *
 * Build (MSVC Developer Command Prompt):
 *   cl /std:c++17 /EHsc /O2 /utf-8 main.cpp /Fe:retaining_wall.exe
 *
 * Run from the folder containing the CSV inputs:
 *   retaining_wall.exe --all .
 *
 * Files:
 *   model.hpp: CSV parsing, distributions, scenario and design generation.
 *   engineering.hpp: transient hydraulics, geotechnical and structural checks, costs.
 *   optimization.hpp: Monte Carlo reliability, verification, CVaR, and selection.
 *   main.cpp: command-line entry point and output writing.
 *
 * The C++ implementation uses mt19937_64, not NumPy PCG64.
 * The resulting random samples are not expected to match a PCG64-based reference.
 */

namespace {

int mainImpl(int argc, char** argv) {
    // The command-line interface runs only the complete simulation.
    if (argc > 3 || (argc > 1 && std::string(argv[1]) != "--all")) {
        std::cerr << "Usage: retaining_wall.exe [--all [csv-folder]]\n";
        return 2;
    }
    const std::string folder = argc > 2 ? argv[2] : ".";

    const rw::Config config = rw::loadConfig(
        folder + "/retaining_wall_inputs_v2.csv",
        folder + "/retaining_wall_settings_v2.csv"
    );
    const auto designs = rw::makeDesigns(config);
    if (designs.size() != 2592)
        throw std::runtime_error("Expected 2592 designs");

    std::cerr << "Generating Monte Carlo scenarios...\n";
    const auto scenarios = rw::makeScenarios(config);
    std::cerr << "Evaluating " << designs.size() << " designs x "
              << scenarios.size() << " scenarios...\n";
    const auto run = rw::calculateAll(config, designs, scenarios, true);

    const std::string outJson = folder + "/retaining_wall_optimization_cpp.json";
    const std::string outCsv  = folder + "/all_design_results_cpp.csv";
    rw::writeResultsJson(outJson, run);
    rw::writeResultsCsv(outCsv, run);

    // Print the same ten reported fields that were written to JSON.
    std::ifstream jsonFile(outJson);
    if (!jsonFile)
        throw std::runtime_error("Cannot reopen result JSON: " + outJson);
    std::cout << "\nMONTE CARLO RESULTS (10 OUTPUT FIELDS)\n";
    std::cout << jsonFile.rdbuf();
    std::cout << "\nFiles: " << outJson << " and " << outCsv << "\n";
    std::cout << "RNG: std::mt19937_64 (not NumPy PCG64)\n";
    return 0;
}

} // anonymous namespace

int main(int argc, char** argv) {
    try {
        return mainImpl(argc, argv);
    } catch (const std::exception& ex) {
        std::cerr << "ERROR: " << ex.what() << "\n";
        return 1;
    }
}
