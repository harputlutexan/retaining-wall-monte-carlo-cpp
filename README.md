# Probabilistic Retaining Wall Design Optimization (C++17)

A reproducible, end-to-end **Monte Carlo engineering simulation and robust design optimization** example for a reinforced-concrete retaining wall. The application enumerates 2,592 candidate designs, evaluates each under the same 5,000 uncertain scenarios, applies reliability and tail-risk constraints, and selects a final design using minimax cost regret.

> **Scope:** This project uses synthetic inputs and simplified, explicitly prescribed engineering models. It is a computational demonstration, **not a validated, code-compliant retaining wall design** and not a substitute for project-specific geotechnical investigation, engineering review, or professional design software.

## Overview

| Item | Configuration |
| --- | --- |
| Language | C++17, standard library only |
| Input format | Two CSV files |
| Fixed engineering inputs | 21 |
| Design variables | 7 |
| Uncertain inputs | 15 |
| Reinforcement detailing input sets | 2 |
| Analysis/optimization setting rows | 105 |
| Candidate designs | 2,592 |
| Common Monte Carlo scenarios | 5,000 |
| Design-scenario evaluations | 12,960,000 |
| Random engine | `std::mt19937_64` with seed `28491`, explicit uniform/normal transforms |
| Primary outputs | One JSON summary and one CSV containing all designs |

### Example results

With the explicit uniform and normal sampling algorithms (seed `28491`), the demonstration produces:

| Metric | Result |
| --- | ---: |
| Baseline-feasible designs | 96 |
| Nominal minimum expected-cost design | #2100 |
| Derived verification scenario | #4514 |
| Robust-feasible designs | 78 |
| Final design | **#2134** |
| Final design expected total cost | **$3,172.20** |
| Final design maximum cost regret | **$0.00** |

These are **example simulation outputs**, not field-verified construction estimates. The final design has zero regret here because it has the lowest deterministic cost in all three prescribed price worlds among the robust-feasible alternatives. Minimax regret does not necessarily produce a different choice from minimum cost for every data set.

## Build and run

Run commands from the repository root, where the two input CSV files are located.

### Linux / macOS (GCC or Clang)

```bash
g++ -std=c++17 -O2 -Wall -Wextra -Wpedantic main.cpp -o retaining_wall
./retaining_wall --all .
```

### Windows (Visual Studio Developer Command Prompt)

```bat
build_windows.bat
retaining_wall.exe --all .
```

Equivalent manual compilation:

```bat
cl /std:c++17 /EHsc /O2 /utf-8 main.cpp /Fe:retaining_wall.exe
retaining_wall.exe --all .
```

### Optional CMake build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
./build/retaining_wall --all .
```

For multi-configuration generators such as Visual Studio, the executable is typically at `build/Release/retaining_wall.exe` rather than `build/retaining_wall`.

The program accepts `--all` and an optional path to the input CSV folder. The application writes these files into that same folder:

- `retaining_wall_optimization_cpp.json`: ten summary metrics.
- `all_design_results_cpp.csv`: one row per candidate design, including reliability, expected cost, feasibility, verification utilization, CVaR95, and regret.

Reference outputs from a successful demonstration run are included under `sample_results/`. Generated result files in the repository root are ignored by Git.

## Input data

### `retaining_wall_inputs_v2.csv`

Contains **45 input rows** organized as follows:

- **21 fixed parameters:** wall height, material and geometric constants, unit weights, reinforcement requirements, and capacity factors.
- **7 design variables:** base width, toe ratio, stem thickness, base thickness, shear-key depth, stem reinforcement target, and base reinforcement target.
- **15 uncertain parameters:** soil properties, loading, material strengths, hydraulic inputs, passive-resistance availability, and cost multipliers.
- **2 detailing option sets:** allowed reinforcement diameters and spacings.

### `retaining_wall_settings_v2.csv`

Contains **105 settings and methodological descriptions**, including stage boundary conditions, engineering criteria, sampling configuration, reliability requirements, and cost-selection rules. Some rows describe prescribed formulas rather than executable configuration variables; the C++ code also contains model-specific formulas. Changing arbitrary descriptive rows will therefore not automatically change the implemented equations.

The configured numerical values include a 5,000-scenario sample, a fixed seed, a minimum baseline success probability of 85%, a verification order-statistic rank of 4,400, and a 95% construction sliding tail-risk measure.

## Methods and workflow

1. **CSV parsing and configuration.** Read both CSV files, validate expected categories, and store fixed inputs, allowable design choices, uncertain distribution definitions, and analysis settings.
2. **Design space enumeration.** Construct the Cartesian product of seven design variables, producing 2,592 designs. Select provided reinforcement from the allowed bar combinations and check deterministic geometry/reinforcement validity.
3. **Uncertainty modeling.** Define 15 uncertain variables with four distribution families:
   - **Truncated normal:** normal draws rejected outside prescribed lower and upper bounds.
   - **Truncated lognormal:** positive, right-skewed draws derived from arithmetic mean and coefficient of variation, with rejection outside bounds.
   - **Triangular:** inverse-CDF sampling using minimum, mode, and maximum.
   - **Categorical:** choose among discrete states using their cumulative probabilities.
4. **Monte Carlo sampling.** Convert the top 53 bits of each raw `mt19937_64` output to a uniform variate in `[0,1)`; draw normal variates with the Marsaglia polar rejection algorithm (without caching the second variate). Generate 5,000 scenarios once. Each scenario contains 15 independently sampled model inputs. Apply the **same scenario set to every design** (common random numbers) to make candidate comparisons more consistent.
5. **Transient hydraulic analysis.** Solve a one-dimensional pressure-diffusion model with an 11-node spatial grid. Use **Crank-Nicolson** time stepping and the **Thomas algorithm** for the nine interior unknowns. Propagate the construction-stage pressure state into service and analytically integrate the piecewise-linear profiles for water forces and moments. Cache hydraulic results that share the same base width and shear-key depth.
6. **Geotechnical and structural evaluation.** Calculate **Rankine active/passive earth pressures**, sliding resistance from force equilibrium, overturning from moment equilibrium, allowable bearing stress, and full/partial foundation contact pressure. Calculate stem/toe/heel bending moments and shears by load integration, check reinforced-concrete flexural/shear capacity, and evaluate cracked-section stiffness and service rotation.
7. **Reliability and expected cost.** Evaluate each design in all 5,000 scenarios. A scenario passes only if the required construction/service geotechnical, structural, and rotation conditions are met. Estimate each design's success probability and expected total cost. Keep designs meeting the **85% baseline reliability threshold** (96 in the sample run).
8. **Nominal selection.** Choose the minimum expected-cost baseline design (sample: #2100), using documented tie-breaking rules.
9. **Verification and tail-risk assessment.** For the nominal design, rank its 5,000 scenario maximum-utilization values and select the scenario at rank 4,400 (sample: #4514). Check every baseline design under that same verification scenario. Also require at least 70% conditional success when passive resistance is unavailable and **CVaR95 <= 1.05** for construction sliding utilization, computed as the average of the 250 largest values among 5,000.
10. **Robust-feasible set.** Retain designs meeting all three supplementary criteria (78 in the sample run).
11. **Minimax regret selection.** Calculate design costs under three **deterministic price worlds** with alternative excavation, concrete, steel, formwork, and backfill multipliers. Compare each eligible design against the lowest cost *within the same world*. Select the design minimizing its maximum cost regret (sample: #2134).

### Interpretation of the risk metrics

- **Maximum utilization, V:** the largest of five aggregated check groups: construction stability, service stability, construction structural demand/capacity, service structural demand/capacity, and service rotation. A selected verification case passes if `V <= 1`.
- **Baseline success probability:** proportion of 5,000 scenarios in which all required engineering checks pass; the minimum accepted probability is 0.85.
- **Passive-out reliability:** conditional success probability among the same shared Monte Carlo scenarios whose passive-availability indicator equals zero; the minimum accepted probability is 0.70.
- **Construction sliding CVaR95:** the average of the **largest** 250 of 5,000 sliding utilization ratios. Its upper limit is 1.05. This measures tail risk in *utilization*, not the average of the smallest safety factors.
- **Minimax regret:** an economic decision rule over three prescribed deterministic cost scenarios, not another 5,000-sample Monte Carlo simulation.

## Project structure

```text
.
|-- main.cpp
|-- model.hpp
|-- engineering.hpp
|-- optimization.hpp
|-- retaining_wall_inputs_v2.csv
|-- retaining_wall_settings_v2.csv
|-- README.md
|-- CMakeLists.txt
|-- build_windows.bat
|-- .gitignore
|-- .gitattributes
|-- scripts/verify_outputs.py
|-- .github/workflows/ci.yml
`-- sample_results/
    |-- retaining_wall_optimization_cpp.json
    `-- all_design_results_cpp.csv
```

### Source organization

- `model.hpp`: CSV parser, data structures, reinforced-concrete detailing, design enumeration, probability sampling.
- `engineering.hpp`: hydraulic finite differences, contact/stability/structural checks, quantities, construction costs.
- `optimization.hpp`: hydraulic caching, Monte Carlo evaluation, reliability screening, verification, CVaR95, robust selection, minimax regret, output serialization.
- `main.cpp`: command-line entry point, orchestration, and result reporting.

## Testing and reproducibility

- The complete simulation compiles with a C++17 compiler and runs without third-party dependencies.
- A clean GCC build was executed for all **12,960,000** design-scenario evaluations. Its results matched a prior MSVC run in all 10 JSON fields and all 2,592 CSV records (ignoring CSV line-ending differences). The updated sampling implementation still needs to be compiled and rerun directly with MSVC for a final independent confirmation.
- GitHub Actions builds and runs the full simulation on both Linux (GCC) and Windows (MSVC) and checks the 10 JSON fields and all 2,592 CSV rows against the included reference outputs.
- The fixed seed and explicitly defined sampling transforms substantially improve cross-toolchain repeatability. The C++ engine `std::mt19937_64` differs from NumPy PCG64. Uniform doubles are generated directly from the upper 53 raw engine bits, and normal values use an explicitly implemented Marsaglia polar transform. This removes the dependency on implementation-specific C++ standard distribution algorithms. Small differences in transcendental floating-point math can still occur across toolchains.

### Engineering limitations

The model intentionally uses idealized earth pressure assumptions (including horizontal backfill and simplified wall/soil friction), a simplified one-dimensional transient hydraulic process, and prescribed contact, structural and cost relationships. It is **not** calibrated against field tests or monitoring data, and independent model verification, convergence studies, parameter-correlation assessment, and project-specific code checks would be required before real-world use. Seismic loading, global slope stability, drainage performance, construction sequencing, and site-dependent failure modes may require additional models.

No private reference document is included in this public-source package.
