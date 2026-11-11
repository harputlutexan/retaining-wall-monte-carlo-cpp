#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

/**
 * Retaining wall Monte Carlo model: CSV parsing, configuration, scenario sampling, and design enumeration.
 */

namespace rw {

constexpr double PI = 3.1415926535897932384626433832795;
constexpr double INF = std::numeric_limits<double>::infinity();

/**
 * Remove leading and trailing ASCII whitespace from a string.
 */
inline std::string trim(const std::string& s) {
    const auto i = s.find_first_not_of(" \t\r\n");
    if (i == std::string::npos) return "";
    return s.substr(i, s.find_last_not_of(" \t\r\n") - i + 1);
}

/**
 * Split a string on a single-character delimiter, retaining trailing empty fields.
 */
inline std::vector<std::string> split(const std::string& s, char delim) {
    std::vector<std::string> out;
    std::stringstream stream(s);
    std::string token;
    while (std::getline(stream, token, delim)) out.push_back(trim(token));
    if (!s.empty() && s.back() == delim) out.emplace_back("");
    return out;
}

/**
 * Parse a CSV record, supporting quoted fields, embedded commas, and escaped quotes.
 */
inline std::vector<std::string> parseCsvRow(const std::string& s) {
    std::vector<std::string> fields;
    std::string field;
    bool insideQuotes = false;
    for (std::size_t i = 0; i < s.size(); ++i) {
        const char ch = s[i];
        if (ch == '"') {
            if (insideQuotes && i + 1 < s.size() && s[i + 1] == '"') {
                field.push_back('"'); ++i;
            } else insideQuotes = !insideQuotes;
        } else if (ch == ',' && !insideQuotes) {
            fields.push_back(trim(field)); field.clear();
        } else field.push_back(ch);
    }
    if (insideQuotes) throw std::runtime_error("CSV: unterminated quoted field");
    fields.push_back(trim(field));
    return fields;
}

/**
 * Read a CSV file into rows of strings. The first row contains column headers.
 */
inline std::vector<std::vector<std::string>> readCsv(const std::string& filename) {
    std::ifstream input(filename);
    if (!input) throw std::runtime_error("Cannot open CSV file: " + filename);
    std::vector<std::vector<std::string>> rows;
    std::string line;
    while (std::getline(input, line)) {
        if (trim(line).empty()) continue;
        rows.push_back(parseCsvRow(line));
    }
    if (rows.empty()) throw std::runtime_error("CSV file is empty: " + filename);
    const std::size_t colCount = rows.front().size();
    for (std::size_t i = 1; i < rows.size(); ++i)
        if (rows[i].size() != colCount)
            throw std::runtime_error("CSV column count mismatch: " + filename + " row " + std::to_string(i + 1));
    return rows;
}

/**
 * Parse a required numeric CSV field and reject missing or invalid values.
 */
inline double requireDouble(const std::string& s, const std::string& key) {
    if (s.empty()) throw std::runtime_error("Missing numeric value: " + key);
    std::size_t read = 0;
    const double value = std::stod(s, &read);
    if (read != s.size() || !std::isfinite(value))
        throw std::runtime_error("Invalid number: " + key + "=" + s);
    return value;
}

/**
 * Convert a pipe-separated option list such as 4.0|4.6 into numeric values.
 */
inline std::vector<double> numbers(const std::string& s) {
    std::vector<double> v;
    for (const auto& entry : split(s, '|')) v.push_back(requireDouble(entry, "allowed_values"));
    return v;
}

/**
 * Keep related data in named fields rather than unrelated scalar variables.
 */
struct UncertainSpec {
    std::string name, distribution, unit;
    double mean = 0, sd = 0, cv = 0, mode = 0, minimum = 0, maximum = 0;
    std::vector<double> states, probabilities;
};

struct Config {
    std::map<std::string, double> fixed;
    std::map<std::string, std::vector<double>> design;
    std::map<std::string, std::string> settings;
    std::map<std::string, double> numericSettings; // Parse numeric CSV settings once.
    std::vector<UncertainSpec> uncertain;
    std::vector<double> barDiameters, barSpacings;

    /**
     * Retrieve a required fixed engineering parameter by name.
     */
    double getFixed(const std::string& name) const { return fixed.at(name); }

    /**
     * Retrieve a configuration setting by its category and parameter name.
     */
    std::string getSetting(const std::string& cat, const std::string& name) const {
        return settings.at(cat + "|" + name);
    }

    /**
     * Retrieve a numeric analysis setting by category and parameter name.
     */
    double getSettingDouble(const std::string& cat, const std::string& name) const {
        return numericSettings.at(cat + "|" + name);
    }
};

/**
 * Load the two CSV files into fixed inputs, design alternatives, uncertain distributions, detailing choices, and settings.
 */
inline Config loadConfig(const std::string& inputPath, const std::string& settingsPath) {
    Config c;
    const auto rows = readCsv(inputPath);
    const std::vector<std::string> expected = {
        "input_group","parameter","fixed_value","allowed_values","distribution",
        "mean","standard_deviation","coefficient_of_variation","mode","states",
        "probabilities","minimum","maximum","unit"
    };
    if (rows.front() != expected) throw std::runtime_error("Unexpected column order in inputs CSV");
    for (std::size_t i = 1; i < rows.size(); ++i) {
        const auto& r = rows[i];
        const auto& group = r[0];
        const auto& key = r[1];
        if (group == "fixed") c.fixed[key] = requireDouble(r[2], key);
        else if (group == "design") c.design[key] = numbers(r[3]);
        else if (group == "detailing") {
            if (key == "main_bar_diameter_options") c.barDiameters = numbers(r[3]);
            else if (key == "main_bar_spacing_options") c.barSpacings = numbers(r[3]);
        } else if (group == "uncertain") {
            UncertainSpec u;
            u.name = key; u.distribution = r[4]; u.unit = r[13];
            if (!r[5].empty()) u.mean = requireDouble(r[5], key + "/mean");
            if (!r[6].empty()) u.sd = requireDouble(r[6], key + "/sd");
            if (!r[7].empty()) u.cv = requireDouble(r[7], key + "/cv");
            if (!r[8].empty()) u.mode = requireDouble(r[8], key + "/mode");
            if (!r[9].empty()) u.states = numbers(r[9]);
            if (!r[10].empty()) u.probabilities = numbers(r[10]);
            if (!r[11].empty()) u.minimum = requireDouble(r[11], key + "/min");
            if (!r[12].empty()) u.maximum = requireDouble(r[12], key + "/max");
            c.uncertain.push_back(std::move(u));
        } else throw std::runtime_error("Unknown input_group: " + group);
    }
    const auto srows = readCsv(settingsPath);
    if (srows.front() != std::vector<std::string>({"category","parameter","value","unit"}))
        throw std::runtime_error("Unexpected settings CSV header");
    for (std::size_t i = 1; i < srows.size(); ++i) {
        const std::string key=srows[i][0] + "|" + srows[i][1];
        c.settings[key]=srows[i][2];
        // Parse numeric values and leave descriptive formula settings as text.
        try {
            std::size_t consumed=0;
            const double value=std::stod(srows[i][2],&consumed);
            if(consumed==srows[i][2].size())c.numericSettings[key]=value;
        } catch(const std::exception&) {}
    }
    if (c.uncertain.size() != 15 || c.fixed.size() != 21 || c.design.size() != 7)
        throw std::runtime_error("Expected 15 uncertain, 21 fixed, and 7 design entries in CSV");
    return c;
}

/**
 * One Monte Carlo scenario contains 15 sampled uncertain engineering and cost inputs.
 */
struct Scenario {
    double phi = 0, gamma = 0, surcharge = 0, qult = 0;
    double fc = 0, fy = 0, qConstruction = 0, etaDelta = 0;
    int passive = 0;
    double etaE = 0, etaConcreteCost = 0, etaSteelCost = 0;
    double initialHead = 0, diffusivity = 0, serviceHead = 0;
};

/**
 * Generate stochastic inputs using a fixed mapping from mt19937_64 raw bits.
 * This avoids implementation-dependent standard-library distribution algorithms.
 * Floating-point math functions may still differ slightly across toolchains.
 */
class Sampler {
    std::mt19937_64 engine_;
public:
    /**
     * Initialize the standard-defined 64-bit Mersenne Twister engine.
     */
    explicit Sampler(std::uint64_t seed) : engine_(seed) {}

    /**
     * Map the top 53 raw engine bits to a uniform double in [0,1).
     * The factor 2^-53 is exactly representable in IEEE-754 binary64.
     */
    double uniform01() {
        constexpr double twoToMinus53 = 1.0 / 9007199254740992.0;
        return static_cast<double>(engine_() >> 11) * twoToMinus53;
    }

    /**
     * Draw one standard normal variate with the Marsaglia polar method.
     * No cached second variate: each call consumes whole uniform pairs,
     * making both the draw order and rejection behavior explicit.
     */
    double normal(double mean, double sd) {
        for (;;) {
            const double x = 2.0 * uniform01() - 1.0;
            const double y = 2.0 * uniform01() - 1.0;
            const double r2 = x*x + y*y;
            if (r2 <= 0.0 || r2 >= 1.0) continue;
            const double z = x * std::sqrt(-2.0 * std::log(r2) / r2);
            return mean + sd * z;
        }
    }

    /**
     * Sample a truncated normal distribution by rejection sampling until a value lies within the specified bounds.
     */
    double truncatedNormal(const UncertainSpec& u) {
        for (int attempts = 0; attempts < 1000000; ++attempts) {
            const double x = normal(u.mean, u.sd);
            if (u.minimum <= x && x <= u.maximum) {
                return x;
            }
        }
        throw std::runtime_error("Truncated normal rejection sampling exceeded maximum attempts: " + u.name);
    }

    /**
     * Convert arithmetic mean and coefficient of variation to log-space parameters; sample a lognormal variate and reject out-of-bounds values.
     */
    double truncatedLognormal(const UncertainSpec& u) {
        const double sigma = std::sqrt(std::log(1.0 + u.cv * u.cv));
        const double mu = std::log(u.mean) - 0.5 * sigma * sigma;
        for (int attempts = 0; attempts < 1000000; ++attempts) {
            const double x = std::exp(normal(mu, sigma));
            if (u.minimum <= x && x <= u.maximum) {
                return x;
            }
        }
        throw std::runtime_error("Truncated lognormal rejection sampling exceeded maximum attempts: " + u.name);
    }

    /**
     * Sample a triangular distribution using its inverse cumulative distribution function and a uniform variate.
     */
    double triangular(const UncertainSpec& u) {
        const double a = u.minimum, b = u.maximum, c = u.mode;
        const double r = uniform01();
        const double p = (c - a) / (b - a);
        if (r < p) return a + std::sqrt(r * (b - a) * (c - a));
        return b - std::sqrt((1.0 - r) * (b - a) * (b - c));
    }

    /**
     * Sample a categorical distribution using cumulative state probabilities and one uniform variate.
     */
    double categorical(const UncertainSpec& u) {
        const double r = uniform01();
        double cumulative = 0;
        for (std::size_t i = 0; i < u.states.size(); ++i) {
            cumulative += u.probabilities.at(i);
            if (r < cumulative) return u.states[i];
        }
        if (std::abs(cumulative - 1.0) > 1e-12)
            throw std::runtime_error("Categorical probabilities do not sum to one: " + u.name);
        return u.states.back(); // Handle potential floating-point rounding at the upper boundary.
    }

    /**
     * Dispatch to the sampling method specified for an uncertain input.
     */
    double sampleValue(const UncertainSpec& u) {
        if (u.distribution == "truncated_normal") return truncatedNormal(u);
        if (u.distribution == "truncated_lognormal") return truncatedLognormal(u);
        if (u.distribution == "triangular") return triangular(u);
        if (u.distribution == "categorical") return categorical(u);
        throw std::runtime_error("Unknown distribution: " + u.distribution);
    }

    /**
     * Sample the 15 uncertain inputs in CSV row order and assign their physical values to a scenario.
     */
    Scenario sampleScenario(const Config& c) {
        Scenario s;
        std::array<double, 15> x{};
        for (std::size_t i = 0; i < 15; ++i) {
            const auto& spec=c.uncertain.at(i);
            x[i]=sampleValue(spec);
        }
        s.phi=x[0]; s.gamma=x[1]; s.surcharge=x[2]; s.qult=x[3];
        s.fc=x[4]; s.fy=x[5]; s.qConstruction=x[6]; s.etaDelta=x[7];
        s.passive=static_cast<int>(x[8]); s.etaE=x[9];
        s.etaConcreteCost=x[10]; s.etaSteelCost=x[11];
        s.initialHead=x[12]; s.diffusivity=x[13]; s.serviceHead=x[14];
        return s;
    }
};

/**
 * Generate the 5,000 scenarios once. Reusing this shared sample set across designs applies the common random numbers method.
 */
inline std::vector<Scenario> makeScenarios(const Config& c) {
    const auto seed = static_cast<std::uint64_t>(c.getSettingDouble("sampling", "random_seed"));
    const auto n = static_cast<std::size_t>(c.getSettingDouble("sampling", "iterations"));
    Sampler sampler(seed);
    std::vector<Scenario> all;
    all.reserve(n);
    for (std::size_t i = 0; i < n; ++i) all.push_back(sampler.sampleScenario(c));
    return all;
}

struct Rebar { double target = 0, area = 0, diameter = 0, spacing = 0; };

struct Design {
    int index = 0;
    int widthIndex = 0, keyIndex = 0;  // Hydraulic cache lookup keys.
    double B = 0, rt = 0, ts = 0, tb = 0, dk = 0;
    Rebar stem, base;
    double toe = 0, heel = 0, xBack = 0, xStem = 0, stemArea = 0;
    double rhoStem = 0, rhoBase = 0;
    double dStemMM = 0, dBaseMM = 0;
    bool valid = false;
};

/**
 * Choose a reinforcement diameter and spacing meeting the target steel area with minimum provided area; apply deterministic tie-breaking.
 */
inline Rebar detailRebar(double target, const Config& c) {
    Rebar best;
    best.target = target;
    best.area = INF;
    for (const double diameter : c.barDiameters)
        for (const double spacing : c.barSpacings) {
            const double area = (PI * diameter * diameter / 4.0) * (1000.0 / spacing);
            if (area + 1e-10 < target) continue;
            const bool better = area < best.area - 1e-10 ||
                (std::abs(area - best.area) <= 1e-10 &&
                 (diameter < best.diameter ||
                  (diameter == best.diameter && spacing > best.spacing)));
            if (better) best.area=area, best.diameter=diameter, best.spacing=spacing;
        }
    if (!std::isfinite(best.area)) throw std::runtime_error("No reinforcement bar diameter/spacing pair meets target");
    return best;
}

/**
 * Complete geometric properties and provided reinforcement, then check geometric and reinforcement feasibility.
 */
inline void finishDesign(Design& d, const Config& c) {
    const double H = c.getFixed("wall_height");
    const double tt = c.getFixed("stem_top_thickness");
    d.stem=detailRebar(d.stem.target,c);
    d.base=detailRebar(d.base.target,c);
    d.toe=d.rt*d.B; d.heel=d.B-d.toe-d.ts;
    d.xBack=d.toe+d.ts;
    d.stemArea=H*(d.ts+tt)/2.0;
    const double centroidOffset=(d.ts*d.ts+d.ts*tt+tt*tt)/(3.0*(d.ts+tt));
    d.xStem=d.xBack-centroidOffset;
    d.rhoStem=d.stem.area/(1000.0*(1000.0*d.ts));
    d.rhoBase=d.base.area/(1000.0*(1000.0*d.tb));
    d.dStemMM=1000.0*(d.ts-c.getFixed("stem_concrete_cover"))-d.stem.diameter/2.0;
    d.dBaseMM=1000.0*(d.tb-c.getFixed("base_concrete_cover"))-d.base.diameter/2.0;
    const double minR=c.getFixed("minimum_reinforcement_ratio");
    const double maxR=c.getFixed("maximum_reinforcement_ratio");
    d.valid=d.heel>0 && d.dStemMM>0 && d.dBaseMM>0 &&
            d.rhoStem>=minR && d.rhoStem<=maxR &&
            d.rhoBase>=minR && d.rhoBase<=maxR;
}

/**
 * Enumerate the Cartesian product of all seven design variables, assigning stable 1-based design identifiers.
 */
inline std::vector<Design> makeDesigns(const Config& c) {
    auto& B=c.design.at("base_width");
    auto& rt=c.design.at("toe_length_ratio");
    auto& ts=c.design.at("stem_base_thickness");
    auto& tb=c.design.at("base_slab_thickness");
    auto& dk=c.design.at("shear_key_depth");
    auto& As=c.design.at("stem_main_steel_area");
    auto& Ab=c.design.at("base_main_steel_area");
    std::vector<Design> designs;
    int index=0;
    for (std::size_t i=0;i<B.size();++i)
        for (double r:rt) for(double stem:ts) for(double base:tb)
            for(std::size_t k=0;k<dk.size();++k)
                for(double stemTarget:As) for(double baseTarget:Ab) {
                    Design d;
                    d.index=++index; d.widthIndex=static_cast<int>(i); d.keyIndex=static_cast<int>(k);
                    d.B=B[i];d.rt=r;d.ts=stem;d.tb=base;d.dk=dk[k];
                    d.stem.target=stemTarget;d.base.target=baseTarget;
                    finishDesign(d,c);
                    designs.push_back(d);
                }
    return designs;
}

} // namespace rw