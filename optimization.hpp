#pragma once
#include "engineering.hpp"
#include <cfenv>
#include <iostream>
#include <functional>
#include <numeric>
#include <thread>
#include <tuple>

namespace rw {

/**
 * Summarize Monte Carlo reliability, expected cost, verification utilization, tail risk, and minimax regret for a design.
 */
struct DesignSummary {
    int index=0, successes=0, passiveOutSuccesses=0;
    bool baseline=false,robust=false,verificationPass=false;
    double probability=0,passiveOutProbability=0,meanCost=INF,selectionCost=INF;
    double maxVolume=0,steelKg=0,verificationV=INF,cvar95=INF,maxRegret=INF;
    std::array<double,3> costs{};
};

/**
 * Store all design summaries together with the reported selection and hydraulic diagnostic indices.
 */
struct SimulationRun {
    int constructionHydraulicCaseIndex=0;
    double serviceHydraulicCaseWaterMoment=std::numeric_limits<double>::quiet_NaN();
    int baselineCount=0,nominalIndex=0,verificationCaseIndex=0;
    int robustCount=0,finalIndex=0;
    int passiveOutCount=0;
    std::vector<DesignSummary> summaries;
};

/**
 * Round to a specified number of decimal places using ties-to-even; avoid rounding intermediate engineering results.
 */
inline double roundHalfEven(double x,int digits) {
    const double factor=std::pow(10.0,digits);
    return std::nearbyint(x*factor)/factor;
}

/**
 * Cache hydraulic solutions by base width, shear-key depth, and scenario to avoid repeated diffusion solves.
 */
inline std::vector<Hydraulic> makeHydraulicCache(const Config& c,
                                                  const std::vector<Scenario>& scenarios) {
    const auto& bs=c.design.at("base_width");
    const auto& ks=c.design.at("shear_key_depth");
    const std::size_t n=scenarios.size();
    std::vector<Hydraulic> cache(bs.size()*ks.size()*n);
    for(std::size_t b=0;b<bs.size();++b)
        for(std::size_t k=0;k<ks.size();++k) {
            Design d; d.B=bs[b];d.dk=ks[k];
            for(std::size_t i=0;i<n;++i)
                cache[(b*ks.size()+k)*n+i]=hydraulicStages(d,scenarios[i],c);
        }
    return cache;
}

/**
 * Retrieve a cached hydraulic result for the current design and scenario without copying it.
 */
inline const Hydraulic& getCachedHydraulics(const std::vector<Hydraulic>& cache,
                                              const Design& d,const Config& c,
                                              std::size_t simulationIndex,std::size_t n) {
    const std::size_t keyWidth=c.design.at("shear_key_depth").size();
    return cache[(static_cast<std::size_t>(d.widthIndex)*keyWidth+
                  static_cast<std::size_t>(d.keyIndex))*n+simulationIndex];
}

/**
 * Sort value/scenario-index pairs by ascending value and then scenario index, returning the 1-based scenario index at the specified rank.
 */
inline int orderStatisticIndex(std::vector<std::pair<double,int>>& values,int rank) {
    if(rank<1 || rank>static_cast<int>(values.size()))
        throw std::runtime_error("Invalid order statistic rank");
    std::sort(values.begin(),values.end(),[](const auto& a,const auto& b){
        if(a.first<b.first)return true;
        if(a.first>b.first)return false;
        return a.second<b.second;
    });
    return values[static_cast<std::size_t>(rank-1)].second;
}

/**
 * Compute the prescribed construction and service hydraulic diagnostics using specified designs and order statistics.
 */
inline std::pair<int,double> hydraulicDiagnostics(const Config& c,
                  const std::vector<Scenario>& scenarios,
                  const std::vector<Design>& designs,
                  const std::vector<Hydraulic>& cache) {
    const std::size_t n=scenarios.size();
    std::vector<std::pair<double,int>> pairs;
    const Design& d146=designs.at(145);
    for(std::size_t i=0;i<n;++i) {
        const auto& h=getCachedHydraulics(cache,d146,c,i,n);
        pairs.emplace_back(h.mwC+d146.tb*h.pwC,static_cast<int>(i+1));
    }
    const int caseC=orderStatisticIndex(pairs,4550);
    const Design& d377=designs.at(376);
    pairs.clear();
    for(std::size_t i=0;i<n;++i) if(scenarios[i].serviceHead>=0.35) {
        const auto& h=getCachedHydraulics(cache,d377,c,i,n);
        pairs.emplace_back(h.mwS+d377.tb*h.pwS,static_cast<int>(i+1));
    }
    if(pairs.empty())return {caseC,std::numeric_limits<double>::quiet_NaN()};
    const int rank=static_cast<int>(std::ceil(0.90*pairs.size()));
    const int caseS=orderStatisticIndex(pairs,rank);
    const auto& h=getCachedHydraulics(cache,d377,c,static_cast<std::size_t>(caseS-1),n);
    return {caseC,h.mwS+d377.tb*h.pwS};
}

/**
 * Choose the nominal baseline design by rounded expected cost, then probability, concrete volume, steel mass, and index.
 */
inline bool isBetterNominal(const DesignSummary& a,const DesignSummary& b) {
    if(a.selectionCost!=b.selectionCost)return a.selectionCost<b.selectionCost;
    if(a.probability!=b.probability)return a.probability>b.probability;
    if(a.maxVolume!=b.maxVolume)return a.maxVolume<b.maxVolume;
    if(a.steelKg!=b.steelKg)return a.steelKg<b.steelKg;
    return a.index<b.index;
}

/**
 * Compute the average of the largest requested number of values, using a partial selection algorithm; used for CVaR95.
 */
inline double cvarLargest(std::vector<double>& u,std::size_t count) {
    if(count==0 || count>u.size())return INF;
    if(count<u.size())
        std::nth_element(u.begin(),u.begin()+static_cast<std::ptrdiff_t>(count),u.end(),
                          std::greater<double>());
    const double sum=std::accumulate(u.begin(),u.begin()+static_cast<std::ptrdiff_t>(count),0.0);
    return sum/static_cast<double>(count);
}

/**
 * Evaluate 2,592 designs over 5,000 common scenarios; select baseline, nominal, verification and robust sets, then minimize cost regret.
 */
inline SimulationRun calculateAll(const Config& c,
                                   const std::vector<Design>& designs,
                                   const std::vector<Scenario>& scenarios,
                                   bool showProgress=false) {
    const std::size_t n=scenarios.size();
    SimulationRun run;
    if(n!=5000)throw std::runtime_error("This configuration requires exactly 5000 scenarios");
    auto cache=makeHydraulicCache(c,scenarios);
    const auto diag=hydraulicDiagnostics(c,scenarios,designs,cache);
    run.constructionHydraulicCaseIndex=diag.first;
    run.serviceHydraulicCaseWaterMoment=diag.second;
    for(const auto& s:scenarios)if(s.passive==0)++run.passiveOutCount;
    run.summaries.resize(designs.size());
    const int successLimit=static_cast<int>(std::ceil(
        c.getSettingDouble("constraint","target_success_probability")*n));
    // Assign disjoint design indices to each worker thread.
    // Parallel tasks share read-only scenario and cache data.
    // Each worker writes only its own design summary slots.
    const unsigned cores=std::max(1u,std::min(4u,std::thread::hardware_concurrency()));
    std::vector<std::thread> workers;
    for(unsigned worker=0;worker<cores;++worker) {
        workers.emplace_back([&,worker]() {
            for(std::size_t j=worker;j<designs.size();j+=cores) {
                const Design& d=designs[j];
                DesignSummary ds;
                ds.index=d.index;
                if(!d.valid) {run.summaries[j]=ds;continue;}
                const Quantities qty=computeQuantities(d,c);
                ds.maxVolume=qty.concrete;
                ds.steelKg=qty.steelKg;
                ds.costs=deterministicCosts(qty,c);
                double costSum=0;
                for(std::size_t i=0;i<n;++i) {
                    const auto& h=getCachedHydraulics(cache,d,c,i,n);
                    Outcome r=evaluate(d,scenarios[i],c,h,qty);
                    if(r.success) {
                        ++ds.successes;
                        if(scenarios[i].passive==0)++ds.passiveOutSuccesses;
                    }
                    costSum+=r.cost;
                }
                ds.probability=ds.successes/static_cast<double>(n);
                ds.passiveOutProbability=run.passiveOutCount?
                    ds.passiveOutSuccesses/static_cast<double>(run.passiveOutCount):
                    std::numeric_limits<double>::quiet_NaN();
                ds.meanCost=costSum/static_cast<double>(n);
                ds.selectionCost=roundHalfEven(ds.meanCost,2);
                ds.baseline=ds.successes>=successLimit;
                run.summaries[j]=ds;
            }
        });
    }
    for(auto& t:workers)t.join(); // Wait for all workers.
    for(const auto& ds:run.summaries)if(ds.baseline)++run.baselineCount;
    if(showProgress)std::cerr<<"  Design evaluations complete ("
                             <<cores<<" parallel C++ workers).\n";
    // Select nominal design; fall back to a valid minimum-cost design if baseline is empty.
    int nominal=-1;
    for(std::size_t j=0;j<designs.size();++j) {
        const auto& ds=run.summaries[j];
        if(!designs[j].valid || (run.baselineCount>0 && !ds.baseline))continue;
        if(nominal<0 || isBetterNominal(ds,run.summaries[static_cast<std::size_t>(nominal)]))
            nominal=static_cast<int>(j);
    }
    if(nominal<0)throw std::runtime_error("No valid design was generated");
    run.nominalIndex=designs[static_cast<std::size_t>(nominal)].index;
    if(run.baselineCount==0) {run.finalIndex=run.nominalIndex;return run;}
    const Design& nd=designs[static_cast<std::size_t>(nominal)];
    const Quantities nQty=computeQuantities(nd,c);
    std::vector<std::pair<double,int>> scored;
    scored.reserve(n);
    for(std::size_t i=0;i<n;++i) {
        const auto& h=getCachedHydraulics(cache,nd,c,i,n);
        const auto result=evaluate(nd,scenarios[i],c,h,nQty);
        scored.emplace_back(result.V,static_cast<int>(i+1));
    }
    run.verificationCaseIndex=orderStatisticIndex(scored,4400);
    const std::size_t verificationI=static_cast<std::size_t>(run.verificationCaseIndex-1);
    // Apply additional risk constraints to baseline-feasible designs only.
    for(std::size_t j=0;j<designs.size();++j) {
        auto& ds=run.summaries[j];
        if(!ds.baseline)continue;
        const Design& d=designs[j];
        const auto qty=computeQuantities(d,c);
        const auto& verificationH=getCachedHydraulics(cache,d,c,verificationI,n);
        ds.verificationV=evaluate(d,scenarios[verificationI],c,verificationH,qty).V;
        ds.verificationPass=ds.verificationV<=1;
        std::vector<double> sliding;
        sliding.reserve(n);
        for(std::size_t i=0;i<n;++i) {
            const auto& h=getCachedHydraulics(cache,d,c,i,n);
            const auto response=evaluateStage(d,scenarios[i],c,h,false);
            sliding.push_back(response.slidingUtilization);
        }
        ds.cvar95=cvarLargest(sliding,250);
        ds.robust=ds.verificationPass &&
                  ds.passiveOutProbability>=0.70 && ds.cvar95<=1.05;
        if(ds.robust)++run.robustCount;
    }
    if(run.robustCount==0) { run.finalIndex=run.nominalIndex;return run; }
    // Determine the minimum cost within each deterministic price world.
    std::array<double,3> bestCosts{{INF,INF,INF}};
    for(const auto& d:run.summaries) if(d.robust)
        for(int w=0;w<3;++w)bestCosts[w]=std::min(bestCosts[w],d.costs[w]);
    int winner=-1;
    for(std::size_t j=0;j<run.summaries.size();++j) {
        auto& d=run.summaries[j];
        if(!d.robust)continue;
        d.maxRegret=std::max({d.costs[0]-bestCosts[0],
                              d.costs[1]-bestCosts[1],
                              d.costs[2]-bestCosts[2]});
        if(winner<0)winner=static_cast<int>(j);
        else {
            const auto& b=run.summaries[static_cast<std::size_t>(winner)];
            const bool better=d.maxRegret<b.maxRegret ||
                (d.maxRegret==b.maxRegret &&
                  (d.meanCost<b.meanCost ||
                  (d.meanCost==b.meanCost &&
                  (d.probability>b.probability ||
                  (d.probability==b.probability && d.index<b.index)))));
            if(better)winner=static_cast<int>(j);
        }
    }
    run.finalIndex=run.summaries[static_cast<std::size_t>(winner)].index;
    return run;
}

/**
 * Retrieve the stored summary for a 1-based design index.
 */
inline const DesignSummary& findSummary(const SimulationRun& run,int index) {
    return run.summaries.at(static_cast<std::size_t>(index-1));
}

/**
 * Write the ten reported metrics as JSON.
 */
inline void writeResultsJson(const std::string& output, const SimulationRun& run) {
    std::ofstream f(output);
    if(!f)throw std::runtime_error("Cannot write JSON file: "+output);
    const auto& selected=findSummary(run,run.finalIndex);
    f<<"{\n";
    f<<"  \"construction_hydraulic_case_index\": "<<run.constructionHydraulicCaseIndex<<",\n";
    f<<"  \"service_hydraulic_case_water_moment\": ";
    if(std::isfinite(run.serviceHydraulicCaseWaterMoment))
        f<<std::fixed<<std::setprecision(3)<<roundHalfEven(run.serviceHydraulicCaseWaterMoment,3);
    else f<<"null";
    f<<",\n  \"baseline_feasible_design_count\": "<<run.baselineCount<<",\n";
    f<<"  \"nominal_selected_design_index\": "<<run.nominalIndex<<",\n";
    f<<"  \"derived_verification_case_index\": "<<run.verificationCaseIndex<<",\n";
    f<<"  \"robust_feasible_design_count\": "<<run.robustCount<<",\n";
    f<<"  \"final_selected_design_index\": "<<run.finalIndex<<",\n";
    f<<"  \"selected_conditional_passive_out_success_probability\": ";
    if(std::isfinite(selected.passiveOutProbability))
        f<<std::fixed<<std::setprecision(4)<<roundHalfEven(selected.passiveOutProbability,4);
    else f<<"null";
    f<<",\n  \"selected_maximum_cost_regret\": ";
    if(std::isfinite(selected.maxRegret))
        f<<std::fixed<<std::setprecision(2)<<roundHalfEven(selected.maxRegret,2);
    else f<<"null";
    f<<",\n  \"expected_total_cost\": "<<std::fixed<<std::setprecision(2)
     <<roundHalfEven(selected.meanCost,2)<<"\n}\n";
}

/**
 * Write success counts, reliability, costs, verification, risk, and regret for every design to CSV.
 */
inline void writeResultsCsv(const std::string& output,const SimulationRun& run) {
    std::ofstream f(output);
    if(!f)throw std::runtime_error("Cannot write summary CSV file: "+output);
    f<<"design_index,success_count,success_probability,mean_cost,baseline,"
        "passive_out_probability,verification_utilization,cvar95,robust,max_regret\n";
    for(const auto& d:run.summaries) {
        f<<d.index<<","<<d.successes<<","<<std::setprecision(12)<<d.probability
         <<","<<d.meanCost<<","<<int(d.baseline)<<",";
        if(std::isfinite(d.passiveOutProbability))f<<d.passiveOutProbability;
        f<<",";
        if(std::isfinite(d.verificationV))f<<d.verificationV;
        f<<",";
        if(std::isfinite(d.cvar95))f<<d.cvar95;
        f<<","<<int(d.robust)<<",";
        if(std::isfinite(d.maxRegret))f<<d.maxRegret;
        f<<"\n";
    }
}

} // namespace rw