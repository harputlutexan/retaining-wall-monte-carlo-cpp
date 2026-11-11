#pragma once
#include "model.hpp"
#include <numeric>

namespace rw {

/**
 * Store the 11-node hydraulic pressure profiles for construction and service, along with integrated water forces and moments.
 */
struct Hydraulic {
    std::array<double,11> initialC{}, firstC{}, finalC{}, initialS{}, finalS{};
    double pwC=0,mwC=0,pwS=0,mwS=0;
};

/**
 * Solve the 9x9 tridiagonal linear system by the Thomas algorithm, with a pivot check and linear time complexity.
 */
inline std::array<double,9> solveTridiagonal(double low, double diag, double high,
                                             std::array<double,9> rhs) {
    std::array<double,9> modifiedUpper{}, solution{};
    double pivot=diag;
    if (std::abs(pivot)<1e-14) throw std::runtime_error("Thomas algorithm encountered a zero pivot");
    modifiedUpper[0]=high/pivot;
    rhs[0]/=pivot;
    for (int j=1;j<9;++j) {
        pivot=diag-low*modifiedUpper[j-1];
        if (std::abs(pivot)<1e-14) throw std::runtime_error("Thomas algorithm encountered a zero pivot");
        modifiedUpper[j]=(j==8?0.0:high/pivot);
        rhs[j]=(rhs[j]-low*rhs[j-1])/pivot;
    }
    solution[8]=rhs[8];
    for(int j=7;j>=0;--j) solution[j]=rhs[j]-modifiedUpper[j]*solution[j+1];
    return solution;
}

/**
 * Advance the one-dimensional diffusion equation by one Crank-Nicolson step with fixed Dirichlet boundary conditions.
 */
inline std::array<double,11> oneCrankNicolsonStep(const std::array<double,11>& old,
                                                   double r) {
    const double off=-0.5*r, diag=1+r;
    std::array<double,9> rhs{};
    for(int j=1;j<=9;++j)
        rhs[j-1]=(1-r)*old[j]+0.5*r*(old[j-1]+old[j+1]);
    rhs[0] += 0.5*r*old[0];  // Add lower Dirichlet boundary contribution.
    rhs[8] += 0.5*r*old[10]; // Add upper Dirichlet boundary contribution.
    const auto inner=solveTridiagonal(off,diag,off,rhs);
    std::array<double,11> next{};
    next[0]=old[0];next[10]=old[10];
    for(int j=1;j<=9;++j) next[j]=inner[j-1];
    return next;
}

/**
 * Integrate a piecewise-linear pressure segment exactly to obtain resultant force and moment.
 */
inline std::pair<double,double> linearProfileIntegral(double z0,double z1,
                                                       double u0,double u1) {
    const double h=z1-z0;
    return {h*(u0+u1)/2.0,
            h*((2*z0+z1)*u0+(z0+2*z1)*u1)/6.0};
}

/**
 * Integrate the 10 linear segments of an 11-node pressure profile to obtain total water force and stem-base moment.
 */
inline std::pair<double,double> integrateWater(const std::array<double,11>& u,
                                                double H) {
    double pw=0,mw=0;
    const double dz=H/10.0;
    for(int j=0;j<10;++j) {
        auto pm=linearProfileIntegral(j*dz,(j+1)*dz,u[j],u[j+1]);
        pw+=pm.first;mw+=pm.second;
    }
    return {pw,mw};
}

/**
 * Integrate water pressure above a specified stem shear section, interpolating where the section cuts a segment.
 */
inline double waterAbove(const std::array<double,11>& u, double H, double section) {
    const double dz=H/10.0;
    double force=0;
    for(int j=0;j<10;++j) {
        const double a=j*dz,b=(j+1)*dz;
        const double left=std::max(a,section);
        if(left>=b)continue;
        const double ul=u[j]+(u[j+1]-u[j])*(left-a)/(b-a);
        force+=linearProfileIntegral(left,b,ul,u[j+1]).first;
    }
    return force;
}

/**
 * Simulate 12 construction and 24 service time steps, passing the final construction pressure profile into the service stage.
 */
inline Hydraulic hydraulicStages(const Design& d,const Scenario& s,const Config& c) {
    Hydraulic h;
    const double H=c.getFixed("wall_height");
    const double gammaW=c.getFixed("water_unit_weight");
    const double dz=H/10.0;
    const double durationC=8.0+1.5*d.B+4.0*d.dk;
    // Fixed stage-step counts: 12 construction steps and 24 service steps over 45 days.
    const int stepsC=12, stepsS=24;
    const double dtC=durationC/stepsC;
    const double dtS=45.0/stepsS;
    const double rC=s.diffusivity*dtC/(dz*dz);
    const double rS=s.diffusivity*dtS/(dz*dz);
    h.initialC[0]=gammaW*H*0.08;
    h.initialC[10]=0;
    for(int j=1;j<=9;++j) h.initialC[j]=gammaW*H*s.initialHead*(1.0-j/10.0);
    auto profile=h.initialC;
    for(int step=0;step<stepsC;++step) {
        profile=oneCrankNicolsonStep(profile,rC);
        if(step==0) h.firstC=profile;
    }
    h.finalC=profile;
    auto pm=integrateWater(h.finalC,H);
    h.pwC=pm.first;h.mwC=pm.second;
    h.initialS=h.finalC;
    h.initialS[0]=gammaW*H*s.serviceHead;
    h.initialS[10]=0;
    profile=h.initialS;
    for(int step=0;step<stepsS;++step)
        profile=oneCrankNicolsonStep(profile,rS);
    h.finalS=profile;
    pm=integrateWater(h.finalS,H);
    h.pwS=pm.first;h.mwS=pm.second;
    return h;
}

/**
 * Represent the supported foundation contact interval and its linear pressure distribution; pressure outside contact is zero.
 */
struct Contact {
    enum class Branch { Full, Toe, Heel, None };
    Branch branch=Branch::None;
    double left=0,right=0,slope=0,intercept=0,qmax=INF,length=0;

    /**
     * Evaluate the soil contact pressure q(x) at a location on the foundation.
     */
    double pressure(double x) const {
        if(x<left || x>right || branch==Branch::None)return 0;
        return slope*x+intercept;
    }
};

/**
 * Choose full, toe-only, heel-only, or no contact based on the resultant location; compute contact pressures from static equilibrium.
 */
inline Contact makeContact(double B,double N,double xR) {
    Contact x;
    if (xR>=B/3 && xR<=2*B/3) {
        x.branch=Contact::Branch::Full; x.left=0;x.right=B;x.length=B;
        x.slope=12*N*(xR-B/2)/(B*B*B);
        x.intercept=N/B-x.slope*B/2;
        x.qmax=std::max(x.pressure(0),x.pressure(B));
    } else if(xR>0 && xR<B/3) {
        x.branch=Contact::Branch::Toe; x.left=0;x.right=3*xR;x.length=x.right;
        x.qmax=2*N/x.length;
        x.intercept=x.qmax;x.slope=-x.qmax/x.length;
    } else if(xR>2*B/3 && xR<B) {
        x.branch=Contact::Branch::Heel;x.length=3*(B-xR);
        x.left=B-x.length;x.right=B;x.qmax=2*N/x.length;
        x.slope=x.qmax/x.length;x.intercept=-x.slope*x.left;
    }
    return x;
}

/**
 * Analytically integrate the net foundation pressure to obtain shear force and bending moment about a specified section.
 */
inline std::pair<double,double> integrateNet(const Contact& contact,
                                              double a,double b,double dead,
                                              double origin) {
    if(b<=a)return {0,0};
    double F=-dead*(b-a);
    double M=-dead*((b*b-a*a)/2-origin*(b-a));
    const double l=std::max(a,contact.left),r=std::min(b,contact.right);
    if(contact.branch!=Contact::Branch::None && r>l) {
        const double A=contact.slope,B=contact.intercept;
        F+=A*(r*r-l*l)/2+B*(r-l);
        M+=A*(r*r*r-l*l*l)/3+(B-A*origin)*(r*r-l*l)/2-B*origin*(r-l);
    }
    return {F,M};
}

/**
 * Hold reinforced-concrete section flexural and shear capacities.
 */
struct Capacity { double moment=0,shear=0,compressionBlock=0; };

/**
 * Compute flexural capacity using a rectangular concrete stress block and concrete shear capacity using the selected provided reinforcement and sampled strengths.
 */
inline Capacity rcCapacity(double As,double d,const Scenario& s,const Config& c) {
    Capacity cap;
    const double b=1000;
    cap.compressionBlock=As*s.fy/(0.85*s.fc*b);
    cap.moment=c.getFixed("flexural_strength_reduction_factor")*As*s.fy*(d-cap.compressionBlock/2)/1e6;
    cap.shear=c.getFixed("shear_strength_reduction_factor")*
        c.getFixed("concrete_shear_coefficient")*std::sqrt(s.fc)*b*d/1e3;
    return cap;
}

/**
 * Collect construction- or service-stage stability and structural demand/capacity results.
 */
struct StageResult {
    bool stabilityPass=false,structuralPass=false;
    double pa=0,pq=0,pw=0,mwStem=0,mwBase=0,normal=0,mResist=0,mOverturn=0;
    double ka=0,kp=0,friction=0,passiveForce=0,fsSliding=0,fsOverturn=0;
    double xResultant=0,qmax=INF,contactRatio=0,stabilityUtilization=INF;
    double slidingUtilization=INF,overturnUtilization=INF,bearingUtilization=INF,contactUtilization=INF;
    Contact contact;
    std::array<double,6> dc{}; // stemM,stemV,toeM,toeV,heelM,heelV
    std::array<double,3> moments{}; // stem,toe,heel kNm/m
    std::array<double,3> shears{};  // stem,toe,heel kN/m
    double dcMax=INF;
};

/**
 * Compute Rankine earth pressures, wall and soil weights, sliding, overturning, bearing and contact, and stem/toe/heel bending and shear.
 */
inline StageResult evaluateStage(const Design& d,const Scenario& s,
                                  const Config& c,const Hydraulic& h,bool service) {
    StageResult r;
    const double H=c.getFixed("wall_height");
    const double gammaC=c.getFixed("concrete_unit_weight");
    const double phiRad=s.phi*PI/180.0;
    const double sinPhi=std::sin(phiRad);
    r.ka=(1-sinPhi)/(1+sinPhi);
    r.kp=(1+sinPhi)/(1-sinPhi);
    const double q=service?s.surcharge:s.qConstruction;
    r.pa=0.5*r.ka*s.gamma*H*H;
    r.pq=r.ka*q*H;
    r.pw=service?h.pwS:h.pwC;
    r.mwStem=service?h.mwS:h.mwC;
    r.mwBase=r.mwStem+d.tb*r.pw;
    const double wBase=gammaC*d.B*d.tb;
    const double wKey=gammaC*c.getFixed("shear_key_width")*d.dk;
    const double wStem=gammaC*d.stemArea;
    const double wHeel=s.gamma*H*d.heel;
    r.normal=wBase+wKey+wStem+wHeel;
    r.mResist=wBase*d.B/2+wKey*(d.toe+d.ts/2)+wStem*d.xStem+
              wHeel*(d.xBack+d.B)/2;
    const double hp=c.getFixed("foundation_embedment_depth")+d.tb+d.dk;
    const double reducedPassive=c.getFixed("passive_resistance_reduction_factor")*
                                0.5*r.kp*s.gamma*hp*hp;
    r.passiveForce=service?s.passive*reducedPassive:0;
    r.friction=r.normal*c.getFixed("interface_friction_factor")*s.etaDelta*
               std::tan(phiRad);
    const double driving=r.pa+r.pq+r.pw;
    r.fsSliding=driving>0?(r.friction+r.passiveForce)/driving:INF;
    r.mOverturn=r.pa*(d.tb+H/3)+r.pq*(d.tb+H/2)+r.mwBase;
    r.fsOverturn=r.mOverturn>0?r.mResist/r.mOverturn:INF;
    r.xResultant=(r.mResist-r.mOverturn)/r.normal;
    r.contact=makeContact(d.B,r.normal,r.xResultant);
    r.qmax=r.contact.qmax;
    r.contactRatio=r.contact.length/d.B;
    // Use the prescribed construction and service contact ratio limits.
    const double minContact=service?0.75:0.60;
    r.slidingUtilization=c.getSettingDouble("constraint","minimum_sliding_factor_of_safety")/r.fsSliding;
    r.overturnUtilization=c.getSettingDouble("constraint","minimum_overturning_factor_of_safety")/r.fsOverturn;
    const double allowable=s.qult/c.getSettingDouble("constraint","bearing_capacity_safety_factor");
    r.bearingUtilization=r.qmax/allowable;
    r.contactUtilization=r.contactRatio>0?minContact/r.contactRatio:INF;
    r.stabilityUtilization=std::max({r.slidingUtilization,r.overturnUtilization,
                                      r.bearingUtilization,r.contactUtilization});
    if(r.contact.branch==Contact::Branch::None)r.stabilityUtilization=INF;
    r.stabilityPass=r.stabilityUtilization<=1;

    const auto stemCap=rcCapacity(d.stem.area,d.dStemMM,s,c);
    const auto baseCap=rcCapacity(d.base.area,d.dBaseMM,s,c);
    r.moments[0]=r.pa*H/3+r.pq*H/2+r.mwStem;
    const double hShear=d.dStemMM/1000.0;
    r.shears[0]=0.5*r.ka*s.gamma*(H-hShear)*(H-hShear)+
                 r.ka*q*(H-hShear)+
                 waterAbove(service?h.finalS:h.finalC,H,hShear);
    const double wSlab=gammaC*d.tb;
    const double wHeelDead=wSlab+s.gamma*H;
    r.moments[1]=std::abs(integrateNet(r.contact,0,d.toe,wSlab,d.toe).second);
    const double shearToeX=std::max(d.toe-d.dBaseMM/1000.0,0.0);
    r.shears[1]=std::abs(integrateNet(r.contact,0,shearToeX,wSlab,d.toe).first);
    r.moments[2]=std::abs(integrateNet(r.contact,d.xBack,d.B,wHeelDead,d.xBack).second);
    const double shearHeelX=std::min(d.xBack+d.dBaseMM/1000.0,d.B);
    r.shears[2]=std::abs(integrateNet(r.contact,shearHeelX,d.B,wHeelDead,d.xBack).first);
    const double demandFactor=c.getSettingDouble("analysis","structural_demand_factor");
    r.dc={demandFactor*std::abs(r.moments[0])/stemCap.moment,
          demandFactor*std::abs(r.shears[0])/stemCap.shear,
          demandFactor*r.moments[1]/baseCap.moment,
          demandFactor*r.shears[1]/baseCap.shear,
          demandFactor*r.moments[2]/baseCap.moment,
          demandFactor*r.shears[2]/baseCap.shear};
    r.dcMax=*std::max_element(r.dc.begin(),r.dc.end());
    r.structuralPass=r.dcMax<=c.getSettingDouble("constraint","maximum_structural_demand_capacity_ratio");
    return r;
}

/**
 * Store the deterministic excavation, concrete, backfill, formwork, and reinforcement quantities.
 */
struct Quantities {
    double excavation=0,concrete=0,backfill=0,formwork=0,steelKg=0;
};

/**
 * Calculate material and construction quantities using the actual selected reinforcement areas.
 */
inline Quantities computeQuantities(const Design& d,const Config& c) {
    Quantities a;
    const double length=c.getFixed("wall_length_basis");
    const double H=c.getFixed("wall_height");
    const double baseV=d.B*d.tb*length;
    const double keyV=c.getFixed("shear_key_width")*d.dk*length;
    const double stemV=d.stemArea*length;
    a.excavation=(c.getFixed("foundation_embedment_depth")+d.tb+d.dk)*
                 (d.B+2*c.getFixed("excavation_working_space"))*length;
    a.concrete=baseV+keyV+stemV;
    a.backfill=std::max(0.0,a.excavation-baseV-keyV-
               d.ts*c.getFixed("foundation_embedment_depth")*length);
    a.formwork=(2*H+2*d.tb+2*d.dk)*length;
    a.steelKg=1e-6*length*(d.stem.area*H+d.base.area*(d.heel+d.toe))*
              c.getFixed("steel_density")*c.getFixed("secondary_steel_mass_factor");
    return a;
}

/**
 * Calculate scenario-specific total construction cost from deterministic quantities and sampled concrete and steel cost multipliers.
 */
inline double totalCost(const Quantities& q,const Scenario& s,const Config& c) {
    const double direct=18*q.excavation+155*s.etaConcreteCost*q.concrete+
         1.35*s.etaSteelCost*q.steelKg+42*q.formwork+28*q.backfill+
         850*c.getFixed("wall_length_basis");
    return 1.12*direct;
}

/**
 * Calculate costs in three prescribed deterministic price worlds for minimax regret; these worlds are separate from Monte Carlo realizations.
 */
inline std::array<double,3> deterministicCosts(const Quantities& q,const Config& c) {
    const double fixed=850*c.getFixed("wall_length_basis");
    return {{
        1.12*(18*q.excavation+155*1.30*q.concrete+1.35*1.00*q.steelKg+
              42*1.00*q.formwork+28*1.00*q.backfill+fixed),
        1.12*(18*1.15*q.excavation+155*1.00*q.concrete+1.35*1.35*q.steelKg+
              42*1.00*q.formwork+28*1.10*q.backfill+fixed),
        1.12*(18*1.40*q.excavation+155*1.05*q.concrete+1.35*1.10*q.steelKg+
              42*1.20*q.formwork+28*1.15*q.backfill+fixed)
    }};
}

/**
 * Compute transformed cracked-section inertia by solving for the neutral-axis depth and then integrating sectional inertia.
 */
inline double solveCrackedInertia(double As,double d,double EConcrete,double ESteel) {
    const double b=1000.0, n=ESteel/EConcrete, nAs=n*As;
    const double root=std::sqrt(nAs*nAs+2*b*nAs*d);
    const double neutralAxis=2*nAs*d/(nAs+root);
    if(!(neutralAxis>0 && neutralAxis<d)) return std::numeric_limits<double>::quiet_NaN();
    const double diff=d-neutralAxis;
    return b*neutralAxis*neutralAxis*neutralAxis/3+nAs*diff*diff;
}

/**
 * Hold the outcome of one design in one scenario: construction/service checks, cracking, stiffness, rotation, cost, and worst utilization.
 */
struct Outcome {
    StageResult construction,service;
    bool crackedC=false,crackedS=false,success=false;
    double EConcrete=0,crackingMoment=0,inertia=0,rotation=INF,rotationUtilization=INF;
    double cost=0, V=INF;
};

/**
 * Evaluate one design-scenario pair across both stages, using the construction crack state for service checks and determining overall success.
 */
inline Outcome evaluate(const Design& d,const Scenario& s,const Config& c,
                        const Hydraulic& h,const Quantities& q) {
    Outcome out;
    out.construction=evaluateStage(d,s,c,h,false);
    out.service=evaluateStage(d,s,c,h,true);
    const double H=c.getFixed("wall_height");
    const double stemMM=d.ts*1000;
    const double b=1000;
    out.EConcrete=4700*std::sqrt(s.fc)*s.etaE;
    const double rupture=0.62*std::sqrt(s.fc);
    const double Ig=b*stemMM*stemMM*stemMM/12;
    out.crackingMoment=(rupture*Ig/(stemMM/2))/1e6;
    out.crackedC=std::abs(out.construction.moments[0])>=out.crackingMoment;
    out.crackedS=out.crackedC || std::abs(out.service.moments[0])>=out.crackingMoment;
    out.inertia=out.crackedS?
        solveCrackedInertia(d.stem.area,d.dStemMM,out.EConcrete,c.getFixed("steel_modulus")):Ig;
    if(std::isfinite(out.inertia) && out.inertia>0)
        out.rotation=std::abs(out.service.moments[0])*1e6*(H*1000)/
            (3*out.EConcrete*out.inertia);
    out.rotationUtilization=out.rotation/c.getFixed("service_stem_rotation_limit");
    out.success=out.construction.stabilityPass && out.construction.structuralPass &&
                out.service.stabilityPass && out.service.structuralPass &&
                out.rotationUtilization<=1;
    out.cost=totalCost(q,s,c);
    out.V=std::max({out.construction.stabilityUtilization,
                    out.service.stabilityUtilization,
                    out.construction.dcMax,out.service.dcMax,out.rotationUtilization});
    return out;
}

} // namespace rw