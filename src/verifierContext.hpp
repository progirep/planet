#ifndef __VERIFIER_CONTEXT__
#define __VERIFIER_CONTEXT__

#include <map>
#include <vector>
#include <list>
#include <set>
#include "lpSolvers.hpp"
#include "supersetdatabase.hpp"


typedef enum { LINEAR, RELU, INPUT, MAXPOOL } NodeType;
typedef enum { LEQ } ComparisonType;


class VerificationProblem {
private:

    // Network
    std::map<std::string,unsigned int> nodeNumberLookup;
    std::vector<std::string> nodeNames;
    std::vector<NodeType> nodeTypes;
    std::vector<std::vector<std::pair<unsigned int, double> > > nodeConnectionIn;
    std::vector<double> nodeBias;
    std::vector<unsigned int> nodeNofPhases;
    unsigned int nofNodesWithPhases;

    // Constraints
    std::vector<std::tuple<ComparisonType,double,std::vector<std::pair<unsigned int, double> > > > constraints;

    // Initial estimates on the neural values
    std::vector<std::pair<double,double> > initialNeuronLimitBounds;

    // Fixed optimization function for network solving
    LPREALTYPE *nullOptimizationFunctionForLPSolving;

    // Data to make the current states of SAT and LP solver during the integrated search process
    // persistent
    SupersetDatabase cachedOptima;
    std::vector<int> lastPartialVariableValuationGivenToCheckPartialNodeFixtureInLPRelaxation;

    // Fixed-throughout-the-solver-run information about the connection of the SAT instance to the NN-verification-problem
    std::vector<int> startingSATVarsPhases;
    unsigned int nofSATVars;
    std::vector<std::pair<int,int> > satVarToNodeAndPhaseMapper;

    // Main LP Problem for checking a partial fixture -- while an assignment is gradually built, it is simply updated
    LPProblem lpPartialFixture;
    std::set<int> literalsAlreadyFixed; // the "lpPartialFixture" needs to be rebuilt if this is empty


protected:
    void addConstraintsToLPInstance(LPProblem &lp, unsigned int nofLPVars);
    void addNodePropagationInequationsToLPInstance(LPProblem &lp, unsigned int nofLPVars);
    void imposeNeuronValueBoundsInLPInstance(LPProblem &lp);
    void addReLUFixture(LPProblem &lp, unsigned int nodeNumber, unsigned int phase, unsigned int nofLPVars);
    void addReLUFixtureWithSlackVar(LPProblem &lp, unsigned int nodeNumber, unsigned int phase, int slackVar, unsigned int nofLPVars);
    void addMaxPoolFixture(LPProblem &lp,unsigned int nofNumber,const std::set<int> knownPhaseComponents, unsigned int nofLPVars);
    void addMaxPoolFixtureWithSlackVar(LPProblem &lp,unsigned int nofNumber,const std::set<int> knownPhaseComponents, int slackVar, unsigned int nofLPVars);


public:

    // Call-back function for the solver to check a partial phase valuation for "sanity".
    bool checkPartialNodeFixtureInLPRelaxation(std::vector<int> partialValuation, int nofLevel0VarsInPartialValuation, std::list<std::vector<int> > &addedClauses);

    // Call-back function for the solver to check data
    bool checkPartialNodeFixtureInSimplePropagation(std::vector<int> currentAssignment, std::list<std::vector<int> > &addedClauses);

    // Call-back function for backtracking.
    void backtrack() { literalsAlreadyFixed.clear(); }

    VerificationProblem() : nullOptimizationFunctionForLPSolving(0), cachedOptima(1) {}
    ~VerificationProblem() {
        if (nullOptimizationFunctionForLPSolving!=0) delete[] nullOptimizationFunctionForLPSolving;
    }

    void load(std::string inFile);
    bool computeInitialNeuronLimitBounds();
    bool solveLPRelavation(std::vector<bool> const &nodesToApproximate, std::vector<int> nodePhases, std::list<std::vector<int> > &newFixationsThatAreAdmissibleInTheRelaxation, std::list<std::vector<int> > &newFixationsThatAreNotAdmissibleInTheRelaxation);
    void verify();
    void printSMTLibInstance(bool includingApproximationConstraints);
    void printILPInstance(bool includingApproximationConstraints);
    inline unsigned int getNofSATVars() { return nofSATVars; }

};





















#endif
