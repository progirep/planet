#include "verifierContext.hpp"
#include <fstream>
#include <iostream>
#include <cassert>
#include <set>
#include <cmath>
#include <limits>
#include <valgrind/callgrind.h>
#include "minisat2/Solver.h"

// LP Solving Parameters
#define AGGREGATED_CHANGE_LIMIT 0.1
#define EPSILON_FOR_CALLING_AN_APPROXIMATED_NODE_VALUE_TO_BE_PRECISE_RELATIVE_TO_BOUNDS 0.00000001



/**
 * @brief Load an input file into the "VerificationProblem" class.
 * @param inFileName The file name
 */
void VerificationProblem::load(std::string inFileName) {

    nofNodesWithPhases = 0;

    std::ifstream inFile(inFileName);
    if (inFile.fail()) throw "Error opening input file.";
    std::string currentLine;
    while (std::getline(inFile,currentLine)) {

        // Trim CurrentLine
        while ((currentLine.size()>0) && (currentLine[0]==' ')) currentLine = currentLine.substr(1,std::string::npos);
        while ((currentLine.size()>0) && (currentLine[currentLine.size()-1]==' ')) currentLine = currentLine.substr(0,currentLine.size()-1);

        if (currentLine.size()>0) {

            // Parse this line
            std::istringstream thisLineStream(currentLine);
            std::string linePrefix;
            thisLineStream >> linePrefix;

            if (linePrefix=="Input") {
                // Input nodes
                std::string nodeName;
                thisLineStream >> nodeName;
                if (nodeNumberLookup.count(nodeName)>0) {
                    std::ostringstream err;
                    err << "Multiply defined node: " << nodeName;
                    throw err.str();
                }
                nodeNumberLookup[nodeName] = nodeNames.size();
                nodeNames.push_back(nodeName);
                nodeTypes.push_back(INPUT);
                nodeBias.push_back(0.0);
                nodeConnectionIn.push_back(std::vector<std::pair<unsigned int, double> >());
                nodeNofPhases.push_back(0);
            } else if (linePrefix=="MaxPool") {

                std::string nodeName;
                thisLineStream >> nodeName;
                if (nodeNumberLookup.count(nodeName)>0) {
                    std::ostringstream err;
                    err << "Multiply defined node: " << nodeName;
                    throw err.str();
                }
                nodeNumberLookup[nodeName] = nodeNames.size();
                nodeNames.push_back(nodeName);
                nodeTypes.push_back(MAXPOOL);
                nofNodesWithPhases++;
                nodeBias.push_back(0.0);
                nodeConnectionIn.push_back(std::vector<std::pair<unsigned int, double> >());

                while (!thisLineStream.eof()) {
                    std::string src;
                    thisLineStream >> src;
                    auto it = nodeNumberLookup.find(src);
                    if (it==nodeNumberLookup.end()) {
                        std::ostringstream err;
                        err << "Error: Did not find previously defined node '" << src << "'";
                        throw err.str();
                    }
                    nodeConnectionIn.back().push_back(std::pair<unsigned int,double>(it->second,1.0));
                }
                nodeNofPhases.push_back(nodeConnectionIn.back().size());

            } else if ((linePrefix=="ReLU") || (linePrefix=="Linear")) {

                // Main nodes
                std::string nodeName;
                thisLineStream >> nodeName;
                if (nodeNumberLookup.count(nodeName)>0) {
                    std::ostringstream err;
                    err << "Multiply defined node: " << nodeName;
                    throw err.str();
                }
                nodeNumberLookup[nodeName] = nodeNames.size();
                nodeNames.push_back(nodeName);
                if (linePrefix=="ReLU") {
                    nodeTypes.push_back(RELU);
                    nofNodesWithPhases++;
                } else if (linePrefix=="Linear")
                    nodeTypes.push_back(LINEAR);
                else
                    throw "Internal Error.";
                nodeConnectionIn.push_back(std::vector<std::pair<unsigned int, double> >());
                double weight;
                thisLineStream >> weight;
                nodeBias.push_back(weight);
                nodeNofPhases.push_back((linePrefix=="ReLU")?2:0);

                // Parse node list.
                while (!thisLineStream.eof()) {
                    thisLineStream >> weight;
                    if (thisLineStream.bad()) throw "Error reading weight.";
                    std::string src;
                    thisLineStream >> src;
                    auto it = nodeNumberLookup.find(src);
                    if (it==nodeNumberLookup.end()) {
                        std::ostringstream err;
                        err << "Error: Did not find previously defined node '" << src << "'";
                        throw err.str();
                    }
                    nodeConnectionIn.back().push_back(std::pair<unsigned int,double>(it->second,weight));
                }
            } else if (linePrefix=="Assert") {

                // Assert statement
                std::string operation;
                thisLineStream >> operation;
                ComparisonType comparisonType;
                bool negateAllWeights;
                if (operation=="<=") {
                    comparisonType = LEQ;
                    negateAllWeights = false;
                } else if (operation==">=") {
                    comparisonType = LEQ;
                    negateAllWeights = true;
                } else {
                    std::ostringstream err;
                    err << "Error: Did not understand comparison operation '" << operation << "'";
                    throw err.str();
                }

                double constant;
                thisLineStream >> constant;
                if (negateAllWeights) constant *= -1.0;
                std::vector<std::pair<unsigned int, double> > factors;

                // Parse node list.
                while (!thisLineStream.eof()) {
                    double weight;
                    thisLineStream >> weight;
                    if (negateAllWeights) weight *= -1.0;
                    if (thisLineStream.bad()) throw "Error reading weight.";
                    std::string src;
                    thisLineStream >> src;
                    auto it = nodeNumberLookup.find(src);
                    if (it==nodeNumberLookup.end()) {
                        std::ostringstream err;
                        err << "Error: Did not find previously defined node '" << src << "'";
                        throw err.str();
                    }
                    factors.push_back(std::pair<unsigned int,double>(it->second,weight));
                }

                constraints.push_back(std::make_tuple(comparisonType,constant,factors));

            } else {
                std::ostringstream err;
                err << "Error: Did not understand line prefix " << linePrefix;
                throw err.str();
            }
        }
    }

    if (inFile.bad()) throw "Error reading input file.";
    }



/**
 * @brief As from the constraints, the limits on the input and output values are not easily recognizable,
 *        we need to compute them using an LP solver once.
 * @return TRUE if successful, and FALSE if the instance is found to be UNSAT.
 */
bool VerificationProblem::computeInitialNeuronLimitBounds() {

    // ====================================================================================
    // First, we need to find upper and lower bounds for all input variables
    // ====================================================================================
    initialNeuronLimitBounds.resize(nodeTypes.size());
    LPProblem lp(nodeTypes.size());
    lp.addRowModeOn();
    addConstraintsToLPInstance(lp,nodeTypes.size());
    lp.addRowModeOff();

    for (unsigned int i=0;i<nodeTypes.size();i++) {
        if (nodeTypes[i]==INPUT) {

            // Lower bound
            LPREALTYPE row[nodeTypes.size()+1];
            for (unsigned int j=0;j<=nodeTypes.size();j++) row[j] = 0.0;
            row[i+1]=1.0;
            lp.setMinim();
            lp.setObjFn(row);

            int lpRes = lp.solve();
            if (lpRes==LPUNBOUNDED) {
                std::ostringstream err;
                err << "Input neuron " << nodeNames[i] << " has an LPUNBOUNDED domain. (A)";
                throw err.str();
            } else if (lpRes==LPINFEASIBLE) {
                return false;
            }
            //lp.writeLP("/tmp/jo.txt");
            double lowerBound = lp.getObjective();

            // Upper bound
            lp.setMaxim();
            lp.setObjFn(row);
            lpRes = lp.solve();
            if (lpRes==LPUNBOUNDED) {
                std::ostringstream err;
                err << "Input neuron " << nodeNames[i] << " has an LPUNBOUNDED domain. (B)";
                throw err.str();
            } else if (lpRes==LPINFEASIBLE) {
                return false;
            } else if ((lpRes!=LPOPTIMUM)) {
                throw "Error solving the LP model. (B)";
            }
            double upperBound = lp.getObjective();

            if ((upperBound==0.0) && (lowerBound>upperBound)) lowerBound = 0.0; // Fix imprecision
            if ((lowerBound==0.0) && (lowerBound>upperBound)) upperBound = 0.0; // Fix imprecision

            initialNeuronLimitBounds[i] = std::pair<double,double>(lowerBound,upperBound);
        }
    }

    // ====================================================================================
    // Now, we propagate the upper and lower bounds through the network
    // (initial estimate)
    // ====================================================================================
    for (unsigned int i=0;i<nodeTypes.size();i++) {
        if (nodeTypes[i]==INPUT) {
            // Nothing
        } else if (nodeTypes[i]==RELU) {
            double max = 0.0;
            double min = 0.0;
            for (auto pair : nodeConnectionIn[i]) {
                if (pair.second<0.0) {
                    min += pair.second*initialNeuronLimitBounds[pair.first].second;
                    max += pair.second*initialNeuronLimitBounds[pair.first].first;
                } else {
                    min += pair.second*initialNeuronLimitBounds[pair.first].first;
                    max += pair.second*initialNeuronLimitBounds[pair.first].second;
                }
            }
            initialNeuronLimitBounds[i] = std::pair<double,double>(std::max(0.0,min+nodeBias[i]),std::max(0.0,max+nodeBias[i]));
        } else if (nodeTypes[i]==LINEAR) {
            double min = 0.0;
            double max = 0.0;
            for (auto pair : nodeConnectionIn[i]) {
                if (pair.second<0.0) {
                    min += pair.second*initialNeuronLimitBounds[pair.first].second;
                    max += pair.second*initialNeuronLimitBounds[pair.first].first;
                } else {
                    min += pair.second*initialNeuronLimitBounds[pair.first].first;
                    max += pair.second*initialNeuronLimitBounds[pair.first].second;
                }
            }
            initialNeuronLimitBounds[i] = std::pair<double,double>(min+nodeBias[i],max+nodeBias[i]);
        } else if (nodeTypes[i]==MAXPOOL) {
            double min = std::numeric_limits<double>::max();
            double max = -1*std::numeric_limits<double>::max();
            for (auto pair : nodeConnectionIn[i]) {
                min = std::min(min,initialNeuronLimitBounds[pair.first].first);
                max = std::max(max,initialNeuronLimitBounds[pair.first].second);
            }
            initialNeuronLimitBounds[i] = std::pair<double,double>(min,max);

        } else {
            throw "Unknown node type.";
        }
    }


    // Print (debugging)
    std::cerr << "Node limits:\n";
    for (unsigned int i=0;i<nodeTypes.size();i++) {
        if (nodeTypes[i]==LINEAR)
            std::cerr << "+ ";
        else
            std::cerr << "- ";
        std::cerr << nodeNames[i] << ", min: " << initialNeuronLimitBounds[i].first << ", max: " << initialNeuronLimitBounds[i].second << std::endl;
    }

    // ====================================================================================
    // Now, we propagate the upper and lower bounds through the network
    // for that, we add the node constraints to the LP instance
    // We do that until there is no change anymore.
    // ====================================================================================
    double aggregatedChange = 1.0;
    unsigned int nofIterations = 0;

    while ((aggregatedChange > AGGREGATED_CHANGE_LIMIT) && ((nofIterations*nodeNames.size()<5000) || (nofIterations<3))) {

        // Remove maxpool connections that are never taken -- this makes the linear approximation tighter
        for (unsigned int i=0;i<nodeTypes.size();i++) {
            if (nodeTypes[i]==MAXPOOL) {
                std::vector<std::pair<unsigned int, double> > localIn;
                for (auto &it : nodeConnectionIn[i]) {
                    if (initialNeuronLimitBounds[it.first].second + EPSILON_FOR_CALLING_AN_APPROXIMATED_NODE_VALUE_TO_BE_PRECISE_RELATIVE_TO_BOUNDS < initialNeuronLimitBounds[i].first) {
                        // out of bounds...
                    } else {
                        localIn.push_back(it);
                    }
                }
                nodeNofPhases[i] = localIn.size();
                nodeConnectionIn[i] = localIn;
                if (nodeConnectionIn[i].size()==0) {
                    // Can happen in UNSAT instances if none of the incoming
                    // nodes supply as much flows as needed.
                    std::cerr << "\n";
                    return false;
                }
            }
        }

        aggregatedChange = 0.0;

        LPProblem lp(nodeTypes.size());
        lp.addRowModeOn();
        addConstraintsToLPInstance(lp,nodeTypes.size());
        addNodePropagationInequationsToLPInstance(lp,nodeTypes.size());
        imposeNeuronValueBoundsInLPInstance(lp);
        lp.addRowModeOff();

        nofIterations++;

        try {

            for (unsigned int i=0;i<nodeTypes.size();i++) {
                if ((nodeTypes[i]==RELU) || (nodeTypes[i]==MAXPOOL) || (nodeTypes[i]==LINEAR) || (nodeTypes[i]==INPUT)) {

                    std::cerr << ".";

                    // Lower bound
                    LPREALTYPE row[nodeTypes.size()+1];
                    double lowerBound;
                    for (unsigned int j=0;j<=nodeTypes.size();j++) row[j] = 0.0;
                    row[i+1]=1.0;
                    lp.setObjFn(row);

                    lp.setMinim();
                    //lp.writeLP("/tmp/jo.txt");
                    int lpRes = lp.solve();
                    if (lpRes==LPUNBOUNDED) {
                        std::ostringstream err;
                        err << "Neuron " << nodeNames[i] << " has an LPUNBOUNDED lower bound.";
                        throw err.str();
                    } else if (lpRes==LPINFEASIBLE) {
                        //lp.writeLP("/tmp/infeasible.txt");
                        std::cerr << "\n";
                        return false;
                    } else if (lpRes!=LPOPTIMUM) {
                        throw "Internal error (336)";
                    }
                    lowerBound = lp.getObjective();

                    // ReLU processing.
                    if (nodeTypes[i]==RELU) {
                        lowerBound = std::max(0.0,lowerBound);
                        if (lowerBound>EPSILON_FOR_CALLING_AN_APPROXIMATED_NODE_VALUE_TO_BE_PRECISE_RELATIVE_TO_BOUNDS) {
                            nodeTypes[i] = LINEAR;
                            nodeNofPhases[i] = 0;
                            std::cerr << "Making " << nodeNames[i] << " linear.\n";
                            addReLUFixture(lp, i, 2, nodeTypes.size());
                        }
                    }

                    // Upper bound
                    lp.setMaxim();
                    lp.setObjFn(row);

                    // ReLU special case: We can assume the second phase here.
                    if (nodeTypes[i]==RELU) {
                        addReLUFixture(lp, i, 2, nodeTypes.size());
                    }

                    lpRes = lp.solve();
                    double upperBound = lp.getObjective();

                    //lp.writeLP("/tmp/jo.txt");
                    if (lpRes==LPUNBOUNDED) {
                        std::ostringstream err;
                        lp.writeLP("/tmp/fuzz.txt");
                        err << "Neuron " << nodeNames[i] << " has an LPUNBOUNDED upper bound.";
                        assert(LPUNBOUNDED==lp.solve());
                        throw err.str();
                    } else if (lpRes==LPINFEASIBLE) {
                        if (nodeTypes[i]!=RELU) {
                            std::cerr << "\n";
                            return false;
                        } else {
                            // Fix the relu to a constant node
                            lp.deleteLastRow();

                            // In an unsatisfiable instance, it may
                            // happen that previously, the lower bound was
                            // >0. But then it's because it's unsat
                            if (initialNeuronLimitBounds[i].first>EPSILON_FOR_CALLING_AN_APPROXIMATED_NODE_VALUE_TO_BE_PRECISE_RELATIVE_TO_BOUNDS) return false;

                            upperBound = 0.0;
                            // Do not make this a linear node as we are losing information then...
                            addReLUFixture(lp, i, 1, nodeTypes.size());
                        }
                    } else if (lpRes!=LPOPTIMUM) {
                        throw "Internal error (351)";
                    }

                    // ReLU special case: Remove the added constraint
                    if (nodeTypes[i]==RELU) {
                        lp.deleteLastRow();
                        upperBound = std::max(upperBound,0.0);
                    }

                    aggregatedChange += std::fabs(lowerBound - initialNeuronLimitBounds[i].first);
                    aggregatedChange += std::fabs(upperBound - initialNeuronLimitBounds[i].second);

                    if (upperBound<lowerBound) {
                        // Can happen due to numerics
                        if ((upperBound==0.0) || (lowerBound==0.0))
                            initialNeuronLimitBounds[i] = std::pair<double,double>(0.0,0.0);
                        else
                            initialNeuronLimitBounds[i] = std::pair<double,double>((upperBound+lowerBound)/2.0,(upperBound+lowerBound)/2.0);
                    } else {
                        initialNeuronLimitBounds[i] = std::pair<double,double>(lowerBound,upperBound);
                    }
                } else {
                    throw "Unsupported (369)";
                }
            }
            std::cerr << "\n"; // The dots need to end here.


        } catch (NumericalFailureException e) {
            // Ok, we have to stop here then.
            aggregatedChange = 0.0;
        }

        // Print (debugging)
        std::cerr << "Node limits:\n";
        for (unsigned int i=0;i<nodeTypes.size();i++) {
            if (nodeTypes[i]==LINEAR)
                std::cerr << "+ ";
            else
                std::cerr << "- ";
            std::cerr << nodeNames[i] << ", min: " << initialNeuronLimitBounds[i].first << ", max: " << initialNeuronLimitBounds[i].second << std::endl;
        }

    }


    // Finally, prepare the optimization function for every LP problem except for bound computation
    nullOptimizationFunctionForLPSolving = new LPREALTYPE[nodeNames.size()+1];
    nullOptimizationFunctionForLPSolving[0] = 0.0;
    // Work in two steps: First, solve with all nodes having a weight, and then remove the weights for linear and input nodes.
    for (unsigned int i=0;i<nodeTypes.size();i++) {
        nullOptimizationFunctionForLPSolving[i+1] = 0.0;
    }
    /*for (int i=nodeTypes.size()-1;i>=0;i--) {
        if (nodeTypes[i]==INPUT) {
            // Terminal
        } else if ((nodeTypes[i]==LINEAR) || (nodeTypes[i]==RELU)) {
            for (auto it : nodeConnectionIn[i]) {
                fixedOptimizationFunctionNetworkSolving[it.first+1] += std::abs(fixedOptimizationFunctionNetworkSolving[i+1]*it.second);
            }
        } else {
            throw "Unsupported";
        }
    }
    for (unsigned int i=0;i<nodeTypes.size();i++) {
        if ((nodeTypes[i]==INPUT) || (nodeTypes[i]==LINEAR)) {
            fixedOptimizationFunctionNetworkSolving[i+1] = 0.0;

        }
    }*/
    return true;
}



void VerificationProblem::addReLUFixture(LPProblem &lp, unsigned int nodeNumber, unsigned int phase, unsigned int nofLPVars) {
    if (phase==1) {
        // Fixed to 0
        // ---> Fix output to 0
        LPREALTYPE row[nofLPVars+1];
        for (unsigned int j=0;j<=nofLPVars;j++) row[j]=0;
        row[nodeNumber+1] = 1.0;
        lp.addEQConstraint(row,0.0);

        // ---> Fix input <= -1*nodeBias
        for (unsigned int j=0;j<=nofLPVars;j++) row[j]=0;
        for (auto pair : nodeConnectionIn[nodeNumber]) {
            row[pair.first+1] = pair.second;
        }
        lp.addLEQConstraint(row, -1*nodeBias[nodeNumber]);

    } else if (phase==2) {
        // Fixed above 0
        LPREALTYPE row[nofLPVars+1];
        for (unsigned int j=0;j<=nofLPVars;j++) row[j]=0;
        for (auto pair : nodeConnectionIn[nodeNumber]) {
            row[pair.first+1] = pair.second;
        }
        row[nodeNumber+1] = -1.0;
        lp.addEQConstraint(row, -1*nodeBias[nodeNumber]);
    } else {
        std::cerr << "Phase: " << phase << std::endl;
        throw "Unknown phase!";
    }
}

void VerificationProblem::addMaxPoolFixture(LPProblem &lp, unsigned int nodeNumber,const std::set<int> knownPhaseComponents, unsigned int nofLPVars) {

    // First constraint: The output must be smaller than the sum of the other nodes
    LPREALTYPE row[nofLPVars+1];
    for (unsigned int j=0;j<=nofLPVars;j++) row[j]=0;
    row[nodeNumber+1] = 1.0;
    for (auto pair : nodeConnectionIn[nodeNumber]) {
        row[pair.first+1] = -1.0;
        //std::cerr << "(IN" << pair.first+1 << ")";
    }
    std::set<int> unselectedNodes;
    for (auto lit : knownPhaseComponents) {
        if (lit<0) {
            row[nodeConnectionIn[nodeNumber][-1*lit-1].first+1] = 0.0;
            unselectedNodes.insert(nodeConnectionIn[nodeNumber][-1*lit-1].first+1);
            //std::cerr << "(R" << nodeConnectionIn[nodeNumber][-1*lit-1].first+1 << ")";
        } else {
            for (unsigned int i=0;i<nodeConnectionIn[nodeNumber].size();i++) {
                if (static_cast<int>(i+1)!=lit) {
                    row[nodeConnectionIn[nodeNumber][i].first+1] = 0.0;
                    unselectedNodes.insert(nodeConnectionIn[nodeNumber][i].first+1);
                    //std::cerr << "(RX" << lit << "," << nodeConnectionIn[nodeNumber][i].first+1 << ")";
                }
            }
        }
        assert(lit!=0);
    }

    // Compute minima
    LPREALTYPE sumMins = 0.0;
    LPREALTYPE maxMin = -1*std::numeric_limits<double>::max();
    for (auto pair : nodeConnectionIn[nodeNumber]) {
        if (row[pair.first+1] != 0.0) {
            double thisMin = initialNeuronLimitBounds[pair.first].first;
            sumMins += thisMin;
            if (thisMin>maxMin) {
                maxMin = thisMin;
            }
        }
    }
#ifndef NDEBUG
    if (maxMin==-1*std::numeric_limits<double>::max()) {
        std::cerr << "Node: " << nodeNames[nodeNumber] << std::endl;
        std::cerr << "Number of phases: " << nodeNofPhases[nodeNumber] << "\n";
        std::cerr << "Incoming edges:";
        for (auto it : nodeConnectionIn[nodeNumber]) {
            std::cerr << " " << nodeNames[it.first];
        }
        std::cerr << "\nLiterals:";
        for (auto it : knownPhaseComponents) std::cerr << " " << it;
        std::cerr << std::endl;
        std::cerr << "nodeConnectionIn.Size: " << nodeConnectionIn[nodeNumber].size() << std::endl;
    }

    assert(maxMin!=-1*std::numeric_limits<double>::max());
#endif

    lp.addLEQConstraint(row, maxMin-sumMins);

    // Second constraint: The output of the node must be >= the unselected values
    for (unsigned int j=0;j<=nofLPVars;j++) row[j]=0;
    row[nodeNumber+1] = 1.0;
    for (auto col : unselectedNodes) {
        row[col] = -1.0;
        lp.addGEQConstraint(row, 0.0);
        row[col] = 0.0;
    }
}

void VerificationProblem::addMaxPoolFixtureWithSlackVar(LPProblem &lp,unsigned int nodeNumber,const std::set<int> knownPhaseComponents, int slackVar, unsigned int nofLPVars) {

    // First constraint: The output must be smaller than the sum of the other nodes
    LPREALTYPE row[nofLPVars+1];
    for (unsigned int j=0;j<=nofLPVars;j++) row[j]=0;
    row[nodeNumber+1] = 1.0;
    row[slackVar+1] = -1.0;
    for (auto pair : nodeConnectionIn[nodeNumber]) {
        row[pair.first+1] = -1.0;
    }
    std::set<int> unselectedNodes;
    for (auto lit : knownPhaseComponents) {
        if (lit<0) {
            row[nodeConnectionIn[nodeNumber][-1*lit-1].first+1] = 0.0;
            unselectedNodes.insert(nodeConnectionIn[nodeNumber][-1*lit-1].first+1);
            //std::cerr << "(R" << nodeConnectionIn[nodeNumber][-1*lit-1].first+1 << ")";
        } else {
            for (unsigned int i=0;i<nodeConnectionIn[nodeNumber].size();i++) {
                if (static_cast<int>(i+1)!=lit) {
                    row[nodeConnectionIn[nodeNumber][i].first+1] = 0.0;
                    unselectedNodes.insert(nodeConnectionIn[nodeNumber][i].first+1);
                    //std::cerr << "(RX" << lit << "," << nodeConnectionIn[nodeNumber][i].first+1 << ")";
                }
            }
        }
        assert(lit!=0);
    }

    // Compute minima
    LPREALTYPE sumMins = 0.0;
    LPREALTYPE maxMin = -1*std::numeric_limits<double>::max();
    for (auto pair : nodeConnectionIn[nodeNumber]) {
        if (row[pair.first+1] != 0.0) {
            double thisMin = initialNeuronLimitBounds[pair.first].first;
            sumMins += thisMin;
            if (thisMin>maxMin) {
                maxMin = thisMin;
            }
        }
    }

    assert(maxMin!=-1*std::numeric_limits<double>::max());
    lp.addLEQConstraint(row, maxMin-sumMins);

    // Second constraint: The output of the node must be >= the unselected values
    for (unsigned int j=0;j<=nofLPVars;j++) row[j]=0;
    row[nodeNumber+1] = 1.0;
    row[slackVar+1] = 1.0;
    for (auto col : unselectedNodes) {
        row[col] = -1.0;
        lp.addGEQConstraint(row, 0.0);
        row[col] = 0.0;
    }
}


void VerificationProblem::addReLUFixtureWithSlackVar(LPProblem &lp, unsigned int nodeNumber, unsigned int phase, int slackVar, unsigned int nofLPVars) {
    if (phase==1) {
        // Fixed to 0
        // ---> Fix output to 0
        LPREALTYPE row[nofLPVars+1];
        for (unsigned int j=0;j<=nofLPVars;j++) row[j]=0;
        row[nodeNumber+1] = 1.0;
        row[slackVar+1] = -1.0;
        lp.addLEQConstraint(row,0.0);

        // ---> Fix input <= -1*nodeBias
        for (unsigned int j=0;j<=nofLPVars;j++) row[j]=0;
        for (auto pair : nodeConnectionIn[nodeNumber]) {
            row[pair.first+1] = pair.second;
        }
        row[slackVar+1] = -1.0;
        lp.addLEQConstraint(row, -1*nodeBias[nodeNumber]);

    } else if (phase==2) {
        // Fixed above 0
        LPREALTYPE row[nofLPVars+1];
        for (unsigned int j=0;j<=nofLPVars;j++) row[j]=0;
        for (auto pair : nodeConnectionIn[nodeNumber]) {
            row[pair.first+1] = pair.second;
        }
        row[nodeNumber+1] = -1.0;
        row[slackVar+1] = 1.0;
        lp.addGEQConstraint(row, -1*nodeBias[nodeNumber]);
        row[slackVar+1] = -1.0;
        lp.addLEQConstraint(row, -1*nodeBias[nodeNumber]);

    } else {
        std::cerr << "Phase: " << phase << std::endl;
        throw "Unknown phase!";
    }
}


/**
 * @brief Takes an LP problem object and adds the constraints on the input and output nodes given in the verification problem.
 * @param lp The LP problem object.
 */
void VerificationProblem::addConstraintsToLPInstance(LPProblem &lp, unsigned int nofLPVars) {
    for (const std::tuple<ComparisonType,double,std::vector<std::pair<unsigned int, double> > > &c : constraints) {
        LPREALTYPE row[nofLPVars+1];
        for (unsigned int i=0;i<=nofLPVars;i++) row[i] = 0.0;
        for (auto l : std::get<2>(c)) {
            //std::cerr << l.second << "*" << nodeNames[l.first] << " ";
            row[l.first+1] = l.second;
        }
        //std::cerr << ">= " << std::get<1>(c) << std::endl;
        lp.addGEQConstraint(row, std::get<1>(c));
        // std::cerr << "Const: row " << row[0] << "," << row[1] << "," <<row[2] << "," << row[3] << ">=" << std::get<1>(c) << std::endl;
    }
}

/**
 * @brief Takes an LP problem object and adds the upper and lower bounds stored in "initialNeuronLimitBounds" to the neurons.
 * @param lp The LP problem object.
 */
void VerificationProblem::imposeNeuronValueBoundsInLPInstance(LPProblem &lp) {
    for (unsigned int i=0;i<nodeTypes.size();i++) {
        lp.setBounds(i,initialNeuronLimitBounds[i].first,initialNeuronLimitBounds[i].second);
    }
}

/**
 * @brief Takes an LP problem object and adds the constraints that restrict the node values by the linear approximation of the verification problem. This function is not incremental.
 * @param lp The LP problem object.
 */
void VerificationProblem::addNodePropagationInequationsToLPInstance(LPProblem &lp, unsigned int nofLPVars) {
    for (unsigned int i=0;i<nodeTypes.size();i++) {
        if (nodeTypes[i]==INPUT) {
            // Nothing
        } else if (nodeTypes[i]==RELU) {
            // RELU -- Make a convex hull

            // GEQ 0
            LPREALTYPE row[nofLPVars+1];
            for (unsigned int j=0;j<=nofLPVars;j++) row[j] = 0.0;
            row[i+1] = 1.0;
            lp.addGEQConstraint(row, 0);

            // GEQ linear combination
            for (auto pair : nodeConnectionIn[i]) {
                row[pair.first+1] = pair.second;
            }
            row[i+1] = -1.0;
            lp.addLEQConstraint(row, -1*nodeBias[i]);

            // Compute theoretical minimum
            double min = nodeBias[i];
            double max = nodeBias[i];
            for (auto pair : nodeConnectionIn[i]) {
                if (pair.second<0.0) {
                    min += pair.second*initialNeuronLimitBounds[pair.first].second;
                    max += pair.second*initialNeuronLimitBounds[pair.first].first;
                } else {
                    min += pair.second*initialNeuronLimitBounds[pair.first].first;
                    max += pair.second*initialNeuronLimitBounds[pair.first].second;
                }
            }
            max = std::min(initialNeuronLimitBounds[i].second,max);

            //std::cerr << "RELUAddNodePropMinMax: " << min << "," << max << std::endl;

            if (max<=0) {
                for (unsigned int j=0;j<=nofLPVars;j++) row[j] = 0.0;
                row[i+1] = 1.0;
                lp.addLEQConstraint(row, 0);

                // Don't declare this node to be linear! The SAT solver should do this.

            } else {
                if (min<0) {
                    // Compute line between (min,0) and (max,max)
                    double factor = max/(max-min);
                    double offset = max-max*factor;

                    for (unsigned int j=0;j<=nofLPVars;j++) row[j] = 0.0;
                    for (auto pair : nodeConnectionIn[i]) {
                        row[pair.first+1] = -1.0*factor*pair.second;
                    }
                    row[i+1] = 1.0;

                    lp.addLEQConstraint(row, nodeBias[i]*factor+offset);


                } else {
                    // Make it linear
                    // std::cerr << "Apply min,max>0 case\n";
                    LPREALTYPE row[nofLPVars+1];
                    for (unsigned int j=0;j<=nofLPVars;j++) row[j] = 0.0;
                    for (auto pair : nodeConnectionIn[i]) {
                        row[pair.first+1] = pair.second;
                    }
                    row[i+1] = -1.0;
                    lp.addLEQConstraint(row, -1*nodeBias[i]);

                    // Don't declare this node to be linear! The SAT solver should do this.
                }
            }

        } else if (nodeTypes[i]==LINEAR) {
            LPREALTYPE row[nofLPVars+1];
            for (unsigned int j=0;j<=nofLPVars;j++) row[j] = 0.0;
            for (auto pair : nodeConnectionIn[i]) {
                row[pair.first+1] = pair.second;
            }
            row[i+1] = -1.0;
            lp.addEQConstraint(row, -1*nodeBias[i]);
        } else if (nodeTypes[i]==MAXPOOL) {

            // First, of all the output must be greater than all inputs
            LPREALTYPE row[nofLPVars+1];
            LPREALTYPE sumMins = 0.0;
            LPREALTYPE maxMin = -1*std::numeric_limits<double>::max();
            for (unsigned int j=0;j<=nofLPVars;j++) row[j] = 0.0;
            row[i+1] = 1.0;
            for (auto pair : nodeConnectionIn[i]) {
                row[pair.first+1] = -1.0;
                lp.addGEQConstraint(row, 0.0);
                row[pair.first+1] = 0.0;
                sumMins += initialNeuronLimitBounds[pair.first].first;
                if (initialNeuronLimitBounds[pair.first].first>maxMin) {
                    maxMin = initialNeuronLimitBounds[pair.first].first;
                }
            }

            // But then, the output must also be smaller than the sum of all inputs
            for (auto pair : nodeConnectionIn[i]) {
                row[pair.first+1] = -1.0;
            }
            lp.addLEQConstraint(row, maxMin-sumMins);

        } else {
            throw "Unknown node type (594).";
        }
    }
}


/**
 * @brief Main verification function -- essentially gives control to Minisat, which then calls "checkPartialSATValuationInLPRelaxation" from time to time.
 */
void VerificationProblem::verify() {

    // Allocate minisat object
    Minisat::Solver minisat;

    // Allocate Minisat vars and initialize "satVarToNodeAndPhaseMapper" -- to make sure that we can work with SAT var numbers starting with "1", we have
    // to initialize "satVarToNodeAndPhaseMapper" with a Pseudo-Elements
    satVarToNodeAndPhaseMapper.push_back(std::pair<int,int>(-1,-1));

    unsigned int nextFreeVar = 1;
    nofNodesWithPhases = 0;
    for (unsigned int i=0;i<nodeTypes.size();i++) {
        std::cerr << "c Starting var phases of " << nodeNames[i] << ": " << nextFreeVar << std::endl;
        startingSATVarsPhases.push_back(nextFreeVar);
        for (unsigned int j=0;j<nodeNofPhases[i];j++) {
            assert(satVarToNodeAndPhaseMapper.size() == nextFreeVar+j);
            satVarToNodeAndPhaseMapper.push_back(std::pair<int,int>(i,j));
        }
        nextFreeVar += nodeNofPhases[i];
        if (nodeNofPhases[i]>0) nofNodesWithPhases++;
    }
    for (unsigned int i=0;i<nextFreeVar-1;i++) minisat.newVar();
    nofSATVars = nextFreeVar-1;
    cachedOptima = SupersetDatabase(nofSATVars+1);

    // Add constraints that enforce that exactly one phase is selected
    for (unsigned int i=0;i<nodeTypes.size();i++) {
        // One is selected
        if (nodeNofPhases[i]>0) {
            std::vector<int> clauseOneIsSelected;
            //std::cerr << "CLAUSE for " << nodeNames[i] << ": ";
            for (unsigned int j=0;j<nodeNofPhases[i];j++) {
                //std::cerr << startingSATVarsPhases[i]+j << " ";
                clauseOneIsSelected.push_back(startingSATVarsPhases[i]+j);
            }
            //std::cerr << std::endl;
            minisat.addClause(clauseOneIsSelected);
        }
        // At most one is selected
        for (unsigned int j=0;j<nodeNofPhases[i];j++) {
            for (unsigned int k=j+1;k<nodeNofPhases[i];k++) {
                std::vector<int> clause;
                //std::cerr << "CLAUSE: " << -1*(int)(startingSATVarsPhases[i]+j) << " " << -1*(int)(startingSATVarsPhases[i]+k) << std::endl;
                clause.push_back(-1*(startingSATVarsPhases[i]+j));
                clause.push_back(-1*(startingSATVarsPhases[i]+k));
                minisat.addClause(clause);
            }
        }
    }

    // For ReLU nodes that have a maximum value of 0, we already know that the phase should be 0
    for (unsigned int j=0;j<nodeNames.size();j++) {
        if (nodeTypes[j]==RELU) {
            if (initialNeuronLimitBounds[j].second<=0.0) {
                std::vector<int> clause;
                clause.push_back(startingSATVarsPhases[j]);
                std::cerr << "CLAUSE RELUPHASE1: " << startingSATVarsPhases[j] << std::endl;
                minisat.addClause(clause);
            }
        }
    }


    CALLGRIND_START_INSTRUMENTATION;

    // Main solving function
    minisat.setVerificationProblem(this);
    bool isSAT = minisat.solve();

    if (!isSAT) {
        std::cout << "\nUNSAT\n";
    } else {
        std::cout << "\nSAT\n";
    }

}

/**
 * @brief Here, we try to find out new upper and lower bounds for nodes dynamically (without using the LP solver), so that
 *        we can detect implications between variables. We don't use the LP solver yet, and if this function
 *        returns new clauses, they should be added to the valuation before running the LP solver.
 * @param currentAssignment
 * @param clausesToPropagate
 * @returns True if no conflict clause has been found.
 */
bool VerificationProblem::checkPartialNodeFixtureInSimplePropagation(std::vector<int> currentAssignment, std::list<std::vector<int> > &clausesToAdd) {
    std::vector<std::pair<double,double> > refinedNeuronLimitBounds(nodeTypes.size());
    std::vector<std::pair<double,double> > howMuchInflowIsAdmissible(nodeTypes.size());
    std::vector<std::set<int> > reasons(nodeTypes.size()); // NEGATED Literals connected to this one.
    std::set<int> currentAssignmentAsSet(currentAssignment.begin(),currentAssignment.end());

    //return true;

    // Iterate over all nodes and propagate DOWNWARDS
    for (unsigned int i=0;i<nodeTypes.size();i++) {
        if (nodeTypes[i]==INPUT) {
            // Nothing to do. Just copy bounds
            refinedNeuronLimitBounds[i] = initialNeuronLimitBounds[i];
            howMuchInflowIsAdmissible[i] = initialNeuronLimitBounds[i];
        } else if (nodeTypes[i]==LINEAR) {
            // Compute new upper and lower limits.
            double min = 0.0;
            double max = 0.0;
            for (auto pair : nodeConnectionIn[i]) {
                if (pair.second<0.0) {
                    min += pair.second*refinedNeuronLimitBounds[pair.first].second;
                    max += pair.second*refinedNeuronLimitBounds[pair.first].first;
                } else {
                    min += pair.second*refinedNeuronLimitBounds[pair.first].first;
                    max += pair.second*refinedNeuronLimitBounds[pair.first].second;
                }
                reasons[i].insert(reasons[pair.first].begin(),reasons[pair.first].end());
            }
            min += nodeBias[i];
            max += nodeBias[i];
            min = std::max(min,initialNeuronLimitBounds[i].first);
            max = std::min(max,initialNeuronLimitBounds[i].second);
            refinedNeuronLimitBounds[i] = std::pair<double,double>(min,max);
            howMuchInflowIsAdmissible[i] = std::pair<double,double>(min,max);

        } else if (nodeTypes[i]==RELU) {
            // Compute new upper and lower limits.
            double min = 0.0;
            double max = 0.0;
            for (auto pair : nodeConnectionIn[i]) {
                if (pair.second<0.0) {
                    min += pair.second*refinedNeuronLimitBounds[pair.first].second;
                    max += pair.second*refinedNeuronLimitBounds[pair.first].first;
                } else {
                    min += pair.second*refinedNeuronLimitBounds[pair.first].first;
                    max += pair.second*refinedNeuronLimitBounds[pair.first].second;
                }
                reasons[i].insert(reasons[pair.first].begin(),reasons[pair.first].end());
            }
            min += nodeBias[i];
            max += nodeBias[i];

            // If min and max are only inherited from initialNeuronLimitBounds, we don't really need the reasons
            if ((min <= initialNeuronLimitBounds[i].first) && (max >= initialNeuronLimitBounds[i].second)) {
                min = initialNeuronLimitBounds[i].first;
                max = initialNeuronLimitBounds[i].second;
                reasons[i].clear();
            }

            // Is is not yet known that this is >0?
            if ((min>EPSILON_FOR_CALLING_AN_APPROXIMATED_NODE_VALUE_TO_BE_PRECISE_RELATIVE_TO_BOUNDS) /*&& (max>min)*/) {
                if (currentAssignmentAsSet.count(startingSATVarsPhases[i]+1)==0) {
                    // Found new necessary phase!
                    // ---> Build implied clause
                    //std::cerr << "RELU MUST BE POSITIVE: " << nodeNames[i] << std::endl;
                    //std::cerr << "RELU-Enforcing clause: " << nodeNames[i] << "---> " << startingSATVarsPhases[i]+1 << "(" << min << "," << max << "," << currentAssignment.size() << ")" << std::endl;
                    std::set<int> sortedSetReasons(reasons[i].begin(),reasons[i].end());
                    std::vector<int> newClause;
                    newClause.push_back(startingSATVarsPhases[i]+1);
                    newClause.insert(newClause.end(),sortedSetReasons.begin(),sortedSetReasons.end());
#ifndef NDEBUG
                    for (unsigned int i=1;i<newClause.size();i++) assert(currentAssignmentAsSet.count(-1*newClause[i])>0);
#endif
                    // Conflict clause or regular clause?
                    if (currentAssignmentAsSet.count(-1*(startingSATVarsPhases[i]+1))>0) {
                        //std::cerr << "Added RELU Conflict clause!";
                        clausesToAdd.push_back(newClause);
                    } else {
                        clausesToAdd.push_back(newClause);
                    }
                }
            } else if (max<=0) {// -1*EPSILON_FOR_CALLING_AN_APPROXIMATED_NODE_VALUE_TO_BE_PRECISE_RELATIVE_TO_BOUNDS) {
                // Ok, this must be in phase 1
                if (currentAssignmentAsSet.count(startingSATVarsPhases[i])==0) {
                    // Found new necessary phase!
                    // ---> Build implied clause
                    //std::cerr << "RELU MUST BE POSITIVE: " << nodeNames[i] << std::endl;
                    //std::cerr << "RELU-FALSE-Enforcing clause: " << nodeNames[i] << "---> " << startingSATVarsPhases[i] << "(" << min << "," << max << "," << currentAssignment.size() << ")" << std::endl;
                    std::set<int> sortedSetReasons(reasons[i].begin(),reasons[i].end());
                    std::vector<int> newClause;
                    newClause.push_back(startingSATVarsPhases[i]);
                    newClause.insert(newClause.end(),sortedSetReasons.begin(),sortedSetReasons.end());

                    // Conflict clause or regular clause?
                    if (currentAssignmentAsSet.count(-1*(startingSATVarsPhases[i]))>0) {
                        //std::cerr << "Added RELU-False Conflict clause!";
                        clausesToAdd.push_back(newClause);
                    } else {
                        clausesToAdd.push_back(newClause);
                    }
                }
            } else {
                // Can it be that this phase is already fixed to be 0? Then fix it as well
                if (currentAssignmentAsSet.count(startingSATVarsPhases[i]+0)>0) {
                    max = 0.0;
                    reasons[i].insert(-1*startingSATVarsPhases[i]+0);
                }
            }

            // Apply RELU logic
            min = std::max(min,initialNeuronLimitBounds[i].first);
            max = std::min(max,initialNeuronLimitBounds[i].second);
            min = std::max(min,0.0);
            max = std::max(max,0.0);

            refinedNeuronLimitBounds[i] = std::pair<double,double>(min,max);
            if (min>EPSILON_FOR_CALLING_AN_APPROXIMATED_NODE_VALUE_TO_BE_PRECISE_RELATIVE_TO_BOUNDS) {
                howMuchInflowIsAdmissible[i] = std::pair<double,double>(min,max);
            } else {
                howMuchInflowIsAdmissible[i] = std::pair<double,double>(-1*std::numeric_limits<double>::max(),max);
            }
        } else if (nodeTypes[i]==MAXPOOL) {

            double min = refinedNeuronLimitBounds[nodeConnectionIn[i][0].first].first;
            double max = refinedNeuronLimitBounds[nodeConnectionIn[i][0].first].second;
            int minLimitingIncomingEdge=0;
            for (unsigned int j=1;j<nodeConnectionIn[i].size();j++) {
                if (refinedNeuronLimitBounds[nodeConnectionIn[i][j].first].first>min) {
                    min = refinedNeuronLimitBounds[nodeConnectionIn[i][j].first].first;
                    minLimitingIncomingEdge = j;
                }
                max = std::max(max,refinedNeuronLimitBounds[nodeConnectionIn[i][j].first].second);
            }

            /*std::cerr << "Refined neuron limits for maxpool node " << nodeNames[i] << ":\n";
            for (unsigned int j=0;j<nodeConnectionIn[i].size();j++) {
                std::cerr << "- " << refinedNeuronLimitBounds[nodeConnectionIn[i][j].first].first << " -- " << refinedNeuronLimitBounds[nodeConnectionIn[i][j].first].second << std::endl;
            }*/

            // Compute reasons
            for (auto pair : nodeConnectionIn[i]) {
                reasons[i].insert(reasons[pair.first].begin(),reasons[pair.first].end());
            }

            // Now which phases can already be deactivated?
            for (unsigned int j=0;j<nodeConnectionIn[i].size();j++) {
                if (refinedNeuronLimitBounds[nodeConnectionIn[i][j].first].second<min-EPSILON_FOR_CALLING_AN_APPROXIMATED_NODE_VALUE_TO_BE_PRECISE_RELATIVE_TO_BOUNDS) {
                    // This cannot be the phase then
                    if (currentAssignmentAsSet.count(-1*startingSATVarsPhases[i]-j)==0) {
                        std::vector<int> newClause;
                        std::set<int> sortedSetReasons(reasons[minLimitingIncomingEdge].begin(),reasons[minLimitingIncomingEdge].end());
                        sortedSetReasons.insert(reasons[i].begin(),reasons[i].end());
                        newClause.push_back(-1*startingSATVarsPhases[i]-j);
                        newClause.insert(newClause.end(),sortedSetReasons.begin(),sortedSetReasons.end());
                        //std::cerr << "Adding new MAXPOOL Forward Propagating Clause: " << newClause[0] << ", length: " << newClause.size() << "\n";
                        clausesToAdd.push_back(newClause);
                    }
                }
            }

            refinedNeuronLimitBounds[i] = std::pair<double,double>(min,max);
            howMuchInflowIsAdmissible[i] = std::pair<double,double>(-1*std::numeric_limits<double>::max(),max);

        } else {
            throw "Unsupported node type....604";
            // TODO: Right now, we don't really reduce the amount of dependencies in the clauses -- perhaps we also want sets everywhere? Or Bitvectors?
        }
    }

    // Iterate over the nodes and propagate *UPWARDS*
    for (int i=nodeTypes.size()-1;i>=0;i--) {
        if (nodeTypes[i]==INPUT) {
            // Nothing to do.
        } else if (nodeTypes[i]!=MAXPOOL) {
            // All other node have incoming edges that uniquely determine their output values

            // Compute lower and upper bounds on the non-fixed ingoing nodes
            double valueFixedMin = nodeBias[i];
            double valueFixedMax = nodeBias[i];
            std::vector<int> theseReasons;

            for (auto pair : nodeConnectionIn[i]) {
                if (pair.second<0.0) {
                    valueFixedMin += pair.second*refinedNeuronLimitBounds[pair.first].second;
                    valueFixedMax += pair.second*refinedNeuronLimitBounds[pair.first].first;
                } else {
                    valueFixedMax += pair.second*refinedNeuronLimitBounds[pair.first].second;
                    valueFixedMin += pair.second*refinedNeuronLimitBounds[pair.first].first;
                }
                theseReasons.insert(theseReasons.end(),reasons[pair.first].begin(),reasons[pair.first].end());
            }

            // How much flow do the other nodes need to contribute?
            //double restMin = refinedNeuronLimitBounds[i].first - valueFixedMin;
            //double restMax = refinedNeuronLimitBounds[i].second - valueFixedMax;
            //std::cerr << "Backward prop node " << nodeNames[i] << ": " << valueFixedMin << " < " << valueFixedMax << std::endl;

            // Now iterate over the incoming nodes and observe which values are forced by the fact that
            // the rest would not be able to supply enough flow if the phase of one incoming node is set "incorrectly".
            for (auto pair : nodeConnectionIn[i]) {
                if ((nodeTypes[pair.first]==INPUT) || (nodeTypes[pair.first]==LINEAR)) {
                    // Nothing to do.
                } else if (nodeTypes[pair.first]==RELU) {
                    // Only look at ReLUs if their values are not fixed already.
                    if ((currentAssignmentAsSet.count(-1*(startingSATVarsPhases[pair.first]))==0)
                            && (currentAssignmentAsSet.count(-1*(startingSATVarsPhases[pair.first]+1))==0)) {

                        // Does the RELU have to be >0?
                        if (valueFixedMax-refinedNeuronLimitBounds[pair.first].second*pair.second+EPSILON_FOR_CALLING_AN_APPROXIMATED_NODE_VALUE_TO_BE_PRECISE_RELATIVE_TO_BOUNDS<howMuchInflowIsAdmissible[i].first) {
                            // ---> Then propagate this fact.

                            //std::cerr << "RELU>0-Enforcing clause: " << nodeNames[pair.first] << std::endl;
                            //std::cerr << valueFixedMax << "from " << refinedNeuronLimitBounds[pair.first].second*pair.second << " to " << howMuchInflowIsAdmissible[i].first << std::endl;
                            //std::cerr << "for incoming node: " << nodeNames[i] << std::endl;
                            std::set<int> sortedSetReasons(reasons[i].begin(),reasons[i].end());
                            std::vector<int> newClause;
                            newClause.push_back(startingSATVarsPhases[pair.first]+1);
                            //std::cerr << "Phase: " << startingSATVarsPhases[pair.first]+1 << std::endl;
                            newClause.insert(newClause.end(),sortedSetReasons.begin(),sortedSetReasons.end());

                            clausesToAdd.push_back(newClause);
                        }
                    }
                } else if (nodeTypes[pair.first]==MAXPOOL) {
                    // For a maxpool layer, check if there is any incoming phase whose selection would lead to insufficient
                    // flow to node "i".
                    for (unsigned int phase=0;phase<nodeNofPhases[pair.first];phase++) {
                        double phaseMin = refinedNeuronLimitBounds[nodeConnectionIn[pair.first][phase].first].first;
                        double phaseMax = refinedNeuronLimitBounds[nodeConnectionIn[pair.first][phase].first].second;
                        bool isBadPhase = false;

                        // The flow would be too high?
                        if (valueFixedMax-pair.second*(phaseMax-phaseMin)+EPSILON_FOR_CALLING_AN_APPROXIMATED_NODE_VALUE_TO_BE_PRECISE_RELATIVE_TO_BOUNDS<howMuchInflowIsAdmissible[i].first) {
                            isBadPhase = true;
                        }

                        // The flow would be too low?
                        if (valueFixedMin+pair.second*(phaseMax-phaseMin)-EPSILON_FOR_CALLING_AN_APPROXIMATED_NODE_VALUE_TO_BE_PRECISE_RELATIVE_TO_BOUNDS>howMuchInflowIsAdmissible[i].second) {
                            isBadPhase = true;
                        }

                        // Something to propagate?
                        if (isBadPhase
                                && (currentAssignmentAsSet.count(-1*(startingSATVarsPhases[pair.first])-phase)==0)
                                && (currentAssignmentAsSet.count(startingSATVarsPhases[pair.first]+phase)==0)) {

                            //std::cerr << "MaxPool Phase Off switcher: " << nodeNames[pair.first] << "," << phase << "," << phaseMin << "," << phaseMax << "," << valueFixedMax << "," << pair.second << "," << howMuchInflowIsAdmissible[i].first << "," << howMuchInflowIsAdmissible[i].second << " due to " << nodeNames[i] << std::endl;
                            std::set<int> sortedSetReasons(reasons[i].begin(),reasons[i].end());
                            std::vector<int> newClause;
                            newClause.push_back(-1*startingSATVarsPhases[pair.first]-phase);
                            newClause.insert(newClause.end(),sortedSetReasons.begin(),sortedSetReasons.end());
                            clausesToAdd.push_back(newClause);
                        }
                    }
                } else {
                        throw "Unsupported node type....715";
                }
            }
        }
    }

    // Debugging output

    /*std::cerr << "- Propagation -- new clauses to propagate:\n";
    for (auto clause : clausesToAdd) {
        std::cerr << " -";
        for (auto it : clause) std::cerr << " " << it;
        std::cerr << std::endl;
    }*/


    return true;
}


/**
 * @brief This is the call-back function for Minisat -- Whenever a new partial valuation is found that is not in conflict with any existing clauses, Minisat calls this function
 *        to check the LP relaxation.
 * @param partialValuation The current partial assignment. This is a vector or -1, 0,
 * @param addedClauses
 * @return
 */

bool VerificationProblem::checkPartialNodeFixtureInLPRelaxation(std::vector<int> varFixture, int nofLevel0VarsInFixture, std::list<std::vector<int> > &clausesToAdd) {

    // Step 1: Ask the "Superset database" if the partial fixture is fine.
    //std::cerr << "Unfixed vars: " << varFixture.size() << "," << nofSATVars << std::endl;

    if (varFixture.size()<nofSATVars) {
        unsigned int index = 0;
        while ((index<varFixture.size())
               && (index<lastPartialVariableValuationGivenToCheckPartialNodeFixtureInLPRelaxation.size())
               && (varFixture[index]==lastPartialVariableValuationGivenToCheckPartialNodeFixtureInLPRelaxation[index])) index++;
        while (lastPartialVariableValuationGivenToCheckPartialNodeFixtureInLPRelaxation.size()>index) {
            if (lastPartialVariableValuationGivenToCheckPartialNodeFixtureInLPRelaxation.back()>0) {
                int data = cachedOptima.removeElementFromTrace();
                (void)data;
                assert(data==lastPartialVariableValuationGivenToCheckPartialNodeFixtureInLPRelaxation.back());
                lastPartialVariableValuationGivenToCheckPartialNodeFixtureInLPRelaxation.pop_back();
            } else {
                lastPartialVariableValuationGivenToCheckPartialNodeFixtureInLPRelaxation.pop_back();
            }
        }
        bool found = true;
        while (index<varFixture.size()) {
            lastPartialVariableValuationGivenToCheckPartialNodeFixtureInLPRelaxation.push_back(varFixture[index]);
            if (varFixture[index]>0) {
                found = cachedOptima.addTraceElement(varFixture[index]);
            }
            index++;
        }
        if (found) return true; // We found a superset -- so no need to check!
    }

    // Print
    /*std::cerr << "Current partial valuation to checkPartialSATValuation: ";
    for (unsigned int i=0;i<varFixture.size();i++) {
            std::cerr << varFixture[i];
            if (i>=varFixture.size()) std::cerr << "*";
            std::cerr << " ";
    }
    std::cerr << "\n";*/

    // Remove negative varFixture elements that are directly implied by positive ones (due to a phase having been selected)
    std::set<int> litsGreaterThanLevel0(varFixture.begin()+nofLevel0VarsInFixture,varFixture.end());
    std::set<int> allLits(varFixture.begin(),varFixture.end());

    {
        std::list<int> newVarFixture;
        for (unsigned int i=nofLevel0VarsInFixture;i<varFixture.size();i++) {
            if (varFixture[i]<0) {
                std::pair<int,int> affectedNode = satVarToNodeAndPhaseMapper[-1*varFixture[i]];
                if (nodeTypes[affectedNode.first]==RELU) {
                    assert((litsGreaterThanLevel0.count(startingSATVarsPhases[affectedNode.first])>0) || (litsGreaterThanLevel0.count(startingSATVarsPhases[affectedNode.first]+1)>0));
                } else if (nodeTypes[affectedNode.first]==MAXPOOL) {
                    // Add anyway -- we may make use of this one!
                    newVarFixture.push_back(varFixture[i]);
                } else {
                    throw "Unsupported (1209)";
                }
            } else {
                newVarFixture.push_back(varFixture[i]);
            }
        }
        varFixture.resize(nofLevel0VarsInFixture);
        varFixture.insert(varFixture.end(),newVarFixture.begin(),newVarFixture.end());
    }


    // Build a new LP instance and check if it is satisfiable
    if (literalsAlreadyFixed.size()==0) {
        lpPartialFixture = LPProblem(nodeTypes.size());
        addConstraintsToLPInstance(lpPartialFixture,nodeTypes.size());
        addNodePropagationInequationsToLPInstance(lpPartialFixture,nodeTypes.size());
        imposeNeuronValueBoundsInLPInstance(lpPartialFixture);
    }

    // Build fixed optimization function, consisting of all phased nodes not yet fixed
    LPREALTYPE fixedOptimizationFunctionNetworkSolving[nodeNames.size()+1];
    fixedOptimizationFunctionNetworkSolving[0] = 0.0;
    for (unsigned int i=0;i<nodeTypes.size();i++) {
        if ((nodeTypes[i]==INPUT) || (nodeTypes[i]==LINEAR)) {
            fixedOptimizationFunctionNetworkSolving[i+1] = 0.0;
        } else if (nodeTypes[i]==RELU) {
            fixedOptimizationFunctionNetworkSolving[i+1] = 1.0;
        } else if (nodeTypes[i]==MAXPOOL) {
            // Optimize a little bit.
            fixedOptimizationFunctionNetworkSolving[i+1] = 0.01;
        } else {
            throw "Unsupported node type. (992)";
        }
    }
    for (auto it : varFixture) {
        int thisVar = std::abs(it);
        // For now, only RELUs and MAXPOOLs are supported
        assert((nodeTypes[satVarToNodeAndPhaseMapper[thisVar].first]==RELU) || (nodeTypes[satVarToNodeAndPhaseMapper[thisVar].first]==MAXPOOL));
        fixedOptimizationFunctionNetworkSolving[satVarToNodeAndPhaseMapper[thisVar].first+1] = 0.0;
    }

    // Set optimization function
    lpPartialFixture.setObjFn(fixedOptimizationFunctionNetworkSolving);
    lpPartialFixture.setMinim();

    // Fix Phases
    unsigned int nofUnfixedPhases = nofNodesWithPhases;

    std::map<int,std::set<int> > maxPoolFixtures; // For buffering what phases have already been fixed.
    std::map<int,std::set<int> > maxPoolFixturesLevel0; // For buffering what phases have already been fixed.

    for (unsigned int i=0;i<varFixture.size();i++) {

        int literal = varFixture[i];

        // Currently, only look at positive fixtures
        if (literal>0) {

            std::pair<int,int> affectedNode = satVarToNodeAndPhaseMapper[literal];
            if (nodeTypes[affectedNode.first]==RELU) {
                if (literalsAlreadyFixed.count(literal)==0)
                    addReLUFixture(lpPartialFixture,affectedNode.first,affectedNode.second+1,nodeTypes.size()); // Phase numbering starts with "1"
                nofUnfixedPhases -= 1;
                literalsAlreadyFixed.insert(literal);
            } else if (nodeTypes[affectedNode.first]==MAXPOOL) {
                maxPoolFixtures[affectedNode.first].insert(affectedNode.second+1);
                if (i<(unsigned int)nofLevel0VarsInFixture)
                    maxPoolFixturesLevel0[affectedNode.first].insert(affectedNode.second+1);
            } else {
                std::cerr << "Literal: " << literal << std::endl;
                std::cerr << "Node: " << nodeNames[affectedNode.first] << std::endl;
                throw "Unsupported node types (625March).";
            }

        } else {
            std::pair<int,int> affectedNode = satVarToNodeAndPhaseMapper[-1*literal];
            if (nodeTypes[affectedNode.first]==MAXPOOL) {
                maxPoolFixtures[affectedNode.first].insert(-1*affectedNode.second-1);
                if (i<(unsigned int)nofLevel0VarsInFixture)
                    maxPoolFixturesLevel0[affectedNode.first].insert(-1*affectedNode.second-1);
            }
        }
    }

    // MaxPoolFixture
    for (auto it : maxPoolFixtures) {
        bool allAreAlreadyFixed = true;
        // std::cerr << "PoolFixtures " << nodeNames[it.first] << ": " << it.second.size() << "->" << startingSATVarsPhases[it.first] << std::endl;
        for (auto it2 : it.second) {
            if (it2<0) {
                if (literalsAlreadyFixed.count(-1*(startingSATVarsPhases[it.first]-it2-1))==0) allAreAlreadyFixed = false;
            } else {
                if (literalsAlreadyFixed.count(startingSATVarsPhases[it.first]+it2-1)==0) allAreAlreadyFixed = false;
            }

        }
        if (!allAreAlreadyFixed) {
            addMaxPoolFixture(lpPartialFixture,it.first,it.second,nodeTypes.size());
            for (auto it2 : it.second) {
                if (it2<0) {
                    literalsAlreadyFixed.insert(-1*(startingSATVarsPhases[it.first]-it2-1));
                } else {
                    literalsAlreadyFixed.insert(startingSATVarsPhases[it.first]+it2-1);
                }
            }
        }
        if (it.second.size()>=nodeNofPhases[it.first]-1) {
            nofUnfixedPhases -= 1;
        } else {
            for (auto it2 : it.second) {
                if (it2>0) {
                    nofUnfixedPhases -= 1;
                }
            }
        }
    }

    //std::cerr << std::endl;

    int result = lpPartialFixture.solve();
    if (result==LPINFEASIBLE) {
        std::cerr << "LPINFEASIBLE\n";

        // Apply the method from the paper "Locating Minimal Infeasible Constraint Sets in Linear Programs."

        // Trvial return
        //std::vector<int> backClause;
        //for (auto lit : varFixture) backClause.push_back(-1*lit);
        //clausesToAdd.push_back(backClause);
        //return false;


        // Prepare cached optima lookup
        //lastPartialVariableValuationGivenToCheckPartialNodeFixtureInLPRelaxation.clear();

        //std::vector<bool> literalsToKeepInConflictClause;
        //lp.setObjFn(nullOptimizationFunctionForLPSolving);

        // Level-0-variables are not taken inco account...
        //for (int i=0;i<nofLevel0VarsInFixture;i++) {
        //    literalsToKeepInConflictClause.push_back(false);
        //}

        // Allocate additional RELU slack variables for minimal infeasible set finding.
        // All variables except for the level-0 ones and the very last one can be unneeded



        int lastPositiveRELULiteral = startingSATVarsPhases.at(satVarToNodeAndPhaseMapper.at(std::abs(varFixture[varFixture.size()-1])).first);
        //std::cerr << "RELU Last Positive RELU Literal: " << lastPositiveRELULiteral << std::endl;
        int nextLPVar = nodeTypes.size();
        std::map<unsigned int /* node number */, std::pair<int /* lp var */, int /* phase */> > reluslackVars;
        for (unsigned int i=nofLevel0VarsInFixture;i<varFixture.size();i++) {
            if (varFixture[i]>0) {
                auto mappedNodeAndPhase = satVarToNodeAndPhaseMapper[varFixture[i]];
                if ((nodeTypes[mappedNodeAndPhase.first]==RELU) && (lastPositiveRELULiteral!=varFixture[i]) && (lastPositiveRELULiteral+1!=varFixture[i])) {
                    reluslackVars[mappedNodeAndPhase.first] = std::pair<int,int>(nextLPVar++,mappedNodeAndPhase.second);
                    //std::cerr << "Added RELU Slack Var for Node " << mappedNodeAndPhase.first << "," << nodeNames[mappedNodeAndPhase.first] << std::endl;
                } else {
                    //std::cerr << "Did not add RELU Slack Node for Lit: " << varFixture[i] << std::endl;
                }
            }
        }

        // Allocate MAXPOOL slack vars. If the node is not fully fixed, we have one variable per phase that is
        // already ruled out. Otherwise, we have one variable for the fixture.
        std::map<unsigned int /* node number */, int /* lp var */> maxpoolSlackVars;
        for (auto it : maxPoolFixtures) {
            unsigned int node = it.first;
            maxpoolSlackVars[node] = nextLPVar++;
        }


        // Prepare optimization function
        LPREALTYPE fixedOptimizationFunctionNetworkSolving[nextLPVar];
        for (unsigned int i=0;i<=(unsigned int)nextLPVar;i++) {
            fixedOptimizationFunctionNetworkSolving[i] = (i<=nodeTypes.size()?0.0:1.0);
        }
        // Let's try out if giing lower weight to the later literals in the trace helps
        for (unsigned int i=nofLevel0VarsInFixture;i<varFixture.size();i++) {
            if (varFixture[i]>0) {
                auto mappedNodeAndPhase = satVarToNodeAndPhaseMapper[varFixture[i]];
                if (nodeTypes[mappedNodeAndPhase.first]==RELU) {
                    if ((lastPositiveRELULiteral!=varFixture[i]) && (lastPositiveRELULiteral+1!=varFixture[i]))
                        fixedOptimizationFunctionNetworkSolving[reluslackVars.at(mappedNodeAndPhase.first).first+1] = varFixture.size()-i+1;
                } else if (nodeTypes[mappedNodeAndPhase.first]==MAXPOOL) {
                    fixedOptimizationFunctionNetworkSolving[maxpoolSlackVars.at(mappedNodeAndPhase.first)+1] = varFixture.size()-i+1;
                } else {
                    throw 12345;
                }
            } else {
                auto mappedNodeAndPhase = satVarToNodeAndPhaseMapper[-1*varFixture[i]];
                if (nodeTypes[mappedNodeAndPhase.first]==MAXPOOL) {
                    if (maxpoolSlackVars.count(varFixture[i])>0) {
                        fixedOptimizationFunctionNetworkSolving[maxpoolSlackVars.at(mappedNodeAndPhase.first)+1] = varFixture.size()-i+1;
                    }
                }

            }
        }

        LPProblem lpInfeasible(nextLPVar);
        addConstraintsToLPInstance(lpInfeasible,nextLPVar);
        addNodePropagationInequationsToLPInstance(lpInfeasible,nextLPVar);
        imposeNeuronValueBoundsInLPInstance(lpInfeasible);
        lpInfeasible.setObjFn(fixedOptimizationFunctionNetworkSolving);
        lpInfeasible.setMinim();

        // Fix bounds of the slack variables
        for (auto it : reluslackVars) {
            lpInfeasible.setBoundVariableOnlyFromBelow(it.second.first,0.0);
        }
        for (auto it : maxpoolSlackVars) {
            lpInfeasible.setBoundVariableOnlyFromBelow(it.second,0.0);
        }

        // Now fix the RELU fixed nodes (modulo slack)
        for (unsigned int i=0;i<varFixture.size();i++) {
            int thisLiteral = varFixture[i];
            if (thisLiteral>0) {
                auto mappedNodeAndPhase = satVarToNodeAndPhaseMapper[thisLiteral];
                if ((nodeTypes[mappedNodeAndPhase.first]==INPUT) || (nodeTypes[mappedNodeAndPhase.first]==LINEAR)) {
                    throw "Should not happen.";
                } else if (nodeTypes[mappedNodeAndPhase.first]==RELU) {

                    // Is the node "open for discussion" in the minimal infeasible subset finding part?
                    if (reluslackVars.count(mappedNodeAndPhase.first)==0) {
                        //std::cerr << "(fix" << thisLiteral << ")";
                        addReLUFixture(lpInfeasible, mappedNodeAndPhase.first, mappedNodeAndPhase.second+1,nextLPVar); // Phase numbering starts with "1".
                    } else {
                        //std::cerr << "(slack)";
                        addReLUFixtureWithSlackVar(lpInfeasible, mappedNodeAndPhase.first, reluslackVars[mappedNodeAndPhase.first].second+1, reluslackVars[mappedNodeAndPhase.first].first,nextLPVar); // Phase numbering starts with "1".
                    }
                } else if (nodeTypes[mappedNodeAndPhase.first]==MAXPOOL) {
                    // Is fixed elsewhere
                } else {
                    throw "Unsupported node type 979";
                }
            }
        }

        // Now fix the Slack nodes
        for (auto it : maxPoolFixtures) {
            //std::cerr << "(slackPool " << it.first << ")";
            addMaxPoolFixtureWithSlackVar(lpInfeasible,it.first,it.second, maxpoolSlackVars[it.first],nextLPVar);
        }

        for (auto it : maxPoolFixturesLevel0) {
            //std::cerr << "(slackPool " << it.first << ")";
            addMaxPoolFixture(lpInfeasible,it.first,it.second,nextLPVar);
        }

        // But all fixtures that result from the last literal (that we hard-add to the resulting clause) bust be *tight*:
        if (nodeTypes[satVarToNodeAndPhaseMapper[std::abs(varFixture[varFixture.size()-1])].first]==MAXPOOL) {
            //std::cerr << "Pooling with:" << varFixture[varFixture.size()-1] << "\n";
            std::set<int> fix;
            if (varFixture[varFixture.size()-1]<0) {
                fix.insert(-1*satVarToNodeAndPhaseMapper[std::abs(varFixture[varFixture.size()-1])].second-1);
                //std::cerr << "Phasse: " << satVarToNodeAndPhaseMapper[std::abs(varFixture[varFixture.size()-1])].second+1 << std::endl;
            } else
                fix.insert(satVarToNodeAndPhaseMapper[std::abs(varFixture[varFixture.size()-1])].second+1);
            addMaxPoolFixture(lpInfeasible,satVarToNodeAndPhaseMapper[std::abs(varFixture[varFixture.size()-1])].first,fix,nextLPVar);
        }

        std::vector<int> resultingClause;
        resultingClause.push_back(-1*varFixture[varFixture.size()-1]);
        std::set<int> nodesAlreadyFixed;
        //std::cerr << "Enter loop";
        while (true) {
            int result2 = lpInfeasible.solve();
            //std::cerr << "Loop res: " << result2 << std::endl;
            if (result2==LPINFEASIBLE) {
                //std::cerr << "Adding clause starting with " << resultingClause[0] << std::endl;
                clausesToAdd.push_back(resultingClause);
                /*if ((resultingClause[0]==-4) && (resultingClause.size()==1)) {
                    lpPartialFixture.writeLP("/tmp/favorite.lp");
                    throw 333;
                }*/

#ifndef NDEBUG
                // Test the clause just added for sanity
                std::cerr << "Check: ";
                for (auto it : resultingClause) {
                    std::cerr << it << " ";
                    assert(allLits.count(-1*it)>0);
                }
                std::cerr << std::endl;
                std::cerr << "#c2: " << clausesToAdd.size() << std::endl;
#endif

                return false;
            } else {
                int bestNode = -1;
                std::vector<int> candidateLits;
                double bestNodeValue = 0.0;
                for (auto it : reluslackVars) {
                    if (nodesAlreadyFixed.count(it.first)==0) {
                        double solutionVal = lpInfeasible.getSolutionVarValue(it.second.first);
                        //std::cerr << "Slack: " << solutionVal;
                        if (solutionVal > bestNodeValue) {
                            bestNodeValue = solutionVal;
                            bestNode = it.first;
                            candidateLits.clear();
                            //std::cerr << "(candidateLit: " << -1*startingSATVarsPhases[it.first]-it.second.second << ")";
                            candidateLits.push_back(-1*startingSATVarsPhases[it.first]-it.second.second);
                        }
                    }
                }
                for (auto it : maxpoolSlackVars) {
                    if (nodesAlreadyFixed.count(it.first)==0) {
                        double solutionVal = lpInfeasible.getSolutionVarValue(it.second);
                        // std::cerr << "SlackPool: " << solutionVal;
                        if (solutionVal > bestNodeValue) {
                            bestNodeValue = solutionVal;
                            bestNode = it.first;
                            candidateLits.clear();
                            if (maxPoolFixtures[it.first].size()==nodeNofPhases[it.first]) {
                                // Push only positive lits
                                //std::cerr << "Positive!\n";
                                for (auto it2 : maxPoolFixtures[it.first]) {
                                    if (it2>0) {
                                        //std::cerr << "(PUSH" << -1*(startingSATVarsPhases[it.first]+it2-1) << ")";
                                        candidateLits.push_back(-1*(startingSATVarsPhases[it.first]+it2-1));
                                    }
                                }
                            } else {
                                // If there is no positive literal, then all others are the candidate literals
                                //std::cerr << "Negative: ";
                                //for (auto it2 : litsGreaterThanLevel0) std::cerr << it2 << " ";
                                //std::cerr << std::endl;
                                for (auto it2 : maxPoolFixtures[it.first]) {
                                    candidateLits.push_back(startingSATVarsPhases[it.first]-it2-1);
                                    assert(it2<0);
                                    //std::cerr << "(push" << it2 << "," << startingSATVarsPhases[it.first]-it2-1 << ")";
                                }
                            }
                        }
                    }
                }

                if (bestNodeValue>0) {
                    assert(bestNode>-1);
                    resultingClause.insert(resultingClause.end(),candidateLits.begin(),candidateLits.end());
                    nodesAlreadyFixed.insert(bestNode);
                } else {
                    // We already found an optimum -- then the last fixed literal suffices. Just return in this case.
                    clausesToAdd.push_back(resultingClause);
#ifndef NDEBUG
                    // Test the clause just added for sanity
                    std::cerr << "Check: ";
                    for (auto it : resultingClause) {
                        std::cerr << it << " ";
                        assert(allLits.count(-1*it)>0);
                    }
                    std::cerr << std::endl;
                    std::cerr << "#c1617: " << clausesToAdd.size() << std::endl;
#endif

                    return false;
                }
                for (auto it : reluslackVars) {
                    if (static_cast<unsigned int>(bestNode)==it.first) {
                        lpInfeasible.setBounds(it.second.first,0.0,0.0);
                    }
                }
                for (auto it : maxpoolSlackVars) {
                    if (static_cast<unsigned int>(bestNode)==it.first) {
                        lpInfeasible.setBounds(it.second,0.0,0.0);
                    }
                }
            }
        }


    } else if (result==LPOPTIMUM) {

        std::set<int> lits(varFixture.begin(),varFixture.end());

        //std::cerr << "LPOPTIMUM " << nofUnfixedPhases << "\n";
        // Ok, is fine.
        if (nofUnfixedPhases==0) {

            lpPartialFixture.writeLP("/tmp/finalLP.lp");

            // Compute Error margin
            double values[nodeNames.size()];
            double error = 0;
            for (unsigned int i=0;i<nodeTypes.size();i++) {
                if (nodeTypes[i]==INPUT) {
                    values[i] = lpPartialFixture.getSolutionVarValue(i);
                } else if ((nodeTypes[i]==RELU) || (nodeTypes[i]==LINEAR)) {
                    double sum=0;
                    for (auto it : nodeConnectionIn[i]) {
                        sum += it.second*values[it.first];
                    }
                    sum += nodeBias[i];
                    if (nodeTypes[i]==RELU) sum = std::max(0.0,sum);
                    values[i] = sum;

                } else if (nodeTypes[i]==MAXPOOL) {
                    double max = -1*std::numeric_limits<double>::max();
                    for (auto it : nodeConnectionIn[i]) {
                        max = std::max(max,values[it.first]);
                    }
                    values[i] = max;
                } else {
                    throw "Unsupported node type.";
                }
                //std::cerr << "DIFF for value " << i << ": " << values[i] - lp.getSolutionVarValue(i) << std::endl;

                error += std::abs(values[i] - lpPartialFixture.getSolutionVarValue(i));
            }

            // SAT -- Compute model
            std::cout << "SAT\n\nValuation:\n";
            for (unsigned int i=0;i<nodeTypes.size();i++) {
                std::cout << "- " << nodeNames[i] << ": " << lpPartialFixture.getSolutionVarValue(i) << " / " << values[i] << std::endl;
            }
            lpPartialFixture.writeLP("/tmp/final.lp");

            // Write SAT literals
            std::cout << "Literals:";
            for (auto it : lits) std::cout << " " << it;
            std::cout << std::endl;

            std::cout << "Total error: " << error << std::endl;

            return true;

        } else {

            // So this partial valuation is fine. See if we find some new propagated (disjunction of) values
            // --- abort when we don't see anything more.
            bool foundUnfixedMaxPoolNode = false;
            //bool numericallyUnstable = false;

            std::vector<double> realNetValuation(nodeTypes.size());
            std::vector<int> newClause;
            bool foundGreater0Node = false;
            std::vector<int> newCachedOptimaEntry;
            for (unsigned int i=0;i<nodeTypes.size();i++) {
                realNetValuation[i] = lpPartialFixture.getSolutionVarValue(i);
                if (nodeTypes[i]==INPUT) {
                    // Nothing to do.
                } else if (nodeTypes[i]==LINEAR) {
                    // Nothing to do.
                } else if (nodeTypes[i]==RELU) {
                    if ((lits.count(startingSATVarsPhases[i]+1)==0) && (lits.count(startingSATVarsPhases[i])==0)) {

                        newClause.push_back(startingSATVarsPhases[i]+1);

                        // The following cannot be >0, as then a value between 0 and epsilon value would be detected as a necessary consequence of the partial valuation, which may be wrong on the Boolean level.
                        if (lpPartialFixture.getSolutionVarValue(i)>0) {

                            foundGreater0Node = true;

                            // If this one is smaller than EPSILON_FOR_CALLING_AN_APPROXIMATED_NODE_VALUE_TO_BE_PRECISE_RELATIVE_TO_BOUNDS) { //EPSILON_FOR_CALLING_AN_APPROXIMATED_NODE_VALUE_TO_BE_PRECISE_RELATIVE_TO_BOUNDS,
                            // then we need to be on the safe side and not learn anything, as we may run into numerical problems
                            if (lpPartialFixture.getSolutionVarValue(i)<EPSILON_FOR_CALLING_AN_APPROXIMATED_NODE_VALUE_TO_BE_PRECISE_RELATIVE_TO_BOUNDS) return true; // Numerically unstable -- can't make use of the result

                            // Check if the value of this node is "tight".
                            double sumIn = nodeBias[i];
                            for (auto it : nodeConnectionIn[i]) {
                                sumIn += it.second*realNetValuation[it.first];
                            }
                            bool isTight = (std::abs(realNetValuation[i]-sumIn)<EPSILON_FOR_CALLING_AN_APPROXIMATED_NODE_VALUE_TO_BE_PRECISE_RELATIVE_TO_BOUNDS);

                            if (isTight)
                                newCachedOptimaEntry.push_back(startingSATVarsPhases[i]+1);
                        } else {
                            // If this one is "0", then mark this in the cacheOptimaEntry
                            newCachedOptimaEntry.push_back(startingSATVarsPhases[i]);
                            //newCachedOptimaEntry.push_back(startingSATVarsPhases[i]+1); <--- IS THIS SAFE?
                        }
                    }
                } else if (nodeTypes[i]==MAXPOOL) {
                    // Ignoring the MAXPOOL layers does not lead to errors as
                    // the fixed maxpool literals will still be added later to the clause
                    bool foundAFixture = false;
                    for (unsigned int j=0;j<nodeNofPhases[j];j++) {
                        if (lits.count(startingSATVarsPhases[i]+j)>0) {
                            foundAFixture = true;
                        }
                    }

                    bool foundTooLarge = false;
                    for (unsigned int j=0;j<nodeNofPhases[j];j++) {
                        if (lits.count(-1*startingSATVarsPhases[i]-j)==0) {
                            if (lpPartialFixture.getSolutionVarValue(i)>lpPartialFixture.getSolutionVarValue(nodeConnectionIn[i][j].first)) foundTooLarge = true;
                        }
                    }

                    foundUnfixedMaxPoolNode |= !foundAFixture | !foundTooLarge;
                } else {
                    throw "NodeType not supported. 981";
                }
            }

            // If we found any phase>0, then add a clause to "propagatedClauses".
            // Otherwise, we should be done with the search process.
            if (foundGreater0Node) {
                for (auto it : litsGreaterThanLevel0) newClause.push_back(-1*it);
                clausesToAdd.push_back(newClause);
            }

            for (auto it : lits) {
                if (it>0) newCachedOptimaEntry.push_back(it);
            }
            cachedOptima.addSet(newCachedOptimaEntry);

            if (!foundGreater0Node && !foundUnfixedMaxPoolNode) {
                // If every open variable had a minimal value, then we actually found a solution and only need to prod the SAT solver towards it.
                // ...but only if all MaxPOOL variables have values!
                //
                // Check if any MaxPool Node is not fully fixed --
                // 'cause then we can't propagate the solution
                for (auto it : maxPoolFixtures) {
                    bool foundPos = false;
                    for (auto it2 : it.second) {
                        if (it2>0) foundPos = true;
                    }
                    if (!foundPos) return true; // ok, no propagation in this case.
                }

                for (unsigned int i=0;i<nodeTypes.size();i++) {
                    if (nodeTypes[i]==INPUT) {
                        // Nothing to do.
                    } else if (nodeTypes[i]==LINEAR) {
                        // Nothing to do.
                    } else if (nodeTypes[i]==RELU) {
                        if ((lits.count(startingSATVarsPhases[i]+1)==0) && (lits.count(startingSATVarsPhases[i])==0)) {
                            newClause.clear();
                            if (realNetValuation[i]>EPSILON_FOR_CALLING_AN_APPROXIMATED_NODE_VALUE_TO_BE_PRECISE_RELATIVE_TO_BOUNDS) {
                                newClause.push_back(startingSATVarsPhases[i]+1);
                                for (auto it : litsGreaterThanLevel0) newClause.push_back(-1*it);
                                clausesToAdd.push_back(newClause);
                            } else if (realNetValuation[i]==0.0) {
                                newClause.push_back(startingSATVarsPhases[i]);
                                for (auto it : litsGreaterThanLevel0) newClause.push_back(-1*it);
                                clausesToAdd.push_back(newClause);
                            }
                            //std::cerr << "_CL_ADD_" << newClause[0] << "__LINE" << __LINE__ << std::endl;
                        }
                    } else if (nodeTypes[i]==MAXPOOL) {
                        // MaxPOOLs can be ignored here.
                    } else {
                        throw "NodeType not supported. 1016";
                    }
                }
            }
        }

        return clausesToAdd.size()==0;
    } else {
        throw "Unexpected LP Solution Return Code";
    }

}

std::string doubleToSMTLib2Expression(double in) {
    std::ostringstream ostream;
    ostream.precision(std::numeric_limits< double >::max_digits10);
    ostream << std::fixed << in;
    std::string a = "";
    std::string b = "1";
    std::string from = ostream.str();
    bool foundDot = false;
    bool minus = false;
    for (unsigned int i=0;i<from.size();i++) {
        char tc = from[i];
        if (tc=='.') {
            foundDot = true;
        } else if (tc=='-') {
          minus = true;
        } else {
            if ((a!="") || (tc!='0')) {
                a = a + tc;
            }
            if (foundDot) b = b + "0";
        }
    }
    if (a.size()==0) return "0";
    if (minus) {
        return "(- 0 (/ "+a+" "+b+"))";
    } else {
        return "(/ "+a+" "+b+")";
    }
}


void VerificationProblem::printILPInstance(bool includingApproximationConstraints) {

    std::cout.precision(std::numeric_limits< double >::max_digits10);
    std::cout << std::scientific;

    std::cout << "Minimize\n";
    std::cout << "  " << nodeNames[0] << std::endl;
    std::cout << "Subject To\n";
    for (unsigned int i=0;i<nodeNames.size();i++) {
        if (nodeTypes[i]==LINEAR) {
            std::cout << "  ";
            bool first = false;
            for (auto it : nodeConnectionIn[i]) {
                if (!first) std::cout << "+ ";
                std::cout << it.second << " " << nodeNames[it.first] << " ";
            }
            std::cout << "- " << nodeNames[i] << " = " << -1*nodeBias[i] << std::endl;
        } else if (nodeTypes[i]==RELU) {

            // Phase 0
            std::cout << "  "+nodeNames[i]+"phase = 0 -> ";
            bool first = true;
            for (auto it : nodeConnectionIn[i]) {
                if (!first) std::cout << "+ ";
                std::cout << it.second << " " << nodeNames[it.first] << " ";
                first = false;
            }
            std::cout << "<= " << -1*nodeBias[i] << std::endl;
            std::cout << "  "+nodeNames[i]+"phase = 0 -> "+nodeNames[i]+" = 0\n";

            // Phase 1
            std::cout << "  "+nodeNames[i]+"phase = 1 -> ";
            first = true;
            for (auto it : nodeConnectionIn[i]) {
                if (!first) std::cout << "+ ";
                std::cout << it.second << " " << nodeNames[it.first] << " ";
                first = false;
            }
            std::cout << "- " << nodeNames[i] << " = " << -1*nodeBias[i] << std::endl;
            std::cout << "  "+nodeNames[i]+"phase = 1 -> "+nodeNames[i]+" >= 0\n";
        } else if (nodeTypes[i]==MAXPOOL) {

            // Maxpool
            unsigned int phase = 0;
            for (auto it : nodeConnectionIn[i]) {
                std::cout << "  " << nodeNames[i] << "phase" << phase << " = 1 -> " << nodeNames[i] << " - " << nodeNames[it.first] << " = 0\n";
                for (auto it2 : nodeConnectionIn[i]) {
                    if (it.first!=it2.first) {
                        std::cout << "  " << nodeNames[i] << "phase" << phase << " = 1 -> " << nodeNames[i] << " - " << nodeNames[it2.first] << " >= 0\n";
                    }
                }
                phase++;
            }
        }
    }

    // Phase selection for MaxPool
    for (unsigned int i=0;i<nodeNames.size();i++) {
        if (nodeTypes[i]==MAXPOOL) {
            std::cout << "  ";
            for (unsigned int j=0;j<nodeConnectionIn[i].size();j++) {
                if (j>0) std::cout << " + ";
                std::cout << nodeNames[i] << "phase" << j;
            }
            std::cout << " = 1\n";
        }
    }

    if (includingApproximationConstraints) {
        for (unsigned int i=0;i<nodeTypes.size();i++) {
            double lower = initialNeuronLimitBounds[i].first;
            double upper = initialNeuronLimitBounds[i].second;
            if (lower<0) lower = lower*(1.0+EPSILON_FOR_CALLING_AN_APPROXIMATED_NODE_VALUE_TO_BE_PRECISE_RELATIVE_TO_BOUNDS)-EPSILON_FOR_CALLING_AN_APPROXIMATED_NODE_VALUE_TO_BE_PRECISE_RELATIVE_TO_BOUNDS;
            if (lower>0) lower = lower*(1.0-EPSILON_FOR_CALLING_AN_APPROXIMATED_NODE_VALUE_TO_BE_PRECISE_RELATIVE_TO_BOUNDS)-EPSILON_FOR_CALLING_AN_APPROXIMATED_NODE_VALUE_TO_BE_PRECISE_RELATIVE_TO_BOUNDS;
            if (upper<0) upper = upper*(1.0-EPSILON_FOR_CALLING_AN_APPROXIMATED_NODE_VALUE_TO_BE_PRECISE_RELATIVE_TO_BOUNDS)+EPSILON_FOR_CALLING_AN_APPROXIMATED_NODE_VALUE_TO_BE_PRECISE_RELATIVE_TO_BOUNDS;
            if (upper>0) upper = upper*(1.0+EPSILON_FOR_CALLING_AN_APPROXIMATED_NODE_VALUE_TO_BE_PRECISE_RELATIVE_TO_BOUNDS)+EPSILON_FOR_CALLING_AN_APPROXIMATED_NODE_VALUE_TO_BE_PRECISE_RELATIVE_TO_BOUNDS;
            std::cout << "  " << nodeNames[i] << " >= " << lower << std::endl;
            std::cout << "  " << nodeNames[i] << " <= " << upper << std::endl;
        }


        for (unsigned int i=0;i<nodeTypes.size();i++) {

            if ((nodeTypes[i]==INPUT) || (nodeTypes[i]==LINEAR)) {
                // Nothing to do, as there are no additional constraints
            } else if (nodeTypes[i]==RELU) {

                // Add "AT LEAST" constraint
                bool first = true;
                std::cout << "  ";
                for (auto it : nodeConnectionIn[i]) {
                    if (!first) std::cout << "+ ";
                    std::cout << it.second << " " << nodeNames[it.first] << " ";
                    first = false;
                }
                std::cout << "- " << nodeNames[i] << " <= " << -1*nodeBias[i] << std::endl;

                // Add ">=0" constraint
                std::cout << "  " << nodeNames[i] << " >= 0\n";

                // Compute theoretical minimum
                double min = nodeBias[i];
                double max = nodeBias[i];
                for (auto pair : nodeConnectionIn[i]) {
                    if (pair.second<0.0) {
                        min += pair.second*initialNeuronLimitBounds[pair.first].second;
                        max += pair.second*initialNeuronLimitBounds[pair.first].first;
                    } else {
                        min += pair.second*initialNeuronLimitBounds[pair.first].first;
                        max += pair.second*initialNeuronLimitBounds[pair.first].second;
                    }
                }
                max = std::min(initialNeuronLimitBounds[i].second,max);

                if (max<=0) {
                    std::cout << "  " << nodeNames[i] << " = 0\n";
                } else {
                    if (min<0) {
                        // Compute line between (min,0) and (max,max)
                        double factor = max/(max-min);
                        double offset = max-max*factor;


                        // Add ">=0" constraint
                        std::cout << "  " << nodeNames[i];
                        for (auto pair : nodeConnectionIn[i]) {
                            std::cout << " - " << 1*factor*pair.second << " " << nodeNames[pair.first];
                        }
                        std::cout << " <= " << nodeBias[i]*factor+offset << "\n";

                    } else {
                        // Make it linear
                        first = true;
                        for (auto it : nodeConnectionIn[i]) {
                            if (!first) std::cout << "+ ";
                            std::cout << it.second << " " << nodeNames[it.first] << " ";
                            first = false;
                        }
                        std::cout << "- " << nodeNames[i] << " = " << -1*nodeBias[i] << std::endl;
                    }
                }
            } else if (nodeTypes[i]==MAXPOOL) {

                // Approx from below
                for (auto it : nodeConnectionIn[i]) {
                    std::cout << "  " << nodeNames[i] << " - " << nodeNames[it.first] << " >= 0\n";
                }

                // Approx from above
                std::cout << "  " << nodeNames[i];
                for (auto it : nodeConnectionIn[i]) {
                    std::cout << " - " << nodeNames[it.first];
                }
                std::cout << " <= 0\n";

            } else {
                throw "Unimplemented 1927!";
            }
        }
    }



    // Constraints
    for (auto constraint : constraints) {
        assert(std::get<0>(constraint)==LEQ);
        std::cout << "  ";
        bool first = true;
        for (auto it : std::get<2>(constraint)) {
            if (!first) std::cout << " + ";
            std::cout << it.second << " " << nodeNames[it.first];
            first = false;
        }
        std::cout << " >= " << std::get<1>(constraint) << std::endl;
    }

    std::cout << "Bounds\n";
    for (unsigned int i=0;i<nodeNames.size();i++) {
        std::cout << "  -infty <= " << nodeNames[i] << " <= infty\n";
    }

    // Phase variables
    std::cout << "Binary\n";
    for (unsigned int i=0;i<nodeNames.size();i++) {
        if (nodeTypes[i]==MAXPOOL) {
            for (unsigned int j=0;j<nodeConnectionIn[i].size();j++)
                std::cout << "  " << nodeNames[i] << "phase" << j << std::endl;
        }
        if (nodeTypes[i]==RELU) {
            std::cout << "  " << nodeNames[i] << "phase" << std::endl;
        }
    }


    std::cout << "End\n";
}

/**
 * @brief Print an SMT instance
 * @param includingApproximationConstraints
 */
void VerificationProblem::printSMTLibInstance(bool includingApproximationConstraints) {
    std::cout << "(set-logic QF_LRA)\n";
    for (unsigned int i=0;i<nodeNames.size();i++) {
        std::cout << "(declare-fun " << nodeNames[i] << " () Real)\n";
    }
    std::cout << "(define-fun maxRELU ((a Real)) Real (ite (<= 0 a) a 0))\n";
    // Define max functions
    size_t maxPoolInputLength = 2;
    for (unsigned int i=0;i<nodeTypes.size();i++) {
        if (nodeTypes[i]==MAXPOOL) maxPoolInputLength = std::max(maxPoolInputLength,nodeConnectionIn[i].size());
    }
    std::cout << "(define-fun maxi2 ((a Real) (b Real)) Real (ite (<= a b) b a))\n";
    for (unsigned int i=3;i<=maxPoolInputLength;i++) {
        std::cout << "(define-fun maxi" << i << " (";
        for (unsigned int j=0;j<i;j++) std::cout << "(x" << j << " Real)";
        std::cout << ") Real (ite (<= x0 x1) (maxi" << i-1 << " x1";
        for (unsigned int j=2;j<i;j++) std::cout << " x" << j;
        std::cout << ") (maxi" << i-1 << " x0";
        for (unsigned int j=2;j<i;j++) std::cout << " x" << j;
        std::cout << ")))\n";
    }



    // Define propagation

    for (unsigned int i=0;i<nodeNames.size();i++) {
        if (nodeTypes[i]==INPUT) {
            // Nothing to be done
        } else if (nodeTypes[i]==LINEAR) {
            std::cout << "(assert (= " << nodeNames[i] << " (+ ";
            std::cout << doubleToSMTLib2Expression(nodeBias[i]) << " ";
            for (auto it : nodeConnectionIn[i]) {
                std::cout << "(* " << doubleToSMTLib2Expression(it.second) << " " << nodeNames[it.first] << ")";
            }
            std::cout << ")))\n";
        } else if (nodeTypes[i]==RELU) {
            std::cout << "(assert (= " << nodeNames[i] << "( maxRELU (+ ";
            std::cout << doubleToSMTLib2Expression(nodeBias[i]) << " ";
            for (auto it : nodeConnectionIn[i]) {
                std::cout << "(* " << doubleToSMTLib2Expression(it.second) << " " << nodeNames[it.first] << ")";
            }
            std::cout << "))))\n";
        } else if (nodeTypes[i]==MAXPOOL) {
            if (nodeConnectionIn[i].size()==1) {
                std::cout << "(assert (= " << nodeNames[i] << " " << nodeNames[nodeConnectionIn[i][0].first] << "))\n";
            } else {
                std::cout << "(assert (= " << nodeNames[i] << "(maxi" << nodeConnectionIn[i].size();
                for (auto it : nodeConnectionIn[i]) {
                    std::cout << " " << nodeNames[it.first];
                }
                std::cout << ")))\n";
            }
        } else {
            throw "Unsupported.";
        }
    }

    if (includingApproximationConstraints) {

        // Upper and lower bound hints
        for (unsigned int i=0;i<nodeTypes.size();i++) {
            double lower = initialNeuronLimitBounds[i].first;
            double upper = initialNeuronLimitBounds[i].second;
            if (lower<0) lower = lower*(1.0+EPSILON_FOR_CALLING_AN_APPROXIMATED_NODE_VALUE_TO_BE_PRECISE_RELATIVE_TO_BOUNDS)-EPSILON_FOR_CALLING_AN_APPROXIMATED_NODE_VALUE_TO_BE_PRECISE_RELATIVE_TO_BOUNDS;
            if (lower>0) lower = lower*(1.0-EPSILON_FOR_CALLING_AN_APPROXIMATED_NODE_VALUE_TO_BE_PRECISE_RELATIVE_TO_BOUNDS)-EPSILON_FOR_CALLING_AN_APPROXIMATED_NODE_VALUE_TO_BE_PRECISE_RELATIVE_TO_BOUNDS;
            if (upper<0) upper = upper*(1.0-EPSILON_FOR_CALLING_AN_APPROXIMATED_NODE_VALUE_TO_BE_PRECISE_RELATIVE_TO_BOUNDS)+EPSILON_FOR_CALLING_AN_APPROXIMATED_NODE_VALUE_TO_BE_PRECISE_RELATIVE_TO_BOUNDS;
            if (upper>0) upper = upper*(1.0+EPSILON_FOR_CALLING_AN_APPROXIMATED_NODE_VALUE_TO_BE_PRECISE_RELATIVE_TO_BOUNDS)+EPSILON_FOR_CALLING_AN_APPROXIMATED_NODE_VALUE_TO_BE_PRECISE_RELATIVE_TO_BOUNDS;
            std::cout << "(assert (>= " << nodeNames[i] << " " << doubleToSMTLib2Expression(lower) << "))\n";
            std::cout << "(assert (<= " << nodeNames[i] << " " << doubleToSMTLib2Expression(upper) << "))\n";
        }

        for (unsigned int i=0;i<nodeTypes.size();i++) {

            if ((nodeTypes[i]==INPUT) || (nodeTypes[i]==LINEAR)) {
                // Nothing to do, as there are no additional constraints
            } else if (nodeTypes[i]==RELU) {
                // Add "AT LEAST" constraint
                std::cout << "(assert (>= " << nodeNames[i] << "(+ ";
                std::cout << doubleToSMTLib2Expression(nodeBias[i]) << " ";
                for (auto it : nodeConnectionIn[i]) {
                    std::cout << "(* " << doubleToSMTLib2Expression(it.second) << " " << nodeNames[it.first] << ")";
                }
                std::cout << ")))\n";

                // Add LEQ constraint

                // Compute theoretical minimum
                double min = nodeBias[i];
                double max = nodeBias[i];
                for (auto pair : nodeConnectionIn[i]) {
                    if (pair.second<0.0) {
                        min += pair.second*initialNeuronLimitBounds[pair.first].second;
                        max += pair.second*initialNeuronLimitBounds[pair.first].first;
                    } else {
                        min += pair.second*initialNeuronLimitBounds[pair.first].first;
                        max += pair.second*initialNeuronLimitBounds[pair.first].second;
                    }
                }
                max = std::min(initialNeuronLimitBounds[i].second,max);

                if (max<=0) {
                    std::cout << "(assert (= " << nodeNames[i] << " 0))\n";
                } else {
                    if (min<0) {
                        // Compute line between (min,0) and (max,max)
                        double factor = max/(max-min);
                        double offset = max-max*factor;

                        std::cout << "(assert (>= " << doubleToSMTLib2Expression(nodeBias[i]*factor+offset) << " " << " (+ " << nodeNames[i];

                        for (auto pair : nodeConnectionIn[i]) {
                            std::cout << " (* " << doubleToSMTLib2Expression(-1.0*factor*pair.second) << " " << nodeNames[pair.first] << ")";
                        }
                        std::cout << ")))\n";

                    } else {
                        // Make it linear
                        std::cout << "(assert (= " << nodeNames[i] << " (+ ";
                        std::cout << doubleToSMTLib2Expression(nodeBias[i]) << " ";
                        for (auto it : nodeConnectionIn[i]) {
                            std::cout << "(* " << doubleToSMTLib2Expression(it.second) << " " << nodeNames[it.first] << ")";
                        }
                        std::cout << ")))\n";
                    }
                }
            } else if (nodeTypes[i]==MAXPOOL) {

                // Approx from below
                for (auto it : nodeConnectionIn[i]) {
                    std::cout << "(assert (>= " << nodeNames[i] << " " << nodeNames[it.first] << "))\n";
                }

                // Approx from above
                if (nodeConnectionIn[i].size()>1) {
                    std::cout << "(assert (<= " << nodeNames[i] << " (+";
                    for (auto it : nodeConnectionIn[i]) {
                        std::cout << " " << nodeNames[it.first];
                    }
                    std::cout << ")))\n";
                }


            } else {
                throw "Unimplemented 1497!";
            }
        }
    }


    // Constraints
    for (auto constraint : constraints) {
        assert(std::get<0>(constraint)==LEQ);
        std::cout << "(assert (<= " << doubleToSMTLib2Expression(std::get<1>(constraint)) << "(+ ";
        for (auto it : std::get<2>(constraint)) {
            std::cout << "(* " << doubleToSMTLib2Expression(it.second) << " " << nodeNames[it.first] << ")";
        }
        std::cout << ")))\n";

    }
    std::cout << "(check-sat)\n";

}
