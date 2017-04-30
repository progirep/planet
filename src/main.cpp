#include "verifierContext.hpp"
#include <iostream>

void printSMTLibInstance(bool includingApproximationConstraints);


int main(int argc, char* argv[])
{
    try {
        VerificationProblem problem;
        int printSMTLibInstance = 0;
        int printILPInstance = 0;
        std::string inFileName = "";
        for (unsigned int i=1;i<static_cast<unsigned int>(argc);i++) {
            if (std::string(argv[i]).substr(0,1)=="-") {
                // Parameters
                if (std::string(argv[i])=="--smtlib") {
                    printSMTLibInstance = 1;
                } else if (std::string(argv[i])=="--smtlibApprox") {
                    printSMTLibInstance = 2;
                } else if (std::string(argv[i])=="--ilp") {
                    printILPInstance = 1;
                } else if (std::string(argv[i])=="--ilpApprox") {
                    printILPInstance = 2;
                } else {
                    std::cerr << "Error: Did not understand parameter '" << argv[i] << "'\n";
                    return 1;
                }
            } else {
                if (inFileName=="") {
                    inFileName = argv[i];
                } else {
                    throw "Error: More than one file name given.";
                }
            }
        }
        problem.load(inFileName);

        if (printSMTLibInstance==1) { // No Approx
            problem.printSMTLibInstance(false);
        } else if (printILPInstance==1) {
            problem.printILPInstance(false);
        } else {
            if (!(problem.computeInitialNeuronLimitBounds())) {
                if (printSMTLibInstance>0)
                    std::cout << "(set-logic QF_LRA) (declare-fun a () Bool) (assert a) (assert (not a)) (check-sat)\n";
                else if (printILPInstance>0)
                    std::cout << "Subject To\n  x <= 0\n  x >= 1\nEnd\n";
                else
                    std::cout << "UNSAT\n";

            } else {
                if (printILPInstance>0) {
                    problem.printILPInstance(true);
                } else if (printSMTLibInstance>0) {
                    problem.printSMTLibInstance(true);
                } else {
                    problem.verify();
                }
            }
        }

    } catch (const char *c) {
        std::cerr << "Error: " << c << std::endl;
        return 1;
    } catch (std::string c) {
        std::cerr << "Error: " << c << std::endl;
        return 1;
    }

}
