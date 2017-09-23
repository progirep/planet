#ifndef __LP_SOLVERS_HPP____
#define __LP_SOLVERS_HPP____

#include <sstream>
#include <iostream>

//=====================================================================================================
// LP Solving Library for the ReLu Verifier
//
// Written for the case that all variables in the LP problem have natural upper and lower bounds
//
// Usage notes:
// - When adding constraints and an objective function, the first element of the array
//   supplied that is looked at is element no. 1, and not 0. However, when giving variable
//   numbers to the library, the first variable has number 0.
//=====================================================================================================

#ifdef USE_GLPK
extern "C" {
    #include "glpk.h"
}
#else
extern "C" {
    #include "lp_lib.h"
}
#endif

class NumericalFailureException {
public:
    NumericalFailureException() {}
};

/**
 * @brief LPSolve Instance encapsulator
 */
#ifndef USE_GLPK
#define LPREALTYPE REAL
#define LPINFEASIBLE INFEASIBLE
#define LPUNBOUNDED UNBOUNDED
#define LPOPTIMUM OPTIMAL

class LPProblem {
private:
    lprec *lp;
public:
    LPProblem(unsigned int nofColumns) {
        lp = make_lp(0,nofColumns);
        if (lp==NULL) throw "Out of memory for LP allocation.";
        for (unsigned int i=0;i<nofColumns;i++) {
            if (set_bounds(lp, i+1, -1*get_infinite(lp), get_infinite(lp))==0) throw "Error setting variable bounds";
        }
        set_verbose(lp, 2);
    }
    ~LPProblem() { if (lp!=NULL) delete_lp(lp); }
    void addRowModeOn() {
        set_add_rowmode(lp, TRUE);
    }
    void addRowModeOff() {
        set_add_rowmode(lp, FALSE);
    }
    /**
     * @brief Adds a constraint
     * @param row pointer to an array of all factors. The array starts at index 1 and the value at position 0 is ignored.
     * @param constraint_type constraint_type Should be ROWTYPE_EQ, ROWTYPE_GEQ, or so on.
     * @param rh the right-hand-side (a constant)
     */
    void addEQConstraint(LPREALTYPE *row, LPREALTYPE rh) {
        if (add_constraint(lp, row, EQ, rh)==0)
            throw "Error adding constraint to LP problem.";
    }
    void addGEQConstraint(LPREALTYPE *row, LPREALTYPE rh) {
        if (add_constraint(lp, row, GE, rh)==0)
            throw "Error adding constraint to LP problem.";
    }
    void addLEQConstraint(LPREALTYPE *row, LPREALTYPE rh) {
        if (add_constraint(lp, row, LE, rh)==0)
            throw "Error adding constraint to LP problem.";
    }

    void setMinim() {
        set_minim(lp);
    }
    void setMaxim() {
        set_maxim(lp);
    }
    void setObjFn(LPREALTYPE *row) {
        if (set_obj_fn(lp,row)==FALSE) throw "Error setting the LP objective function.";
    }
    int solve() {
        int lpRes = ::solve(lp);
        if (lpRes==PRESOLVED) {
            set_presolve(lp, PRESOLVE_NONE, 0);
            lpRes = ::solve(lp);
        }
        if (lpRes==NUMFAILURE) {
            throw NumericalFailureException();
        }
        if ((lpRes!=LPUNBOUNDED) && (lpRes!=LPINFEASIBLE) && (lpRes!=OPTIMAL) ) {
            std::ostringstream err;
            err << "LPSolve returned error " << lpRes;
            throw err.str();
        }
        return lpRes;
    }
    double getObjective() {
        double data = get_objective(lp);

        if (data <= -1*getInfinite()) throw "Error: Objective has an open lower bound";
        if (data >= getInfinite()) throw "Error: Objective has an open upper bound";
        return data;
    }
    double getInfinite() {
        return get_infinite(lp);
    }
    void writeLP(const char *filename) {
        if (write_lp(lp,const_cast<char*>(filename))==0) throw "Error writing LP instance to file.";
    }
    void setBounds(unsigned int varNumberWithoutPlusOne,double lower, double upper) {
        if (set_bounds(lp, varNumberWithoutPlusOne+1, lower, upper)==0) throw "Error setting variable bounds";
    }
    double getSolutionVarValue(unsigned int varNumWithoutPlusOne) {
        return get_var_primalresult(lp,varNumWithoutPlusOne+1+get_Nrows(lp));
    }
};
#endif



/**
 * @brief LPSolve Instance encapsulator
 */
#ifdef USE_GLPK
#define LPREALTYPE double
#define LPINFEASIBLE 12345
#define LPUNBOUNDED 12346
#define LPOPTIMUM 12347
class LPProblem {
private:
    glp_prob *lp;
    unsigned int nofVars;
    int *allVariables;
public:
    LPProblem(unsigned int nofColumns) : allVariables(NULL) {
        lp = glp_create_prob();
        glp_term_out(GLP_OFF);
        if (lp==NULL) throw "Out of memory for LP allocation.";
        glp_add_cols(lp, nofColumns);
        for (unsigned int i=1;i<=nofColumns;i++) {
            glp_set_col_bnds(lp, i, GLP_FR, 0.0, 0.0);
        }
        nofVars = nofColumns;
        allVariables = new int[nofVars+1];
        for (unsigned int i=0;i<=nofVars;i++) {
            allVariables[i] = i;
        }
    }

    // Default constructor -- Yields an unsable LPProblem.
    LPProblem() {
        lp = 0;
        allVariables = 0;
    }

    // Copying of this object is not allowed. However, a move is ok.
    LPProblem(const LPProblem &other) = delete;
    LPProblem(LPProblem &&other) {
        lp = other.lp;
        nofVars = other.nofVars;
        allVariables = other.allVariables;
        other.lp = 0;
        other.allVariables = 0;
    }
    LPProblem& operator=(const LPProblem &other) = delete;
    LPProblem& operator=(LPProblem &&other) {
        if (lp!=NULL) glp_delete_prob(lp);
        if (allVariables!=NULL) delete[] allVariables;
        lp = other.lp;
        nofVars = other.nofVars;
        allVariables = other.allVariables;
        other.lp = 0;
        other.allVariables = 0;
        return *this;
    }

    ~LPProblem() {
        if (lp!=NULL) glp_delete_prob(lp);
        if (allVariables!=NULL) delete[] allVariables;
    }
    void addRowModeOn() {
        // Nothing to do
    }
    void addRowModeOff() {
        // Nothing to do
    }
    /**
     * @brief Adds a constraint
     * @param row pointer to an array of all factors. The array starts at index 1 and the value at position 0 is ignored.
     * @param constraint_type constraint_type Should be ROWTYPE_EQ, ROWTYPE_GEQ, or so on.
     * @param rh the right-hand-side (a constant)
     */
    void addLEQConstraint(LPREALTYPE *coefficients, LPREALTYPE rh) {
        glp_add_rows(lp, 1);
        glp_set_row_bnds(lp, glp_get_num_rows(lp), GLP_UP, 0.0, rh);
        glp_set_mat_row(lp, glp_get_num_rows(lp), nofVars,allVariables,coefficients);
    }
    void addEQConstraint(LPREALTYPE *coefficients, LPREALTYPE rh) {
        glp_add_rows(lp, 1);
        glp_set_row_bnds(lp, glp_get_num_rows(lp), GLP_FX, rh, rh);
        glp_set_mat_row(lp, glp_get_num_rows(lp), nofVars,allVariables,coefficients);
    }
    void addGEQConstraint(LPREALTYPE *coefficients, LPREALTYPE rh) {
        glp_add_rows(lp, 1);
        glp_set_row_bnds(lp, glp_get_num_rows(lp), GLP_LO, rh, 0.0);
        glp_set_mat_row(lp, glp_get_num_rows(lp), nofVars,allVariables,coefficients);
    }
    void setMinim() {
        glp_set_obj_dir(lp, GLP_MIN);
    }
    void setMaxim() {
        glp_set_obj_dir(lp, GLP_MAX);
    }
    void setObjFn(LPREALTYPE *row) {
        for (unsigned int i=0;i<nofVars;i++) {
            glp_set_obj_coef(lp, i+1, row[i+1]);
        }
    }
    int getNofRows() {
        return glp_get_num_rows(lp);
    }
    void resetBasis() {
        glp_adv_basis(lp,0);
    }

    int solve() {
        //glp_cpx_basis(lp);
        int lpRes = glp_simplex(lp, NULL);

        // If it failed, try again immediately.
        if (lpRes==GLP_EFAIL) {
            std::cerr << "GLP Failed. Trying again once.\n";
            lpRes = glp_simplex(lp, NULL);
        }

        // If the solution failed with a EBADB error, then an added constraint made the initial basis incorrect. Then we have to recompute it and try again.
        if (lpRes==GLP_EBADB) {
            //glp_adv_basis(lp, 0);
            glp_adv_basis(lp,0);
            lpRes = glp_simplex(lp, NULL);
        }

        switch (lpRes) {
        case GLP_EFAIL:
        case GLP_EBADB:
        case GLP_ESING:
        case GLP_ECOND:
        case GLP_EOBJLL:
        case GLP_EOBJUL:
        case GLP_EITLIM:
        case GLP_ETMLIM:
        case GLP_ENOPFS:
        case GLP_ENODFS:
        {
                std::ostringstream err;
                err << "GLP returned error " << lpRes;
                if (lpRes==GLP_EFAIL) err << " (GLP_EFAIL)";
                else if (lpRes==GLP_EBADB) err << " (GLP_EFAIL)";
                else if (lpRes==GLP_ESING) err << " (GLP_ESING)";
                else if (lpRes==GLP_ECOND) err << " (GLP_ECOND)";
                else if (lpRes==GLP_EOBJLL) err << " (GLP_EOBJLL)";
                else if (lpRes==GLP_EOBJUL) err << " (GLP_EOBJUL)";
                else if (lpRes==GLP_EITLIM) err << " (GLP_EITLIM)";
                else if (lpRes==GLP_ENOPFS) err << " (GLP_ENOPFS)";
                writeLP("/tmp/glp_error.lp");
                throw err.str();
        }
        case 0:
        {
            int status = glp_get_status(lp);
            if (status==GLP_OPT) return LPOPTIMUM;
            if (status==GLP_FEAS) return LPOPTIMUM; // TODO: Check if optimization function was really unset. Otherwise this should be an error
            if (status==GLP_INFEAS) return LPINFEASIBLE;
            if (status==GLP_NOFEAS) return LPINFEASIBLE;
            if (status==GLP_UNBND) return LPUNBOUNDED;
            if (status==GLP_UNDEF) throw "LP problem has an undefined solution";
            return LPOPTIMUM;
        }
        case GLP_EBOUND:
            writeLP("/tmp/glp_bounds.lp");
            throw "Incorrectly set bounds.";

        }
        throw "Error: Unknown GLPK Return status.";
    }
    double getObjective() {
        return glp_get_obj_val(lp);
    }
    void writeLP(const char *filename) {
        if (glp_write_lp(lp, NULL,filename)!=0) throw "Error writing LP instance to file.";
    }
    void setBounds(unsigned int varNumberWithoutPlusOne,double lower, double upper) {
        if (lower==upper) {
            glp_set_col_bnds(lp, varNumberWithoutPlusOne+1, GLP_FX, lower, upper);
        } else {
            glp_set_col_bnds(lp, varNumberWithoutPlusOne+1, GLP_DB, lower, upper);
        }
    }
    void setBoundVariableOnlyFromBelow(unsigned int varNumberWithoutPlusOne, double lower) {
        glp_set_col_bnds(lp,varNumberWithoutPlusOne+1,GLP_LO,lower,0.0);
    }
    double getSolutionVarValue(unsigned int varNumWithoutPlusOne) {
        return glp_get_col_prim(lp, varNumWithoutPlusOne+1);
    }
    void deleteLastRow() {
        int rowNum[2];
        rowNum[1] = glp_get_num_rows(lp);
        glp_del_rows(lp,1,rowNum);


    }

    /**
     * @brief Delete rows
     * @param Row number. The first row has number "1".
     */
    void deleteRow(unsigned int _rowNum) {
        int rowNum[2];
        rowNum[1] = _rowNum;
        glp_del_rows(lp,1,rowNum);
    }

    void deleteCloseToLastRow(unsigned int lastRowsToKeep) {
        int rowNum[2];
        rowNum[1] = glp_get_num_rows(lp)-lastRowsToKeep;
        glp_del_rows(lp,1,rowNum);


    }
};
#endif












#endif
