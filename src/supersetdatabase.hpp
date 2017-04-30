#ifndef __SUPERSETDATABASE_HPP__
#define __SUPERSETDATABASE_HPP__

#include <vector>
#include <cstdlib>

/**
 * @brief A special database class that stores subsets of {0, ..., n} for some a-priori defined number n.
 * The main query for this database is to slowly add numbers from {0, ..., n} to a trace, and to always return
 * whether the database (still) contains a superset. Traces are like stacks -- only the respective last element
 * can be removed.
 *
 * Adding elements can only be done if for the current trace there is no match.
 */
class SupersetDatabase {
private:
    std::vector<int*> setsInDatabase; // The first element of every int* is the number of elements
    int nofObjects;
    unsigned int currentSetPointer;
    int *trace;
    unsigned int *traceDatabasePointers;
    unsigned int nofElementsInTrace;
public:
    SupersetDatabase(int _nofObjects) : nofObjects(_nofObjects), currentSetPointer(0),nofElementsInTrace(0) {
        trace = new int[_nofObjects];
        traceDatabasePointers = new unsigned int[_nofObjects];
    }
    ~SupersetDatabase() {
        for (auto it : setsInDatabase) delete[] it;
        delete[] traceDatabasePointers;
        delete[] trace;
    }

    // Copy constructor - Not actually quick, as it rebuilds the trace here.
    SupersetDatabase(const SupersetDatabase &other) {
        nofObjects = other.nofObjects;
        trace = new int[nofObjects];
        traceDatabasePointers = new unsigned int[nofObjects];
        setsInDatabase.resize(other.setsInDatabase.size());
        for (unsigned int i=0;i<other.setsInDatabase.size();i++) {
            int *otherSet = other.setsInDatabase[i];
            int *newOne = new int[otherSet[0]+1];
            for (int j=0;j<=otherSet[0];j++) newOne[j] = otherSet[j];
            setsInDatabase[i] = newOne;
        }
        currentSetPointer = 0;
        nofElementsInTrace = 0;
        for (unsigned int i=0;i<other.nofElementsInTrace;i++) addTraceElement(other.trace[i]);
    }

    // Assignment operator. Not actually quick, as it copies the whole database
    SupersetDatabase& operator=(const SupersetDatabase &other) {
        for (auto it : setsInDatabase) delete[] it;
        delete[] traceDatabasePointers;
        delete[] trace;
        nofObjects = other.nofObjects;
        trace = new int[nofObjects];
        traceDatabasePointers = new unsigned int[nofObjects];
        setsInDatabase.resize(other.setsInDatabase.size());
        for (unsigned int i=0;i<other.setsInDatabase.size();i++) {
            int *otherSet = other.setsInDatabase[i];
            int *newOne = new int[otherSet[0]+1];
            for (int j=0;j<=otherSet[0];j++) newOne[j] = otherSet[j];
            setsInDatabase[i] = newOne;
        }
        currentSetPointer = 0;
        nofElementsInTrace = 0;
        for (unsigned int i=0;i<other.nofElementsInTrace;i++) addTraceElement(other.trace[i]);
        return *this;
    }



    void addSet(std::vector<int> const &elements);
    bool addTraceElement(int element);
    int removeElementFromTrace();
    void resetTrace();

};



#endif
