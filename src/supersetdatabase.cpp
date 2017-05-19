#include "supersetdatabase.hpp"
#include <algorithm>
#include <array>
#include <iostream>

void SupersetDatabase::addSet(std::vector<int> const &elements) {
    int *cpy = new int[elements.size()+1];
    for (unsigned int i=0;i<elements.size();i++) { cpy[i+1] = elements.at(i); }
    std::sort(cpy+1,cpy + elements.size()+1);
    cpy[0] = elements.size();
    setsInDatabase.push_back(cpy);

    // If we are beyond the last element, we have to check the remaining one with the trace.
    if (setsInDatabase.size()!=currentSetPointer+1) return;

    int *currentPtr = setsInDatabase[currentSetPointer];
    int currentSize = currentPtr[0];

    // Search for a match
    int *backPtrSetInDB = currentPtr+currentSize;
    int *backPtrInTrace = trace+nofElementsInTrace-1;
    while ((backPtrInTrace>trace) && (backPtrSetInDB>currentPtr)) {
        int diff = (*backPtrSetInDB);
        diff -= *backPtrInTrace;
        if (diff==0) {
            backPtrSetInDB--;
            backPtrInTrace--;
        } else if (diff<0) {
            backPtrSetInDB = currentPtr;
        } else {
            backPtrSetInDB--;
        }
    }

    if (backPtrInTrace<trace) return;
    currentSetPointer++;

}

bool SupersetDatabase::addTraceElement(int element) {
    trace[nofElementsInTrace] = element;
    traceDatabasePointers[nofElementsInTrace] = currentSetPointer;
    nofElementsInTrace++;

    // Case 1: We are at the end already...
    if (currentSetPointer>=setsInDatabase.size()) return false;

    // Case 2: The current pointer is valid! Then check if the current set also contains this element.
    int *currentPtr = setsInDatabase[currentSetPointer];
    int currentSize = currentPtr[0];
    if (std::binary_search(currentPtr+1,currentPtr+currentSize+1,element)) return true;

    // Case 3: We have to search for another set containing all other elements
    std::vector<int> sortedElementsInTrace(nofElementsInTrace);
    for (unsigned int i=0;i<nofElementsInTrace;i++) sortedElementsInTrace[i] = trace[i];
    std::sort(sortedElementsInTrace.begin(),sortedElementsInTrace.end());
    currentSetPointer++;
    while (currentSetPointer < setsInDatabase.size()) {

        int *currentPtr = setsInDatabase[currentSetPointer];
        int currentSize = currentPtr[0];

        // Search for a match
        int *backPtrSetInDB = currentPtr+currentSize;
        int *backPtrInTrace = sortedElementsInTrace.data()+nofElementsInTrace-1;
        while ((backPtrInTrace>=sortedElementsInTrace.data()) && (backPtrSetInDB>currentPtr)) {
            int diff = (*backPtrSetInDB) - *backPtrInTrace;
            if (diff==0) {
                backPtrSetInDB--;
                backPtrInTrace--;
            } else if (diff<0) {
                backPtrSetInDB = currentPtr;
            } else {
                backPtrSetInDB--;
            }
        }

        if (backPtrInTrace<sortedElementsInTrace.data()) return true;
        currentSetPointer++;
    }
    return false;
}

int SupersetDatabase::removeElementFromTrace() {
    nofElementsInTrace--;
    currentSetPointer = traceDatabasePointers[nofElementsInTrace];
    return trace[nofElementsInTrace];
}

void SupersetDatabase::resetTrace() {
    nofElementsInTrace = 0;
    currentSetPointer = 0;
}

