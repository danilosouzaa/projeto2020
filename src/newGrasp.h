#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <omp.h>


//#include "constraintsManipulation.h"
#include "coverLifted.h"

constraintsReal *graspLciAdam(constraintsReal *constraintsFull, int precision, char **nameConstraints, char **nameVariables, int szPoolCutsMax, int nIterationGrasp, float alpha, int minimal);

int *coverLCIAdamPerConstraints(cutSmall *knapsackConstraints, int precision, int constraint, int nIterationGrasp, int szPoolCutsMax, float alpha, int minimal, int *nPoolUsed);

void initialLCIAdamGRASP(int *solution, int sz,constraintMin *cMin, cutSmall *constraintsSmall, int precision, int constraint, float alpha);

constraintMin *remove_X_0_and_1_one_Constraints(cutSmall *knapsackConstraints, int constraint, int precision);

void copyPosGraspCMin(int *solution, int *solutionTemp, constraintMin *cMin, cutSmall *knapsackConstraints, int constraint, int precision);
