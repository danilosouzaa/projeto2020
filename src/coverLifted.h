#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <omp.h>

#include "constraintsManipulation.h"


int *LCIAdam(int *coverSolution, cutCover *constraintsCover, TNumberConstraints constraint);

int *LCIBallas(int *coverSolution, cutCover *constraintsCover, TNumberConstraints constraint);

constraintsReal *runCCwithGrasp(constraintsReal *constraintsFull, int precision, char **nameConstraints, char **nameVariables, int szPoolCutsMax, int nIterationGrasp, float alpha, int minimal, int typeLifted);

int *createCoverGraspIndividual(cutSmall *constraintsSmall, int precision, TNumberConstraints constraint, int numberIteration, int szPoolCutsMax, float alpha, int minimal, int typeLifted, int *nPoolTemp);

cutSmall *reduceCutFullForCutSmall(constraintsReal *constraints, int *typeIntOrFloat, int precision);

void createInitialCoverGRASP(int *solution, int sz, cutSmall *constraintsSmall, int precision, int constraint, float alpha, int typeLifted);

void copyAndVerifyPoolSolution(int *solution, int sz, int *poolSolution, int *numberSolutionAtual);

int verifySolutionCoverMinimal(int *solution, cutSmall *constraintsSmall, TNumberConstraints constraint);

constraintsReal *createCutsCoverGrasp(cutCover *cutsCover, constraintsReal *constraintsOriginal, cutSmall *constraintsSmall, int *idc_Cover, int constraint, int nCuts, int precision);

int *localSearch(int *solution, int sz, cutSmall *constraintsSmall, int precision, TNumberConstraints constraint, int minimal);

int verifySolutionCover(int *solution, cutSmall *constraintsSmall, int precision, TNumberConstraints constraint);

void shuffleVectorInt(int *vec, int sz);

int verifyViolationGreedy(int *solutionFinal, cutSmall *constraitsUsed, int constraint, int precision);

constraintsReal *runCCGreedy(constraintsReal *constraintsFull, int precision, char **nameConstraints, char **nameVariables, int typeLifted);

int *greedyMethodForCC(cutSmall *constraintsUsed, TNumberConstraints constraint, int precision, int typeLift);