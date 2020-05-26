#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <assert.h>
#include <math.h>

#include "preProcessing.h"

constraintsReal *AllocStrConstraintsReal(TCont cont, TNumberConstraints nConstraints, TNumberVariables nVariables);

cutSmall *AllocStrCutSmall(TCont cont, TNumberConstraints nConstraints, TNumberVariables nVariables);

void freeStrConstraintsReal(constraintsReal *cut);

constraintsReal *fillStructPerLP(LinearProgram *lp, char **nameConstraints, char **nameVariables, int *typeVariables, double *lbVariables, double *ubVariables);

void showStructFull(constraintsReal *constraintsFull, char **nameConstraints, char **nameVariables);

int verifyCutsValidatedPerSolutionInteger(constraintsReal *constraintsOriginal, int cut, double *sol, char **nameVariables);

int *returnBinaryConstraints(constraintsReal *constraintsOriginal, int *typeVariables);

TNumberConstraints countConstraintsBinaryUsed(int *binaryConstraints, int sz);

constraintsReal *convertBinaryConstraints(constraintsReal *constraintsOriginal, int *BinaryConstraints, int *typeVariables, double *lbVariables, double *ubVariables);

constraintsReal *removeNegativeCoefficientsAndSort(constraintsReal *constraintsOriginal, int *convertVector);

void SortByCoefficients(constraintsReal *h_cut);

void quicksortCof(TCoefficientsFull *values, int *idc, int began, int end);

void quicksortDouble(double *values, int began, int end);

int *returnVectorTypeContraintsIntOrFloat(constraintsReal *constraints);

int verifyOfFloatIsInteger(TCoefficientsFull coef);

cutCover *CopyCutToCover(cutSmall *h_cut);

cutCover *AllocStrCover(TCont cont, TNumberConstraints nConstraints);

int verifyRepeated(constraintsReal *originalConstraints, int posCover);

int verifyRepeatedIncidency(int **matrizIncidencia, constraintsReal *originalConstraints, int posCover);

double valueViolation(cutCover *cCover, cutSmall *constraintsSmall, TNumberConstraints idCover, TNumberConstraints ogConstraint, int precision);

void quicksortTParameters(TParametersGroup *values, TNumberConstraints *idc, TNumberConstraints began, TNumberConstraints end);

constraintsReal *returnVariablesOriginals(constraintsReal *constraintsOriginal, int *convertVector, int precision, int nVariablesInitial);

constraintsReal *convertBinaryOfOriginalConstraints(constraintsReal *constraintsOriginal, constraintsReal *constraintsBinary, int nInitialBinary);