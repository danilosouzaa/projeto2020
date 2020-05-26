#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "lp.h"
#include "defs.h"
#include <omp.h>
#include <sys/time.h>

char **createNameVariablesInitial(LinearProgramPtr lp); //read name variables and create struct with name

TNumberConstraints countConstraintsValided(LinearProgramPtr lp); // count number of constraints validated for create strutc

char **createStructNameConstraintsInitial(LinearProgramPtr lp);// create struct name Constraints - no fill

double *initialSolutionForValidation(LinearProgramPtr lp);

double *readSolFile(const char *name, int nVariables);

void freeStructName(char **name, TNumberVariables sz);

int verifyOfFloatIsInteger(TCoefficientsFull coef);

char **renamedNameConstraints(char **nameConstraints, int typeContraints, TNumberConstraints szNewConstraints, TNumberConstraints szCuts, TNumberConstraints lastCut);

void insertConstraintsLPDebug(LinearProgramPtr lp, constraintsReal *constraintsOriginal, int nConstrainsInitial, char **nameConstraints, int *verifyBug);

void insertConstraintsLP(LinearProgramPtr lp, constraintsReal *constraintsOriginal, int nConstrainsInitial, char **nameConstraints);

int *vectorNonRepeteadNonDominated(constraintsReal *constraintsOriginal, int nConstraintsInitial);

TCoefficients cutMaxDivisorCommonVector(TCoefficients coefs[], TNumberVariables nElem);

TCoefficients cutMaxDivisorCommonRec(TCoefficients m, TCoefficients n);

int verifyRepeatCuts(constraintsReal *constraintsOriginal, int cutOriginal, int cutCreate);

double fRand(double fMin, double fMax);

