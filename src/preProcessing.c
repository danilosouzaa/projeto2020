#include "preProcessing.h"

char **createNameVariablesInitial(LinearProgramPtr lp){
    int i;
    int nVar = lp_cols(lp);
    char **nameVariables = (char **)malloc(nVar * sizeof(char *));
    for (i = 0; i < nVar; i++)
    {
        nameVariables[i] = (char *)malloc(255 * sizeof(char));
        lp_col_name(lp, i, nameVariables[i]);
    }
    return nameVariables;
}

TNumberConstraints countConstraintsValided(LinearProgramPtr lp)
{
    TNumberConstraints nConstraints, nConstraintsValided;
    TNumberVariables nVariables;
    nVariables = lp_cols(lp);
    nConstraints = lp_rows(lp);
    nConstraintsValided = 0;
    int szConst, i, j;
    TElements *idx = (TElements*)malloc(sizeof(TElements) * nVariables);
    TCoefficientsFull *coef = (TCoefficientsFull *)malloc(sizeof(TCoefficientsFull) * nVariables);
    for (i = 0; i < nConstraints; i++)
    {
        if ((lp_sense(lp, i) == 'L') || (lp_sense(lp, i) == 'G'))
        {
            nConstraintsValided++;
            for (j = 0; j < nVariables; j++)
            {
                coef[j] = 0.0;
                idx[j] = 0;
            }
            lp_row(lp, i, idx, coef);
            szConst = 0;
            for (j = 0; j < nVariables; j++)
            {
                if (coef[j] != 0)
                {
                    szConst++;
                }
            }
            if (szConst < 2)
            {
                nConstraintsValided--;
            }
        }
        else
        {
            nConstraintsValided += 2;
            for (j = 0; j < nVariables; j++)
            {
                coef[j] = 0.0;
                idx[j] = 0;
            }
            lp_row(lp, i, idx, coef);
            szConst = 0;
            for (j = 0; j < nVariables; j++)
            {
                if (coef[j] != 0)
                {
                    szConst++;
                }
            }
            if (szConst < 2)
            {
                nConstraintsValided -= 2;
            }
        }
    }
    free(idx);
    free(coef);
    return nConstraintsValided;
}

char **createStructNameConstraintsInitial(LinearProgramPtr lp){
    int i, nConstraintsValid = countConstraintsValided(lp);
    char **nameConstraints = (char **)malloc(nConstraintsValid * sizeof(char *));
    for (i = 0; i < nConstraintsValid; i++)
    {
        nameConstraints[i] = (char *)malloc(255 * sizeof(char));
    }
    return nameConstraints;

}
double *initialSolutionForValidation(LinearProgramPtr lp){
    int nVariables = lp_cols(lp);
    lp_set_max_solutions(lp, 1);
    lp_optimize(lp);
    lp_write_sol(lp, "../sol/solutionTest.sol");
    double *sol = readSolFile("../sol/solutionTest.sol", nVariables);
    return sol;
}

double *readSolFile(const char *name, int nVariables)
{
    double *sol = (double *)malloc(sizeof(double) * nVariables);
    int idx = -1, i;
    double value, aux;
    char n[255];
    for (i = 0; i < nVariables; i++)
    {
        sol[i] = 0;
    }
    FILE *f;
    f = fopen(name, "r");
    if (f != NULL)
    {
        int flag = 0;
        while (!feof(f))
        {
            if (flag == 0)
            {
                fscanf(f, "%[^\n]", n);
                flag = 1;
                continue;
            }
            fscanf(f, "%d %s %lf %lf\n", &idx, n, &value, &aux);
            if(idx !=-1){
                    sol[idx] = value;
            }
        }
    }
    return sol;
}

void freeStructName(char **name, TNumberVariables sz){
    int i;
    for (i = 0; i < sz; i++)
    {
        free(name[i]);
    }
    free(name);
}


char **renamedNameConstraints(char **nameConstraints, int typeContraints, TNumberConstraints szNewConstraints, TNumberConstraints szCuts, TNumberConstraints lastCut)
{
    int i;
    char **newNameConstraints = (char **)malloc(szNewConstraints * sizeof(char *));
    for (i = 0; i < szNewConstraints; i++)
    {
        newNameConstraints[i] = (char *)malloc(255 * sizeof(char));
        if (i < szNewConstraints - szCuts)
        {
            strcpy(newNameConstraints[i], nameConstraints[i]);
        }
    }
    for (i = szNewConstraints - szCuts; i < szNewConstraints; i++)
    {
        if (typeContraints == 1)
        {
            char name[255];
            sprintf(name, "CG1(%d)", lastCut);
            lastCut++;
            strcpy(newNameConstraints[i], name);
        }
        if (typeContraints == 2)
        {
            char name[255];
            sprintf(name, "CG2(%d)", lastCut);
            lastCut++;
            strcpy(newNameConstraints[i], name);
        }
        if (typeContraints == 3)
        {
            char name[255];
            sprintf(name, "CC(%d)", lastCut);
            lastCut++;
            strcpy(newNameConstraints[i], name);
        }
    }
    for (i = 0; i < szNewConstraints - szCuts; i++)
    {
        free(nameConstraints[i]);
    }
    free(nameConstraints);
    return newNameConstraints;
}



void insertConstraintsLPDebug(LinearProgramPtr lp, constraintsReal *constraintsOriginal, int nConstrainsInitial, char **nameConstraints, int *verifyBug)
{
    int i, j, w;
    int *idx;
    double *Coef;
    int *cutsNonInserted = vectorNonRepeteadNonDominated(constraintsOriginal, nConstrainsInitial);
    for (i = nConstrainsInitial; i < constraintsOriginal->numberConstraints; i++)
    {
        if ((cutsNonInserted[i] == 1)||(verifyBug[i]==0))
        {
            continue;
        }
        int sz = constraintsOriginal->ElementsConstraints[i + 1] - constraintsOriginal->ElementsConstraints[i];
        idx = (int *)malloc(sizeof(int) * sz);
        Coef = (double *)malloc(sizeof(double) * sz);
        for (j = constraintsOriginal->ElementsConstraints[i]; j < constraintsOriginal->ElementsConstraints[i + 1]; j++)
        {
            idx[w] = constraintsOriginal->Elements[j];
            Coef[w] = constraintsOriginal->Coefficients[j];
            w++;
        }
        double rhs = constraintsOriginal->rightSide[i];
        lp_add_row(lp, sz, idx, Coef, nameConstraints[i], 'L', rhs);
        //}
        w = 0;
        free(idx);
        free(Coef);
    }

    free(cutsNonInserted);
}


void insertConstraintsLP(LinearProgramPtr lp, constraintsReal *constraintsOriginal, int nConstrainsInitial, char **nameConstraints)
{
    int i, j, w;
    int *idx;
    double *Coef;
    int *cutsNonInserted = vectorNonRepeteadNonDominated(constraintsOriginal, nConstrainsInitial);
    for (i = nConstrainsInitial; i < constraintsOriginal->numberConstraints; i++)
    {
        if (cutsNonInserted[i] == 1)
        {
            continue;
        }
        int sz = constraintsOriginal->ElementsConstraints[i + 1] - constraintsOriginal->ElementsConstraints[i];
        idx = (int *)malloc(sizeof(int) * sz);
        Coef = (double *)malloc(sizeof(double) * sz);
        for (j = constraintsOriginal->ElementsConstraints[i]; j < constraintsOriginal->ElementsConstraints[i + 1]; j++)
        {
            idx[w] = constraintsOriginal->Elements[j];
            Coef[w] = constraintsOriginal->Coefficients[j];
            w++;
        }
        double rhs = constraintsOriginal->rightSide[i];
        lp_add_row(lp, sz, idx, Coef, nameConstraints[i], 'L', rhs);
        //}
        w = 0;
        free(idx);
        free(Coef);
    }

    free(cutsNonInserted);
}


int *vectorNonRepeteadNonDominated(constraintsReal *constraintsOriginal, int nConstraintsInitial)
{
    int i, j, k = 0, qntRepeat = 0;
    int *isRepeat = (int *)malloc(sizeof(int) * constraintsOriginal->numberConstraints);
    for (i = 0; i < nConstraintsInitial; i++)
    {
        isRepeat[i] = 0;
    }
    int *v_aux = (int *)malloc(sizeof(int) * (constraintsOriginal->numberVariables+1) );
    for (i = nConstraintsInitial; i < constraintsOriginal->numberConstraints; i++)
    {
        k = 0;
        for (j = constraintsOriginal->ElementsConstraints[i]; j < constraintsOriginal->ElementsConstraints[i + 1]; j++)
        {
            v_aux[k] = constraintsOriginal->Coefficients[j];
            k++;
        }
        v_aux[k] = constraintsOriginal->rightSide[i];
        int mdc = cutMaxDivisorCommonVector(v_aux, k);
        for (j = constraintsOriginal->ElementsConstraints[i]; j < constraintsOriginal->ElementsConstraints[i + 1]; j++)
        {
            constraintsOriginal->Coefficients[j] = constraintsOriginal->Coefficients[j] / mdc;
        }
        constraintsOriginal->rightSide[i] = constraintsOriginal->rightSide[i] / mdc;
        for (j = 0; j < i; j++)
        {
            isRepeat[i] = verifyRepeatCuts(constraintsOriginal, j, i);
            qntRepeat += isRepeat[i];
            if (isRepeat[i] == 1)
            {
                break;
            }
        }
    }
    //printf("ncRpt: %d\t", qntRepeat);
    free(v_aux);
    return isRepeat;
}

TCoefficients cutMaxDivisorCommonVector(TCoefficients coefs[], TNumberVariables nElem)
{
    TCoefficients n = coefs[nElem];
    TCoefficients mdc = 1;
    while (nElem > 0)
    {
        TCoefficients m = coefs[nElem - 1];
        mdc = cutMaxDivisorCommonRec(m, n);
        n = mdc;
        nElem--;
    }
    TCoefficients m = coefs[nElem];
    mdc = cutMaxDivisorCommonRec(m, n);
    n = mdc;
    return n;
}

/*calculates the maximum common divisor for two integers.*/
TCoefficients cutMaxDivisorCommonRec(TCoefficients m, TCoefficients n)
{

    TCoefficients t = 0;
    m = m < 0 ? -m : m; /* abs(u) */
    n = n < 0 ? -n : n;
    if (m < n)
    {
        t = m;
        m = n;
        n = t;
    }

    if (n == 0)
        return m;
    else
    {
        TCoefficients resto = m % n;
        return cutMaxDivisorCommonRec(n, resto);
    }
}

int verifyRepeatCuts(constraintsReal *constraintsOriginal, int cutOriginal, int cutCreate)
{
    int i, j, aux = 1;
    int szOri = constraintsOriginal->ElementsConstraints[cutOriginal + 1] - constraintsOriginal->ElementsConstraints[cutOriginal];
    int szCre = constraintsOriginal->ElementsConstraints[cutCreate + 1] - constraintsOriginal->ElementsConstraints[cutCreate];
    if ((szOri != szCre) || (constraintsOriginal->rightSide[cutOriginal] != constraintsOriginal->rightSide[cutCreate]))
    {
        return 0;
    }
    for (i = constraintsOriginal->ElementsConstraints[cutOriginal]; i < constraintsOriginal->ElementsConstraints[cutOriginal + 1]; i++)
    {
        aux = 0;
        for (j = constraintsOriginal->ElementsConstraints[cutCreate]; j < constraintsOriginal->ElementsConstraints[cutCreate + 1]; j++)
        {
            if ((constraintsOriginal->Coefficients[i] == constraintsOriginal->Coefficients[j]) && (constraintsOriginal->Elements[i] == constraintsOriginal->Elements[j]))
            {
                aux = 1;
            }
        }
        if (aux == 0)
        {
            return 0;
        }
    }
    return 1;
}


double fRand(double fMin, double fMax)
{
    struct timeval time;
    gettimeofday(&time, NULL);
    srand((time.tv_sec * 1000) + (time.tv_usec / 1000));
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}