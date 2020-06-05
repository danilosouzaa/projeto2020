

#include "constraintsManipulation.h"

constraintsReal *AllocStrConstraintsReal(TCont cont, TNumberConstraints nConstraints, TNumberVariables nVariables)
{

    constraintsReal *cut = (constraintsReal *)malloc(sizeof(constraintsReal));
    cut->Coefficients = (TCoefficientsFull *)malloc(sizeof(TCoefficientsFull) * cont);
    cut->Elements = (TElements *)malloc(sizeof(TElements) * cont);
    cut->ElementsConstraints = (TElementsConstraints *)malloc(sizeof(TElementsConstraints) * (nConstraints + 1));
    cut->rightSide = (TRightSideFull *)malloc(sizeof(TRightSideFull) * nConstraints);
    cut->xAsterisc = (TXAsteriscFull *)malloc(sizeof(TXAsteriscFull) * nVariables);
    cut->numberVariables = nVariables;
    cut->numberConstraints = nConstraints;
    cut->cont = cont;
    return cut;
}

cutSmall *AllocStrCutSmall(TCont cont, TNumberConstraints nConstraints, TNumberVariables nVariables)
{
    size_t size_cut = sizeof(cutSmall) +
                      sizeof(TCoefficients) * (cont) +
                      sizeof(TElements) * (cont) +
                      sizeof(TElementsConstraints) * (nConstraints + 1) +
                      sizeof(TRightSide) * (nConstraints) +
                      sizeof(TXAsterisc) * (nVariables) +
                      sizeof(TNumberConstraints) * (nConstraints);

    cutSmall *cut = (cutSmall *)malloc(size_cut);
    assert(cut != NULL);
    memset(cut, 0, size_cut);
    cut->Coefficients = (TCoefficients *)(cut + 1);
    cut->Elements = (TElements *)(cut->Coefficients + cont);
    cut->ElementsConstraints = (TElementsConstraints *)(cut->Elements + cont);
    cut->rightSide = (TRightSide *)(cut->ElementsConstraints + (nConstraints + 1));
    cut->xAsterisc = (TXAsterisc *)(cut->rightSide + (nConstraints));
    cut->originalConstraints = (TNumberConstraints *)(cut->xAsterisc + (nVariables));
    cut->numberVariables = nVariables;
    cut->numberConstraints = nConstraints;
    cut->cont = cont;
    return cut;
}


constraintMin *AllocStrConstraintsMin(TCont cont){
      size_t size_c= sizeof(constraintMin) +
                      sizeof(TCoefficients) * (cont) +
                      sizeof(TPosList) * (cont) +
                      sizeof(TActivedCut) * (cont);

    constraintMin *c = (constraintMin *)malloc(size_c);
    assert(c != NULL);
    memset(c, 0, size_c);
    c->Coefficients = (TCoefficients *)(c + 1);
    c->positionOriginal = (TPosList *)(c->Coefficients + cont);
    c->Cover = (TActivedCut *)(c->positionOriginal + cont);
    c->cont = cont;
    return c;
}

void freeStrConstraintsReal(constraintsReal *cut)
{
    free(cut->Coefficients);
    free(cut->Elements);
    free(cut->ElementsConstraints);
    free(cut->rightSide);
    free(cut->xAsterisc);
    free(cut);
}



constraintsReal *fillStructPerLP(LinearProgram *lp, char **nameConstraints, char **nameVariables, int *typeVariables, double *lbVariables, double *ubVariables)
{
    TNumberConstraints nConstraints, nConstraintsValided;
    TNumberVariables nVariables, nNonZeroValided = 0;
    nVariables = lp_cols(lp);
    nConstraints = lp_rows(lp);

    int szConst, i, j;
    for (i = 0; i < nVariables; i++)
    {
        typeVariables[i] = lp_is_binary(lp, i);
        lbVariables[i] = lp_col_lb(lp, i);
        ubVariables[i] = lp_col_ub(lp, i);
        //printf("Teste: %d %f %f \n", typeVariables[i],lbVariables[i], ubVariables[i]);
    }

    nConstraintsValided = 0;

    TElements *idx = (TElements *)malloc(sizeof(TElements) * nVariables);
    double *coef = (double *)malloc(sizeof(double) * nVariables);
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
                    nNonZeroValided++;
                }
            }
            if (szConst < 2)
            {
                nConstraintsValided--;
                nNonZeroValided -= szConst;
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
                    nNonZeroValided += 2;
                }
            }
            if (szConst < 2)
            {
                nConstraintsValided -= 2;
                nNonZeroValided -= 2 * szConst;
            }
        }
    }
    double rhs;
    TXAsteriscFull *xTemp;
    lp_set_print_messages(lp, 0);
    lp_optimize_as_continuous(lp);
    xTemp = lp_x(lp);
    //saveSoluctionFrac(xTemp, numberVariables,nameInst,lp,0);

    constraintsReal *h_cut = AllocStrConstraintsReal(nNonZeroValided, nConstraintsValided, nVariables);

    //TCoefficientsFull *v_aux = (TCoefficientsFull*)malloc(sizeof(TCoefficientsFull)*(nVariables + 1));

    for (i = 0; i < nVariables; i++)
    {
        coef[i] = 0.0;
        idx[i] = 0;
        h_cut->xAsterisc[i] = xTemp[i];
    }

    int aux = 0;
    h_cut->ElementsConstraints[0] = 0;
    int contador = 0;
    int a = 1;
    for (i = 0; i < nConstraints; i++)
    {
        if (lp_sense(lp, i) == 'L')
        {
            lp_row(lp, i, idx, coef);
            rhs = lp_rhs(lp, i);
            if (rhs == 1)
            {
                a = 2;
            }
            else
            {
                a = 1;
            }
            szConst = 0;
            for (j = 0; j < nVariables; j++)
            {
                if (coef[j] != 0.0)
                {
                    szConst++;
                }
            }
            if (szConst < 2)
            {
                continue;
            }
            for (j = 0; j < nVariables; j++)
            {
                if (coef[j] != 0.0)
                {

                    h_cut->Coefficients[aux] = a * coef[j];
                    h_cut->Elements[aux] = idx[j];
                    aux++;
                }

                coef[j] = 0.0;
                idx[j] = 0;
            }

            lp_row_name(lp, i, nameConstraints[contador]);

            h_cut->ElementsConstraints[contador + 1] = aux;
            h_cut->rightSide[contador] = a * rhs;
            contador++;
        }
        else if (lp_sense(lp, i) == 'G')
        {
            lp_row(lp, i, idx, coef);
            rhs = lp_rhs(lp, i);
            if (rhs == 1)
            {
                a = 2;
            }
            else
            {
                a = 1;
            }
            szConst = 0;
            for (j = 0; j < nVariables; j++)
            {
                if (coef[j] != 0.0)
                {
                    szConst++;
                }
            }
            if (szConst < 2)
            {
                continue;
            }

            for (j = 0; j < nVariables; j++)
            {
                if (coef[j] != 0.0)
                {
                    h_cut->Coefficients[aux] = -1 * coef[j] * a;
                    h_cut->Elements[aux] = idx[j];
                    aux++;
                }

                coef[j] = 0.0;
                idx[j] = 0;
            }

            lp_row_name(lp, i, nameConstraints[contador]);

            h_cut->ElementsConstraints[contador + 1] = aux;
            h_cut->rightSide[contador] = -1 * rhs * a;
            contador++;
        }
        else
        {

            lp_row(lp, i, idx, coef);
            rhs = lp_rhs(lp, i);
            if (rhs == 1)
            {
                a = 2;
            }
            else
            {
                a = 1;
            }
            szConst = 0;
            for (j = 0; j < nVariables; j++)
            {
                if (coef[j] != 0.0)
                {
                    szConst++;
                }
            }
            if (szConst < 2)
            {
                continue;
            }
            for (j = 0; j < nVariables; j++)
            {
                if (coef[j] != 0.0)
                {
                    h_cut->Coefficients[aux] = a * coef[j];
                    h_cut->Elements[aux] = idx[j];
                    aux++;
                }

                //coef[j] = 0.0;
                //idx[j] = 0;
            }
            char nameTemp[255];
            lp_row_name(lp, i, nameTemp);
            strcat(nameTemp, "_1");
            strcpy(nameConstraints[contador], nameTemp);
            h_cut->ElementsConstraints[contador + 1] = aux;
            h_cut->rightSide[contador] = a * rhs;
            contador++;

            for (j = 0; j < nVariables; j++)
            {
                if (coef[j] != 0.0)
                {
                    h_cut->Coefficients[aux] = -1 * coef[j] * a;
                    h_cut->Elements[aux] = idx[j];
                    aux++;
                }

                coef[j] = 0.0;
                idx[j] = 0;
            }
            strcpy(nameTemp, "");
            lp_row_name(lp, i, nameTemp);
            strcat(nameTemp, "_2");
            strcpy(nameConstraints[contador], nameTemp);
            h_cut->ElementsConstraints[contador + 1] = aux;
            h_cut->rightSide[contador] = -1 * rhs * a;
            contador++;
        }
    }

    //free(xTemp);
    free(idx);
    free(coef);
    return h_cut;
}

void showStructFull(constraintsReal *constraintsFull, char **nameConstraints, char **nameVariables)
{
    int i, j, el, flag;
    for (i = 0; i < constraintsFull->numberConstraints; i++)
    {
        flag = 0;
        printf("%s:\t ", nameConstraints[i]);
        for (j = constraintsFull->ElementsConstraints[i]; j < constraintsFull->ElementsConstraints[i + 1]; j++)
        {
            el = constraintsFull->Elements[j];
            if (constraintsFull->Coefficients[j] > 0)
            {
                if (flag == 0)
                {
                    printf("%.2f %s  = %lf ", constraintsFull->Coefficients[j], nameVariables[el], constraintsFull->xAsterisc[el]);
                    flag = 1;
                }
                else
                {
                    printf("+ %.2f %s  = %lf ", constraintsFull->Coefficients[j], nameVariables[el], constraintsFull->xAsterisc[el]);
                }
            }
            else
            {
                flag = 1;
                printf("%.2f %s = %lf ", constraintsFull->Coefficients[j], nameVariables[el], constraintsFull->xAsterisc[el]);
            }
        }
        printf("<= %.2f\n", constraintsFull->rightSide[i]);
    }
}

int verifyCutsValidatedPerSolutionInteger(constraintsReal *constraintsOriginal, int cut, double *sol, char **nameVariables)
{
    double lhs = 0;
    int j, el;
    for (j = constraintsOriginal->ElementsConstraints[cut]; j < constraintsOriginal->ElementsConstraints[cut + 1]; j++)
    {
        el = constraintsOriginal->Elements[j];
        lhs += constraintsOriginal->Coefficients[j] * (sol[el]);
    }
    if (lhs - 1e-5 <= constraintsOriginal->rightSide[cut])
    {
        return 1;
    }
    else
    {
        //printf("lhs: %lf rhs: %lf\n", lhs, constraintsOriginal->rightSide[cut]);
        return 0;
    }
}

int *returnBinaryConstraints(constraintsReal *constraintsOriginal, int *typeVariables)
{
    int *BinaryConstraints = (int *)malloc(sizeof(int) * constraintsOriginal->numberConstraints);
    int i, j, qntBin, qntNBin, el;
    for (i = 0; i < constraintsOriginal->numberConstraints; i++)
    {
        qntBin = 0;
        qntNBin = 0;
        for (j = constraintsOriginal->ElementsConstraints[i]; j < constraintsOriginal->ElementsConstraints[i + 1]; j++)
        {
            el = constraintsOriginal->Elements[j];
            //printf("%d \t", typeVariables[el]);
            if (typeVariables[el] == 0)
            {
                qntNBin++;
            }
            else
            {
                qntBin++;
            }
        }
        if (qntBin < 1)
        {
            BinaryConstraints[i] = 0; // restrição que não vale a pena transformar
        }
        else if (qntNBin == 0)
        {
            BinaryConstraints[i] = 1; // restrição com todos os elementos binários
        }
        else
        {
            BinaryConstraints[i] = 2; // restrição vale a pena transformar em binária
        }
    }
    return BinaryConstraints;
}

TNumberConstraints countConstraintsBinaryUsed(int *binaryConstraints, int sz)
{
    int constraintsUsed = 0, i;
    for (i = 0; i < sz; i++)
    {
        if (binaryConstraints[i] != 0)
        {
            constraintsUsed++;
        }
    }
    return constraintsUsed;
}

constraintsReal *convertBinaryConstraints(constraintsReal *constraintsOriginal, int *BinaryConstraints, int *typeVariables, double *lbVariables, double *ubVariables)
{
    int i, j, el;
    int contConstraints, newCont, contAux;
    ;
    int *constraintsAvalible = (int *)malloc(sizeof(int) * constraintsOriginal->numberConstraints);

    contConstraints = 0;
    newCont = 0;
    for (i = 0; i < constraintsOriginal->numberConstraints; i++)
    {
        constraintsAvalible[i] = 0;
        if (BinaryConstraints[i] != 0)
        {
            contAux = 0;
            for (j = constraintsOriginal->ElementsConstraints[i]; j < constraintsOriginal->ElementsConstraints[i + 1]; j++)
            {
                el = constraintsOriginal->Elements[j];
                if (typeVariables[el] != 1)
                {
                    if (((constraintsOriginal->Coefficients[j] < 0) && (ubVariables[el] == DBL_MAX)) || ((constraintsOriginal->Coefficients[j] > 0) && (lbVariables[el] == -DBL_MAX)))
                    {
                        contAux = 0;
                        break;
                    }
                }
                else
                {
                    contAux++;
                }
            }
            if (contAux >= 2)
            {
                newCont += contAux;
                contConstraints++;
                constraintsAvalible[i] = 1;
            }
        }
    }
    constraintsReal *constraintsBinary = AllocStrConstraintsReal(newCont, contConstraints, constraintsOriginal->numberVariables);
    newCont = 0;
    contConstraints = 0;
    constraintsBinary->ElementsConstraints[0] = 0;
    int rhsTemp = 0;
    for (i = 0; i < constraintsOriginal->numberConstraints; i++)
    {
        if (constraintsAvalible[i] == 1)
        {
            rhsTemp = 0;
            for (j = constraintsOriginal->ElementsConstraints[i]; j < constraintsOriginal->ElementsConstraints[i + 1]; j++)
            {
                el = constraintsOriginal->Elements[j];
                if (typeVariables[el] == 1)
                {
                    constraintsBinary->Coefficients[newCont] = constraintsOriginal->Coefficients[j];
                    constraintsBinary->Elements[newCont] = el;
                    newCont++;
                }
                else
                {
                    if (constraintsOriginal->Coefficients[j] < 0)
                    {
                        if ((ubVariables[el] >= 0) && (lbVariables[el] >= 0))
                        {
                            rhsTemp += (constraintsOriginal->Coefficients[j] * -1) * ubVariables[el];
                        }
                        else if ((ubVariables[el] <= 0) && (lbVariables[el] <= 0))
                        {
                            rhsTemp += 0;
                        }
                        else
                        {
                            rhsTemp += (constraintsOriginal->Coefficients[j] * -1) * ubVariables[el];
                        }
                    }
                    else
                    {
                        if ((ubVariables[el] >= 0) && (lbVariables[el] >= 0))
                        {
                            rhsTemp += 0;
                        }
                        else if ((ubVariables[el] <= 0) && (lbVariables[el] <= 0))
                        {
                            rhsTemp += constraintsOriginal->Coefficients[j] * (lbVariables[el] * (-1));
                        }
                        else
                        {
                            rhsTemp += constraintsOriginal->Coefficients[j] * (lbVariables[el] * -1);
                        }
                    }
                }
            }
            constraintsBinary->ElementsConstraints[contConstraints + 1] = newCont;
            constraintsBinary->rightSide[contConstraints] = constraintsOriginal->rightSide[i] + rhsTemp;
            contConstraints++;
        }
    }
    for (i = 0; i < constraintsOriginal->numberVariables; i++)
    {
        constraintsBinary->xAsterisc[i] = constraintsOriginal->xAsterisc[i];
    }
    free(constraintsAvalible);
    return constraintsBinary;
}

constraintsReal *removeNegativeCoefficientsAndSort(constraintsReal *constraintsOriginal, int *convertVector)
{
    int i, j;
    int qntX = constraintsOriginal->numberVariables;
    int qntNegative = 0;
    TRightSideFull rhs = 0;
    int el = 0;
    for (i = 0; i < constraintsOriginal->cont; i++)
    {
        if (constraintsOriginal->Coefficients[i] < 0)
        {
            qntNegative++;
        }
    }
    qntNegative += constraintsOriginal->numberVariables;
    constraintsReal *newConstraints = AllocStrConstraintsReal(constraintsOriginal->cont, constraintsOriginal->numberConstraints, qntNegative);

    for (i = 0; i < constraintsOriginal->numberVariables; i++)
    {
        newConstraints->xAsterisc[i] = constraintsOriginal->xAsterisc[i];
    }

    newConstraints->ElementsConstraints[0] = 0;
    for (i = 0; i < constraintsOriginal->numberConstraints; i++)
    {
        rhs = constraintsOriginal->rightSide[i];
        for (j = constraintsOriginal->ElementsConstraints[i]; j < constraintsOriginal->ElementsConstraints[i + 1]; j++)
        {
            if (constraintsOriginal->Coefficients[j] < 0)
            {
                newConstraints->Coefficients[j] = constraintsOriginal->Coefficients[j] * (-1);
                rhs += newConstraints->Coefficients[j];

                newConstraints->Elements[j] = qntX;
                el = constraintsOriginal->Elements[j];
                newConstraints->xAsterisc[qntX] = 1 - constraintsOriginal->xAsterisc[el];
                convertVector[qntX - constraintsOriginal->numberVariables] = constraintsOriginal->Elements[j];
                qntX++;
            }
            else
            {
                newConstraints->Coefficients[j] = constraintsOriginal->Coefficients[j];
                newConstraints->Elements[j] = constraintsOriginal->Elements[j];
            }
        }
        newConstraints->rightSide[i] = rhs;
        newConstraints->ElementsConstraints[i + 1] = constraintsOriginal->ElementsConstraints[i + 1];
    }

    SortByCoefficients(newConstraints);
    freeStrConstraintsReal(constraintsOriginal);
    return newConstraints;
}

void SortByCoefficients(constraintsReal *h_cut)
{
    int i = 0;
    for (i = 0; i < h_cut->numberConstraints; i++)
    {
        quicksortCof(h_cut->Coefficients, h_cut->Elements, h_cut->ElementsConstraints[i], h_cut->ElementsConstraints[i + 1]);
    }
}

void quicksortCof(TCoefficientsFull *values, int *idc, int began, int end)
{
    int i, j;
    TCoefficientsFull pivo, aux;
    i = began;
    j = end - 1;
    pivo = values[(began + end) / 2];
    while (i <= j)
    {
        while (values[i] > pivo && i < end)
        {
            i++;
        }
        while (values[j] < pivo && j > began)
        {
            j--;
        }
        if (i <= j)
        {
            aux = values[i];
            values[i] = values[j];
            values[j] = aux;

            aux = idc[i];
            idc[i] = idc[j];
            idc[j] = aux;

            i++;
            j--;
        }
    }
    if (j > began)
        quicksortCof(values, idc, began, j + 1);
    if (i < end)
        quicksortCof(values, idc, i, end);
}

void quicksortDouble(double *values, int began, int end)
{
    int i, j;
    double pivo, aux;
    i = began;
    j = end - 1;
    pivo = values[(began + end) / 2];
    while (i <= j)
    {
        while (values[i] > pivo && i < end)
        {
            i++;
        }
        while (values[j] < pivo && j > began)
        {
            j--;
        }
        if (i <= j)
        {
            aux = values[i];
            values[i] = values[j];
            values[j] = aux;
            i++;
            j--;
        }
    }
    if (j > began)
        quicksortDouble(values, began, j + 1);
    if (i < end)
        quicksortDouble(values, i, end);
}

int verifyOfFloatIsInteger(TCoefficientsFull coef)
{
    //printf("%f %f", ceil(coef), coef);
    if (ceil(coef) == coef)
    {
        return 1;
    }

    return 0;
}

int *returnVectorTypeContraintsIntOrFloat(constraintsReal *constraints)
{
    int flag = 1, i, j, aux = 0;
    TCoefficientsFull coeffAux, lhsCoefs = 0;
    int *vectorTypeInt = (int *)malloc(sizeof(int) * constraints->numberConstraints);
    for (i = 0; i < constraints->numberConstraints; i++)
    {
        flag = 1;
        lhsCoefs = 0;
        for (j = constraints->ElementsConstraints[i]; j < constraints->ElementsConstraints[i + 1]; j++)
        {
            coeffAux = constraints->Coefficients[j];
            lhsCoefs += constraints->Coefficients[j];
            aux = verifyOfFloatIsInteger(coeffAux);
            if (aux == 0)
            {
                flag = 0;
                break;
            }
        }
        aux = verifyOfFloatIsInteger(constraints->rightSide[i]);
        if ((aux == 0) && (lhsCoefs <= constraints->rightSide[i]))
        {
            flag = 0;
        }
        if (flag)
        {
            vectorTypeInt[i] = 1;
        }
        else
        {
            vectorTypeInt[i] = 0;
        }
    }
    return vectorTypeInt;
}

cutCover *CopyCutToCover(cutSmall *h_cut)
{
    //Cover_gpu *cover = AllocationStructCover(h_cut->cont,h_cut->numberConstrains);
    cutCover *cover = AllocStrCover(h_cut->cont, h_cut->numberConstraints);
    int i;
    for (i = 0; i < h_cut->cont; i++) //VERIFY WHY H_CUT->CONT AND NOT NCONSTRAINTSINI
    {
        cover->Coefficients[i] = h_cut->Coefficients[i];
    }
    for (i = 0; i < h_cut->numberConstraints; i++)
    {
        cover->ElementsConstraints[i] = h_cut->ElementsConstraints[i];
        cover->rightSide[i] = h_cut->rightSide[i];
    }
    cover->ElementsConstraints[i] = h_cut->ElementsConstraints[i];
    return cover;
}

cutCover *AllocStrCover(TCont cont, TNumberConstraints nConstraints)
{
    size_t size_cover = sizeof(cutCover) +
                        sizeof(TCoefficients) * (cont) +
                        sizeof(TElementsConstraints) * (nConstraints + 1) +
                        sizeof(TRightSide) * (nConstraints);
    cutCover *h_cover = (cutCover *)malloc(size_cover);
    assert(h_cover != NULL);
    memset(h_cover, 0, size_cover);
    h_cover->Coefficients = (TCoefficients *)(h_cover + 1);
    h_cover->ElementsConstraints = (TElementsConstraints *)(h_cover->Coefficients + cont);
    h_cover->rightSide = (TRightSide *)(h_cover->ElementsConstraints + (nConstraints + 1));
    h_cover->numberConstraints = nConstraints;
    h_cover->cont = cont;
    return h_cover;
}

int verifyRepeated(constraintsReal *originalConstraints, int posCover)
{
    int i, j, k, cont = 0;
    int sz = originalConstraints->ElementsConstraints[posCover + 1] - originalConstraints->ElementsConstraints[posCover];
    for (i = 0; i < posCover; i++)
    {
        for (j = originalConstraints->ElementsConstraints[posCover]; j < originalConstraints->ElementsConstraints[posCover + 1]; j++)
        {
            if ((sz == (originalConstraints->ElementsConstraints[i + 1] - originalConstraints->ElementsConstraints[i])) && (originalConstraints->rightSide[posCover] == originalConstraints->rightSide[i]))
            {
                for (k = originalConstraints->ElementsConstraints[i]; k < originalConstraints->ElementsConstraints[i + 1]; k++)
                {
                    if ((originalConstraints->Elements[j] == originalConstraints->Elements[k]) && (originalConstraints->Coefficients[j] == originalConstraints->Coefficients[k]))
                    {
                        cont++;
                        break;
                    }
                }
            }
        }
        if (cont == sz)
        {
            return 0;
        }
        else
        {
            cont = 0;
        }
    }
    return 1;
}

int verifyRepeatedIncidency(int **matrizIncidencia, constraintsReal *originalConstraints, int posCover)
{
    int i, j, cont = 0, el, aux;
    int sz = originalConstraints->ElementsConstraints[posCover + 1] - originalConstraints->ElementsConstraints[posCover];
    for (i = 0; i < posCover; i++)
    {
        if ((sz == (originalConstraints->ElementsConstraints[i + 1] - originalConstraints->ElementsConstraints[i])) && (originalConstraints->rightSide[posCover] == originalConstraints->rightSide[i]))
        {
            for (j = originalConstraints->ElementsConstraints[posCover]; j < originalConstraints->ElementsConstraints[posCover + 1]; j++)
            {
                el = originalConstraints->Elements[j];
                if (matrizIncidencia[i][el] != -1)
                {
                    aux = matrizIncidencia[i][el];
                    if (originalConstraints->Coefficients[j] == originalConstraints->Coefficients[aux])
                    {
                        cont++;
                    }
                }
                else
                {
                    break;
                }
            }
        }
        if (cont == sz)
        {
            return 0;
        }
        else
        {
            cont = 0;
        }
    }
    return 1;
}

double valueViolation(cutCover *cCover, cutSmall *constraintsSmall, TNumberConstraints idCover, TNumberConstraints ogConstraint, int precision)
{
    int i = 0, el;
    double lhs = 0, violation = 0;
    int aux = constraintsSmall->ElementsConstraints[ogConstraint];
    for (i = cCover->ElementsConstraints[idCover]; i < cCover->ElementsConstraints[idCover + 1]; i++)
    {
        el = constraintsSmall->Elements[aux];
        lhs += ((double)cCover->Coefficients[i]) * ((double)constraintsSmall->xAsterisc[el] / (double)precision);
        //       printf("x* %d= %d coef = %d \n",el, constraintsSmall->xAsterisc[el], constraintsSmall->Coefficients[i]);
        aux++;
        //    }
    }
    // printf("lhs %f rhs %d\n",lhs, cCover->rightSide[idCover]);
    violation = lhs - (double)cCover->rightSide[idCover];
    if (violation + 1e-5 < 0.0)
    {
        violation = 0.0;
    }
    return violation;
}

void quicksortTParameters(TParametersGroup *values, TNumberConstraints *idc, TNumberConstraints began, TNumberConstraints end)
{
    TNumberConstraints i, j;
    TParametersGroup pivo, aux;
    TNumberConstraints auxIdc;
    i = began;
    j = end - 1;
    pivo = values[(began + end) / 2];
    while (i <= j)
    {
        while (values[i] > pivo && i < end)
        {
            i++;
        }
        while (values[j] < pivo && j > began)
        {
            j--;
        }
        if (i <= j)
        {
            aux = values[i];
            values[i] = values[j];
            values[j] = aux;

            auxIdc = idc[i];
            idc[i] = idc[j];
            idc[j] = auxIdc;

            i++;
            j--;
        }
    }
    if (j > began)
        quicksortTParameters(values, idc, began, j + 1);
    if (i < end)
        quicksortTParameters(values, idc, i, end);
}

constraintsReal *returnVariablesOriginals(constraintsReal *constraintsOriginal, int *convertVector, int precision, int nVariablesInitial)
{
    //printf("%d - nVariablesInitial\n", nVariablesInitial);
    constraintsReal *newConstraints = AllocStrConstraintsReal(constraintsOriginal->cont, constraintsOriginal->numberConstraints, nVariablesInitial);
    int i, j, el;
    TRightSideFull rhs;
    for (i = 0; i < nVariablesInitial; i++)
    {
        newConstraints->xAsterisc[i] = constraintsOriginal->xAsterisc[i];
    }
    for (i = 0; i < constraintsOriginal->numberConstraints; i++)
    {
        rhs = constraintsOriginal->rightSide[i];
        for (j = constraintsOriginal->ElementsConstraints[i]; j < constraintsOriginal->ElementsConstraints[i + 1]; j++)
        {
            el = constraintsOriginal->Elements[j];
            //printf("valor de el: %d\n",el);
            if (el >= nVariablesInitial)
            {
                newConstraints->Coefficients[j] = constraintsOriginal->Coefficients[j] * (-1);
                rhs -= constraintsOriginal->Coefficients[j];
                newConstraints->Elements[j] = convertVector[el - nVariablesInitial];
                newConstraints->xAsterisc[newConstraints->Elements[j]] = 1 - constraintsOriginal->xAsterisc[el]; //erro aqui
            }
            else
            {
                newConstraints->Coefficients[j] = constraintsOriginal->Coefficients[j];
                newConstraints->Elements[j] = constraintsOriginal->Elements[j];
            }
        }

        newConstraints->rightSide[i] = rhs;
        newConstraints->ElementsConstraints[i] = constraintsOriginal->ElementsConstraints[i];
    }
    newConstraints->ElementsConstraints[i] = constraintsOriginal->ElementsConstraints[i];
    freeStrConstraintsReal(constraintsOriginal);
    return newConstraints;
}

constraintsReal *convertBinaryOfOriginalConstraints(constraintsReal *constraintsOriginal, constraintsReal *constraintsBinary, int nInitialBinary)
{
    TNumberConstraints sizeNew = constraintsOriginal->numberConstraints + (constraintsBinary->numberConstraints - nInitialBinary);
    TCont contNew = constraintsOriginal->cont;
    int i, j, jTemp;
    ;
    for (i = nInitialBinary; i < constraintsBinary->numberConstraints; i++)
    {
        contNew += constraintsBinary->ElementsConstraints[i + 1] - constraintsBinary->ElementsConstraints[i];
    }
    if (constraintsBinary->numberConstraints - nInitialBinary == 0)
    {
        return constraintsOriginal;
    }
    constraintsReal *newConstraintsOriginal = AllocStrConstraintsReal(contNew, sizeNew, constraintsOriginal->numberVariables);

    for (j = 0; j < constraintsOriginal->numberVariables; j++)
    {
        newConstraintsOriginal->xAsterisc[j] = constraintsOriginal->xAsterisc[j];
    }
    newConstraintsOriginal->ElementsConstraints[0] = constraintsOriginal->ElementsConstraints[0];
    for (i = 0; i < constraintsOriginal->numberConstraints; i++)
    {
        newConstraintsOriginal->ElementsConstraints[i + 1] = constraintsOriginal->ElementsConstraints[i + 1];
        newConstraintsOriginal->rightSide[i] = constraintsOriginal->rightSide[i];
        for (j = constraintsOriginal->ElementsConstraints[i]; j < constraintsOriginal->ElementsConstraints[i + 1]; j++)
        {
            newConstraintsOriginal->Coefficients[j] = constraintsOriginal->Coefficients[j];
            newConstraintsOriginal->Elements[j] = constraintsOriginal->Elements[j];
        }
    }
    contNew = newConstraintsOriginal->ElementsConstraints[i];
    jTemp = constraintsOriginal->numberConstraints;
    for (i = nInitialBinary; i < constraintsBinary->numberConstraints; i++)
    {
        for (j = constraintsBinary->ElementsConstraints[i]; j < constraintsBinary->ElementsConstraints[i + 1]; j++)
        {
            newConstraintsOriginal->Coefficients[contNew] = constraintsBinary->Coefficients[j];
            newConstraintsOriginal->Elements[contNew] = constraintsBinary->Elements[j];
            contNew++;
        }
        newConstraintsOriginal->ElementsConstraints[jTemp + 1] = contNew;
        newConstraintsOriginal->rightSide[jTemp] = constraintsBinary->rightSide[i];
        jTemp++;
    }

    freeStrConstraintsReal(constraintsOriginal);
    return newConstraintsOriginal;
}