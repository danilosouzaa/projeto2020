#include "newGrasp.h"

constraintsReal *graspLciAdam(constraintsReal *constraintsFull, int precision, char **nameConstraints, char **nameVariables, int szPoolCutsMax, int nIterationGrasp, float alpha, int minimal)
{
    int i, j, k, itePool;
    int iteConstraints = 0, iteSzPool = 0;
    cutCover *coverCuts;
    int c_AuxSolution = 0;
    int c_XSolution = 0;
    int *intOrFloat = returnVectorTypeContraintsIntOrFloat(constraintsFull); //verifica se as restrições podem ser transformadas em restrições da mochila
    cutSmall *knapsackConstraints = reduceCutFullForCutSmall(constraintsFull, intOrFloat, precision);
    coverCuts = CopyCutToCover(knapsackConstraints);
    for (i = 0; i < knapsackConstraints->numberConstraints; i++)
    {
        int lhs = 0;
        if (coverCuts->rightSide[i] <= 1)
        {
            continue;
        }
        // printf("TESTE: %d\n!!", i);
        c_AuxSolution = 0;
        c_XSolution = 0;
        int szConstraint = coverCuts->ElementsConstraints[i + 1] - coverCuts->ElementsConstraints[i];
        cutCover *cutsCoverSolution = AllocStrCover(szConstraint * szPoolCutsMax, szPoolCutsMax);
        int nPoolUsed = 0;
        int *poolSolution = coverLCIAdamPerConstraints(knapsackConstraints, precision, i, nIterationGrasp, szPoolCutsMax, alpha, minimal, &nPoolUsed);
        cutsCoverSolution->ElementsConstraints[0] = 0;
        //printf("Pool used:%d\n", nPoolUsed);
        for (itePool = 0; itePool < nPoolUsed; itePool++)
        {
            //qnt = 0;
            int caux = 0;
            lhs = 0;
            int *solutionTemp = (int *)malloc(sizeof(int) * szConstraint);
            for (k = coverCuts->ElementsConstraints[i]; k < coverCuts->ElementsConstraints[i + 1]; k++)
            {
                //  qnt += poolSolution[caux + itePool * szConstraint];

                solutionTemp[k - coverCuts->ElementsConstraints[i]] = poolSolution[(k - coverCuts->ElementsConstraints[i]) + itePool * szConstraint];
                lhs += coverCuts->Coefficients[k] * poolSolution[caux + itePool * szConstraint];
                caux++;
            }
            if (lhs <= coverCuts->rightSide[i])
            {
                free(solutionTemp);
                continue;
            }

            int *cutCoverLifted = LCIAdam(solutionTemp, coverCuts, i);
            for (k = 0; k < szConstraint; k++)
            {
                cutsCoverSolution->Coefficients[c_XSolution] = cutCoverLifted[k];
                c_XSolution++;
            }
            cutsCoverSolution->ElementsConstraints[c_AuxSolution + 1] = c_XSolution;
            cutsCoverSolution->rightSide[c_AuxSolution] = cutCoverLifted[szConstraint];
            c_AuxSolution++;
            free(cutCoverLifted);

            free(solutionTemp);
        }
        int *idc_cover = (int *)malloc(sizeof(int) * c_AuxSolution);
        constraintsFull = createCutsCoverGrasp(cutsCoverSolution, constraintsFull, knapsackConstraints, idc_cover, i, c_AuxSolution, precision);
        free(cutsCoverSolution);
        free(idc_cover);
        free(poolSolution);
    }
    free(coverCuts);
    free(knapsackConstraints);
    free(intOrFloat);
    return constraintsFull;
}

int *coverLCIAdamPerConstraints(cutSmall *knapsackConstraints, int precision, int constraint, int nIterationGrasp, int szPoolCutsMax, float alpha, int minimal, int *nPoolUsed)
{
    int sz = knapsackConstraints->ElementsConstraints[constraint + 1] - knapsackConstraints->ElementsConstraints[constraint];
    int *solution = (int *)malloc(sizeof(int) * sz);
    int *poolSolution = (int *)malloc(sizeof(int) * sz * szPoolCutsMax);
    int nPoolSolution = 0;
    int ite = 0;
    double testAlpha = alpha;
    int i = 0;
    memset(poolSolution, 0, sizeof(int) * sz * szPoolCutsMax);
    // for (ite = 0; ite < sz * szPoolCutsMax; ite++)
    // {
    //     poolSolution[ite] = 0;
    // }
    if (testAlpha == -1)
    {
        //printf("entre 0.0 e 0.1");
        alpha = fRand(0.0, 0.1);
        //printf("alpha: %f\n", alpha);
    }

    ite = 0;
    // if (numberIteration == 0)
    // {
    //     copyAndVerifyPoolSolution(solution, sz, poolSolution, &nPoolSolution);
    // }
    constraintMin *cMin = remove_X_0_and_1_one_Constraints(knapsackConstraints, constraint, precision);
    if(cMin->Cover[0]==-1){
        free(solution);
        free(cMin);
        return poolSolution;
    }
    //printf("cMIn: %d %d \n",cMin->Cover[0], cMin->cont);
    while (ite <= nIterationGrasp)
    {
        int *solutionTemp = (int *)malloc(sizeof(int) * cMin->cont);
       // memset(solution,0,sizeof(int)*sz);
        initialLCIAdamGRASP(solutionTemp, cMin->cont, cMin, knapsackConstraints, precision, constraint, alpha);
        // for(i = 0;i< cMin->cont;i++){
        //     printf("x = %d \t", solutionTemp[i]);
        // }
        // getchar();
        copyPosGraspCMin(solution, solutionTemp, cMin, knapsackConstraints, constraint, precision);
        if (verifySolutionCoverMinimal(solution, knapsackConstraints, constraint) == 1)
        {
            int *solFinalTemp;
            cutCover *coverTemp = CopyCutToCover(knapsackConstraints);
            solFinalTemp = LCIAdam(solution, coverTemp, constraint);
            if (verifyViolationGreedy(solFinalTemp, knapsackConstraints, constraint, precision) == 1)
            {

                copyAndVerifyPoolSolution(solution, sz, poolSolution, &nPoolSolution);
                //printf("teste aqui");
                (*nPoolUsed)++;
            }

            free(solFinalTemp);
            free(coverTemp);
        }
        else
        {
            if (testAlpha == -1)
            {
                alpha = fRand(0.0, 0.1);
                //printf("alpha: %f\n", alpha);
            }
            ite++;
            free(solutionTemp);
            continue;
        }

        solution = localSearch(solution, sz, knapsackConstraints, precision, constraint, minimal);
        if (verifySolutionCoverMinimal(solution, knapsackConstraints, constraint) == 1)
        {
            int *solFinalTemp;
            cutCover *coverTemp = CopyCutToCover(knapsackConstraints);
            solFinalTemp = LCIAdam(solution, coverTemp, constraint);
            if (verifyViolationGreedy(solFinalTemp, knapsackConstraints, constraint, precision) == 1)
            {
                copyAndVerifyPoolSolution(solution, sz, poolSolution, &nPoolSolution);
                (*nPoolUsed)++;
            }
            free(coverTemp);
            free(solFinalTemp);
        }
        if (nPoolSolution == szPoolCutsMax)
        {
            break;
        }
        ite++;
        if (testAlpha == -1)
        {
            alpha = fRand(0.0, 0.1);
            //printf("alpha: %f\n", alpha);
        }
        free(solutionTemp);
    }
    free(cMin);
    free(solution);
    return poolSolution;
}

void initialLCIAdamGRASP(int *solution, int sz, constraintMin *cMin, cutSmall *constraintsSmall, int precision, int constraint, float alpha)
{
    int verifyCompleteSolution = 0;
    int *setId = (int *)malloc(sz * sizeof(int));
    double c_min = 0, c_max = 0, value, value_best;
    int i, aux = 0, el, szAux, lhs;
    struct timeval time;
    gettimeofday(&time, NULL);
    srand((time.tv_sec * 1000) + (time.tv_usec / 1000));
    for (i = 0; i < sz; i++)
    {
        cMin->Cover[i] = 0;
        solution[i] = 0;
    }
    for (i = 0; i < cMin->cont; i++)
    {
        setId[i] = i;
    }
    szAux = cMin->cont;
    lhs = 0;
    int test = 0;
    int contSizeNewSolution = 0;
    while (verifyCompleteSolution == 0)
    {
        aux = 0;
        for (i = 0; i < szAux; i++)
        {
            el = constraintsSmall->Elements[cMin->positionOriginal[setId[i]]];
            value = (double)constraintsSmall->xAsterisc[el] / (double)cMin->Coefficients[setId[i]];
            if (aux == 0)
            {
                c_min = value;
                c_max = value;
                aux = 1;
            }
            if (value < c_min)
            {
                c_min = value;
            }
            if (value > c_max)
            {
                c_max = value;
            }
        }
        value_best = c_min - alpha * (c_min - c_max);
        test++;
        int *setTemp = (int *)malloc(szAux * sizeof(int));
        int *posTemp = (int *)malloc(szAux * sizeof(int));
        int aux_t = 0;
        for (i = 0; i < szAux; i++)
        {
            el = constraintsSmall->Elements[cMin->positionOriginal[setId[i]]];
            value = (double)constraintsSmall->xAsterisc[el] / (double)cMin->Coefficients[setId[i]];
            if (value >= value_best)
            {
                setTemp[aux_t] = setId[i];
                posTemp[aux_t] = i;
                aux_t++;
            }
        }
        int itemAdd = rand() % aux_t;
        el = setTemp[itemAdd]; // - constraintsSmall->ElementsConstraints[constraint];
        solution[el] = 1;
        //cMin->Cover[el] = 1;
        contSizeNewSolution++;
        szAux--;
        lhs += cMin->Coefficients[setTemp[itemAdd]];
        if (lhs - 1e-5 > cMin->rightSide)
        {
            verifyCompleteSolution = 1;
            free(setTemp);
            free(posTemp);
            break;
        }

        free(setTemp);
        int *newSet = (int *)malloc(sizeof(int) * szAux);
        aux_t = 0;
        for (i = 0; i < szAux + 1; i++)
        {
            if (i != posTemp[itemAdd])
            {
                newSet[aux_t] = setId[i];
                aux_t++;
            }
        }
        free(setId);
        setId = newSet;
        free(posTemp);
        if (contSizeNewSolution == sz)
        {
            verifyCompleteSolution = 1;
        }
    }
    free(setId);
}
/*
constraintMin *remove_X_0_and_1_one_Constraints(cutSmall *knapsackConstraints, int constraint, int precision)
{
    int i = constraint, j, el, cont = 0, contAux = 0;
    int flag = 0;
    for (j = knapsackConstraints->ElementsConstraints[constraint]; j < knapsackConstraints->ElementsConstraints[constraint + 1]; j++)
    {
        el = knapsackConstraints->Elements[j];
        if ((knapsackConstraints->xAsterisc[el] != 1 * precision) && (knapsackConstraints->xAsterisc[el] != 0))
        {
            contAux++;
        }
    }
    if (contAux < 2)
    {
        flag = 1;
    }
    if (flag == 0)
    {
        constraintMin *cTemp = AllocStrConstraintsMin(contAux);
        TRightSide rhsTemp = knapsackConstraints->rightSide[constraint];
        cont = 0;
        for (j = knapsackConstraints->ElementsConstraints[i]; j < knapsackConstraints->ElementsConstraints[i + 1]; j++)
        {
            el = knapsackConstraints->Elements[j];
            if ((knapsackConstraints->xAsterisc[el] != 1 * precision) && (knapsackConstraints->xAsterisc[el] != 0))
            {
                cTemp->Coefficients[cont] = knapsackConstraints->Coefficients[j];
                cTemp->positionOriginal[cont] = j;
                cTemp->Cover[cont] = 0;
                cont++;
            }
            else if (knapsackConstraints->xAsterisc[el] == 1 * precision)
            {
                rhsTemp -= knapsackConstraints->Coefficients[j];
            }
        }
        cTemp->rightSide = rhsTemp;
        return cTemp;
    }
    else
    {
        cont = knapsackConstraints->ElementsConstraints[constraint + 1] - knapsackConstraints->ElementsConstraints[constraint];
        constraintMin *cTemp = AllocStrConstraintsMin(cont);
        cont = 0;
        for (j = knapsackConstraints->ElementsConstraints[i]; j < knapsackConstraints->ElementsConstraints[i + 1]; j++)
        {
            cTemp->Coefficients[cont] = knapsackConstraints->Coefficients[j];
            cTemp->positionOriginal[cont] = j;
            cTemp->Cover[cont] = -1;
            cont++;
        }
        cTemp->rightSide = knapsackConstraints->rightSide[constraint];
        return cTemp;
    }
}
*/
constraintMin *remove_X_0_and_1_one_Constraints(cutSmall *knapsackConstraints, int constraint, int precision)
{
    int i = constraint, j, el, cont = 0, contAux = 0;
    int sz = knapsackConstraints->ElementsConstraints[i + 1] - knapsackConstraints->ElementsConstraints[i];
    constraintMin *cTemp = AllocStrConstraintsMin(sz);
    TRightSide rhsTemp = knapsackConstraints->rightSide[constraint];
    cont = 0;
    int lhs = 0;
    for (j = knapsackConstraints->ElementsConstraints[i]; j < knapsackConstraints->ElementsConstraints[i + 1]; j++)
    {
        el = knapsackConstraints->Elements[j];
        //if ((knapsackConstraints->xAsterisc[el] != 1 * precision)&&(knapsackConstraints->xAsterisc[el] != 0))
        if (knapsackConstraints->xAsterisc[el] != 0)
        {
            cTemp->Coefficients[cont] = knapsackConstraints->Coefficients[j];
            lhs+=cTemp->Coefficients[cont];
            cTemp->positionOriginal[cont] = j;
            cTemp->Cover[cont] = 0;
            cont++;
        }
        //else if (knapsackConstraints->xAsterisc[el] == 1 * precision)
        //{
        //    rhsTemp -= knapsackConstraints->Coefficients[j];
        //}
    }
    cTemp->rightSide = rhsTemp;
    cTemp->cont = cont;
    if(lhs<cTemp->rightSide){
        cTemp->Cover[0]=-1;
    }
    //if(cont<2){
        //printf("AQUI!");
    //    cTemp->Cover[0]=-1;
    //}
    return cTemp;
}


void copyPosGraspCMin(int *solution, int *solutionTemp, constraintMin *cMin, cutSmall *knapsackConstraints, int constraint, int precision)
{
    int el, i = 0;
    // int sz = knapsackConstraints->ElementsConstraints[constraint + 1] - knapsackConstraints->ElementsConstraints[constraint];
    // if (cMin->Cover[0] != -1)
    // {
         for (i = knapsackConstraints->ElementsConstraints[constraint]; i < knapsackConstraints->ElementsConstraints[constraint + 1]; i++)
         {
    //         el = knapsackConstraints->Elements[i];
    //         if (knapsackConstraints->xAsterisc[el] == 1 * precision)
    //         {
    //            solution[i - knapsackConstraints->ElementsConstraints[constraint]] = 1;
    //        }
    //        else
    //        {
                 solution[i - knapsackConstraints->ElementsConstraints[constraint]] = 0;
    //        }
         }
     //}
    for (i = 0; i < cMin->cont; i++)
    {
        el = cMin->positionOriginal[i] - knapsackConstraints->ElementsConstraints[constraint];
        solution[el] = solutionTemp[i];
    }
}