#include "newGrasp.h"

constraintsReal *graspLciAdam(constraintsReal *constraintsFull, int precision, char **nameConstraints, char **nameVariables, int szPoolCutsMax, int nIterationGrasp, float alpha, int minimal)
{
    int i, j;
    int iteConstraints = 0, iteSzPool = 0;
    cutCover *coverCuts;
    int *intOrFloat = returnVectorTypeContraintsIntOrFloat(constraintsFull); //verifica se as restrições podem ser transformadas em restrições da mochila
    cutSmall *knapsackConstraints = reduceCutFullForCutSmall(constraintsFull, intOrFloat, precision);
    for (i = 0; i < knapsackConstraints->numberConstraints; i++)
    {
        if (knapsackConstraints->rightSide <= 1)
        {
            continue;
        }
        int szConstraint = knapsackConstraints->ElementsConstraints[i + 1] - knapsackConstraints->ElementsConstraints[i];
        cutCover *cutsCoverSolution = AllocStrCover(szConstraint * szPoolCutsMax, szPoolCutsMax);
        int nPoolUsed = 0;
        int *poolSolution = coverLCIAdamPerConstraints(knapsackConstraints, precision, i, nIterationGrasp, szPoolCutsMax, alpha, minimal, &nPoolUsed);
    }
}

int *coverLCIAdamPerConstraints(cutSmall *knapsackConstraints, int precision, int constraint, int nIterationGrasp, int szPoolCutsMax, double alpha, int minimal, int *nPoolUsed)
{
    int sz = knapsackConstraints->ElementsConstraints[constraint + 1] - knapsackConstraints->ElementsConstraints[constraint];
    int *solution = (int *)malloc(sizeof(int) * sz);
    int *poolSolution = (int *)malloc(sizeof(int) * sz * szPoolCutsMax);
    int nPoolSolution = 0;
    int ite = 0;
    double testAlpha = alpha;
    int i = 0;
    for (ite = 0; ite < sz * szPoolCutsMax; ite++)
    {
        poolSolution[ite] = 0;
    }
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
    if (minimal == 0)
    {
        while (ite <= nIterationGrasp)
        {
            memset(solution,0, sizeof(int)*sz);
            
            
            
            initialLCIAdamGRASP(solutionNew, szNew, knapsackConstraintsTemp, precision, constraint, alpha);
            for(i=knapsackConstraintsTemp->ElementsConstraints[constraint];i<knapsackConstraintsTemp->ElementsConstraints[constraint+1]; i++){
                if(solutionNew[ i - knapsackConstraintsTemp->ElementsConstraints[constraint] ] == 1){
                    int cAux = knapsackConstraintsTemp->originalConstraints[constr]
                    solution[ i - ]
                }
            }
            knapsackConstraints->originalConstraints
            int *solFinalTemp;
            cutCover *coverTemp = CopyCutToCover(constraintsSmall);
            solFinalTemp = LCIAdam(solution, coverTemp, constraint);
            if (verifyViolationGreedy(solFinalTemp, constraintsSmall, constraint, precision) == 1)
            {
                copyAndVerifyPoolSolution(solution, sz, poolSolution, &nPoolSolution);
                //printf("teste aqui");
                (*nPoolTemp)++;
            }

            free(solFinalTemp);
            solution = localSearch(solution, sz, constraintsSmall, precision, constraint, minimal);
            solFinalTemp = LCIAdam(solution, coverTemp, constraint);

            if (verifyViolationGreedy(solFinalTemp, constraintsSmall, constraint, precision) == 1)
            {
                copyAndVerifyPoolSolution(solution, sz, poolSolution, &nPoolSolution);
                // printf("teste aqui");
                (*nPoolTemp)++;
            }
            free(coverTemp);

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
        }
    }
    else
    {

        while (ite < nIterationGrasp)
        {
            initialLCIAdamGRASP(solution, sz, constraintsSmall, precision, constraint, alpha);
            if (verifySolutionCoverMinimal(solution, constraintsSmall, constraint) == 1)
            {
                int *solFinalTemp;
                cutCover *coverTemp = CopyCutToCover(constraintsSmall);
                solFinalTemp = LCIAdam(solution, coverTemp, constraint);
                if (verifyViolationGreedy(solFinalTemp, constraintsSmall, constraint, precision) == 1)
                {
                    copyAndVerifyPoolSolution(solution, sz, poolSolution, &nPoolSolution);
                    //     printf("teste aqui");
                    (*nPoolTemp)++;
                }
                free(solFinalTemp);
                free(coverTemp);
            }
            else
            {
                ite++;
                continue;
            }
            solution = localSearch(solution, sz, constraintsSmall, precision, constraint, minimal);
            if (verifySolutionCoverMinimal(solution, constraintsSmall, constraint) == 1)
            {
                int *solFinalTemp;
                cutCover *coverTemp = CopyCutToCover(constraintsSmall);

                solFinalTemp = LCIAdam(solution, coverTemp, constraint);

                if (verifyViolationGreedy(solFinalTemp, constraintsSmall, constraint, precision) == 1)
                {
                    copyAndVerifyPoolSolution(solution, sz, poolSolution, &nPoolSolution);
                    //printf("teste aqui");
                    (*nPoolTemp)++;
                }
                free(solFinalTemp);
                free(coverTemp);
            }
            if (nPoolSolution == szPoolCutsMax)
            {
                break;
            }

            if (testAlpha == -1)
            {
                alpha = fRand(0.0, 0.1);
                //printf("alpha: %f\n", alpha);
            }
            ite++;
        }
    }
    free(solution);
    return poolSolution;
}





void initialLCIAdamGRASP(int *solution, int sz, cutSmall *constraintsSmall, int precision, int constraint, float alpha)
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
        solution[i] = 0;
    }
    for (i = constraintsSmall->ElementsConstraints[constraint]; i < constraintsSmall->ElementsConstraints[constraint + 1]; i++)
    {
        setId[aux] = i;
        aux++;
    }
    szAux = aux;
    lhs = 0;
    int test = 0;
    int contSizeNewSolution = 0;
    while (verifyCompleteSolution == 0)
    {
        aux = 0;
        for (i = 0; i < szAux; i++)
        {
            el = constraintsSmall->Elements[setId[i]];
            value = (double)constraintsSmall->xAsterisc[el] / (double)constraintsSmall->Coefficients[setId[i]];
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
            el = constraintsSmall->Elements[setId[i]];
            value = (double)constraintsSmall->xAsterisc[el] / (double)constraintsSmall->Coefficients[setId[i]];
            if (value >= value_best)
            {
                setTemp[aux_t] = setId[i];
                posTemp[aux_t] = i;
                aux_t++;
            }
        }
        int itemAdd = rand() % aux_t;
        el = setTemp[itemAdd] - constraintsSmall->ElementsConstraints[constraint];
        solution[el] = 1;
        contSizeNewSolution++;
        szAux--;
        lhs += constraintsSmall->Coefficients[setTemp[itemAdd]];
        if (lhs - 1e-5 > constraintsSmall->rightSide[constraint])
        {
            //int *solTemp;
            //cutCover *coverTemp = CopyCutToCover(constraintsSmall);
            //if(typeLifted==0){
            //    solTemp = LCIAdam(solution, coverTemp,constraint);
            //}else{
            //    solTemp = LCIBallas(solution, coverTemp,constraint);
            //}
            //if(verifyViolationGreedy(solTemp,constraintsSmall,constraint,precision)==1){

            //  free(coverTemp);
            //free(solTemp);
            verifyCompleteSolution = 1;
            free(setTemp);
            free(posTemp);
            break;
            //}
            //free(coverTemp);
            //free(solTemp);
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


cutSmall *remove_X_0_and_1(cutSmall *knapsackConstraints, int precision)
{
    int i, j, el, cont = 0, contAux, nConstraintsTemp = knapsackConstraints->numberConstraints;
    int *testConstraints = (int *)malloc(sizeof(int) * nConstraintsTemp);
    for (i = 0; i < knapsackConstraints->numberConstraints; i++)
    {
        testConstraints[i] = 1;
        contAux = 0;
        for (j = knapsackConstraints->ElementsConstraints[i]; j < knapsackConstraints->ElementsConstraints[i + 1]; j++)
        {
            el = knapsackConstraints->Elements[j];
            if ((knapsackConstraints->xAsterisc[el] != 1 * precision) && (knapsackConstraints->xAsterisc[el] != 0))
            {
                contAux++;
            }
        }
        if (contAux < 2)
        {
            testConstraints[i] = 0;
            nConstraintsTemp--;
        }
        else
        {
            cont += contAux;
        }
    }
    cutSmall *knapTemp = AllocStrCutSmall(cont, nConstraintsTemp, knapsackConstraints->numberVariables);
    TRightSide rhsTemp;
    cont = 0;
    contAux = 0;
    knapTemp->ElementsConstraints[0] = 0;
    for (i = 0; i < knapsackConstraints->numberConstraints; i++)
    {
        if (testConstraints[i])
        {
            rhsTemp = knapsackConstraints->rightSide[i];
            for (j = knapsackConstraints->ElementsConstraints[i]; j < knapsackConstraints->ElementsConstraints[i + 1]; j++)
            {
                el = knapsackConstraints->Elements[j];
                if ((knapsackConstraints->xAsterisc[el] != 1 * precision) && (knapsackConstraints->xAsterisc[el] != 0))
                {
                    knapTemp->Coefficients[cont] = knapsackConstraints->Coefficients[j];
                    knapTemp->Elements[cont] = el;
                    cont++;
                }
                else if (knapsackConstraints->xAsterisc[el] == 1 * precision)
                {
                    rhsTemp -= knapsackConstraints->Coefficients[j];
                }
            }
            knapTemp->originalConstraints[contAux] = i;
            contAux++;
            knapTemp->ElementsConstraints[contAux] = cont;
        }
    }
    for (i = 0; i < knapsackConstraints->numberVariables; i++)
    {
        knapTemp->xAsterisc[i] = knapsackConstraints->xAsterisc[i];
    }
    return knapTemp;
}
