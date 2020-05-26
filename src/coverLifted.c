#include "coverLifted.h"

int *LCIAdam(int *coverSolution, cutCover *constraintsCover, TNumberConstraints constraint)
{
    double fillBag = 0;
    int i, szCover = 0;
    int szConstraints = constraintsCover->ElementsConstraints[constraint + 1] - constraintsCover->ElementsConstraints[constraint];
    for (i = constraintsCover->ElementsConstraints[constraint]; i < constraintsCover->ElementsConstraints[constraint + 1]; i++)
    {
        fillBag += coverSolution[i - constraintsCover->ElementsConstraints[constraint]] * constraintsCover->Coefficients[i];
        szCover += coverSolution[i - constraintsCover->ElementsConstraints[constraint]];
    }
    int *n_coef = (int *)malloc(sizeof(int) * szCover);
    int *n_el = (int *)malloc(sizeof(int) * szCover);
    int caux = 0, qnt = 0;
    for (i = constraintsCover->ElementsConstraints[constraint]; i < constraintsCover->ElementsConstraints[constraint + 1]; i++)
    {
        if (coverSolution[caux] == 1)
        {
            n_coef[qnt] = constraintsCover->Coefficients[i];
            n_el[qnt] = i;
            qnt++;
        }
        caux++;
    }
    double b = (double)constraintsCover->rightSide[constraint];
    double delta = 0;
    double phi = ((double)fillBag - (double)constraintsCover->rightSide[constraint]);
    int k = 1, w;
    double a_barra = (double)n_coef[0];
    for (w = 1; w < qnt; w++)
    {

        delta = a_barra - (double)n_coef[w];
        if (((double)k * delta) < phi)
        {
            a_barra = (double)n_coef[w];
            phi = phi - ((double)k * delta);
        }
        else
        {
            a_barra = a_barra - (phi / (double)k);
            phi = 0;
            break;
        }
        k++;
    }
    if (phi > 0)
    {
        a_barra = (double)b / (double)szCover;
    }
    //printf("a_barra: %f\n", a_barra);
    int *c_menus = (int *)malloc(sizeof(int) * szCover);
    int *c_mais = (int *)malloc(sizeof(int) * szCover);
    double *S_barra = (double *)malloc(sizeof(double) * (szCover + 1));
    double *aux = (double *)malloc(sizeof(double) * szCover);
    int id1 = 0, id2 = 0;

    for (w = 0; w < qnt; w++)
    {
        if ((double)n_coef[w] <= a_barra)
        {
            c_menus[id1] = w;
            id1++;
        }
        else
        {
            c_mais[id2] = w;
            id2++;
        }
    }
    for (w = 0; w < qnt; w++)
    {
        if (n_coef[w] < a_barra)
        {
            aux[w] = n_coef[w];
        }
        else
        {
            aux[w] = a_barra;
        }
    }

    quicksortDouble(aux, 0, qnt);
    S_barra[0] = 0;
    for (w = 1; w < qnt; w++)
    {
        S_barra[w] = S_barra[w - 1] + aux[w - 1];
    }
    S_barra[qnt] = constraintsCover->rightSide[constraint];

    int ini = 0;
    int *cutCoverLifted = (int *)malloc(sizeof(int) * (szConstraints + 1));
    for (w = constraintsCover->ElementsConstraints[constraint]; w < constraintsCover->ElementsConstraints[constraint + 1]; w++)
    {
        int flag = 0;
        //cutsCoverSolution->Coefficients[c_XSolution] = cutsCover->Coefficients[w];
        for (ini = 0; ini < qnt; ini++)
        {
            if ((constraintsCover->Coefficients[w] > S_barra[ini]) && (constraintsCover->Coefficients[w] <= S_barra[ini + 1]))
            {
                cutCoverLifted[w - constraintsCover->ElementsConstraints[constraint]] = ini;
                flag = 1;
                break;
            }
        }
        if (flag == 0)
        {
            cutCoverLifted[w - constraintsCover->ElementsConstraints[constraint]] = constraintsCover->Coefficients[w];
        }
    }
    int el;
    for (w = 0; w < id1; w++)
    {
        el = n_el[c_menus[w]];
        cutCoverLifted[el - constraintsCover->ElementsConstraints[constraint]] = 1;
    }
    cutCoverLifted[szConstraints] = szCover - 1;
    free(c_menus);
    free(c_mais);
    free(S_barra);
    free(aux);
    free(n_coef);
    free(n_el);
    return cutCoverLifted;
}

int *LCIBallas(int *coverSolution, cutCover *constraintsCover, TNumberConstraints constraint)
{
    int i, w, szCover = 0;
    int szConstraints = constraintsCover->ElementsConstraints[constraint + 1] - constraintsCover->ElementsConstraints[constraint];
    for (i = constraintsCover->ElementsConstraints[constraint]; i < constraintsCover->ElementsConstraints[constraint + 1]; i++)
    {
        szCover += coverSolution[i - constraintsCover->ElementsConstraints[constraint]];
    }
    double *S_barra = (double *)malloc(sizeof(double) * (szCover + 1));
    S_barra[0] = 0;
    w = 1;
    for (i = constraintsCover->ElementsConstraints[constraint]; i < constraintsCover->ElementsConstraints[constraint + 1]; i++)
    {
        if (coverSolution[i - constraintsCover->ElementsConstraints[constraint]] == 1)
        {
            S_barra[w] = S_barra[w - 1] + constraintsCover->Coefficients[i];
            w++;
        }
    }

    int ini = 0;
    int *cutCoverLifted = (int *)malloc(sizeof(int) * (szConstraints + 1));
    for (w = constraintsCover->ElementsConstraints[constraint]; w < constraintsCover->ElementsConstraints[constraint + 1]; w++)
    {

        //cutsCoverSolution->Coefficients[c_XSolution] = cutsCover->Coefficients[w];
        if (coverSolution[w - constraintsCover->ElementsConstraints[constraint]] == 1)
        {
            cutCoverLifted[w - constraintsCover->ElementsConstraints[constraint]] = 1;
        }
        else
        {
            int flag = 0;
            for (ini = 0; ini < szCover; ini++)
            {
                if ((constraintsCover->Coefficients[w] > S_barra[ini]) && (constraintsCover->Coefficients[w] <= S_barra[ini + 1]))
                {
                    cutCoverLifted[w - constraintsCover->ElementsConstraints[constraint]] = ini;
                    flag = 1;
                    break;
                }
            }
            if (flag == 0)
            {
                cutCoverLifted[w - constraintsCover->ElementsConstraints[constraint]] = constraintsCover->Coefficients[w];
            }
        }
    }

    cutCoverLifted[szConstraints] = szCover - 1;
    free(S_barra);
    return cutCoverLifted;
}

cutSmall *reduceCutFullForCutSmall(constraintsReal *constraints, int *typeIntOrFloat, int precision)
{
    TNumberConstraints numberConstraintsSmall = 0;
    TCont cont = 0;
    int i, j, cont_aux;
    fflush(stdin);
    for (i = 0; i < constraints->numberConstraints; i++)
    {
        if (typeIntOrFloat[i])
        {
            numberConstraintsSmall++;
            cont += (constraints->ElementsConstraints[i + 1] - constraints->ElementsConstraints[i]);
        }
    }

    cutSmall *cSmall = AllocStrCutSmall(cont, numberConstraintsSmall, constraints->numberVariables);
    cont = 0;
    cont_aux = 0;
    cSmall->ElementsConstraints[0] = 0;
    for (i = 0; i < constraints->numberConstraints; i++)
    {
        if (typeIntOrFloat[i])
        {
            for (j = constraints->ElementsConstraints[i]; j < constraints->ElementsConstraints[i + 1]; j++)
            {
                cSmall->Coefficients[cont] = (TCoefficients)constraints->Coefficients[j];
                cSmall->Elements[cont] = constraints->Elements[j];
                cont++;
            }
            cSmall->rightSide[cont_aux] = (TRightSide)constraints->rightSide[i];
            cSmall->ElementsConstraints[cont_aux + 1] = cont;
            cSmall->originalConstraints[cont_aux] = i;
            cont_aux++;
        }
    }
    fflush(stdout);
    //printf("teste: %d %d \n", cont, numberConstraintsSmall);
    for (i = 0; i < constraints->numberVariables; i++)
    {
        cSmall->xAsterisc[i] = (int)(constraints->xAsterisc[i] * (double)precision);
    }
    return cSmall;
}

constraintsReal *runCCwithGrasp(constraintsReal *constraintsFull, int precision, char **nameConstraints, char **nameVariables, int szPoolCutsMax, int nIterationGrasp, float alpha, int minimal, int typeLifted)
{
    int i, k, itePool = 0;
    //int qnt,j;
    int c_XSolution = 0;
    int c_AuxSolution = 0;
    cutCover *cutsCover;
    int *intOrFloat = returnVectorTypeContraintsIntOrFloat(constraintsFull); //verifica se as restrições podem ser transformadas em restrições da mochila
    cutSmall *newConstraintsSmall = reduceCutFullForCutSmall(constraintsFull, intOrFloat, precision);
    cutsCover = CopyCutToCover(newConstraintsSmall);
    for (i = 0; i < cutsCover->numberConstraints; i++)
    {
        int lhs = 0;
        // for (j = newConstraintsSmall->ElementsConstraints[i]; j < newConstraintsSmall->ElementsConstraints[i + 1]; j++)
        // {
        //     lhs += newConstraintsSmall->Coefficients[j];
        // }
        if ((cutsCover->rightSide[i] <= 1))
        {
            continue;
        }
        c_AuxSolution = 0;
        c_XSolution = 0;
        int szConstraint = cutsCover->ElementsConstraints[i + 1] - cutsCover->ElementsConstraints[i];
        cutCover *cutsCoverSolution = AllocStrCover(szConstraint * szPoolCutsMax, szPoolCutsMax);
        int nPoolUsed = 0;
        int *poolSolution = createCoverGraspIndividual(newConstraintsSmall, precision, i, nIterationGrasp, szPoolCutsMax, alpha, minimal, typeLifted, &nPoolUsed);
        //printf("pooooolllll::%d\n", nPoolUsed);
        cutsCoverSolution->ElementsConstraints[0] = 0;
        for (itePool = 0; itePool < nPoolUsed; itePool++)
        {
            //qnt = 0;
            int caux = 0;
            lhs = 0;
            int *solutionTemp = (int *)malloc(sizeof(int) * szConstraint);
            for (k = cutsCover->ElementsConstraints[i]; k < cutsCover->ElementsConstraints[i + 1]; k++)
            {
                //  qnt += poolSolution[caux + itePool * szConstraint];

                solutionTemp[k - cutsCover->ElementsConstraints[i]] = poolSolution[(k - cutsCover->ElementsConstraints[i]) + itePool * szConstraint];
                lhs += cutsCover->Coefficients[k] * poolSolution[caux + itePool * szConstraint];
                caux++;
            }
            if (lhs <= cutsCover->rightSide[i])
            {
                free(solutionTemp);
                continue;
            }

            // for (k = 0; k < szConstraint; k++)
            // {
            //     solutionTemp[k] = poolSolution[k + itePool * szConstraint];
            // }
            if (typeLifted == 0)
            {
                // printf("LCI Adam!\n");
                int *cutCoverLifted = LCIAdam(solutionTemp, cutsCover, i);
                for (k = 0; k < szConstraint; k++)
                {
                    cutsCoverSolution->Coefficients[c_XSolution] = cutCoverLifted[k];
                    c_XSolution++;
                }
                cutsCoverSolution->ElementsConstraints[c_AuxSolution + 1] = c_XSolution;
                cutsCoverSolution->rightSide[c_AuxSolution] = cutCoverLifted[szConstraint];
                c_AuxSolution++;
                free(cutCoverLifted);
            }
            else
            {
                // printf("LCI Ballas!\n");
                minimal = 1;
                int *cutCoverLifted = LCIBallas(solutionTemp, cutsCover, i);
                for (k = 0; k < szConstraint; k++)
                {
                    cutsCoverSolution->Coefficients[c_XSolution] = cutCoverLifted[k];
                    c_XSolution++;
                }
                cutsCoverSolution->ElementsConstraints[c_AuxSolution + 1] = c_XSolution;
                cutsCoverSolution->rightSide[c_AuxSolution] = cutCoverLifted[szConstraint];
                c_AuxSolution++;
                free(cutCoverLifted);
            }

            free(solutionTemp);
        }
        //int qnt_cuts_cover = c_AuxSolution;
        int *idc_cover = (int *)malloc(sizeof(int) * c_AuxSolution);
        //memset(idc_cover,1, sizeof(int) * c_AuxSolution);
        // for (j = 0; j < c_AuxSolution; j++)
        // {
        //     idc_cover[j] = 1;
        //     qnt_cuts_cover++;
        // }
        constraintsFull = createCutsCoverGrasp(cutsCoverSolution, constraintsFull, newConstraintsSmall, idc_cover, i, c_AuxSolution, precision);
        free(cutsCoverSolution);
        free(idc_cover);
        free(poolSolution);
    }
    free(cutsCover);
    free(newConstraintsSmall);
    free(intOrFloat);
    return constraintsFull;
}

int *createCoverGraspIndividual(cutSmall *constraintsSmall, int precision, TNumberConstraints constraint, int numberIteration, int szPoolCutsMax, float alpha, int minimal, int typeLifted, int *nPoolTemp)
{
    int sz = constraintsSmall->ElementsConstraints[constraint + 1] - constraintsSmall->ElementsConstraints[constraint];
    int *solution = (int *)malloc(sizeof(int) * sz);
    int *poolSolution = (int *)malloc(sizeof(int) * sz * szPoolCutsMax);
    int nPoolSolution = 0;
    int ite = 0;
    double testAlpha = alpha;
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
        while (ite <= numberIteration)
        {

            createInitialCoverGRASP(solution, sz, constraintsSmall, precision, constraint, alpha, typeLifted);
            int *solFinalTemp;
            cutCover *coverTemp = CopyCutToCover(constraintsSmall);
            if (typeLifted == 0)
            {
                solFinalTemp = LCIAdam(solution, coverTemp, constraint);
            }
            else
            {
                solFinalTemp = LCIBallas(solution, coverTemp, constraint);
            }
            if (verifyViolationGreedy(solFinalTemp, constraintsSmall, constraint, precision) == 1)
            {

                copyAndVerifyPoolSolution(solution, sz, poolSolution, &nPoolSolution);
                //printf("teste aqui");
                (*nPoolTemp)++;
            }

            free(solFinalTemp);
            solution = localSearch(solution, sz, constraintsSmall, precision, constraint, minimal);
            if (typeLifted == 0)
            {
                solFinalTemp = LCIAdam(solution, coverTemp, constraint);
            }
            else
            {
                solFinalTemp = LCIBallas(solution, coverTemp, constraint);
            }
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

        while (ite < numberIteration)
        {
            createInitialCoverGRASP(solution, sz, constraintsSmall, precision, constraint, alpha, typeLifted);
            if (verifySolutionCoverMinimal(solution, constraintsSmall, constraint) == 1)
            {
                int *solFinalTemp;
                cutCover *coverTemp = CopyCutToCover(constraintsSmall);
                if (typeLifted == 0)
                {
                    solFinalTemp = LCIAdam(solution, coverTemp, constraint);
                }
                else
                {
                    solFinalTemp = LCIBallas(solution, coverTemp, constraint);
                }
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
                if (typeLifted == 0)
                {
                    solFinalTemp = LCIAdam(solution, coverTemp, constraint);
                }
                else
                {
                    solFinalTemp = LCIBallas(solution, coverTemp, constraint);
                }
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

void createInitialCoverGRASP(int *solution, int sz, cutSmall *constraintsSmall, int precision, int constraint, float alpha, int typeLifted)
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

void copyAndVerifyPoolSolution(int *solution, int sz, int *poolSolution, int *numberSolutionAtual)
{
    int i, j, atual;
    int flag = *numberSolutionAtual;
    atual = flag;
    for (i = 0; i < *numberSolutionAtual; i++)
    {
        for (j = 0; j < sz; j++)
        {
            if (solution[j] != poolSolution[j + i * sz])
            {
                flag--;
                break;
            }
        }
    }
    if (flag == 0)
    {
        for (j = 0; j < sz; j++)
        {
            poolSolution[j + atual * sz] = solution[j];
        }
        (*numberSolutionAtual)++;
    }
}

constraintsReal *createCutsCoverGrasp(cutCover *cutsCover, constraintsReal *constraintsOriginal, cutSmall *constraintsSmall, int *idc_Cover, int constraint, int nCuts, int precision)
{
    if (nCuts <= 0)
    {
        return constraintsOriginal;
    }
    long int i = 0, j = 0, cont = 0, contConstraints = 0;
    for (i = 0; i < nCuts; i++)
    {
        // if (idc_Cover[i] == 0)
        // {
        //     continue;
        // }
        // double violation = valueViolation(cutsCover, constraintsSmall, i, constraint, precision);
        //printf("%f\n",violation);
        //if (violation == 0)
        //{
        //    idc_Cover[i] = 0;
        //}
        //if (idc_Cover[i] == 1)
        //{
        contConstraints++;
        for (j = cutsCover->ElementsConstraints[i]; j < cutsCover->ElementsConstraints[i + 1]; j++)
        {
            if (cutsCover->Coefficients[j] != 0)
            {
                cont++;
            }
        }
        //}
    }
    if (contConstraints == 0)
    {
        return constraintsOriginal;
    }

    constraintsReal *outCutsNew = AllocStrConstraintsReal(constraintsOriginal->cont + cont, constraintsOriginal->numberConstraints + contConstraints, constraintsOriginal->numberVariables);
    outCutsNew->ElementsConstraints[0] = constraintsOriginal->ElementsConstraints[0];

    for (i = 0; i < constraintsOriginal->numberConstraints; i++)
    {
        outCutsNew->rightSide[i] = constraintsOriginal->rightSide[i];
        outCutsNew->ElementsConstraints[i + 1] = constraintsOriginal->ElementsConstraints[i + 1];
    }
    for (i = 0; i < constraintsOriginal->cont; i++)
    {
        outCutsNew->Elements[i] = constraintsOriginal->Elements[i];
        outCutsNew->Coefficients[i] = constraintsOriginal->Coefficients[i];
    }
    for (i = 0; i < constraintsOriginal->numberVariables; i++)
    {
        outCutsNew->xAsterisc[i] = constraintsOriginal->xAsterisc[i];
    }
    int aux = constraintsOriginal->numberConstraints;
    int c_aux = constraintsOriginal->cont;
    for (i = 0; i < nCuts; i++)
    {
        //if (idc_Cover[i] == 1)
        //{
        int c_XSolution = constraintsSmall->ElementsConstraints[constraint];
        outCutsNew->rightSide[aux] = cutsCover->rightSide[i];
        for (j = cutsCover->ElementsConstraints[i]; j < cutsCover->ElementsConstraints[i + 1]; j++)
        {
            if (cutsCover->Coefficients[j] != 0)
            {
                outCutsNew->Elements[c_aux] = constraintsSmall->Elements[c_XSolution];
                outCutsNew->Coefficients[c_aux] = cutsCover->Coefficients[j];
                c_aux++;
            }
            c_XSolution++;
        }
        outCutsNew->ElementsConstraints[aux + 1] = c_aux;
        aux++;
        //}
    }
    int **matrizIncidencia = (int **)malloc(sizeof(int *) * outCutsNew->numberConstraints);
    for (i = 0; i < outCutsNew->numberConstraints; i++)
    {
        matrizIncidencia[i] = (int *)malloc(sizeof(int) * outCutsNew->numberVariables);
        for (j = 0; j < outCutsNew->numberVariables; j++)
        {
            matrizIncidencia[i][j] = -1;
        }
        for (j = outCutsNew->ElementsConstraints[i]; j < outCutsNew->ElementsConstraints[i + 1]; j++)
        {       
                aux = outCutsNew->Elements[j];
                matrizIncidencia[i][aux] = outCutsNew->ElementsConstraints[i];
        }
    }

    int *validated = (int *)malloc(sizeof(int) * outCutsNew->numberConstraints);
    validated[0] = 1;
    for (i = 1; i < outCutsNew->numberConstraints; i++)
    {
        validated[i] = verifyRepeatedIncidency(matrizIncidencia, outCutsNew, i);
    }
    cont = 0, contConstraints = 0;
    for (i = 0; i < outCutsNew->numberConstraints; i++)
    {
        if (validated[i] == 1)
        {
            contConstraints++;
            for (j = outCutsNew->ElementsConstraints[i]; j < outCutsNew->ElementsConstraints[i + 1]; j++)
            {
                cont++;
            }
        }
    }
    if (contConstraints == outCutsNew->numberConstraints)
    {
        freeStrConstraintsReal(constraintsOriginal);
        free(validated);
        return outCutsNew;
    }
    constraintsReal *cutsNewNoRepetead = AllocStrConstraintsReal(cont, contConstraints, outCutsNew->numberVariables);
    cont = 0, contConstraints = 0;
    cutsNewNoRepetead->ElementsConstraints[0] = 0;
    for (i = 0; i < outCutsNew->numberConstraints; i++)
    {
        if (validated[i] == 1)
        {
            cutsNewNoRepetead->rightSide[contConstraints] = outCutsNew->rightSide[i];
            for (j = outCutsNew->ElementsConstraints[i]; j < outCutsNew->ElementsConstraints[i + 1]; j++)
            {
                cutsNewNoRepetead->Elements[cont] = outCutsNew->Elements[j];
                cutsNewNoRepetead->Coefficients[cont] = outCutsNew->Coefficients[j];
                cont++;
            }
            cutsNewNoRepetead->ElementsConstraints[contConstraints + 1] = cont;
            contConstraints++;
        }
    }
    for (i = 0; i < cutsNewNoRepetead->numberVariables; i++)
    {
        cutsNewNoRepetead->xAsterisc[i] = outCutsNew->xAsterisc[i];
    }
    for (i = 0; i < outCutsNew->numberConstraints; i++)
    {
        free(matrizIncidencia[i]);
    }
    free(matrizIncidencia);

    free(validated);
    freeStrConstraintsReal(constraintsOriginal);
    freeStrConstraintsReal(outCutsNew);
    return cutsNewNoRepetead;
}

int verifySolutionCoverMinimal(int *solution, cutSmall *constraintsSmall, TNumberConstraints constraint)
{
    int i = 0, lhs = 0, menor;
    menor = constraintsSmall->Coefficients[constraintsSmall->ElementsConstraints[constraint]];
    for (i = constraintsSmall->ElementsConstraints[constraint]; i < constraintsSmall->ElementsConstraints[constraint + 1]; i++)
    {
        if (solution[i - constraintsSmall->ElementsConstraints[constraint]] == 1)
        {
            lhs += constraintsSmall->Coefficients[i];
            if (constraintsSmall->Coefficients[i] < menor)
            {
                menor = constraintsSmall->Coefficients[i];
            }
        }
    }
    if (lhs - menor <= constraintsSmall->rightSide[constraint])
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

int *localSearch(int *solution, int sz, cutSmall *constraintsSmall, int precision, TNumberConstraints constraint, int minimal)
{
    int *newSolution = (int *)malloc(sizeof(int) * sz);
    int *vNActived = (int *)malloc(sizeof(int) * sz);
    int *vActived = (int *)malloc(sizeof(int) * sz);
    int i, szActive = 0, szNActive = 0;
    for (i = 0; i < sz; i++)
    {
        newSolution[i] = solution[i];
        if (solution[i] == 1)
        {
            vActived[szActive] = i;
            szActive++;
        }
        else
        {
            vNActived[szNActive] = i;
            szNActive++;
        }
    }
    shuffleVectorInt(vActived, szActive);
    newSolution[vActived[0]] = 0;
    while (szNActive > 0)
    {
        if (szNActive > 1)
        {
            shuffleVectorInt(vNActived, szNActive);
        }
        int k_t = vNActived[szNActive - 1];
        newSolution[k_t] = 1;
        szNActive--;
        if (minimal == 0)
        {
            if (verifySolutionCover(newSolution, constraintsSmall, precision, constraint) == 1)
            {
                free(solution);
                free(vNActived);
                free(vActived);
                return newSolution;
            }
        }
        else
        {
            if (verifySolutionCover(newSolution, constraintsSmall, precision, constraint) == 1)
            {
                if (verifySolutionCoverMinimal(newSolution, constraintsSmall, constraint) == 1)
                {
                    free(solution);
                    free(vNActived);
                    free(vActived);
                    return newSolution;
                }
            }
        }
    }
    free(newSolution);
    free(vActived);
    free(vNActived);
    return solution;
}

int verifySolutionCover(int *solution, cutSmall *constraintsSmall, int precision, TNumberConstraints constraint)
{
    double lhs = 0;
    int i, aux = 0;
    for (i = constraintsSmall->ElementsConstraints[constraint]; i < constraintsSmall->ElementsConstraints[constraint + 1]; i++)
    {
        //el = constraintsSmall->Elements[i];
        lhs += constraintsSmall->Coefficients[i] * solution[aux];
        aux++;
    }
    if ((lhs - 1e-5) > constraintsSmall->rightSide[constraint])
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

void shuffleVectorInt(int *vec, int sz)
{
    struct timeval time;
    gettimeofday(&time, NULL);
    srand((time.tv_sec * 1000) + (time.tv_usec / 1000));
    int j;
    float *num_temp = (float *)malloc(sizeof(float) * sz); //free
    for (j = 0; j < sz; j++)
    {
        num_temp[j] = (float)(rand() % RAND_MAX);
    }
    quicksortTParameters(num_temp, vec, 0, sz);
    free(num_temp);
}

int *greedyMethodForCC(cutSmall *constraintsUsed, TNumberConstraints constraint, int precision, int typeLift)
{
    cutCover *constraintsCover = CopyCutToCover(constraintsUsed);
    int sz = constraintsUsed->ElementsConstraints[constraint + 1] - constraintsUsed->ElementsConstraints[constraint];
    int *solutionCover = (int *)malloc(sizeof(int) * sz);
    double *xTemp = (double *)malloc(sizeof(double) * sz);
    int *idc = (int *)malloc(sizeof(int) * sz);
    int i, el, aux = 0;
    for (i = constraintsUsed->ElementsConstraints[constraint]; i < constraintsUsed->ElementsConstraints[constraint + 1]; i++)
    {
        el = constraintsUsed->Elements[i];
        xTemp[aux] = (double)constraintsUsed->xAsterisc[el] / (double)precision;
        idc[aux] = i;
        solutionCover[aux] = 0;
        aux++;
    }
    //getchar();
    quicksortCof(xTemp, idc, 0, sz);
    int lhs = 0, posRef = 0;
    for (i = 0; i < sz; i++)
    {
        el = idc[i];
        lhs += constraintsUsed->Coefficients[el];
        if (lhs > constraintsUsed->rightSide[constraint])
        {
            lhs -= constraintsUsed->Coefficients[el];
            posRef = i;
            break;
        }
        else
        {
            solutionCover[el - constraintsUsed->ElementsConstraints[constraint]] = 1;
        }
    }
    int k = 0;
    for (k = posRef; k < sz; k++)
    {
        int *solutionFinal;
        el = idc[k];
        solutionCover[el - constraintsUsed->ElementsConstraints[constraint]] = 1;
        lhs += constraintsUsed->Coefficients[el];
        if (typeLift == 0)
        {
            solutionFinal = LCIAdam(solutionCover, constraintsCover, constraint);
        }
        else
        {
            solutionFinal = LCIBallas(solutionCover, constraintsCover, constraint);
        }
        if (verifyViolationGreedy(solutionFinal, constraintsUsed, constraint, precision) == 1)
        {
            free(xTemp);
            free(idc);
            free(solutionCover);
            free(constraintsCover);
            return solutionFinal;
        }
        free(solutionFinal);
    }
    free(xTemp);
    free(idc);
    free(solutionCover);
    free(constraintsCover);
    return NULL;
}

int verifyViolationGreedy(int *solutionFinal, cutSmall *constraitsUsed, int constraint, int precision)
{
    double lhs = 0.0;
    int i, el;
    int sz = constraitsUsed->ElementsConstraints[constraint + 1] - constraitsUsed->ElementsConstraints[constraint];
    for (i = constraitsUsed->ElementsConstraints[constraint]; i < constraitsUsed->ElementsConstraints[constraint + 1]; i++)
    {
        el = constraitsUsed->Elements[i];
        int aux = (double)solutionFinal[i - constraitsUsed->ElementsConstraints[constraint]];
        aux *= (double)constraitsUsed->xAsterisc[el];
        lhs += aux;
    }
    lhs = lhs / (double)precision;
    if (lhs + 1e-5 > solutionFinal[sz])
    {
        //printf("%f %d", lhs, solutionFinal[sz]);
        return 1;
    }
    return 0;
}

constraintsReal *runCCGreedy(constraintsReal *constraintsFull, int precision, char **nameConstraints, char **nameVariables, int typeLifted)
{
    int i, j, k;
    cutCover *cutsCover;
    int *intOrFloat = returnVectorTypeContraintsIntOrFloat(constraintsFull);
    cutSmall *newConstraintsSmall = reduceCutFullForCutSmall(constraintsFull, intOrFloat, precision);
    cutsCover = CopyCutToCover(newConstraintsSmall);

    for (i = 0; i < cutsCover->numberConstraints; i++)
    {
        int lhs = 0;
        int szConstraint = newConstraintsSmall->ElementsConstraints[i + 1] - newConstraintsSmall->ElementsConstraints[i];
        cutCover *cutsCoverSolution = AllocStrCover(szConstraint + 1, 1);
        cutsCoverSolution->ElementsConstraints[0] = 0;
        int c_XSolution = 0;
        int c_AuxSolution = 0;

        for (j = newConstraintsSmall->ElementsConstraints[i]; j < newConstraintsSmall->ElementsConstraints[i + 1]; j++)
        {
            lhs += newConstraintsSmall->Coefficients[j];
        }
        if ((cutsCover->rightSide[i] <= 1) || (lhs <= cutsCover->rightSide[i]))
        {
            continue;
        }
        int *coverResult = greedyMethodForCC(newConstraintsSmall, i, precision, typeLifted);
        if (coverResult != NULL)
        {
            for (k = 0; k < szConstraint; k++)
            {
                //printf("Passou aqui!: %d\n", coverResult[k]);
                cutsCoverSolution->Coefficients[c_XSolution] = coverResult[k];
                c_XSolution++;
            }
            cutsCoverSolution->ElementsConstraints[c_AuxSolution + 1] = c_XSolution;
            cutsCoverSolution->rightSide[c_AuxSolution] = coverResult[szConstraint];
            c_AuxSolution++;
            free(coverResult);
            int *idc_cover = (int *)malloc(sizeof(int) * c_AuxSolution);
            idc_cover[0] = 1;
            constraintsFull = createCutsCoverGrasp(cutsCoverSolution, constraintsFull, newConstraintsSmall, idc_cover, i, c_AuxSolution, precision);
            free(idc_cover);
        }
        free(cutsCoverSolution);
    }
    //int qnt_cuts_cover = 0;

    //constraintsFull = createCutsCoverGrasp(cutsCoverSolution, constraintsFull, newConstraintsSmall, idc_cover, i, c_AuxSolution, precision);
    //free(cutsCoverSolution);

    free(cutsCover);
    free(newConstraintsSmall);
    free(intOrFloat);
    return constraintsFull;
}
