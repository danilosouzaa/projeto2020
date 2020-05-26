#include "newGrasp.h"

constraintsReal *graspLciAdam(constraintsReal *constraintsFull, int precision, char **nameConstraints, char **nameVariables, int szPoolCutsMax, int nIterationGrasp, float alpha, int minimal){
    int i, j;
    int iteConstraints = 0, iteSzPool = 0;
    cutCover *coverCuts;
    int *intOrFloat = returnVectorTypeContraintsIntOrFloat(constraintsFull);  //verifica se as restrições podem ser transformadas em restrições da mochila
    cutSmall *knapsackConstraints = reduceCutFullForCutSmall(constraintsFull, intOrFloat, precision);
    cutSmall *knapsackConstraintsTemp = remove_X_0_and_1(knapsackConstraints,precision);
    for(i=0;i< knapsackConstraintsTemp->numberConstraints;i++){
        if(knapsackConstraintsTemp->rightSide<=1){
            continue;
        }
        






    }


}

cutSmall *remove_X_0_and_1(cutSmall *knapsackConstraints, int precision){
    int i, j, el, cont = 0, contAux, nConstraintsTemp = knapsackConstraints->numberConstraints;
    int *testConstraints = (int*)malloc(sizeof(int)*nConstraintsTemp);
    for(i=0;i<knapsackConstraints->numberConstraints;i++){
        testConstraints[i] = 1;
        contAux = 0;
        for(j = knapsackConstraints->ElementsConstraints[i];j < knapsackConstraints->ElementsConstraints[i+1];j++){
            el = knapsackConstraints->Elements[j];
            if ((knapsackConstraints->xAsterisc[el] != 1*precision)&&(knapsackConstraints->xAsterisc[el]!=0)){
                contAux ++;
            }
        }
        if(contAux<2){
            testConstraints[i] = 0;
            nConstraintsTemp--;
        }else{
            cont += contAux;
        }

    }
    cutSmall *knapTemp = AllocStrCutSmall(cont,nConstraintsTemp,knapsackConstraints->numberVariables);
    TRightSide rhsTemp;
    cont = 0;
    contAux = 0;
    knapTemp->ElementsConstraints[0] = 0;
    for(i=0;i<knapsackConstraints->numberConstraints;i++){
        if(testConstraints[i]){
            rhsTemp = knapsackConstraints->rightSide[i];
            for(j = knapsackConstraints->ElementsConstraints[i];j < knapsackConstraints->ElementsConstraints[i+1];j++){
                el = knapsackConstraints->Elements[j];
                if ((knapsackConstraints->xAsterisc[el] != 1*precision)&&(knapsackConstraints->xAsterisc[el]!=0)){
                    knapTemp->Coefficients[cont] = knapsackConstraints->Coefficients[j];
                    knapTemp->Elements[cont]= el;
                    cont ++;
                }else if(knapsackConstraints->xAsterisc[el] == 1*precision){
                    rhsTemp -= knapsackConstraints->Coefficients[j];
                }
            }
            knapTemp->originalConstraints[contAux] = i;
            contAux++;
            knapTemp->ElementsConstraints[contAux] = cont;
        }
    }
    for(i=0;i<knapsackConstraints->numberVariables;i++){
        knapTemp->xAsterisc[i] = knapsackConstraints->xAsterisc[i];
    }
    return knapTemp;
}
