#include "SA.h"






constraintsReal *runCCwithSA(constraintsReal *constraintsFull, int precision, char **nameConstraints, char **nameVariables, int szPoolCutsMax, int nIterationSA, float TempInitial, float FactorResf){
    
    int *intOrFloat = returnVectorTypeContraintsIntOrFloat(constraintsFull); //verifica se as restrições podem ser transformadas em restrições da mochila
    cutSmall *knapsackConstraints = reduceCutFullForCutSmall(constraintsFull, intOrFloat, precision);
    cutCover *coverCuts = CopyCutToCover(knapsackConstraints);
    int i, j, k, itePool;

    for(i = 0; i< knapsackConstraints->numberConstraints;i++){
        int lhs = 0;
        if (coverCuts->rightSide[i]<=1){
            continue;
        }
        





    

    //Criar Solução Inicial



    //Loop do Resfriamento

    //Loop da condição de parada




    //Atulizar a solução


    }
    
    //retornar a solução






}

