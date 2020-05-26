#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <omp.h>


extern "C"
{
    
    //#include "preProcessing.h"
    #include "coverLifted.h"
    //#include "lp.h"
}



int main(int argc, const char *argv[]){
    char nameFileInstance[255] = "../input/";
    char nameInst[255] = "";
    strcat(nameInst, argv[1]);
    strcat(nameFileInstance, argv[1]);
    double timeMax = atof(argv[2]);
    int precision = atoi(argv[3]);
    int typeCC = 0; // 0 - GRASP, 1 -Guloso
    int typeLift = 0;//  0 - Adam Lethford 1 - Bala1s
    int minimal = 0; // 0 - no minimal 1- minimal
    int szPoolCutsMaxCC = 0;
    int nIterationCCGrasp = 0;
    float alpha = 0.0;
    int i;
    int cg1, cg2, ck, cc ; 
    cg1 = 0;
    cg2 = 0;
    ck = 0;
    ck = 0;
    for (i = 0; i < argc; i++)
    {
        if (strcmp(argv[i], "-CG1") == 0)
        {
            cg1 = 1;
        }
        if (strcmp(argv[i], "-CG2") == 0)
        {
            cg2 = 1;
        }
        if (strcmp(argv[i], "-CC") == 0)
        {
            cc = 1;
            typeCC = atoi(argv[i+1]); //0 - grasp or 1 - greedy
            typeLift = atoi(argv[i+2]); //0 - Adam e  1 - ballas or adam
            minimal  = atoi(argv[i+3]); // 1 - cover minimal or 0 - not minimal
            szPoolCutsMaxCC = atoi(argv[i+4]); //size pool iterGrasp
            nIterationCCGrasp = atof(argv[i+5]); //
            alpha = atof(argv[i+6]);
        }
        if (strcmp(argv[i], "-CK") == 0)
        {
            ck = 1;
        }
    }
    int numberAux = 0, numberCutsCC = 0;
    LinearProgram *lp = lp_create();
    lp_read(lp, nameFileInstance);
    lp_set_print_messages(lp,0);
    lp_write_lp(lp,"daniloSe.lp");
    char **nameVariables = createNameVariablesInitial(lp);// struct and name
    char **nameConstraints = createStructNameConstraintsInitial(lp);// struct
    // TNumberConstraints nConstraints = countConstraintsValided(lp);
    int nVariables = lp_cols(lp);
    #ifdef DEBUG
        double *sol = initialSolutionForValidation(lp);
    #endif // DEBUG
    int *typeVariables = (int *)malloc(sizeof(int) * nVariables);
    double *lbVariables = (double *)malloc(sizeof(double) * nVariables);
    double *ubVariables = (double *)malloc(sizeof(double) * nVariables);
    constraintsReal *constraintsFull = fillStructPerLP(lp,nameConstraints, nameVariables,typeVariables,lbVariables, ubVariables);
    //showStructFull(constraintsFull, nameConstraints, nameVariables);
    lp_write_lp(lp,"../output/temp.lp");
  
    lp_optimize_as_continuous(lp);
    //getchar();
    double iniObjSol = lp_obj_value(lp);
    printf("Solution initial: %f \n", iniObjSol);
    lp_set_cuts(lp,'0');

    //------------------------------------------------------------------------------------
    // Initial Time Counting
    //------------------------------------------------------------------------------------
    double startT = omp_get_wtime();
    double _time = 0;
    _time = ((double)timeMax - (omp_get_wtime() - startT));
    int contNoImprovement = 0;
    int contSize = 1;
    if(typeCC==0){
        printf("Grasp\n");
    }else{
        printf("Greedy\n");
    }
    do{
    //------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------

        int numberAuxConstraints = constraintsFull->numberConstraints;

        
        //execute Separation 
        if(cc==1){  
            int *binaryConstraints = returnBinaryConstraints(constraintsFull, typeVariables); // verify constraints can be used fo method
            TNumberConstraints nConstraintsUsed = countConstraintsBinaryUsed(binaryConstraints,constraintsFull->numberConstraints); 
            constraintsReal *constraintsBinary = convertBinaryConstraints(constraintsFull, binaryConstraints, typeVariables, lbVariables, ubVariables);
            int nInitialBinary = constraintsBinary->numberConstraints;
            int numberVariablesInitial = constraintsBinary->numberVariables;
            int *convertVariables = (int *)malloc(sizeof(int) * constraintsBinary->cont);
            constraintsBinary = removeNegativeCoefficientsAndSort(constraintsBinary, convertVariables);
            //printf("ncUsed:%d\t", nConstraintsUsed);
            numberAux = constraintsBinary->numberConstraints;
            if(typeCC==0){
                
                constraintsBinary = runCCwithGrasp(constraintsBinary,precision,nameConstraints,nameVariables,szPoolCutsMaxCC, nIterationCCGrasp, alpha,minimal,typeLift);
            }else{
                
                constraintsBinary = runCCGreedy(constraintsBinary,precision,nameConstraints,nameVariables,typeLift);
            }

            numberAux = constraintsBinary->numberConstraints - numberAux;
            constraintsBinary = returnVariablesOriginals(constraintsBinary, convertVariables, precision, numberVariablesInitial);
            constraintsFull = convertBinaryOfOriginalConstraints(constraintsFull, constraintsBinary, nInitialBinary);
            if(numberAux != 0){
                nameConstraints = renamedNameConstraints(nameConstraints, 3, constraintsFull->numberConstraints, numberAux, numberCutsCC);
            }
            numberCutsCC += numberAux;
            freeStrConstraintsReal(constraintsBinary);
            free(binaryConstraints);
            free(convertVariables);
        }
        #ifdef DEBUG //validation of feasible of constraints
            int *verifyTest = (int *)malloc(sizeof(int) * constraintsFull->numberConstraints);
            for (i = 0; i < constraintsFull->numberConstraints; i++)
            {
                verifyTest[i] = verifyCutsValidatedPerSolutionInteger(constraintsFull, i, sol, nameVariables);
                if (verifyTest[i] == 0)
                {
                    printf("validado depois: %d %s\n", verifyTest[i], nameConstraints[i]);
                }
            }
        #endif
        int totalCuts = constraintsFull->numberConstraints - numberAuxConstraints;
        // printf("Depois: %d \n", constraintsOriginal->numberConstraints);
        if (totalCuts > 0)
        {
            printf("Round : %d LCI: %d\t",contSize, totalCuts);
            #ifdef DEBUG
                insertConstraintsLPDebug(lp, constraintsFull, numberAuxConstraints, nameConstraints, verifyTest);
            #else
                insertConstraintsLP(lp, constraintsFull, numberAuxConstraints, nameConstraints);
            #endif // DEBUG
            
            lp_write_lp(lp, "../output/temp.lp");
            lp_optimize_as_continuous(lp);
            double *xTemp = lp_x(lp);
            for (i = 0; i < constraintsFull->numberVariables; i++)
            {   //printf("%f\t", constraintsFull->xAsterisc[i]);
                constraintsFull->xAsterisc[i] = xTemp[i];
                //printf("%f\n", constraintsFull->xAsterisc[i]);
            }   
            if( (typeCC==1) && (iniObjSol==lp_obj_value(lp)) ){
                timeMax = 0;    
            }
            if( iniObjSol==lp_obj_value(lp) ){
                contNoImprovement ++;
            }else{
                _time = ((double)timeMax - (omp_get_wtime() - startT));
                printf("tImp: %f\t", (omp_get_wtime() - startT));
                contNoImprovement = 0;
            }
        
            //getchar();
            iniObjSol =  lp_obj_value(lp);
            printf("obj: %f\n", iniObjSol);
            
        }else{
            //printf("\n");
            contNoImprovement++;
        }
        if(contNoImprovement>1000){
            timeMax = 0;
        }
        #ifdef DEBUG
            free(verifyTest);
        #endif
    //------------------------------------------------------------------------------------
    // Final Time Counting
    //------------------------------------------------------------------------------------
        _time = ((double)timeMax - (omp_get_wtime() - startT));
    //}while(_time > 1);
        contSize++;
    }while((contSize <= 50)&&(_time>0));
    printf("ncF: %d tF: %f \n", numberCutsCC, (omp_get_wtime() - startT));
    //------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------
    #ifdef DEBUG    
        free(sol);
    #endif
    free(typeVariables);
    free(lbVariables);
    free(ubVariables);
    freeStructName(nameVariables, lp_cols(lp));
    freeStructName(nameConstraints,constraintsFull->numberConstraints);
    freeStrConstraintsReal(constraintsFull);
    lp_free(&lp);
    



}