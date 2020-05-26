#define INT 1

#ifdef INT
    typedef int TCoefficients;
    typedef int TElements;
    typedef int TElementsConstraints;
    typedef int TRightSide;
    typedef int TXAsterisc;
    typedef int TNumberVariables;
    typedef int TNumberConstraints;
    typedef int TCont;
    typedef short TInterval;
    typedef int TList;
    typedef int TPosList;

#endif 

typedef float TParametersGroup;
typedef double TCoefficientsFull;
typedef double TRightSideFull;
typedef double TXAsteriscFull;
typedef char TActivedCut;

typedef struct {
    TNumberVariables numberVariables;
    TNumberConstraints numberConstraints;
    TCont cont;
    TCoefficients *Coefficients;
    TElements *Elements;
    TElementsConstraints *ElementsConstraints;
    TRightSide *rightSide;
    TXAsterisc *xAsterisc;
    TNumberConstraints *originalConstraints;
}cutSmall;

typedef struct{
    TNumberVariables numberVariables;
    TNumberConstraints numberConstraints;
    TCont cont;
    TCoefficientsFull *Coefficients;
    TElements *Elements;
    TElementsConstraints *ElementsConstraints;
    TRightSideFull *rightSide;
    TXAsteriscFull *xAsterisc;
}constraintsReal;

typedef struct{
    TNumberConstraints numberCuts;
    TActivedCut activedCut;
    constraintsReal *pCuts;
}cutPool;

typedef struct{
    TNumberConstraints numberConstraints;
    TCont cont;
    TCoefficients *Coefficients;
    TElementsConstraints *ElementsConstraints;
    TRightSide *rightSide;
}cutCover;

