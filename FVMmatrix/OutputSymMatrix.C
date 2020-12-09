#include <cstdio>
#include <stdio.h>
#include <stdlib.h>

#include "OutputSymMatrix.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(OutputSymMatrix, 0);

    lduMatrix::solver::addsymMatrixConstructorToTable<OutputSymMatrix>
        addOutputSymMatrixSymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::OutputSymMatrix::OutputSymMatrix
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const FieldField<Field, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const dictionary& solverControls
)
:
    lduMatrix::solver
    (
        fieldName,
        matrix,
        interfaceBouCoeffs,
        interfaceIntCoeffs,
        interfaces,
        solverControls
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::solverPerformance Foam::OutputSymMatrix::solve
(
    scalarField& psi,
    const scalarField& source,
    const direction cmpt
) const
{
    // --- Setup class containing solver performance data
    solverPerformance solverPerf
    (
        typeName,
        fieldName_
    );

    int n = psi.size();
    Info << "  *********** SymMatrix **********  " <<endl;
    Info << " sizeof(scalar) = " << sizeof(scalar) << endl;
    scalar* a = (scalar *)malloc(sizeof(scalar)*n*n);
    //scalar *a = new scalar(n^3);
    

    /*for(int i = 0; i < n; i++){
	    for(int j = 0; j < n; j++){
		    a[i*n + j] = 0.;
	    }
    }*/

    printf("upper\n");
    for(int i = 0; i < matrix_.lower().size(); i++){
	    printf("%d %d\n", i, matrix_.lduAddr().upperAddr()[i]);
    }
    printf("lower\n");
    for(int i = 0; i < matrix_.upper().size(); i++){
	    printf("%d %d\n", i, matrix_.lduAddr().lowerAddr()[i]);
    }

    for(int i = 0; i < matrix_.lower().size(); i++){
	    int r = matrix_.lduAddr().upperAddr()[i];
	    int c = matrix_.lduAddr().lowerAddr()[i];
	    a[r*n + c] = matrix_.lower()[i];
    }

    for(int i = 0; i < matrix_.diag().size(); i++){
	    a[i*n + i] = matrix_.diag()[i];
    }

    for(int i = 0; i < matrix_.upper().size(); i++){
	    int r = matrix_.lduAddr().lowerAddr()[i];
	    int c = matrix_.lduAddr().upperAddr()[i];
	    a[r*n + c] = matrix_.upper()[i];
    }

    char filename[32];
    sprintf(filename, "matrix-%s.txt", fieldName_.c_str());

    FILE *fp = fopen(filename, "w");
    if(fp == NULL){
        printf("error: Can't open file: OutputSymMatrix.txt\n");
	std::exit(0);
    }

    fprintf(fp, "%d\n", n);

    for(int i = 0; i < n; i++){
	    for(int j = 0; j < n; j++){
		    fprintf(fp, "%15.7e", a[i*n + j]);
	    }
	    fprintf(fp, "\n");
    }

    for(int i = 0; i < n; i++) fprintf(fp, "%15.7e\n", source[i]);

    fclose(fp);

    free(a);
    
    return solverPerf;
}


// ************************************************************************* //
