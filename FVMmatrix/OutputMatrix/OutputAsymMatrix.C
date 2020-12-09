#include <cstdio>
#include <stdlib.h>
#include <fstream>
#include <iostream>

#include "OutputAsymMatrix.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(OutputAsymMatrix, 0);

    lduMatrix::solver::addasymMatrixConstructorToTable<OutputAsymMatrix>
        addOutputAsymMatrixAsymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::OutputAsymMatrix::OutputAsymMatrix
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

Foam::solverPerformance Foam::OutputAsymMatrix::solve
(
    scalarField& psi,
    const scalarField& source,
    const direction cmpt
) const
{
    // --- Setup class containing solver performance data
    solverPerformance solverPerf
    (
        lduMatrix::preconditioner::getName(controlDict_) + typeName,
        fieldName_
    );

    int n = psi.size();
    //scalar* a = (scalar *)malloc(sizeof(scalar)*n*n);
    scalar* a = new scalar[n*n];
    scalar* b = new scalar[n][n];
 
    
    for(int i = 0; i < n; i++){
	    for(int j = 0; j < n; j++){
		    a[i*n + j] = 0.;
	    }
    }

    forAll(matrix_.lower(), i){
	    int r = matrix_.lduAddr().upperAddr()[i];
	    int c = matrix_.lduAddr().lowerAddr()[i];
	    a[r*n + c] = matrix_.lower()[i];
    }

    forAll(matrix_.diag(), i){
	    a[i*n + i] = matrix_.diag()[i];
    }

    forAll(matrix_.upper(), i){
	    int r = matrix_.lduAddr().lowerAddr()[i];
	    int c = matrix_.lduAddr().upperAddr()[i];
	    a[r*n + c] = matrix_.upper()[i];
    }

    char filename[32];
    sprintf(filename, "matrix-%s.txt", fieldName_.c_str());
    std::ofstream coefFile(filename);
    for(int i=0; i < n; i++)
    {
        for(int j=0; j<n; j++) 
        {
           coefFile << a[i*n+j] << "    ";
        } coefFile<< " | " << source[i] << std::endl;
        
    } 
    coefFile.close();

    delete(a);


    return solverPerf;
}


// ************************************************************************* //
