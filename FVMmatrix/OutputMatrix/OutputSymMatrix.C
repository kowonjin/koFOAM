/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/
#include <cstdio>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include "OutputSymMatrix.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(OutputSymMatrix, 0);

    lduMatrix::solver::addsymMatrixConstructorToTable<OutputSymMatrix>
        addPCGSymMatrixConstructorToTable_;
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
        lduMatrix::preconditioner::getName(controlDict_) + typeName,
        fieldName_
    );

    int n = psi.size();
    //scalar* a = (scalar *)malloc(sizeof(scalar)*n*n);
    scalar* a = new scalar[n*n];
    

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


    std::ofstream coefFile("P.txt");
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
