#ifndef OutputAsymMatrix_H
#define OutputAsymMatrix_H

#include "lduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class OutputAsymMatrix Declaration
\*---------------------------------------------------------------------------*/

class OutputAsymMatrix
:
    public lduMatrix::solver
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        OutputAsymMatrix(const OutputAsymMatrix&);

        //- Disallow default bitwise assignment
        void operator=(const OutputAsymMatrix&);


public:

    //- Runtime type information
    TypeName("OutputAsymMatrix");


    // Constructors

        //- Construct from matrix components and solver data stream
        OutputAsymMatrix
        (
            const word& fieldName,
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& interfaceBouCoeffs,
            const FieldField<Field, scalar>& interfaceIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaces,
            const dictionary& solverControls
        );


    // Destructor

        virtual ~OutputAsymMatrix()
        {}


    // Member Functions

        //- Solve the matrix with this solver
        virtual solverPerformance solve
        (
            scalarField& psi,
            const scalarField& source,
            const direction cmpt=0
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
