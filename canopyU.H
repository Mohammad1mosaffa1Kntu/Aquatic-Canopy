/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2019 OpenFOAM Foundation
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

Class
    Foam::fv::canopyU

Description
    Calculates and applies the buoyancy force rho*g to the momentum equation
    corresponding to the specified velocity field.

Usage
    Example usage:
    \verbatim
    fields          (U);                    // Name of velocity field
    \endverbatim

SourceFiles
    canopyU.C

\*---------------------------------------------------------------------------*/

#ifndef canopyU_H
#define canopyU_H

#include "fvModel.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
               Class canopyU Declaration
\*---------------------------------------------------------------------------*/

class canopyU
:
    public fvModel
{
    // Private Data

        //- Name of density field
        const word rhoName_;
        
         //- Optional phase name
        word phaseName_;
        
        //- Name of velocity field
        word UName_;

        // Fields
       //- Plant canopy drag coefficient field [-]
       volScalarField plantCd_;

       //- Leaf area density field [1/m]
       volScalarField leafAreaDensity_;
       
       //- Non-virtual read
        void readCoeffs();

public:

    //- Runtime type information
    TypeName("canopyU");


    // Constructors

        //- Construct from explicit source name and mesh
        canopyU
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Disallow default bitwise copy construction
        canopyU(const canopyU&) = delete;


    // Member Functions
    
    	    //- Return the list of fields for which the fvModel adds source term
            //  to the transport equation
            virtual wordList addSupFields() const;	

        // Evaluate

            //- Add explicit contribution to incompressible momentum equation
            virtual void addSup
            (
                fvMatrix<vector>& eqn,
                const word& fieldName
            ) const;

            //- Add explicit contribution to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const word& fieldName
            ) const;

	    //- Add explicit contribution to phase momentum equation
            virtual void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const word& fieldName
            ) const;	
        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const canopyU&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
