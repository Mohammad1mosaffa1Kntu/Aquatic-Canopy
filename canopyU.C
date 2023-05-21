/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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

#include "canopyU.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "fvmSup.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(canopyU, 0);

    addToRunTimeSelectionTable
    (
        fvModel,
        canopyU,
        dictionary
    );
}
}
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::canopyU::readCoeffs()
{
    phaseName_ = coeffs().lookupOrDefault<word>("phase", word::null);

    UName_ =
        coeffs().lookupOrDefault<word>
        (
            "U",
            IOobject::groupName("U", phaseName_)
        );
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::canopyU::canopyU
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fvModel(name, modelType, dict, mesh),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    phaseName_(word::null),
    UName_(word::null),
    plantCd_
    (
        IOobject
        (
            "plantCd",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    leafAreaDensity_
    (
        IOobject
        (
            "leafAreaDensity",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    )

{
        readCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::wordList Foam::fv::canopyU::addSupFields() const
{
    return wordList(1, UName_);
}

void Foam::fv::canopyU::addSup 
(
    fvMatrix<vector>& eqn,
    const word& fieldName
)const
{
    //Info << "add1" << endl;
    const volVectorField& U = eqn.psi();// access to the field
    //const volVectorField& U = mesh_.lookupObject<volVectorField>("U");
    //eqn -= (plantCd_*leafAreaDensity_*mag(U))*U;
    eqn += fvm::Sp(plantCd_*leafAreaDensity_*mag(U),U);//TODO: check + or -
}


void Foam::fv::canopyU::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const word& fieldName
)const
{
    //Info << "add2" << endl; 	
    const volVectorField& U = eqn.psi();
    //const volVectorField& U = mesh_.lookupObject<volVectorField>("U");	
    //eqn -= rho*(plantCd_*leafAreaDensity_*mag(U))*U;
    eqn -= fvm::Sp(rho*plantCd_*leafAreaDensity_*mag(U),U);//TODO: check + or -
}

void Foam::fv::canopyU::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const word& fieldName
) const
{
    const volVectorField& U = eqn.psi();
    //const volVectorField& U = mesh_.lookupObject<volVectorField>("U");	
    //eqn -= rho*(plantCd_*leafAreaDensity_*mag(U))*U;
    eqn -= fvm::Sp(alpha*rho*plantCd_*leafAreaDensity_*mag(U),U);//TODO: check + or -
    //eqn += alpha*rho*g_;
}

bool Foam::fv::canopyU::read(const dictionary& dict)
{
    if (fvModel::read(dict))
    {
        readCoeffs();
        return true;
    }
    else
    {
        return false;
    }
}




// ************************************************************************* //
