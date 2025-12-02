/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "OpticallyThin.H"
#include "fvmLaplacian.H"
#include "fvmSup.H"
#include "absorptionEmissionModel.H"
#include "scatterModel.H"
#include "constants.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(OpticallyThin, 0);
        addToRadiationRunTimeSelectionTables(OpticallyThin);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::OpticallyThin::OpticallyThin(const volScalarField& T)
:
    radiationModel(typeName, T),
    a_
    (
        IOobject
        (
            "a",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimLength, Zero)
    ),
    delqr_
    (
        IOobject
        (
            "delqr",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/dimLength/pow3(dimTime), Zero)
    ),
    Rp_
    (
        IOobject
        (
            "RpCalc",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/dimLength/pow3(dimTime)/pow4(dimTemperature), Zero)
    )
{}


Foam::radiation::OpticallyThin::OpticallyThin(const dictionary& dict, const volScalarField& T)
:
    radiationModel(typeName, dict, T),
    a_
    (
        IOobject
        (
            "a",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimLength, Zero)
    ),
    delqr_
    (
        IOobject
        (
            "delqr",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/dimLength/pow3(dimTime), Zero)
    ),
    Rp_
    (
        IOobject
        (
            "RpCalc",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/dimLength/pow3(dimTime)/pow4(dimTemperature), Zero)
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::radiation::OpticallyThin::read()
{
    if (radiationModel::read())
    {
        return true;
    }

    return false;
}


void Foam::radiation::OpticallyThin::calculate()
{
    // Only works for kplanckFsckAbsorptionEmission model
    a_ = absorptionEmission_->a();

    // Evaluate the divergence of radiative heat flux
    Rp_ = a_*4.0*physicoChemical::sigma;
    delqr_= Rp_*pow4(T_); //a_*(4.0*physicoChemical::sigma*pow4(T_));
}


Foam::tmp<Foam::volScalarField> Foam::radiation::OpticallyThin::Rp() const
{
    return Rp_;
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::radiation::OpticallyThin::Ru() const
{
    const volScalarField::Internal E = absorptionEmission_->ECont()()();
    // Will return zero for kplanckFsckAbsorptionEmission model (other models are invalid)
    // Emission is not contributing to radiation source term in optically thin model
    return E;
}


Foam::label Foam::radiation::OpticallyThin::nBands() const
{
    return absorptionEmission_->nBands();
}


// ************************************************************************* //
