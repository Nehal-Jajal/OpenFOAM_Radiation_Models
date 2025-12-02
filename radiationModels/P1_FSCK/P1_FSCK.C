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

#include "P1_FSCK.H"
#include "fvmLaplacian.H"
#include "fvmSup.H"
#include "absorptionEmissionModel.H"
#include "scatterModel.H"
#include "constants.H"
#include "addToRunTimeSelectionTable.H"
#include "fsckContainer.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(P1_FSCK, 0);
        addToRadiationRunTimeSelectionTables(P1_FSCK);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::P1_FSCK::P1_FSCK(const volScalarField& T)
:
    radiationModel(typeName, T),
    G_
    (
        IOobject
        (
            "G",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    qr_
    (
        IOobject
        (
            "qr",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/pow3(dimTime), Zero)
    ),
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
    e_
    (
        IOobject
        (
            "e",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimLength, Zero)
    ),
    E_
    (
        IOobject
        (
            "E",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/dimLength/pow3(dimTime), Zero)
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
    ),
    Ru_
    (
        IOobject
        (
            "RuCalc",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            false
        ),
        mesh_,
        dimensionedScalar(dimMass/dimLength/pow3(dimTime), Zero)
    )
{}


Foam::radiation::P1_FSCK::P1_FSCK(const dictionary& dict, const volScalarField& T)
:
    radiationModel(typeName, dict, T),
    G_
    (
        IOobject
        (
            "G",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    qr_
    (
        IOobject
        (
            "qr",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/pow3(dimTime), Zero)
    ),
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
    e_
    (
        IOobject
        (
            "e",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimLength, Zero)
    ),
    E_
    (
        IOobject
        (
            "E",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/dimLength/pow3(dimTime), Zero)
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
    ),
    Ru_
    (
        IOobject
        (
            "RuCalc",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            false
        ),
        mesh_,
        dimensionedScalar(dimMass/dimLength/pow3(dimTime), Zero)
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::radiation::P1_FSCK::read()
{
    if (radiationModel::read())
    {
        // nothing to read

        return true;
    }

    return false;
}


void Foam::radiation::P1_FSCK::calculate()
{
    //a_ = absorptionEmission_->a();
    //e_ = absorptionEmission_->e();
    E_ = absorptionEmission_->E();
    const volScalarField sigmaEff(scatter_->sigmaEff());

    const dimensionedScalar a0("a0", a_.dimensions(), ROOTVSMALL);

    // Construct diffusion
    volScalarField gamma
    (
        IOobject
        (
            "gammaRad",
            G_.mesh().time().timeName(),
            G_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        1.0/(3.0*a_ + sigmaEff + a0)
    );

    volScalarField aScale
    (
        IOobject
        (
            "aScaleFsck",
            G_.mesh().time().timeName(),
            G_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        G_.mesh(),
        dimensionedScalar(dimless, Zero)
    );

    Info<<"\n Running FSCK version of P1 \n"<<endl;
    //scalar kfsck[8] = {9.767899e-7,5.360105e-6,8.152979e-5,1.340214e-3,2.387425e-2,0.301935,0.986868,3.687812};
    //scalar wq[8] = {0.333888,0.261700,0.180514,0.120834,6.253695e-2,3.095662e-2,8.205659e-3,1.364635e-3};
    //scalar hundred = 1e2;

    std::vector<std::vector<double>> kcell;
    std::vector<std::vector<double>> acell;
    Foam::radiation::fsckContainer* fsck = Foam::radiation::fsckContainer::getInstance();
    //absorptionEmission_->fsckGetCellka(kcell,acell);
    fsck->fsckGetCellka(kcell,acell);
    std::vector<double> weight;
    //weight = absorptionEmission_->getWeight();
    weight = fsck->getWeight();

    // Convert units to m-1 from cm-1
    //kfsck = hundred*kfsck;

    // Calculate radiative heat flux on boundaries.
    volScalarField::Boundary& qrBf = qr_.boundaryFieldRef();
    qrBf = 0;

    delqr_ = dimensionedScalar(dimMass/dimLength/pow3(dimTime), Zero);
    Rp_ = dimensionedScalar(dimMass/dimLength/pow3(dimTime)/pow4(dimTemperature), Zero);
    Ru_ = delqr_;
    const label points = weight.size(); // Need to find out how many quadrature points are being used!!!
    for (label j=0; j < points; j++)
    {
        //a_ = dimensionedScalar(dimless/dimLength,1e2*kfsck[j]);
        //e_ = a_;
        forAll(G_.mesh().cells(),celli)
        {
            a_[celli]=kcell[celli][j];//dimensionedScalar(dimless/dimLength,kcell[celli][j]);
            e_[celli]=a_[celli]*acell[celli][j];
            aScale[celli]=acell[celli][j];
        }

        a_.correctBoundaryConditions();
        e_.correctBoundaryConditions();
        aScale.correctBoundaryConditions();
        
        gamma = 1.0/(3.0*a_ + sigmaEff + a0);

        //absorptionEmission_->bandSet(j);
        fsck->bandSet(j);

        G_ = 0.0*G_;
        Info<<"Solving for quadrature point "<<j<<endl;
        // Solve G transport equation
        solve
        (
            fvm::laplacian(gamma, G_)
          - fvm::Sp(a_, G_)
         ==
          - 4.0*(e_*physicoChemical::sigma*pow4(T_)) - E_
        );

        delqr_= delqr_ + weight[j]*a_*(4.0*aScale*physicoChemical::sigma*pow4(T_) - G_); // FIX THE ACELL HERE!

        // Compute the Rp and Ru terms here
        Rp_ = Rp_ + weight[j]*a_*aScale*4.0*physicoChemical::sigma;
        Ru_ = Ru_ + weight[j]*a_*G_;

        const volScalarField::Boundary& GBf = G_.boundaryField();
        const volScalarField::Boundary& gammaBf = gamma.boundaryField();    

        forAll(mesh_.boundaryMesh(), patchi)
        {
            // All patches are coupled - Work out this qr function properly
            //Info << "Patch " << patchi << endl;
            if (!GBf[patchi].coupled())
            {
                //Info << "This patch was coupled " << patchi << endl;
                qrBf[patchi] = qrBf[patchi] - gammaBf[patchi]*GBf[patchi].snGrad();
            }
        }
        // Sync the solution over quadrature points over all processors
        label sync = j;
        reduce(sync, sumOp<label>());
     }
}


Foam::tmp<Foam::volScalarField> Foam::radiation::P1_FSCK::Rp() const
{
    /*return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Rp",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            4.0*absorptionEmission_->eCont()*physicoChemical::sigma
        )
    );*/

    return Rp_;
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::radiation::P1_FSCK::Ru() const
{
    /*const volScalarField::Internal& G = G_();
    const volScalarField::Internal E = absorptionEmission_->ECont()()();
    const volScalarField::Internal a = absorptionEmission_->aCont()()();

    return a*G - E;*/
    const volScalarField::Internal& Ru_int = Ru_;
    return Ru_int;
}


Foam::label Foam::radiation::P1_FSCK::nBands() const
{
    return absorptionEmission_->nBands();
} 
// ************************************************************************* //
