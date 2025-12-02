/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "kplanckAbsorptionEmission.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(kplanckAbsorptionEmission, 0);

        addToRunTimeSelectionTable
        (
            absorptionEmissionModel,
            kplanckAbsorptionEmission,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::kplanckAbsorptionEmission::kplanckAbsorptionEmission
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    absorptionEmissionModel(dict, mesh),
    coeffsDict_(dict.optionalSubDict(typeName + "Coeffs")),
    kappaBg("backgroundKappa",coeffsDict_),
    a_
    (
        IOobject
        (
            "a",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless/dimLength, Zero)
    ),
    e_
    (
        IOobject
        (
            "e",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless/dimLength, Zero)
    ),
    E_
    (
        IOobject
        (
            "E",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimMass/dimLength/pow3(dimTime), Zero)
    )
{
    /*const volScalarField& T=this->mesh().objectRegistry::lookupObject<volScalarField>("T");
    const volScalarField& CO=this->mesh().objectRegistry::lookupObject<volScalarField>("CO");
    const volScalarField& CO2=this->mesh().objectRegistry::lookupObject<volScalarField>("CO2");
    const volScalarField& H2O=this->mesh().objectRegistry::lookupObject<volScalarField>("H2O");

    scalar kp_CO;
    scalar kp_CO2;
    scalar kp_H2O;
    forAll(mesh.cells(),celli)
    {
        kp_CO2=(-7.90115071277856e-29*pow(T[celli],9)
            +1.24707742855733e-24*pow(T[celli],8)
            -8.45678888740049e-21*pow(T[celli],7)
            +3.21699070979214e-17*pow(T[celli],6)
            -7.51143134718651e-14*pow(T[celli],5)
            +1.10486914624208e-10*pow(T[celli],4)
            -1.00741978549599e-07*pow(T[celli],3)
            +5.35226089320996e-05*pow(T[celli],2)
            -0.01463453045937304*T[celli]
            +1.82109714454140);
        kp_CO=(-5.41774420621837e-30*pow(T[celli],9)
            +8.27249165538230e-26*pow(T[celli],8)
            -5.37799559368833e-22*pow(T[celli],7)
            +1.93490005574305e-18*pow(T[celli],6)
            -4.18121036975696e-15*pow(T[celli],5)
            +5.47645428973631e-12*pow(T[celli],4)
            -4.10211417397809e-09*pow(T[celli],3)
            +1.43055573284301e-06*pow(T[celli],2)
            -4.58276209080535e-05*T[celli]
            -0.0337678342397966);
        kp_H2O=(-2.01015011815984e-29*pow(T[celli],9)
            +3.20110711055158e-25*pow(T[celli],8)
            -2.20265716090849e-21*pow(T[celli],7)
            +8.58129731742980e-18*pow(T[celli],6)
            -2.08423549337768e-14*pow(T[celli],5)
            +3.27608028430091e-11*pow(T[celli],4)
            -3.34940928455180e-08*pow(T[celli],3)
            +2.17801398577614e-05*pow(T[celli],2)
            -0.00846009525729881*T[celli]
            +1.62725753686949);
        
        a_[celli]=CO2[celli]*kp_CO2+CO[celli]*kp_CO+H2O[celli]*kp_H2O;
    }

    // Convert units from cm-1 to m-1
    a_=100.0*a_;
    e_=a_;
    */

    // Vertical Gradient in sq. cavity
    /*const scalarField& C = mesh.C().component(1);
    forAll(mesh.cells(),celli)
    {
        a_[celli]=C[celli]+0.5;
    }
    e_=a_;
    */
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiation::kplanckAbsorptionEmission::aCont(const label bandI) const
{
    tmp<volScalarField> ta
    (
        new volScalarField
        (
            IOobject
            (
                "a",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar(dimless/dimLength, Zero)
        )
    );

    scalarField& a = ta.ref().primitiveFieldRef();
    const volScalarField& T=this->mesh().objectRegistry::lookupObject<volScalarField>("T");
    const volScalarField& CO2=this->mesh().objectRegistry::lookupObject<volScalarField>("CO2");
    const volScalarField& H2O=this->mesh().objectRegistry::lookupObject<volScalarField>("H2O");

    scalar kp_CO2;
    scalar kp_H2O;
    forAll(this->mesh().cells(),celli)
    //forAll(T,celli)
    //forAll(a_, celli)
    {
        kp_CO2=(-7.90115071277856e-29*pow(T[celli],9)
            +1.24707742855733e-24*pow(T[celli],8)
            -8.45678888740049e-21*pow(T[celli],7)
            +3.21699070979214e-17*pow(T[celli],6)
            -7.51143134718651e-14*pow(T[celli],5)
            +1.10486914624208e-10*pow(T[celli],4)
            -1.00741978549599e-07*pow(T[celli],3)
            +5.35226089320996e-05*pow(T[celli],2)
            -0.01463453045937304*T[celli]
            +1.82109714454140);
        kp_H2O=(-2.01015011815984e-29*pow(T[celli],9)
            +3.20110711055158e-25*pow(T[celli],8)
            -2.20265716090849e-21*pow(T[celli],7)
            +8.58129731742980e-18*pow(T[celli],6)
            -2.08423549337768e-14*pow(T[celli],5)
            +3.27608028430091e-11*pow(T[celli],4)
            -3.34940928455180e-08*pow(T[celli],3)
            +2.17801398577614e-05*pow(T[celli],2)
            -0.00846009525729881*T[celli]
            +1.62725753686949);

        a[celli] = CO2[celli]*kp_CO2+H2O[celli]*kp_H2O;
        a[celli] = 100.0*a[celli]+kappaBg.value();
    }

    // Convert units from cm-1 to m-1
    //a_=100.0*a_;
    //a_=ta;
    //e_=a_;
    return ta;
    //return a_;
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::kplanckAbsorptionEmission::eCont(const label bandI) const
{
    /*tmp<volScalarField> te
    (
        new volScalarField
        (
            IOobject
            (
                "e",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            e_
        )
    );*/

    return aCont(bandI);
    //return e_;
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::kplanckAbsorptionEmission::ECont(const label bandI) const
{
    /*tmp<volScalarField> tE
    (
        new volScalarField
        (
            IOobject
            (
                "E",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            E_
        )
    );*/

    return E_;
}


// ************************************************************************* //
