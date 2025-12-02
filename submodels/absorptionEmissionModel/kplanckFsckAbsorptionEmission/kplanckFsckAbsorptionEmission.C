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

#include "kplanckFsckAbsorptionEmission.H"
#include "addToRunTimeSelectionTable.H"
#include "fskTableFunction.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(kplanckFsckAbsorptionEmission, 0);

        addToRunTimeSelectionTable
        (
            absorptionEmissionModel,
            kplanckFsckAbsorptionEmission,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::kplanckFsckAbsorptionEmission::kplanckFsckAbsorptionEmission
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
    ),
    calcTref(coeffsDict_.getOrDefault<bool>("calcTref", false)),
    isV2(coeffsDict_.getOrDefault<bool>("isV2", false)),
    mixingModel(coeffsDict_.getOrDefault<bool>("mixingModel", false)),
    mix_scheme(coeffsDict_.getOrDefault<label>("mixingScheme", 1)),
    Tref_input("Tref", coeffsDict_)
{
    Tref = Tref_input.value(); // Units - K
    // Load the FSK Table Data into memory (variable kq)
    std::string fpath;
    word dictPath = word(coeffsDict_.lookup("DataPath"));
    if (dictPath == "OSC")
    {
        if (isV2)
        {
            fpath = "/users/PAA0008/nehaljajal/NongrayData/fsck/";
        }
        else
        {
            fpath = "/users/PAA0008/nehaljajal/fskNewDB/DataBank/";
        }
    }
    else
    {
        fpath = "/home/cfd-user/fskNewDB/DataBank/";
    }
    //std::string fpath = "/home/cfd-user/fskNewDB/DataBank/";
    if (mixingModel)
    {
        /*loadIntoMem(kq, fpath);
        fpath = "/users/PAA0008/nehaljajal/NongrayDataFSCK/fuelFsck/";
        loadIntoMemV2(kq2, fpath, Tref);*/
        
        fpath = "/users/PAA0008/nehaljajal/NongrayDataFSCK/productFsck/";
        loadIntoMemProduct(kq, fpath, Tref);
        fpath = "/users/PAA0008/nehaljajal/NongrayDataFSCK/fuelFsck/";
        loadIntoMemFuel(kq2, fpath, Tref);
        Info << "FSCK Tables with fuel and products loaded into memory." << endl;
    }
    else
    {
        if (isV2)
        {
            loadIntoMemV2(kq, fpath, Tref);
        }
        else
        {
            loadIntoMem(kq, fpath);
        }
    }

    // Generate the quadrature points
    const bool cheb2 = false;
    const int points = 32;
    const float alpha = 2.0;
    g.assign(points,0.0);
    weight.assign(points,0.0);
    quadgen2(cheb2,g,weight,points,alpha);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiation::kplanckFsckAbsorptionEmission::aCont(const label bandI) const
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
    
    // Get the volScalar fields for computing Planck Mean from FSCK table
    const volScalarField& T=this->mesh().objectRegistry::lookupObject<volScalarField>("T");
    const volScalarField& CO2=this->mesh().objectRegistry::lookupObject<volScalarField>("X_CO2");
    const volScalarField& H2O=this->mesh().objectRegistry::lookupObject<volScalarField>("X_H2O");
    const volScalarField& CO=this->mesh().objectRegistry::lookupObject<volScalarField>("X_CO");
    const volScalarField& NH3=this->mesh().objectRegistry::lookupObject<volScalarField>("X_NH3");
    const volScalarField& CH4=this->mesh().objectRegistry::lookupObject<volScalarField>("X_CH4");

    // Compute Reference temperature (Volume average of the domain)
    // Reference Temperature and vector with to load all data into memory
    double Tref_compute;
    if (calcTref)
    {
        
        Tref_compute = 0.0;
        scalar totalVolume = 0.0;
        forAll(T,celli)
        {
            Tref_compute = Tref_compute + T[celli]*mesh_.V()[celli];
            totalVolume = totalVolume + mesh_.V()[celli];
        }
        reduce(Tref_compute, sumOp<double>());
        reduce(totalVolume, sumOp<scalar>());
        Tref_compute = Tref_compute/totalVolume;
    }
    else
    {
        Tref_compute = Tref;
    }

    Info << "Reference temperature is set to - " << Tref_compute << endl;
    Info << "Background kappa value is - " << kappaBg.value() << endl;
    
    // Compute Planck Mean abs coeff values for all the cells
    forAll(this->mesh().cells(),celli)
    {
        
        GasMixInfo Mix_Info;
        Mix_Info.T = T[celli];
        Mix_Info.xCO = CO[celli];
        Mix_Info.xCO2 = CO2[celli];
        Mix_Info.xH2O = H2O[celli];
        Mix_Info.xNH3 = NH3[celli];
        Mix_Info.xCH4 = CH4[celli];
        
        std::vector<double> k;
        std::vector<double> a;
        if (mixingModel)
        {
            get_ka_mix(Mix_Info,Tref,kq,kq2,g,weight,k,a,mix_scheme);
        }
        else
        {
            get_ka(Mix_Info,Tref,kq,g,weight,isV2,k,a);
        }

        int points = k.size();
        for (int i = 1; i <= points; ++i)
        {
            a[celli] = a[celli] + weight[i-1]*a[i-1]*k[i-1];
        }
    }

    // Convert units from cm-1 to m-1
    a=100.0*a+kappaBg.value();

    return ta;
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::kplanckFsckAbsorptionEmission::eCont(const label bandI) const
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
Foam::radiation::kplanckFsckAbsorptionEmission::ECont(const label bandI) const
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
