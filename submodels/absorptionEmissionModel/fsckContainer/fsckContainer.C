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

#include "fsckContainer.H"
#include "addToRunTimeSelectionTable.H"
#include "rhoThermo.H"
#include "fskTableFunction.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(fsckContainer, 0);
    }
}

Foam::radiation::fsckContainer* Foam::radiation::fsckContainer::thisContainer = nullptr;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::fsckContainer::fsckContainer
(
    const dictionary& dict,
    const fvMesh& mesh,
    double Tref_in
)
:
    mesh_(mesh),
    coeffsDict_(dict.optionalSubDict(typeName + "Coeffs")),
    kappaBg("backgroundKappa",coeffsDict_),
    quadPoints(coeffsDict_.getOrDefault<label>("quadPoints", 32)),
    calcTref(coeffsDict_.getOrDefault<bool>("calcTref", false)),
    isV2(coeffsDict_.getOrDefault<bool>("isV2", false)),
    mixingModel(coeffsDict_.getOrDefault<bool>("mixingModel", false)),
    mix_scheme(coeffsDict_.getOrDefault<label>("mixingScheme", 1))
{
    // If mixing model is true, set isV2 to false
    if (mixingModel)
    {
        isV2 = false;
    }
    // Set default to first quadrature point
    // Set to zero for computation of k values across the domain
    bandID = 0;

    // Set reference Temperature for the fsck Table
    // CODE THIS PROPERLY TO ALTER BASE ON DOMAIN OR FIX
    Tref = Tref_in;

    // Load the FSK Table Data into memory (variable kq)
    std::string fpath;
    word dictPath = word(coeffsDict_.lookup("DataPath"));
    if (dictPath == "OSC")
    {
        if (isV2)
        {
            fpath = "/users/PAA0008/nehaljajal/NongrayDataFSCK/productFsck/";
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
    // Write out the time taken to load the data
    Info << "Time taken to load the data - " << mesh_.time().elapsedCpuTime() << endl;
    // Reset the time to zero
    //Info << "Resetting the time to zero " << endl;
    //resetCpuTime();

    // Generate the quadrature points
    const bool cheb2 = false;
    const int points = 32;
    const float alpha = 2.0;
    gBase.assign(points,0.0);
    weightBase.assign(points,0.0);
    quadgen2(cheb2,gBase,weightBase,points,alpha);
    if (quadPoints != points)
    {
        g.assign(quadPoints,0.0);
        weight.assign(quadPoints,0.0);
        quadgen2(cheb2,g,weight,quadPoints,alpha);
    }
    else
    {
        g = gBase;
        weight = weightBase;
    }
    
    /*
    // Write out the quadrature points and weights
    Info << "Base quadrature points and weights - " << endl;
    for (int i = 1; i <= points; ++i)
    {
        Info << i << "  " << gBase[i-1] << "  " << weightBase[i-1] << endl;
    }
    // Write out the quadrature points and weights
    Info << "Quadrature points and weights - " << endl;
    for (int i = 1; i <= quadPoints; ++i)
    {
        Info << i << "  " << g[i-1] << "  " << weight[i-1] << endl;
    }*/

    //TEST THE USAGE OF THE KDIST FUNCTIONS FROM THE CONSTRUCTOR
    /*GasMixInfo Mix_Info;
    Mix_Info.T = 1850.0;
    Mix_Info.P = 1.0;
    Mix_Info.fv = 0.0;
    Mix_Info.xCO = 0.0;
    Mix_Info.xCO2 = 0.05;
    Mix_Info.xH2O = 0.0;
    double Tref = 1000.0;
    Info << "Gas Mix info - " << endl << "T - " << Mix_Info.T << "\nP - " << Mix_Info.P
         << "\nfv - " << Mix_Info.fv << "\nxCO - " << Mix_Info.xCO << "\nxCO2 - " << Mix_Info.xCO2
         << "\nxH2O - " << Mix_Info.xH2O << "\nTref - " << Tref << endl;
    //std::vector<double> k;
    //std::vector<double> a;
    get_ka(Mix_Info,Tref,kq,gBase,weightBase,isV2,kBase,aBase);
    Info << endl << "k and a values for above inputs " << endl;
    for (int i = 1; i <= points; ++i)
    {
        Info << i << "  " << gBase[i-1] << "  " << kBase[i-1] << "  " << aBase[i-1] << endl;
    }*/
    

    // Implement code for fewer quadrature points
    /*k.assign(quadPoints,0.0);
    a.assign(quadPoints,0.0);
    lin_interp(gBase,kBase,g,k);
    lin_interp(gBase,aBase,g,a);

    Info << endl << "k and a values for fewer quadrature inputs " << endl;
    for (int i = 1; i <= quadPoints; ++i)
    {
        Info << i << "  " << g[i-1] << "  " << k[i-1] << "  " << a[i-1] << endl;
    }*/

    /*label bFaceSize = 0;
    forAll(mesh.boundaryMesh(),patchi)
    {
        const volScalarField& G = mesh.objectRegistry::lookupObject<volScalarField>("T");
        
        const scalarField& boundaryG = G.boundaryField()[patchi];
        
        Info<<"The boundary condition for this patch is "<<G.boundaryField()[patchi].type()<<endl;

        bFaceSize = bFaceSize + mesh.boundaryMesh()[patchi].size();
        forAll(mesh.boundaryMesh()[patchi],facei)
        {
            Info<<"The boundary face index is "
            <<mesh.boundaryMesh()[patchi].start()+facei-mesh.nInternalFaces()<<endl;
        }
    }*/

    //Info<<"Number of total faces "<<mesh.faces().size()<<endl;
    //Info<<"Number of internal faces "<<mesh.nInternalFaces()<<endl;
    //Info<<"Number of boundary faces "<<mesh.faces().size()-mesh.nInternalFaces()<<endl;
    //Info<<"Calculated boundary faces "<<bFaceSize<<endl;

    /*forAll(mesh.faces(),facei)
    {
        if(!mesh.isInternalFace(facei))
        {
            if(mesh.isInternalFace(facei-1))
            {
                Info<<"Boundary patch started at "<<facei<<endl;
            }
        }
    }*/

    // Initialize the 2D vectors for holding all quadrature point data
    kcellInternal.assign(mesh.cells().size(),std::vector<double>(quadPoints,0.0));
    acellInternal.assign(mesh.cells().size(),std::vector<double>(quadPoints,0.0));
    //kfaceInternal.assign(mesh_.boundaryMesh().size(),std::vector<double>(32,0.0));
    const label bFaceSize = mesh.faces().size()-mesh.nInternalFaces();
    afaceInternal.assign(bFaceSize,std::vector<double>(quadPoints,0.0));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::radiation::fsckContainer* Foam::radiation::fsckContainer::getInstance()
{
    return thisContainer;
}

Foam::radiation::fsckContainer* Foam::radiation::fsckContainer::getInstance(const dictionary& dict, const fvMesh& mesh, double Tref_in)
{
    if (thisContainer == nullptr)
    {
        thisContainer = new fsckContainer(dict, mesh, Tref_in);
    }
    return thisContainer;
}

void Foam::radiation::fsckContainer::bandSet(label bandI)
{
    bandID = bandI;
    Info<<"Band set to "<<bandID<<endl;
}

Foam::label Foam::radiation::fsckContainer::bandGet() const
{
    Info<<"Band retrieved "<<bandID<<endl;
    return bandID;
}

void Foam::radiation::fsckContainer::fsckGetCellka
(
    std::vector<std::vector<double>>& kcell,
    std::vector<std::vector<double>>& acell
)
{
    kcell.assign(mesh_.cells().size(),std::vector<double>(quadPoints,0.0));
    acell.assign(mesh_.cells().size(),std::vector<double>(quadPoints,0.0));

    const volScalarField& T=mesh_.objectRegistry::lookupObject<volScalarField>("T");
    
    //const volScalarField& CO2=mesh_.objectRegistry::lookupObject<volScalarField>("CO2");
    //const volScalarField& H2O=mesh_.objectRegistry::lookupObject<volScalarField>("H2O");
    //const volScalarField& CO=mesh_.objectRegistry::lookupObject<volScalarField>("CO");
    //const volScalarField& CH4=mesh_.objectRegistry::lookupObject<volScalarField>("CH4");
    //const volScalarField& NH3=mesh_.objectRegistry::lookupObject<volScalarField>("NH3");
    // Load the mole fractions with format "X_speciesName"
    // Requires the activation of moleFraction functionObject in controlDict
    const volScalarField& CO2=mesh_.objectRegistry::lookupObject<volScalarField>("X_CO2");
    const volScalarField& H2O=mesh_.objectRegistry::lookupObject<volScalarField>("X_H2O");
    const volScalarField& CO=mesh_.objectRegistry::lookupObject<volScalarField>("X_CO");
    const volScalarField& CH4=mesh_.objectRegistry::lookupObject<volScalarField>("X_CH4");
    const volScalarField& NH3=mesh_.objectRegistry::lookupObject<volScalarField>("X_NH3");
    /*bool presentCO = false;
    if (mesh_.objectRegistry::foundObject<volScalarField>("CO"))
    {
        presentCO = true;
    }*/

    // Compute Reference temperature (Volume average of the domain)
    // Reference Temperature and vector with to load all data into memory
    if (calcTref)
    {
        Tref = 0.0;
        scalar totalVolume = 0.0;
        forAll(T,celli)
        {
            Tref = Tref + T[celli]*mesh_.V()[celli];
            totalVolume = totalVolume + mesh_.V()[celli];
        }
        reduce(Tref, sumOp<double>());
        reduce(totalVolume, sumOp<scalar>());
        Tref = Tref/totalVolume;
    }

    Info << "Reference temperature is set to - " << Tref << endl;
    Info << "Background kappa value is - " << kappaBg.value() << endl;
    
    // Compute k values and scaling factor a for all the cells
    forAll(mesh_.cells(),celli)
    {
        
        GasMixInfo Mix_Info;
        Mix_Info.T = T[celli];
        Mix_Info.xCO = CO[celli]; //0.0;
        Mix_Info.xCO2 = CO2[celli];
        Mix_Info.xH2O = H2O[celli];
        // Added for Ammonia Database
        Mix_Info.xCH4 = CH4[celli];//0.0; // Not used in this model
        Mix_Info.xNH3 = NH3[celli];//0.0; // Not used in this model
        
        //std::vector<double> k;
        //std::vector<double> a;
        if (mixingModel)
        {
            get_ka_mix(Mix_Info,Tref,kq,kq2,gBase,weightBase,kBase,aBase,mix_scheme);
        }
        else
        {
            get_ka(Mix_Info,Tref,kq,gBase,weightBase,isV2,kBase,aBase);
        }

        int points = kBase.size();
        if (points != quadPoints)
        {
            k.assign(quadPoints,0.0);
            a.assign(quadPoints,0.0);
            lin_interp(gBase,kBase,g,k);
            lin_interp(gBase,aBase,g,a);
            points = quadPoints;
        }
        else
        {
            k = kBase;
            a = aBase;
        }
        
            /*Info << endl << "k values for above inputs " << endl;
            for (int i = 1; i <= points; ++i)
            {
                Info << i << "  " << g[i-1] << "  " << k[i-1] << endl;
            }
            Info << endl << "a values for above inputs " << endl;
            for (int i = 1; i <= points; ++i)
            {
                Info << i << "  " << g[i-1] << "  " << a[i-1] << endl;
            }
            Info << "Working till here " << endl;*/
        for (int i = 1; i <= points; ++i)
        {
            // Convert cm-1 to m-1;
            kcell[celli][i-1] = k[i-1]*100.0+kappaBg.value();
            acell[celli][i-1] = a[i-1];
            //acellInternal[celli][i-1] = a[i-1];
        }
    }

    // Compute the scaling factor for all the boundary faces
    forAll(mesh_.boundaryMesh(),patchi)
    {
        // skip empty patches that reduce dimensionality of the problem
        // Temperature field on empty patches will lead to seg fault!
        if (mesh_.boundaryMesh()[patchi].type()!="empty")
        {
            //Pout << "Boundary type is " << mesh_.boundaryMesh()[patchi].type() << endl;
            
            const scalarField& boundaryT = T.boundaryField()[patchi];

            forAll(mesh_.boundaryMesh()[patchi],facei)
            {
                GasMixInfo Mix_Info;
                const label bFaceIdx = mesh_.boundaryMesh()[patchi].start()+facei-mesh_.nInternalFaces();
                //Info<<"The boundary face index is "<<bFaceIdx<<endl;
                Mix_Info.T = boundaryT[facei];

                const label bCell = mesh_.boundaryMesh()[patchi].faceCells()[facei];
                Mix_Info.xCO = CO[bCell]; //0.0;
                Mix_Info.xCO2 = CO2[bCell];
                Mix_Info.xH2O = H2O[bCell];
                // Added for Ammonia Database
                Mix_Info.xCH4 = CH4[bCell];//0.0; // Not used in this model
                Mix_Info.xNH3 = NH3[bCell];//0.0; // Not used in this model

                //Pout << "Boundary index " << bFaceIdx << endl;
                //Pout << "T - " << Mix_Info.T << " ; xCO -" << Mix_Info.xCO << " ; xCO2 - " << Mix_Info.xCO2 << " ; xH2O - " << Mix_Info.xH2O << endl;
                //Pout << "Scaling factors are - ";
                if (mixingModel)
                {
                    get_ka_mix(Mix_Info,Tref,kq,kq2,gBase,weightBase,kBase,aBase,mix_scheme);
                }
                else
                {
                    get_ka(Mix_Info,Tref,kq,gBase,weightBase,isV2,kBase,aBase);
                }

                int points = kBase.size();
                if (points != quadPoints)
                {
                    k.assign(quadPoints,0.0);
                    a.assign(quadPoints,0.0);
                    lin_interp(gBase,kBase,g,k);
                    lin_interp(gBase,aBase,g,a);
                    points = quadPoints;
                }
                else
                {
                    k = kBase;
                    a = aBase;
                }

                for (int i = 1; i <= points; ++i)
                {
                    afaceInternal[bFaceIdx][i-1] = a[i-1];
                    //acellInternal[celli][i-1] = a[i-1];
                    //Pout << a[i-1] << "; ";
                }
                //Pout << endl;
            }
        }
    }
}

void Foam::radiation::fsckContainer::fsckGetBFka
(
    std::vector<double>& aBF
) const
{
    const label bFaceSize = mesh_.faces().size()-mesh_.nInternalFaces();
    aBF.assign(bFaceSize,0.0);
    
    forAll(mesh_.boundaryMesh(),patchi)
    {
        forAll(mesh_.boundaryMesh()[patchi],facei)
        {
            const label bFaceIdx = mesh_.boundaryMesh()[patchi].start()+facei-mesh_.nInternalFaces();
            aBF[bFaceIdx] = afaceInternal[bFaceIdx][bandID]; // for bandID -> C++ index start from 0
        }
    }
}

std::vector<double> Foam::radiation::fsckContainer::getWeight() const
{
    return weight;
}

void Foam::radiation::fsckContainer::lin_interp
(
    const std::vector<double>& x_in,
    const std::vector<double>& y_in,
    std::vector<double>& x_out,
    std::vector<double>& y_out
)
{
    int n_in = x_in.size();
    int n_out = x_out.size();
    int l_indx, u_indx;
    // Linear interpolation function for quadrature points
    /*for ii=1,n_out
        do jj=1,n_data-1
            if (x_out(ii)<x_data(jj+1) .AND. x_out(ii)>x_data(jj)) then
                u_indx=jj+1
                l_indx=jj
                exit
            end if
        end do
        y_out(ii)=y_data(l_indx)+(y_data(u_indx)-y_data(l_indx)) &
            *(x_out(ii)-x_data(l_indx)) &
            /(x_data(u_indx)-x_data(l_indx))
    end do*/
    // Converted above Fortran to C++ code
    for (int ii = 0; ii < n_out; ++ii)
    {
        for (int jj = 0; jj < n_in-1; ++jj)
        {
            if (x_out[ii] < x_in[jj+1] && x_out[ii] > x_in[jj])
            {
                u_indx = jj+1;
                l_indx = jj;
                break;
            }
        }
        y_out[ii] = y_in[l_indx] + (y_in[u_indx]-y_in[l_indx])*(x_out[ii]-x_in[l_indx])/(x_in[u_indx]-x_in[l_indx]);
    }
}
// ************************************************************************* //
