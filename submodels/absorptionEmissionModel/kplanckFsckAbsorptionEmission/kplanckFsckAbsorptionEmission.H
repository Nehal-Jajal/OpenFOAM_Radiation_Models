/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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
    Foam::radiation::constantAbsorptionEmission

Group
    grpRadiationAbsorptionEmissionSubModels

Description
    Constant radiation absorption and emission coefficients for continuous
    phase

SourceFiles
    constantAbsorptionEmission.C

\*---------------------------------------------------------------------------*/

#ifndef radiation_kplanckFsckAbsorptionEmission_H
#define radiation_kplanckFsckAbsorptionEmission_H

#include "absorptionEmissionModel.H"
#include "volFields.H"
#include "fskTableFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                  Class kplanckAbsorptionEmission Declaration
\*---------------------------------------------------------------------------*/

class kplanckFsckAbsorptionEmission
:
    public absorptionEmissionModel
{
    // Private data

        //- Absorption model dictionary
        dictionary coeffsDict_;

        //- Absorption coefficient / [1/m]
        volScalarField a_;

        //- Emission coefficient / [1/m]
        volScalarField e_;

        //- Emission contribution / [kg/(m s^3)]
        volScalarField E_;

        //- Input background absorption coefficient (for stabilizing P1 model)
        dimensionedScalar kappaBg;

        //// FSCK Table v3 Data
        std::vector<float> kq;
        // Quadrature points and weights
        std::vector<double> g;
        std::vector<double> weight;
        // Vectors to extract pointwise k and a values from getka function call
        std::vector<double> k;
        std::vector<double> a;

        //// Switch for using Ammonia Database
        // Set flag to true for using ammonia database
        bool isV2;
        // Reference Temperature for the fsck Table
        dimensionedScalar Tref_input;
        double Tref;
        bool calcTref;

        // Set flag for using mixing models
        bool mixingModel;
        int mix_scheme; // 1 for multiplicative, 2 for superposition, 3 for ModestRiazzi
        std::vector<float> kq2;


public:

    //- Runtime type information
    TypeName("kplanckFsckAbsorptionEmission");


    // Constructors

        //- Construct from components
        kplanckFsckAbsorptionEmission(const dictionary& dict, const fvMesh& mesh);


    //- Destructor
    virtual ~kplanckFsckAbsorptionEmission() = default;


    // Member Functions

        // Access

            // Absorption coefficient

                //- Absorption coefficient for continuous phase
                tmp<volScalarField> aCont(const label bandI = 0) const;


            // Emission coefficient

                //- Emission coefficient for continuous phase
                tmp<volScalarField> eCont(const label bandI = 0) const;


            // Emission contribution

                //- Emission contribution for continuous phase
                tmp<volScalarField> ECont(const label bandI = 0) const;


    // Member Functions

        inline bool isGrey() const
        {
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
