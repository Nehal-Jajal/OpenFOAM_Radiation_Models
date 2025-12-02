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
    Foam::radiation::fsckAbsorptionEmission

Group
    grpRadiationAbsorptionEmissionSubModels

Description
    FSCK radiation absorption and emission coefficients for continuous
    phase

SourceFiles
    fsckAbsorptionEmission.C

\*---------------------------------------------------------------------------*/

#ifndef radiation_fsckAbsorptionEmission_H
#define radiation_fsckAbsorptionEmission_H

#include "absorptionEmissionModel.H"
#include "fskTableFunction.H"
#include "fsckContainer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                  Class fsckAbsorptionEmission Declaration
\*---------------------------------------------------------------------------*/

class fsckAbsorptionEmission
:
    public absorptionEmissionModel
{
    // Private data

        //- Absorption model dictionary
        dictionary coeffsDict_;

        //- Absorption coefficient / [1/m]
        dimensionedScalar a_;

        //- Emission coefficient / [1/m]
        dimensionedScalar e_;

        //- Emission contribution / [kg/(m s^3)]
        dimensionedScalar E_;

        //- Input Reference Temperature / Enter without dimensions
        dimensionedScalar Tref_input;

        //// FSCK Table v3 Data
        // Reference Temperature and vector with to load all data into memory
        double Tref;       
        

public:

    //- Runtime type information
    TypeName("fsckAbsorptionEmission");


    // Constructors

        //- Construct from components
        fsckAbsorptionEmission(const dictionary& dict, const fvMesh& mesh);


    //- Destructor
    virtual ~fsckAbsorptionEmission() = default;


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
