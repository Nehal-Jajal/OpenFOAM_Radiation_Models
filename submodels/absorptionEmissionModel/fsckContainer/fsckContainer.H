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
    Foam::radiation::fsckContainer

Group
    grpRadiationAbsorptionEmissionSubModels

Description
    FSCK radiation absorption and emission coefficients for continuous
    phase

SourceFiles
    fsckContainer.C

\*---------------------------------------------------------------------------*/

#ifndef radiation_fsckContainer_H
#define radiation_fsckContainer_H

#include "fskTableFunction.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                  Class fsckContainer Declaration
\*---------------------------------------------------------------------------*/

class fsckContainer
{
    // Private data

        //- Use singleton approach for generating one object of this class from fsckAbsorptionEmission
        static fsckContainer* thisContainer;

        //- Absorption model dictionary
        const dictionary coeffsDict_;

        //- Mesh
        const fvMesh& mesh_;

        //- Input Reference Temperature / Enter without dimensions
        dimensionedScalar Tref_input;

        //- Input Background Absorption Coefficient (for stabilizing P1 model)
        dimensionedScalar kappaBg;

        //- Current Band
        label bandID;

        //- Number of quadrature points
        label quadPoints;

        //// FSCK Table v3 Data
        // Reference Temperature and vector with to load all data into memory
        double Tref;
        std::vector<float> kq;
        // Calculate Tref each iteration as the volume average of the domain
        bool calcTref;
        // Quadrature points and weights baseline
        std::vector<double> gBase;
        std::vector<double> weightBase;
        // Quadrature points and weights
        std::vector<double> g;
        std::vector<double> weight;
        // Vectors to extract pointwise k and a values from getka function call
        // Baseline
        std::vector<double> kBase;
        std::vector<double> aBase;
        // For input quadrature points
        std::vector<double> k;
        std::vector<double> a;
        // 2D Vectors for storing k and a values for all quadrature points at all cells
        std::vector<std::vector<double>> kcellInternal;
        std::vector<std::vector<double>> acellInternal;
        // 2D Vectors for storing k and a values for all quadrature points at all faces
        //std::vector<std::vector<double>> kfaceInternal;   // k-values not needed at boundary faces
        std::vector<std::vector<double>> afaceInternal;

        //// FSCK Table with Ammonia Database (used with V2 functions)
        // Set flag to true for using ammonia database
        bool isV2;

        // Set flag for using mixing models
        bool mixingModel;
        int mix_scheme; // 1 for multiplicative, 2 for superposition, 3 for ModestRiazzi
        std::vector<float> kq2;

        // Constructors
        //- Construct from components (Nade private for singleton approach)
        fsckContainer(const dictionary& dict, const fvMesh& mesh, double Tref_in);

public:

    //- Runtime type information
    TypeName("fsckContainer");

    //- Destructor
    virtual ~fsckContainer() = default;


    // Member Functions
    // Create one instance through 
    static fsckContainer* getInstance();
    static fsckContainer* getInstance(const dictionary& dict, const fvMesh& mesh, double Tref_in);

    // FSCK Setter and Getter

    void bandSet(label bandI);
    label bandGet() const;

    // Extract k and a values for all cells values
    void fsckGetCellka
    (
        std::vector<std::vector<double>>& kcell,
        std::vector<std::vector<double>>& acell
    );

    // Extract k and a values for all boundary faces (BF)
    void fsckGetBFka
    (
        std::vector<double>& aBF
    ) const;

    std::vector<double> getWeight() const;

    void lin_interp
    (
        const std::vector<double>& x_in,
        const std::vector<double>& y_in,
        std::vector<double>& x_out,
        std::vector<double>& y_out
    );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
