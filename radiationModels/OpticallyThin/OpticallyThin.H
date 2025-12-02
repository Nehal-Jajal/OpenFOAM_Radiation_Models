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

Class
    Foam::radiation::OpticallyThin

Group
    grpRadiationModels

Description
    The optically thin radiation model that ignores reabsorption of gas
    Only the emission part is considered.
    Works well for combustion applications where optical thickness, tau is
    small

SourceFiles
    OpticallyThin.C

\*---------------------------------------------------------------------------*/

#ifndef radiation_OpticallyThin_H
#define radiation_OpticallyThin_H

#include "radiationModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                           Class P1Custom Declaration
\*---------------------------------------------------------------------------*/

class OpticallyThin
:
    public radiationModel
{
    // Private data

        //- Absorption coefficient
        volScalarField a_;
        
        //- Radiation Souce Term
        volScalarField delqr_;

        //- Rp term computed within calculate
        volScalarField Rp_;


    // Private Member Functions

        //- No copy construct
        OpticallyThin(const OpticallyThin&) = delete;

        //- No copy assignment
        void operator=(const OpticallyThin&) = delete;


public:

    //- Runtime type information
    TypeName("OpticallyThin");


    // Constructors

        //- Construct from components
        OpticallyThin(const volScalarField& T);

        //- Construct from components
        OpticallyThin(const dictionary& dict, const volScalarField& T);


    //- Destructor
    virtual ~OpticallyThin() = default;


    // Member functions

        // Edit

            //- Solve radiation equation(s)
            void calculate();

            //- Read radiation properties dictionary
            bool read();


        // Access

            //- Source term component (for power of T^4)
            virtual tmp<volScalarField> Rp() const;

            //- Source term component (constant)
            virtual tmp<volScalarField::Internal> Ru() const;

            //- Number of bands
            label nBands() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
