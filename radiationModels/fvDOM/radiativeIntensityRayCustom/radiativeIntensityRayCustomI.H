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

inline const Foam::volScalarField&
Foam::radiation::radiativeIntensityRayCustom::I() const
{
    return I_;
}


inline const Foam::volScalarField&
Foam::radiation::radiativeIntensityRayCustom::qr() const
{
    return qr_;
}


inline Foam::volScalarField& Foam::radiation::radiativeIntensityRayCustom::qr()
{
    return qr_;
}

inline const Foam::volScalarField& Foam::radiation::
radiativeIntensityRayCustom::qin() const
{
    return qin_;
}


inline Foam::volScalarField& Foam::radiation::radiativeIntensityRayCustom::qin()
{
    return qin_;
}


inline const Foam::volScalarField& Foam::radiation::
radiativeIntensityRayCustom::qem() const
{
    return qem_;
}


inline Foam::volScalarField& Foam::radiation::radiativeIntensityRayCustom::qem()
{
    return qem_;
}


inline const Foam::vector& Foam::radiation::radiativeIntensityRayCustom::d() const
{
    return d_;
}


inline const Foam::vector& Foam::radiation::radiativeIntensityRayCustom::dAve() const
{
    return dAve_;
}


inline Foam::vector& Foam::radiation::radiativeIntensityRayCustom::d()
{
    return d_;
}


inline Foam::vector& Foam::radiation::radiativeIntensityRayCustom::dAve()
{
    return dAve_;
}


inline Foam::scalar Foam::radiation::radiativeIntensityRayCustom::nLambda() const
{
    return nLambda_;
}


inline Foam::scalar Foam::radiation::radiativeIntensityRayCustom::phi() const
{
    return phi_;
}


inline Foam::scalar Foam::radiation::radiativeIntensityRayCustom::theta() const
{
    return theta_;
}


inline Foam::scalar Foam::radiation::radiativeIntensityRayCustom::omega() const
{
    return omega_;
}


inline const Foam::volScalarField&
Foam::radiation::radiativeIntensityRayCustom::ILambda
(
    const label lambdaI
) const
{
    return ILambda_[lambdaI];
}


// ************************************************************************* //
