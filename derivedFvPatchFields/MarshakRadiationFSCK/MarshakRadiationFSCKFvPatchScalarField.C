/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2016,2020 OpenCFD Ltd.
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

#include "MarshakRadiationFSCKFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "radiationModel.H"
#include "physicoChemicalConstants.H"
#include "boundaryRadiationProperties.H"
#include "fsckContainer.H"

// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * *//

/*void Foam::radiation::MarshakRadiationFvPatchScalarField::initialise()
{
    if (this->found("absorptionEmissionModel"))
    {
        absorptionEmission_.reset
        (
            absorptionEmissionModel::New(*this,mesh_).ptr()
        );
    }
}*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::MarshakRadiationFSCKFvPatchScalarField::
MarshakRadiationFSCKFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    TName_("T")
    //absorptionEmission_(nullptr)
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 0.0;
}


Foam::radiation::MarshakRadiationFSCKFvPatchScalarField::
MarshakRadiationFSCKFvPatchScalarField
(
    const MarshakRadiationFSCKFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    TName_(ptf.TName_)
    //absorptionEmission_(nullptr)
{}


Foam::radiation::MarshakRadiationFSCKFvPatchScalarField::
MarshakRadiationFSCKFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    TName_(dict.getOrDefault<word>("T", "T"))
    //absorptionEmission_(nullptr)
{
    if (dict.found("value"))
    {
        refValue() = scalarField("value", dict, p.size());
    }
    else
    {
        refValue() = 0.0;
    }

    // zero gradient
    refGrad() = 0.0;

    valueFraction() = 1.0;

    fvPatchScalarField::operator=(refValue());
}


Foam::radiation::MarshakRadiationFSCKFvPatchScalarField::
MarshakRadiationFSCKFvPatchScalarField
(
    const MarshakRadiationFSCKFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    TName_(ptf.TName_)
    //absorptionEmission_(nullptr)
{}


Foam::radiation::MarshakRadiationFSCKFvPatchScalarField::
MarshakRadiationFSCKFvPatchScalarField
(
    const MarshakRadiationFSCKFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    TName_(ptf.TName_)
    //absorptionEmission_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::radiation::MarshakRadiationFSCKFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

    // Temperature field
    const scalarField& Tp =
        patch().lookupPatchField<volScalarField, scalar>(TName_);

    //const label fieldSize = Tp.size();
    //Info<<"Number of boundary patch faces is "<< fieldSize << endl;
    //Info<<"Number of boundary faces is "<< patch().boundaryMesh().size() <<endl;

    // Temperature Vol field
    //const volScalarField& T = patch().boundaryMesh().mesh().objectRegistry::lookupObject<volScalarField>("T");

    // *FSCK* Modification
    // Wall @ 300K
    //scalar p1_af[8] = {1.410031,0.359589,0.960889,1.946702,3.036567,0.155868,9.230135e-2,4.649323e-2};
    // Wall @ 1500K
    //scalar p1_af[8] = {0.831306,0.882217,0.920929,0.756753,0.686080,0.635293,0.665711,0.618451};
    //const volScalarField& T=this->mesh().objectRegistry::lookupObject<volScalarField>("T");
    //autoPtr<radiation::radiationModel> radmodel(radiation::radiationModel::New(T));
    const radiationModel& radmodel = patch().boundaryMesh().mesh().objectRegistry::lookupObject<radiationModel>("radiationProperties");
    Foam::radiation::fsckContainer* fsck = Foam::radiation::fsckContainer::getInstance();
    /*absorptionEmission_.reset
    (
        absorptionEmissionModel::New(radmodel,patch().boundaryMesh().mesh()).ptr()
    );*/
    //const scalar a_f = radmodel->af();
    //const label bandI = radmodel->absorptionEmission().bandGet();
    //const label bandI = radmodel.absorptionEmission().bandGet();
    //const label bandI = 0;
    std::vector<double> aBF;
    //radmodel.absorptionEmission().fsckGetBFka(aBF);
    fsck->fsckGetBFka(aBF);

    // Re-calc reference value
    //refValue() = 4.0*constant::physicoChemical::sigma.value()*pow4(Tp)*p1_af[bandI];
    //Info<<"The scaling factor for band is - "<<p1_af[bandI]<<endl;
    const label patchStart = patch().start()-patch().boundaryMesh().mesh().nInternalFaces();
    forAll(patch(),facei)
    {
        //const label bFaceIdx = mesh_.boundaryMesh()[patchi].start()+facei-mesh_.nInternalFaces();
        //const label bFaceIdx = patch().start()+facei-patch().boundaryMesh().mesh().nInternalFaces();
        const label bFaceIdx = patchStart+facei;
        refValue()[facei] = 4.0*constant::physicoChemical::sigma.value()*pow4(Tp[facei])*aBF[bFaceIdx];
        //Pout << "Boundary face index is " << bFaceIdx << endl;
        //Info << "Temperature is " << Tp[facei] << endl;
        //Info << "Scaling factor is " << aBF[bFaceIdx] << endl;
        //Pout << "Ref value is " << refValue()[facei] << endl;
        //Pout << "Patch name is " <<patch().name() << " and patch type is " << patch().type() << endl;
    }
    

    // Diffusion coefficient - created by radiation model's ::updateCoeffs()
    const scalarField& gamma =
        patch().lookupPatchField<volScalarField, scalar>("gammaRad");

    const boundaryRadiationProperties& boundaryRadiation =
        boundaryRadiationProperties::New(internalField().mesh());

    const tmp<scalarField> temissivity
    (
        boundaryRadiation.emissivity(patch().index())
    );

    const scalarField& emissivity = temissivity();

    const scalarField Ep(emissivity/(2.0*(scalar(2) - emissivity)));

    // Set value fraction
    valueFraction() = 1.0/(1.0 + gamma*patch().deltaCoeffs()/Ep);

    // Restore tag
    UPstream::msgType() = oldTag;

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::radiation::MarshakRadiationFSCKFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    os.writeEntryIfDifferent<word>("T", "T", TName_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{
    makePatchTypeField
    (
        fvPatchScalarField,
        MarshakRadiationFSCKFvPatchScalarField
    );
}
}

// ************************************************************************* //
