/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "fvDOM_FSCK.H"
#include "absorptionEmissionModel.H"
#include "scatterModel.H"
#include "constants.H"
#include "unitConversion.H"
#include "fvm.H"
#include "addToRunTimeSelectionTable.H"
#include "fsckContainer.H"

using namespace Foam::constant;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(fvDOM_FSCK, 0);
        addToRadiationRunTimeSelectionTables(fvDOM_FSCK);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::radiation::fvDOM_FSCK::rotateInitialRays(const vector& sunDir)
{
    // Rotate Y spherical cordinates to Sun direction.
    // Solid angles on the equator are better fit for planar radiation
    const tensor coordRot = rotationTensor(vector(0, 1, 0), sunDir);

    forAll(IRay_, rayId)
    {
        IRay_[rayId].dAve() = coordRot & IRay_[rayId].dAve();
        IRay_[rayId].d() = coordRot & IRay_[rayId].d();
    }
}


void Foam::radiation::fvDOM_FSCK:: alignClosestRayToSun(const vector& sunDir)
{
    label SunRayId(-1);
    scalar maxSunRay = -GREAT;

    // Looking for the ray closest to the Sun direction
    forAll(IRay_, rayId)
    {
        const vector& iD = IRay_[rayId].d();
        scalar dir = sunDir & iD;
        if (dir > maxSunRay)
        {
            maxSunRay = dir;
            SunRayId = rayId;
        }
    }

    // Second rotation to align colimated radiation with the closest ray
    const tensor coordRot = rotationTensor(IRay_[SunRayId].d(), sunDir);

    forAll(IRay_, rayId)
    {
        IRay_[rayId].dAve() = coordRot & IRay_[rayId].dAve();
        IRay_[rayId].d() = coordRot & IRay_[rayId].d();
    }

    Info << "Sun direction : " << sunDir << nl << endl;
    Info << "Sun ray ID : " << SunRayId << nl << endl;
}


void Foam::radiation::fvDOM_FSCK::updateRaysDir()
{
    solarCalculator_->correctSunDirection();
    const vector sunDir = solarCalculator_->direction();

    // First iteration
    if (updateTimeIndex_ == 0)
    {
        rotateInitialRays(sunDir);
        alignClosestRayToSun(sunDir);
    }
    else if (updateTimeIndex_ > 0)
    {
        alignClosestRayToSun(sunDir);
    }
}


void Foam::radiation::fvDOM_FSCK::initialise()
{
    coeffs_.readIfPresent("useExternalBeam", useExternalBeam_);

    if (useExternalBeam_)
    {
        spectralDistributions_.reset
        (
            Function1<scalarField>::New
            (
                "spectralDistribution",
                coeffs_,
                &mesh_
            )
        );

        spectralDistribution_ =
            spectralDistributions_->value(mesh_.time().value());

        spectralDistribution_ =
            spectralDistribution_/sum(spectralDistribution_);

        const dictionary& solarDict = this->subDict("solarCalculatorCoeffs");
        solarCalculator_.reset(new solarCalculator(solarDict, mesh_));

        if (mesh_.nSolutionD() != 3)
        {
            FatalErrorInFunction
                << "External beam model only available in 3D meshes "
                << abort(FatalError);
        }

        if (solarCalculator_->diffuseSolarRad() > 0)
        {
            FatalErrorInFunction
                << "External beam model does not support Diffuse "
                << "Solar Radiation. Set diffuseSolarRad to zero"
                << abort(FatalError);
        }
        if (spectralDistribution_.size() != nLambda_)
        {
            FatalErrorInFunction
                << "The epectral energy distribution has different bands "
                << "than the absoprtivity model "
                << abort(FatalError);
        }
    }

    // 3D
    if (mesh_.nSolutionD() == 3)
    {
        nRay_ = 4*nPhi_*nTheta_;

        IRay_.setSize(nRay_);

        const scalar deltaPhi = pi/(2*nPhi_);
        const scalar deltaTheta = pi/nTheta_;

        label i = 0;

        for (label n = 1; n <= nTheta_; n++)
        {
            for (label m = 1; m <= 4*nPhi_; m++)
            {
                scalar thetai = (2*n - 1)*deltaTheta/2.0;
                scalar phii = (2*m - 1)*deltaPhi/2.0;

                IRay_.set
                (
                    i,
                    new radiativeIntensityRay_FSCK
                    (
                        *this,
                        mesh_,
                        phii,
                        thetai,
                        deltaPhi,
                        deltaTheta,
                        nLambda_,
                        *absorptionEmission_,
                        blackBody_,
                        i
                    )
                );
                i++;
            }
        }
    }
    // 2D
    else if (mesh_.nSolutionD() == 2)
    {
        const scalar thetai = piByTwo;
        const scalar deltaTheta = pi;
        nRay_ = 4*nPhi_;
        IRay_.setSize(nRay_);
        const scalar deltaPhi = pi/(2.0*nPhi_);
        label i = 0;
        for (label m = 1; m <= 4*nPhi_; m++)
        {
            const scalar phii = (2*m - 1)*deltaPhi/2.0;
            IRay_.set
            (
                i,
                new radiativeIntensityRay_FSCK
                (
                    *this,
                    mesh_,
                    phii,
                    thetai,
                    deltaPhi,
                    deltaTheta,
                    nLambda_,
                    *absorptionEmission_,
                    blackBody_,
                    i
                )
            );
            i++;
        }
    }
    // 1D
    else
    {
        const scalar thetai = piByTwo;
        const scalar deltaTheta = pi;
        nRay_ = 2;
        IRay_.setSize(nRay_);
        const scalar deltaPhi = pi;
        label i = 0;
        for (label m = 1; m <= 2; m++)
        {
            const scalar phii = (2*m - 1)*deltaPhi/2.0;
            IRay_.set
            (
                i,
                new radiativeIntensityRay_FSCK
                (
                    *this,
                    mesh_,
                    phii,
                    thetai,
                    deltaPhi,
                    deltaTheta,
                    nLambda_,
                    *absorptionEmission_,
                    blackBody_,
                    i
                )
            );
            i++;
        }
    }


    // Construct absorption field for each wavelength
    forAll(aLambda_, lambdaI)
    {
        aLambda_.set
        (
            lambdaI,
            new volScalarField
            (
                IOobject
                (
                    "aLambda_" + Foam::name(lambdaI) ,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                a_
            )
        );
    }

    Info<< "fvDOM : Allocated " << IRay_.size()
        << " rays with average orientation:" << nl;

    if (useExternalBeam_)
    {
        // Rotate rays for Sun direction
        updateRaysDir();
    }

    scalar totalOmega = 0;
    forAll(IRay_, rayId)
    {
        if (omegaMax_ <  IRay_[rayId].omega())
        {
            omegaMax_ = IRay_[rayId].omega();
        }
        totalOmega += IRay_[rayId].omega();
        Info<< '\t' << IRay_[rayId].I().name() << " : " << "dAve : "
            << '\t' << IRay_[rayId].dAve() << " : " << "omega : "
            << '\t' << IRay_[rayId].omega() << " : " << "d : "
            << '\t' << IRay_[rayId].d() << nl;
    }

    Info << "Total omega : " << totalOmega << endl;

    Info<< endl;

    coeffs_.readIfPresent("useSolarLoad", useSolarLoad_);

    if (useSolarLoad_)
    {
        if (useExternalBeam_)
        {
            FatalErrorInFunction
                << "External beam with fvDOM can not be used "
                << "with the solar load model"
                << abort(FatalError);
        }
        const dictionary& solarDict = this->subDict("solarLoadCoeffs");
        solarLoad_.reset(new solarLoad(solarDict, T_));

        if (solarLoad_->nBands() != this->nBands())
        {
            FatalErrorInFunction
                << "Requested solar radiation with fvDOM. Using "
                << "different number of bands for the solar load is not allowed"
                << abort(FatalError);
        }

        Info<< "Creating Solar Load Model " << nl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::fvDOM_FSCK::fvDOM_FSCK(const volScalarField& T)
:
    radiationModel(typeName, T),
    G_
    (
        IOobject
        (
            "G",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,//NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_//,
        //dimensionedScalar(dimMass/pow3(dimTime), Zero)  ! MODIFIED FOR HYBRID SOLVER
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
    qem_
    (
        IOobject
        (
            "qem",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/pow3(dimTime), Zero)
    ),
    qin_
    (
        IOobject
        (
            "qin",
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
    ),
    nTheta_(coeffs_.get<label>("nTheta")),
    nPhi_(coeffs_.get<label>("nPhi")),
    nRay_(0),
    nLambda_(absorptionEmission_->nBands()),
    aLambda_(nLambda_),
    blackBody_(nLambda_, T),
    IRay_(0),
    tolerance_
    (
        coeffs_.getOrDefaultCompat<scalar>
        (
            "tolerance",
            {{"convergence", 1712}},
            0
        )
    ),
    maxIter_(coeffs_.getOrDefault<label>("maxIter", 50)),
    omegaMax_(0),
    useSolarLoad_(false),
    solarLoad_(),
    meshOrientation_
    (
        coeffs_.getOrDefault<vector>("meshOrientation", Zero)
    ),
    useExternalBeam_(false),
    spectralDistribution_(),
    spectralDistributions_(),
    solarCalculator_(),
    updateTimeIndex_(0),
    useHybridSolver_(coeffs_.getOrDefault<bool>("Hybrid", false)),
    hyb_L_(coeffs_.getOrDefault<dimensionedScalar>("lengthScale", 0.0)),
    tau_u_(coeffs_.getOrDefault<dimensionedScalar>("tauU", 0.0)),
    tau_l_(coeffs_.getOrDefault<dimensionedScalar>("tauL", 0.0))
{
    initialise();
}


Foam::radiation::fvDOM_FSCK::fvDOM_FSCK
(
    const dictionary& dict,
    const volScalarField& T
)
:
    radiationModel(typeName, dict, T),
    G_
    (
        IOobject
        (
            "G",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/pow3(dimTime), Zero)
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
    qem_
    (
        IOobject
        (
            "qem",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/pow3(dimTime), Zero)
    ),
    qin_
    (
        IOobject
        (
            "qin",
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
    ),
    nTheta_(coeffs_.get<label>("nTheta")),
    nPhi_(coeffs_.get<label>("nPhi")),
    nRay_(0),
    nLambda_(absorptionEmission_->nBands()),
    aLambda_(nLambda_),
    blackBody_(nLambda_, T),
    IRay_(0),
    tolerance_
    (
        coeffs_.getOrDefaultCompat<scalar>
        (
            "tolerance",
            {{"convergence", 1712}},
            0
        )
    ),
    maxIter_(coeffs_.getOrDefault<label>("maxIter", 50)),
    omegaMax_(0),
    useSolarLoad_(false),
    solarLoad_(),
    meshOrientation_
    (
        coeffs_.getOrDefault<vector>("meshOrientation", Zero)
    ),
    useExternalBeam_(false),
    spectralDistribution_(),
    spectralDistributions_(),
    solarCalculator_(),
    updateTimeIndex_(0),
    useHybridSolver_(coeffs_.getOrDefault<bool>("Hybrid", false)),
    hyb_L_(coeffs_.getOrDefault<dimensionedScalar>("lengthScale", 0.0)),
    tau_u_(coeffs_.getOrDefault<dimensionedScalar>("tauU", 1.0)),
    tau_l_(coeffs_.getOrDefault<dimensionedScalar>("tauL", -1.0))
{
    initialise();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::radiation::fvDOM_FSCK::read()
{
    if (radiationModel::read())
    {
        // Only reading solution parameters - not changing ray geometry
        coeffs_.readIfPresentCompat
        (
            "tolerance", {{"convergence", 1712}}, tolerance_
        );
        coeffs_.readIfPresent("maxIter", maxIter_);

        return true;
    }

    return false;
}


void Foam::radiation::fvDOM_FSCK::calculate()
{
    absorptionEmission_->correct(a_, aLambda_);

    updateBlackBodyEmission();

    if (useSolarLoad_)
    {
        solarLoad_->calculate();
    }

    if (useExternalBeam_)
    {
        switch (solarCalculator_->sunDirectionModel())
        {
            case solarCalculator::mSunDirConstant:
            {
                break;
            }
            case solarCalculator::mSunDirTracking:
            {
                label updateIndex = label
                (
                    mesh_.time().value()
                   /solarCalculator_->sunTrackingUpdateInterval()
                );

                if (updateIndex > updateTimeIndex_)
                {
                    Info << "Updating Sun position..." << endl;
                    updateTimeIndex_ = updateIndex;
                    updateRaysDir();
                }
                break;
            }
        }
    }

    Info<<"\n Running FSCK version of fvDOM (hardcoded) \n"<<endl;
    //scalar kfsck[8] = {9.767899e-7,5.360105e-6,8.152979e-5,1.340214e-3,2.387425e-2,0.301935,0.986868,3.687812};
    //scalar wq[8] = {0.333888,0.261700,0.180514,0.120834,6.253695e-2,3.095662e-2,8.205659e-3,1.364635e-3};

    std::vector<std::vector<double>> kcell;
    std::vector<std::vector<double>> acell;
    Foam::radiation::fsckContainer* fsck = Foam::radiation::fsckContainer::getInstance();
    //absorptionEmission_->fsckGetCellka(kcell,acell);
    fsck->fsckGetCellka(kcell,acell);
    std::vector<double> weight;
    //weight = absorptionEmission_->getWeight();
    weight = fsck->getWeight();

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

    delqr_ = dimensionedScalar(dimMass/dimLength/pow3(dimTime), Zero);
    Rp_ = dimensionedScalar(dimMass/dimLength/pow3(dimTime)/pow4(dimTemperature), Zero);
    Ru_ = delqr_;
    const label points = weight.size(); // Need to find out how many quadrature points are being used!!!

    // For Hybrid solver
    // Write out whether the solver is hybrid or not
    if (useHybridSolver_)
    {
        Info << "Hybrid solver is being used." << endl;
    }
    else
    {
        Info << "Hybrid solver is not being used." << endl;
    }
    // Labels for counting the number of iterations for each solver
    label fvDomIter = 0;
    label p1Iter = 0;
    // Find out hot and cold spots in domain for Hybrid solver implementation
    // Find the maximum temperature index in the domain
    scalar maxT = max(T_).value();
    //Info << "Maximum temperature in the domain is " << maxT << endl;
    Pout << "Maximum temperature in the domain is " << maxT << endl;
    // Index of the maximum temperature
    //label maxTIdx = T_.findIndex(maxT);
    //Info << "Index of the maximum temperature is " << maxTIdx << " and temperature is (bug check)" << T_[maxTIdx] << endl;

    // Find the minimum temperature index in the domain
    scalar minT = min(T_).value();
    //Info << "Minimum temperature in the domain is " << minT << endl;
    Pout << "Minimum temperature in the domain is " << minT << endl;
    // Index of the minimum temperature
    //label minTIdx = T_.findIndex(minT);
    //Info << "Index of the minimum temperature is " << minTIdx << " and temperature is (bug check)" << T_[minTIdx] << endl;
    labelPair minMaxIDs = findMinMax(T_);
    label minTIdx = minMaxIDs.first();
    label maxTIdx = minMaxIDs.second();
    //Info << "Index of the minimum temperature is " << minTIdx << " and temperature is (bug check)" << T_[minTIdx] << endl;
    //Info << "Index of the maximum temperature is " << maxTIdx << " and temperature is (bug check)" << T_[maxTIdx] << endl;
    Pout << "Index of the minimum temperature is " << minTIdx << " and temperature is (bug check)" << T_[minTIdx] << endl;
    Pout << "Index of the maximum temperature is " << maxTIdx << " and temperature is (bug check)" << T_[maxTIdx] << endl;
    // Global maximum and minimum temperature (needed for parallelization)
    maxT= T_[maxTIdx];
    minT= T_[minTIdx];
    scalar globalMaxT = maxT;
    scalar globalMinT = minT;
    reduce(globalMaxT, maxOp<scalar>());
    reduce(globalMinT, minOp<scalar>());
    // Set up scalars that are set to 1 if the maximum and minimum temperatures are found on the current processor or 0 otherwise
    scalar maxTempOnThisProc = (maxT == globalMaxT) ? 1.0 : 0.0;
    scalar minTempOnThisProc = (minT == globalMinT) ? 1.0 : 0.0;

    for (label j=0; j < points; j++)
    {
        forAll(G_.mesh().cells(),celli)
        {
            a_[celli]=kcell[celli][j];//dimensionedScalar(dimless/dimLength,kcell[celli][j]);
            //e_[celli]=a_[celli]*acell[celli][j];
            aScale[celli]=acell[celli][j];
        }
        a_.correctBoundaryConditions();
        aScale.correctBoundaryConditions();
        
        // Use the index of maximum and minimum temperature to find the maximum and minimum absorption coefficient
        scalar kappa_hyb_h = a_[maxTIdx]*maxTempOnThisProc;
        scalar kappa_hyb_c = a_[minTIdx]*minTempOnThisProc;
        // Find the maximum and minimum absorption coefficient in the domain across all paralllel processors
        reduce(kappa_hyb_h, sumOp<scalar>());
        reduce(kappa_hyb_c, sumOp<scalar>());

        //if ((kappa_hyb_h*hyb_L<crit_u .and. kappa_hyb_h*hyb_L>crit_l) .or. &
        //            (kappa_hyb_c*hyb_L<crit_u .and. kappa_hyb_c*hyb_L>crit_l)) then
        // Converting the above condition to C++ code
        if ((kappa_hyb_h*hyb_L_<tau_u_ && kappa_hyb_h*hyb_L_>tau_l_) || (kappa_hyb_c*hyb_L_<tau_u_ && kappa_hyb_c*hyb_L_>tau_l_) || useHybridSolver_==false)
        {
            // Set rays convergence false
            List<bool> rayIdConv(nRay_, false);

            scalar maxResidual = 0;
            label radIter = 0;
            //absorptionEmission_->bandSet(j);
            fsck->bandSet(j);

            Info<<"Solving using fvDOM for quadrature point "<<j<<endl;
            do
            {
                Info<< "Radiation solver iter: " << radIter << endl;

                radIter++;
                maxResidual = 0;
                forAll(IRay_, rayI)
                {
                    if (!rayIdConv[rayI])
                    {
                        scalar maxBandResidual = IRay_[rayI].correct_FSCK(a_,aScale);
                        maxResidual = max(maxBandResidual, maxResidual);

                        if (maxBandResidual < tolerance_)
                        {
                            rayIdConv[rayI] = true;
                        }
                    }
                }

            } while (maxResidual > tolerance_ && radIter < maxIter_);

            updateG();

            //delqr_= delqr_ + weight[j]*dimensionedScalar(dimless/dimLength,1e2*kfsck[j])*(4.0*physicoChemical::sigma*pow4(T_) - G_);
            delqr_= delqr_ + weight[j]*a_*(4.0*aScale*physicoChemical::sigma*pow4(T_) - G_);

            // Compute the Rp and Ru terms here
            Rp_ = Rp_ + weight[j]*a_*aScale*4.0*physicoChemical::sigma;
            Ru_ = Ru_ + weight[j]*a_*G_;

            // Update the fvDOM iteration count
            fvDomIter++;
        }
        else
        {
            //a_ = dimensionedScalar(dimless/dimLength,1e2*kfsck[j]);
            //e_ = a_;
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
                1.0/(3.0*a_)
            );

            //absorptionEmission_->bandSet(j);
            fsck->bandSet(j);

            Info<<"Solving using P1 for quadrature point "<<j<<endl;
            // Solve G transport equation
            solve
            (
                fvm::laplacian(gamma, G_)
            - fvm::Sp(a_, G_)
            ==
            - 4.0*(a_*aScale*physicoChemical::sigma*pow4(T_))
            );

            //delqr_= delqr_ + weight[j]*a_*(4.0*physicoChemical::sigma*pow4(T_) - G_);
            delqr_= delqr_ + weight[j]*a_*(4.0*aScale*physicoChemical::sigma*pow4(T_) - G_);

            // Compute the Rp and Ru terms here
            Rp_ = Rp_ + weight[j]*a_*aScale*4.0*physicoChemical::sigma;
            Ru_ = Ru_ + weight[j]*a_*G_;

            // Update the P1 iteration count
            p1Iter++;
        }
        // Sync the solution over quadrature points over all processors
        label sync = j;
        reduce(sync, sumOp<label>());
    }
    // Write out the number of iterations for each solver
    Info << "Number of iterations for fvDOM solver: " << fvDomIter << endl;
    Info << "Number of iterations for P1 solver: " << p1Iter << endl;
}


Foam::tmp<Foam::volScalarField> Foam::radiation::fvDOM_FSCK::Rp() const
{
    // Construct using contribution from first frequency band
    /*tmp<volScalarField> tRp
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
            (
                4
               *physicoChemical::sigma
               *(aLambda_[0] - absorptionEmission_->aDisp(0)())
               *blackBody_.deltaLambdaT(T_, absorptionEmission_->bands(0))
            )
        )
    );

    volScalarField& Rp=tRp.ref();

    // Add contributions over remaining frequency bands
    for (label j=1; j < nLambda_; j++)
    {
        Rp +=
        (
            4
           *physicoChemical::sigma
           *(aLambda_[j] - absorptionEmission_->aDisp(j)())
           *blackBody_.deltaLambdaT(T_, absorptionEmission_->bands(j))
        );
    }*/

    //return tRp;
    return Rp_;
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::radiation::fvDOM_FSCK::Ru() const
{
    /*tmp<DimensionedField<scalar, volMesh>> tRu
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "Ru",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar(dimensionSet(1, -1, -3, 0, 0), Zero)
        )
    );

    DimensionedField<scalar, volMesh>& Ru=tRu.ref();

    // Sum contributions over all frequency bands
    for (label j=0; j < nLambda_; j++)
    {
        // Compute total incident radiation within frequency band
        tmp<DimensionedField<scalar, volMesh>> Gj
        (
            IRay_[0].ILambda(j)()*IRay_[0].omega()
        );

        for (label rayI=1; rayI < nRay_; rayI++)
        {
            Gj.ref() += IRay_[rayI].ILambda(j)()*IRay_[rayI].omega();
        }

        Ru += (aLambda_[j]() - absorptionEmission_->aDisp(j)()())*Gj
             - absorptionEmission_->ECont(j)()();
    }

    return tRu;*/
    const volScalarField::Internal& Ru_int = Ru_;
    return Ru_int;
}


void Foam::radiation::fvDOM_FSCK::updateBlackBodyEmission()
{
    for (label j=0; j < nLambda_; j++)
    {
        blackBody_.correct(j, absorptionEmission_->bands(j));
    }
}


void Foam::radiation::fvDOM_FSCK::updateG()
{
    G_ = dimensionedScalar(dimMass/pow3(dimTime), Zero);
    qr_ = dimensionedScalar(dimMass/pow3(dimTime), Zero);
    qem_ = dimensionedScalar(dimMass/pow3(dimTime), Zero);
    qin_ = dimensionedScalar(dimMass/pow3(dimTime), Zero);

    forAll(IRay_, rayI)
    {
        IRay_[rayI].addIntensity();
        G_ += IRay_[rayI].I()*IRay_[rayI].omega();
        qr_.boundaryFieldRef() += IRay_[rayI].qr().boundaryField();
        qem_.boundaryFieldRef() += IRay_[rayI].qem().boundaryField();
        qin_.boundaryFieldRef() += IRay_[rayI].qin().boundaryField();
    }
}


void Foam::radiation::fvDOM_FSCK::setRayIdLambdaId
(
    const word& name,
    label& rayId,
    label& lambdaId
) const
{
    // Assuming name is in the form: CHARS_rayId_lambdaId
    const auto i1 = name.find('_');
    const auto i2 = name.find('_', i1+1);

    rayId    = readLabel(name.substr(i1+1, i2-i1-1));
    lambdaId = readLabel(name.substr(i2+1));
}


const Foam::solarCalculator& Foam::radiation::fvDOM_FSCK::solarCalc() const
{
    return solarCalculator_();
}


// ************************************************************************* //
