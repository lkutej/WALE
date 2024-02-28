/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "WALEmod1.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volSymmTensorField> WALEmod1<BasicTurbulenceModel>::Sd
(
    const volTensorField& gradU
) const
{
    return dev(symm(gradU & gradU));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> WALEmod1<BasicTurbulenceModel>::k
(
    const volTensorField& gradU
) const
{
    volScalarField magSqrSd(magSqr(Sd(gradU)));

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("k", this->alphaRhoPhi_.group()),
                this->runTime_.timeName(),
                this->mesh_
            ),
            sqr(sqr(Cw_)*this->delta()/Ck_)*
            (
                pow3(magSqrSd)
               /(
                   sqr
                   (
                       pow(magSqr(symm(gradU)), 5.0/2.0)
                     + pow(magSqrSd, 5.0/4.0)
                   )
                 + dimensionedScalar
                   (
                       "SMALL",
                       dimensionSet(0, 0, -10, 0, 0),
                       SMALL
                   )
               )
            )
        )
    );
}


template<class BasicTurbulenceModel>
void WALEmod1<BasicTurbulenceModel>::correctNut()
{
    this->nut_ = Ck_*this->delta()*sqrt(this->k(fvc::grad(this->U_)));
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
WALEmod1<BasicTurbulenceModel>::WALEmod1
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    LESeddyViscosity<BasicTurbulenceModel>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    Ck_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ck",
            this->coeffDict_,
            0.094
        )
    ),

    Cw_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cw",
            this->coeffDict_,
            0.325
        )
    ),
/*
    startAveraging_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "startAveraging",
            //&const_cast<dictionary>(this->runTime_.controlDict()),  
            this->coeffDict_,
            0
        )
    ),
*/
/*
    startAveraging_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "startAveraging",
            (runTime.controlDict()),  
            //this->coeffDict_,
            0
        )
    ),
*/
  startAveraging_( readScalar(this->runTime_.controlDict().lookup("startAveraging"))),
//    startAveraging_ (this->runTime_.controlDict().lookupOrDefault<scalar>("startAveraging", 0)),


    kRES_
    (
        IOobject
        (
            IOobject::groupName("kRES", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("kRES", dimVelocity*dimVelocity, 0)
    ),

    kSGS_
    (
        IOobject
        (
            IOobject::groupName("kSGS", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("kSGS", dimVelocity*dimVelocity, 0)
    ),

    kRatio_
    (
        IOobject
        (
            IOobject::groupName("kRatio", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("kRatio", dimless, 0)
    ),

    epsilonSGS_
    (
        IOobject
        (
            IOobject::groupName("epsilonSGS", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("epsilonSGS", dimVelocity*dimVelocity/dimTime, 0)
    ),

    epsilonRES_
    (
        IOobject
        (
            IOobject::groupName("epsilonRES", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("epsilonRES", dimVelocity*dimVelocity/dimTime, 0)
    ),

    deltaOverEta_
    (
        IOobject
        (
            IOobject::groupName("deltaOverEta", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("deltaOverEta", dimless, 0)
    ),

    Uavg_
    (
        IOobject
        (
            IOobject::groupName("Uavg", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedVector("Uavg", dimVelocity, Zero)
    )

{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool WALEmod1<BasicTurbulenceModel>::read()
{
    if (LESeddyViscosity<BasicTurbulenceModel>::read())
    {
        Ck_.readIfPresent(this->coeffDict());
        Cw_.readIfPresent(this->coeffDict());
        startAveraging_.readIfPresent(this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
tmp<volScalarField> WALEmod1<BasicTurbulenceModel>::epsilon() const
{
    volScalarField k(this->k(fvc::grad(this->U_)));

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("epsilon", this->alphaRhoPhi_.group()),
                this->runTime_.timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->Ce_*k*sqrt(k)/this->delta()
        )
    );
}


template<class BasicTurbulenceModel>
void WALEmod1<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    LESeddyViscosity<BasicTurbulenceModel>::correct();


    const volVectorField& U = this->U_;
    tmp<volTensorField> tgradU(fvc::grad(U));
    const volTensorField& gradU = tgradU();

    volSymmTensorField Sd= dev(symm(gradU & gradU));
    volScalarField magSqrSd = magSqr(Sd);

    kSGS_ = sqr(sqr(Cw_)*this->delta()/Ck_)*(pow3(magSqrSd)/(sqr(pow(magSqr(symm(gradU)), 5.0/2.0)+ pow(magSqrSd, 5.0/4.0))+ dimensionedScalar("SMALL",dimensionSet(0, 0, -10, 0, 0),SMALL)));

    epsilonSGS_ = this->Ce_*kSGS_*sqrt(kSGS_)/this->delta();

    // *******************************************************************************************************
    //
    //
    //
    // *******************************************************************************************************
    scalar time = this->runTime_.value();
    scalar dt = this->runTime_.deltaT().value();

    scalar startAvg = startAveraging_.value();

    if (time > startAvg )
    {
        if ( (time - dt) < startAvg )
        {
            Uavg_ = U;
        }
        else
        {
            Uavg_ = ( Uavg_*(time-startAvg) + U*dt ) / (time-startAvg + dt);
        }
    }

    // k resolved
    kRES_ =  0.5*magSqr(U-Uavg_);
    
    // kRatio
    kRatio_ = kSGS_/(kSGS_ + kRES_);

    // epsilon resolved
    tmp<volTensorField> tgradUp(fvc::grad(U-Uavg_));
    const volTensorField& gradUp = tgradUp();
    epsilonRES_ = this->nu()*(gradUp && gradUp);

    // deltaOverEta
//    deltaOverEta = 1;
//    deltaOverEta_ = pow(this->delta(),1.0/3.0)/pow((pow3(this->nu())/(epsilonSGS_+epsilonRES_)),0.25);
    deltaOverEta_ = this->delta()/pow((pow3(this->nu())/(epsilonSGS_+epsilonRES_)),0.25);

    //

    // *******************************************************************************************************
    //
    //
    //
    // *******************************************************************************************************

    //        Uavg_ = ( Uavg_*(time) + U*dt ) / (time + dt);
    



    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
