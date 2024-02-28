/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2016 OpenFOAM Foundation
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

Class
    Foam::LESModels::WALEmod2

Group
    grpLESTurbulence

Description
    The Wall-adapting local eddy-viscosity (WALEmod2) SGS model.

    Reference:
    \verbatim
        Nicoud, F., & Ducros, F. (1999).
        Subgrid-scale stress modelling based on the square of the velocity
        gradient tensor.
        Flow, Turbulence and Combustion, 62(3), 183-200.
    \endverbatim

    The default model coefficients are
    \verbatim
        WALEmod2Coeffs
        {
            Ck                  0.094;
            Ce                  1.048;e
            Cw                  0.325;
        }
    \endverbatim

See also
    Foam::LESModels::Smagorinsky

SourceFiles
    WALEmod2.C

\*---------------------------------------------------------------------------*/

#ifndef WALEmod2_H
#define WALEmod2_H

#include "LESModel.H"
#include "LESeddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                           Class WALEmod2 Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class WALEmod2
:
    public LESeddyViscosity<BasicTurbulenceModel>
{
    // Private Member Functions

        //- No copy construct
        WALEmod2(const WALEmod2&) = delete;

        //- No copy assignment
        void operator=(const WALEmod2&) = delete;


protected:

    // Protected data

        dimensionedScalar Ck_;
        dimensionedScalar Cw_;
    dimensionedScalar startAveraging_;

    volScalarField kRES_;
    volScalarField kSGS_;
    volScalarField kRatio_;
    volScalarField epsilonSGS_;
    volScalarField epsilonRES_;
    volScalarField deltaOverEta_;
    volVectorField Uavg_;

    // Protected Member Functions

        //- Return the deviatoric symmetric part of the square of the given
        //  velocity gradient field
        tmp<volSymmTensorField> Sd(const volTensorField& gradU) const;

        //- Return SGS kinetic energy
        //  calculated from the given velocity gradient
        tmp<volScalarField> k(const volTensorField& gradU) const;

        //- Update the SGS eddy-viscosity
        virtual void correctNut();


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("WALEmod2");


    // Constructors

        //- Construct from components
        WALEmod2
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~WALEmod2() = default;


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read();

        //- Return SGS kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k(fvc::grad(this->U_));
        }

        //- Return sub-grid dissipation rate
        virtual tmp<volScalarField> epsilon() const;

        //- Correct Eddy-Viscosity and related properties
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "WALEmod2.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
