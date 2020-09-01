/*---------------------------------------------------------------------------*\
   =========                 |
   \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
    \\    /   O peration     |
     \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
      \\/     M anipulation  |
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

#include "farFieldSonicBoom.H"
#include "FFSBtest.H"

#include "volFields.H"
#include "dictionary.H"
#include "Time.H"
#include "argList.H"
#include "IOmanip.H"
#include "volPointInterpolation.H"
#include "PatchTools.H"
#include "ListListOps.H"
#include "stringListOps.H"
#include "addToRunTimeSelectionTable.H"
#include "uniformSet.H"
#include "meshSearch.H"
#include "pointFieldFwd.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(farFieldPressure, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        farFieldPressure,
        dictionary
    );
}

Foam::scalar Foam::farFieldPressure::mergeTol_ = 1e-10;

// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::FFPtests::FFPtests
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    mesh_
    (
        refCast<const fvMesh>
        (
            t.lookupObject<objectRegistry>
            (
                dict.lookupOrDefault("region", polyMesh::defaultRegion)
            )
        )
    ),

//    name_(name),
//    interpolationScheme_(word::null)

{
//    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::FFPtests::~FFPtests()
{}

// * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

//testing constructor for Whitham functionObject
bool Foam::FFPtest::test_constructor()
{

    return true;
}

// testing the module for reading the dictionry for Whitham functionObject
bool Foam::FFPtest::test_read(const dictionary& dict)
{

    return true;
}

// testing the module for checking initial values
bool Foam::FFPtest::test_initialValues(const dictionary& dict)
{

    return true;
}


// testing the module for calculating velocity magnitude in the distance from the body
bool Foam::FFPtest::test_infVelMag()
{

    return true;
}

// testing the module for calculating beta coefficent
bool Foam::FFPtest::test_beta(const scalar, const scalar, const scalar, const scalar )
{

    return true;
}

// testing the module for calculating equivalent body
bool Foam::FFPtest::test_dSeqBydX(const vector , const scalar , const vector , const vector, label)
{

    return true;
}

// testing the module for calculating equivalent body in new coordinate system
bool Foam::FFPtest::test_new_dSeqBydX(const scalar, const scalar, const scalar, const scalar, const scalar, const scalar)
{

    return true;
}

// testing the module for calculating part of Whitham function
bool Foam::FFPtest::test_PhiOfEta( point, point, label, List<scalar>)
{

    return true;
}

// testing the module for calculating Whitham function
bool Foam::FFPtest::test_FofEta(label, List<scalar>, point, point)
{

    return true;
}

// testing the function for calculating attenuation coefficient k for pressure in the far field
bool Foam::FFPtest::test_kCoef()
{

    return true;
}

// testing the function for calculating attenuation coefficient k1 for pressure in the far field
bool Foam::FFPtest::test_k1Coef()
{

    return true;
}

// testing the function for calculating attenuation coefficient k2 for pressure in the far field
bool Foam::FFPtest::test_k2Coef()
{

    return true;
}

// testing the function for calculating pressure in the far field
bool Foam::FFPtest::test_deltaP()
{

    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
