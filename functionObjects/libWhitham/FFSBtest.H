/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

Class
    Foam::FFPtests

Description

    Module for testing Whitham library.

    Input data:

    Output:
    -   file with results of testing Whitham functions.

SourceFiles
    FFPtests.C

\*---------------------------------------------------------------------------*/

#ifndef FFSBtest_H
#define FFSBtest_H

#include "coordinateSystem.H"
#include "coordinateSystems.H"
#include "primitiveFieldsFwd.H"
#include "volFieldsFwd.H"
#include "HashSet.H"
#include "Tuple2.H"
#include "OFstream.H"
#include "Switch.H"

#include "probes.H"
#include "pointFieldFwd.H"

#include "polyMesh.H"
#include "mapPolyMesh.H"

#include "volFields.H"
#include "passiveParticleCloud.H"

#include "functionObject.H"
#include "sampledSet.H"

#include "wordReList.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace Foam
{
//Forward declaration of classes

class fvMesh;
class dictionary;

/*---------------------------------------------------------------------------*\
                     Class FFPtests Declaration
\*---------------------------------------------------------------------------*/

class FFPtests
:
    public functionObject
{

private:

    public:

//

protected:

    // Private data

private:

    //- Disallow default bitwise copy construct
    FFPtests(const FFPtests&);

    //- Disallow default bitwise assignment
    void operator=(const FFPtests&);

    //- Write geometry only
    void writeGeometry() const;

    //- If the soundPressureSampler file has not been created create it
//    void makeFile();

public:

    //- Runtime type information
    TypeName("FFPtests");

    // Constructors

    //- Construct from Time and dictionary
    FFPtests
    (
        const word& name,
        const Time& time,
        const dictionary& dict
    );

    //- Destructor

    virtual ~FFPtests();

    // Member Functions


    //- Read the farFieldPressure data
    virtual bool read(const Foam::dictionary&);

    //- Execute, currently does nothing
    virtual bool execute();

    //- Execute at the final time-loop, currently does nothing
//    virtual bool end();

    // Write farFieldPessure
    virtual bool write();



    //Place for test function



    //testing constructor for Whitham functionObject
    bool Foam::FFPtest::test_constructor();

    // testing the module for reading the dictionry for Whitham functionObject
    bool Foam::FFPtest::test_read(const dictionary& dict);

    // testing the module for checking initial values
    bool Foam::FFPtest::test_initialValues(const dictionary& dict);

    // testing the module for calculating velocity magnitude in the distance from the body
    bool Foam::FFPtest::test_infVelMag();

    // testing the module for calculating beta coefficent
    bool Foam::FFPtest::test_beta(const scalar, const scalar, const scalar, const scalar );

    // testing the module for calculating equivalent body
    bool Foam::FFPtest::test_dSeqBydX(const vector , const scalar , const vector , const vector, label);

    // testing the module for calculating equivalent body in new coordinate system
    bool Foam::FFPtest::test_new_dSeqBydX(const scalar, const scalar, const scalar, const scalar, const scalar, const scalar);

    // testing the module for calculating part of Whitham function
    bool Foam::FFPtest::test_PhiOfEta( point, point, label, List<scalar>);

    // testing the module for calculating Whitham function
    bool Foam::FFPtest::test_FofEta(label, List<scalar>, point, point);

    // testing the function for calculating attenuation coefficient k for pressure in the far field
    bool Foam::FFPtest::test_kCoef();

    // testing the function for calculating attenuation coefficient k1 for pressure in the far field
    bool Foam::FFPtest::test_k1Coef();

    // testing the function for calculating attenuation coefficient k2 for pressure in the far field
    bool Foam::FFPtest::test_k2Coef();

    // testing the function for calculating pressure in the far field
    bool Foam::FFPtest::test_deltaP();

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
