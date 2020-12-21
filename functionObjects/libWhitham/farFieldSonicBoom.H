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
    Foam::farFieldSonicBoom

Description

    Return sound amplitudes of sound pressure in each point of STL-surface.
    Data are sorted by frequencies.

    Input data:

    Output:
    -   file with pressure history.

SourceFiles
    farFieldSonicBoom.C

\*---------------------------------------------------------------------------*/

#ifndef farFieldSonicBoom_H
#define farFieldSonicBoom_H

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


#include "forces.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// Forward declaration of classes

class fvMesh;
class dictionary;

/*---------------------------------------------------------------------------*\
                    Class farFieldPressure Declaration
\*---------------------------------------------------------------------------*/

class farFieldPressure
:
    public functionObject
{

   //- Center of surface for sampling
//   vector c_;
   
   //- Force acting on patch
//   tdPrimitive<vector> F_;

private:

    //- Class used for sampling volFields
    template<class Type>
    class volFieldSampler
    :
        public List<Field<Type>>
    {
        //- Name of this collection of values
        const word name_;

    public:

        //- Construct interpolating field to the sampleSets
        volFieldSampler
        (
            const word& interpolationScheme,
            const GeometricField<Type, fvPatchField, volMesh>& field,
            const PtrList<sampledSet>&
        );

        //- Return the field name
        const word& name() const
        {
            return name_;
        }
    };

    //- Tolerance for merging points (fraction of mesh bounding box)
    static scalar mergeTol_;

protected:

    // Private data

    //- Const reference to fvMesh
    const fvMesh& mesh_;

    //
//    const bool loadFromFiles;

    //- Interpolation scheme to use for sampling
    word interpolationScheme_;

    // Also used as the name of the probes directory.
    word name_;

    //- glogal string to keep directory name for any output
    std::string dirName_;

    //- Standard OpenFOAM object registry
//    const objectRegistry& obr_;

    //- Switch to send output to Info as well as to file
//    Switch log_;

    //
    label nLines;

    //
    label nPoints;

    //
    point xs;

    //
    point xe;

    //
    point xn;

    //
    scalar yn;

    //
    scalar H;

    //
    vector infVel;

    //
    vector TotalForce;

    //
    scalar Wgamma;

    //
    scalar WR;

    //
    scalar infT;

    //
    scalar length;

    //
    scalar AoA;
    //- Name of pressure field
//    word fieldName_;

    //- Velocity magnitude value at infinity
    scalar infVelMagVal;

    //
    scalar betaVal;

    //
    scalar M_ef;

    //
    scalar beta_ef;

    //
    scalar rho_00;

    //
    scalar rho_0;

    //
    scalar rho_inf;

    //
    scalar gamma_bar;



private:

    //- Disallow default bitwise copy construct
    farFieldPressure(const farFieldPressure&);

    //- Disallow default bitwise assignment
    void operator=(const farFieldPressure&);

    //- Write geometry only
    void writeGeometry() const;

    //- If the farFieldPressure file has not been created create it
    void makeFile(List<scalar> res_dp);

public:

    //- Runtime type information
    TypeName("farFieldPressure");

    // Constructors

    //- Construct from Time and dictionary
    farFieldPressure
    (
        const word& name,
        const Time& time,
        const dictionary& dict
    );

    //- Destructor

    virtual ~farFieldPressure();

    // Member Functions

    //
    virtual const word& name() const
    {
        return name_;
    }

    //- Read the farFieldPressure data
    virtual bool read(const Foam::dictionary&);

    //- Execute, currently does nothing
    virtual bool execute();

    //- Execute at the final time-loop, currently does nothing
    //    virtual bool end();

    // Write farFieldPessure
    virtual bool write();

    //
    scalar infVelMag ();

    //
    scalar beta( );

    //
    scalar dSeqBydX(const vector , const scalar , const vector , const vector, label);

    //
    scalar new_dSeqBydX( const scalar);

    //
    List<Foam::scalar> PhiOfEta( point, point, label, List<scalar>);

    //
    List<Foam::scalar> FofEta(label, List<scalar>, point, point);

    //
    Foam::scalar getL( scalar , scalar , scalar , scalar );

    //
    Foam::scalar getYforce();

    //
    Foam::scalar k1Coef( scalar);

    //
    List<Foam::scalar> deltaP( label , scalar , List<scalar> );

/*
    //
    scalar kCoef();

    //
    scalar k2Coef();
*/

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#ifdef NoRepository
    #include "farFieldSonicBoomTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
