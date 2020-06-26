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
#include "volFields.H"
#include "dictionary.H"
#include "Time.H"
#include "argList.H"

//sampledSurfaces stuff
//#include "meshedSurfRef.H"
#include "IOmanip.H"
//#include "volPointInterpolation.H"
#include "PatchTools.H"
#include "ListListOps.H"
//#include "stringListOps.H"
#include "surfaceFields.H"
//#include "addToRunTimeSelectionTable.H"

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

Foam::word Foam::farFieldPressure::fieldName(const word& name) const
{
    return this->name() + ":" + name;
}


void Foam::farFieldPressure::createFiles()


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::farFieldPressure::farFieldPressure
(
    const word& name,
    const Time& t,
    const volVectorField& U,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    log_(false),
    outputFormat_(word::null),
    fieldName_(word::null),
    pointA_(vector::zero),
    pointB_ (vector::zero),
    nose_ (vector::zero),
    infVel_ (vector::zero),
    infVelMag (scalar 0),
    Wgamma_ (scalar 0),
    WR_ (scalar 0),
    length_ (scalar 0),
    dSeqByX_ (scalar 0),
    new_dSeqByX_ (scalar 0)

{
    Info << "in constructor" << nl;

    read(dict);
}

Foam::farFieldPressure::farFieldPressure
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    log_(false),
    outputFormat_(word::null),
    fieldName_(word::null),
    pointA_(vector::zero),
    pointB_ (vector::zero),
    nose_ (vector::zero),
    infVel_ (vector::zero),
    infVelMag_ (scalar 0),
    Wgamma_ (scalar 0),
    WR_ (scalar 0),
    length_ (scalar 0),
    dSeqByX_ (scalar 0),
    new_dSeqByX_ (scalar 0)
{
    read(dict);
}

// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::farFieldPressure::~farFieldPressure()
{}

// * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

bool Foam::soundPressureSampler::read(const dictionary& dict)
{
    //where should be log info 
    log_ = dict.lookupOrDefault<Switch>("log", false);

    if (!log_)
    {
        Info << "Direct logging to stdio disabled" << endl
        << " to enable, please insert string:" << endl
        << "log\t\t true;" << endl
        << "in dictionary" << endl;
    }

    dict.lookup("pointA") >> pointA_;
    dict.lookup("pointB") >> pointB_;
    dict.lookup("nose") >> nose_;
    dict.lookup("infVel") >> infVel_;
    dict.lookup("Wgamma") >> Wgamma_;
    dict.lookup("WR") >> WR_;
    dict.lookup("length") >> length_;
    //read mesh
    //const fvMesh& mesh_ = refCast<const fvMesh>(obr_);

    // Ensure all surfaces and merge information are expired
    expire();

    dict.lookup("fieldName") >> fieldName_;

    return true;
}

void Foam::soundPressureSampler::makeFile()
{
    fileName soundPressureSamplerDir;

    if (Pstream::master() && Pstream::parRun())
    {
        //farFieldSonicBoomDir = obr_.time().rootPath() + "/" + obr_.time().caseName().path() + "/data" ;
        farFieldPressureDir = mesh_.time().rootPath() + "/" + mesh_.time().caseName().path() + "/data" ;

        mkDir(farFieldPressure);
    }
    else if (!Pstream::parRun())
    {
        //farFieldSonicBoomDir = obr_.time().rootPath() + "/" + obr_.time().caseName() + "/data" ;
        farFieldPressureDir = mesh_.time().rootPath() + "/" + mesh_.time().caseName() + "/data" ;

        mkDir(farFieldPressureDir);
    }
}
Foam::scalar infVelMag( infVel_, const scalar Wgamma_, const scalar WR_, const scalar infT_) const
(
    scalar infVelMag_;
)
{
    infVelMag_ = sqrt(infVel_&infVel_);
    return infVelMag_;
}

Foam::scalar Foam::farFieldPressure::beta(const scalar infVelMag_, const scalar Wgamma_, const scalar WR_, const scalar infT_) const
{
    beta_ = infVelMag_*infVelMag_/(Wgamma_*WR_*infT_);
    return beta_;
}

Foam::scalar Foam::farFieldPressure::dSeqByX(const vector infVel_, const scalar beta_, const vector new_pointA_, const vector new_pointB_) const
{
    vector BA = ((new_pointA[0] - new_pointB[0]) (new_pointA[1]-new_pointB[1]) (new_pointA_[2]-new_pointB_[2]));
    scalar magBA = mag(BA);
    const pointField & pp = mesh.points();

    forAll ( mesh.C().internalField(), celli)
    {
        const cell & cc = mesh.cells()[celli];
        labelList pLabels(cc.labels(ff));
        pointField pLocal(pLabels.size(), vector::zero);
        vector vz= U.component(2);
        vector BC = ((mech.C()[0]-new_pointB[0])  (mech.C()[2]-new_pointB[2]) 0);
        scalar distC = mag(BC&BA/magBA);
        scalar sumVz = 0;
        forAll (pLabels, pointi)
            pLobal[pointi] = pp[pLabels[pointi]];

        if(Foam::min(pLocal & vector(0,1,0)) < (new_pointB_[1]+new_pointA_[1])/2 && (new_pointB_[1]+new_pointA_[1])/2 < Foam::max(pLocal & vector(0,1,0)))
        {
            scalar ZdimC = mag(Foam::max(pLocal & vector(0,0,1))-Foam::min(pLocal & vector(0,0,1)));
            scalar XdimC = mag(Foam::max(pLocal & vector(1,0,0))-Foam::min(pLocal & vector(1,0,0)));
            scalar XZdimC = sqrt(ZdimC*ZdimC+XdimC*XdimC)/2;

            if(XZdimC<distC)
            {
               sumVz = sumVz + vx*ZdimC; 
            }
        }
    }
    scalar dSeqByX_= -2/infVelMag_sumVz;
    return dSeqByX_;
}

Foam::scalar  Foam::farFieldPressure::new_dSeqByX(const scalar length_, const scalar WR_, const scalar infT_, const scalar beta_, const scalar Wgamma_, const scalar dSeqByX_) const
{
     new_dSeqByX_ = Wgamma_*WR_*Wgamma_*infT_/(2*length_*constant::mathematical::pi*sqrt(2*beta_))*dSeqByX_;
     return new_dSeqByX_;
}

Foam::volScalarField Foam::farFieldPressure::PhiOfEta() const
{
    scalar PhiOfEta_=0;
    return PhiOfEta_;
}

Foam::scalar Foam::farFieldPressure::FoEta() const
{
    scalar FoEta_=0;
    return FoEta_;
}

Foam::scalar Foam::farFieldPressure::kCoef() const
{
    scalar kCoef_=0;
    return kCoef_;
}

Foam::scalar Foam::farFieldPressure::k1Coef() const
{
    scalar k1Coef_=0;
    return k1Coef_;
}

Foam::scalar Foam::farFieldPressure::deltaP() const
{
    scalar deltaP_ = 0;
    return deltaP_;
}
// ************************************************************************* //
