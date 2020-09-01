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

#include "IOmanip.H"
#include "volPointInterpolation.H"
#include "PatchTools.H"
#include "ListListOps.H"
#include "stringListOps.H"

#include "addToRunTimeSelectionTable.H"

#include "uniformSet.H"

#include "meshSearch.H"
#include "pointFieldFwd.H"

#include "forces.H"


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

/*
Foam::word Foam::farFieldPressure::fieldName(const word& name)
{
    return this->name() + ":" + name;
}
*/

//void Foam::farFieldPressure::createFiles()


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::farFieldPressure::farFieldPressure
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
//    loadFromFiles = false,

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

    interpolationScheme_(word::null)
    
{
    read(dict);
}

// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::farFieldPressure::~farFieldPressure()
{}

// * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

bool Foam::farFieldPressure::read(const dictionary& dict)
{
//    name_(name);
//    log_(false);

    dict.lookup("interpolationScheme") >> interpolationScheme_;

    infVel_ = ( dict.get<vector>("infVel") );
    infT_   = ( dict.get<scalar>("infT")   );
    Wgamma_ = ( dict.get<scalar>("Wgamma") );
    WR_     = ( dict.get<scalar>("WR")     );
    length_ = ( dict.get<scalar>("length") );

    scalar rho_00    = ( dict.get<scalar>( "rho_00")   );
    scalar rho_0     = ( dict.get<scalar>( "rho_0" )   );
    scalar rho_inf   = ( dict.get<scalar>( "rho_00")   );
    scalar gamma_bar = ( dict.get<scalar>("gamma_bar") );


    point newPointA;
    point newPointB;

    uniformSet baseLine ("baseLine", mesh_, meshSearch(mesh_), dict.subDict("baseLine"));

    Info << "dict:     "<< dict << nl;
    Info << "baseLine  "<< dict.subDict("baseLine") << nl;

    nLines= ( readLabel( dict.lookup("nLines")) );
    nPoints=( readLabel( dict.subDict("baseLine").lookup("nPoints")) );


    xs=( dict.subDict("baseLine").lookup("start") );
    xe=( dict.subDict("baseLine").lookup("end")   );
    xn=( dict.lookup("nose")  );

    point xm( (xe + xs)*0.5 );
    point xmn( xm-xn );
    xmn[1]= 0.0;

    scalar infVelMag_( farFieldPressure::infVelMag( infVel_) );
    point newXmn( farFieldPressure::beta(infVelMag_,Wgamma_,WR_,infT_)*abs(xm[1]-xn[1])/sqrt(mag(xm-xn)*mag(xm-xn)-(xm[1]-xn[1])*(xm[1]-xn[1]))*xmn );
    newXmn[1] = xm[1]-xn[1];
    point newXn (xn+newXmn);


    Info << "nose:      " << xn                   << nl;
    Info << "newNose:   " << newXn                << nl;
    Info << "xm:        " << xm                   << nl;

    Info << "infVelMag: " << infVelMag_                                           << nl;
    Info << "beta:      " << farFieldPressure::beta(infVelMag_,Wgamma_,WR_,infT_) << nl;
    Info << "           " << nl;

    point dy( (newXn - xm) / static_cast<scalar>(nLines) );

    PtrList<sampledSet>  probeSets_(nLines);

    for (int i=0; i < nLines; i++)
    {
        newPointA = xs + i*dy;
        newPointB = xe + i*dy;
        word iLine = "";
        char numLine [6];
        sprintf (numLine, "%d", i+1);
        word inLine (numLine);
        iLine = "Line" + inLine;

        probeSets_.set (i, new uniformSet(iLine, mesh_, meshSearch(mesh_), "xyz", newPointA, newPointB, nPoints));

        Info << "line:         " << iLine                << nl;
        Info << "              " << probeSets_[i]        << nl;
        Info << "start point:  " << newPointA            << nl;
        Info << "end point:    " << newPointB            << nl;
        Info << "line size:    " << probeSets_[i].size() << nl;

        if( probeSets_[i].size()< nPoints)
        {
            FatalErrorInFunction
            << "Error bad line set parameters, lines should be inside of domain"
            << nl << exit(FatalError);
        }
    }

    dict.lookup("fieldName") >> fieldName_;

    volScalarField field_ (mesh_.thisDb().lookupObject<volVectorField>("U").component(1));

    volFieldSampler<Foam::scalar>  volFieldSampler_(interpolationScheme_, field_, probeSets_);

    Info << "volFieldSampler:          " << volFieldSampler_       << nl;
    Info << "volFieldSampler[1][2]:    " << volFieldSampler_[1][2] << nl;

/*
    //where should be log info 
    log_ = dict.lookupOrDefault<Switch>("log", false);

    if (!log_)
    {
        Info << "Direct logging to stdio disabled" << endl
        << " to enable, please insert string:" << endl
        << "log\t\t true;" << endl
        << "in dictionary" << endl;
    }

        Info<< "Function object "<< name_<<":" << nl;

        Info<< endl;

        writeGeometry();
        Info << "writeGeometry OK" << nl;
    }
*/

    //Setting default adress for files to save.
    if (Pstream::parRun() )
    {
        dirName_ = mesh_.time().rootPath() + "/" + mesh_.time().caseName().path() + "/libWhitham/";
    }
    else
    {
        dirName_ = mesh_.time().rootPath() + "/" + mesh_.time().caseName() + "/libWhitham/";
    }
    return true;
}

bool Foam::farFieldPressure::write()
{
//    update();
//    dict.lookup("fieldName") >> fieldName_;

    point newPointA;
    point newPointB;

    point xm( (xe + xs)*0.5 );

    Info << "write xm:  " << xm << nl;

    point xmn( xm-xn );
    xmn[1]= 0.0;

    scalar infVelMag_ = farFieldPressure::infVelMag( infVel_);
    point newXmn( farFieldPressure::beta(infVelMag_,Wgamma_,WR_,infT_)*abs(xm[1]-xn[1])/sqrt(mag(xm-xn)*mag(xm-xn)-(xm[1]-xn[1])*(xm[1]-xn[1]))*xmn );
    newXmn[1] = xm[1]-xn[1];
    point newXn (xn+newXmn);

    point dy( (newXn - xm) / static_cast<scalar>(nLines) );

    PtrList<sampledSet>  probeSets_(nLines);

    for (int i=0; i < nLines; i++)
    {
        newPointA = xs + i*dy;
        newPointB = xe + i*dy;
        word iLine = "";
        char numLine [6];
        sprintf (numLine, "%d", i+1);
        word inLine (numLine);
        iLine = "Line" + inLine;

        probeSets_.set (i, new uniformSet(iLine, mesh_, meshSearch(mesh_), "xyz", newPointA, newPointB, nPoints));

        Info << "line:         " << iLine                << nl;
        Info << "              " << probeSets_[i]        << nl;
        Info << "start point:  " << newPointA            << nl;
        Info << "end point:    " << newPointB            << nl;
        Info << "line size:    " << probeSets_[i].size() << nl;

        if( probeSets_[i].size()< nPoints)
        {
            FatalErrorInFunction
            << "Error bad line set parameters, lines should be inside of domain"
            << nl << exit(FatalError);
        }
    }

    volScalarField field_ (mesh_.thisDb().lookupObject<volVectorField>("U").component(0));

    volFieldSampler<Foam::scalar>  volFieldSampler_(interpolationScheme_, field_, probeSets_);

    Info << "volFieldSampler:          " << volFieldSampler_       << nl;
    Info << "volFieldSampler[1][2]:    " << volFieldSampler_[1][2] << nl;

//    scalar infVelMag_ = infVelMag( infVel_);


    //
    Foam::List<Foam::scalar> dSeqBydX_(nLines);
    for (int i=0; i < nLines; i++)
    {
        dSeqBydX_[i]= 0.0;
        Info << "volFieldSampler[i]:    " << volFieldSampler_[i] << nl;

        for (int j=0; j < volFieldSampler_[i].size() ; j++)
        {
            dSeqBydX_[i] += volFieldSampler_[i][j]* mag(newPointB[2]-newPointA[2]) / volFieldSampler_[i].size();
        }
        dSeqBydX_[i]= -2/infVelMag_*dSeqBydX_[i];
    }


    //
    beta_ = beta( infVelMag_, Wgamma_, WR_, infT_ );


    //
    Foam::List<Foam::scalar> new_dSeqBydX_(nLines);
    for (int i=0; i < nLines; i++)
    {
        new_dSeqBydX_[i] = 0.0;
        new_dSeqBydX_[i] = new_dSeqBydX(length_, WR_, infT_, beta_, infVelMag_, dSeqBydX_[i]);
    }

    Info << "new_dSeqBydX_1:    " << new_dSeqBydX_[1] << nl;


    Foam::List<Foam::scalar> PhiOfEta_(nLines);
    PhiOfEta_ = PhiOfEta( xn, xm, nLines, new_dSeqBydX_);

    //
    List<scalar> FofEta_(nLines-1);
    FofEta_ = FofEta(nLines, PhiOfEta_, xn, xm);


    return true;
}

bool Foam::farFieldPressure::execute()
{
    return true;
}


/*
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
*/

/*
void Foam::farFieldPressure::sampleAndCalc
(
    const GeometricField<scalar, fvPatchField, volMesh>& vField,
    const label surfI
)
{
    Field<scalar> values;
    const sampledSurface& surface = controlSurfaces_.operator[](surfI);

    getValues(vField, surface, values);

    if (!values.empty())
    {
//        gainCurrentValues(vField, values);
    }

//    refreshPreviousValues(values, vField.time().timeOutputValue());

//        makeSnapshot(vField, values);

//            prepareAndGatherData();

    if ( Pstream::master() )
    {
        calcWhithamFunc(); //this function also writes results

        // Saving coordinates
//        writeCoordinates(dirName_, surface, surfI);

        Info << nl << "All files in field "<< fieldName_ <<" were saved."<< nl
        <<"Whitham function for field "<< fieldName_<<" is finished"<< nl;
    }
}
*/

//
Foam::scalar Foam::farFieldPressure::infVelMag(vector infVel_)
{
    scalar infVelMag_ = 0;
    infVelMag_=sqrt(infVel_&infVel_);

    return infVelMag_;
}

//
Foam::scalar Foam::farFieldPressure::beta(const scalar infVelMag_, const scalar Wgamma_, const scalar WR_, const scalar infT_)
{
    scalar beta_ = infVelMag_*infVelMag_/(Wgamma_*WR_*infT_);
    return beta_;
}

//
Foam::scalar Foam::farFieldPressure::dSeqBydX(const vector infVel_, const scalar beta_, const vector new_pointA_, const vector new_pointB_, const label nPoints)
{
    scalar dSeqByX_ =0;

    return dSeqByX_;
}

//
Foam::scalar  Foam::farFieldPressure::new_dSeqBydX(const scalar length_, const scalar WR_, const scalar infT_, const scalar beta_, const scalar infVelMag_, const scalar dSeqByX_)
{
     scalar new_dSeqByX_ = infVelMag_*infVelMag_/WR_/infT_/(2*length_*constant::mathematical::pi*sqrt(2*beta_))*dSeqByX_;

     return new_dSeqByX_;
}

//
Foam::List<Foam::scalar> Foam::farFieldPressure::PhiOfEta( point xn_, point xm_, label nLines, List<scalar> new_dSeqByX_)
{
    Foam::List<Foam::scalar> PhiOfEta_(nLines);

    PhiOfEta_[nLines-1] = 0;
    for (int i= (nLines -2); i > -1; i--)
    {
        PhiOfEta_[i] = PhiOfEta_[i+1] +  new_dSeqByX_[i]*mag(xm_[0]-xn_[0])/nLines;
    }

    return PhiOfEta_;
}

//
Foam::List<Foam::scalar> Foam::farFieldPressure::FofEta(label nLines, List<scalar> PhiOfEta_, point xn, point xm)
{
    Foam::List<Foam::scalar> FofEta_(nLines-1);

//    FofEta_[nLines-1] = (PhiOfEta_[nLines-2]-PhiofEta_[nLines-1])*nLines/(xm[0]-xn[0]);
//    FofEta_[0] = (PhiOfEta_[0]-PhiofEta_[1])*nLines/(xm[0]-xn[0]);

    for (int i= 0; i < nLines-1; i++) 
    {
        FofEta_[i] = (PhiOfEta_[i] - PhiOfEta_[i+1])/(xm[0]-xn[0])*nLines;
    }

    return FofEta_;
}


//
Foam::scalar Foam::farFieldPressure::getL( scalar yn_, scalar y_, scalar beta_ef_, scalar M_ef_ )
{
    scalar getL_ = M_ef_ / beta_ef_ * mag( yn_ - y_ );

    return getL_;
}

//
Foam::scalar Foam::farFieldPressure::getYforce()
{
    scalar Yforce_ = 1;

    return Yforce_;
}

//
Foam::scalar Foam::farFieldPressure::k1Coef( scalar M_ef_, scalar M, scalar beta_ef_, scalar yn_, scalar yEarth_, scalar y_, scalar rho_inf)
{
  //List of params

    scalar rho_00 = rho_inf;
    scalar rho_0  = rho_inf;

//    scalar gamma_bar = 0;
    scalar l_        = getL ( yn_, y_,      beta_ef_, M_ef_ );
    scalar l_0_      = getL ( yn_, yEarth_, beta_ef_, M_ef_ );

    scalar Yforce_   = getYforce();

  //end list of params

    scalar k1Coef_ = 0;

    k1Coef_ = 1/constant::mathematical::pi*sqrt(M*rho_00/(2*rho_0*l_0_*mag(1-gamma_bar)))*Yforce_/(l_*sqrt(l_)) ;     

    return k1Coef_;
}

//
Foam::List<Foam::scalar> Foam::farFieldPressure::deltaP(label nLines_, scalar k1_, List<scalar> FofEta_)
{
    Foam::List<Foam::scalar> deltaP_ (nLines);

    for (int i= 0; i < nLines-1; i++)
    {
        deltaP_[i] = k1_*FofEta_[i];
    }

    return deltaP_;
}


/*

//
//Additional part to farFieldPressure lib
//

//
Foam::scalar Foam::farFieldPressure::kCoef()
{
    scalar kCoef_ = 0;

    return kCoef_;
}

//
Foam::scalar Foam::farFieldPressure::k2Coef()
{
    scalar k2Coef_ = 0;

    return k2Coef_;
}

//
Foam::scalar Foam::farFieldPressure::TbyTwo( scalar eta_, scalar kCoef_, scalar k2Coef_)
{
    scalar TbyTwo_ = 0;

    TbyTwo_ = ;

    return TbyTwo_;
}
*/

//************************************************************************** //
