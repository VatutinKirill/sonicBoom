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

#include "IFstream.H"
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



// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::farFieldPressure::farFieldPressure
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

    dict.lookup("interpolationScheme") >> interpolationScheme_;

    infVel     = ( dict.get<vector>("infVel") );
    TotalForce = ( dict.get<vector>("TotalForce") );
    infT       = ( dict.get<scalar>("infT")   );
    Wgamma     = ( dict.get<scalar>("Wgamma") );
    M_ef       = ( dict.get<scalar>("M_ef")   );
    AoA        = ( dict.get<scalar>("AoA")    );
    WR         = ( dict.get<scalar>("WR")     );
    length     = ( dict.get<scalar>("length") );
    
    rho_00    = ( dict.get<scalar>( "rho_00")   );
    rho_0     = ( dict.get<scalar>( "rho_0" )   );
    rho_inf   = ( dict.get<scalar>( "rho_inf")  );
    gamma_bar = ( dict.get<scalar>("gamma_bar") );

    H         = ( dict.get<scalar>("H") );

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

    yn = xn[vector::Y];

    point xm( (xe + xs)*0.5 );
    point xmn( xm-xn );
    xmn[vector::Y]= 0.0;

    infVelMagVal = farFieldPressure::infVelMag();

    point newXmn( farFieldPressure::beta() * abs( xm[vector::Y]-xn[vector::Y] )/sqrt(sqr(mag(xm-xn)) - sqr(xm[vector::Y]-xn[vector::Y]))*xmn );
    
    newXmn[vector::Y] = xm[vector::Y]-xn[vector::Y];
    newXmn[vector::X] = newXmn[vector::X] + mag(xm[vector::Y]-xn[vector::Y])*betaVal;
    point newXn (xn+newXmn);

    Info << "nose:      " << xn                   << nl;
    Info << "newNose:   " << newXn                << nl;
    Info << "xm:        " << xm                   << nl;

    Info << "infVelMag: " << infVelMagVal << nl;
    Info << "beta:      " << farFieldPressure::beta() << nl;
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

//        Info << "line:         " << iLine                << nl;
//        Info << "              " << probeSets_[i]        << nl;
//        Info << "start point:  " << newPointA            << nl;
//        Info << "end point:    " << newPointB            << nl;
//        Info << "line size:    " << probeSets_[i].size() << nl;

        if( probeSets_[i].size()< nPoints)
        {
            FatalErrorInFunction
            << "Error bad line set parameters, lines should be inside of domain"
            << nl << exit(FatalError);
        }
    }

    return true;
}

bool Foam::farFieldPressure::write()
{

    point newPointA;
    point newPointB;

    point xm( (xe + xs)*0.5 );

    scalar y = xm[vector::Y];

    Info << "write xm:  " << xm << nl;

    point xmn( xm-xn );
    xmn[vector::Y]= 0.0;

    point newXmn( farFieldPressure::beta() * abs(xm[vector::Y]-xn[vector::Y])/sqrt( sqr(mag(xm-xn)) - sqr(xm[vector::Y]-xn[vector::Y]))*xmn );
    newXmn[vector::Y] = xm[vector::Y]-xn[vector::Y];
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

//        Info << "line:         " << iLine                << nl;
//        Info << "              " << probeSets_[i]        << nl;
//        Info << "start point:  " << newPointA            << nl;
//        Info << "end point:    " << newPointB            << nl;
//        Info << "line size:    " << probeSets_[i].size() << nl;

        if( probeSets_[i].size()< nPoints)
        {
            FatalErrorInFunction
            << "Error bad line set parameters, lines should be inside of domain"
            << nl << exit(FatalError);
        }
    }


//    volScalarField velField (mesh_.thisDb().lookupObject<volVectorField>("U").component(vector::Y));
    
//    scalar infVelY = infVel[vector::Y];
    
    volScalarField fieldY (mesh_.thisDb().lookupObject<volVectorField>("U").component(vector::Y)); 
    volScalarField fieldX (mesh_.thisDb().lookupObject<volVectorField>("U").component(vector::X));
    volScalarField field2 = fieldY * Foam::cos ( AoA * constant::mathematical::pi/180.0 ) - fieldX* Foam::sin ( AoA * constant::mathematical::pi/180.0);

//    velField -= dimensionedScalar("Uinf",dimLength/dimTime,infVelY);
    


    volFieldSampler<Foam::scalar>  volFieldSampler_(interpolationScheme_, field2, probeSets_);


    //
    Foam::List<Foam::scalar> dSeqBydXVal (nLines);
    for (int i=0; i < nLines; i++)
    {
        dSeqBydXVal[i]= 0.0;

        for (int j=0; j < volFieldSampler_[i].size() ; j++)
        {
            dSeqBydXVal[i] += volFieldSampler_[i][j]* mag(newPointB[vector::Z]-newPointA[vector::Z]) / volFieldSampler_[i].size();
        }
        dSeqBydXVal[i]= -2.0/infVelMagVal*dSeqBydXVal[i];
     Info <<  "Line1" << nl; 
     Info <<  "Line:"<< i  << "  dSeqBydXVal[i]:    " <<  dSeqBydXVal[i] << nl;
    }


    //
    betaVal = beta( );


    //
    Foam::List<Foam::scalar> new_dSeqBydXVal(nLines);
    for (int i=0; i < nLines; i++)
    {
        new_dSeqBydXVal[i] = 0.0;
        new_dSeqBydXVal[i] = new_dSeqBydX( dSeqBydXVal[i]);
    //    write results  
    Info <<  "Line:"<< i  << " new_ dSeqBydXVal[i]:    " <<  new_dSeqBydXVal[i] << nl;
    }

    //
    Foam::List<Foam::scalar> PhiOfEtaVal(nLines);
    PhiOfEtaVal = PhiOfEta( xn, xm, nLines, new_dSeqBydXVal);

    //    write results  
          Info << "PhiOfEtaVal:    " << PhiOfEtaVal << nl;

    List<scalar> FofEtaVal(nLines-1);
    FofEtaVal = FofEta(nLines, PhiOfEtaVal, xn, xm);

    //    write results  
        Info << "FofEtaVal:    " << FofEtaVal << nl;

    beta_ef   = beta();

    //
    scalar k1 = k1Coef( y );

    //
    List<scalar> res_dp = deltaP(nLines - 1, k1, FofEtaVal);

    //    write results  
    Info << "res_dp:    " << res_dp << nl;

    makeFile(res_dp);
//    farFieldPressureDir
//                       /time
//                            /coord( x y z) dp

    return true;
}

bool Foam::farFieldPressure::execute()
{
    return true;
}


void Foam::farFieldPressure::makeFile(List<scalar> res_dp)
{
    word farFieldPressureDir;
    std::string fileName;

    fileName = "/farFieldDp ";

    if (Pstream::master() && Pstream::parRun())
    {
        //farFieldSonicBoomDir = obr_.time().rootPath() + "/" + obr_.time().caseName().path() + "/data" ;
        farFieldPressureDir = mesh_.time().rootPath() + "/" + mesh_.time().caseName().path() + "/data" ;

        mkDir(farFieldPressureDir);
    }
    else if (!Pstream::parRun())
    {
        //farFieldSonicBoomDir = obr_.time().rootPath() + "/" + obr_.time().caseName() + "/data" ;
        farFieldPressureDir = mesh_.time().rootPath() + "/" + mesh_.time().caseName() + "/data" ;

        mkDir(farFieldPressureDir);
    }

//    std::string dirName = farFieldPressureDir + fileName.c_str();

//    Info << nl <<"Field: "<< fieldName_ << ": Saving "<< fileName << " to file.";

//    Foam::OFstream theFile(dirName);
  //  theFile << Foam::scientific << Foam::setprecision(10);


   // Info    << "res_dp " << res_dp  << nl;
   // theFile <<  res_dp              << nl;
}

//
Foam::scalar Foam::farFieldPressure::infVelMag()
{
    scalar infVelMagnitude = 0.0;
    infVelMagnitude        = sqrt(infVel&infVel);

    return infVelMagnitude;
}

//
Foam::scalar Foam::farFieldPressure::beta()
{
    scalar cInfSqr  = Wgamma*WR*infT;

//    Info << "infVelMag_ : " << infVelMagVal << nl;
//    Info << "cInfSqr    : " << cInfSqr      << nl;

    if( mag(infVelMagVal) < sqrt(cInfSqr))
    {
        Info << "initial velocity less than sound velocity. Beta not defiend. " << nl;
    }

    scalar final_beta    = sqrt(sqr(infVelMagVal)/cInfSqr - 1.0);
    return final_beta;
}

/*
Foam::scalar Foam::farFieldPressure::dSeqBydX(const vector infVel_, const scalar beta_, const vector new_pointA_, const vector new_pointB_, const label nPoints)
{
    
}
*/

//
Foam::scalar  Foam::farFieldPressure::new_dSeqBydX(const scalar dSeqByX_)
{

//     scalar new_dSeqByX_ = sqr(infVelMagVal)/WR/infT/(2*length*constant::mathematical::pi*sqrt(2*betaVal))*dSeqByX_;

     scalar new_dSeqByX2 = dSeqByX_*Wgamma*M_ef*M_ef/(2.0*length*constant::mathematical::pi*sqrt(2.0*betaVal));

//    Info << " new_dSeqByX_    :  " <<  new_dSeqByX_    << nl;

    Info << " new_dSeqByX2    :  " <<  new_dSeqByX2    << nl;

     return new_dSeqByX2;
}

//
Foam::List<Foam::scalar> Foam::farFieldPressure::PhiOfEta( point xn_, point xm_, label nLines, List<scalar> new_dSeqByX_)
{
    Foam::List<Foam::scalar> PhiOfEta_(nLines);

    PhiOfEta_ = 0.0;

    scalar dx  = (xm_.x() - xn_.x()) / scalar(nLines) / length;

    for (int i=1; i<nLines; i++)
    {
        scalar eta = dx * static_cast<scalar>(i);
        
        for (int k=0; k<i; k++)
        {
            scalar xi  = dx * static_cast<scalar>(k);        
            PhiOfEta_[i] += 
                new_dSeqByX_[nLines-1-k]*dx/sqrt(eta - xi);
        }
    }

    return PhiOfEta_;
}

//
Foam::List<Foam::scalar> Foam::farFieldPressure::FofEta(label nLines, List<scalar> PhiOfEta_, point xn, point xm)
{
    Foam::List<Foam::scalar> FofEta_(nLines-1,0.0);

    const scalar dx = (xm.x() - xn.x())/static_cast<scalar>(nLines) / length ;

    for (int i= 0; i < nLines-1; i++)
    {
        FofEta_[i] = -(PhiOfEta_[i] - PhiOfEta_[i+1])/dx;
    }

    return FofEta_;
}


//
Foam::scalar Foam::farFieldPressure::getL( scalar yn_, scalar y_, scalar beta_ef_, scalar M_ef_ )
{ 
    scalar LVal = 0.0;

    Info << " yn_      :  " <<  yn_      << nl;
    Info << " y_       :  " <<  y_       << nl;
    Info << " beta_ef_ :  " <<  beta_ef_ << nl;
    Info << " M_ef_    :  " <<  M_ef_    << nl;

    LVal = M_ef_ / beta_ef_ * mag( yn_ - y_ );

    Info << " L(yn - y) :    " <<  LVal << nl;

    return LVal;
}

//
Foam::scalar Foam::farFieldPressure::getYforce()
{

//    Foam::functionObjects::forces::calcForcesMoment();
    //calcForcesMoment();

//    scalar Yforce_ = Foam::functionObjects::forces::forceEff();
//    vector dotF = F_.dot(obr_.time().value(), 0);
//
//scalar Yforce_ = 0.15;

//    scalar TotalForceX = 76.36467358;
//    scalar TotalForceY = 53.33356928;

    scalar TotalForceX = TotalForce[vector::X];
    scalar TotalForceY = TotalForce[vector::Y];

//    scalar TotalForceY = 53.33356928;

    scalar Yforce_  = (TotalForceY * Foam::cos( AoA*constant::mathematical::pi/180.0) - TotalForceX*Foam::sin( AoA*constant::mathematical::pi/180.0)) * 2.0;

    Info << " Yforce__ :    " <<  Yforce_ << nl;

    return Yforce_;
}

//
Foam::scalar Foam::farFieldPressure::k1Coef( scalar y_)
{
  //List of params

    scalar rho_00 = rho_inf;
    scalar rho_0  = rho_inf;

    scalar l    = getL ( H, y_, beta_ef, M_ef );
    scalar l_0  = getL ( H,  0, beta_ef, M_ef );

    scalar YforceVal   = getYforce();

    Info << " l   : " << l   << nl;
    Info << " l_0 : " << l_0 << nl;


  //end list of params

    scalar k1CoefVal = 0.0;

    Info << " rho_0  : " << rho_0  << nl;
    Info << " rho_00 : " << rho_00 << nl;


 //   scalar k0Coef_ = sqrt(M_ef_*rho_00/(2.0*rho_0*l_0_*mag(1.0-gamma_bar)))*Yforce_/(l_*sqrt(l_)) ;

//    Info << " k_Coef_ :    " <<  k0Coef_ << nl;


    k1CoefVal = 1.0/constant::mathematical::pi*sqrt(M_ef*rho_00/(2.0*rho_0*l_0*mag(1.0-gamma_bar)))*YforceVal/(l*sqrt(l)) ;     

    Info << " k1:    " <<  k1CoefVal << nl;

    return k1CoefVal;
}

//
Foam::List<Foam::scalar> Foam::farFieldPressure::deltaP(label nLines_, scalar k1_, List<scalar> FofEta_)
{
    Foam::List<Foam::scalar> deltaP_ (nLines_);

    for (int i= 0; i < nLines_; i++)
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
*////************************************************************************////
  //**************************************************************************//
 ////************************************************************************////
