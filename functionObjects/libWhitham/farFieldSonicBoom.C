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


#include "soundPressureSampler.H"
#include "volFields.H"
#include "dictionary.H"
#include "Time.H"
#include "argList.H"

//sampledSurfaces stuff
//#include "meshedSurfRef.H"
//#include "IOmanip.H"
//#include "volPointInterpolation.H"
#include "PatchTools.H"
#include "ListListOps.H"
//#include "stringListOps.H"
//#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(farFieldSonicBoom, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        farFieldSonicBoom,
        dictionary
    );
}


Foam::scalar Foam::farFieldSonicBoom::mergeTol_ = 1e-10;

// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * //


void Foam::farfieldSonicBoom::writeGeometry() const
{
    const fileName outputDir = "surfaceGeometryData";

    Info << "write geometry \n";

    Info << controlSurfaces_.size() << nl;


    forAll(controlSurfaces_, surfI)
    {
        Info << "in surfaces list\n";

        const sampledSurface& s = controlSurfaces_.operator[](surfI);


        if (Pstream::parRun())
        {
            Info << "surfI size = " << mergeList_[surfI].points.size() << nl;

            if (Pstream::master() && mergeList_[surfI].faces.size())
            {
                formatter_->write
                (
                    outputDir,
                    s.name(),
                    meshedSurfRef
                    (
                        mergeList_[surfI].points,
                        mergeList_[surfI].faces
                    )
                );
            }
        }
        else if (s.points().size())
        {
            Info << s.points().size() << nl;

            formatter_->write
            (
                outputDir,
                s.name(),
                meshedSurfRef
                (
                    s.points(),
                    s.faces()
                )
            );
        }
        else
        {
            Info << "empty surface info" << nl;
        }
    }
}



// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::farFieldSonicBoom::farFieldSonicBoom
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
    loadFromFiles_(false),
    name_(name),
    log_(false),
    controlSurfaces_(0),
    outputFormat_(word::null),
    fieldName_(word::null)
{
    Info << "in constructor" << nl;

    read(dict);
}

Foam::farFieldSonicBoom::farFieldSonicBoom
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    functionObject(name),
    name_(name),
    log_(false),
    controlSurfaces_(0),
    outputFormat_(word::null),
    fieldName_(word::null)
{
    read(dict);
}

// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::farFieldSonicBoom::~farFieldSonicBoom()
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

   	dict.lookup("interpolationScheme") >> interpolationScheme_;

	dict.lookup("timeIndexRead") >> timeIndexRead_;
	dict.lookup("nSnapshots") >> nSnapshots_;
	dict.lookup("nmodes") >> nmodes_;
	dict.lookup("poddebug") >> poddebug_;
	dict.lookup("PODStartTime") >> PODStartTime_;

    const word writeType (dict.lookup("outputGeometryFormat"));

    //read mesh
    //const fvMesh& mesh_ = refCast<const fvMesh>(obr_);

        // Define the surface formatter
        // Optionally defined extra controls for the output formats
        formatter_ = surfaceWriter::New
        (
            writeType,
            dict.subOrEmptyDict("formatOptions").subOrEmptyDict(writeType)
        );



    //read surfaces for sampling
    PtrList<sampledSurface> newList
    (
        dict.lookup("surfaces"),
        sampledSurface::iNew(mesh_)
    );

    controlSurfaces_.transfer(newList);

    //Parallel fix as it was implemented in sampledSurfaces class
    if (Pstream::parRun())
    {
        mergeList_.setSize(controlSurfaces_.size());
    }
        
    // Ensure all surfaces and merge information are expired
    expire();

    if (controlSurfaces_.size())
    {

        Info<< "Function object "<< name_<<":" << nl;

        Info<< " Reading control surface description:" << nl;
            
        forAll(controlSurfaces_, surfI)
        {
            Info<< " " << controlSurfaces_.operator[](surfI).name() << nl;
            Info<< "needsUpdate: " << controlSurfaces_.operator[](surfI).needsUpdate() << nl;
            Info<< "size: " << controlSurfaces_.operator[](surfI).points().size() << nl;
        }

        Info<< endl;

        //pHistory_.resize(controlSurfaces_.size());

        update();

        writeGeometry();

        Info << "writeGeometry OK" << nl;
    }

    dict.lookup("fieldName") >> fieldName_;

    return true;
}

void Foam::soundPressureSampler::makeFile()
{
    fileName soundPressureSamplerDir;

    if (Pstream::master() && Pstream::parRun())
    {
        //farFieldSonicBoomDir = obr_.time().rootPath() + "/" + obr_.time().caseName().path() + "/data" ;
        farFieldSonicBoomDir = mesh_.time().rootPath() + "/" + mesh_.time().caseName().path() + "/data" ;

        mkDir(farFieldSonicBoomDir);
    }
    else if (!Pstream::parRun())
    {
        //farFieldSonicBoomDir = obr_.time().rootPath() + "/" + obr_.time().caseName() + "/data" ;
        farFieldSonicBoomDir = mesh_.time().rootPath() + "/" + mesh_.time().caseName() + "/data" ;

        mkDir(farFieldSonicBoomDir);
    }

    // File update
    if (Pstream::master() || !Pstream::parRun())
    {
        // Create the soundPressureSampler file if not already created
        farFieldSonicBoomFilePtr_.resize(controlSurfaces_.size());

        forAll(controlSurfaces_, surfI)
        {
            sampledSurface& s = controlSurfaces_.operator[](surfI);

            fileName pFileName = soundPressureSamplerDir + "/" + (s.name() + "_" + name_ + ".txt");

            if (isFile(pFileName))
            {
                Info << "File " << pFileName << " is already exists" << nl;
            }

            if (farFieldSonicBoomFilePtr_[surfI].empty())
            {
                // Open new file at start up
                Info << pFileName << nl;

                farFieldSonicBoomFilePtr_[surfI].reset
                (
                    new OFstream
                    (
                        pFileName
                    )
                );
            }
        }
    }

    //write();
}

Foam::scalar Foam::farFieldSonicBoom::beta(const scalar Mref)
{
    
}

Foam::scalar Foam::farFieldSonicBoom::dSeqByX()
{

}

Foam::scalar Foam::farFieldSonicBoom::new_dSeqByX()
{

}

Foam::scalar Foam::farFieldSonicBoom::PhiOfEta()
{

}

Foam::scalar Foam::farFieldSonicBoom::FoEta()
{

}

Foam::scalar Foam::farFieldSonicBoom::kCoef()
{

}

Foam::scalar Foam::farFieldSonicBoom::k1Coef()
{

}

Foam::scalar Foam::farFieldSonicBoom::deltaP()
{

}
// ************************************************************************* //
