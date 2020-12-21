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

//#include "sampledSets.H"
#include "volFields.H"
#include "globalIndex.H"


// * * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * //

// Function from sampledSets.H

template<class Type>
Foam::farFieldPressure::volFieldSampler<Type>::volFieldSampler
(
     const word& interpolationScheme,
     const GeometricField<Type, fvPatchField, volMesh>& field,
     const PtrList<sampledSet>& samplers
 )
 :
     List<Field<Type>>(samplers.size()),
     name_(field.name())
 {
     autoPtr<interpolation<Type>> interpolator
     (
         interpolation<Type>::New(interpolationScheme, field)

     );

     forAll(samplers, setI)
     {
         Field<Type>& values = this->operator[](setI);
         const sampledSet& samples = samplers[setI];
  
         values.setSize(samples.size());
         forAll(samples, sampleI)
         {
             const point& samplePt = samples[sampleI];
             label celli = samples.cells()[sampleI];
             label facei = samples.faces()[sampleI];
  
             if (celli == -1 && facei == -1)
             {

                 values[sampleI] = pTraits<Type>::max;
             }
             else
             {
                 values[sampleI] = interpolator().interpolate
                 (
                     samplePt,
                     celli,
                     facei
                 );
             }
         }
     }
 }

/*
(
    const GeometricField<Type, fvPatchField, volMesh>& field,
    const PtrList<sampledSet>& samplers
)
:
    List<Field<Type>>(samplers.size()),
    name_(field.name())
{
    forAll(samplers, setI)
    {
        Field<Type>& values = this->operator[](setI);
        const sampledSet& samples = samplers[setI];

        values.setSize(samples.size());
        forAll(samples, sampleI)
        {
            label celli = samples.cells()[sampleI];

            if (celli ==-1)
            {
                values[sampleI] = pTraits<Type>::max;
            }
            else
            {
                values[sampleI] = field[celli];
            }
        }
    }
}
*/
// ************************************************************************* //
