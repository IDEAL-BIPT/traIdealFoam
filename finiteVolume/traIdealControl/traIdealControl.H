/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::traIdealControl

Description
    traIdeal control class to supply convergence information/checks for
    the traIdeal loop.

\*---------------------------------------------------------------------------*/

#ifndef traIdealControl_H
#define traIdealControl_H

#include "solutionControl.H"

//- Declare that pimpleControl will be used
#define TRAIDEAL_CONTROL

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class traIdealControl Declaration
\*---------------------------------------------------------------------------*/

class traIdealControl
:
    public solutionControl
{

protected:

    // Protected Data

        // Solution controls

            //- Flag to indicate whether to solve for the flow
            bool solveFlow_;

            //- Maximum number of PIMPLE correctors
            label nCorrPIMPLE_;//nCorrTRAIDEAL_;

            //- Maximum number of PISO correctors
            label nCorrPISO_;

            //- Current PISO corrector
            label corrPISO_;

            //- Flag to indicate whether to update density in SIMPLE
            //  rather than PISO mode
            bool SIMPLErho_;

            //- Flag to indicate whether to only solve turbulence on final iter
            bool turbOnFinalIterOnly_;

            //- Converged flag
            bool converged_;

        //- Initialised flag  //this is not in pimpleControl
        //bool initialised_;

        //- The number of first inner iteration process for pressure equation
        label N1_;

        //- The number of second inner iteration process for pressure equation
        label N2_;


    // Protected Member Functions

        //- Read controls from fvSolution dictionary
        virtual void read();

        //- Return true if all convergence checks are satisfied
        virtual bool criteriaSatisfied();


private:

    // Private member functions

        //- Disallow default bitwise copy construct
        traIdealControl(const traIdealControl&);

        //- Disallow default bitwise assignment
        void operator=(const traIdealControl&);


public:

    // Static Data Members

        //- Run-time type information
        TypeName("traIdealControl");


    // Constructors

        //- Construct from mesh
        traIdealControl(fvMesh& mesh, const word& dictName="TRAIDEAL");
        //traIdealControl(fvMesh& mesh);


    //- Destructor
    virtual ~traIdealControl();


   // Member Functions

        // Access

            //- Maximum number of PIMPLE correctors
            inline label nCorrPIMPLE() const;

            //- Maximum number of PISO correctors
            inline label nCorrPISO() const;

            //- Current PISO corrector index
            inline label corrPISO() const;

            //- Flag to indicate whether to update density in SIMPLE
            //  rather than PISO mode
            inline bool SIMPLErho() const;

            //- The number of first inner iteration process for pressure equation
            inline label N1() const;

            //- The number of second inner iteration process for pressure equation
            inline label N2() const;


        // Solution control

            //- PIMPLE loop
            virtual bool loop();

            //- Pressure corrector loop control
            inline bool correct();

            //- Return true to store the intial residuals
            inline bool storeInitialResiduals() const;

            //- Return true for first PIMPLE (outer) iteration
            inline bool firstIter() const;

            //- Return true fore final PIMPLE (outer) iteration
            inline bool finalIter() const;

            //- Return true for final inner iteration
            inline bool finalInnerIter() const;

            //- Return true to solve for flow
            inline bool solveFlow() const;

            //- Return true to solve for turbulence
            inline bool turbCorr() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "traIdealControlI.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
