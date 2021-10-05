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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::label Foam::traIdealControl::nCorrPIMPLE() const
{
    return nCorrPIMPLE_;
}

inline Foam::label Foam::traIdealControl::nCorrPISO() const
{
    return nCorrPISO_;
}

inline Foam::label Foam::traIdealControl::N1() const
{
    return N1_;
}


inline Foam::label Foam::traIdealControl::N2() const
{
    return N2_;
}


inline Foam::label Foam::traIdealControl::corrPISO() const
{
    return corrPISO_;
}


inline bool Foam::traIdealControl::SIMPLErho() const
{
    return SIMPLErho_;
}


inline bool Foam::traIdealControl::correct()
{
    corrPISO_++;

    if (debug)
    {
        Info<< algorithmName_ << " correct: corrPISO = " << corrPISO_ << endl;
    }

    if (corrPISO_ <= nCorrPISO_)
    {
        return true;
    }
    else
    {
        corrPISO_ = 0;
        return false;
    }
}


inline bool Foam::traIdealControl::storeInitialResiduals() const
{
    // Start from second PIMPLE iteration
    return (corr_ == 2) && (corrPISO_ == 0) && (corrNonOrtho_ == 0);
}


inline bool Foam::traIdealControl::firstIter() const
{
    return corr_ == 1;
}


inline bool Foam::traIdealControl::finalIter() const
{
    return converged_ || (corr_ == nCorrPIMPLE_);
}


inline bool Foam::traIdealControl::finalInnerIter() const
{
    return
       corrPISO_ == nCorrPISO_
    && corrNonOrtho_ == nNonOrthCorr_ + 1;
}


inline bool Foam::traIdealControl::solveFlow() const
{
    return solveFlow_;
}


inline bool Foam::traIdealControl::turbCorr() const
{
    return !turbOnFinalIterOnly_ || finalIter();
}


// ************************************************************************* //
