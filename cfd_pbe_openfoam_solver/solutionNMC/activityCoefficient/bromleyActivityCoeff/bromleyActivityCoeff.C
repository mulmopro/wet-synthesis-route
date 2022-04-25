/* ----------------------------------------------------------------------------
Copyright © 2021 Politecnico di Torino

This file is part of WetSynthRoute.

WetSynthRoute is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

WetSynthRoute is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with WetSynthRoute.  If not, see <https://www.gnu.org/licenses/>.

This offering is not approved or endorsed by OpenCFD Limited (ESI Group)
and the OpenFOAM Foundation.
---------------------------------------------------------------------------- */

#include "bromleyActivityCoeff.H"
#include "solutionNMC.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace activityCoeffModels
{
    defineTypeNameAndDebug(bromleyActivityCoeff, 0);
    addToRunTimeSelectionTable(activityCoeffModel, bromleyActivityCoeff, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar
Foam::activityCoeffModels::bromleyActivityCoeff::pure_pair_activity_coeff
(
    label id_cation, label id_anion, scalar z_c, scalar z_a, scalar I_s
) const
{
    scalar B = B_(id_cation, id_anion);

    const scalar& nu_c = z_a;
    const scalar& nu_a = z_c;

    scalar zByZ = (nu_c*z_c*z_c + nu_a*z_a*z_a) / (nu_c + nu_a);

    scalar sqrt_I = Foam::sqrt(I_s);

    scalar log_gamma = 
        -1.0 * A_gamma_ * zByZ * sqrt_I / (1 + sqrt_I)
      + (0.06 + 0.6*B) * zByZ * I_s / Foam::sqr(1.0 + 1.5 * I_s / zByZ)
      + B * I_s;

    scalar E = E_(id_cation, id_anion);

    if (E)
    {
        if (nu_c == 2 && nu_a == 2)
        {
            log_gamma -= 70.0 * E * sqrt_I * (1.0 - Foam::exp(-70.0 * sqrt_I));
        }
    }

    return log_gamma;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::activityCoeffModels::bromleyActivityCoeff::bromleyActivityCoeff
(
    const dictionary& dict,
    const solutionNMC& solution
)
:
    activityCoeffModel(dict, solution),

    A_gamma_(0.511),

    nCation_(solution.nMetals() + 3),
    nAnion_(2),

    z_a_({2, 1})
{
    double B_allMetals[3][2];

    B_allMetals[0][0] = 0.1056; B_allMetals[0][1] = -0.080;
    B_allMetals[1][0] = 0.1226; B_allMetals[1][1] = -0.097;
    B_allMetals[2][0] = 0.1244; B_allMetals[2][1] = -0.085;

    double B_nonMetals[3][2];
    B_nonMetals[0][0] = -0.0204; B_nonMetals[0][1] = 0.0747;
    B_nonMetals[1][0] = -0.0287; B_nonMetals[1][1] = 0.0540;
    B_nonMetals[2][0] = 0.0606; B_nonMetals[2][1] = 0.0605;

    double E_allMetals[3][2];
    E_allMetals[0][0] = 0.00524; E_allMetals[0][1] = 0.0;
    E_allMetals[1][0] = 0.00599; E_allMetals[1][1] = 0.0;
    E_allMetals[2][0] = 0.00498; E_allMetals[2][1] = 0.0;

    label nMetals(solution.nMetals());

    B_.setSize(nCation_, nAnion_);
    E_.setSize(nCation_, nAnion_);

    for (int j=0; j < nAnion_; j++)
    {
        int i;
        for (i=0; i < nMetals; i++)
        {
            B_(i, j) = B_allMetals[i][j];
            E_(i, j) = E_allMetals[i][j];
        }

        for (i=nMetals; i < nCation_; i++)
        {
            B_(i, j) = B_nonMetals[i - nMetals][j];
            E_(i, j) = 0.0;
        }
    }

    z_c_.setSize(nCation_);

    int i;
    for (i=0; i < nMetals; i++)
    {
        z_c_[i] = 2;
    }

    for (i=nMetals; i < nCation_; i++)
    {
        z_c_[i] = 1;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::activityCoeffModels::bromleyActivityCoeff::~bromleyActivityCoeff()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::activityCoeffModels::bromleyActivityCoeff::ionic_strength
(
    const Foam::List<Foam::scalar>& cationMolalConc,
    const Foam::List<Foam::scalar>& anionMolalConc
) const
{
    scalar I_s = 0.0;

    int i;
    for (i=0; i < nCation_; i++)
    {
        I_s += cationMolalConc[i] * Foam::sqr(z_c_[i]);
    }

    for (i=0; i < nAnion_; i++)
    {
        I_s += anionMolalConc[i] * Foam::sqr(z_a_[i]);
    }

    return 0.5*I_s;
}


Foam::scalar Foam::activityCoeffModels::bromleyActivityCoeff
::pair_activity_coeff
(
    Foam::label id_metal,
    const Foam::List<Foam::scalar>& cationMolalConc,
    const Foam::List<Foam::scalar>& anionMolalConc,
    Foam::scalar I_s
) const
{
    label id_OH = 1;

    scalar z_c = z_c_[id_metal];
    scalar z_a = z_a_[id_OH];

    scalar log_gamma_pure_ac = pure_pair_activity_coeff(
        id_metal, id_OH, z_c, z_a, I_s);

    scalar sqrt_I = Foam::sqrt(I_s);

    scalar coeff = A_gamma_ * sqrt_I / (1 + sqrt_I);

    scalar z_i;
    scalar log_gamma_pure;

    scalar y_ac;
    scalar f_c = 0.0;
    int i;
    for (i=0; i < nAnion_; i++)
    {
        if (i != id_OH)
        {
            z_i = z_a_[i];

            log_gamma_pure = pure_pair_activity_coeff(
                id_metal, i, z_c, z_i, I_s);

            y_ac = Foam::sqr(z_c + z_i) * anionMolalConc[i] / (4 * I_s);

            f_c += y_ac * (log_gamma_pure + coeff*z_c*z_i);
        }
        else
        {
            y_ac = Foam::sqr(z_c + z_a) * anionMolalConc[i] / (4 * I_s);

            f_c += y_ac * (log_gamma_pure_ac + coeff*z_c*z_a);
        }
    }

    scalar x_ac;
    scalar f_a = 0.0;
    for (i=0; i < nCation_; i++)
    {
        if (i != id_metal)
        {
            z_i = z_c_[i];

            log_gamma_pure = pure_pair_activity_coeff(
                i, id_OH, z_i, z_a, I_s);

            x_ac = Foam::sqr(z_i + z_a) * cationMolalConc[i] / (4 * I_s);

            f_a += x_ac * (log_gamma_pure + coeff*z_i*z_a);
        }
        else
        {
            x_ac = Foam::sqr(z_c + z_a) * cationMolalConc[i] / (4 * I_s);

            f_a += x_ac * (log_gamma_pure_ac + coeff*z_c*z_a);
        }
    }

    const scalar& nu_c = z_a;
    const scalar& nu_a = z_c;

    scalar log_gamma =
            (nu_c*f_c + nu_a*f_a - coeff * (nu_c*z_c*z_c + nu_a*z_a*z_a))
          / (nu_c + nu_a);

    scalar gamma_ac =
        Foam::pow
        (
            10,
            log_gamma
        );

    return gamma_ac;
}


// ************************************************************************* //
