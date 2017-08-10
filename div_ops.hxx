/*
  Finite volume discretisations of divergence operators
 
  ***********

    Copyright B.Dudson, J.Leddy, University of York, September 2016
              email: benjamin.dudson@york.ac.uk

    This file is part of Hermes.

    Hermes is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Hermes is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Hermes.  If not, see <http://www.gnu.org/licenses/>.
  
 */

#ifndef __DIV_OPS_H__
#define __DIV_OPS_H__

#include <field3d.hxx>
#include <vector3d.hxx>

// Uses flux-conservative form with ZIP interpolation to cell faces
const Field3D Div_n_a_bxGrad_f_B(const Field3D &n, const Field3D &a, const Field3D &f, bool xflux=false, bool yderiv=true);

// This uses a combination of two flux-conservative discretisations to get an operator
// which is almost, but not quite, anti-symmetric
const Field3D Div_n_bxGrad_f_B(const Field3D &n, const Field3D &f);

// Div ( a * Grad_perp(f) )
const Field3D Div_Perp_Lap(const Field3D &a, const Field3D &f);

const Field3D Div_Perp_Lap_x3(const Field3D &a, const Field3D &f, bool xflux=false);

const Field3D Div_Perp_Lap_XYZ(const Field3D &a, const Field3D &f, bool bndryflux=false);


/*!
 * Parallel diffusion (in y)
 *
 * Calculated in terms of fluxes through cell faces. Takes the average
 * coefficient K from the cells either side, and the gradient across the 
 * boundary. Calculated flux is added to one cell, subtracted from the other.
 *
 * Div_par( K Grad_par(f) )
 *
 * @param[in] K The diffusion coefficient
 * @param[in] f The variable to be differentiated
 * @param[in] bndry_flux  Are fluxes calculated through Y boundaries?
 *
 */
const Field3D Div_par_diffusion(const Field3D &k, const Field3D &f, bool bndry_flux=true);

/*!
 * Parallel heat conduction, assuming a heat conduction coefficient
 * K which depends on the temperature Te^2.5
 * 
 * Div_par( K0 Te^2.5 Grad_par(Te) )
 *
 * To calculate K0*Te^2.5 the temperature is averaged from cell centre
 * to cell boundary.
 *
 * @param[in] K0  Constant coefficient in the conductivity
 * @param[in] Te  Temperature
 * @param[in] bndry_flux  Are fluxes through the boundary calculated?
 */
const Field3D Div_par_spitzer(BoutReal K0, const Field3D &Te, bool bndry_flux=true);

/*!
 * Diffusion using upwinding of the conduction coefficient
 *
 * Depending on the sign of the gradient, the value of K from the 
 * "upwind" side is used in calculating the flux, rather than taking 
 * the average of upstream and downstream sides. 
 *
 * Div_par( K Grad_par(f) )
 * 
 * @param[in] K  The diffusion coefficient
 * @param[in] f  The variable which is differentiated
 * @param[in] bndry_flux   Are boundary fluxes calculated?
 */
const Field3D Div_par_diffusion_upwind(const Field3D &K, const Field3D &f, bool bndry_flux=true);

/*!
 * Diffusion in index space
 * 
 * Similar to using Div_par_diffusion(SQ(mesh->dy)*mesh->g_22, f)
 *
 * @param[in] The field to be differentiated
 * @param[in] bndry_flux  Are fluxes through the boundary calculated?
 */
const Field3D Div_par_diffusion_index(const Field3D &f, bool bndry_flux=true);

/*!
 * Added Dissipation scheme (related to Momentum Interpolation)
 *
 * This uses a 3rd-order derivative of the pressure as
 * a correction to the velocity. 
 *
 * This should appear in the form
 * 
 * df/dt = ... + AddedDissipation(N, P, f);
 */
const Field3D AddedDissipation(const Field3D &N, const Field3D &P, const Field3D f, bool bndry_flux=true);

const Field3D Div_n_bxGrad_f_B_XPPM(const Field3D &n, const Field3D &f, bool bndry_flux=true, bool poloidal=false, bool positive=false);

const Field3D Div_f_v_XPPM(const Field3D &n, const Vector3D &v, bool bndry_flux=true, bool positive=false);

void communicateFluxes(Field3D &f);

const Field3D Div_Perp_Lap_FV(const Field3D &n, const Field3D &f, bool xflux);
const Field3D Div_Perp_Lap_FV_Index(const Field3D &a, const Field3D &f, bool xflux);

/*!
 * Finite volume parallel divergence
 *
 * Assumes there are (at least) two guard cells (MYG >= 2)
 * 
 * @param[in] f  The field being advected
 * @param[in] v  The advection velocity
 */
const Field3D Div_par_FV(const Field3D &f, const Field3D &v);

/*!
 * Parallel divergence, flux splitting version
 *
 * @param[in] f   The field being advected
 * @param[in] v   The advection velocity
 * @param[in] a   Maximum wave speed. Used to determine the amount of upwinding
 *
 * Split into fluxes with speed v+a and v-a
 */
const Field3D Div_par_FV_FS(const Field3D &f, const Field3D &v, const Field3D &a, bool fixflux = true);

// Finite volume parallel divergence of a flow velocity
const Field3D Div_parV_FV(const Field3D &v);

/*!
 * 4th-order derivative
 *
 * Implemented as a flux through cell boundaries, calculated
 * using one-sided 3rd derivative at the boundary.
 *
 * @param[in]  d  Coefficient, averaged from neighbouring cells
 * @param[in]  f  The field being differentiated
 */
const Field3D D4DY4_FV(const Field3D &d, const Field3D &f, bool bndry_flux=false);

/*
 * 4th-order dissipation term
 * 
 *
 * A one-sided 3rd-order derivative, given a value
 * at a boundary is:
 * 
 * d3f/dx3 ~= 16/5 f_b - 6 f_0 + 4 f_1 - 6/5 f_2
 * 
 * where f_b is the value on the boundary; f_0 is the cell 
 * to the left of the boundary; f_1 to the left of f_0 and f_2
 * to the left of f_1
 *
 *    f_2 | f_1 | f_0 |
 *                   f_b
 */
const Field3D D4DY4_FV_Index(const Field3D &f, bool bndry_flux=false);

// 4th-order flux conserving term, in index space
const Field3D D4DX4_FV_Index(const Field3D &f, bool bndry_flux=false);

// Div ( k * Grad(f) )
const Field2D Laplace_FV(const Field2D &k, const Field2D &f);

#endif //  __DIV_OPS_H__
