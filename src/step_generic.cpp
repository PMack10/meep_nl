/* Copyright (C) 2005-2025 Massachusetts Institute of Technology
%
%  This program is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation; either version 2, or (at your option)
%  any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software Foundation,
%  Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
*/

#include "meep.hpp"
#include "meep_internals.hpp"
#include "config.h"
#include "newton_raphson.hpp"
#include <iostream>
#include <unistd.h>
#include <iterator>

#define RPR realnum *restrict

using namespace std;

namespace meep {

#define SWAP(t, a, b)                                                                              \
  {                                                                                                \
    t xxxx = a;                                                                                    \
    a = b;                                                                                         \
    b = xxxx;                                                                                      \
  }

/* update step for df/dt = curl g,
   i.e. f += dt curl g = dt/dx (dg1 - dg2)
   where dgk = gk[i] - gk[i+sk].

   g = (g1,g2), where g1 or g2 may be NULL.  Note that dt/dx and/or s1
   and s2 may be negative to flip signs of derivatives.

   PML: sig[k] = sigma[k]*dt/2, siginv[k] = 1 / (kap[k] + sigma[k]*dt/2).
   Here, k is the index in the dsig direction.  if dsig ==
   NO_DIRECTION, then PML is not used.  (dsig is the sigma direction.)

   if non-NULL, then cnd is an array of conductivity values, changing
   the underlying PDE to:
       df/dt = curl g - cnd f
   which is updated as:
       f = [ dt * curl g + (1 - dt cnd/2) f ] / (1 + dt cnd/2)
   cndinv should be an array of 1 / (1 + dt cnd/2).  In the case
   of PML, cndinv should contain 1 / (1 + dt (cnd + sigma)/2).

   fcnd is an auxiliary field used ONLY when we simultaneously have
   PML (dsig != NO_DIR) and conductivity, in which case fcnd solves
       dfcnd/dt = curl g - cnd*fcnd
   and f satisfies
       df/dt = dfcnd/dt - sigma*f.

   fu is another auxiliary field used only in PML (dsigu != NO_DIR),
   in which case f solves:
       df/dt = dfu/dt - sigma_u * f
   and fu replaces f in the equations above (fu += dt curl g etcetera).
*/
void step_curl(RPR f, component c, const RPR g1, const RPR g2, ptrdiff_t s1,
               ptrdiff_t s2, // strides for g1/g2 shift
               const grid_volume &gv, const ivec is, const ivec ie, realnum dtdx, direction dsig,
               const RPR sig, const RPR kap, const RPR siginv, RPR fu, direction dsigu,
               const RPR sigu, const RPR kapu, const RPR siginvu, realnum dt, const RPR cnd,
               const RPR cndinv, RPR fcnd) {
  (void)c;   // currently unused
  if (!g1) { // swap g1 and g2
    SWAP(const RPR, g1, g2);
    SWAP(ptrdiff_t, s1, s2);
    dtdx = -dtdx; // need to flip derivative sign
  }

  /* The following are a bunch of special cases of the "MOST GENERAL CASE"
     loop below.  We make copies of the loop for each special case in
     order to keep the innermost loop efficient.  This is especially
     important because the non-PML cases are actually more common.
     (The "right" way to do this is by partial evaluation of the
      most general case, but that would require a code generator.) */

  if (dsig == NO_DIRECTION) {    // no PML in f update
    if (dsigu == NO_DIRECTION) { // no fu update
      if (cnd) {
        realnum dt2 = dt * 0.5;
        if (g2) {
          PLOOP_OVER_IVECS(gv, is, ie, i) {
            f[i] = ((1 - dt2 * cnd[i]) * f[i] - dtdx * (g1[i + s1] - g1[i] + g2[i] - g2[i + s2])) *
                   cndinv[i];
          }
        }
        else {
          PLOOP_OVER_IVECS(gv, is, ie, i) {
            f[i] = ((1 - dt2 * cnd[i]) * f[i] - dtdx * (g1[i + s1] - g1[i])) * cndinv[i];
          }
        }
      }
      else { // no conductivity
        if (g2) {
          PLOOP_OVER_IVECS(gv, is, ie, i) {
            f[i] -= dtdx * (g1[i + s1] - g1[i] + g2[i] - g2[i + s2]);
          }
        }
        else {
          PLOOP_OVER_IVECS(gv, is, ie, i) { f[i] -= dtdx * (g1[i + s1] - g1[i]); }
        }
      }
    }
    else { // fu update, no PML in f update
      KSTRIDE_DEF(dsigu, ku, is, gv);
      if (cnd) {
        realnum dt2 = dt * 0.5;
        if (g2) {
          PLOOP_OVER_IVECS(gv, is, ie, i) {
            DEF_ku;
            realnum fprev = fu[i];
            fu[i] =
                ((1 - dt2 * cnd[i]) * fprev - dtdx * (g1[i + s1] - g1[i] + g2[i] - g2[i + s2])) *
                cndinv[i];
            f[i] = siginvu[ku] * ((kapu[ku] - sigu[ku]) * f[i] + fu[i] - fprev);
          }
        }
        else {
          PLOOP_OVER_IVECS(gv, is, ie, i) {
            DEF_ku;
            realnum fprev = fu[i];
            fu[i] = ((1 - dt2 * cnd[i]) * fprev - dtdx * (g1[i + s1] - g1[i])) * cndinv[i];
            f[i] = siginvu[ku] * ((kapu[ku] - sigu[ku]) * f[i] + fu[i] - fprev);
          }
        }
      }
      else { // no conductivity
        if (g2) {
          PLOOP_OVER_IVECS(gv, is, ie, i) {
            DEF_ku;
            realnum fprev = fu[i];
            fu[i] -= dtdx * (g1[i + s1] - g1[i] + g2[i] - g2[i + s2]);
            f[i] = siginvu[ku] * ((kapu[ku] - sigu[ku]) * f[i] + fu[i] - fprev);
          }
        }
        else {
          PLOOP_OVER_IVECS(gv, is, ie, i) {
            DEF_ku;
            realnum fprev = fu[i];
            fu[i] -= dtdx * (g1[i + s1] - g1[i]);
            f[i] = siginvu[ku] * ((kapu[ku] - sigu[ku]) * f[i] + fu[i] - fprev);
          }
        }
      }
    }
  }
  else { /* PML in f update */
    KSTRIDE_DEF(dsig, k, is, gv);
    if (dsigu == NO_DIRECTION) { // no fu update
      if (cnd) {
        realnum dt2 = dt * 0.5;
        if (g2) {
          PLOOP_OVER_IVECS(gv, is, ie, i) {
            DEF_k;
            realnum fcnd_prev = fcnd[i];
            fcnd[i] =
                ((1 - dt2 * cnd[i]) * fcnd[i] - dtdx * (g1[i + s1] - g1[i] + g2[i] - g2[i + s2])) *
                cndinv[i];
            f[i] = ((kap[k] - sig[k]) * f[i] + (fcnd[i] - fcnd_prev)) * siginv[k];
          }
        }
        else {
          PLOOP_OVER_IVECS(gv, is, ie, i) {
            DEF_k;
            realnum fcnd_prev = fcnd[i];
            fcnd[i] = ((1 - dt2 * cnd[i]) * fcnd[i] - dtdx * (g1[i + s1] - g1[i])) * cndinv[i];
            f[i] = ((kap[k] - sig[k]) * f[i] + (fcnd[i] - fcnd_prev)) * siginv[k];
          }
        }
      }
      else { // no conductivity (other than PML conductivity)
        if (g2) {
          PLOOP_OVER_IVECS(gv, is, ie, i) {
            DEF_k;
            f[i] = ((kap[k] - sig[k]) * f[i] - dtdx * (g1[i + s1] - g1[i] + g2[i] - g2[i + s2])) *
                   siginv[k];
          }
        }
        else {
          PLOOP_OVER_IVECS(gv, is, ie, i) {
            DEF_k;
            f[i] = ((kap[k] - sig[k]) * f[i] - dtdx * (g1[i + s1] - g1[i])) * siginv[k];
          }
        }
      }
    }
    else { // fu update + PML in f update
      KSTRIDE_DEF(dsigu, ku, is, gv);
      if (cnd) {
        realnum dt2 = dt * 0.5;
        if (g2) {
          //////////////////// MOST GENERAL CASE //////////////////////
          PLOOP_OVER_IVECS(gv, is, ie, i) {
            DEF_k;
            DEF_ku;
            realnum fprev = fu[i];
            realnum fcnd_prev = fcnd[i];
            fcnd[i] =
                ((1 - dt2 * cnd[i]) * fcnd[i] - dtdx * (g1[i + s1] - g1[i] + g2[i] - g2[i + s2])) *
                cndinv[i];
            fu[i] = ((kap[k] - sig[k]) * fu[i] + (fcnd[i] - fcnd_prev)) * siginv[k];
            f[i] = siginvu[ku] * ((kapu[ku] - sigu[ku]) * f[i] + fu[i] - fprev);
          }
          /////////////////////////////////////////////////////////////
        }
        else {
          PLOOP_OVER_IVECS(gv, is, ie, i) {
            DEF_k;
            DEF_ku;
            realnum fprev = fu[i];
            realnum fcnd_prev = fcnd[i];
            fcnd[i] = ((1 - dt2 * cnd[i]) * fcnd[i] - dtdx * (g1[i + s1] - g1[i])) * cndinv[i];
            fu[i] = ((kap[k] - sig[k]) * fu[i] + (fcnd[i] - fcnd_prev)) * siginv[k];
            f[i] = siginvu[ku] * ((kapu[ku] - sigu[ku]) * f[i] + fu[i] - fprev);
          }
        }
      }
      else { // no conductivity (other than PML conductivity)
        if (g2) {
          PLOOP_OVER_IVECS(gv, is, ie, i) {
            DEF_k;
            DEF_ku;
            realnum fprev = fu[i];
            fu[i] = ((kap[k] - sig[k]) * fu[i] - dtdx * (g1[i + s1] - g1[i] + g2[i] - g2[i + s2])) *
                    siginv[k];
            f[i] = siginvu[ku] * ((kapu[ku] - sigu[ku]) * f[i] + fu[i] - fprev);
          }
        }
        else {
          PLOOP_OVER_IVECS(gv, is, ie, i) {
            DEF_k;
            DEF_ku;
            realnum fprev = fu[i];
            fu[i] = ((kap[k] - sig[k]) * fu[i] - dtdx * (g1[i + s1] - g1[i])) * siginv[k];
            f[i] = siginvu[ku] * ((kapu[ku] - sigu[ku]) * f[i] + fu[i] - fprev);
          }
        }
      }
    }
  }
}

/* field-update equation f += betadt * g (plus variants for conductivity
   and/or PML).  This is used in 2d calculations to add an exp(i beta z)
   time dependence, which gives an additional i \beta \hat{z} \times
   cross-product in the curl equations. */
void step_beta(RPR f, component c, const RPR g, const grid_volume &gv, const ivec is, const ivec ie,
               realnum betadt, direction dsig, const RPR siginv, RPR fu, direction dsigu,
               const RPR siginvu, const RPR cndinv, RPR fcnd) {
  (void)c; // currently unused
  if (!g) return;
  if (dsig != NO_DIRECTION) { // PML in f update
    KSTRIDE_DEF(dsig, k, is, gv);
    if (dsigu != NO_DIRECTION) { // PML in f + fu
      KSTRIDE_DEF(dsigu, ku, is, gv);
      if (cndinv) { // conductivity + PML
        //////////////////// MOST GENERAL CASE //////////////////////
        PLOOP_OVER_IVECS(gv, is, ie, i) {
          DEF_k;
          DEF_ku;
          realnum df;
          realnum dfcnd = betadt * g[i] * cndinv[i];
          fcnd[i] += dfcnd;
          fu[i] += (df = dfcnd * siginv[k]);
          f[i] += siginvu[ku] * df;
        }
        /////////////////////////////////////////////////////////////
      }
      else { // PML only
        PLOOP_OVER_IVECS(gv, is, ie, i) {
          DEF_k;
          DEF_ku;
          realnum df;
          fu[i] += (df = betadt * g[i] * siginv[k]);
          f[i] += siginvu[ku] * df;
        }
      }
    }
    else {          // PML in f, no fu
      if (cndinv) { // conductivity + PML
        PLOOP_OVER_IVECS(gv, is, ie, i) {
          DEF_k;
          realnum dfcnd = betadt * g[i] * cndinv[i];
          fcnd[i] += dfcnd;
          f[i] += dfcnd * siginv[k];
        }
      }
      else { // PML only
        PLOOP_OVER_IVECS(gv, is, ie, i) {
          DEF_k;
          f[i] += betadt * g[i] * siginv[k];
        }
      }
    }
  }
  else {                         // no PML in f update
    if (dsigu != NO_DIRECTION) { // fu, no PML in f
      KSTRIDE_DEF(dsigu, ku, is, gv);
      if (cndinv) { // conductivity, no PML
        PLOOP_OVER_IVECS(gv, is, ie, i) {
          DEF_ku;
          realnum df;
          fu[i] += (df = betadt * g[i] * cndinv[i]);
          f[i] += siginvu[ku] * df;
        }
      }
      else { // no conductivity or PML
        PLOOP_OVER_IVECS(gv, is, ie, i) {
          DEF_ku;
          realnum df;
          fu[i] += (df = betadt * g[i]);
          f[i] += siginvu[ku] * df;
        }
      }
    }
    else {          // no PML, no fu
      if (cndinv) { // conductivity, no PML
        PLOOP_OVER_IVECS(gv, is, ie, i) { f[i] += betadt * g[i] * cndinv[i]; }
      }
      else { // no conductivity or PML
        PLOOP_OVER_IVECS(gv, is, ie, i) { f[i] += betadt * g[i]; }
      }
    }
  }
}
// allows fixed angle broadband simulations
void step_bfast(RPR f, component c, const RPR g1, const RPR g2, ptrdiff_t s1,
                ptrdiff_t s2, // strides for g1/g2 shift
                const grid_volume &gv, const ivec is, const ivec ie, realnum dtdx, direction dsig,
                const RPR sig, const RPR kap, const RPR siginv, RPR fu, direction dsigu,
                const RPR sigu, const RPR kapu, const RPR siginvu, realnum dt, const RPR cnd,
                const RPR cndinv, RPR fcnd, RPR F, realnum k1, realnum k2) {
  (void)c;   // currently unused
  if (!g1) { // swap g1 and g2
    SWAP(const RPR, g1, g2);
    SWAP(ptrdiff_t, s1, s2);
    SWAP(realnum, k1, k2); // need to swap in cross product
  }
  if (dsig == NO_DIRECTION) {    // no PML in f update
    if (dsigu == NO_DIRECTION) { // no fu update
      if (cnd) {
        if (g2) {
          PLOOP_OVER_IVECS(gv, is, ie, i) {
            realnum F_prev = F[i];
            F[i] = (k1 * (g1[i + s1] + g1[i]) - k2 * (g2[i + s2] + g2[i])) - F[i];
            f[i] += (F[i] - F_prev) * cndinv[i];
          }
        }
        else {
          PLOOP_OVER_IVECS(gv, is, ie, i) {
            realnum F_prev = F[i];
            F[i] = k1 * (g1[i + s1] + g1[i]) - F[i];
            f[i] += (F[i] - F_prev) * cndinv[i];
          }
        }
      }
      else { // no conductivity
        if (g2) {
          PLOOP_OVER_IVECS(gv, is, ie, i) {
            realnum F_prev = F[i];
            F[i] = (k1 * (g1[i + s1] + g1[i]) - k2 * (g2[i + s2] + g2[i])) - F[i];
            f[i] += (F[i] - F_prev);
          }
        }
        else {
          PLOOP_OVER_IVECS(gv, is, ie, i) {
            realnum F_prev = F[i];
            F[i] = k1 * (g1[i + s1] + g1[i]);
            f[i] += (F[i] - F_prev);
          }
        }
      }
    }
    else { // fu update, no PML in f update
      KSTRIDE_DEF(dsigu, ku, is, gv);
      if (cnd) {
        if (g2) {
          PLOOP_OVER_IVECS(gv, is, ie, i) {
            DEF_ku;
            realnum df;
            realnum F_prev = F[i];
            F[i] = (k1 * (g1[i + s1] + g1[i]) - k2 * (g2[i + s2] + g2[i])) - F[i];
            fu[i] += (df = (F[i] - F_prev) * cndinv[i]);
            f[i] += siginvu[ku] * df;
          }
        }
        else {
          PLOOP_OVER_IVECS(gv, is, ie, i) {
            DEF_ku;
            realnum df;
            realnum F_prev = F[i];
            F[i] = k1 * (g1[i + s1] + g1[i]) - F[i];
            fu[i] += (df = (F[i] - F_prev) * cndinv[i]);
            f[i] += siginvu[ku] * df;
          }
        }
      }
      else { // no conductivity
        if (g2) {
          PLOOP_OVER_IVECS(gv, is, ie, i) {
            DEF_ku;
            realnum df;
            realnum F_prev = F[i];
            F[i] = (k1 * (g1[i + s1] + g1[i]) - k2 * (g2[i + s2] + g2[i])) - F[i];
            fu[i] += (df = (F[i] - F_prev));
            f[i] += siginvu[ku] * df;
          }
        }
        else {
          PLOOP_OVER_IVECS(gv, is, ie, i) {
            DEF_ku;
            realnum df;
            realnum F_prev = F[i];
            F[i] = k1 * (g1[i + s1] + g1[i]) - F[i];
            fu[i] += (df = (F[i] - F_prev));
            f[i] += siginvu[ku] * df;
          }
        }
      }
    }
  }
  else { // PML in f update
    KSTRIDE_DEF(dsig, k, is, gv);
    if (dsigu == NO_DIRECTION) { // no fu update
      if (cnd) {
        realnum dt2 = dt * 0.5;
        if (g2) {
          PLOOP_OVER_IVECS(gv, is, ie, i) {
            DEF_k;
            realnum F_prev = F[i];
            F[i] = (k1 * (g1[i + s1] + g1[i]) - k2 * (g2[i + s2] + g2[i])) - F[i];
            realnum dfcnd = (F[i] - F_prev) * cndinv[i];
            fcnd[i] += dfcnd;
            f[i] += dfcnd * siginv[k];
          }
        }
        else {
          PLOOP_OVER_IVECS(gv, is, ie, i) {
            DEF_k;
            realnum F_prev = F[i];
            F[i] = k1 * (g1[i + s1] + g1[i]) - F[i];
            realnum dfcnd = (F[i] - F_prev) * cndinv[i];
            fcnd[i] += dfcnd;
            f[i] += dfcnd * siginv[k];
          }
        }
      }
      else { // no conductivity (other than PML conductivity)
        if (g2) {
          PLOOP_OVER_IVECS(gv, is, ie, i) {
            DEF_k;
            realnum F_prev = F[i];
            F[i] = (k1 * (g1[i + s1] + g1[i]) - k2 * (g2[i + s2] + g2[i])) - F[i];
            f[i] += (F[i] - F_prev) * siginv[k];
          }
        }
        else {
          PLOOP_OVER_IVECS(gv, is, ie, i) {
            DEF_k;
            realnum F_prev = F[i];
            F[i] = k1 * (g1[i + s1] + g1[i]) - F[i];
            f[i] += (F[i] - F_prev) * siginv[k];
          }
        }
      }
    }
    else { // fu update + PML in f update
      KSTRIDE_DEF(dsigu, ku, is, gv);
      if (cnd) {
        if (g2) {
          //////////////////// MOST GENERAL CASE //////////////////////
          PLOOP_OVER_IVECS(gv, is, ie, i) {
            DEF_k;
            DEF_ku;
            realnum df;
            realnum F_prev = F[i];
            F[i] = (k1 * (g1[i + s1] + g1[i]) - k2 * (g2[i + s2] + g2[i])) - F[i];
            realnum dfcnd = (F[i] - F_prev) * cndinv[i];
            fcnd[i] += dfcnd;
            fu[i] += (df = dfcnd * siginv[k]);
            f[i] += siginvu[ku] * df;
          }
          /////////////////////////////////////////////////////////////
        }
        else {
          PLOOP_OVER_IVECS(gv, is, ie, i) {
            DEF_k;
            DEF_ku;
            realnum df;
            realnum F_prev = F[i];
            F[i] = k1 * (g1[i + s1] + g1[i]) - F[i];
            realnum dfcnd = (F[i] - F_prev) * cndinv[i];
            fcnd[i] += dfcnd;
            fu[i] += (df = dfcnd * siginv[k]);
            f[i] += siginvu[ku] * df;
          }
        }
      }
      else { // no conductivity (other than PML conductivity)
        if (g2) {
          PLOOP_OVER_IVECS(gv, is, ie, i) {
            DEF_k;
            DEF_ku;
            realnum df;
            realnum F_prev = F[i];
            F[i] = (k1 * (g1[i + s1] + g1[i]) - k2 * (g2[i + s2] + g2[i])) - F[i];
            fu[i] += (df = (F[i] - F_prev) * siginv[k]);
            f[i] += siginvu[ku] * df;
          }
        }
        else {
          PLOOP_OVER_IVECS(gv, is, ie, i) {
            DEF_k;
            DEF_ku;
            realnum df;
            realnum F_prev = F[i];
            F[i] = k1 * (g1[i + s1] + g1[i]) - F[i];
            fu[i] += (df = (F[i] - F_prev) * siginv[k]);
            f[i] += siginvu[ku] * df;
          }
        }
      }
    }
  }
}

/* Given Dsqr = |D|^2 and Di = component of D, compute the factor f so
   that Ei = chi1inv * f * Di.   In principle, this would involve solving
   a cubic equation, but instead we use a Pade approximant that is
   accurate to several orders.  This is inaccurate if the nonlinear
   index change is large, of course, but in that case the chi2/chi3
   power-series expansion isn't accurate anyway, so the cubic isn't
   physical there either. */
inline realnum calc_nonlinear_u(const realnum Dsqr, const realnum Di, const realnum chi1inv,
                                const realnum chi2, const realnum chi3) {
  realnum c2 = Di * chi2 * (chi1inv * chi1inv);
  realnum c3 = Dsqr * chi3 * (chi1inv * chi1inv * chi1inv);
  return (1 + c2 + 2 * c3) / (1 + 2 * c2 + 3 * c3);
}

/* Update E from D using epsilon and PML, *or* update H from B using
   mu and PML.

   To be generic, here we set f = u * g, where u may
   be a tensor, and we also have a nonlinear susceptibility chi.
   Here, g = (g,g1,g2) where g1 and g2 are the off-diagonal
   components, if any (g2 may be NULL).

   In PML (dsigw != NO_DIR), we have an additional auxiliary field fw,
   which is updated by the equations:
          fw = u * g
          df/dt = kappaw dfw/dt - sigmaw * fw
   That is, fw is updated like the non-PML f, and f is updated from
   fw by a little ODE.  Here, sigw[k] = sigmaw[k]*dt/2, kappaw[k] = kapw[k]

*/

/// TODO add old macro and fn defs for this
void step_update_EDHB(RPR f, component fc, const grid_volume &gv, const ivec is, const ivec ie,
                      const RPR g, const RPR g1, const RPR g2, const RPR u, const RPR u1,
                      const RPR u2, ptrdiff_t s, ptrdiff_t s1, ptrdiff_t s2, const RPR chi2,
                      const RPR chi3, RPR fw, direction dsigw, const RPR sigw, const RPR kapw) {
  (void)fc; // currently unused
  if (!f) return;
  // f= E field,  fc = field component (ec).
  // g, g1 g2 are the D components (Dx, Dy..) of the fminusp data
  //  u, u1, u2 are the E components of ' s->chi1inv' (whatever that is...)
  // s, s1, s2 are the strides in the relevant directions
  // chi2 is s->chi2[ec], suggesting only has components for the three E directions exclusivley (but
  // need to check what chi2[] actually looks like) fw is W auxiliary field. Last three are sigma
  // components.

  if ((!g1 && g2) || (g1 && g2 && !u1 && u2)) { /* swap g1 and g2 */
    SWAP(const RPR, g1, g2);
    SWAP(const RPR, u1, u2);
    SWAP(ptrdiff_t, s1, s2);
  }
///cout << "in linear loop  "  << endl;
// stable averaging of offdiagonal components
#define OFFDIAG(u, g, sx)                                                                          \
  (0.25 * ((g[i] + g[i - sx]) * u[i] + (g[i + s] + g[(i + s) - sx]) * u[i + s]))

  /* As with step_curl, these loops are all essentially copies
     of the "MOST GENERAL CASE" loop with various terms thrown out. */
  if (chi3) { cout << "IN CHI3!! " << endl;
  }

  if (dsigw != NO_DIRECTION) { //////// PML case (with fw) /////////////
  ///cout << "in linear loop PML" << endl;
    KSTRIDE_DEF(dsigw, kw, is, gv);
    if (u1 && u2) { // 3x3 off-diagonal u
      if (chi3) {
        //////////////////// MOST GENERAL CASE //////////////////////
        PLOOP_OVER_IVECS(gv, is, ie, i) {
          realnum g1s = g1[i] + g1[i + s] + g1[i - s1] + g1[i + (s - s1)];
          realnum g2s = g2[i] + g2[i + s] + g2[i - s2] + g2[i + (s - s2)];
          realnum gs = g[i];
          realnum us = u[i];
          DEF_kw;
          realnum fwprev = fw[i], kapwkw = kapw[kw], sigwkw = sigw[kw];
          fw[i] = (gs * us + OFFDIAG(u1, g1, s1) + OFFDIAG(u2, g2, s2)) *
                  calc_nonlinear_u(gs * gs + 0.0625 * (g1s * g1s + g2s * g2s), gs, us, chi2[i],
                                   chi3[i]);
          f[i] += (kapwkw + sigwkw) * fw[i] - (kapwkw - sigwkw) * fwprev;
        }
        /////////////////////////////////////////////////////////////
      }
      else {
        PLOOP_OVER_IVECS(gv, is, ie, i) {
          realnum gs = g[i];
          realnum us = u[i];
          DEF_kw;
          realnum fwprev = fw[i], kapwkw = kapw[kw], sigwkw = sigw[kw];
          fw[i] = (gs * us + OFFDIAG(u1, g1, s1) + OFFDIAG(u2, g2, s2));
          f[i] += (kapwkw + sigwkw) * fw[i] - (kapwkw - sigwkw) * fwprev;
        }
      }
    }
    else if (u1) { // 2x2 off-diagonal u
      if (chi3) {
        PLOOP_OVER_IVECS(gv, is, ie, i) {
          realnum g1s = g1[i] + g1[i + s] + g1[i - s1] + g1[i + (s - s1)];
          realnum gs = g[i];
          realnum us = u[i];
          DEF_kw;
          realnum fwprev = fw[i], kapwkw = kapw[kw], sigwkw = sigw[kw];
          fw[i] = (gs * us + OFFDIAG(u1, g1, s1)) *
                  calc_nonlinear_u(gs * gs + 0.0625 * (g1s * g1s), gs, us, chi2[i], chi3[i]);
          f[i] += (kapwkw + sigwkw) * fw[i] - (kapwkw - sigwkw) * fwprev;
        }
      }
      else {
        PLOOP_OVER_IVECS(gv, is, ie, i) {
          realnum gs = g[i];
          realnum us = u[i];
          DEF_kw;
          realnum fwprev = fw[i], kapwkw = kapw[kw], sigwkw = sigw[kw];
          fw[i] = (gs * us + OFFDIAG(u1, g1, s1));
          f[i] += (kapwkw + sigwkw) * fw[i] - (kapwkw - sigwkw) * fwprev;
        }
      }
    }
    else if (u2) { // 2x2 off-diagonal u
      meep::abort("bug - didn't swap off-diagonal terms!?");
    }
    else { // diagonal u
      if (chi3) {
        if (g1 && g2) {
          PLOOP_OVER_IVECS(gv, is, ie, i) {
            realnum g1s = g1[i] + g1[i + s] + g1[i - s1] + g1[i + (s - s1)];
            realnum g2s = g2[i] + g2[i + s] + g2[i - s2] + g2[i + (s - s2)];
            realnum gs = g[i];
            realnum us = u[i];
            DEF_kw;
            realnum fwprev = fw[i], kapwkw = kapw[kw], sigwkw = sigw[kw];
            fw[i] = (gs * us) * calc_nonlinear_u(gs * gs + 0.0625 * (g1s * g1s + g2s * g2s), gs, us,
                                                 chi2[i], chi3[i]);
            f[i] += (kapwkw + sigwkw) * fw[i] - (kapwkw - sigwkw) * fwprev;
          }
        }
        else if (g1) {
          PLOOP_OVER_IVECS(gv, is, ie, i) {
            realnum g1s = g1[i] + g1[i + s] + g1[i - s1] + g1[i + (s - s1)];
            realnum gs = g[i];
            realnum us = u[i];
            DEF_kw;
            realnum fwprev = fw[i], kapwkw = kapw[kw], sigwkw = sigw[kw];
            fw[i] = (gs * us) *
                    calc_nonlinear_u(gs * gs + 0.0625 * (g1s * g1s), gs, us, chi2[i], chi3[i]);
            f[i] += (kapwkw + sigwkw) * fw[i] - (kapwkw - sigwkw) * fwprev;
          }
        }
        else if (g2) { meep::abort("bug - didn't swap off-diagonal terms!?"); }
        else {
          PLOOP_OVER_IVECS(gv, is, ie, i) {
            realnum gs = g[i];
            realnum us = u[i];
            DEF_kw;
            realnum fwprev = fw[i], kapwkw = kapw[kw], sigwkw = sigw[kw];
            fw[i] = (gs * us) * calc_nonlinear_u(gs * gs, gs, us, chi2[i], chi3[i]);
            f[i] += (kapwkw + sigwkw) * fw[i] - (kapwkw - sigwkw) * fwprev;
          }
        }
      }
      else if (u) {
        PLOOP_OVER_IVECS(gv, is, ie, i) {
          realnum gs = g[i];
          realnum us = u[i];
          DEF_kw;
          realnum fwprev = fw[i], kapwkw = kapw[kw], sigwkw = sigw[kw];
          fw[i] = (gs * us);
          f[i] += (kapwkw + sigwkw) * fw[i] - (kapwkw - sigwkw) * fwprev;
        }
      }
      else {
        PLOOP_OVER_IVECS(gv, is, ie, i) {
          DEF_kw;
          realnum fwprev = fw[i], kapwkw = kapw[kw], sigwkw = sigw[kw];
          fw[i] = g[i];
          f[i] += (kapwkw + sigwkw) * fw[i] - (kapwkw - sigwkw) * fwprev;
        }
      }
    }
  }
  else {            /////////////// no PML (no fw) ///////////////////
  ///cout << "in linear loop Non pml" << endl;

    if (u1 && u2) { // 3x3 off-diagonal u
      ///cout << "at chi3" << endl;
      if (chi3) {
        PLOOP_OVER_IVECS(gv, is, ie, i) {
          realnum g1s = g1[i] + g1[i + s] + g1[i - s1] + g1[i + (s - s1)];
          realnum g2s = g2[i] + g2[i + s] + g2[i - s2] + g2[i + (s - s2)];
          realnum gs = g[i];
          realnum us = u[i];
          f[i] = (gs * us + OFFDIAG(u1, g1, s1) + OFFDIAG(u2, g2, s2)) *
                 calc_nonlinear_u(gs * gs + 0.0625 * (g1s * g1s + g2s * g2s), gs, us, chi2[i],
                                  chi3[i]);
        }
      }
      

      else {
        ///cout << "at chi3 else" << endl;
        PLOOP_OVER_IVECS(gv, is, ie, i) {
          realnum gs = g[i];
          realnum us = u[i];
          f[i] = (gs * us + OFFDIAG(u1, g1, s1) + OFFDIAG(u2, g2, s2));
        }
      }
    }
    else if (u1) { // 2x2 off-diagonal u
      ///cout << "at chi3 u1" << endl;

      if (chi3) {
        PLOOP_OVER_IVECS(gv, is, ie, i) {
          realnum g1s = g1[i] + g1[i + s] + g1[i - s1] + g1[i + (s - s1)];
          realnum gs = g[i];
          realnum us = u[i];
          f[i] = (gs * us + OFFDIAG(u1, g1, s1)) *
                 calc_nonlinear_u(gs * gs + 0.0625 * (g1s * g1s), gs, us, chi2[i], chi3[i]);
        }
      }
      

      else {
        ///cout << "at chi3 u1 else" << endl;
        PLOOP_OVER_IVECS(gv, is, ie, i) {
          realnum gs = g[i];
          realnum us = u[i];
          f[i] = (gs * us + OFFDIAG(u1, g1, s1));
        }
      }
    }
    else if (u2) { // 2x2 off-diagonal u
      ///cout << "abrte" << endl;

      meep::abort("bug - didn't swap off-diagonal terms!?");
    }
    else { // diagonal u
      ///cout << "at chi3 diag " << endl;

      if (chi3) {
        ///cout << "rzjrejz " << endl;

        if (g1 && g2) {
          PLOOP_OVER_IVECS(gv, is, ie, i) {
            realnum g1s = g1[i] + g1[i + s] + g1[i - s1] + g1[i + (s - s1)];
            realnum g2s = g2[i] + g2[i + s] + g2[i - s2] + g2[i + (s - s2)];
            realnum gs = g[i];
            realnum us = u[i];
            f[i] = (gs * us) * calc_nonlinear_u(gs * gs + 0.0625 * (g1s * g1s + g2s * g2s), gs, us,
                                                chi2[i], chi3[i]);
          }
        }
        else if (g1) {
          ///cout << "ahrewjr " << endl;

          PLOOP_OVER_IVECS(gv, is, ie, i) {
            realnum g1s = g1[i] + g1[i + s] + g1[i - s1] + g1[i + (s - s1)];
            realnum gs = g[i];
            realnum us = u[i];
            f[i] = (gs * us) *
                   calc_nonlinear_u(gs * gs + 0.0625 * (g1s * g1s), gs, us, chi2[i], chi3[i]);
          }
        }
        else if (g2) {
          ///cout << "ykdrs " << endl;

          meep::abort("bug - didn't swap off-diagonal terms!?");
        }

        else {
        ///cout << "lydtser " << endl;
          PLOOP_OVER_IVECS(gv, is, ie, i) {
            realnum gs = g[i];
            realnum us = u[i];
            f[i] = (gs * us) * calc_nonlinear_u(gs * gs, gs, us, chi2[i], chi3[i]);
          }
        }
      }

      else if (u) {
        ///cout << "most basic case  " << endl;
        PLOOP_OVER_IVECS(gv, is, ie, i) {
       //   ///cout << u[i] << "  " << g[i] << endl;

          realnum gs = g[i];
          realnum us = u[i];
          f[i] = (gs * us);
        }
      }
      else {
        cout << " SHOULDN'T be here!" << endl;
        PLOOP_OVER_IVECS(gv, is, ie, i)
        {
        f[i] = g[i];
      }
    }
    }
  }
}


/// TODO rename macro and fn defs for this
void step_update_EDHB_NL(RPR f, RPR f_2, RPR f_3, component fc, const grid_volume &gv, const ivec is, const ivec is_2, const ivec is_3, const ivec ie, const RPR g, const RPR g1, const RPR g2, const RPR u, const RPR u_2, const RPR u_3, const RPR u1, const RPR u2, ptrdiff_t s, ptrdiff_t s1, ptrdiff_t s2, const realnum* chi2new, const RPR chi3, RPR fw, RPR fw_2_atZ, RPR fw_3_atZ, RPR fw_2, RPR fw_3, direction dsigw, direction dsigw_2, direction dsigw_3, const RPR sigw, const RPR sigw_2, const RPR sigw_3, const RPR kapw, const RPR kapw_2, const RPR kapw_3) {
  (void)fc; // currently unused
  if (!f) {
  ///cout << "returning for bad reasons" << endl;
      return; }

  ///cout << "NL ..wafewfaef" << endl;

// stable averaging of offdiagonal components
#define OFFDIAG(u, g, sx)                                                                          \
  (0.25 * ((g[i] + g[i - sx]) * u[i] + (g[i + s] + g[(i + s) - sx]) * u[i + s]))

  if (dsigw != NO_DIRECTION) { //////// PML case (with fw) /////////////  // TODO need to also
                               ///implement for non pml case, since this fn operates chunk-wise, and
                               ///some chunks won't have pml even if some do.
    KSTRIDE_DEF(dsigw, kw, is, gv); /// Used in DEF_kw.  dsigw, is, and gv come from fn. kw0,
                                    /// skw1/2/3 out. TODO check wrt position
  ///cout << "in Pml NL case 1 " << endl;
    // TODO implement off-diagonal epsilon terms into NewtonRaphson solver
  ///  if (chi3) { 
      /// Build NR solver into this section here. **should add similar copy to the section further
      /// down without PML... Will be callable by defining the 2nd order NL material with a full 3x3
      /// offdiagonal epsilon as a flag (i.e. u1 && u2 must be non-zero, however they will NOT be
      /// used in the calc) and also setting chi3 as nonzero as a flag (chi3 will not be used, but
      /// setting it to non-zero is a hack to make the code enter this 'if' statement for the NL
      /// material
      ///
      /// f, f_2 and f_3 (i.e, Z, X, and Y respectively) are (in non-PML case) i-th cell electric
      /// field the output from this fn. In PML case they are updated as per original meep (but not
      /// clear what exactly the fields represent in this case)... fc = field cmpnt gv = grid volume
      /// (need to ploopoverivecs) is, ie, start and end of chunk. is, is_2, is_3 are start ivecs
      /// for Z, X, and Y directions respectively g, g1, g2 = 'D - P'  field components z, x, y u is
      /// inverse epsilon Z direction. u1, u2 = offdiagonal inverse epsilon(NOT used - except as
      /// flag). u_2 and u_3 are inverse epsilon in the X and Y directions s, s1, s2 = strides (for
      /// getting avged orthogonal fields over adjacent cells (see yeecell diag it makes sense). For
      /// this NL fn, these will have been set in the most recent 'FOR_FT_COMPONENTS' loop, in which
      /// 'ec' would have been Ez, therefore s is Ez stride, s1 and s2 are x and y respectively...
      /// TODO - DOUBLE CHECK THIS IS TRUE from update_eh! chi2 = actual simulation value of chi2
      /// (in meep units..?), will be passed into NR chi3 = a flag (set to nonzero float to indicate
      /// NL material) fw[i] = PML case: i-th cell Ez field the output from this fn. fw_2_atZ and
      /// fw_3_atZ are X and Y E fields AT Z LOCATIONS (see below). fw_2/3 are for the final
      /// interpolated Ex and Ey fields. dsigw =           /// sigw =           /// kapw =

      ///  fw_2_atZ and fw_3_atZ need to have the same dim as fw. If fw[i] is i'th Ex field, fw_2/3
      ///  are the Ex and Ey fields at the SAME location as the Ez field. According
      /// to the Yee cell, the Ex and Ey fields are not in fact at the same location as Ez or each
      /// other. The Ex and Ey fields at the correct yee cell locations must therefore subsequently
      /// calculated by interpolation (same principle as for gs_2 below), in subsequent
      /// ploopoverivecs

      PLOOP_OVER_IVECS(gv, is, ie, i) {

          if (chi2new[i] == 0.0) { 
            ///cout << "Continuing" << endl;
              continue; 
          }
          cout << "IN PML CASE!"<<endl; // in theory it should never enter this case anyway, since the PML should not be located within the chi2 material...
          sleep(2);
        realnum gs = g[i]; // dmpZ
        // avg orthogonal D-P fields over adjacent cells (see yee cell diag to understand why...)
        realnum gs_2 = (g1[i] + g1[i + s] + g1[i - s1] + g1[i + (s - s1)]) * 0.25; // dmpX
        realnum gs_3 = (g2[i] + g2[i + s] + g2[i - s2] + g2[i + (s - s2)]) * 0.25; // dmpY

        /// taking inverse of chi1inverse is easiest way to access epsilon...
        realnum us = 1 / u[i];
        realnum us_2 = 1 / u_2[i];
        realnum us_3 = 1 / u_3[i];

       // if ( (us - 11.5)*(us-10.5) >0){}
       // cout << "us " << us << " us2 " << us_2 << endl;


        /// #will be format Parameters p1 = {prevF D-P_X, eps, 0, 0, 0, chi2, 0, 0 } etc;
        Parameters p1 = {gs_2, us_2, 0.0, 0.0, 0.0, chi2new[i], 0.0, 0.0}; // X
        Parameters p2 = {gs_3, us_3, 0.0, 0.0, 0.0, 0.0, chi2new[i], 0.0}; // Y
        Parameters p3 = {gs,  us,  0.0,       0.0, 0.0,
                         0.0, 0.0, chi2new[i]}; // Z. currently using all chi2 tensor components
                                                // equal (as zinc blende)

        realnum seed1 = fw[i]; // Z
        realnum seed2 = fw_2_atZ[i]; // X
        realnum seed3 = fw_3_atZ[i]; // Y
        realnum fwprev = fw[i];
  ///cout << "in Pml NL case 1 at NR " << endl;

        /// Newton Raphson for calculating Ez, Ex and Ey fields, (AT Z LOCATIONS):
        ///  Seeded with previous field vals. Passing in field array pointers to be assigned new
        ///  vals.
        runNR(seed2, seed3, seed1, &fw_2_atZ[i], &fw_3_atZ[i], &fw[i], p1, p2, p3);
  ///cout << "in Pml NL case 1 done NR" << endl;

        // Do the other fields for PML (whatever they do exactly..)
        DEF_kw; 
        /// uses kw0, skw1/2/3. not sure what it does exactly...
       realnum kapwkw = kapw[kw], sigwkw = sigw[kw]; // Ez
        f[i] += (kapwkw + sigwkw) * fw[i] - (kapwkw - sigwkw) * fwprev;
      }

      // now do the other two PLOOPs to interpolate the X and Y fields to their correct positions,
      // and then calculate f_2 and f_3..
      KSTRIDE_DEF(dsigw_2, kw_2, is_2, gv);
      PLOOP_OVER_IVECS(gv, is_2, ie, i) { /// Round two for interpolating X

        if (chi2new[i] == 0.0) { continue; }// TODO should this be in these two interpolation loops??

        DEF_kw_2; // these will throw if the necessary values are null, but they shouldn't be
                  // because this ploop should only run if we're doing chi3 flagged NL stuff
        realnum fwprev_2 = fw_2[i];
        realnum kapwkw_2 = kapw_2[kw_2], sigwkw_2 = sigw_2[kw_2];

        fw_2[i] = (fw_2_atZ[i] + fw_2_atZ[i + s] + fw_2_atZ[i - s1] + fw_2_atZ[i + (s - s1)]) *
                  0.25; // interpolation here.
        //(Gets 'Ex fields at X cell locations' from 'Ex fields at Z cell locations')
        f_2[i] += (kapwkw_2 + sigwkw_2) * fw_2[i] - (kapwkw_2 - sigwkw_2) * fwprev_2; // update f_2
      }

      KSTRIDE_DEF(dsigw_3, kw_3, is_3, gv);
      PLOOP_OVER_IVECS(gv, is_3, ie, i) { /// Round three for interpolating Y

                  if (chi2new[i] == 0.0) { continue; } // TODO should this be in these two interpolation loops??

        DEF_kw_3;
        realnum fwprev_3 = fw_3[i];
        realnum kapwkw_3 = kapw_3[kw_3], sigwkw_3 = sigw_3[kw_3];

        fw_3[i] =
            (fw_3_atZ[i] + fw_3_atZ[i + s] + fw_3_atZ[i - s2] + fw_3_atZ[i + (s - s2)]) *
            0.25; // interpolation here
                  //(Gets 'Ey fields at Y cell locations' from 'Ey fields at Z cell locations')
        f_3[i] += (kapwkw_3 + sigwkw_3) * fw_3[i] - (kapwkw_3 - sigwkw_3) * fwprev_3;
      }

      //                                   z            x            z     x
      // realnum g1sZatX = g1Z[i] + g1[i + s] + g1[i - s1] + g1[i + (s - s1)];
  ///  }
  }
 
  else {            /////////////// no PML (no fw) ///////////////////
  ///cout << "in NONPml NL case 1 " << endl;
    
    ///  if (chi3) { /// // TODO delete if
        
          /// NR Solver. Similar (simplified as no separate fw and f field updates) to the version above for the PML case(which probably won't be used but who knows))...
          /// Again, callable by defining the 2nd order NL material with a full 3x3 offdiagonal epsilon as a flag (i.e. u1 && u2 must be non-zero, however they will NOT be used in the calc) and also
          /// setting chi3 as nonzero as a flag (chi3 will not be used, but setting it to non-zero is a hack to make the code enter this 'if' statement for the NL material
         
          ///  fw_2_atZ and fw_3_atZ need to have the same dim as f. If f[i] is i'th Ex field, f_2/3 are the Ex and Ey fields at the SAME location as the Ez field. According
          /// to the Yee cell, the Ex and Ey fields are not in fact at the same location as Ez or each other. The Ex and Ey fields at the correct yee cell locations 
          /// must therefore subsequently calculated by interpolation (same principle as for gs_2 below), in subsequent ploopoverivecs
  const size_t g1size = gv.ntot();
  const size_t g2size = gv.ntot();
  const size_t fw2zsize = gv.ntot();
  const size_t fw3zsize = gv.ntot();

        PLOOP_OVER_IVECS(gv, is, ie, i) {
                                   
            if (chi2new[i] == 0.0) { continue; }// 

            if (i + s >= g1size || i - s2 < 0 || i + (s - s2) >= g1size) {
              std::cerr << "1: g1 out of bounds at i = " << i << "\n";
              sleep(15);
            }
            if (i + s >= g2size || i - s2 < 0 || i + (s - s2) >= g2size) {
              std::cerr << "1: g2 out of bounds at i = " << i << "\n";
              sleep(15);
            }
            if (i + s >= fw2zsize || i - s2 < 0 || i + (s - s2) >= fw2zsize) {
              std::cerr << "1: fw2z out of bounds at i = " << i << "\n";
              sleep(15);
            }
            if (i + s >= fw3zsize || i - s2 < 0 || i + (s - s2) >= fw3zsize) {
              std::cerr << "1: fw3z out of bounds at i = " << i << "\n";
              sleep(15);
            }


            realnum gs = g[i]; //dmpZ
          // avg orthogonal D-P fields over adjacent cells (see yee cell diag to understand why...):
            realnum gs_2 = (g1[i] + g1[i + s] + g1[i - s1] + g1[i + (s - s1)]) * 0.25; //dmpX at Z locations
            realnum gs_3 = (g2[i] + g2[i + s] + g2[i - s2] + g2[i + (s - s2)]) * 0.25; // dmpY at Z locations

                      if (i % 150 == 13) {
                          cout << "at entry  " << fw_2_atZ[i] << "  " << fw_3_atZ[i] << "  "<<  f_3[i]  << endl;
                          cout << "at entry dmp  " << g[i] << "  " << g1[i] << "  "<<  g2[i]  << endl;
                          cout << "dmp interp  " << gs_2 << "  " << gs_3  << endl;
                          
                    //  cout << u[i] << "    " << u_2[i] << "     " << u_2[i + s] << "   " << u_3[i]                           << endl;
                    }


            /// taking inverse of chi1inverse is easiest way to access epsilon...
            //realnum us = 1 / u[i]; 
            //realnum us_2 = 1 / u_2[i]; 
            //realnum us_3 = 1 / u_3[i];          
           

            if (u[i] == 0 || u_2[i] == 0 || u_3[i] == 0) {
              cout << "u is zero!! " << u[i] << "  " << u_2[i] << "  "<< u_3[i]<<endl;
              sleep(5);

            }
            realnum us = 1 / u[i]   ; 
            realnum us_2 = 1 / (   u_2[i]   );
            realnum us_3 = 1 / u_3[i];

          ///cout << "lgsfdg " << endl;
          //  cout << "us " << us << " us2 " << us_2 << endl;

            // will be format Parameters p1 = {prevF D-P_X, eps, 0, 0, 0, chi2new, 0, 0 } etc;
            Parameters p1 = {gs_2, us_2, 0.0, 0.0, 0.0, chi2new[i], 0.0, 0.0}; // X 
            Parameters p2 = {gs_3, us_3, 0.0, 0.0, 0.0, 0.0, chi2new[i], 0.0}; // Y
            Parameters p3 = {gs, us,     0.0, 0.0, 0.0, 0.0, 0.0, chi2new[i]}; // Z. currently using all chi2 tensor components equal (as zinc blende)

            realnum seed1 = f[i];
            realnum seed2 = fw_2_atZ[i];  
            realnum seed3 = fw_3_atZ[i];

          /*        cout << "PRENR s1" << endl;
            cout << seed1 << endl;
                  cout << " fw_2_atZ[i]" << seed2 << endl;
            cout << " s3" << seed3 << endl;
                  cout << " chi2:" << chi2new[i] << endl;
            cout << "chi3:" << chi3[i] << endl;
                  cout << " us " << us << endl;
            cout << " us_2 " << us_2 << endl;
                  cout << " us_3 " << us_3 << endl;
                  cout << " gs_2" << gs_2 << endl;*/

           /* fw_2_atZ[i] = gs_2 * u_2[i];
            fw_3_atZ[i] = gs_3 * u_3[i]; /// TODO THIS IS JUST A CHECK TEMPORARIOLY
            f[i] = gs * u[i];*/
            fw_2_atZ[i] = (  (g1[i] + g1[i - s1]) * u_2[i]   +    (g1[i + s] + g1[i + s - s1])*u_2[i+s]   )* 0.25;
            fw_3_atZ[i] = (  (g2[i] + g2[i - s2]) * u_3[i]   +    (g2[i + s] + g2[i + s - s2])*u_3[i+s]   )* 0.25; /// TODO THIS IS JUST A CHECK TEMPORARIOLY
            f[i] = gs * u[i];

            ///Newton Raphson for calculating Ez, Ex and Ey fields, (AT Z LOCATIONS):
            /// Seeded with previous field vals. Passing in field array pointers to be assigned new vals.
        //    runNR(seed2, seed3, seed1, &fw_2_atZ[i], &fw_3_atZ[i], &f[i], p1, p2, p3); // note fw_2_atZ variable is named 'fw' but is used for 'f' here
          //  if (i % 300 == 0) {
            if (isnan(fw_2_atZ[i]) ||  isnan(fw_3_atZ[i]) || isnan(f[i]) ||isinf(fw_2_atZ[i]) ||  isinf(fw_3_atZ[i]) || isinf(f[i]) ){
        cout << "err " << fw_2_atZ[i] << "    " << fw_3_atZ[i] << "     " << f[i] << endl;
              cout << u[i] << "    " << u_2[i] << "     " << u_2[i + s] << "   " << u_3[i]<< endl;
        sleep(5);
            }
           /*   if (i % 500 == 0) {
        cout<< "general print"<< fw_2_atZ[i] << "    " << fw_3_atZ[i] << "     " << f[i] << endl;
        cout << u[i] << "    " << u_2[i] << "     " << u_2[i + s] << "   " << u_3[i] << endl;

            }*/

        }

        // now do the other two PLOOPs to interpolate the X and Y fields to their correct positions, and then calculate f_2 and f_3..
        PLOOP_OVER_IVECS(gv, is_2, ie, i) { /// Round two for interpolating X

                                if (chi2new[i] == 0.0) { continue; }// TODO should this be in these two interpolation loops??

             //(Gets 'Ex fields at X cell locations' from 'Ex fields at Z cell locations')
                     //   f_2[i] = (fw_2_atZ[i] + fw_2_atZ[i + s] + fw_2_atZ[i - s1] + fw_2_atZ[i + (s - s1)]) * 0.25; // interpolation here. //TODO THIS IS ERROR?
        }

        PLOOP_OVER_IVECS(gv, is_3, ie, i) { /// Round three for interpolating Y

                                if (chi2new[i] == 0.0) { continue; }// TODO should this be in these two interpolation loops??

          //(Gets 'Ey fields at y cell locations' from 'Ey fields at Z cell locations')
                     f_3[i] = (   fw_3_atZ[i] + fw_3_atZ[i + s] + fw_3_atZ[i - s2] + fw_3_atZ[i + (s - s2)]    )*0.25; // interpolation here  //TODO THIS IS ERROR?
                   }
        //                                   z            x            z     x
        // realnum g1sZatX = g1Z[i] + g1[i + s] + g1[i - s1] + g1[i + (s - s1)];

      

     /// }

  }
}

//#ifndef DONT_CALL_NL
//void __dummy_call_to_step_update_EDHB_NL() {
//  step_update_EDHB_NL(nullptr, nullptr, nullptr, NO_COMPONENT, *(grid_volume *)nullptr, ivec(),
//                      ivec(), ivec(), ivec(), nullptr, nullptr, nullptr, nullptr, nullptr, nullptr,
//                      nullptr, nullptr, 0, 0, 0, nullptr, nullptr, nullptr, nullptr, nullptr,
//                      nullptr, NO_DIRECTION, NO_DIRECTION, NO_DIRECTION, nullptr, nullptr, nullptr,
//                      nullptr, nullptr, nullptr);
//}
//#endif

} // namespace meep




///// TODO rename macro and fn defs for this
//void step_update_EDHB_NL(RPR f, RPR f_2, RPR f_3, component fc, const grid_volume &gv,
//                         const ivec is, const ivec is_2, const ivec is_3, const ivec ie,
//                         const RPR g, const RPR g1, const RPR g2, const RPR u, const RPR u_2,
//                         const RPR u_3, const RPR u1, const RPR u2, ptrdiff_t s, ptrdiff_t s1,
//                         ptrdiff_t s2, const realnum *chi2new, const RPR chi3, RPR fw, RPR fw_2_atZ,
//                         RPR fw_3_atZ, RPR fw_2, RPR fw_3, direction dsigw, direction dsigw_2,
//                         direction dsigw_3, const RPR sigw, const RPR sigw_2, const RPR sigw_3,
//                         const RPR kapw, const RPR kapw_2, const RPR kapw_3) {
//  (void)fc; // currently unused
//  if (!f) return;
//
//// stable averaging of offdiagonal components
//#define OFFDIAG(u, g, sx)                                                                          \
//  (0.25 * ((g[i] + g[i - sx]) * u[i] + (g[i + s] + g[(i + s) - sx]) * u[i + s]))
//
//  if (dsigw != NO_DIRECTION) { //////// PML case (with fw) /////////////  // TODO need to also
//                               /// implement for non pml case, since this fn operates chunk-wise,
//                               /// and some chunks won't have pml even if some do.
//    KSTRIDE_DEF(dsigw, kw, is, gv); /// Used in DEF_kw.  dsigw, is, and gv come from fn. kw0,
//                                    /// skw1/2/3 out. TODO check wrt position
//
//    // TODO implement off-diagonal epsilon terms into NewtonRaphson solver
//    ///  if (chi3) {
//    /// Build NR solver into this section here. **should add similar copy to the section further
//    /// down without PML... Will be callable by defining the 2nd order NL material with a full 3x3
//    /// offdiagonal epsilon as a flag (i.e. u1 && u2 must be non-zero, however they will NOT be
//    /// used in the calc) and also setting chi3 as nonzero as a flag (chi3 will not be used, but
//    /// setting it to non-zero is a hack to make the code enter this 'if' statement for the NL
//    /// material
//    ///
//    /// f, f_2 and f_3 (i.e, Z, X, and Y respectively) are (in non-PML case) i-th cell electric
//    /// field the output from this fn. In PML case they are updated as per original meep (but not
//    /// clear what exactly the fields represent in this case)... fc = field cmpnt gv = grid volume
//    /// (need to ploopoverivecs) is, ie, start and end of chunk. is, is_2, is_3 are start ivecs
//    /// for Z, X, and Y directions respectively g, g1, g2 = 'D - P'  field components z, x, y u is
//    /// inverse epsilon Z direction. u1, u2 = offdiagonal inverse epsilon(NOT used - except as
//    /// flag). u_2 and u_3 are inverse epsilon in the X and Y directions s, s1, s2 = strides (for
//    /// getting avged orthogonal fields over adjacent cells (see yeecell diag it makes sense). For
//    /// this NL fn, these will have been set in the most recent 'FOR_FT_COMPONENTS' loop, in which
//    /// 'ec' would have been Ez, therefore s is Ez stride, s1 and s2 are x and y respectively...
//    /// TODO - DOUBLE CHECK THIS IS TRUE from update_eh! chi2 = actual simulation value of chi2
//    /// (in meep units..?), will be passed into NR chi3 = a flag (set to nonzero float to indicate
//    /// NL material) fw[i] = PML case: i-th cell Ez field the output from this fn. fw_2_atZ and
//    /// fw_3_atZ are X and Y E fields AT Z LOCATIONS (see below). fw_2/3 are for the final
//    /// interpolated Ex and Ey fields. dsigw =           /// sigw =           /// kapw =
//
//    ///  fw_2_atZ and fw_3_atZ need to have the same dim as fw. If fw[i] is i'th Ex field, fw_2/3
//    ///  are the Ex and Ey fields at the SAME location as the Ez field. According
//    /// to the Yee cell, the Ex and Ey fields are not in fact at the same location as Ez or each
//    /// other. The Ex and Ey fields at the correct yee cell locations must therefore subsequently
//    /// calculated by interpolation (same principle as for gs_2 below), in subsequent
//    /// ploopoverivecs
//
//    PLOOP_OVER_IVECS(gv, is, ie, i) {
//
//      if (chi2new[i] == 0.0) { continue; }
//      realnum gs = g[i]; // dmpZ
//      // avg orthogonal D-P fields over adjacent cells (see yee cell diag to understand why...)
//      realnum gs_2 = (g1[i] + g1[i + s] + g1[i - s1] + g1[i + (s - s1)]) * 0.25; // dmpX
//      realnum gs_3 = (g2[i] + g2[i + s] + g2[i - s2] + g2[i + (s - s2)]) * 0.25; // dmpY
//
//      /// taking inverse of chi1inverse is easiest way to access epsilon...
//      realnum us = 1 / u[i];
//      realnum us_2 = 1 / u_2[i];
//      realnum us_3 = 1 / u_3[i];
//
//      /// #will be format Parameters p1 = {prevF D-P_X, eps, 0, 0, 0, chi2, 0, 0 } etc;
//      Parameters p1 = {gs_2, us_2, 0.0, 0.0, 0.0, chi2new[i], 0.0, 0.0}; // X
//      Parameters p2 = {gs_3, us_3, 0.0, 0.0, 0.0, 0.0, chi2new[i], 0.0}; // Y
//      Parameters p3 = {gs,  us,  0.0, 0.0,
//                       0.0, 0.0, 0.0, chi2new[i]}; // Z. currently using all chi2 tensor components
//                                                   // equal (as zinc blende)
//
//      realnum seed1 = fw[i];
//      realnum seed2 =
//          fw_2_atZ[i]; // TODO THIS MIGHT FAIL BECAUSE FW FIELDS MAY NOT YET HAVE BEEN INITIALISED
//                       // SO MAY NOT BE ABLE TO BE USED AS A SEED NUMBER ON FIRST LOOP...
//      realnum seed3 = fw_3_atZ[i];
//
//      /// Newton Raphson for calculating Ez, Ex and Ey fields, (AT Z LOCATIONS):
//      ///  Seeded with previous field vals. Passing in field array pointers to be assigned new
//      ///  vals.
//      runNR(seed2, seed3, seed1, &fw_2_atZ[i], &fw_3_atZ[i], &fw[i], p1, p2, p3);
//
//      // Do the other fields for PML (whatever they do exactly..)
//      DEF_kw; /// TODO  - might not be relevant because its for PML and PML won't use this
//              /// implementation so long as we don't have a chi3(flag) defined within a PML
//              /// layer?!
//      /// uses kw0, skw1/2/3. not sure what it does exactly...
//      realnum fwprev = fw[i], kapwkw = kapw[kw], sigwkw = sigw[kw]; // Ez
//      f[i] += (kapwkw + sigwkw) * fw[i] - (kapwkw - sigwkw) * fwprev;
//    }
//
//    // now do the other two PLOOPs to interpolate the X and Y fields to their correct positions,
//    // and then calculate f_2 and f_3..
//    KSTRIDE_DEF(dsigw_2, kw_2, is_2, gv);
//    PLOOP_OVER_IVECS(gv, is_2, ie, i) { /// Round two for interpolating X
//
//      if (chi2new[i] == 0.0) { continue; } // TODO should this be in these two interpolation loops??
//
//      DEF_kw_2; // these will throw if the necessary values are null, but they shouldn't be
//                // because this ploop should only run if we're doing chi3 flagged NL stuff
//      realnum fwprev_2 = fw_2[i];
//      realnum kapwkw_2 = kapw_2[kw_2], sigwkw_2 = sigw_2[kw_2];
//
//      fw_2[i] = (fw_2_atZ[i] + fw_2_atZ[i + s] + fw_2_atZ[i - s1] + fw_2_atZ[i + (s - s1)]) *
//                0.25; // interpolation here.
//      //(Gets 'Ex fields at X cell locations' from 'Ex fields at Z cell locations')
//      f_2[i] += (kapwkw_2 + sigwkw_2) * fw_2[i] - (kapwkw_2 - sigwkw_2) * fwprev_2; // update f_2
//    }
//
//    KSTRIDE_DEF(dsigw_3, kw_3, is_3, gv);
//    PLOOP_OVER_IVECS(gv, is_3, ie, i) { /// Round three for interpolating Y
//
//      if (chi2new[i] == 0.0) { continue; } // TODO should this be in these two interpolation loops??
//
//      DEF_kw_3;
//      realnum fwprev_3 = fw_3[i];
//      realnum kapwkw_3 = kapw_3[kw_3], sigwkw_3 = sigw_3[kw_3];
//
//      fw_3[i] = (fw_3_atZ[i] + fw_3_atZ[i + s] + fw_3_atZ[i - s2] + fw_3_atZ[i + (s - s2)]) *
//                0.25; // interpolation here
//                      //(Gets 'Ey fields at Y cell locations' from 'Ey fields at Z cell locations')
//      f_3[i] += (kapwkw_3 + sigwkw_3) * fw_3[i] - (kapwkw_3 - sigwkw_3) * fwprev_3;
//    }
//
//    //                                   z            x            z     x
//    // realnum g1sZatX = g1Z[i] + g1[i + s] + g1[i - s1] + g1[i + (s - s1)];
//    ///  }
//  }
//
//  else { /////////////// no PML (no fw) ///////////////////
//
//    ///  if (chi3) { /// // TODO delete if
//
//    /// NR Solver. Similar (simplified as no separate fw and f field updates) to the version above
//    /// for the PML case(which probably won't be used but who knows))... Again, callable by defining
//    /// the 2nd order NL material with a full 3x3 offdiagonal epsilon as a flag (i.e. u1 && u2 must
//    /// be non-zero, however they will NOT be used in the calc) and also setting chi3 as nonzero as
//    /// a flag (chi3 will not be used, but setting it to non-zero is a hack to make the code enter
//    /// this 'if' statement for the NL material
//
//    ///  fw_2_atZ and fw_3_atZ need to have the same dim as f. If f[i] is i'th Ex field, f_2/3 are
//    ///  the Ex and Ey fields at the SAME location as the Ez field. According
//    /// to the Yee cell, the Ex and Ey fields are not in fact at the same location as Ez or each
//    /// other. The Ex and Ey fields at the correct yee cell locations must therefore subsequently
//    /// calculated by interpolation (same principle as for gs_2 below), in subsequent ploopoverivecs
//
//    PLOOP_OVER_IVECS(gv, is, ie, i) {
//
//      if (chi2new[i] == 0.0) { continue; } //
//
//      realnum gs = g[i]; // dmpZ
//      // avg orthogonal D-P fields over adjacent cells (see yee cell diag to understand why...):
//      realnum gs_2 = (g1[i] + g1[i + s] + g1[i - s1] + g1[i + (s - s1)]) * 0.25; // dmpX
//      realnum gs_3 = (g2[i] + g2[i + s] + g2[i - s2] + g2[i + (s - s2)]) * 0.25; // dmpY
//
//      /// taking inverse of chi1inverse is easiest way to access epsilon...
//      realnum us = 1 / u[i];
//      realnum us_2 = 1 / u_2[i];
//      realnum us_3 = 1 / u_3[i];
//
//      // will be format Parameters p1 = {prevF D-P_X, eps, 0, 0, 0, chi2new, 0, 0 } etc;
//      Parameters p1 = {gs_2, us_2, 0.0, 0.0, 0.0, chi2new[i], 0.0, 0.0}; // X
//      Parameters p2 = {gs_3, us_3, 0.0, 0.0, 0.0, 0.0, chi2new[i], 0.0}; // Y
//      Parameters p3 = {gs,  us,  0.0,       0.0, 0.0,
//                       0.0, 0.0, chi2new[i]}; // Z. currently using all chi2 tensor components equal
//                                              // (as zinc blende)
//
//      realnum seed1 = f[i];
//      realnum seed2 =
//          fw_2_atZ[i]; // TODO THIS MIGHT FAIL BECAUSE FW FIELDS MAY NOT YET HAVE BEEN INITIALISED
//                       // SO MAY NOT BE ABLE TO BE USED AS A SEED NUMBER ON FIRST LOOP...
//      realnum seed3 = fw_3_atZ[i];
//
//      cout << "PRENR s1" << seed1 << " fw_2_atZ[i]" << fw_2_atZ[i] << " s3" << seed3
//           << " chi2:" << chi2new[i] << "chi3:" << chi3[i] << " us " << us << " us_2 " << us_2
//           << " us_3 " << us_3 << " u1 " << u1[i] << " gs_2" << gs_2 << endl;
//
//      /// Newton Raphson for calculating Ez, Ex and Ey fields, (AT Z LOCATIONS):
//      ///  Seeded with previous field vals. Passing in field array pointers to be assigned new vals.
//      runNR(seed2, seed3, seed1, &fw_3_atZ[i], &fw_2_atZ[i], &f[i], p1, p2,
//            p3); // note fw_2_atZ variable is named 'fw' but is used for 'f' here
//    }
//
//    // now do the other two PLOOPs to interpolate the X and Y fields to their correct positions, and
//    // then calculate f_2 and f_3..
//    PLOOP_OVER_IVECS(gv, is_2, ie, i) { /// Round two for interpolating X
//
//      if (chi2new[i] == 0.0) { continue; } // TODO should this be in these two interpolation loops??
//
//      //(Gets 'Ex fields at X cell locations' from 'Ex fields at Z cell locations')
//      f_2[i] = (fw_2_atZ[i] + fw_2_atZ[i + s] + fw_2_atZ[i - s1] + fw_2_atZ[i + (s - s1)]) *
//               0.25; // interpolation here.
//    }
//
//    PLOOP_OVER_IVECS(gv, is_3, ie, i) { /// Round three for interpolating Y
//
//      if (chi2new[i] == 0.0) { continue; } // TODO should this be in these two interpolation loops??
//
//      //(Gets 'Ey fields at y cell locations' from 'Ey fields at Z cell locations')
//      f_3[i] = (fw_3_atZ[i] + fw_3_atZ[i + s] + fw_3_atZ[i - s2] + fw_3_atZ[i + (s - s2)]) *
//               0.25; // interpolation here
//    }
//    //                                   z            x            z     x
//    // realnum g1sZatX = g1Z[i] + g1[i + s] + g1[i - s1] + g1[i + (s - s1)];
//
//    /// }
//  }
//}