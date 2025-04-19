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

#include <string.h>
#include <assert.h>
#include "meep.hpp"
#include "meep_internals.hpp"
#include <iostream>


using namespace std;

namespace meep {

void fields::update_eh(field_type ft, bool skip_w_components) {
  if (ft != E_stuff && ft != H_stuff) meep::abort("update_eh only works with E/H");

  // split the chunks' volume into subdomains for tiled execution of update_eh loop
  for (int i = 0; i < num_chunks; i++)
    if (chunks[i]->is_mine() && changed_materials) {
      bool is_aniso = false;
      FOR_FT_COMPONENTS(ft, cc) {
        const direction d_c = component_direction(cc);
        const direction d_1 = cycle_direction(chunks[i]->gv.dim, d_c, 1);
        const direction d_2 = cycle_direction(chunks[i]->gv.dim, d_c, 2);
        if (chunks[i]->s->chi1inv[cc][d_1] && chunks[i]->s->chi1inv[cc][d_2]) {
          is_aniso = true;
          break;
        }
      }
      if (!chunks[i]->gvs_eh[ft].empty()) chunks[i]->gvs_eh[ft].clear();
      if (loop_tile_base_eh > 0 && is_aniso) {
        split_into_tiles(chunks[i]->gv, &chunks[i]->gvs_eh[ft], loop_tile_base_eh);
        check_tiles(chunks[i]->gv, chunks[i]->gvs_eh[ft]);
      }
      else { chunks[i]->gvs_eh[ft].push_back(chunks[i]->gv); }
    }

  for (int i = 0; i < num_chunks; i++)
    if (chunks[i]->is_mine())
      if (chunks[i]->update_eh(ft, skip_w_components)) {
        chunk_connections_valid = false; // E/H allocated - reconnect chunks
        assert(changed_materials);
      }
}

bool fields_chunk::needs_W_prev(component c) const {
  for (susceptibility *chiP = s->chiP[type(c)]; chiP; chiP = chiP->next)
    if (chiP->needs_W_prev()) return true;
  return false;
}

bool fields_chunk::update_eh(field_type ft, bool skip_w_components) {
  field_type ft2 = ft == E_stuff ? D_stuff : B_stuff; // for sources etc.
  bool allocated_eh = false;

  bool have_int_sources = false;
  if (!doing_solve_cw) {
    for (const src_vol &sv : sources[ft2]) {
      if (sv.t()->is_integrated) {
        have_int_sources = true;
        break;
      }
    }
  }


  //########## ALL CODE FROM HERE TO THE BOTTOM OF THE PAGE HAS BEEN COPY-PASTED FROM MY OLD ANNOTATED COPY OF MEEP - THE CODE SHOULD BE THE SAME, JUST HAS USEFUL COMMENTS...

FOR_FT_COMPONENTS(ft,ec) { // Iter thro field type components, i.e., for Estuff; Ex, Ey, Er, Ep Ez...
    component dc =
        field_type_component(ft2, ec); // ft2 is D_stuff. dc = Dx, Dy, etc for ec = Ex, Ey,...
    DOCMP {                            // Re and Im
      bool need_fmp = false;
      if (f[ec][cmp]) { // if E field component is nonzero...
        need_fmp = have_int_sources;
        for (polarization_state *p = pol[ft]; p && !need_fmp; p = p->next)
          need_fmp = need_fmp || p->s->needs_P(ec, cmp, f); // allocate memory?
      }
      if (need_fmp) { // allocate memory
        if (!f_minus_p[dc][cmp])
          f_minus_p[dc][cmp] =
              new realnum[gv.ntot()]; // set each field and component an empty array of size
                                      // sufficient for the gv (i.e. points in the chunk)...
      }
      else if (f_minus_p[dc][cmp]) { // remove unneeded f_minus_p
        delete[] f_minus_p[dc][cmp];
        f_minus_p[dc][cmp] = 0; // bind to 0 so as to avoid having a dangling pointer
      }
    }
  }
  bool have_f_minus_p = false;
  FOR_FT_COMPONENTS(ft2, dc) { // check if any d field components are actually required
    if (f_minus_p[dc][0]) {
      have_f_minus_p = true;
      break;
    }
  }

  const size_t ntot = s->gv.ntot();

  if (have_f_minus_p && doing_solve_cw)
    meep::abort("dispersive materials are not yet implemented for solve_cw");

  //////////////////////////////////////////////////////////////////////////
  // First, initialize f_minus_p to D - P, if necessary

  FOR_FT_COMPONENTS(ft, ec) if (f[ec][0]) {
    component dc = field_type_component(ft2, ec);
    DOCMP if (f_minus_p[dc][cmp]) {
      realnum *fmp = f_minus_p[dc][cmp];
      memcpy(fmp, f[dc][cmp], sizeof(realnum) * ntot);
    }
  }

  for (polarization_state *p = pol[ft]; p; p = p->next)
    if (p->data) p->s->subtract_P(ft, f_minus_p, p->data);

  //////////////////////////////////////////////////////////////////////////
  // Next, subtract time-integrated sources (i.e. polarizations, not currents)

  if (have_f_minus_p && !doing_solve_cw) {
    for (const src_vol &sv : sources[ft2]) {
      if (sv.t()->is_integrated && f[sv.c][0] && ft == type(sv.c)) {
        component c = field_type_component(ft2, sv.c);
        for (size_t j = 0; j < sv.num_points(); ++j) {
          const complex<double> A = sv.dipole(j);
          DOCMP { f_minus_p[c][cmp][sv.index_at(j)] -= (cmp) ? imag(A) : real(A); }
        }
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////
  // Finally, compute E = chi1inv * D

  realnum *dmp[NUM_FIELD_COMPONENTS][2];
  FOR_FT_COMPONENTS(ft2, dc) DOCMP2 {
    dmp[dc][cmp] = f_minus_p[dc][cmp] ? f_minus_p[dc][cmp] : f[dc][cmp];
  }


  /// New local field component defs just for convenience:
     const component ex = Ex;
  const component ey = Ey;
  const component ez = Ez; // <<<TODO NEED to verify this, it's unclear how meep is actually
                           // allocating the fields in the f[component] arrays, because Ez is after
                           // the cylindrical coordinate cpmnts in the switch
  //, which confuses things... (Is it accidentally just storing Ez in the Er array when Dcyl is
  ///false...? I don't think so but it's possible...

  /// ############### THIS is the loop that needs modifying for Newton Raphson implementation
  for (size_t i = 0; i < gvs_eh[ft].size(); ++i) {
    DOCMP { /// added { here to split the two loop macros
      /// ADDED variable declarations here (and placeholder inits), so they are in-scope both in AND
      /// after FOR_FT_COMPONENTS
      component ecInLoop = ex;
      component dc = field_type_component(D_stuff, ex);
      direction d_ec = component_direction(ex); // Direction of main field component  (e.g. X)
      ptrdiff_t s_ec = gv.stride(d_ec);         // stride of main field component
      direction d_1 =
          cycle_direction(gv.dim, d_ec, 1);          // Direction of next field component (e.g. Y)
      component dc_1 = direction_component(dc, d_1); // next field component
      ptrdiff_t s_1 = gv.stride(d_1);                // stride of next field component
      direction d_2 = cycle_direction(gv.dim, d_ec, 2); // '' #2   (e.g. Z)
      component dc_2 = direction_component(dc, d_2);    // '' #2
      ptrdiff_t s_2 = gv.stride(d_2);                   // '' #2
      direction dsigw0 = d_ec;
      direction dsigw = NO_DIRECTION;

      FOR_FT_COMPONENTS(ft, ec) { // For E_Stuff, loop over ec = Ex, Ey, etc
        if (f[ec][cmp]) {
          if (type(ec) != ft) meep::abort("bug in FOR_FT_COMPONENTS");
          /// TODO check if these are the right thing when they go into the NL version
          dc = field_type_component(ft2, ec);
          d_ec = component_direction(ec); // Direction of main field component  (e.g. X)
          s_ec = gv.stride(d_ec) * (ft == H_stuff ? -1 : +1); // stride of main field component
          d_1 = cycle_direction(gv.dim, d_ec, 1); // Direction of next field component (e.g. Y)
          dc_1 = direction_component(dc, d_1);    // next field component
          s_1 = gv.stride(d_1) * (ft == H_stuff ? -1 : +1); // stride of next field component
          d_2 = cycle_direction(gv.dim, d_ec, 2);           // '' #2   (e.g. Z)
          dc_2 = direction_component(dc, d_2);              // '' #2
          s_2 = gv.stride(d_2) * (ft == H_stuff ? -1 : +1); // '' #2
          dsigw0 = d_ec;
          dsigw = s->sigsize[dsigw0] > 1 ? dsigw0 : NO_DIRECTION;
          ecInLoop = ec;

          // lazily allocate any E/H fields that are needed (H==B initially)
          if (i == 0 && f[ec][cmp] == f[dc][cmp] &&
              (s->chi1inv[ec][d_ec] || have_f_minus_p || dsigw != NO_DIRECTION)) {
            f[ec][cmp] = new realnum[gv.ntot()];
            memcpy(f[ec][cmp], f[dc][cmp], gv.ntot() * sizeof(realnum));
            allocated_eh = true;
          }

          // lazily allocate W auxiliary field
          if (i == 0 && !f_w[ec][cmp] && dsigw != NO_DIRECTION) {
            f_w[ec][cmp] = new realnum[gv.ntot()];
            memcpy(f_w[ec][cmp], f[ec][cmp], gv.ntot() * sizeof(realnum));
            if (needs_W_notowned(ec)) allocated_eh = true; // communication needed
          }

          // for solve_cw, when W exists we get W and E from special variables
          if (f_w[ec][cmp] && skip_w_components) continue;

          // save W field from this timestep in f_w_prev if needed by pols
          if (i == 0 && needs_W_prev(ec)) {
            if (!f_w_prev[ec][cmp]) f_w_prev[ec][cmp] = new realnum[gv.ntot()];
            memcpy(f_w_prev[ec][cmp], f_w[ec][cmp] ? f_w[ec][cmp] : f[ec][cmp],
                   sizeof(realnum) * gv.ntot());
          }

          /// if (!s->chi3[ec] || ft == H_stuff ) { /// Add this 'if not chi3' (hack, using chi3 as
          /// a flag, actual value not relevant so long as it is non-zero for the 2nd order NL
          /// material) statement wrapper around this STEP_UPDATE_EDHB, so that the non-nonlinear
          /// chunks
          /// run field-component-piecewise. TODO need to check s->chi3[ec] does what I want tho.
          /// TODO may be worth making a CHECKPOINT here print ft to check H_stuff isn't going into
          /// NR loop...

          if (f[ec][cmp] != f[dc][cmp]) {
            STEP_UPDATE_EDHB(
                f[ec][cmp], ec, gv, gvs_eh[ft][i].little_owned_corner0(ec),
                gvs_eh[ft][i].big_corner(), dmp[dc][cmp], dmp[dc_1][cmp], dmp[dc_2][cmp],
                s->chi1inv[ec][d_ec], dmp[dc_1][cmp] ? s->chi1inv[ec][d_1] : NULL,
                dmp[dc_2][cmp] ? s->chi1inv[ec][d_2] : NULL, s_ec, s_1, s_2,
                NULL, /// CHI2 set to NULL here as will be applied in subsequent STEP_UPDATE_EDHB_NL
                s->chi3[ec], f_w[ec][cmp], dsigw, s->sig[dsigw], s->kap[dsigw]);

            // if (gv.dim == Dcyl) {
            //  ivec is = gvs_eh[ft][i].little_owned_corner(ec);
            //  if (is.r() == 0) {
            //    ivec ie = gvs_eh[ft][i].big_corner();
            //    ie.set_direction(R, 0);
            //    /* pass NULL for off-diagonal terms since they must be
            //       zero at r=0 for an axisymmetric structure: */
            //    STEP_UPDATE_EDHB(f[ec][cmp], ec, gv, is, ie, dmp[dc][cmp], NULL, NULL,
            //                     s->chi1inv[ec][d_ec], NULL, NULL, s_ec, s_1, s_2, s->chi2[ec],
            //                     s->chi3[ec], f_w[ec][cmp], dsigw, s->sig[dsigw], s->kap[dsigw]);
            //  }
            //}
          }
        }
        cout << "Done linear" << ft << endl;
      } /// The FOR_FT_COMPONENTS loop should close out on ec = Ez, therefore for convenience, start
        /// the NL STEP_UPDATE_EDHB with the main field and field locations as Z, and X and Y as
      /// the 'auxiliaries' which are to be calculated by interpolation...

        cout << "Done linear 2 " <<  endl;  
      cout << "Done linear 2b " << typeid(s->chi2[ez]).name() << endl; // Pd
        cout << "Done linear 2b " << typeid(s->chi2[ez][0]).name() << endl; // 
      cout << "Done linear 2c "  << ecInLoop << endl;   // 4


      /// START OF NL VERSION OF STEP_UPDATE_EDHB>>>>>>>>> This is outside the FOR_FT_COMPONENTS
      /// loop becasue we're doing all 3 field components at once in the NR solver!
      /// TODO this bit is only for the PML case [is it??], need another one where it's not pml
      /// case..? handle the NL chunks and pass in all xyz field components in at once.
      if (s->chi2[ez] && ft == E_stuff) { /// if chi2 is non-zero (only z direction checked, but assumes defined for ALL axes)
          cout << "Doing Nonlinear 1: " << ft << "  " << s->chi2[ez]
               << "  " << ecInLoop << endl; /// of chi3... TODO check s->chi3[ec] 

        if (f[ez][cmp]) { // added this as it's also wrapping the stuffin FOR_FT_COMPONENTS
          cout << "Doing Nonlinear 2" << endl;
            //if (ec != ez) { // CHECKPOINT - TODO doesn't currently work because ec isn't accessible here anyway..
            //  std::cout << "ec != ez!! ec:" << ec
            //            << std::endl; /// TODO, ec might not print as a variable I suppose...
            //                          /// depends on it's type.
            //  meep::abort("ec != ez!! ec:");
            //}

            /// Now need to create the two extra temp Z dimensioned field arrays for storing Ex and
            /// Ey at the Z positions (for subsequent interpolation to their correct positions)
            // lazily allocate the two temp extra Z-dimensioned W auxiliary fields:
            if (i == 0 && !fTempNlFieldsForInterpolation[0][cmp]) {
              fTempNlFieldsForInterpolation[0][cmp] =
                  new realnum[gv.ntot()]; // for temp Ex at Z positions fields
              fTempNlFieldsForInterpolation[1][cmp] =
                  new realnum[gv.ntot()]; // for temp Ey at Z positions fields
            }

            direction dsigw0_2 = d_1; /// X  Additional dsigw terms for other two field directions
            direction dsigw_2 = s->sigsize[dsigw0] > 1 ? dsigw0_2 : NO_DIRECTION;
            direction dsigw0_3 = d_2; /// Y
            direction dsigw_3 = s->sigsize[dsigw0] > 1 ? dsigw0_3 : NO_DIRECTION;

            /// Now do nonlinear xyz e field step update:  /// TODO! need to ensure correct 'ec'
            /// components go into all these (i.e, z, x, y)
       //  if (f[ecInLoop][cmp] != f[dc][cmp]) { // not sure if this 'if' is still needed - might cause
                                            // probs? TODO - leave for now, see what happens
          STEP_UPDATE_EDHB_NL(f[ez][cmp], f[ex][cmp], f[ey][cmp],
              ez, gv, 
              gvs_eh[ft][i].little_owned_corner0(ez), 
              gvs_eh[ft][i].little_owned_corner0(ex),
              gvs_eh[ft][i].little_owned_corner0(ey), 
              gvs_eh[ft][i].big_corner(), 
              dmp[dc][cmp], dmp[dc_1][cmp], dmp[dc_2][cmp], //TODO check these dmp field compoennts
              s->chi1inv[ez][component_direction(ez)], s->chi1inv[ex][component_direction(ex)], s->chi1inv[ey][component_direction(ey)], // principal epsilon (inverse) components
              dmp[dc_1][cmp] ? s->chi1inv[ez][d_1] : NULL, dmp[dc_2][cmp] ? s->chi1inv[ez][d_2] : NULL, //TODO check - think these should be fine as they are the offdiag ones
              s_ec, s_1, s_2, /// strides of Z, X, and Y directions respectively
              s->chi2[ez], s->chi3[ez], 
              f_w[ez][cmp], fTempNlFieldsForInterpolation[0][cmp],
                    fTempNlFieldsForInterpolation[1][cmp], f_w[ex][cmp], f_w[ey][cmp],                    
              dsigw, dsigw_2, dsigw_3, 
              s->sig[dsigw], s->sig[dsigw_2], s->sig[dsigw_3],
              s->kap[dsigw], s->kap[dsigw_2], s->kap[dsigw_3]);

          //}
        }
      } /// end of new NL stuff
      cout << "next DOCMP" << endl;
    } /// added } here to split the two loop macros
  }

  return allocated_eh;
}

} // namespace meep