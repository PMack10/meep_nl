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
#include <vector>
#include <string.h>
#include "meep.hpp"

using namespace std;

namespace meep {

#define DOCMP for (int cmp = 0; cmp < 2 - is_real; cmp++)
#define DOCMP2 for (int cmp = 0; cmp < 2; cmp++)

// note that C99 has a round() function, but I don't want to rely on it
static inline int my_round(double x) { return int(floor(fabs(x) + 0.5) * (x < 0 ? -1 : 1)); }

/* implement mirror boundary conditions for i outside 0..n-1: */
int mirrorindex(int i, int n);

/* map the cell coordinates into the range [0,1] */
void map_coordinates(double rx, double ry, double rz, int nx, int ny, int nz, int &x1, int &y1,
                     int &z1, int &x2, int &y2, int &z2, double &dx, double &dy, double &dz,
                     bool do_fabs = true);

inline int small_r_metal(int m) { return m - 1; }

inline int rmin_bulk(int m) {
  int r = 1 + small_r_metal(m);
  if (r < 1) r = 1;
  return r;
}

// A source volume
// Moveable and copyable
class src_vol {
public:
  // Constructs a new source volume. Takes ownership of `ind` and `amps`.
  // Requirement: ind.size() == amps.size()
  src_vol(component cc, src_time *st, std::vector<ptrdiff_t> &&ind,
          std::vector<std::complex<double> > &&amps, bool fix_boundaries = false);

  // Checks whether `a` and `b` are combinable, i.e. have the same indices and point to the same
  // `src_time` instance, but have potentially different amplitudes.
  static bool combinable(const src_vol &a, const src_vol &b);

  ptrdiff_t index_at(size_t pos) const { return index[pos]; }
  std::complex<double> amplitude(size_t j) const { return amp[j]; };
  void set_amplitude(size_t j, std::complex<double> a) { amp[j] = a; };
  void set_amplitude(size_t j, double a) { amp[j] = a; };
  const std::complex<double> &amplitude_at(size_t pos) const { return amp[pos]; }
  size_t num_points() const { return index.size(); };
  const src_time *t() const { return src_t; };

  std::complex<double> dipole(size_t j) const { return amp[j] * src_t->dipole(); };
  std::complex<double> current(size_t j) const { return amp[j] * src_t->current(); };
  void update(double time, double dt) { src_t->update(time, dt); };
  // Merges the amplitudes from volume source `other` into `this`.
  // Requirement: other.num_points() == this->num_points().
  // It is recommended to use `combinable` before calling this method.
  void add_amplitudes_from(const src_vol &other);

  const component c;       // field component the source applies to
  bool needs_boundary_fix; // whether fix_boundary_sources needs calling
private:
  src_time *src_t;                        // Not owned by us.
  std::vector<ptrdiff_t> index;           // locations of sources in grid (indices)
  std::vector<std::complex<double> > amp; // amplitudes
};

const int num_bandpts = 32;

symmetry r_to_minus_r_symmetry(int m);

// functions in step_generic.cpp:

void step_curl(realnum *f, component c, const realnum *g1, const realnum *g2, ptrdiff_t s1,
               ptrdiff_t s2, // strides for g1/g2 shift
               const grid_volume &gv, const ivec is, const ivec ie, realnum dtdx, direction dsig,
               const realnum *sig, const realnum *kap, const realnum *siginv, realnum *fu,
               direction dsigu, const realnum *sigu, const realnum *kapu, const realnum *siginvu,
               realnum dt, const realnum *cnd, const realnum *cndinv, realnum *fcnd);


void step_update_EDHB(realnum *f, component fc, const grid_volume &gv, const ivec is, const ivec ie,
                      const realnum *g, const realnum *g1, const realnum *g2, const realnum *u,
                      const realnum *u1, const realnum *u2, ptrdiff_t s, ptrdiff_t s1, ptrdiff_t s2,
                      const realnum *chi2, const realnum *chi3, realnum *fw, direction dsigw,
                      const realnum *sigw, const realnum *kapw);



void step_update_EDHB_NL(realnum *f, realnum *f_2, realnum *f_3, component fc,const grid_volume &gv, const ivec is, const ivec is_2, const ivec is_3, const ivec ie,const realnum *g, const realnum *g1, const realnum *g2, const realnum *u,const realnum *u_2, const realnum *u_3, const realnum *u1, const realnum *u2, ptrdiff_t s, ptrdiff_t s1, ptrdiff_t s2,const realnum *chi2, const realnum *chi3, realnum *fw, realnum *fw_2_atZ, realnum *f2_3_atZ, realnum *fw_2, realnum *fw_3,direction dsigw, direction dsigw_2,direction dsigw_3,const realnum *sigw,const realnum *sigw_2, const realnum *sigw_3, const realnum *kapw, const realnum *kapw_2, const realnum *kapw_3);

void step_beta(realnum *f, component c, const realnum *g, const grid_volume &gv, const ivec is,
               const ivec ie, realnum betadt, direction dsig, const realnum *siginv, realnum *fu,
               direction dsigu, const realnum *siginvu, const realnum *cndinv, realnum *fcnd);

void step_bfast(realnum *f, component c, const realnum *g1, const realnum *g2, ptrdiff_t s1,
                ptrdiff_t s2, // strides for g1/g2 shift
                const grid_volume &gv, const ivec is, const ivec ie, realnum dtdx, direction dsig,
                const realnum *sig, const realnum *kap, const realnum *siginv, realnum *fu,
                direction dsigu, const realnum *sigu, const realnum *kapu, const realnum *siginvu,
                realnum dt, const realnum *cnd, const realnum *cndinv, realnum *fcnd, realnum *F,
                realnum k1, realnum k2);

// functions in step_generic_stride1.cpp, generated from step_generic.cpp:

void step_curl_stride1(realnum *f, component c, const realnum *g1, const realnum *g2, ptrdiff_t s1,
                       ptrdiff_t s2, // strides for g1/g2 shift
                       const grid_volume &gv, const ivec is, const ivec ie, realnum dtdx,
                       direction dsig, const realnum *sig, const realnum *kap,
                       const realnum *siginv, realnum *fu, direction dsigu, const realnum *sigu,
                       const realnum *kapu, const realnum *siginvu, realnum dt, const realnum *cnd,
                       const realnum *cndinv, realnum *fcnd);

void step_update_EDHB_stride1(realnum *f, component fc, const grid_volume &gv, const ivec is,
                              const ivec ie, const realnum *g, const realnum *g1, const realnum *g2,
                              const realnum *u, const realnum *u1, const realnum *u2, ptrdiff_t s,
                              ptrdiff_t s1, ptrdiff_t s2, const realnum *chi2, const realnum *chi3,
                              realnum *fw, direction dsigw, const realnum *sigw,
                              const realnum *kapw);

void step_update_EDHB_NL_stride1(realnum *f, realnum *f_2, realnum *f_3, component fc,const grid_volume &gv, const ivec is, const ivec is_2,const ivec is_3, const ivec ie, const realnum *g, const realnum *g1,const realnum *g2, const realnum *u, const realnum *u_2,const realnum *u_3, const realnum *u1, const realnum *u2, ptrdiff_t s,ptrdiff_t s1, ptrdiff_t s2, const realnum *chi2, const realnum *chi3,realnum *fw, realnum *fw_2_atZ, realnum *f2_3_atZ, realnum *fw_2,realnum *fw_3, direction dsigw, direction dsigw_2, direction dsigw_3,const realnum *sigw, const realnum *sigw_2, const realnum *sigw_3,const realnum *kapw, const realnum *kapw_2, const realnum *kapw_3);

void step_beta_stride1(realnum *f, component c, const realnum *g, const grid_volume &gv,
                       const ivec is, const ivec ie, realnum betadt, direction dsig,
                       const realnum *siginv, realnum *fu, direction dsigu, const realnum *siginvu,
                       const realnum *cndinv, realnum *fcnd);

void step_bfast_stride1(realnum *f, component c, const realnum *g1, const realnum *g2, ptrdiff_t s1,
                        ptrdiff_t s2, // strides for g1/g2 shift
                        const grid_volume &gv, const ivec is, const ivec ie, realnum dtdx,
                        direction dsig, const realnum *sig, const realnum *kap,
                        const realnum *siginv, realnum *fu, direction dsigu, const realnum *sigu,
                        const realnum *kapu, const realnum *siginvu, realnum dt, const realnum *cnd,
                        const realnum *cndinv, realnum *fcnd, realnum *F, realnum k1, realnum k2);

/* macro wrappers around time-stepping functions: for performance reasons,
   if the inner loop is stride-1 then we use the stride-1 versions,
   which allow gcc (and possibly other compilers) to do additional
   optimizations, especially loop vectorization */

#define STEP_CURL(f, c, g1, g2, s1, s2, gv, is, ie, dtdx, dsig, sig, kap, siginv, fu, dsigu, sigu, \
                  kapu, siginvu, dt, cnd, cndinv, fcnd)                                            \
  do {                                                                                             \
    if (LOOPS_ARE_STRIDE1(gv))                                                                     \
      step_curl_stride1(f, c, g1, g2, s1, s2, gv, is, ie, dtdx, dsig, sig, kap, siginv, fu, dsigu, \
                        sigu, kapu, siginvu, dt, cnd, cndinv, fcnd);                               \
    else                                                                                           \
      step_curl(f, c, g1, g2, s1, s2, gv, is, ie, dtdx, dsig, sig, kap, siginv, fu, dsigu, sigu,   \
                kapu, siginvu, dt, cnd, cndinv, fcnd);                                             \
  } while (0)


#define STEP_UPDATE_EDHB(f, fc, gv, is, ie, g, g1, g2, u, u1, u2, s, s1, s2, chi2, chi3, fw,       \
                         dsigw, sigw, kapw)                                                        \
  do {                                                                                             \
    if (LOOPS_ARE_STRIDE1(gv))                                                                     \
      step_update_EDHB_stride1(f, fc, gv, is, ie, g, g1, g2, u, u1, u2, s, s1, s2, chi2, chi3, fw, \
                               dsigw, sigw, kapw);                                                 \
    else                                                                                           \
      step_update_EDHB(f, fc, gv, is, ie, g, g1, g2, u, u1, u2, s, s1, s2, chi2, chi3, fw, dsigw,  \
                       sigw, kapw);                                                                \
  } while (0)




// TODO don't think if(LOOPS_ARE_STRIDE1(gv)) returns as true as can't see stride1 vn defined anywhere - confirm
#define STEP_UPDATE_EDHB_NL(f, f2, f3, fc, gv, is, is2, is3, ie, g, g1, g2, u, u_2, u_3, u1, u2, s, s1, s2, chi2, chi3, fw, fTemp2, fTemp3, fw2, fw3,       \
                         dsigw, dsigw2, dsigw3, sigw, sigw2, sigw3, kapw, kapw2, kapw3)                                                        \
  do {                                                                                             \
    if (LOOPS_ARE_STRIDE1(gv))                                                                     \ 
      step_update_EDHB_NL_stride1(f, f2, f3, fc, gv, is, is2, is3, ie, g, g1, g2, u, u_2, u_3, u1, u2, s, s1, s2, chi2, chi3, fw, fTemp2, fTemp3, fw2, fw3, dsigw,  \
                       dsigw2, dsigw3, sigw, sigw2, sigw3, kapw, kapw2, kapw3);                                                                \
    else                                                                                           \
      step_update_EDHB_NL(f, f2, f3, fc, gv, is, is2, is3, ie, g, g1, g2, u, u_2, u_3, u1, u2, s,  \
                          s1, s2, chi2, chi3, fw, fTemp2, fTemp3, fw2, fw3, dsigw,  \
                       dsigw2, dsigw3, sigw, sigw2, sigw3, kapw, kapw2, kapw3);                                                                \
  } while (0)

#define STEP_BETA(f, c, g, gv, is, ie, betadt, dsig, siginv, fu, dsigu, siginvu, cndinv, fcnd)     \
  do {                                                                                             \
    if (LOOPS_ARE_STRIDE1(gv))                                                                     \
      step_beta_stride1(f, c, g, gv, is, ie, betadt, dsig, siginv, fu, dsigu, siginvu, cndinv,     \
                        fcnd);                                                                     \
    else                                                                                           \
      step_beta(f, c, g, gv, is, ie, betadt, dsig, siginv, fu, dsigu, siginvu, cndinv, fcnd);      \
  } while (0)

#define STEP_BFAST(f, c, g1, g2, s1, s2, gv, is, ie, dtdx, dsig, sig, kap, siginv, fu, dsigu,      \
                   sigu, kapu, siginvu, dt, cnd, cndinv, fcnd, F, k1, k2)                          \
  do {                                                                                             \
    if (LOOPS_ARE_STRIDE1(gv))                                                                     \
      step_bfast_stride1(f, c, g1, g2, s1, s2, gv, is, ie, dtdx, dsig, sig, kap, siginv, fu,       \
                         dsigu, sigu, kapu, siginvu, dt, cnd, cndinv, fcnd, F, k1, k2);            \
    else                                                                                           \
      step_bfast(f, c, g1, g2, s1, s2, gv, is, ie, dtdx, dsig, sig, kap, siginv, fu, dsigu, sigu,  \
                 kapu, siginvu, dt, cnd, cndinv, fcnd, F, k1, k2);                                 \
  } while (0)

// analytical Green's functions from near2far.cpp, which we might want to expose someday
void green3d(std::complex<double> *EH, const vec &x, double freq, double eps, double mu,
             const vec &x0, component c0, std::complex<double> f0);
void green2d(std::complex<double> *EH, const vec &x, double freq, double eps, double mu,
             const vec &x0, component c0, std::complex<double> f0);
void greencyl(std::complex<double> *EH, const vec &x, double freq, double eps, double mu,
              const vec &x0, component c0, std::complex<double> f0, double m, double tol);

// functions in array_slice.cpp:

complex<realnum> *collapse_array(complex<realnum> *array, int *rank, size_t dims[3],
                                 direction dirs[3], volume where);

void reduce_array_dimensions(volume where, int full_rank, size_t dims[3], direction dirs[3],
                             size_t stride[3], int &reduced_rank, size_t reduced_dims[3],
                             direction reduced_dirs[3], size_t reduced_stride[3]);

/**** macros used internally by step_generic and friends: ****/
/* These macros get into the guts of the PLOOP_OVER_VOL loops to
   efficiently construct the index k into a PML sigma array.
   Basically, k needs to increment by 2 for each increment of one of
   LOOP's for-loops, starting at the appropriate corner of the grid_volume,
   and these macros define the relevant strides etc. for each loop.
   KSTRIDE_DEF defines the relevant strides etc. and goes outside the
   LOOP, wheras KDEF defines the k index and goes inside the LOOP. */
#define KSTRIDE_DEF(dsig, k, is, gv)                                                               \
  const int k##0 = is.in_direction(dsig) - gv.little_corner().in_direction(dsig);                  \
  const int s##k##1 = gv.yucky_direction(0) == dsig ? 2 : 0;                                       \
  const int s##k##2 = gv.yucky_direction(1) == dsig ? 2 : 0;                                       \
  const int s##k##3 = gv.yucky_direction(2) == dsig ? 2 : 0
#define KDEF(k, dsig)                                                                              \
  const int k = ((k##0 + s##k##1 * loop_i1) + s##k##2 * loop_i2) + s##k##3 * loop_i3
#define DEF_k KDEF(k, dsig)
#define DEF_ku KDEF(ku, dsigu)
#define DEF_kw KDEF(kw, dsigw)
#define DEF_kw_2 KDEF(kw_2, dsigw_2)
#define DEF_kw_3 KDEF(kw_3, dsigw_3)

} // namespace meep
