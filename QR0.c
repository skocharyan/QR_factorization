/*
 * QR0.c
 *
 * Code generation for model "QR0".
 *
 * Model version              : 1.312
 * Simulink Coder version : 23.2 (R2023b) 01-Aug-2023
 * C source code generated on : Tue Nov 21 14:00:25 2023
 *
 * Target selection: grt.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: ASIC/FPGA->ASIC/FPGA
 * Emulation hardware selection:
 *    Differs from embedded hardware (Custom Processor->MATLAB Host Computer)
 * Code generation objective: Execution efficiency
 * Validation result: Not run
 */

#include "QR0.h"
#include "rtwtypes.h"
#include <string.h>
#include <math.h>
#include "QR0_private.h"
#include <emmintrin.h>

/* Block signals (default storage) */
B_QR0_T QR0_B;

/* Block states (default storage) */
DW_QR0_T QR0_DW;

/* External inputs (root inport signals with default storage) */
ExtU_QR0_T QR0_U;

/* External outputs (root outports fed by signals with default storage) */
ExtY_QR0_T QR0_Y;

/* Real-time model */
static RT_MODEL_QR0_T QR0_M_;
RT_MODEL_QR0_T *const QR0_M = &QR0_M_;

/* Forward declaration for local functions */
static real_T QR0_xnrm2(int32_T n, const creal_T x[1040], int32_T ix0);
static real_T QR0_xdlapy3(real_T x1, real_T x2, real_T x3);
static creal_T QR0_recip(const creal_T y);
static void QR0_xzlarf_c(int32_T m, int32_T n, int32_T iv0, const creal_T tau,
  const creal_T C[1040], int32_T ic0, creal_T work[4], creal_T b_C[1040]);
static void QR0_xgeqp3(const creal_T A[1040], creal_T b_A[1040], creal_T tau[4],
  int32_T jpvt[4]);
int32_T div_nde_s32_floor(int32_T numerator, int32_T denominator)
{
  return (((numerator < 0) != (denominator < 0)) && (numerator % denominator !=
           0) ? -1 : 0) + numerator / denominator;
}

real_T rt_hypotd(real_T u0, real_T u1)
{
  real_T a;
  real_T b;
  real_T y;
  a = fabs(u0);
  b = fabs(u1);
  if (a < b) {
    a /= b;
    y = sqrt(a * a + 1.0) * b;
  } else if (a > b) {
    b /= a;
    y = sqrt(b * b + 1.0) * a;
  } else {
    y = a * 1.4142135623730951;
  }

  return y;
}

static real_T QR0_xnrm2(int32_T n, const creal_T x[1040], int32_T ix0)
{
  real_T absxk;
  real_T scale;
  real_T t;
  real_T y;
  int32_T k;
  int32_T kend;
  y = 0.0;
  if (n < 1) {
  } else if (n == 1) {
    y = rt_hypotd(x[ix0 - 1].re, x[ix0 - 1].im);
  } else {
    scale = 3.3121686421112381E-170;
    kend = ix0 + n;
    for (k = ix0; k < kend; k++) {
      absxk = fabs(x[k - 1].re);
      if (absxk > scale) {
        t = scale / absxk;
        y = y * t * t + 1.0;
        scale = absxk;
      } else {
        t = absxk / scale;
        y += t * t;
      }

      absxk = fabs(x[k - 1].im);
      if (absxk > scale) {
        t = scale / absxk;
        y = y * t * t + 1.0;
        scale = absxk;
      } else {
        t = absxk / scale;
        y += t * t;
      }
    }

    y = scale * sqrt(y);
  }

  return y;
}

static real_T QR0_xdlapy3(real_T x1, real_T x2, real_T x3)
{
  real_T a;
  real_T b;
  real_T c;
  real_T y;
  a = fabs(x1);
  b = fabs(x2);
  c = fabs(x3);
  y = fmax(a, b);
  if (c > y) {
    y = c;
  }

  if (y > 0.0) {
    a /= y;
    b /= y;
    c /= y;
    y *= sqrt((a * a + c * c) + b * b);
  } else {
    y = (a + b) + c;
  }

  return y;
}

static creal_T QR0_recip(const creal_T y)
{
  creal_T z;
  real_T bim;
  real_T brm;
  brm = fabs(y.re);
  bim = fabs(y.im);
  if (y.im == 0.0) {
    z.re = 1.0 / y.re;
    z.im = 0.0;
  } else if (y.re == 0.0) {
    z.re = 0.0;
    z.im = -1.0 / y.im;
  } else if (brm > bim) {
    brm = y.im / y.re;
    _mm_storeu_pd((real_T *)&z, _mm_div_pd(_mm_set_pd(-brm, 1.0), _mm_set1_pd
      (brm * y.im + y.re)));
  } else if (brm == bim) {
    real_T br;
    bim = 0.5;
    if (y.re < 0.0) {
      bim = -0.5;
    }

    br = 0.5;
    if (y.im < 0.0) {
      br = -0.5;
    }

    _mm_storeu_pd((real_T *)&z, _mm_div_pd(_mm_set_pd(-br, bim), _mm_set1_pd(brm)));
  } else {
    brm = y.re / y.im;
    _mm_storeu_pd((real_T *)&z, _mm_div_pd(_mm_set_pd(-1.0, brm), _mm_set1_pd
      (brm * y.re + y.im)));
  }

  return z;
}

static void QR0_xzlarf_c(int32_T m, int32_T n, int32_T iv0, const creal_T tau,
  const creal_T C[1040], int32_T ic0, creal_T work[4], creal_T b_C[1040])
{
  real_T tmp_0[2];
  int32_T coltop;
  int32_T iac;
  int32_T ijA;
  int32_T jA;
  int32_T lastc;
  int32_T lastv;
  memcpy(&b_C[0], &C[0], 1040U * sizeof(creal_T));
  if ((tau.re != 0.0) || (tau.im != 0.0)) {
    boolean_T exitg2;
    lastv = m;
    lastc = (iv0 + m) - 2;
    while ((lastv > 0) && ((C[lastc].re == 0.0) && (C[lastc].im == 0.0))) {
      lastv--;
      lastc--;
    }

    lastc = n - 1;
    exitg2 = false;
    while ((!exitg2) && (lastc + 1 > 0)) {
      int32_T exitg1;
      coltop = lastc * 260 + ic0;
      jA = coltop;
      do {
        exitg1 = 0;
        if (jA <= (coltop + lastv) - 1) {
          if ((C[jA - 1].re != 0.0) || (C[jA - 1].im != 0.0)) {
            exitg1 = 1;
          } else {
            jA++;
          }
        } else {
          lastc--;
          exitg1 = 2;
        }
      } while (exitg1 == 0);

      if (exitg1 == 1) {
        exitg2 = true;
      }
    }
  } else {
    lastv = 0;
    lastc = -1;
  }

  if (lastv > 0) {
    real_T alpha1_im;
    real_T alpha1_re;
    real_T c_im;
    real_T c_re;
    int32_T d;
    if (lastc + 1 != 0) {
      if (lastc >= 0) {
        memset(&work[0], 0, (uint32_T)(lastc + 1) * sizeof(creal_T));
      }

      coltop = 260 * lastc + ic0;
      for (iac = ic0; iac <= coltop; iac += 260) {
        __m128d tmp;
        c_re = 0.0;
        c_im = 0.0;
        d = iac + lastv;
        for (jA = iac; jA < d; jA++) {
          ijA = ((iv0 + jA) - iac) - 1;
          tmp = _mm_add_pd(_mm_add_pd(_mm_mul_pd(_mm_loadu_pd((const real_T *)
            &b_C[ijA]), _mm_set1_pd(b_C[jA - 1].re)), _mm_mul_pd(_mm_mul_pd
            (_mm_shuffle_pd(_mm_loadu_pd((const real_T *)&b_C[ijA]),
                            _mm_loadu_pd((const real_T *)&b_C[ijA]), 1),
             _mm_set1_pd(b_C[jA - 1].im)), _mm_set_pd(-1.0, 1.0))), _mm_set_pd
                           (c_im, c_re));
          _mm_storeu_pd(&tmp_0[0], tmp);
          c_re = tmp_0[0];
          c_im = tmp_0[1];
        }

        ijA = div_nde_s32_floor(iac - ic0, 260);
        tmp = _mm_add_pd(_mm_loadu_pd((const real_T *)&work[ijA]), _mm_set_pd
                         (c_im, c_re));
        _mm_storeu_pd((real_T *)&work[ijA], tmp);
      }
    }

    alpha1_re = -tau.re;
    alpha1_im = -tau.im;
    if ((-tau.re != 0.0) || (-tau.im != 0.0)) {
      jA = ic0;
      for (iac = 0; iac <= lastc; iac++) {
        real_T work_im;
        c_im = work[iac].re;
        work_im = work[iac].im;
        if ((c_im != 0.0) || (work_im != 0.0)) {
          c_re = c_im * alpha1_re + work_im * alpha1_im;
          c_im = c_im * alpha1_im - work_im * alpha1_re;
          d = (lastv + jA) - 1;
          for (ijA = jA; ijA <= d; ijA++) {
            real_T b_C_im_tmp;
            coltop = ((iv0 + ijA) - jA) - 1;
            work_im = b_C[coltop].re;
            b_C_im_tmp = b_C[coltop].im;
            b_C[ijA - 1].re += work_im * c_re - b_C_im_tmp * c_im;
            b_C[ijA - 1].im += work_im * c_im + b_C_im_tmp * c_re;
          }
        }

        jA += 260;
      }
    }
  }
}

static void QR0_xgeqp3(const creal_T A[1040], creal_T b_A[1040], creal_T tau[4],
  int32_T jpvt[4])
{
  __m128d tmp_0;
  creal_T work[4];
  creal_T b_A_0;
  creal_T c_atmp;
  creal_T tau_0;
  creal_T tmp;
  real_T vn1[4];
  real_T vn2[4];
  real_T absxk;
  real_T b_A_1;
  real_T beta1;
  real_T scale;
  real_T t;
  int32_T b_A_tmp;
  int32_T b_k;
  int32_T itemp;
  int32_T ix;
  int32_T iy;
  int32_T j;
  int32_T kend;
  int32_T nmip1;
  int32_T temp_re_tmp;
  memcpy(&b_A[0], &A[0], 1040U * sizeof(creal_T));
  memset(&tau[0], 0, sizeof(creal_T) << 2U);
  work[0].re = 0.0;
  work[0].im = 0.0;
  work[1].re = 0.0;
  work[1].im = 0.0;
  work[2].re = 0.0;
  work[2].im = 0.0;
  work[3].re = 0.0;
  work[3].im = 0.0;
  for (b_k = 0; b_k < 4; b_k++) {
    jpvt[b_k] = b_k + 1;
    ix = b_k * 260;
    beta1 = 0.0;
    scale = 3.3121686421112381E-170;
    for (itemp = ix + 1; itemp <= ix + 260; itemp++) {
      absxk = fabs(A[itemp - 1].re);
      if (absxk > scale) {
        t = scale / absxk;
        beta1 = beta1 * t * t + 1.0;
        scale = absxk;
      } else {
        t = absxk / scale;
        beta1 += t * t;
      }

      absxk = fabs(A[itemp - 1].im);
      if (absxk > scale) {
        t = scale / absxk;
        beta1 = beta1 * t * t + 1.0;
        scale = absxk;
      } else {
        t = absxk / scale;
        beta1 += t * t;
      }
    }

    absxk = scale * sqrt(beta1);
    vn1[b_k] = absxk;
    vn2[b_k] = absxk;
  }

  for (b_k = 0; b_k < 4; b_k++) {
    j = b_k + 1;
    kend = b_k * 260 + b_k;
    nmip1 = 4 - b_k;
    iy = 0;
    if (4 - b_k > 1) {
      scale = vn1[b_k];
      for (itemp = 2; itemp <= nmip1; itemp++) {
        absxk = vn1[(b_k + itemp) - 1];
        if (absxk > scale) {
          iy = itemp - 1;
          scale = absxk;
        }
      }
    }

    nmip1 = b_k + iy;
    if (nmip1 != b_k) {
      ix = nmip1 * 260;
      iy = b_k * 260;
      for (itemp = 0; itemp < 260; itemp++) {
        temp_re_tmp = ix + itemp;
        scale = b_A[temp_re_tmp].re;
        absxk = b_A[temp_re_tmp].im;
        b_A_tmp = iy + itemp;
        b_A[temp_re_tmp] = b_A[b_A_tmp];
        b_A[b_A_tmp].re = scale;
        b_A[b_A_tmp].im = absxk;
      }

      itemp = jpvt[nmip1];
      jpvt[nmip1] = jpvt[b_k];
      jpvt[b_k] = itemp;
      vn1[nmip1] = vn1[b_k];
      vn2[nmip1] = vn2[b_k];
    }

    ix = kend + 2;
    scale = b_A[kend].re;
    absxk = b_A[kend].im;
    tau[b_k].re = 0.0;
    tau[b_k].im = 0.0;
    beta1 = QR0_xnrm2(259 - b_k, b_A, kend + 2);
    if ((beta1 != 0.0) || (b_A[kend].im != 0.0)) {
      beta1 = QR0_xdlapy3(scale, absxk, beta1);
      if (scale >= 0.0) {
        beta1 = -beta1;
      }

      if (fabs(beta1) < 1.0020841800044864E-292) {
        iy = -1;
        do {
          iy++;
          nmip1 = (kend - b_k) + 260;
          for (itemp = ix; itemp <= nmip1; itemp++) {
            tmp_0 = _mm_mul_pd(_mm_loadu_pd((const real_T *)&b_A[itemp - 1]),
                               _mm_set1_pd(9.9792015476736E+291));
            _mm_storeu_pd((real_T *)&b_A[itemp - 1], tmp_0);
          }

          beta1 *= 9.9792015476736E+291;
          scale *= 9.9792015476736E+291;
          absxk *= 9.9792015476736E+291;
        } while ((fabs(beta1) < 1.0020841800044864E-292) && (iy + 1 < 20));

        beta1 = QR0_xdlapy3(scale, absxk, QR0_xnrm2(259 - b_k, b_A, kend + 2));
        if (scale >= 0.0) {
          beta1 = -beta1;
        }

        t = beta1 - scale;
        if (0.0 - absxk == 0.0) {
          tau[b_k].re = t / beta1;
          tau[b_k].im = 0.0;
        } else if (t == 0.0) {
          tau[b_k].re = 0.0;
          tau[b_k].im = (0.0 - absxk) / beta1;
        } else {
          tau[b_k].re = t / beta1;
          tau[b_k].im = (0.0 - absxk) / beta1;
        }

        c_atmp.re = scale - beta1;
        c_atmp.im = absxk;
        tmp = QR0_recip(c_atmp);
        scale = tmp.re;
        absxk = tmp.im;
        for (itemp = ix; itemp <= nmip1; itemp++) {
          t = b_A[itemp - 1].im;
          b_A_1 = b_A[itemp - 1].re;
          b_A[itemp - 1].re = b_A_1 * scale - t * absxk;
          b_A[itemp - 1].im = t * scale + b_A_1 * absxk;
        }

        for (itemp = 0; itemp <= iy; itemp++) {
          beta1 *= 1.0020841800044864E-292;
        }

        scale = beta1;
        absxk = 0.0;
      } else {
        t = beta1 - scale;
        if (0.0 - absxk == 0.0) {
          tau[b_k].re = t / beta1;
          tau[b_k].im = 0.0;
        } else if (t == 0.0) {
          tau[b_k].re = 0.0;
          tau[b_k].im = (0.0 - absxk) / beta1;
        } else {
          tau[b_k].re = t / beta1;
          tau[b_k].im = (0.0 - absxk) / beta1;
        }

        b_A_0.re = scale - beta1;
        b_A_0.im = absxk;
        tmp = QR0_recip(b_A_0);
        scale = tmp.re;
        absxk = tmp.im;
        nmip1 = (kend - b_k) - 1;
        for (itemp = ix; itemp <= nmip1 + 261; itemp++) {
          t = b_A[itemp - 1].im;
          b_A_1 = b_A[itemp - 1].re;
          b_A[itemp - 1].re = b_A_1 * scale - t * absxk;
          b_A[itemp - 1].im = t * scale + b_A_1 * absxk;
        }

        scale = beta1;
        absxk = 0.0;
      }
    }

    b_A[kend].re = scale;
    b_A[kend].im = absxk;
    if (b_k + 1 < 4) {
      b_A[kend].re = 1.0;
      b_A[kend].im = 0.0;
      tau_0.re = tau[b_k].re;
      tau_0.im = -tau[b_k].im;
      memcpy(&QR0_B.b_A_m[0], &b_A[0], 1040U * sizeof(creal_T));
      QR0_xzlarf_c(260 - b_k, 3 - b_k, kend + 1, tau_0, QR0_B.b_A_m, kend + 261,
                   work, b_A);
      b_A[kend].re = scale;
      b_A[kend].im = absxk;
    }

    for (nmip1 = j + 1; nmip1 < 5; nmip1++) {
      itemp = ((nmip1 - 1) * 260 + b_k) + 1;
      absxk = vn1[nmip1 - 1];
      if (absxk != 0.0) {
        beta1 = rt_hypotd(b_A[itemp - 1].re, b_A[itemp - 1].im) / absxk;
        beta1 = 1.0 - beta1 * beta1;
        if (beta1 < 0.0) {
          beta1 = 0.0;
        }

        scale = absxk / vn2[nmip1 - 1];
        scale = scale * scale * beta1;
        if (scale <= 1.4901161193847656E-8) {
          ix = itemp;
          beta1 = 0.0;
          scale = 3.3121686421112381E-170;
          kend = (itemp - b_k) - 1;
          for (itemp = ix + 1; itemp <= kend + 260; itemp++) {
            absxk = fabs(b_A[itemp - 1].re);
            if (absxk > scale) {
              t = scale / absxk;
              beta1 = beta1 * t * t + 1.0;
              scale = absxk;
            } else {
              t = absxk / scale;
              beta1 += t * t;
            }

            absxk = fabs(b_A[itemp - 1].im);
            if (absxk > scale) {
              t = scale / absxk;
              beta1 = beta1 * t * t + 1.0;
              scale = absxk;
            } else {
              t = absxk / scale;
              beta1 += t * t;
            }
          }

          beta1 = scale * sqrt(beta1);
          vn1[nmip1 - 1] = beta1;
          vn2[nmip1 - 1] = beta1;
        } else {
          vn1[nmip1 - 1] = absxk * sqrt(beta1);
        }
      }
    }
  }
}

/* Model step function */
void QR0_step(void)
{
  creal_T b_R[16];
  creal_T tau[4];
  int32_T b_jpvt[4];
  int32_T b_i;
  int32_T b_j;
  int32_T tmp;

  /* MATLABSystem: '<S1>/QR Factor with pivot' incorporates:
   *  Inport: '<Root>/A'
   */
  QR0_xgeqp3(QR0_U.A, QR0_B.b_A, tau, b_jpvt);
  for (b_j = 0; b_j < 4; b_j++) {
    for (b_i = 0; b_i <= b_j; b_i++) {
      b_R[b_i + (b_j << 2)] = QR0_B.b_A[260 * b_j + b_i];
    }

    for (b_i = b_j + 2; b_i < 5; b_i++) {
      tmp = ((b_j << 2) + b_i) - 1;
      b_R[tmp].re = 0.0;
      b_R[tmp].im = 0.0;
    }
  }

  /* Outport: '<Root>/Out1' incorporates:
   *  MATLABSystem: '<S1>/QR Factor with pivot'
   */
  memcpy(&QR0_Y.Out1[0], &b_R[0], sizeof(creal_T) << 4U);

  /* Outport: '<Root>/PIVOTPERM' incorporates:
   *  MATLABSystem: '<S1>/QR Factor with pivot'
   */
  QR0_Y.PIVOTPERM[0] = b_jpvt[0];
  QR0_Y.PIVOTPERM[1] = b_jpvt[1];
  QR0_Y.PIVOTPERM[2] = b_jpvt[2];
  QR0_Y.PIVOTPERM[3] = b_jpvt[3];
}

/* Model initialize function */
void QR0_initialize(void)
{
  /* Registration code */

  /* initialize error status */
  rtmSetErrorStatus(QR0_M, (NULL));

  /* states (dwork) */
  (void) memset((void *)&QR0_DW, 0,
                sizeof(DW_QR0_T));

  /* external inputs */
  (void)memset(&QR0_U, 0, sizeof(ExtU_QR0_T));

  /* external outputs */
  (void)memset(&QR0_Y, 0, sizeof(ExtY_QR0_T));

  /* Start for MATLABSystem: '<S1>/QR Factor with pivot' */
  QR0_DW.objisempty = true;
  QR0_DW.obj.isInitialized = 1;
}

/* Model terminate function */
void QR0_terminate(void)
{
  /* (no terminate code required) */
}
