/*
 * QR0_private.h
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

#ifndef RTW_HEADER_QR0_private_h_
#define RTW_HEADER_QR0_private_h_
#include "rtwtypes.h"
#include "multiword_types.h"
#include "QR0_types.h"
#include "QR0.h"
#include "rtw_continuous.h"
#include "rtw_solver.h"

extern real_T rt_hypotd(real_T u0, real_T u1);
extern int32_T div_nde_s32_floor(int32_T numerator, int32_T denominator);

#endif                                 /* RTW_HEADER_QR0_private_h_ */
