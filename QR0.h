/*
 * QR0.h
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

#ifndef RTW_HEADER_QR0_h_
#define RTW_HEADER_QR0_h_
#ifndef QR0_COMMON_INCLUDES_
#define QR0_COMMON_INCLUDES_
#include "rtwtypes.h"
#include "rtw_continuous.h"
#include "rtw_solver.h"
#endif                                 /* QR0_COMMON_INCLUDES_ */

#include "QR0_types.h"
#include <stddef.h>
#include <string.h>

/* Macros for accessing real-time model data structure */
#ifndef rtmGetErrorStatus
#define rtmGetErrorStatus(rtm)         ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
#define rtmSetErrorStatus(rtm, val)    ((rtm)->errorStatus = (val))
#endif

/* Block signals (default storage) */
typedef struct {
  creal_T b_A[1040];
  creal_T b_A_m[1040];
} B_QR0_T;

/* Block states (default storage) for system '<Root>' */
typedef struct {
  dsp_simulink_QRFactorization__T obj; /* '<S1>/QR Factor with pivot' */
  boolean_T objisempty;                /* '<S1>/QR Factor with pivot' */
} DW_QR0_T;

/* External inputs (root inport signals with default storage) */
typedef struct {
  creal_T A[1040];                     /* '<Root>/A' */
} ExtU_QR0_T;

/* External outputs (root outports fed by signals with default storage) */
typedef struct {
  creal_T Out1[16];                    /* '<Root>/Out1' */
  real_T PIVOTPERM[4];                 /* '<Root>/PIVOTPERM' */
} ExtY_QR0_T;

/* Real-time Model Data Structure */
struct tag_RTM_QR0_T {
  const char_T *errorStatus;
};

/* Block signals (default storage) */
extern B_QR0_T QR0_B;

/* Block states (default storage) */
extern DW_QR0_T QR0_DW;

/* External inputs (root inport signals with default storage) */
extern ExtU_QR0_T QR0_U;

/* External outputs (root outports fed by signals with default storage) */
extern ExtY_QR0_T QR0_Y;

/* Model entry point functions */
extern void QR0_initialize(void);
extern void QR0_step(void);
extern void QR0_terminate(void);

/* Real-time Model object */
extern RT_MODEL_QR0_T *const QR0_M;

/*-
 * The generated code includes comments that allow you to trace directly
 * back to the appropriate location in the model.  The basic format
 * is <system>/block_name, where system is the system number (uniquely
 * assigned by Simulink) and block_name is the name of the block.
 *
 * Note that this particular code originates from a subsystem build,
 * and has its own system numbers different from the parent model.
 * Refer to the system hierarchy for this subsystem below, and use the
 * MATLAB hilite_system command to trace the generated code back
 * to the parent model.  For example,
 *
 * hilite_system('Fixed_point_conveted_LCMV_2/Time Delayed LCMV Beamformer/lcmv weigth calculator simplified/qrlinsolve/QR Factorization ')    - opens subsystem Fixed_point_conveted_LCMV_2/Time Delayed LCMV Beamformer/lcmv weigth calculator simplified/qrlinsolve/QR Factorization
 * hilite_system('Fixed_point_conveted_LCMV_2/Time Delayed LCMV Beamformer/lcmv weigth calculator simplified/qrlinsolve/QR Factorization /Kp') - opens and selects block Kp
 *
 * Here is the system hierarchy for this model
 *
 * '<Root>' : 'Fixed_point_conveted_LCMV_2/Time Delayed LCMV Beamformer/lcmv weigth calculator simplified/qrlinsolve'
 * '<S1>'   : 'Fixed_point_conveted_LCMV_2/Time Delayed LCMV Beamformer/lcmv weigth calculator simplified/qrlinsolve/QR Factorization '
 */
#endif                                 /* RTW_HEADER_QR0_h_ */
