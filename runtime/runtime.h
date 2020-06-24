#ifndef _PLASMA_RUNTIME_H_
#define _PLASMA_RUNTIME_H_

#include "common.h"
#include "lapacke.h"
#include <omp.h>

#if defined(PLASMA_WITH_CUDA_HYBRID) || defined(PLASMA_WITH_CUDA_PURE)
#include "cublas_v2.h"
#endif

#define RT_DEBUG 1
#define OUTFILE stderr

#endif //_PLASMA_RUNTIME_H_
