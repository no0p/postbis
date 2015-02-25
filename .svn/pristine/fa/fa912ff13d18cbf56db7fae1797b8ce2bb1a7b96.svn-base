/*-------------------------------------------------------------------------
*
* Copyright (c) 2013, Max Planck Institute for Marine Microbiology
*
* This software is released under the PostgreSQL License
*
* Author: Michael Schneider <mschneid@mpi-bremen.de>
*
* IDENTIFICATION
*   include/utils/debug.h
*
*-------------------------------------------------------------------------
*/
#ifndef UTILS_DEBUG_H_
#define UTILS_DEBUG_H_

#include "postgres.h"

/*
 * To compile debug messages, add this to makefile:
 * PG_CPPFLAGS+=-D DEBUG
 */
#ifdef DEBUG

/*
 * Important summarizing information
 * Rule of thumb: everything outside loops
 */
#define PB_DEBUG1(error_message) \
	ereport(DEBUG1,(error_message));

/*
 * Detailed step-by-step information
 * Rule of thumb: everything inside loops
 */
#define PB_DEBUG2(error_message) \
	ereport(DEBUG2,(error_message));

/*
 * En-/decoding step-by-step information
 * Rule of thum: Everything inside the coding loops
 */
#define PB_DEBUG3(error_message) \
	ereport(DEBUG3,(error_message));

/*
 * Function call trace
 */
#define PB_TRACE(error_message) \
	ereport(DEBUG2,(error_message));

#else

#define PB_DEBUG1(error_message)

#define PB_DEBUG2(error_message)

#define PB_DEBUG3(error_message)

#define PB_TRACE(error_message)

#endif

#endif /* UTILS_DEBUG_H_ */
