/*-------------------------------------------------------------------------
*
* Copyright (c) 2013, Max Planck Institute for Marine Microbiology
*
* This software is released under the PostgreSQL License
*
* Author: Michael Schneider <mschneid@mpi-bremen.de>
*
* IDENTIFICATION
*   include/sequence/code_set_creation.h
*
*-------------------------------------------------------------------------
*/

#ifndef SEQUENCE_CODE_SET_CREATION_H_
#define SEQUENCE_CODE_SET_CREATION_H_

#include "sequence/sequence.h"

/**
 * get_equal_lengths_code()
 * 		Creates a code set where all codewords have equal lengths.
 *
 * 	RLE statistics will be ignored.
 *
 * 	PB_SequenceInfo* info : info about sequence at hand
 */
PB_CodeSet* get_equal_lengths_code(const PB_SequenceInfo* info);

/**
 * get_huffman_code()
 * 		Build huffman code from given sequence info.
 *
 * If the huffman tree, that is created during this process, is deeper
 * than the size of PB_PrefixCode allows the function returns NULL.
 *
 * 	RLE statistics will be ignored.
 *
 * 	PB_SequenceInfo* info : stats of the sequence
 */
PB_CodeSet* get_huffman_code(const PB_SequenceInfo* info);

/**
 * get_huffman_code_rle()
 * 		Build huffman code from given rle sequence info.
 *
 * If the huffman tree, that is created during this process, is deeper
 * than the size of PB_PrefixCode allows the function returns NULL.
 *
 * 	PB_SequenceInfo* info : stats of the sequence
 */
PB_CodeSet* get_huffman_code_rle(const PB_SequenceInfo* info);

/**
 * truncate_huffman_tree()
 * 		Truncate huffman code if useful.
 *
 * 	If it is advantageous a truncated huffman tree will returned. If
 * 	not the result is NULL.
 *
 * 	PB_CodeSet* codeset : original huffman code created from stats
 * 	PB_SequenceInfo* info : stats of input sequence
 */
PB_CodeSet* truncate_huffman_code(const PB_CodeSet* codeset,
								  const PB_SequenceInfo* info);

/**
 * get_optimal_code()
 * 		Creates an optimal code for a given sequence.
 *
 * 	PB_SequenceInfo* info : info about sequence at hand
 */
PB_CodeSet* get_optimal_code(const PB_SequenceInfo* info);

#endif /* SEQUENCE_CODE_SET_CREATION_H_ */
