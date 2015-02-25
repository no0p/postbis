/*-------------------------------------------------------------------------
*
* Copyright (c) 2013, Max Planck Institute for Marine Microbiology
*
* This software is released under the PostgreSQL License
*
* Author: Michael Schneider <mschneid@mpi-bremen.de>
*
* IDENTIFICATION
*   include/sequence/compression.h
*
*-------------------------------------------------------------------------
*/
#ifndef SEQUENCE_COMPRESSION_H_
#define SEQUENCE_COMPRESSION_H_

#include "postgres.h"

#include "sequence/sequence.h"

/**
 * get_compressed_size()
 * 		Compute the size of a compressed sequence.
 *
 * 	PB_SequenceInfo* info : stats of input sequence
 * 	PB_CodeSet* codeset : code to test
 */
uint32 get_compressed_size(const PB_SequenceInfo* info,
						   PB_CodeSet* codeset);

/**
 * encode()
 * 		Encode a sequence.
 *
 * 	uint8* input : input sequence, not null-terminated
 * 					length must be in info
 * 	uint32 compressed_size : length of the compressed sequence, calculated
 * 							with get_compressed_size()
 * 	PB_CodeSet* codeset : codeset for encoding
 * 	PB_SequenceInfo* info : info about the sequence to compress
 */
PB_CompressedSequence* encode(uint8* input,
							  uint32 compressed_size,
							  PB_CodeSet* codeset,
							  PB_SequenceInfo* info);

/**
 * decode()
 * 		Decode a compressed sequence.
 *
 * 	Varlena* input : pointer to non-detoasted compressed sequence.
 * 	uint8* output : pointer to sufficient space to store the decoded sequence
 * 	uint32 start_position : position to start decoding from, first is 0
 * 	uint32 out_length : number of characters to decode
 * 	PB_CodeSet** codeset : list of fixed codesets
 */
void decode(Varlena* input,
			uint8* output,
			uint32 start_position,
			uint32 out_length,
			PB_CodeSet** fixed_codesets);

#endif /* SEQUENCE_COMPRESSION_H_ */
