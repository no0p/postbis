/*-------------------------------------------------------------------------
*
* Copyright (c) 2013, Max Planck Institute for Marine Microbiology
*
* This software is released under the PostgreSQL License
*
* Author: Michael Schneider <mschneid@mpi-bremen.de>
*
* IDENTIFICATION
*   include/sequence/functions.h
*
*-------------------------------------------------------------------------
*/
#ifndef SEQUENCE_FUNCTIONS_H_
#define SEQUENCE_FUNCTIONS_H_

#include "sequence/sequence.h"

/**
 * reverse()
 * 		Reverses a compressed sequence.
 */
PB_CompressedSequence* reverse(PB_CompressedSequence* sequence, PB_CodeSet** fixed_codesets);

/**
 * sequence_equal()
 * 		Compares two compressed sequences. Returns (-1) if equal, 0 if not.
 *
 * 	Varlena* seq1 : first possibly toasted sequence
 * 	Varlena* seq2 : second possibly toasted sequence
 */
bool sequence_equal(Varlena* raw_seq1, Varlena* raw_seq2, PB_CodeSet** fixed_codesets);

/**
 * sequence_compare()
 * 		Compares two compressed sequences. Returns (-1) if first is smaller,
 * 		0 if both are equal, 1 if second sequence is smaller (lexicographically).
 *
 * 	Varlena* seq_a : first possibly toasted sequence
 * 	Varlena* seq_b : second possibly toasted sequence
 */
int sequence_compare(Varlena* seq_a, Varlena* seq_b, PB_CodeSet** fixed_codesets);

/*
 * sequence_crc32()
 *		Compute CRC32 of a compressed sequence.
 * 		This code implements the AUTODIN II polynomial
 * 		Original code by Spencer Garrett <srg@quick.com>
 * 		Taken into PostBIS from pgsql/src/contrib/hstore/crc32.c
 *
 * 	PB_CompressedSequence* seq: detoasted input sequence
 * 	PB_CodeSet** fixed_codesets : fixed codes
 */
uint32 sequence_crc32(PB_CompressedSequence* seq, PB_CodeSet** fixed_codesets);

/*
 * sequence_strpos()
 * 		Find position of given string.
 *
 * 	PB_CompressedSequence* seq : detoasted sequence to search in
 * 	Text* search : detoasted text to search for
 * 	PB_CodeSet** fixed_codesets : fixed codes
 */
uint32 sequence_strpos(PB_CompressedSequence* seq, text* search, PB_CodeSet** fixed_codesets);

#endif /* SEQUENCE_FUNCTIONS_H_ */
