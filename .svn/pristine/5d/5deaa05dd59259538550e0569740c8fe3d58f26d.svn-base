/*-------------------------------------------------------------------------
*
* Copyright (c) 2013, Max Planck Institute for Marine Microbiology
*
* This software is released under the PostgreSQL License
*
* Author: Michael Schneider <mschneid@mpi-bremen.de>
*
* IDENTIFICATION
*   include/sequence/stats.h
*
*-------------------------------------------------------------------------
*/

#ifndef SEQUENCE_STATS_H_
#define SEQUENCE_STATS_H_

#include "postgres.h"

#include "sequence/sequence.h"

/*
 * collect_alphabet()
 * 		Retrieve the alphabet from the frequencies, order symbols according
 * 		to their frequency and builds bitmaps.
 *
 * 	uint32* frequencies : input parameter symbol frequencies
 * 	uint8* n_symbols : output parameter for the total number of symbols
 * 	uint8** symbols : output parameter for the symbols
 * 	uint64* ascii_bitmap_low : output parameter for ASCII bitmap from 0 to 63
 * 	uint64* ascii_bitmap_high : output parameter for ASCII bitmap from 64 to 127
 */
void collect_alphabet(uint32* frequencies,
					  uint8* n_symbols,
					  uint8** symbols,
					  uint64* ascii_bitmap_low,
					  uint64* ascii_bitmap_high);

/*
 * check_ascii()
 * 		Check if any non-ASCII symbols were observed in UTF-8 cstring.
 * 		The first character of a multi-byte UTF-8 character must be
 * 		somewhere between 194 and 244 by definition. Single-byte UTF-8
 * 		characters are fully ASCII-compatible, hence always below 128.
 * 		So a maximum of 128 different symbols can be in an alphabet.
 *
 * 	PB_SequenceInfo* sequence_info : info to check for ASCII compliance
 */
void check_ascii(PB_SequenceInfo* sequence_info);

/**
 * get_sequence_info_??? will treat upper and lower case characters as different.
 */
#define PB_SEQUENCE_INFO_CASE_SENSITIVE		1

/**
 * get_sequence_info_??? will treat upper and lower case characters as equal.
 */
#define PB_SEQUENCE_INFO_CASE_INSENSITIVE	0

/**
 * get_sequence_info_??? will collect information for RLE.
 */
#define PB_SEQUENCE_INFO_WITH_RLE			2

/**
 * get_sequence_info_??? will not collect information for RLE.
 */
#define PB_SEQUENCE_INFO_WITHOUT_RLE		0

/**
 * get_sequence_info_cstring()
 * 		Obtain sequence length, symbol frequencies and alphabet.
 *
 * 		char* input : a null-terminated input string
 * 		int mode : determines case-sensitivity and whether
 * 				   run-length stats should be collected.
 * 				   Possible modes:
 * 					PB_SEQUENCE_INFO_CASE_SENSITIVE
 * 					PB_SEQUENCE_INFO_CASE_INSENSITIVE
 * 					PB_SEQUENCE_INFO_WITH_RLE
 * 					PB_SEQUENCE_INFO_WITHOUT_RLE
 * 				   Both pairs are mutually exclusive.
 */
PB_SequenceInfo* get_sequence_info_cstring(uint8* input,
										   int mode);

/**
 * get_sequence_info_text()
 * 		Obtain sequence length, symbol frequencies and alphabet.
 *
 * 		text* input : text input string
 * 		int mode : determines case-sensitivity and whether
 * 				   run-length stats should be collected.
 * 				   Possible modes:
 * 					PB_SEQUENCE_INFO_CASE_SENSITIVE
 * 					PB_SEQUENCE_INFO_CASE_INSENSITIVE
 * 					PB_SEQUENCE_INFO_WITH_RLE
 * 					PB_SEQUENCE_INFO_WITHOUT_RLE
 * 				   Both pairs are mutually exclusive.
 */
PB_SequenceInfo* get_sequence_info_text(text* input,
										int mode);


#endif /* SEQUENCE_STATS_H_ */
