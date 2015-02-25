/*-------------------------------------------------------------------------
*
* Copyright (c) 2013, Max Planck Institute for Marine Microbiology
*
* This software is released under the PostgreSQL License
*
* Author: Michael Schneider <mschneid@mpi-bremen.de>
*
* IDENTIFICATION
*   include/types/aa_sequence.h
*-------------------------------------------------------------------------
*/

#ifndef TYPES_AA_SEQUENCE_H_
#define TYPES_AA_SEQUENCE_H_

#include "postgres.h"
#include "fmgr.h"

#include "sequence/sequence.h"

/*
 * Section 1 - the type modifier type
 *
 * Type modifiers must be represented by a single
 * integer value. This is the interface to handle
 * the aa_sequence-type modifiers as a single
 * integer value.
 */

typedef struct {
	uint32 case_sensitive : 1;
	uint32 restricting_alphabet : 2;
} PB_AaSequenceTypMod;

#define PB_AA_TYPMOD_CASE_INSENSITIVE 0
#define PB_AA_TYPMOD_CASE_SENSITIVE	1

#define PB_AA_TYPMOD_IUPAC			0
#define PB_AA_TYPMOD_ASCII			1

/*
 * Section 2 - public functions
 */

/**
 * aa_sequence_typmod_to_int()
 * 		Convert from PB_AaSequenceTypMod to int
 *
 * 	PB_AaSequenceTypMod tm : typmod
 */
int aa_sequence_typmod_to_int(PB_AaSequenceTypMod tm);

/**
 * int_to_aa_sequence_typmod()
 * 		Convert from int to PB_AaSequenceTypMod
 *
 * 	int tm : typmod
 */
PB_AaSequenceTypMod int_to_aa_sequence_typmod(int tm);

/*
 * Section 2 - Interface functions for pgsql
 *
 * The Version 1 Calling Conventions require each function
 * that will be called from pgsql to be declared like:
 *  * Datum function(PG_FUNCTION_ARGS)
 * Parameters are hidden from the function declaration, yet
 * described in the comment.
 */

/**
 * get_fixed_aa_code()
 * 		Returns a fixed code for the specified id.
 *
 * 	Id	|	Code Description
 * -----------------------------------------------------------
 * 	0	|	AA IUPAC code
 * 	1	|	AA IUPAC code, case sensitive
 *
 */
PB_CodeSet* get_fixed_aa_code(unsigned int fixed_code_id);

/**
 * get_fixed_aa_codes()
 * 		Returns pointer to fixed AA codes.
 */
PB_CodeSet** get_fixed_aa_codes(void);

/**
 * compress_aa_sequence()
 * 		Compress an AA sequence.
 *
 *	First it will be checked whether the given sequence matches
 *	the restricting alphabet specified with the type modifiers.
 *	Then one of three types of compression will be performed,
 *	depending on which one is most efficient.
 *
 *	A) Simple Huffman-Codes
 *	B) Equal lengths codes
 *	C) Fixed pre-built codes
 *
 * 	uint8* input : unterminated or null-terminated input sequence
 * 	PB_DnaSequenceTypMod typmod : target type modifier
 * 	PB_SequenceInfo info : given by input function
 */
PB_CompressedSequence* compress_aa_sequence(uint8* input,
											PB_AaSequenceTypMod typmod,
											PB_SequenceInfo* info);

/**
 * decompress_aa_sequence()
 * 		Decompress a AA sequence
 *
 * 	PB_CompressedSequence* input : compressed input sequence
 * 	uint8* output : already allocated appropriately sized target memory area
 * 	uint32 from_position : position to start decompressing, first = 0
 * 	uint32 length : length of sequence to decompress
 */
void decompress_aa_sequence(PB_CompressedSequence* input,
							uint8* output,
							uint32 from_position,
							uint32 length);

/*
 * Section 3 - Interface functions for pgsql
 *
 * The Version 1 Calling Conventions require each function
 * that will be called from pgsql to be declared like:
 *  * Datum function(PG_FUNCTION_ARGS)
 * Parameters are hidden from the function declaration, yet
 * described in the comment.
 *
 */

/**
 *	aa_sequence_typmod_in()
 *		Condense type modifier keywords into single integer value.
 *
 *	cstring[] input : lower-case keywords separated into array
 */
Datum aa_sequence_typmod_in(PG_FUNCTION_ARGS);

/**
 * 	aa_sequence_typmod_out()
 * 		Restore type modifier keywords from single integer value.
 *
 * 	int typmod : single value representing the type modifiers
 */
Datum aa_sequence_typmod_out(PG_FUNCTION_ARGS);

/**
 * aa_sequence_in()
 * 		Compress a given input sequence.
 *
 * 	Due to "bizarrely inconsistent rules" (Tom Lane) in the
 * 	SQL standard, pgsql will always set the typmod parameter
 * 	to (-1). See coerce_type comment in parser/parse_coerce.h.
 * 	If type modifiers were specified the cast function will be
 * 	called afterwards. That means that the whole compression
 * 	process will inevitably be performed twice.
 *
 * 	uint8* input : null-terminated input sequence (cstring)
 * 	Oid oid : oid of the sequence type
 * 	int typmod : single value representing target type modifier
 */
Datum aa_sequence_in (PG_FUNCTION_ARGS);

/**
 * aa_sequence_in_varlena()
 * 		Compress a given input sequence.
 *
 * 	This function expects a varlena input sequence, that is
 * 	text, varchar or char. It will be called by the respective
 * 	cast functions.
 *
 * 	Varlena* input : input sequence
 * 	int typmod : single value representing target type modifier
 */
Datum aa_sequence_in_varlena (PG_FUNCTION_ARGS);

/**
 * aa_sequence_cast()
 * 		Decompress a given sequence and compress it again
 * 		using a different compression
 *
 * 	Varlena* input : input sequence
 * 	int typmod : single value representing target type modifier
 */
Datum aa_sequence_cast (PG_FUNCTION_ARGS);

/**
 * aa_sequence_out()
 * 		Decompress a sequence.
 *
 * 	PB_CompressedSequence* input : compressed input sequence
 */
Datum aa_sequence_out (PG_FUNCTION_ARGS);

/**
 * aa_sequence_out_varlena()
 * 		Decompress a sequence into varlena.
 *
 * 	PB_CompressedSequence* input : compressed input sequence
 */
Datum aa_sequence_out_varlena (PG_FUNCTION_ARGS);

/**
 * aa_sequence_substring()
 * 		Decompress a substring of a sequence.
 *
 * 	This function mimics the originals substr function's
 * 	behaviour. The first position is 1.
 *
 * 	Varlena* input : compressed input sequence
 * 	int start : position to start from
 * 	int len : length of substring
 */
Datum aa_sequence_substring (PG_FUNCTION_ARGS);

/**
 * aa_sequence_compression_ratio()
 * 		Get compression ratio.
 *
 * 	The ratio between the size of the sequence as pgsql text
 * 	type and the size of the compressed sequence including all
 * 	required meta-data, such as the substring-index.
 *
 * 	PB_CompressedSequence* input : compressed input sequence
 */
Datum aa_sequence_compression_ratio (PG_FUNCTION_ARGS);

/**
 * aa_sequence_char_length()
 * 		Get length of sequence.
 *
 * 	CompressedSequence* input : compressed input sequence
 */
Datum aa_sequence_char_length (PG_FUNCTION_ARGS);

/**
 * aa_sequence_reverse()
 * 		Returns the reverse of a AA sequence.
 *
 * 	PB_CompressedSequence* input : compressed input sequence
 */
Datum aa_sequence_reverse(PG_FUNCTION_ARGS);

/**
 * get_alphabet_aa_sequence()
 * 		Calculates alphabet from AA sequence.
 *
 * 	PB_CompressedSequence* input : compressed input sequence
 */
Datum get_alphabet_aa_sequence(PG_FUNCTION_ARGS);

/**
 * equal_aa()
 * 		Compares two AA sequence for equality
 *
 * 	Varlena* seq1 : first sequence
 * 	Varlena* seq2 : second sequence
 */
Datum equal_aa(PG_FUNCTION_ARGS);

/**
 * compare_aa_lt()
 * 		Compares two AA sequence for equality Less-than.
 *
 * 	Varlena* seq1 : first sequence
 * 	Varlena* seq2 : second sequence
 */
Datum compare_aa_lt(PG_FUNCTION_ARGS);

/**
 * compare_aa_le()
 * 		Compares two AA sequence for equality Less or equal.
 *
 * 	Varlena* seq1 : first sequence
 * 	Varlena* seq2 : second sequence
 */
Datum compare_aa_le(PG_FUNCTION_ARGS);

/**
 * compare_aa_gt()
 * 		Compares two AA sequence for equality. greater than
 *
 * 	Varlena* seq1 : first sequence
 * 	Varlena* seq2 : second sequence
 */
Datum compare_aa_gt(PG_FUNCTION_ARGS);

/**
 * compare_aa_ge()
 * 		Compares two AA sequence for equality. greater or equal
 *
 * 	Varlena* seq1 : first sequence
 * 	Varlena* seq2 : second sequence
 */
Datum compare_aa_ge(PG_FUNCTION_ARGS);

/**
 * compare_aa()
 * 		Compares two AA sequence for equality.
 *
 * 	Varlena* seq1 : first possibly toasted sequence
 * 	Varlena* seq2 : second possibly toasted sequence
 */
Datum compare_aa(PG_FUNCTION_ARGS);

/**
 * hash_aa()
 * 		Returns a CRC32 for a AA sequence.
 *
 * 	PB_CompressedSequence* seq1 : input sequence
 */
Datum hash_aa(PG_FUNCTION_ARGS);

/**
 * strpos_aa()
 * 		Finds the first occurrence of a pattern in a sequence.
 */
Datum strpos_aa(PG_FUNCTION_ARGS);

/**
 * octet_length_aa()
 * 		Returns byte size of datum.
 */
Datum octet_length_aa(PG_FUNCTION_ARGS);

#endif /* TYPES_AA_SEQUENCE_H_ */
