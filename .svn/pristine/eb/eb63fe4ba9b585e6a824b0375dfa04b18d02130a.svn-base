/*-------------------------------------------------------------------------
*
* Copyright (c) 2013, Max Planck Institute for Marine Microbiology
*
* This software is released under the PostgreSQL License
*
* Author: Michael Schneider <mschneid@mpi-bremen.de>
*
* IDENTIFICATION
*   include/types/aligned_dna_sequence.h
*
*-------------------------------------------------------------------------
*/
#ifndef TYPES_ALIGNED_DNA_SEQUENCE_H_
#define TYPES_ALIGNED_DNA_SEQUENCE_H_

#include "postgres.h"

/*
 * Section 1 - the type modifier type
 *
 * Type modifiers must be represented by a single
 * integer value. This is the interface to handle
 * the aligned_rna_sequence-type modifiers as a single
 * integer value.
 */

typedef struct {
	uint32 case_sensitive : 1;
	uint32 restricting_alphabet : 2;
} PB_AlignedDnaSequenceTypMod;

#define PB_ALIGNED_DNA_TYPMOD_CASE_INSENSITIVE	0
#define PB_ALIGNED_DNA_TYPMOD_CASE_SENSITIVE	1

#define PB_ALIGNED_DNA_TYPMOD_IUPAC			0
#define PB_ALIGNED_DNA_TYPMOD_FLC			1
#define PB_ALIGNED_DNA_TYPMOD_ASCII			2

/*
 * Section 2 - public functions
 */

/**
 * aligned_dna_sequence_typmod_to_int()
 * 		Convert from PB_AlignedDnaSequenceTypMod to int
 *
 * 	PB_AlignedDnaSequenceTypMod typmod : type modifier
 */
int aligned_dna_sequence_typmod_to_int(PB_AlignedDnaSequenceTypMod typmod);

/**
 * int_to_aligned_dna_sequence_typmod()
 * 		Convert from int to PB_AlignedDnaSequenceTypMod
 *
 * 	int typmod : type modifier
 */
PB_AlignedDnaSequenceTypMod int_to_aligned_dna_sequence_typmod(int typmod);

/**
 * get_fixed_aligned_dna_code()
 * 		Returns a fixed code for the specified id.
 *
 * 	Id	|	Code Description
 * -----------------------------------------------------------
 * 	0	|	DNA four-letter code
 * 	1	|	DNA four-letter code, case sensitive
 * 	2	|	DNA IUPAC code
 * 	3	|	DNA IUPAC code, case sensitive
 * 	4	|	DNA four-letter code complement
 * 	5	|	DNA four-letter code complement, case sensitive
 * 	6	|	DNA IUPAC code complement
 * 	7	|	DNA IUPAC code complement, case sensitive
 *
 */
PB_CodeSet* get_fixed_aligned_dna_code(unsigned int fixed_code_id);

/**
 * get_fixed_aligned_dna_codes()
 * 		Returns pointer to fixed aligned DNA codes.
 */
PB_CodeSet** get_fixed_aligned_dna_codes(void);

/**
 * compress_aligned_dna_sequence()
 * 		Compress an aligned DNA sequence.
 *
 *	First it will be checked whether the given sequence matches
 *	the restricting alphabet specified with the type modifiers.
 *	Then the sequence will be optimally compressed.
 *
 * 	uint8* input : unterminated or null-terminated input sequence
 * 	PB_DnaSequenceTypMod typmod : target type modifier
 * 	PB_SequenceInfo info : given by input function
 */
PB_CompressedSequence* compress_aligned_dna_sequence(uint8* input,
													 PB_AlignedDnaSequenceTypMod typmod,
													 PB_SequenceInfo* info);

/**
 * decompress_aligned_dna_sequence()
 * 		Decompress an aligned DNA sequence
 *
 * 	PB_CompressedSequence* input : compressed input sequence
 * 	uint8* output : already allocated appropriately sized target memory area
 * 	uint32 from_position : position to start decompressing, first = 0
 * 	uint32 length : length of sequence to decompress
 */
void decompress_aligned_dna_sequence(PB_CompressedSequence* input,
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
 *	aligned_dna_sequence_typmod_in()
 *		Condense type modifier keywords into single integer value.
 *
 *	cstring[] input : lower-case keywords separated into array
 */
Datum aligned_dna_sequence_typmod_in(PG_FUNCTION_ARGS);

/**
 * 	aligned_dna_sequence_typmod_out()
 * 		Restore type modifier keywords from single integer value.
 *
 * 	int typmod : single value representing the type modifiers
 */
Datum aligned_dna_sequence_typmod_out(PG_FUNCTION_ARGS);

/**
 * aligned_dna_sequence_in()
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
Datum aligned_dna_sequence_in (PG_FUNCTION_ARGS);

/**
 * aligned_dna_sequence_in_varlena()
 * 		Compress a given input sequence.
 *
 * 	This function expects a varlena input sequence, that is
 * 	text, varchar or char. It will be called by the respective
 * 	cast functions.
 *
 * 	Varlena* input : input sequence
 * 	int typmod : single value representing target type modifier
 */
Datum aligned_dna_sequence_in_varlena (PG_FUNCTION_ARGS);

/**
 * aligned_dna_sequence_cast()
 * 		Decompress a given sequence and compress it again
 * 		using a different compression
 *
 * 	Varlena* input : input sequence
 * 	int typmod : single value representing target type modifier
 */
Datum aligned_dna_sequence_cast (PG_FUNCTION_ARGS);

/**
 * aligned_dna_sequence_out()
 * 		Decompress a sequence.
 *
 * 	PB_CompressedSequence* input : compressed input sequence
 */
Datum aligned_dna_sequence_out (PG_FUNCTION_ARGS);

/**
 * aligned_dna_sequence_out_varlena()
 * 		Decompress a sequence into varlena.
 *
 * 	PB_CompressedSequence* input : compressed input sequence
 */
Datum aligned_dna_sequence_out_varlena (PG_FUNCTION_ARGS);

/**
 * aligned_dna_sequence_substring()
 * 		Decompress a substring of a sequence.
 *
 * 	This function mimics the behaviour of the originals substr function.
 * 	The first position is 1.
 *
 * 	Varlena* input : compressed input sequence
 * 	int start : position to start from
 * 	int len : length of substring
 */
Datum aligned_dna_sequence_substring (PG_FUNCTION_ARGS);

/**
 * aligned_dna_sequence_compression_ratio()
 * 		Get compression ratio.
 *
 * 	The ratio between the size of the sequence as pgsql text
 * 	type and the size of the compressed sequence including all
 * 	required meta-data, such as the substring-index.
 *
 * 	PB_CompressedSequence* input : compressed input sequence
 */
Datum aligned_dna_sequence_compression_ratio (PG_FUNCTION_ARGS);

/**
 * aligned_dna_sequence_char_length()
 * 		Get length of sequence.
 *
 * 	PB_CompressedSequence* input : compressed input sequence
 */
Datum aligned_dna_sequence_char_length (PG_FUNCTION_ARGS);

/**
 * aligned_dna_sequence_complement()
 * 		Returns the complement of an aligned DNA sequence.
 *
 * 	PB_CompressedSequence* input : compressed input sequence
 */
Datum aligned_dna_sequence_complement(PG_FUNCTION_ARGS);

/**
 * aligned_dna_sequence_reverse()
 * 		Returns the reverse of an aligned DNA sequence.
 *
 * 	PB_CompressedSequence* input : compressed input sequence
 */
Datum aligned_dna_sequence_reverse(PG_FUNCTION_ARGS);

/**
 * aligned_dna_sequence_reverse_complement()
 * 		Returns the reverse of an aligned DNA sequence.
 *
 * 	PB_CompressedSequence* input : compressed input sequence
 */
Datum aligned_dna_sequence_reverse_complement(PG_FUNCTION_ARGS);

/**
 * get_alphabet_aligned_dna_sequence()
 * 		Calculates alphabet from an aligned DNA sequence.
 *
 * 	PB_CompressedSequence* input : compressed input sequence
 */
Datum get_alphabet_aligned_dna_sequence(PG_FUNCTION_ARGS);
/**
 * equal_aligned_dna()
 * 		Compares two aligned DNA sequence for equality
 *
 * 	Varlena* seq1 : first sequence
 * 	Varlena* seq2 : second sequence
 */
Datum equal_aligned_dna(PG_FUNCTION_ARGS);

/**
 * compare_aligned_dna_lt()
 * 		Compares two aligned DNA sequence for equality Less-than.
 *
 * 	Varlena* seq1 : first sequence
 * 	Varlena* seq2 : second sequence
 */
Datum compare_aligned_dna_lt(PG_FUNCTION_ARGS);

/**
 * compare_aligned_dna_le()
 * 		Compares two aligned DNA sequence for equality Less or equal.
 *
 * 	Varlena* seq1 : first sequence
 * 	Varlena* seq2 : second sequence
 */
Datum compare_aligned_dna_le(PG_FUNCTION_ARGS);

/**
 * compare_aligned_dna_gt()
 * 		Compares two aligned DNA sequence for equality. greater than
 *
 * 	Varlena* seq1 : first sequence
 * 	Varlena* seq2 : second sequence
 */
Datum compare_aligned_dna_gt(PG_FUNCTION_ARGS);

/**
 * compare_aligned_dna_ge()
 * 		Compares two aligned DNA sequence for equality. greater or equal
 *
 * 	Varlena* seq1 : first sequence
 * 	Varlena* seq2 : second sequence
 */
Datum compare_aligned_dna_ge(PG_FUNCTION_ARGS);

/**
 * compare_aligned_dna()
 * 		Compares two aligned DNA sequence for equality.
 *
 * 	Varlena* seq1 : first possibly toasted sequence
 * 	Varlena* seq2 : second possibly toasted sequence
 */
Datum compare_aligned_dna(PG_FUNCTION_ARGS);

/**
 * hash_aligned_dna()
 * 		Returns a CRC32 for a aligned DNA sequence.
 *
 * 	PB_CompressedSequence* seq1 : input sequence
 */
Datum hash_aligned_dna(PG_FUNCTION_ARGS);

/**
 * strpos_aligned_dna()
 * 		Finds the first occurrence of a pattern in a sequence.
 */
Datum strpos_aligned_dna(PG_FUNCTION_ARGS);

/**
 * octet_length_aligned_dna()
 * 		Returns byte size of datum.
 */
Datum octet_length_aligned_dna(PG_FUNCTION_ARGS);

#endif /* TYPES_ALIGNED_DNA_SEQUENCE_H_ */
