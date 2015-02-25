/*-------------------------------------------------------------------------
*
* Copyright (c) 2013, Max Planck Institute for Marine Microbiology
*
* This software is released under the PostgreSQL License
*
* Author: Michael Schneider <mschneid@mpi-bremen.de>
*
* IDENTIFICATION
*   include/types/alphabet.h
*-------------------------------------------------------------------------
*/

#ifndef TYPES_ALPHABET_H_
#define TYPES_ALPHABET_H_

#include "postgres.h"
#include "fmgr.h"

#include "sequence/sequence.h"

/**
 * Dummy type. To facilitate future changes.
 */
typedef uint8	PB_AlphabetHeader;

/**
 * Dummy type. To facilitate future changes.
 */
typedef uint8	PB_Symbol;

/**
 * Dummy type. To facilitate future changes.
 */
typedef float4	PB_SymbolProbability;

/**
 * Holds an alphabet, possibly including probabilities.
 *
 * To create an alphabet:
 * 		int size = 10;
 * 		PB_Alphabet* alpha;
 * 		PB_CREATE_ALPHABET_WITH_PROBABILITIES(alpha, size);
 * 	or
 * 		PB_CREATE_ALPHABET_WITHOUT_PROBABILITIES(alpha, size);
 *
 *  To get the size:
 *  	size = PB_ALPHABET_SIZE(alpha);
 *
 *  To check whether probabilities are included:
 *  	if (PB_ALPHABET_TYPE(alpha) == PB_ALPHABET_WITH_PROBABILITIES)
 *  	...
 *
 *  To access the symbols:
 *  	PB_Symbol* symbols = PB_ALPHABET_SYMBOL_POINTER(alpha)
 *
 *  To access the frequencies:
 *  	PB_SymbolProbability* probabilities =
 *  		PB_ALPHABET_SYMBOL_PROBABILITY_POINTER(alpha)
 *
 *  To delete it:
 *  	pfree(alpha);
 */
typedef struct
{
	int32				_vl_len;
	PB_AlphabetHeader	alphabet_header;
	uint8				data[];
} PB_Alphabet;

/**
 * Returns the number of elements in the alphabet.
 */
#define PB_ALPHABET_SIZE(a)	((a)->alphabet_header & 0x7F)

/**
 * Returns the type of alphabet.
 * 	PB_ALPHABET_TYPE(alpha) == PB_ALPHABET_WITH_PROBABILITIES
 * 	PB_ALPHABET_TYPE(alpha) == PB_ALPHABET_WITHOUT_PROBABILITIES
 */
#define PB_ALPHABET_TYPE(a)	(((a)->alphabet_header & 0x80) >> 7)

/**
 * Changes the number of elements. Should not be used after initialization.
 */
#define PB_ALPHABET_SET_SIZE(a,s)	\
	((a)->alphabet_header = ((a)->alphabet_header & 0x80) | (s & 0x7F))

/**
 * Changes the type of alphabet. Should not be used after initialization.
 */
#define PB_ALPHABET_SET_TYPE(a,t)	\
	((a)->alphabet_header = ((a)->alphabet_header & 0x7F) | ((t << 7) & 0x80))

#define PB_ALPHABET_WITHOUT_PROBABILITIES		0
#define PB_ALPHABET_WITH_PROBABILITIES			1

/**
 * Calculates how much memory is required for an alphabet
 * with n elements and no probabilities.
 */
#define PB_ALPHABET_CALC_MEMSIZE_WITHOUT_PROBABILITIES(n)	\
	(sizeof(int32) + \
	sizeof(PB_AlphabetHeader) + \
	sizeof(PB_Symbol) * (n))

/**
 * Calculates how much memory is required for an alphabet
 * with n elements including probabilities.
 */
#define PB_ALPHABET_CALC_MEMSIZE_WITH_PROBABILITIES(n)	\
	(PB_ALPHABET_CALC_MEMSIZE_WITHOUT_PROBABILITIES((n)) + \
	sizeof(PB_SymbolProbability) * (n))

/**
 * Returns a pointer to the symbols in the alphabet.
 */
#define PB_ALPHABET_SYMBOL_POINTER(a) \
	((PB_Symbol*)((a)->data))

/**
 * Returns a pointer to the symbol probabilities in
 * the alphabet.
 */
#define PB_ALPHABET_SYMBOL_PROBABILITY_POINTER(a) \
	( (PB_SymbolProbability*) \
		( \
			PB_ALPHABET_TYPE((a)) == PB_ALPHABET_WITH_PROBABILITIES \
			? ( \
				((uint8*)((a))->data) + \
					PB_ALPHABET_SIZE((a)) * sizeof(PB_Symbol) \
			) : 0 \
		) \
	)

/**
 * Allocates memory for an alphabet of given size and
 * initializes the header.
 */
#define PB_ALPHABET_CREATE_WITH_PROBABILITIES(alpha, size) \
	alpha = palloc0(PB_ALPHABET_CALC_MEMSIZE_WITH_PROBABILITIES(size)); \
	SET_VARSIZE(alpha, PB_ALPHABET_CALC_MEMSIZE_WITH_PROBABILITIES(size)); \
	PB_ALPHABET_SET_SIZE(alpha, size); \
	PB_ALPHABET_SET_TYPE(alpha, PB_ALPHABET_WITH_PROBABILITIES);

/**
 * Allocates memory for an alphabet of given size and
 * initializes the header.
 */
#define PB_ALPHABET_CREATE_WITHOUT_PROBABILITIES(alpha, size) \
	alpha = palloc0(PB_ALPHABET_CALC_MEMSIZE_WITHOUT_PROBABILITIES(size)); \
	SET_VARSIZE(alpha, PB_ALPHABET_CALC_MEMSIZE_WITHOUT_PROBABILITIES(size)); \
	PB_ALPHABET_SET_SIZE(alpha, size); \
	PB_ALPHABET_SET_TYPE(alpha, PB_ALPHABET_WITHOUT_PROBABILITIES);

/**
 * parse_alphabet_from_text()
 * 		Internal conversion function.
 *
 * 	text* input : input text
 */
PB_Alphabet* parse_alphabet_from_text(text* input);

/**
 * alphabet_to_text()
 * 		Internal conversion function.
 *
 * 	PB_Alphabet* alphabet : input alphabet
 */
text* alphabet_to_text(PB_Alphabet* alphabet);

/**
 * alphabet_in()
 * 		Creates an alphabet from a cstring.
 *
 * 	uint8* input : null-terminated input (cstring)
 */
Datum alphabet_in (PG_FUNCTION_ARGS);

/**
 * alphabet_in_text()
 * 		Creates an alphabet from a text.
 *
 * 	text* input : text input
 */
Datum alphabet_in_text (PG_FUNCTION_ARGS);

/**
 * alphabet_out()
 * 		Alphabet to cstring.
 *
 * 	PB_Alphabet* input : Alphabet
 */
Datum alphabet_out (PG_FUNCTION_ARGS);

/**
 * alphabet_out_text()
 * 		Alphabet to text.
 *
 * 	PB_Alphabet* input : Alphabet
 */
Datum alphabet_out_text (PG_FUNCTION_ARGS);

/**
 * Creates an alphabet from a text sequence.
 *
 * text* input : the input sequence
 */
Datum get_alphabet_text_sequence (PG_FUNCTION_ARGS);

/**
 * Creates an alphabet from compressed sequence.
 *
 * PB_CompressedSequence* input : the input sequence
 */
PB_Alphabet* get_alphabet_compressed_sequence (PB_CompressedSequence* input, PB_CodeSet** fixed_codesets);

#endif /* TYPES_ALPHABET_H_ */
