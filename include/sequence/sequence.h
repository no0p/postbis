/*-------------------------------------------------------------------------
*
* Copyright (c) 2013, Max Planck Institute for Marine Microbiology
*
* This software is released under the PostgreSQL License
*
* Author: Michael Schneider <mschneid@mpi-bremen.de>
*
* IDENTIFICATION
*   include/sequence/sequence.h
*
*-------------------------------------------------------------------------
*/
#ifndef SEQUENCE_SEQUENCE_H_
#define SEQUENCE_SEQUENCE_H_

#include "postgres.h"

/**
 * Size of input string alphabet
 */
#define PB_SOURCE_ALPHABET_SIZE		(1 << (sizeof(uint8) * 8))

/**
 * Size of ASCII alphabet
 */
#define PB_ASCII_SIZE		(128)

/**
 * Limit for input sequence length
 */
#define PB_MAX_INPUT_SEQUENCE_LENGTH		4294967295

/**
 * Limit for compressed sequence size
 */
#define PB_MAX_COMPRESSED_SEQUENCE_SIZE		1073741823

/**
 * Number of characters between substring-index entries
 * Smaller number -> denser index -> faster substring access
 * 									+ larger sequence in total
 * Larger number -> sparser index -> slower substring access
 * 										+ smaller sequence in total
 */
#define PB_INDEX_PART_SIZE 65536

/**
 * Type to buffer compressed data
 */
typedef uint64 PB_CompressionBuffer;

/**
 * Buffer size in bits
 */
#define PB_COMPRESSION_BUFFER_BIT_SIZE		(sizeof(PB_CompressionBuffer) * 8)

/**
 * Buffer size in bytes
 */
#define PB_COMPRESSION_BUFFER_BYTE_SIZE		(sizeof(PB_CompressionBuffer))

/**
 * Determines the data alignment of a data type.
 */
#define PB_ALIGNOF(type) offsetof (struct { char c; type member; }, member)

/**
 * Aligns bit or byte size to match PB_CompressionBuffer.
 */
#define PB_ALIGN_SIZE(x,y)		((x / y + (x % y > 0 ? 1 : 0)) * y)
#define PB_ALIGN_BIT_SIZE(x)	PB_ALIGN_SIZE(x, PB_COMPRESSION_BUFFER_BIT_SIZE)
#define PB_ALIGN_BYTE_SIZE(x)	PB_ALIGN_SIZE(x, sizeof(PB_CompressionBuffer))

/**
 * Determines the maximum of two numbers.
 */
#define PB_MAX(a,b) (a > b) ? a : b;

/*
 *  about 3% faster than toupper()
 */
#define TO_UPPER(c)	((c >= 97 && c <= 122) ? c - 32 : c)
#define TO_LOWER(c)	((c >= 65 && c <= 90) ? c + 32 : c)

/**
 * Type for prefix codes
 *
 * uint8 is strongly recommended as uint16 will result in huge decoding
 * maps, which can slow down compression significantly.
 */
typedef uint8 PB_PrefixCode;

/**
 * Size of a PB_PrefixCode
 */
#define PB_PREFIX_CODE_BIT_SIZE		(sizeof(PB_PrefixCode) * 8)
#define PB_PREFIX_CODE_BYTE_SIZE	(sizeof(PB_PrefixCode))

/*
 * Type to store number of consecutive equal characters.
 */
typedef uint8 PB_RunLength;

/**
 * Size of a PB_RunLength
 */
#define PB_RUN_LENGTH_BIT_SIZE		(sizeof(PB_RunLength) * 8)
#define PB_RUN_LENGTH_BYE_SIZE		(sizeof(PB_RunLength))

/**
 * Minimum number of consecutive equal characters to trigger
 * run-length encoding
 */
#define PB_MIN_RUN_LENGTH			8

/**
 * Maximum number of consecutive equal characters that can be
 * expressed by a run-length
 */
#define PB_MAX_RUN_LENGTH			(PB_MIN_RUN_LENGTH + (1 << PB_RUN_LENGTH_BIT_SIZE))

/**
 * Run-length symbol 26 0x1A SUB
 */
#define PB_RUN_LENGTH_SYMBOL		'\x1a'

/**
 * Type to store a swap run-length
 */
typedef uint16 PB_SwapRunLength;

/**
 * Size of PB_SwapRunLength
 */
#define PB_SWAP_RUN_LENGTH_BIT_SIZE		(sizeof(PB_SwapRunLength) * 8)
#define PB_SWAP_RUN_LENGTH_BYTE_SIZE	(sizeof(PB_SwapRunLength))

/**
 * Maximal number of characters that can be expressed with a swap run-length.
 */
#define PB_MAX_SWAP_RUN_LENGTH		((1 << PB_SWAP_RUN_LENGTH_BIT_SIZE) - 1)

/**
 * Minimal length of a sequence for the optimal code function
 * to actually try swapping.
 */
#define PB_MIN_LENGTH_FOR_SWAPPING 	(32768)

/**
 * Generic type for variable length data
 */
typedef struct varlena Varlena;

/**
 * PB_RleInfo stores the frequencies of run-length elements.
 */
typedef struct {
	uint32 rle_frequencies[PB_SOURCE_ALPHABET_SIZE];
	uint8 n_symbols;
	uint8* symbols;
} PB_RleInfo;

/**
 * PB_SequenceInfo stores information about a sequence such as:
 * 		-its length
 * 		-symbol frequencies
 * 		-bitmaps of occuring ASCII symbols
 * 		-list of symbols in order of frequency
 * 		-rle info
 */
typedef struct {
	uint32 sequence_length;
	uint32 frequencies[PB_SOURCE_ALPHABET_SIZE];
	PB_RleInfo* rle_info;
	uint64 ascii_bitmap_low;
	uint64 ascii_bitmap_high;
	uint8 n_symbols;
	bool ignore_case : 1;
	uint8* symbols;
} PB_SequenceInfo;

/**
 * pfrees all elements of a PB_SequenceInfo instance.
 */
#define PB_SEQUENCE_INFO_PFREE(sequence_info) \
	if (sequence_info != NULL) { \
		if (((PB_SequenceInfo*) sequence_info)->rle_info != NULL) { \
			if (((PB_RleInfo*) ((PB_SequenceInfo*) sequence_info)->rle_info)->symbols != NULL) \
				pfree(((PB_RleInfo*) ((PB_SequenceInfo*) sequence_info)->rle_info)->symbols); \
			\
			pfree(((PB_SequenceInfo*) sequence_info)->rle_info) ; \
		} \
		if (((PB_SequenceInfo*) sequence_info)->symbols != NULL) \
			pfree(((PB_SequenceInfo*) sequence_info)->symbols); \
		\
		pfree(sequence_info); \
	}

/**
 * A codeword for a symbol.
 */
typedef struct {
	PB_PrefixCode code;
	uint8 code_length;
	uint8 symbol;
} PB_Codeword;

/**
 * Holds a set of prefix codes.
 */
typedef struct {
	uint8 n_symbols;
	uint8 max_codeword_length;
	uint8 n_swapped_symbols;
	uint8 max_swapped_codeword_length;
	bool has_equal_length : 1;
	bool is_fixed : 1;
	bool uses_rle : 1;
	bool ignore_case : 1;
	uint8 fixed_id;
	uint64 swap_savings;
	uint64 ascii_bitmap_low;
	uint64 ascii_bitmap_high;
	PB_Codeword words[];
} PB_CodeSet;

/*
 * Stores a position in a compressed stream.
 */
typedef struct {
	uint32 block;
	uint8 bit;
	uint16 rle_shift;
	PB_SwapRunLength swap_shift;
} PB_IndexEntry;

/*
 * Stores a compressed sequence.
 *
 *	uint32 _vl_len			:	pgsql specific 4-byte length field; must only be set and get
 *								with pgsqls macros SET_VARSIZE() and VARSIZE()
 *	uint32 sequence_length	:	number of characters in the original sequence
 *	uint8 n_symbols			:	total number of symbols in code
 *	uint8 n_swapped_symbols	:	number of swapped symbols
 *	bool has_equal_length	:	true if all codewords have same length
 *	bool has_index			:	true if index is included
 *	bool is_fixed			:	true if fixed code was used
 *	bool uses_rle			:	true if rle was used
 *	uint8 _align2			:	data alignment
 *
 * The layout of the variable part in 'data' member of this struct is:
 * 	Variable member					|	size
 * ----------------------------------------------------------------------------
 * 	PB_Codeword symbols[];				|	a = sizeof(PB_Codeword) * (n_symbols - n_swapped_symbols)
 *	PB_Codeword swapped_symbols[];		|	b = sizeof(PB_Codeword) * (n_swapped_symbols)
 *	PB_IndexEntry index[];				|	c = has_index == true ? sizeof(PB_IndexEntry) * (sequence_length / PB_INDEX_PART_SIZE) : 0
 *	PB_CompressionBuffer stream[];		|	d = VARSIZE(_vl_len) - roundupto8(12 + a + b + c)
 *
 * Pointers to the variable members can be obtained by the following functions:
 * 	Variable member					|	function
 * ----------------------------------------------------------------------------
 * 	PB_Codeword symbols[];				|	PB_COMPRESSED_SEQUENCE_SYMBOL_POINTER(seq)
 *	PB_Codeword swapped_symbols[];		|	PB_COMPRESSED_SEQUENCE_SWAPPED_SYMBOL_POINTER(seq)
 *	PB_IndexEntry index[];				|	PB_COMPRESSED_SEQUENCE_INDEX_POINTER(seq)
 *	PB_CompressionBuffer stream[];		|	PB_COMPRESSED_SEQUENCE_STREAM_POINTER(seq)
 *
 */
typedef struct {
	uint32 _vl_len;
	uint32 sequence_length;
	uint8 n_symbols;
	uint8 n_swapped_symbols;
	bool has_equal_length : 1;
	bool has_index : 1;
	bool is_fixed : 1;
	bool uses_rle : 1;
	uint8 _align2;
	uint8 data[];
} PB_CompressedSequence;

/**
 * Number of elements in index table
 */
#define PB_COMPRESSED_SEQUENCE_INDEX_N_ELEMENTS(seq) \
	(((PB_CompressedSequence*)seq)->has_index ? \
	(((PB_CompressedSequence*)seq)->sequence_length / PB_INDEX_PART_SIZE) : 0)

/**
 * Returns the identifier of the employed fixed code. If
 * a sequence specific code is included it returns (-1).
 */
#define PB_COMPRESSED_SEQUENCE_FIXED_CODE_ID(seq) \
	(((PB_CompressedSequence*)seq)->is_fixed ? \
	((PB_CompressedSequence*)seq)->n_swapped_symbols : \
	-1)

/**
 * Returns a (PB_Codeword*) pointer to sequence specific
 * codewords. Returns NULL if a fixed code was used.
 */
#define PB_COMPRESSED_SEQUENCE_SYMBOL_POINTER(seq) \
	((((PB_CompressedSequence*)seq)->is_fixed) ? \
	NULL : \
	((PB_Codeword*)(((PB_CompressedSequence*)seq)->data)))

/**
 * Returns a (PB_Codeword*) pointer to sequence specific
 * swap codewords. Returns NULL if there are no swapped symbols.
 */
#define PB_COMPRESSED_SEQUENCE_SWAPPED_SYMBOL_POINTER(seq) \
	(((((PB_CompressedSequence*)seq)->is_fixed) || \
	((((PB_CompressedSequence*)seq)->n_swapped_symbols) == 0)) ? \
	NULL : \
	(((PB_Codeword*)(((PB_CompressedSequence*)seq)->data)) \
	+ (((PB_CompressedSequence*)seq)->n_symbols - ((PB_CompressedSequence*)seq)->n_swapped_symbols)))

/**
 * Returns a (PB_IndexEntry*) pointer to the first index
 * entry. Returns NULL if there is no index.
 */
#define PB_COMPRESSED_SEQUENCE_INDEX_POINTER(seq) \
	(((((PB_CompressedSequence*)seq)->has_index) == false) ? \
	NULL : \
	((PB_IndexEntry*)(((PB_Codeword*)(((PB_CompressedSequence*)seq)->data)) \
	+ ((PB_CompressedSequence*)seq)->n_symbols)))

/**
 * Returns the offset of the stream.
 */
#define PB_COMPRESSED_SEQUENCE_STREAM_OFFSET(seq) \
	(PB_ALIGN_BYTE_SIZE(( \
	sizeof(PB_CompressedSequence) + \
	((PB_CompressedSequence*)seq)->n_symbols * sizeof(PB_Codeword) + \
	PB_COMPRESSED_SEQUENCE_INDEX_N_ELEMENTS(seq) * sizeof(PB_IndexEntry))))

/**
 * Returns a (PB_CompressionBuffer*) to the beginning of
 * the compressed stream.
 */
#define PB_COMPRESSED_SEQUENCE_STREAM_POINTER(seq) \
	((PB_CompressionBuffer*) \
	(((uint8*) seq) + \
	PB_COMPRESSED_SEQUENCE_STREAM_OFFSET(seq)))

/**
 * Checks if a codeset can encode a given sequence.
 */
#define PB_CHECK_CODESET(codeset,info) \
	(!(info->n_symbols > codeset->n_symbols || \
	((info->ascii_bitmap_high ^ codeset->ascii_bitmap_high) & info->ascii_bitmap_high) || \
	((info->ascii_bitmap_low ^ codeset->ascii_bitmap_low) & info->ascii_bitmap_low)))

#endif /* SEQUENCE_SEQUENCE_H_ */
