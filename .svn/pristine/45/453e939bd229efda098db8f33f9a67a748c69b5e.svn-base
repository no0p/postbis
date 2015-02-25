/*-------------------------------------------------------------------------
*
* Copyright (c) 2013, Max Planck Institute for Marine Microbiology
*
* This software is released under the PostgreSQL License
*
* Author: Michael Schneider <mschneid@mpi-bremen.de>
*
* IDENTIFICATION
*   include/sequence/decompression_iteration.h
*
*-------------------------------------------------------------------------
*/

#ifndef SEQUENCE_DECOMPRESSION_ITERATION_H_
#define SEQUENCE_DECOMPRESSION_ITERATION_H_

#include "postgres.h"
#include "access/tuptoaster.h"
#include "c.h"

#include "sequence/sequence.h"
#include "sequence/compression.h"

#include "utils/debug.h"

typedef struct {
	uint8 symbol;
	uint8 code_length;
} PB_DecodingMap;

#define PB_DECODE_MAP_SIZE (1 << PB_PREFIX_CODE_BIT_SIZE)

#define DECODE(input_pointer, buffer, bits_in_buffer, val, length, map) { \
	val = buffer >> (PB_COMPRESSION_BUFFER_BIT_SIZE - PB_PREFIX_CODE_BIT_SIZE); \
	length = map[val].code_length; \
	if (length <= bits_in_buffer) { \
		bits_in_buffer -= length; \
		buffer = buffer << length; \
	} else { \
		PB_CompressionBuffer next = *input_pointer; \
		val = (buffer | (next >> bits_in_buffer)) >> (PB_COMPRESSION_BUFFER_BIT_SIZE - PB_PREFIX_CODE_BIT_SIZE); \
		length = map[val].code_length; \
		bits_in_buffer -= length; \
		buffer = next << (-bits_in_buffer); \
		bits_in_buffer += PB_COMPRESSION_BUFFER_BIT_SIZE; \
		input_pointer++; \
	} \
}

/**
 * This macro reads a given number of bits into a target variable.
 *
 * Parameters:
 * 	PB_CompressionBuffer* input_pointer : compressed sequence
 * 	PB_CompressionBuffer buffer : compression buffer
 * 	int bits_in_buffer : bits in compression buffer
 *  (??) target : where to write to
 * 	int length : length of next code
 */
#define READ_N_BITS(input_pointer, buffer, bits_in_buffer, target, no_bits) { \
	int len = (int) no_bits; \
	if (len > bits_in_buffer) { \
		PB_CompressionBuffer next = *input_pointer; \
		target = (buffer | (next >> bits_in_buffer)) >> (PB_COMPRESSION_BUFFER_BIT_SIZE - len); \
		bits_in_buffer -= len; \
		buffer = next << (-bits_in_buffer); \
		bits_in_buffer += PB_COMPRESSION_BUFFER_BIT_SIZE; \
		input_pointer++; \
	} else { \
		target = buffer >> (PB_COMPRESSION_BUFFER_BIT_SIZE - len); \
		bits_in_buffer -= len; \
		buffer = buffer << len; \
	} \
}

#define PB_NO_SWAP_MAP		0
#define PB_SWAP_MAP			1

/**
 * get_decoding_map()
 * 		Creates a decoding map for a prefix code set.
 */
static PB_DecodingMap* get_decoding_map(const PB_CodeSet* codeset, int mode)
{
	PB_DecodingMap* map;
	int i,j;
	int from = 0;
	int to = 0;

	if (mode == PB_NO_SWAP_MAP)
	{
		from = 0;
		to = codeset->n_symbols - codeset->n_swapped_symbols;
	}
	else if (mode == PB_SWAP_MAP)
	{
		from = codeset->n_symbols - codeset->n_swapped_symbols;
		to = codeset->n_symbols;
	}

	map = (PB_DecodingMap*) palloc(PB_DECODE_MAP_SIZE * sizeof(PB_DecodingMap));
	memset(map, 0xFF, PB_DECODE_MAP_SIZE * sizeof(PB_DecodingMap));

	for (i = from; i < to; i++)
	{
		const int lower_bound = codeset->words[i].code;
		const int upper_bound = lower_bound + (1 << (PB_PREFIX_CODE_BIT_SIZE - codeset->words[i].code_length));

		for (j = lower_bound; j < upper_bound; j++)
		{
			map[j].symbol = codeset->words[i].symbol;
			map[j].code_length = codeset->words[i].code_length;
		}

		PB_DEBUG3(errmsg("get_decoding_map(): i:%d lower_bound:%u upper_bound:%u symbol:%c code:%u length:%u/%ld",i,lower_bound,upper_bound,codeset->words[i].symbol,codeset->words[i].code,codeset->words[i].code_length, PB_PREFIX_CODE_BIT_SIZE));
	}

	return map;
}

/*
 * Iterate over a compressed sequence.
 *
 * I know this macro is a big ugly monster but it is 4 to 5-fold faster than
 * any other approach I tried.
 *
 * If you have suggestions to improve this, please contact me at
 * 		mschneid@mpi-bremen.de
 *
 * 	Usage:
 * 		uint8 c;
 *		PB_BEGIN_DECODE(input, from_position, length, fixed_dna_codes, c) {
 *			// do something with c
 *		} PB_END_DECODE
 */
#define PB_BEGIN_DECODE(__pb_decode_input, __pb_decode_start_position, __pb_decode_output_length, __pb_decode_fixed_codesets, __pb_decode_output) {\
	PB_CompressedSequence* __pb_decode_input_header;\
	PB_CodeSet* __pb_decode_codeset;\
	PB_DecodingMap* __pb_decode_map;\
	PB_DecodingMap* __pb_decode_swap_map;\
	PB_IndexEntry* __pb_decode_start_entry = NULL;\
	PB_CompressionBuffer __pb_decode_buffer;\
	int __pb_decode_bits_in_buffer;\
	int __pb_decode_i;\
	int __pb_decode_swap_counter;\
	int __pb_decode_stream_offset;\
	Varlena* __pb_decode_input_slice;\
	PB_CompressionBuffer* __pb_decode_input_pointer;\
	uint8 __pb_decode_current = 0;\
	int __pb_decode_n_rle_out = 0;\
	uint8 __pb_decode_master_symbol = 0;\
	int __pb_decode_max_codeword_length ;\
	const int __pb_decode_raw_size = toast_raw_datum_size((Datum)__pb_decode_input);\
\
	PB_TRACE(errmsg("BEGIN_DECODE(%u,%u)", __pb_decode_start_position, __pb_decode_output_length));\
\
	__pb_decode_input_header = (PB_CompressedSequence*) PG_DETOAST_DATUM_SLICE(__pb_decode_input,\
																			   0,\
																			   sizeof(PB_CompressedSequence) - VARHDRSZ);\
\
	PB_DEBUG1(errmsg("PB_BEGIN_DECODE(): input header detoasted\n\tsequence_length:%u\n\tn_symbols:%u\n\tn_swapped_symbols:%u\n\thas_equal_length:%d\n\thas_index:%d\n\tis_fixed:%d\n\tuses_rle:%d",\
			__pb_decode_input_header->sequence_length, __pb_decode_input_header->n_symbols, __pb_decode_input_header->n_swapped_symbols, __pb_decode_input_header->has_equal_length,\
			__pb_decode_input_header->has_index, __pb_decode_input_header->is_fixed, __pb_decode_input_header->uses_rle));\
\
	if (__pb_decode_input_header->is_fixed) {\
		__pb_decode_codeset = __pb_decode_fixed_codesets[__pb_decode_input_header->n_swapped_symbols];\
	} else {\
		int __pb_decode_code_size = sizeof(PB_Codeword) * __pb_decode_input_header->n_symbols;\
		PB_Codeword* __pb_decode_code;\
\
		__pb_decode_codeset = palloc0(sizeof(PB_CodeSet) + __pb_decode_code_size);\
		__pb_decode_codeset->n_symbols = __pb_decode_input_header->n_symbols;\
		__pb_decode_codeset->n_swapped_symbols = __pb_decode_input_header->n_swapped_symbols;\
		__pb_decode_codeset->is_fixed = FALSE;\
		__pb_decode_codeset->has_equal_length = __pb_decode_input_header->has_equal_length;\
		__pb_decode_codeset->uses_rle = __pb_decode_input_header->uses_rle;\
\
		pfree(__pb_decode_input_header);\
		__pb_decode_input_header = (PB_CompressedSequence*)\
									PG_DETOAST_DATUM_SLICE(__pb_decode_input,0,sizeof(PB_CompressedSequence) - VARHDRSZ +\
														   __pb_decode_code_size);\
\
		__pb_decode_code = PB_COMPRESSED_SEQUENCE_SYMBOL_POINTER(__pb_decode_input_header);\
		memcpy(__pb_decode_codeset->words, __pb_decode_code, __pb_decode_code_size);\
\
		PB_DEBUG1(errmsg("PB_BEGIN_DECODE():Sequence specific code copied"));\
	}\
\
	__pb_decode_map = get_decoding_map(__pb_decode_codeset, PB_NO_SWAP_MAP);\
	__pb_decode_swap_map = get_decoding_map(__pb_decode_codeset, PB_SWAP_MAP);\
\
	if (__pb_decode_codeset->n_swapped_symbols > 0) {\
		__pb_decode_master_symbol = __pb_decode_codeset->words[__pb_decode_codeset->n_symbols - __pb_decode_codeset->n_swapped_symbols].symbol;\
		__pb_decode_max_codeword_length = __pb_decode_codeset->words[__pb_decode_codeset->n_symbols - 1].code_length +\
										  __pb_decode_codeset->words[__pb_decode_master_symbol].code_length +\
										  PB_SWAP_RUN_LENGTH_BIT_SIZE;\
	} else {\
		__pb_decode_max_codeword_length = __pb_decode_codeset->words[__pb_decode_codeset->n_symbols - 1].code_length;\
	}\
\
	if (__pb_decode_input_header->has_index) {\
		int __pb_decode_start_entry_no;\
\
		__pb_decode_start_entry_no = (__pb_decode_start_position + 1) / PB_INDEX_PART_SIZE - 1;\
		if (__pb_decode_start_entry_no >= 0) {\
			Varlena* __pb_decode_data_slice = (Varlena*)\
				PG_DETOAST_DATUM_SLICE(__pb_decode_input,\
									   sizeof(PB_CompressedSequence) - VARHDRSZ +\
									   sizeof(PB_Codeword) * __pb_decode_input_header->n_symbols +\
									   sizeof(PB_IndexEntry) * __pb_decode_start_entry_no,\
									   sizeof(PB_IndexEntry));\
\
			__pb_decode_start_entry = palloc0(sizeof(PB_IndexEntry));\
			memcpy(__pb_decode_start_entry, VARDATA_ANY(__pb_decode_data_slice), sizeof(PB_IndexEntry));\
			pfree(__pb_decode_data_slice);\
\
			PB_DEBUG1(errmsg("PB_BEGIN_DECODE(): index found, uses entry no %d", __pb_decode_start_entry_no));\
		}\
	}\
\
	__pb_decode_stream_offset = PB_COMPRESSED_SEQUENCE_STREAM_OFFSET(__pb_decode_input_header) - VARHDRSZ;\
\
	PB_DEBUG1(errmsg("PB_BEGIN_DECODE():calculated stream offset:%u", __pb_decode_stream_offset));\
\
	if (__pb_decode_start_entry == NULL) {\
		int __pb_decode_slice_size = (__pb_decode_start_position + __pb_decode_output_length) *\
									  __pb_decode_max_codeword_length;\
\
		if (__pb_decode_codeset->uses_rle)\
			__pb_decode_slice_size += __pb_decode_max_codeword_length + 8;\
\
		__pb_decode_slice_size = __pb_decode_slice_size / PB_COMPRESSION_BUFFER_BIT_SIZE + 1;\
		__pb_decode_slice_size *= PB_COMPRESSION_BUFFER_BYTE_SIZE;\
\
		if (__pb_decode_slice_size + __pb_decode_stream_offset > __pb_decode_raw_size)\
			__pb_decode_slice_size = __pb_decode_raw_size - __pb_decode_stream_offset;\
\
		__pb_decode_input_slice = (Varlena*)\
				PG_DETOAST_DATUM_SLICE(__pb_decode_input,\
									   __pb_decode_stream_offset,\
									   __pb_decode_slice_size);\
\
		PB_DEBUG1(errmsg("PB_BEGIN_DECODE(): skipping through sequence\n\tslice size is %d bytes", __pb_decode_slice_size));\
\
		__pb_decode_i = __pb_decode_start_position - 1;\
		__pb_decode_input_pointer = (PB_CompressionBuffer*) VARDATA_ANY(__pb_decode_input_slice);\
\
		if (__pb_decode_codeset->n_swapped_symbols > 0) {\
			__pb_decode_buffer = *__pb_decode_input_pointer;\
			__pb_decode_input_pointer++;\
			__pb_decode_swap_counter = __pb_decode_buffer >> (PB_COMPRESSION_BUFFER_BIT_SIZE - PB_SWAP_RUN_LENGTH_BIT_SIZE);\
			__pb_decode_bits_in_buffer = PB_COMPRESSION_BUFFER_BIT_SIZE - PB_SWAP_RUN_LENGTH_BIT_SIZE;\
			__pb_decode_buffer = __pb_decode_buffer << PB_SWAP_RUN_LENGTH_BIT_SIZE;\
		} else {\
			__pb_decode_buffer = 0;\
			__pb_decode_bits_in_buffer = 0;\
			__pb_decode_swap_counter = __pb_decode_input_header->sequence_length + 1;\
		}\
\
		PB_DEBUG1(errmsg("buf:%08X%08X bib:%d", (unsigned int) (__pb_decode_buffer >> 32), (unsigned int) __pb_decode_buffer, __pb_decode_bits_in_buffer));\
	} else {\
		int __pb_decode_slice_start = __pb_decode_stream_offset +\
									  __pb_decode_start_entry->block * PB_COMPRESSION_BUFFER_BYTE_SIZE;\
		int __pb_decode_slice_size = (__pb_decode_output_length + (__pb_decode_start_position % PB_INDEX_PART_SIZE)) *\
									  __pb_decode_max_codeword_length;\
\
		if (__pb_decode_codeset->uses_rle)\
			__pb_decode_slice_size += __pb_decode_max_codeword_length + 8;\
\
		__pb_decode_slice_size = __pb_decode_slice_size / PB_COMPRESSION_BUFFER_BIT_SIZE + 1;\
		__pb_decode_slice_size *= PB_COMPRESSION_BUFFER_BYTE_SIZE;\
\
		if (__pb_decode_slice_size + __pb_decode_slice_start > __pb_decode_raw_size)\
			__pb_decode_slice_size = __pb_decode_raw_size - __pb_decode_slice_start;\
\
		__pb_decode_input_slice = (Varlena*)\
					PG_DETOAST_DATUM_SLICE(__pb_decode_input,\
										   __pb_decode_slice_start,\
										   __pb_decode_slice_size);\
\
		PB_DEBUG1(errmsg("PB_BEGIN_DECODE(): index entry given\n\tstarting in block %u\n\tslice starts at byte %d\n\tslice size is %d bytes\nrle_shift:%u\nswap_shift:%u",\
						 __pb_decode_start_entry->block, __pb_decode_slice_start, __pb_decode_slice_size, __pb_decode_start_entry->rle_shift, __pb_decode_start_entry->swap_shift));\
\
		__pb_decode_input_pointer = (PB_CompressionBuffer*) VARDATA_ANY(__pb_decode_input_slice);\
		__pb_decode_bits_in_buffer = PB_COMPRESSION_BUFFER_BIT_SIZE - __pb_decode_start_entry->bit;\
		__pb_decode_buffer = *(__pb_decode_input_pointer) << __pb_decode_start_entry->bit;\
		__pb_decode_input_pointer++;\
		__pb_decode_i = ((__pb_decode_start_position + 1) % PB_INDEX_PART_SIZE) - 1 + __pb_decode_start_entry->rle_shift;\
		if (__pb_decode_codeset->n_swapped_symbols > 0)\
			__pb_decode_swap_counter = __pb_decode_start_entry->swap_shift;\
		else\
			__pb_decode_swap_counter = __pb_decode_input_header->sequence_length + 1;\
	}\
\
	if (__pb_decode_input_header->is_fixed == FALSE)\
		pfree(__pb_decode_codeset);\
\
	if (__pb_decode_start_entry != NULL)\
		pfree(__pb_decode_start_entry);\
\
	PB_DEBUG1(errmsg("PB_BEGIN_DECODE(): reading %d chars to skip, swap_counter = %d, bib=%d, rle_shift=%u",\
					  __pb_decode_i + 1, __pb_decode_swap_counter, __pb_decode_bits_in_buffer,\
					  __pb_decode_start_entry == NULL ? -1 : __pb_decode_start_entry->rle_shift));\
\
	while (__pb_decode_i >= 0) {\
		PB_PrefixCode __pb_decode_val;\
		int __pb_decode_length;\
\
		DECODE(__pb_decode_input_pointer, __pb_decode_buffer, __pb_decode_bits_in_buffer, __pb_decode_val, __pb_decode_length, __pb_decode_map);\
		__pb_decode_current = __pb_decode_map[__pb_decode_val].symbol;\
		__pb_decode_i--;\
\
		if (__pb_decode_current == __pb_decode_master_symbol) {\
			__pb_decode_swap_counter--;\
			if (__pb_decode_swap_counter < 0) {\
				DECODE(__pb_decode_input_pointer, __pb_decode_buffer, __pb_decode_bits_in_buffer, __pb_decode_val, __pb_decode_length, __pb_decode_swap_map);\
				__pb_decode_current = __pb_decode_swap_map[__pb_decode_val].symbol;\
				READ_N_BITS(__pb_decode_input_pointer, __pb_decode_buffer, __pb_decode_bits_in_buffer, __pb_decode_swap_counter, PB_SWAP_RUN_LENGTH_BIT_SIZE);\
			}\
		}\
\
		if (__pb_decode_current == PB_RUN_LENGTH_SYMBOL) {\
			int __pb_decode_repeated_chars = 0;\
\
			READ_N_BITS(__pb_decode_input_pointer, __pb_decode_buffer, __pb_decode_bits_in_buffer, __pb_decode_repeated_chars, PB_RUN_LENGTH_BIT_SIZE);\
\
			__pb_decode_i -= __pb_decode_repeated_chars + PB_MIN_RUN_LENGTH - 2;\
		}\
	}\
\
	if (__pb_decode_i < -1) {\
		int __pb_decode_repeated_chars = -__pb_decode_i;\
		PB_PrefixCode __pb_decode_val;\
		int __pb_decode_length;\
\
		if (__pb_decode_repeated_chars >= __pb_decode_output_length)\
		{\
			__pb_decode_repeated_chars = __pb_decode_output_length;\
		}\
\
		PB_DEBUG1(errmsg("PB_BEGIN_DECODE(): %d remaining chars from last rle, i=%d", __pb_decode_repeated_chars, __pb_decode_i));\
\
		DECODE(__pb_decode_input_pointer, __pb_decode_buffer, __pb_decode_bits_in_buffer, __pb_decode_val, __pb_decode_length, __pb_decode_map);\
		__pb_decode_i--;\
		__pb_decode_current = __pb_decode_map[__pb_decode_val].symbol;\
		if (__pb_decode_current == __pb_decode_master_symbol)\
		{\
			__pb_decode_swap_counter--;\
			if (__pb_decode_swap_counter < 0)\
			{\
				PB_PrefixCode __pb_decode_val;\
				int __pb_decode_length;\
\
				DECODE(__pb_decode_input_pointer, __pb_decode_buffer, __pb_decode_bits_in_buffer, __pb_decode_val, __pb_decode_length, __pb_decode_swap_map);\
				__pb_decode_current = __pb_decode_swap_map[__pb_decode_val].symbol;\
\
				READ_N_BITS(__pb_decode_input_pointer, __pb_decode_buffer, __pb_decode_bits_in_buffer, __pb_decode_swap_counter, PB_SWAP_RUN_LENGTH_BIT_SIZE);\
			}\
		}\
\
		PB_DEBUG1(errmsg("PB_BEGIN_DECODE(): write %dx%c", __pb_decode_repeated_chars, __pb_decode_current));\
\
	__pb_decode_n_rle_out = __pb_decode_repeated_chars;\
	}\
\
	__pb_decode_i = __pb_decode_output_length - 1;\
\
	PB_DEBUG1(errmsg("PB_BEGIN_DECODE(): reading %d chars, swap_counter = %d, bib=%d, nrle=%d",\
					 __pb_decode_i + 1, __pb_decode_swap_counter, __pb_decode_bits_in_buffer, __pb_decode_n_rle_out));\
\
	while (__pb_decode_i >= 0) {\
		__pb_decode_i--;\
		__pb_decode_n_rle_out--;\
\
		if (__pb_decode_n_rle_out < 0) {\
			PB_PrefixCode __pb_decode_val;\
			int __pb_decode_length;\
\
			DECODE(__pb_decode_input_pointer, __pb_decode_buffer, __pb_decode_bits_in_buffer, __pb_decode_val, __pb_decode_length, __pb_decode_map);\
			__pb_decode_current = __pb_decode_map[__pb_decode_val].symbol;\
\
			if (__pb_decode_current == __pb_decode_master_symbol) {\
				__pb_decode_swap_counter--;\
				if (__pb_decode_swap_counter < 0) {\
					PB_PrefixCode __pb_decode_val;\
					int __pb_decode_length;\
\
					DECODE(__pb_decode_input_pointer, __pb_decode_buffer, __pb_decode_bits_in_buffer, __pb_decode_val, __pb_decode_length, __pb_decode_swap_map);\
					__pb_decode_current = __pb_decode_swap_map[__pb_decode_val].symbol;\
\
					READ_N_BITS(__pb_decode_input_pointer, __pb_decode_buffer, __pb_decode_bits_in_buffer, __pb_decode_swap_counter, PB_SWAP_RUN_LENGTH_BIT_SIZE);\
				}\
			}\
\
			if (__pb_decode_current == PB_RUN_LENGTH_SYMBOL) {\
				PB_CompressionBuffer __pb_decode_repeated_chars = 0;\
\
				READ_N_BITS(__pb_decode_input_pointer, __pb_decode_buffer, __pb_decode_bits_in_buffer, __pb_decode_repeated_chars, PB_RUN_LENGTH_BIT_SIZE);\
\
				DECODE(__pb_decode_input_pointer, __pb_decode_buffer, __pb_decode_bits_in_buffer, __pb_decode_val, __pb_decode_length, __pb_decode_map);\
				__pb_decode_repeated_chars += PB_MIN_RUN_LENGTH;\
				__pb_decode_current = __pb_decode_map[__pb_decode_val].symbol;\
\
				if (__pb_decode_current == __pb_decode_master_symbol) {\
					__pb_decode_swap_counter--;\
					if (__pb_decode_swap_counter < 0) {\
						PB_PrefixCode __pb_decode_val;\
						int __pb_decode_length;\
\
						DECODE(__pb_decode_input_pointer, __pb_decode_buffer, __pb_decode_bits_in_buffer, __pb_decode_val, __pb_decode_length, __pb_decode_swap_map);\
						__pb_decode_current = __pb_decode_swap_map[__pb_decode_val].symbol;\
\
						READ_N_BITS(__pb_decode_input_pointer, __pb_decode_buffer, __pb_decode_bits_in_buffer, __pb_decode_swap_counter, PB_SWAP_RUN_LENGTH_BIT_SIZE);\
					}\
				}\
\
				PB_DEBUG3(errmsg("SRLE %c (%d) code:%u len:%u i:%d swap_counter:%lu bitsInBuffer:%d buffer:%08X%08X",\
								__pb_decode_current, __pb_decode_current, __pb_decode_val, __pb_decode_length, __pb_decode_i,\
								__pb_decode_repeated_chars, __pb_decode_bits_in_buffer, (uint32) (__pb_decode_buffer >> 32), (uint32) __pb_decode_buffer));\
\
				__pb_decode_n_rle_out = __pb_decode_repeated_chars - 1;\
			}\
		}\
\
		__pb_decode_output = __pb_decode_current;


#define PB_END_DECODE\
	}\
\
	pfree((PB_DecodingMap*) __pb_decode_map);\
	pfree((PB_DecodingMap*) __pb_decode_swap_map);\
\
	PB_TRACE(errmsg("<-PB_BEGIN_DECODE()"))\
}

#endif /* SEQUENCE_DECOMPRESSION_ITERATION_H_ */
