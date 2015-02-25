/*-------------------------------------------------------------------------
*
* Copyright (c) 2013, Max Planck Institute for Marine Microbiology
*
* This software is released under the PostgreSQL License
*
* Author: Michael Schneider <mschneid@mpi-bremen.de>
*
* IDENTIFICATION
*   src/sequence/compression.c
*
*-------------------------------------------------------------------------
*/
#include "postgres.h"
#include "fmgr.h"
#include "access/tuptoaster.h"
#include "c.h"

#include "sequence/sequence.h"
#include "sequence/stats.h"
#include "utils/debug.h"

#include "sequence/compression.h"

/*
 * local types
 */

/**
 * Type for decoding map elements, so the symbol and the length
 * of the code can be stored close to each other in one data
 * structure. An actual decoding map would be an array of
 * PB_DecodingMap of size PB_DECODE_MAP_SIZE (256 if PB_PrefixCode is
 * uint8)
  */
typedef struct {
	uint8 symbol;
	uint8 code_length;
} PB_DecodingMap;

#define PB_DECODE_MAP_SIZE (1 << PB_PREFIX_CODE_BIT_SIZE)

/*
 * Type for encoding map elements, so the code and its length
 * can be stored close to each other in one data structure.
 * An actual encoding map would be an array of PB_EncodingMap of size
 * PB_SOURCE_ALPHABET_SIZE.
 *
 * For encoding the codes will be stored right-aligned!
 */
typedef struct {
	PB_PrefixCode code;
	uint8 code_length;
} PB_EncodingMap;

/*
 * local function declarations
 */

static PB_EncodingMap* get_encoding_map(const PB_CodeSet* codeset, int mode);
static PB_DecodingMap* get_decoding_map(const PB_CodeSet* codeset, int mode);

static void encode_pc(uint8* input,
					  PB_CompressedSequence* output,
					  PB_CodeSet* codeset);
static void encode_pc_idx(uint8* input,
						  PB_CompressedSequence* output,
						  PB_CodeSet* codeset);
static void encode_pc_rle(uint8* input,
						  PB_CompressedSequence* output,
						  PB_CodeSet* codeset);
static void encode_pc_rle_idx(uint8* input,
							  PB_CompressedSequence* output,
							  PB_CodeSet* codeset);
static void encode_pc_swp(uint8* input,
						  PB_CompressedSequence* output,
						  PB_CodeSet* codeset);
static void encode_pc_swp_idx(uint8* input,
							  PB_CompressedSequence* output,
							  PB_CodeSet* codeset);
static void encode_pc_swp_rle(uint8* input,
							  PB_CompressedSequence* output,
							  PB_CodeSet* codeset);
static void encode_pc_swp_rle_idx(uint8* input,
								  PB_CompressedSequence* output,
								  PB_CodeSet* codeset);

static void decode_pc_idx(Varlena* input,
						  uint8* output,
						  uint32 start_position,
						  uint32 output_length,
						  PB_IndexEntry* start_entry,
						  PB_CodeSet* codeset);
static void decode_pc_rle_idx(Varlena* input,
							  uint8* output,
							  uint32 start_position,
							  uint32 output_length,
							  PB_IndexEntry* start_entry,
							  PB_CodeSet* codeset);
static void decode_pc_swp_idx(Varlena* input,
							  uint8* output,
							  uint32 start_position,
							  uint32 output_length,
							  PB_IndexEntry* start_entry,
							  PB_CodeSet* codeset);
static void decode_pc_swp_rle_idx(Varlena* input,
								  uint8* output,
								  uint32 start_position,
								  uint32 output_length,
								  PB_IndexEntry* start_entry,
								  PB_CodeSet* codeset);


/*
 * macros
 */

/**
 * This macro writes a character to a compression buffer. It can
 * only be used if the next character will definitely fit into
 * the buffer. This way there is one jump less than there would be
 * with the ENCODE-macro. Use it to speed up encoding by ca. 5%.
 *
 * Parameters:
 *  PB_PrefixCode code	: the code of the next char
 * 	int code_length : its length
 * 	PB_CompressionBuffer buffer : compression buffer
 * 	int bits_free : free bits in compression buffer
 * 	uint8* input_pointer : points to current position of input sequence
 * 	int i : counts down to end of sequence
 */
#define BURST(code,code_length,buffer,bits_free,input_pointer,i) {	\
	buffer = (buffer << code_length) | code;\
	bits_free -= code_length;\
	input_pointer++;\
	i--;\
}

/**
 * This macro writes a character to the compression buffer. If the
 * compression buffer is full, the code will be split up and the buffer
 * will be written to memory and cleared. Whatever did not fit in the
 * buffer at first will be written to the buffer after it has been
 * cleared.
 *
 * Parameters:
 *  PB_PrefixCode code : code to be written to buffer
 * 	int code_length : length of code to be written to buffer
 * 	PB_CompressionBuffer buffer : compression buffer
 * 	int bits_free : free bits remaining in compression buffer
 * 	PB_CompressionBuffer* output_pointer : compressed sequence
 */
#define ENCODE(code, code_length, buffer, bits_free, output_pointer) { \
	if (code_length <= bits_free) { \
		buffer = (buffer << code_length) | code; \
		bits_free -= code_length; \
	} else { \
		buffer = (buffer << bits_free) | code >> (code_length - bits_free); \
		*output_pointer = buffer; \
		output_pointer++; \
		bits_free = bits_free - code_length + PB_COMPRESSION_BUFFER_BIT_SIZE; \
		buffer = code; \
	} \
}

/**
 * Does the same as ENCODE but binary ORs on output (required to write to stream
 * that has already been written to.
 *
 * Slower than ENCODE.
 */
#define ENCODE_OR(code, code_length, buffer, bits_free, output_pointer) { \
	if (code_length <= bits_free) { \
		buffer = (buffer << code_length) | code; \
		bits_free -= code_length; \
	} else { \
		buffer = (buffer << bits_free) | code >> (code_length - bits_free); \
		*output_pointer |= buffer; \
		output_pointer++; \
		bits_free = bits_free - code_length + PB_COMPRESSION_BUFFER_BIT_SIZE; \
		buffer = code; \
	} \
}

/**
 * This macro reads one character from a compression buffer. If there
 * are not enough bits in buffer to match a prefix code, the next block
 * will be read from memory, and the code will be combined.
 *
 * Parameters:
 * 	PB_CompressionBuffer* input_pointer : compressed sequence
 * 	PB_CompressionBuffer buffer : compression buffer
 * 	int bits_in_buffer : bits in compression buffer
 *  PB_PrefixCode val : next code
 * 	int length : length of next code
 * 	DecodingMap* map : decoding map
 */
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

/*
 * local functions
 */
#define PB_NO_SWAP_MAP		0
#define PB_SWAP_MAP			1

/**
 * get_encoding_map()
 * 		Creates an encoding map for a prefix code set.
 *
 * 	PB_CodeSet* codeset : codeset to create map from
 * 	int mode : PB_ENCODING_MAP_NO_SWAP returns normal map
 * 			   PB_ENCODING_MAP_SWAP returns swap map
 */
static PB_EncodingMap* get_encoding_map(const PB_CodeSet* codeset, int mode)
{
	PB_EncodingMap* map;
	int i;
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

	map = (PB_EncodingMap*) palloc0(PB_ASCII_SIZE * sizeof(PB_EncodingMap));
	memset(map, 0xFF, PB_ASCII_SIZE * sizeof(PB_EncodingMap));

	if (codeset->ignore_case)
	{
		for (i = from; i < to; i++)
		{
			const int symbol = codeset->words[i].symbol;
			const int upper = TO_UPPER(symbol);
			const int lower = TO_LOWER(symbol);

			map[upper].code = codeset->words[i].code >> (PB_PREFIX_CODE_BIT_SIZE - codeset->words[i].code_length);
			map[upper].code_length = codeset->words[i].code_length;
			map[lower].code = codeset->words[i].code >> (PB_PREFIX_CODE_BIT_SIZE - codeset->words[i].code_length);
			map[lower].code_length = codeset->words[i].code_length;

			PB_DEBUG2(errmsg("get_encoding_map(): i:%d symbol:%c (%d) code:%u length:%u/%ld",i,codeset->words[i].symbol,codeset->words[i].symbol,codeset->words[i].code,codeset->words[i].code_length, PB_PREFIX_CODE_BIT_SIZE));
		}
	}
	else
	{
		for (i = from; i < to; i++)
		{
			const int j = codeset->words[i].symbol;

			map[j].code = codeset->words[i].code >> (PB_PREFIX_CODE_BIT_SIZE - codeset->words[i].code_length);
			map[j].code_length = codeset->words[i].code_length;

			PB_DEBUG2(errmsg("get_encoding_map(): i:%d symbol:%c (%d) code:%u length:%u/%ld",i,codeset->words[i].symbol,codeset->words[i].symbol,codeset->words[i].code,codeset->words[i].code_length, PB_PREFIX_CODE_BIT_SIZE));
		}
	}

	PB_TRACE(errmsg("<-get_encoding_map()"));

	return map;
}

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

		PB_DEBUG2(errmsg("get_decoding_map(): i:%d lower_bound:%u upper_bound:%u symbol:%c code:%u length:%u/%ld",i,lower_bound,upper_bound,codeset->words[i].symbol,codeset->words[i].code,codeset->words[i].code_length, PB_PREFIX_CODE_BIT_SIZE));
	}

	return map;
}

/**
 * encode_pc()
 * 		Performs simple prefix code encoding.
 */
static void encode_pc(uint8* input,
					  PB_CompressedSequence* output,
					  PB_CodeSet* codeset)
{
	const int burst_size = codeset->max_codeword_length * 4;
	const PB_EncodingMap* map = get_encoding_map(codeset, PB_NO_SWAP_MAP);

	PB_CompressionBuffer buffer = 0;
	int bits_free = PB_COMPRESSION_BUFFER_BIT_SIZE;
	int i = output->sequence_length - 1;

	PB_CompressionBuffer* output_pointer = PB_COMPRESSED_SEQUENCE_STREAM_POINTER(output);
	uint8* input_pointer = input;

	PB_TRACE(errmsg("->encode_pc(), burst_size=%d, len=%d", burst_size, output->sequence_length));

	while (i >= 0)
	{
		if (burst_size <= bits_free && i > 2)
		{
			BURST(map[(uint8) *input_pointer].code, map[(uint8) *input_pointer].code_length, buffer, bits_free, input_pointer, i);
			BURST(map[(uint8) *input_pointer].code, map[(uint8) *input_pointer].code_length, buffer, bits_free, input_pointer, i);
			BURST(map[(uint8) *input_pointer].code, map[(uint8) *input_pointer].code_length, buffer, bits_free, input_pointer, i);
			BURST(map[(uint8) *input_pointer].code, map[(uint8) *input_pointer].code_length, buffer, bits_free, input_pointer, i);

			PB_DEBUG3(errmsg("BURST i:%d buf:%08X%08X bib:%d", i, (unsigned int) (buffer >> 32), (unsigned int) buffer, bits_free));
		}
		else
		{
			const uint8 current = *input_pointer;
			const PB_PrefixCode code = map[current].code;
			const int code_length = map[current].code_length;

			ENCODE(code, code_length, buffer, bits_free, output_pointer);
			i--;
			input_pointer++;

			PB_DEBUG3(errmsg("SINGLE i:%d val:%u len:%u s:%c buf:%08X%08X bib:%d", i, code, code_length, current, (unsigned int) (buffer >> 32), (unsigned int) buffer, bits_free));
		}
	}

	if (bits_free < PB_COMPRESSION_BUFFER_BIT_SIZE)
		*output_pointer = (buffer << bits_free);

	pfree((PB_EncodingMap*) map);

#ifdef DEBUG
	if (((uint8*)output_pointer) - ((uint8*) output) >= VARSIZE(output)
			&& bits_free < PB_COMPRESSION_BUFFER_BIT_SIZE)
	{
		ereport(ERROR,(errmsg("segmentation fault"),
				errhint("Recognized in %s at line %d.", __FILE__, __LINE__),
				errdetail("Output pointer %p is %ld bytes after a block of size %u.",
				output_pointer,
				(((uint8*)output_pointer) - ((uint8*) output))- VARSIZE(output),
				VARSIZE(output))));
	}
#endif

	PB_TRACE(errmsg("<-encode_pc()"));
 }

/**
 * encode_pc_idx()
 *	 	Encode a sequence with a huffman code and an index.
 */
static void encode_pc_idx(uint8* input,
						  PB_CompressedSequence* output,
						  PB_CodeSet* codeset)
{
	const PB_EncodingMap* map = get_encoding_map(codeset, PB_NO_SWAP_MAP);

	PB_CompressionBuffer buffer = 0;
	int bits_free = PB_COMPRESSION_BUFFER_BIT_SIZE;
	int i = output->sequence_length - 1;
	int index_counter = PB_INDEX_PART_SIZE - 1;

	uint8* input_pointer = input;
	PB_CompressionBuffer* stream_start = PB_COMPRESSED_SEQUENCE_STREAM_POINTER(output);
	PB_CompressionBuffer* output_pointer = stream_start;
	PB_IndexEntry* index_pointer = PB_COMPRESSED_SEQUENCE_INDEX_POINTER(output);

	PB_TRACE(errmsg("->encode_pc_idx()\n\toutput at:%p\n\tstream at:%p\n\tindex at:%p", output, stream_start, index_pointer));

	while (i >= 0)
	{
		const uint8 current = *input_pointer;
		const PB_PrefixCode code = map[current].code;
		const int code_length = map[current].code_length;

		index_counter--;
		if (index_counter < 0)
		{
			index_counter += PB_INDEX_PART_SIZE;
			if (bits_free > 0)
			{
				index_pointer->bit = PB_COMPRESSION_BUFFER_BIT_SIZE - bits_free;
				index_pointer->block = output_pointer - stream_start;
			}
			else
			{
				index_pointer->bit = 0;
				index_pointer->block = (output_pointer + 1) - stream_start;
			}
			PB_DEBUG3(errmsg("IDX block:%u bit:%d", index_pointer->block, index_pointer->bit));
			index_pointer++;
		}

		ENCODE(code, code_length, buffer, bits_free, output_pointer);
		input_pointer++;
		i--;

		PB_DEBUG3(errmsg("OUT i:%d val:%u len:%u s:%c buf:%08X%08X bib:%d", i, code, code_length, current, (unsigned int) (buffer >> 32), (unsigned int) buffer, bits_free));
	}

	/*
	 * Flush buffer.
	 */
	if (bits_free < PB_COMPRESSION_BUFFER_BIT_SIZE)
		*output_pointer = (buffer << bits_free);

	pfree((PB_EncodingMap*) map);

#ifdef DEBUG
	if (((uint8*)output_pointer) - ((uint8*) output) >= VARSIZE(output)
			&& bits_free < PB_COMPRESSION_BUFFER_BIT_SIZE)
	{
		ereport(ERROR,(errmsg("segmentation fault"),
				errhint("Recognized in %s at line %d.", __FILE__, __LINE__),
				errdetail("Output pointer %p is %ld bytes after a block of size %u.",
				output_pointer,
				(((uint8*)output_pointer) - ((uint8*) output))- VARSIZE(output),
				VARSIZE(output))));
	}
#endif

	PB_TRACE(errmsg("<-encode_pc_idx(), output at %p (%d bits free)", output_pointer, bits_free));

}

 /**
  * encode_pc_rle()
  * 		Encode a sequence with a huffman code and rle.

  */
static void encode_pc_rle(uint8* input,
						  PB_CompressedSequence* output,
						  PB_CodeSet* codeset)
{
	const PB_EncodingMap* map = get_encoding_map(codeset, PB_NO_SWAP_MAP);
	const int rlecode_length = map[PB_RUN_LENGTH_SYMBOL].code_length;

	PB_CompressionBuffer buffer = 0;
	int bits_free = PB_COMPRESSION_BUFFER_BIT_SIZE;
	int i = output->sequence_length - 1;
	int repeated_chars = 0;

	uint8 recent = 0;

	uint8* input_pointer = input;
	PB_CompressionBuffer* output_pointer = PB_COMPRESSED_SEQUENCE_STREAM_POINTER(output);

	PB_TRACE(errmsg("->encode_pc_rle()"));

	/*
	 * Read first char.
	 */
	recent = *input_pointer;
	input_pointer++;
	i--;

	while (i >= 0)
	{
		const uint8 current = *input_pointer;
		input_pointer++;
		repeated_chars++;
		i--;

		if (current != recent || repeated_chars >= PB_MAX_RUN_LENGTH - 1)
		{
			encode_pc_rle_out:
			if (repeated_chars < PB_MIN_RUN_LENGTH)
			{
				/*
				 * too few repeated chars to write an RLE word -> write out sequentially
				 */
				const PB_PrefixCode code = map[recent].code;
				const int code_length = map[recent].code_length;
					int k = repeated_chars - 1;

				while (k >= 0)
				{
					k--;
					ENCODE(code, code_length, buffer, bits_free, output_pointer);
					PB_DEBUG3(errmsg("OUT i:%d val:%u len:%u s:%c buf:%08X%08X bib:%d", i, code, code_length, recent, (unsigned int) (buffer >> 32), (unsigned int) buffer, bits_free));
				}
			}
			else
			{
				/*
				 * enough repeated chars to write an RLE word
				 */
				const int recent_length = map[recent].code_length;
				const int code_length = recent_length +
									rlecode_length +
									PB_RUN_LENGTH_BIT_SIZE;
				const PB_CompressionBuffer code = ((PB_CompressionBuffer) map[PB_RUN_LENGTH_SYMBOL].code <<
													(
														PB_RUN_LENGTH_BIT_SIZE +
														recent_length
													)
												) |
												(((PB_CompressionBuffer) (repeated_chars - PB_MIN_RUN_LENGTH)) <<
														recent_length
												) |
												((PB_CompressionBuffer) map[recent].code
												);
				ENCODE(code, code_length, buffer, bits_free, output_pointer);
				PB_DEBUG3(errmsg("RLE i:%d val:%lu len:%u s:%c buf:%08X%08X bib:%d", i, code, code_length, recent, (unsigned int) (buffer >> 32), (unsigned int) buffer, bits_free));
			}
			recent = current;
			repeated_chars = 0;
		}
	}

	/*
	 * Last symbol.
	 *
	 * This is ugly. But the alternatives are:
	 * 		a. another cmp + jmp in the loop -> cost 5% performance
	 * 		b. repeating a bunch of lines -> even uglier
	 */
	if (i == (-1)) {
		i--;
		repeated_chars++;
		PB_DEBUG3(errmsg("last symbol"));
		goto encode_pc_rle_out;
	}

	/*
	 * Flush buffer.
	 */
	if (bits_free < PB_COMPRESSION_BUFFER_BIT_SIZE)
		*output_pointer = (buffer << bits_free);

	pfree((PB_EncodingMap*) map);

#ifdef DEBUG
	if (((uint8*)output_pointer) - ((uint8*) output) >= VARSIZE(output)
			&& bits_free < PB_COMPRESSION_BUFFER_BIT_SIZE) {
		ereport(ERROR,(errmsg("segmentation fault"),
				errhint("Recognized in %s at line %d.", __FILE__, __LINE__),
				errdetail("Output pointer %p is %ld bytes after a block of size %u.",
				output_pointer,
				(((uint8*)output_pointer) - ((uint8*) output))- VARSIZE(output),
				VARSIZE(output))));
	}
#endif

	PB_TRACE(errmsg("<-encode_pc_rle()--%d", bits_free));
 }

 /**
  * encode_pc_rle_idx()
  * 		Encode a sequence with a huffman code, run-length encoding and an index
  *
 */
static void encode_pc_rle_idx(uint8* input,
							  PB_CompressedSequence* output,
							  PB_CodeSet* codeset)
{
	const PB_EncodingMap* map = get_encoding_map(codeset, PB_NO_SWAP_MAP);
	const int rlecode_length = map[PB_RUN_LENGTH_SYMBOL].code_length;

	PB_CompressionBuffer buffer = 0;
	int bits_free = PB_COMPRESSION_BUFFER_BIT_SIZE;
	int i = output->sequence_length - 1;
	int index_counter = PB_INDEX_PART_SIZE - 1;
	int repeated_chars = 0;

	uint8 recent = 0;

	uint8* input_pointer = input;
	PB_CompressionBuffer* stream_start = PB_COMPRESSED_SEQUENCE_STREAM_POINTER(output);
	PB_CompressionBuffer* output_pointer = stream_start;
	PB_IndexEntry* index_pointer = PB_COMPRESSED_SEQUENCE_INDEX_POINTER(output);

	PB_TRACE(errmsg("->encode_pc_rle_idx()"));

	/*
	 * Read first char.
	 */
	recent = *input_pointer;
	input_pointer++;
	i--;

	while (i >= 0)
	{
		const uint8 current = *input_pointer;

		repeated_chars++;
		input_pointer++;
		i--;

		if (current != recent || repeated_chars >= PB_MAX_RUN_LENGTH - 1)
		{
			encode_pc_rle_idx_out:
			if (repeated_chars < PB_MIN_RUN_LENGTH)
			{
				/*
				 * too few repeated chars to write an RLE word -> write out sequentially
				 */
				const int code_length = map[recent].code_length;
				const PB_PrefixCode code = map[recent].code;

				int k = repeated_chars - 1;

				while (k >= 0)
				{
					index_counter--;
					if (index_counter < 0)
					{
						index_counter += PB_INDEX_PART_SIZE;
						if (bits_free > 0)
						{
							index_pointer->bit = PB_COMPRESSION_BUFFER_BIT_SIZE - bits_free;
							index_pointer->block = output_pointer - stream_start;
						}
						else
						{
							index_pointer->bit = 0;
							index_pointer->block = (output_pointer + 1)- stream_start;
						}
						PB_DEBUG3(errmsg("IDX block:%u bit:%d", index_pointer->block, index_pointer->bit));
						index_pointer++;
					}

					k--;
					ENCODE(code, code_length, buffer, bits_free, output_pointer);

					PB_DEBUG3(errmsg("OUT i:%d val:%u len:%u s:%c buf:%08X%08X bib:%d", i, code, code_length, recent, (unsigned int) (buffer >> 32), (unsigned int) buffer, bits_free));
				}
			}
			else
			{
				/*
				 * enough repeated chars to write an RLE word
				 */
				const int recent_length = map[recent].code_length;
				const int code_length = recent_length +
									rlecode_length +
									PB_RUN_LENGTH_BIT_SIZE;
				const PB_CompressionBuffer code = ((PB_CompressionBuffer) map[PB_RUN_LENGTH_SYMBOL].code <<
													(
														PB_RUN_LENGTH_BIT_SIZE +
														recent_length
													)
												) |
												(((PB_CompressionBuffer) (repeated_chars - PB_MIN_RUN_LENGTH)) <<
														recent_length
												) |
												((PB_CompressionBuffer) map[recent].code
												);

				if (index_counter < repeated_chars)
				{
					if (bits_free > 0)
					{
						index_pointer->bit = PB_COMPRESSION_BUFFER_BIT_SIZE - bits_free;
						index_pointer->block = output_pointer - stream_start;
					}
					else
					{
						index_pointer->bit = 0;
						index_pointer->block = (output_pointer + 1)- stream_start;
					}
					index_pointer->rle_shift = (uint16) index_counter;
					index_counter += PB_INDEX_PART_SIZE;
					PB_DEBUG3(errmsg("IDX block:%u bit:%d", index_pointer->block, index_pointer->bit));
					index_pointer++;
				}
				index_counter-=repeated_chars;

				ENCODE(code, code_length, buffer, bits_free, output_pointer);

				PB_DEBUG3(errmsg("RLE i:%d val:%lu len:%u s:%c buf:%08X%08X bib:%d", i, code, code_length, recent, (unsigned int) (buffer >> 32), (unsigned int) buffer, bits_free));
			}

			recent = current;
			repeated_chars = 0;
		}
	}

	/*
	 * Last symbol.
	 *
	 * This is ugly. But the alternatives are:
	 * 		a. another cmp + jmp in the loop -> cost 5% performance
	 * 		b. repeating a bunch of lines -> even uglier
	 */
	if (i == (-1)) {
		i--;
		repeated_chars++;
		PB_DEBUG3(errmsg("last symbol"));
		goto encode_pc_rle_idx_out;
	}

	/*
	 * Flush buffer.
	 */
	if (bits_free < PB_COMPRESSION_BUFFER_BIT_SIZE)
		*output_pointer = (buffer << bits_free);

	pfree((PB_EncodingMap*) map);

#ifdef DEBUG
	if (((uint8*)output_pointer) - ((uint8*) output) >= VARSIZE(output)
			&& bits_free < PB_COMPRESSION_BUFFER_BIT_SIZE) {
		ereport(ERROR,(errmsg("segmentation fault"),
				errhint("Recognized in %s at line %d.", __FILE__, __LINE__),
				errdetail("Output pointer %p is %ld bytes after a block of size %u.",
				output_pointer,
				(((uint8*)output_pointer) - ((uint8*) output))- VARSIZE(output),
				VARSIZE(output))));
	}
#endif

	PB_TRACE(errmsg("<-encode_pc_rle_idx()"));
 }

 /**
  * encode_pc_swp()
  * 		Encode a sequence with a huffman code and rare symbol swapping
  *
  */
static void encode_pc_swp(uint8* input,
						  PB_CompressedSequence* output,
						  PB_CodeSet* codeset)
{
	const PB_EncodingMap* master_map = get_encoding_map(codeset, PB_NO_SWAP_MAP);
	const PB_EncodingMap* swap_map = get_encoding_map(codeset, PB_SWAP_MAP);
	const uint8 master_symbol = codeset->words[codeset->n_symbols - codeset->n_swapped_symbols].symbol;
	const PB_PrefixCode master_symbol_code = master_map[master_symbol].code;
	const int master_symbol_code_length = master_map[master_symbol].code_length;

	PB_CompressionBuffer buffer = 0;
	int bits_free = PB_COMPRESSION_BUFFER_BIT_SIZE;
	int i = output->sequence_length - 1;
	int swap_counter = PB_MAX_SWAP_RUN_LENGTH;

	uint8* input_pointer = input;
	PB_CompressionBuffer * output_pointer = PB_COMPRESSED_SEQUENCE_STREAM_POINTER(output);

	PB_CompressionBuffer* swap_pointer;
	int swap_bits;

	PB_TRACE(errmsg("->encode_pc_swp()"));

	/*
	 * First encode number of master symbols to come.
	 */
	swap_pointer = output_pointer;
	swap_bits = bits_free - PB_SWAP_RUN_LENGTH_BIT_SIZE;
	ENCODE(0, PB_SWAP_RUN_LENGTH_BIT_SIZE, buffer, bits_free, output_pointer);

	while (i >= 0)
	{
		const uint8 current = *input_pointer;

		input_pointer++;
		i--;

		if (swap_map[current].code_length == 0xFF)
		{
			/*
			 * Symbols that do not occur in swap code, will be written to stream directly.
			 */
			const PB_PrefixCode code = master_map[current].code;
			const uint8 code_length = master_map[current].code_length;

			ENCODE_OR(code, code_length, buffer, bits_free, output_pointer);
		}
		else
		{
			/*
			 * Current char is master symbol or swapped symbol
			 */
			ENCODE_OR(master_symbol_code, master_symbol_code_length, buffer, bits_free, output_pointer);

			if (current == master_symbol)
			{
				/*
				 * Current char is master symbol
				 */
				swap_counter--;
			}

			if (current != master_symbol || swap_counter < 0)
			{
				/*
				 * Recent is swapped symbol or swap RLE is full.
				 */
				const PB_CompressionBuffer pos = swap_counter < 0 ? PB_MAX_SWAP_RUN_LENGTH : (PB_MAX_SWAP_RUN_LENGTH - swap_counter);
				const PB_PrefixCode code = swap_map[current].code;
				const uint8 code_length = swap_map[current].code_length;

				ENCODE_OR(code, code_length, buffer, bits_free, output_pointer);

				if (swap_bits < 0)
				{
					*swap_pointer = *swap_pointer | (pos >> (-swap_bits));
					*(swap_pointer+1) = *(swap_pointer+1) | (pos << (swap_bits + PB_COMPRESSION_BUFFER_BIT_SIZE));
				}
				else
				{
					*swap_pointer = *swap_pointer | (pos << swap_bits);
				}

				swap_pointer = output_pointer;
				swap_bits = bits_free - PB_SWAP_RUN_LENGTH_BIT_SIZE;
				ENCODE_OR(0, PB_SWAP_RUN_LENGTH_BIT_SIZE, buffer, bits_free, output_pointer);

				swap_counter = PB_MAX_SWAP_RUN_LENGTH;
			}
		}
	}

	/*
	 * Flush buffer.
	 */
	if (bits_free < PB_COMPRESSION_BUFFER_BIT_SIZE)
		*output_pointer |= (buffer << bits_free);

	if (swap_bits < 0)
	{
		*swap_pointer = *swap_pointer | ((PB_CompressionBuffer) 0xFFFF >> (-swap_bits));
		*(swap_pointer+1) = *(swap_pointer+1) | ((PB_CompressionBuffer)0xFFFF << (swap_bits + PB_COMPRESSION_BUFFER_BIT_SIZE));
	}
	else
	{
		*swap_pointer = *swap_pointer | ((PB_CompressionBuffer) 0xFFFF << swap_bits);
	}

	pfree((PB_EncodingMap*) master_map);
	pfree((PB_EncodingMap*) swap_map);

#ifdef DEBUG
	if (((uint8*)output_pointer) - ((uint8*) output) >= VARSIZE(output)
			&& bits_free < PB_COMPRESSION_BUFFER_BIT_SIZE) {
		ereport(ERROR,(errmsg("segmentation fault"),
				errhint("Recognized in %s at line %d.", __FILE__, __LINE__),
				errdetail("Output pointer %p is %ld bytes after a block of size %u.",
				output_pointer,
				(((uint8*)output_pointer) - ((uint8*) output))- VARSIZE(output),
				VARSIZE(output))));
	}
#endif

	PB_TRACE(errmsg("<-encode_pc_swp()"));
 }

 /**
  * encode_pc_swp_idx()
  * 		Encode a sequence with a huffman code, rare symbol swapping and an index
  *
  */
static void encode_pc_swp_idx(uint8* input,
							  PB_CompressedSequence* output,
							  PB_CodeSet* codeset)
{
	const PB_EncodingMap* master_map = get_encoding_map(codeset, PB_NO_SWAP_MAP);
	const PB_EncodingMap* swap_map = get_encoding_map(codeset, PB_SWAP_MAP);
	const uint8 master_symbol = codeset->words[codeset->n_symbols - codeset->n_swapped_symbols].symbol;
	const PB_PrefixCode master_symbol_code = master_map[master_symbol].code;
	const int master_symbol_code_length = master_map[master_symbol].code_length;

	PB_CompressionBuffer buffer = 0;
	int bits_free = PB_COMPRESSION_BUFFER_BIT_SIZE;
	int i = output->sequence_length - 1;
	int swap_counter = PB_MAX_SWAP_RUN_LENGTH;
	int index_counter = PB_INDEX_PART_SIZE - 1;

	uint8* input_pointer = input;
	PB_CompressionBuffer* stream_start = PB_COMPRESSED_SEQUENCE_STREAM_POINTER(output);
	PB_CompressionBuffer* output_pointer = stream_start;
	PB_IndexEntry* index_pointer = PB_COMPRESSED_SEQUENCE_INDEX_POINTER(output);
	int n_swap_index_pointers = 0;

	PB_CompressionBuffer* swap_pointer;
	int swap_bits;

	PB_TRACE(errmsg("->encode_pc_swp_idx(), out:%p stream:%p", output, stream_start));

	/*
	 * First encode number of master symbols to come.
	 */
	swap_pointer = output_pointer;
	swap_bits = bits_free - PB_SWAP_RUN_LENGTH_BIT_SIZE;
	ENCODE(0, PB_SWAP_RUN_LENGTH_BIT_SIZE, buffer, bits_free, output_pointer);

	while (i >= 0)
	{
		const uint8 current = *input_pointer;

		index_counter--;
		if (index_counter < 0)
		{
			index_counter += PB_INDEX_PART_SIZE;
			if (bits_free > 0)
			{
				index_pointer->bit = PB_COMPRESSION_BUFFER_BIT_SIZE - bits_free;
				index_pointer->block = output_pointer - stream_start;
			}
			else
			{
				index_pointer->bit = 0;
				index_pointer->block = (output_pointer + 1)- stream_start;
			}
			index_pointer->swap_shift = swap_counter;
			n_swap_index_pointers++;
			index_pointer++;
		}

		input_pointer++;
		i--;

		if (swap_map[current].code_length == 0xFF)
		{
			/*
			 * Symbols that do not occur in swap code, will be written to stream directly.
			 */
			const PB_PrefixCode code = master_map[current].code;
			const uint8 code_length = master_map[current].code_length;

			ENCODE_OR(code, code_length, buffer, bits_free, output_pointer);
		}
		else
		{
			/*
			 * Current char is master symbol or swapped symbol
			 */
			ENCODE_OR(master_symbol_code, master_symbol_code_length, buffer, bits_free, output_pointer);

			if (current == master_symbol)
			{
				/*
				 * Current char is master symbol
				 */
				swap_counter--;
			}

			if (current != master_symbol || swap_counter < 0)
			{
				/*
				 * Recent is swapped symbol or swap RLE is full.
				 */
				const PB_CompressionBuffer pos = swap_counter < 0 ? PB_MAX_SWAP_RUN_LENGTH : (PB_MAX_SWAP_RUN_LENGTH - swap_counter);
				const PB_PrefixCode code = swap_map[current].code;
				const uint8 code_length = swap_map[current].code_length;

				ENCODE_OR(code, code_length, buffer, bits_free, output_pointer);

				if (n_swap_index_pointers > 0)
				{
					int index_entry_it;
					for (index_entry_it = 1; index_entry_it <= n_swap_index_pointers; index_entry_it++)
					{
						(index_pointer - index_entry_it)->swap_shift -= swap_counter;
					}
				}
				n_swap_index_pointers = 0;

				if (swap_bits < 0)
				{
					*swap_pointer = *swap_pointer | (pos >> (-swap_bits));
					*(swap_pointer+1) = *(swap_pointer+1) | (pos << (swap_bits + PB_COMPRESSION_BUFFER_BIT_SIZE));
				}
				else
				{
					*swap_pointer = *swap_pointer | (pos << swap_bits);
				}

				swap_pointer = output_pointer;
				swap_bits = bits_free - PB_SWAP_RUN_LENGTH_BIT_SIZE;
				ENCODE_OR(0, PB_SWAP_RUN_LENGTH_BIT_SIZE, buffer, bits_free, output_pointer);

				swap_counter = PB_MAX_SWAP_RUN_LENGTH;
			}
		}
	}

	/*
	 * Flush buffer.
	 */
	if (bits_free < PB_COMPRESSION_BUFFER_BIT_SIZE)
		*output_pointer |= (buffer << bits_free);

	if (swap_bits < 0)
	{
		*swap_pointer = *swap_pointer | ((PB_CompressionBuffer) 0xFFFF >> (-swap_bits));
		*(swap_pointer+1) = *(swap_pointer+1) | ((PB_CompressionBuffer)0xFFFF << (swap_bits + PB_COMPRESSION_BUFFER_BIT_SIZE));
	}
	else
	{
		*swap_pointer = *swap_pointer | ((PB_CompressionBuffer) 0xFFFF << swap_bits);
	}

	if (n_swap_index_pointers > 0)
	{
		int index_entry_it;
		for (index_entry_it = 1; index_entry_it <= n_swap_index_pointers; index_entry_it++)
		{
			(index_pointer - index_entry_it)->swap_shift -= swap_counter;
		}
	}

	pfree((PB_EncodingMap*) master_map);
	pfree((PB_EncodingMap*) swap_map);

#ifdef DEBUG
	if (((uint8*)output_pointer) - ((uint8*) output) >= VARSIZE(output)
			&& bits_free < PB_COMPRESSION_BUFFER_BIT_SIZE) {
		ereport(ERROR,(errmsg("segmentation fault"),
				errhint("Recognized in %s at line %d.", __FILE__, __LINE__),
				errdetail("Output pointer %p is %ld bytes after a block of size %u starting at %p.",
				output_pointer,
				(((uint8*)output_pointer) - ((uint8*) output)) - VARSIZE(output),
				VARSIZE(output),
				output)));
	}
#endif

	PB_TRACE(errmsg("<-encode_pc_swp_idx()"));
 }

/**
 * encode_pc_swp_rle()
 * 		Encode a sequence with a huffman code, run-length encoding and
 * 		rare symbol swapping.
 *
 */
static void encode_pc_swp_rle(uint8* input,
							  PB_CompressedSequence* output,
							  PB_CodeSet* codeset)
{
	const PB_EncodingMap* master_map = get_encoding_map(codeset, PB_NO_SWAP_MAP);
	const PB_EncodingMap* swap_map = get_encoding_map(codeset, PB_SWAP_MAP);

	const uint8 master_symbol = codeset->words[codeset->n_symbols - codeset->n_swapped_symbols].symbol;
	const PB_PrefixCode master_symbol_code = master_map[master_symbol].code;
	const int master_symbol_code_length = master_map[master_symbol].code_length;

	const bool rle_is_swapped = master_map[PB_RUN_LENGTH_SYMBOL].code_length == 0xFF ? true : false;
	const PB_PrefixCode rlecode = rle_is_swapped ? master_map[master_symbol].code : master_map[PB_RUN_LENGTH_SYMBOL].code;
	const int rlecode_length = rle_is_swapped ? master_map[master_symbol].code_length : master_map[PB_RUN_LENGTH_SYMBOL].code_length;

	PB_CompressionBuffer buffer = 0;
	int bits_free = PB_COMPRESSION_BUFFER_BIT_SIZE;
	int i = output->sequence_length - 1;
	int swap_counter = PB_MAX_SWAP_RUN_LENGTH;
	int repeated_chars = 0;

	uint8 recent = 0;

	uint8* input_pointer = input;
	PB_CompressionBuffer* stream_start = PB_COMPRESSED_SEQUENCE_STREAM_POINTER(output);
	PB_CompressionBuffer* output_pointer = stream_start;

	PB_CompressionBuffer* swap_pointer;
	int swap_bits;

	PB_TRACE(errmsg("->encode_pc_swp_rle()"));

	/*
	 * First encode number of master symbols to come.
	 */
	swap_pointer = output_pointer;
	swap_bits = bits_free - PB_SWAP_RUN_LENGTH_BIT_SIZE;
	ENCODE(0, PB_SWAP_RUN_LENGTH_BIT_SIZE, buffer, bits_free, output_pointer);

	/*
	 * Read first char.
	 */
	recent = *input_pointer;
	input_pointer++;
	i--;

	while (i >= 0)
	{
		const uint8 current = *input_pointer;

		repeated_chars++;
		input_pointer++;
		i--;

		if (current != recent || repeated_chars >= PB_MAX_RUN_LENGTH - 1)
		{
			encode_pc_swp_rle_out:
			if (repeated_chars < PB_MIN_RUN_LENGTH)
			{
				/*
				 * too few repeated chars to write an RLE word -> write out sequentially
				 */
				int write_out_chars = repeated_chars - 1;

				if (swap_map[recent].code_length == 0xFF)
				{
					/*
					 * recent is a non-swapped, non-master symbol
					 */
					const int code_length = master_map[recent].code_length;
					const PB_PrefixCode code = master_map[recent].code;

					PB_DEBUG3(errmsg("WOSNS %c (%d) code:%u len:%u i:%d index_counter:%d bits_free:%d outat:%ld buffer:%08X%08X", recent, recent,code, code_length,  i, write_out_chars + 1, bits_free, (output_pointer - stream_start) * 8, (uint32) (buffer >> 32), (uint32) buffer));

					while (write_out_chars >= 0)
					{
						write_out_chars--;
						ENCODE_OR(code, code_length, buffer, bits_free, output_pointer);
					}
				}
				else
				{
					/*
					 * recent is swapped symbol or master symbol
					 */
					const int code_length = swap_map[recent].code_length;
					const PB_PrefixCode code = swap_map[recent].code;

					PB_DEBUG3(errmsg("WOSS1 %c (%d) code:%u len:%u i:%d index_counter:%d bits_free:%d buffer:%08X%08X", master_symbol, master_symbol, master_symbol_code, master_symbol_code_length,  i, write_out_chars + 1, bits_free, (uint32) (buffer >> 32), (uint32) buffer));
					PB_DEBUG3(errmsg("WOSS2 %c (%d) code:%u len:%u i:%d l:%d bits_free:%d buffer:%08X%08X", recent, recent, code, code_length,  i, swap_counter, bits_free, (uint32) (buffer >> 32), (uint32) buffer));

					while (write_out_chars >= 0)
					{
						ENCODE_OR(master_symbol_code, master_symbol_code_length, buffer, bits_free, output_pointer);

						write_out_chars--;

						if (recent == master_symbol)
						{
							swap_counter--;
						}

						if (recent != master_symbol || swap_counter < 0)
						{
							/*
							 * Recent is swapped symbol or swap RLE is full.
							 */
							const PB_CompressionBuffer pos = swap_counter < 0 ? PB_MAX_SWAP_RUN_LENGTH : (PB_MAX_SWAP_RUN_LENGTH - swap_counter);

							ENCODE_OR(code, code_length, buffer, bits_free, output_pointer);

							if (swap_bits < 0)
							{
								*swap_pointer = *swap_pointer | (pos >> (-swap_bits));
								*(swap_pointer+1) = *(swap_pointer+1) | (pos << (swap_bits + PB_COMPRESSION_BUFFER_BIT_SIZE));
							}
							else
							{
								*swap_pointer = *swap_pointer | (pos << swap_bits);
							}

							swap_pointer = output_pointer;
							swap_bits = bits_free - PB_SWAP_RUN_LENGTH_BIT_SIZE;
							ENCODE_OR(0, PB_SWAP_RUN_LENGTH_BIT_SIZE, buffer, bits_free, output_pointer);

							swap_counter = PB_MAX_SWAP_RUN_LENGTH;
						}
					}
				}
			}
			else
			{
				/*
				 * 1. Write out RLE symbol
				 */
				PB_DEBUG3(errmsg("RLE %c (%d) code:%u len:%u i:%d index_counter:%d bits_free:%d buffer:%08X%08X", master_symbol, master_symbol,rlecode, rlecode_length,  i, repeated_chars, bits_free, (uint32) (buffer >> 32), (uint32) buffer));

				ENCODE_OR(rlecode,rlecode_length,buffer,bits_free,output_pointer);

				if (rle_is_swapped)
				{
					/*
					 * RLE-symbol is swapped but not master symbol
					 *  -> write out swap info
					 */
					const PB_CompressionBuffer pos = (PB_MAX_SWAP_RUN_LENGTH - swap_counter);
					const int code_length = swap_map[PB_RUN_LENGTH_SYMBOL].code_length;
					const PB_PrefixCode code = swap_map[PB_RUN_LENGTH_SYMBOL].code;

					ENCODE_OR(code, code_length, buffer, bits_free, output_pointer);

					if (swap_bits < 0)
					{
						*swap_pointer = *swap_pointer | (pos >> (-swap_bits));
						*(swap_pointer+1) = *(swap_pointer+1) | (pos << (swap_bits + PB_COMPRESSION_BUFFER_BIT_SIZE));
					}
					else
					{
						*swap_pointer = *swap_pointer | (pos << swap_bits);
					}

					swap_pointer = output_pointer;
					swap_bits = bits_free - PB_SWAP_RUN_LENGTH_BIT_SIZE;
					ENCODE_OR(0, PB_SWAP_RUN_LENGTH_BIT_SIZE, buffer, bits_free, output_pointer);
					swap_counter = PB_MAX_SWAP_RUN_LENGTH;
				}
				else if (master_symbol == PB_RUN_LENGTH_SYMBOL)
				{
					/*
					 * RLE symbol is master symbol
					 */
					swap_counter--;

					if (swap_counter < 0)
					{
						const PB_CompressionBuffer pos = PB_MAX_SWAP_RUN_LENGTH;
						const PB_PrefixCode code = swap_map[master_symbol].code;
						const uint8 code_length = swap_map[master_symbol].code_length;

						ENCODE_OR(code, code_length, buffer, bits_free, output_pointer);

						if (swap_bits < 0)
						{
							*swap_pointer = *swap_pointer | (pos >> (-swap_bits));
							*(swap_pointer+1) = *(swap_pointer+1) | (pos << (swap_bits + PB_COMPRESSION_BUFFER_BIT_SIZE));
						}
						else
						{
							*swap_pointer = *swap_pointer | (pos << swap_bits);
						}

						swap_pointer = output_pointer;
						swap_bits = bits_free - PB_SWAP_RUN_LENGTH_BIT_SIZE;
						ENCODE_OR(0, PB_SWAP_RUN_LENGTH_BIT_SIZE, buffer, bits_free, output_pointer);
						swap_counter = PB_MAX_SWAP_RUN_LENGTH;
					}
				}

				/*
				 * 2. Write out run-length
				 */
				ENCODE_OR((PB_RunLength) (repeated_chars - PB_MIN_RUN_LENGTH), sizeof(PB_RunLength) * 8, buffer, bits_free, output_pointer);

				/**
				 * 3. Write out symbol
				 */
				if (swap_map[recent].code_length == 0xFF)
				{
					/*
					 * Recent is neither swapped nor master.
					 */
					const PB_PrefixCode code = master_map[recent].code;
					const uint8 code_length = master_map[recent].code_length;

					ENCODE_OR(code, code_length, buffer, bits_free, output_pointer);
				}
				else
				{
					/*
					 * Recent IS either swapped or master.
					 */
					ENCODE_OR(master_symbol_code, master_symbol_code_length, buffer, bits_free, output_pointer);

					if (master_symbol == recent)
					{
						swap_counter--;
					}

					if (recent != master_symbol || swap_counter < 0)
					{
						/*
						 * Recent is swapped symbol or swap RLE is full.
						 */
						const PB_CompressionBuffer pos = swap_counter < 0 ? PB_MAX_SWAP_RUN_LENGTH : (PB_MAX_SWAP_RUN_LENGTH - swap_counter);
						const int code_length = swap_map[recent].code_length;
						const PB_PrefixCode code = swap_map[recent].code;

						ENCODE_OR(code, code_length, buffer, bits_free, output_pointer);

						if (swap_bits < 0)
						{
							*swap_pointer = *swap_pointer | (pos >> (-swap_bits));
							*(swap_pointer+1) = *(swap_pointer+1) | (pos << (swap_bits + PB_COMPRESSION_BUFFER_BIT_SIZE));
						}
						else
						{
							*swap_pointer = *swap_pointer | (pos << swap_bits);
						}

						swap_pointer = output_pointer;
						swap_bits = bits_free - PB_SWAP_RUN_LENGTH_BIT_SIZE;
						ENCODE_OR(0, PB_SWAP_RUN_LENGTH_BIT_SIZE, buffer, bits_free, output_pointer);

						swap_counter = PB_MAX_SWAP_RUN_LENGTH;
					}
				}
			}

			recent = current;
			repeated_chars = 0;
		}
	}

	/*
	 * Last symbol.
	 *
	 * This is ugly. But the alternatives are:
	 * 		a. another cmp + jmp in the loop -> cost 5% performance
	 * 		b. repeating a bunch of lines -> even uglier
	 */
	if (i == (-1)) {
		i--;
		repeated_chars++;
		PB_DEBUG3(errmsg("last symbol"));
		goto encode_pc_swp_rle_out;
	}

	/*
	 * Flush buffer.
	 */
	if (bits_free < PB_COMPRESSION_BUFFER_BIT_SIZE)
		*output_pointer |= (buffer << bits_free);

	if (swap_bits < 0)
	{
		*swap_pointer = *swap_pointer | ((PB_CompressionBuffer) 0xFFFF >> (-swap_bits));
		*(swap_pointer+1) = *(swap_pointer+1) | ((PB_CompressionBuffer) 0xFFFF << (swap_bits + PB_COMPRESSION_BUFFER_BIT_SIZE));
	}
	else
	{
		*swap_pointer = *swap_pointer | ((PB_CompressionBuffer) 0xFFFF << swap_bits);
	}

	PB_DEBUG1(errmsg("encode_pc_swp_rle(): pos:%ld", (output_pointer - PB_COMPRESSED_SEQUENCE_STREAM_POINTER(output))*8));

	pfree((PB_EncodingMap*) master_map);
	pfree((PB_EncodingMap*) swap_map);

#ifdef DEBUG
	if (((uint8*)output_pointer) - ((uint8*) output) >= VARSIZE(output)
			&& bits_free < PB_COMPRESSION_BUFFER_BIT_SIZE) {
		ereport(ERROR,(errmsg("segmentation fault"),
				errhint("Recognized in %s at line %d.", __FILE__, __LINE__),
				errdetail("Output pointer %p is %ld bytes after a block of size %u.",
				output_pointer,
				(((uint8*)output_pointer) - ((uint8*) output))- VARSIZE(output),
				VARSIZE(output))));
	}
#endif

	PB_TRACE(errmsg("<-encode_pc_swp_rle()"));
}

/**
 * encode_pc_swp_rle_idx()
 * 		Encode a sequence with a huffman code, run-length encoding,
 * 		rare symbol swapping and an index.
 *
 */
static void encode_pc_swp_rle_idx(uint8* input,
								  PB_CompressedSequence* output,
								  PB_CodeSet* codeset)
{
	const PB_EncodingMap* master_map = get_encoding_map(codeset, PB_NO_SWAP_MAP);
	const PB_EncodingMap* swap_map = get_encoding_map(codeset, PB_SWAP_MAP);

	const uint8 master_symbol = codeset->words[codeset->n_symbols - codeset->n_swapped_symbols].symbol;
	const PB_PrefixCode master_symbol_code = master_map[master_symbol].code;
	const int master_symbol_code_length = master_map[master_symbol].code_length;

	const bool rle_is_swapped = master_map[PB_RUN_LENGTH_SYMBOL].code_length == 0xFF ? true : false;
	const PB_PrefixCode rlecode = rle_is_swapped ? master_map[master_symbol].code : master_map[PB_RUN_LENGTH_SYMBOL].code;
	const int rlecode_length = rle_is_swapped ? master_map[master_symbol].code_length : master_map[PB_RUN_LENGTH_SYMBOL].code_length;

	PB_CompressionBuffer buffer = 0;
	int bits_free = PB_COMPRESSION_BUFFER_BIT_SIZE;
	int i = output->sequence_length - 1;
	int swap_counter = PB_MAX_SWAP_RUN_LENGTH;
	int repeated_chars = 0;

	uint8 recent = 0;

	uint8* input_pointer = input;
	PB_CompressionBuffer* stream_start = PB_COMPRESSED_SEQUENCE_STREAM_POINTER(output);
	PB_CompressionBuffer* output_pointer = stream_start;

	PB_CompressionBuffer* swap_pointer;
	int swap_bits;

	int index_counter = PB_INDEX_PART_SIZE - 1;
	PB_IndexEntry* index_pointer = PB_COMPRESSED_SEQUENCE_INDEX_POINTER(output);
	int n_swap_index_pointers = 0;

	PB_TRACE(errmsg("->encode_pc_swp_rle_idx(): %d chars", i));

	/*
	 * First encode number of master symbols to come.
	 */
	swap_pointer = output_pointer;
	swap_bits = bits_free - PB_SWAP_RUN_LENGTH_BIT_SIZE;
	ENCODE(0, PB_SWAP_RUN_LENGTH_BIT_SIZE, buffer, bits_free, output_pointer);

	/*
	 * Read first char.
	 */
	recent = *input_pointer;
	input_pointer++;
	i--;

	while (i >= 0)
	{
		const uint8 current = *input_pointer;

		repeated_chars++;
		input_pointer++;
		i--;

		if (current != recent || repeated_chars >= PB_MAX_RUN_LENGTH - 1)
		{
			encode_pc_swp_rle_idx_out:
			if (repeated_chars < PB_MIN_RUN_LENGTH)
			{
				/*
				 * too few repeated chars to write an RLE word -> write out sequentially
				 */
				int write_out_chars = repeated_chars - 1;

				if (swap_map[recent].code_length == 0xFF)
				{
					/*
					 * recent is a non-swapped, non-master symbol
					 */
					const int code_length = master_map[recent].code_length;
					const PB_PrefixCode code = master_map[recent].code;

					PB_DEBUG3(errmsg("WOSNS %c (%d) code:%u len:%u i:%d index_counter:%d bits_free:%d outat:%ld buffer:%08X%08X", recent, recent,code, code_length,  i, write_out_chars + 1, bits_free, (output_pointer - stream_start) * 8, (uint32) (buffer >> 32), (uint32) buffer));

					while (write_out_chars >= 0)
					{
						index_counter--;
						if (index_counter < 0)
						{
							index_counter += PB_INDEX_PART_SIZE;
							if (bits_free > 0)
							{
								index_pointer->bit = PB_COMPRESSION_BUFFER_BIT_SIZE - bits_free;;
								index_pointer->block = output_pointer - stream_start;
							}
							else
							{
								index_pointer->bit = 0;
								index_pointer->block = (output_pointer + 1)- stream_start;
							}
							index_pointer->swap_shift = swap_counter;
							n_swap_index_pointers++;
							index_pointer++;
						}

						write_out_chars--;
						ENCODE_OR(code, code_length, buffer, bits_free, output_pointer);
					}
				}
				else
				{
					/*
					 * recent is swapped symbol or master symbol
					 */
					const int code_length = swap_map[recent].code_length;
					const PB_PrefixCode code = swap_map[recent].code;

					PB_DEBUG3(errmsg("WOSS1 %c (%d) code:%u len:%u i:%d index_counter:%d bits_free:%d buffer:%08X%08X", master_symbol, master_symbol, master_symbol_code, master_symbol_code_length,  i, write_out_chars + 1, bits_free, (uint32) (buffer >> 32), (uint32) buffer));
					PB_DEBUG3(errmsg("WOSS2 %c (%d) code:%u len:%u i:%d l:%d bits_free:%d buffer:%08X%08X", recent, recent,code, code_length,  i, swap_counter, bits_free, (uint32) (buffer >> 32), (uint32) buffer));

					while (write_out_chars >= 0)
					{
						index_counter--;
						if (index_counter < 0)
						{
							index_counter += PB_INDEX_PART_SIZE;
							if (bits_free > 0)
							{
								index_pointer->bit = PB_COMPRESSION_BUFFER_BIT_SIZE - bits_free;;
								index_pointer->block = output_pointer - stream_start;
							}
							else
							{
								index_pointer->bit = 0;
								index_pointer->block = (output_pointer + 1)- stream_start;
							}
							index_pointer->swap_shift = swap_counter;
							n_swap_index_pointers++;
							index_pointer++;
						}

						ENCODE_OR(master_symbol_code, master_symbol_code_length, buffer, bits_free, output_pointer);

						write_out_chars--;

						if (recent == master_symbol) {
							swap_counter--;
						}

						if (recent != master_symbol || swap_counter < 0)
						{
							/*
							 * Recent is swapped symbol or swap RLE is full.
							 */
							const PB_CompressionBuffer pos = swap_counter < 0 ? PB_MAX_SWAP_RUN_LENGTH : (PB_MAX_SWAP_RUN_LENGTH - swap_counter);
							int index_entry_it;

							ENCODE_OR(code, code_length, buffer, bits_free, output_pointer);

							for (index_entry_it = 1; index_entry_it <= n_swap_index_pointers; index_entry_it++)
							{
								(index_pointer - index_entry_it)->swap_shift -= swap_counter;
							}
							n_swap_index_pointers = 0;

							if (swap_bits < 0)
							{
								*swap_pointer = *swap_pointer | (pos >> (-swap_bits));
								*(swap_pointer+1) = *(swap_pointer+1) | (pos << (swap_bits + PB_COMPRESSION_BUFFER_BIT_SIZE));
							}
							else
							{
								*swap_pointer = *swap_pointer | (pos << swap_bits);
							}

							swap_pointer = output_pointer;
							swap_bits = bits_free - PB_SWAP_RUN_LENGTH_BIT_SIZE;
							ENCODE_OR(0, PB_SWAP_RUN_LENGTH_BIT_SIZE, buffer, bits_free, output_pointer);

							swap_counter = PB_MAX_SWAP_RUN_LENGTH;
						}
					}
				}
			}
			else
			{
				/*
				 * 1. Write out index entry
				 */
				if (index_counter < repeated_chars)
				{
					PB_DEBUG3(errmsg("idx - i1:%d c:%d rc:%d", output->sequence_length - i, index_counter, repeated_chars));
					if (bits_free > 0)
					{
						index_pointer->bit = PB_COMPRESSION_BUFFER_BIT_SIZE - bits_free;
						index_pointer->block = output_pointer - stream_start;
					}
					else
					{
						index_pointer->bit = 0;
						index_pointer->block = (output_pointer + 1)- stream_start;
					}
					index_pointer->rle_shift = (uint16) index_counter;
					index_pointer->swap_shift = swap_counter;
					n_swap_index_pointers++;
					index_pointer++;
					index_counter += PB_INDEX_PART_SIZE;
				}
				index_counter -= repeated_chars;

				/*
				 * 2. Write out RLE symbol
				 */
				PB_DEBUG3(errmsg("RLE %c (%d) code:%u len:%u i:%d index_counter:%d bits_free:%d buffer:%08X%08X", master_symbol, master_symbol,rlecode, rlecode_length,  i, repeated_chars, bits_free, (uint32) (buffer >> 32), (uint32) buffer));

				ENCODE_OR(rlecode,rlecode_length,buffer,bits_free,output_pointer);

				if (rle_is_swapped)
				{
					/*
					 * RLE-symbol is swapped but not master symbol
					 *  -> write out swap info
					 */
					const PB_CompressionBuffer pos = (PB_MAX_SWAP_RUN_LENGTH - swap_counter);
					const int code_length = swap_map[PB_RUN_LENGTH_SYMBOL].code_length;
					const PB_PrefixCode code = swap_map[PB_RUN_LENGTH_SYMBOL].code;

					int index_entry_it;

					ENCODE_OR(code, code_length, buffer, bits_free, output_pointer);

					for (index_entry_it = 1; index_entry_it <= n_swap_index_pointers; index_entry_it++)
					{
						(index_pointer - index_entry_it)->swap_shift -= swap_counter;
					}
					n_swap_index_pointers = 0;

					if (swap_bits < 0)
					{
						*swap_pointer = *swap_pointer | (pos >> (-swap_bits));
						*(swap_pointer+1) = *(swap_pointer+1) | (pos << (swap_bits + PB_COMPRESSION_BUFFER_BIT_SIZE));
					}
					else
					{
						*swap_pointer = *swap_pointer | (pos << swap_bits);
					}

					swap_pointer = output_pointer;
					swap_bits = bits_free - PB_SWAP_RUN_LENGTH_BIT_SIZE;
					ENCODE_OR(0, PB_SWAP_RUN_LENGTH_BIT_SIZE, buffer, bits_free, output_pointer);
					swap_counter = PB_MAX_SWAP_RUN_LENGTH;
				}
				else if (master_symbol == PB_RUN_LENGTH_SYMBOL)
				{
					/*
					 * RLE symbol is master symbol
					 */
					swap_counter--;

					if (swap_counter < 0)
					{
						const PB_CompressionBuffer pos = PB_MAX_SWAP_RUN_LENGTH;
						const PB_PrefixCode code = swap_map[master_symbol].code;
						const uint8 code_length = swap_map[master_symbol].code_length;

						int index_entry_it;

						ENCODE_OR(code, code_length, buffer, bits_free, output_pointer);

						for (index_entry_it = 1; index_entry_it <= n_swap_index_pointers; index_entry_it++)
						{
							(index_pointer - index_entry_it)->swap_shift -= swap_counter;
						}
						n_swap_index_pointers = 0;

						if (swap_bits < 0)
						{
							*swap_pointer = *swap_pointer | (pos >> (-swap_bits));
							*(swap_pointer+1) = *(swap_pointer+1) | (pos << (swap_bits + PB_COMPRESSION_BUFFER_BIT_SIZE));
						}
						else
						{
							*swap_pointer = *swap_pointer | (pos << swap_bits);
						}

						swap_pointer = output_pointer;
						swap_bits = bits_free - PB_SWAP_RUN_LENGTH_BIT_SIZE;
						ENCODE_OR(0, PB_SWAP_RUN_LENGTH_BIT_SIZE, buffer, bits_free, output_pointer);
						swap_counter = PB_MAX_SWAP_RUN_LENGTH;
					}
				}

				/*
				 * 3. Write out run-length
				 */
				ENCODE_OR((PB_RunLength) (repeated_chars - PB_MIN_RUN_LENGTH), sizeof(PB_RunLength) * 8, buffer, bits_free, output_pointer);

				/**
				 * 4. Write out symbol
				 */
				if (swap_map[recent].code_length == 0xFF)
				{
					/*
					 * Recent is neither swapped nor master.
					 */
					const PB_PrefixCode code = master_map[recent].code;
					const uint8 code_length = master_map[recent].code_length;

					ENCODE_OR(code, code_length, buffer, bits_free, output_pointer);
				}
				else
				{
					/*
					 * Recent IS either swapped or master.
					 */
					ENCODE_OR(master_symbol_code, master_symbol_code_length, buffer, bits_free, output_pointer);

					if (master_symbol == recent)
					{
						swap_counter--;
					}

					if (recent != master_symbol || swap_counter < 0)
					{
						/*
						 * Recent is swapped symbol or swap RLE is full.
						 */
						const PB_CompressionBuffer pos = swap_counter < 0 ? PB_MAX_SWAP_RUN_LENGTH : (PB_MAX_SWAP_RUN_LENGTH - swap_counter);
						const int code_length = swap_map[recent].code_length;
						const PB_PrefixCode code = swap_map[recent].code;

						int index_entry_it;

						ENCODE_OR(code, code_length, buffer, bits_free, output_pointer);

						for (index_entry_it = 1; index_entry_it <= n_swap_index_pointers; index_entry_it++)
						{
							(index_pointer - index_entry_it)->swap_shift -= swap_counter;
						}
						n_swap_index_pointers = 0;

						if (swap_bits < 0)
						{
							*swap_pointer = *swap_pointer | (pos >> (-swap_bits));
							*(swap_pointer+1) = *(swap_pointer+1) | (pos << (swap_bits + PB_COMPRESSION_BUFFER_BIT_SIZE));
						}
						else
						{
							*swap_pointer = *swap_pointer | (pos << swap_bits);
						}

						swap_pointer = output_pointer;
						swap_bits = bits_free - PB_SWAP_RUN_LENGTH_BIT_SIZE;
						ENCODE_OR(0, PB_SWAP_RUN_LENGTH_BIT_SIZE, buffer, bits_free, output_pointer);

						swap_counter = PB_MAX_SWAP_RUN_LENGTH;
					}
				}
			}

			recent = current;
			repeated_chars = 0;
		}
	}

	/*
	 * Last symbol.
	 *
	 * This is ugly. But the alternatives are:
	 * 		a. another cmp + jmp in the loop -> cost 5% performance
	 * 		b. repeating a bunch of lines -> even uglier
	 *
	 */
	if (i == (-1)) {
		i--;
		repeated_chars++;
		PB_DEBUG3(errmsg("last symbol"));
		goto encode_pc_swp_rle_idx_out;
	}

	/*
	 * Flush buffer.
	 */
	if (bits_free < PB_COMPRESSION_BUFFER_BIT_SIZE)
		*output_pointer |= (buffer << bits_free);

	if (swap_bits < 0)
	{
		*swap_pointer = *swap_pointer | ((PB_CompressionBuffer) 0xFFFF >> (-swap_bits));
		*(swap_pointer+1) = *(swap_pointer+1) | ((PB_CompressionBuffer) 0xFFFF << (swap_bits + PB_COMPRESSION_BUFFER_BIT_SIZE));
	}
	else
	{
		*swap_pointer = *swap_pointer | ((PB_CompressionBuffer) 0xFFFF << swap_bits);
	}

	if (n_swap_index_pointers > 0)
	{
		int index_entry_it;
		for (index_entry_it = 1; index_entry_it <= n_swap_index_pointers; index_entry_it++)
		{
			(index_pointer - index_entry_it)->swap_shift -= swap_counter;
		}
	}

	PB_DEBUG1(errmsg("encode_pc_swp_rle_idx(): pos:%ld", (output_pointer - PB_COMPRESSED_SEQUENCE_STREAM_POINTER(output))*8));

	pfree((PB_EncodingMap*) master_map);
	pfree((PB_EncodingMap*) swap_map);

#ifdef DEBUG
	if (((uint8*)output_pointer) - ((uint8*) output) >= VARSIZE(output)
			&& bits_free < PB_COMPRESSION_BUFFER_BIT_SIZE) {
		ereport(ERROR,(errmsg("segmentation fault"),
				errhint("Recognized in %s at line %d.", __FILE__, __LINE__),
				errdetail("Output pointer %p is %ld bytes after a block of size %u.",
				output_pointer,
				(((uint8*)output_pointer) - ((uint8*) output))- VARSIZE(output),
				VARSIZE(output))));
	}
#endif

	PB_TRACE(errmsg("<-encode_pc_swp_rle_idx()"));
}


/**
 * decode_pc_idx()
 * 		Decode a sequence, possibly indexed
 *
 */
static void decode_pc_idx(Varlena* input,
						  uint8* output,
						  uint32 start_position,
						  uint32 output_length,
						  PB_IndexEntry* start_entry,
						  PB_CodeSet* codeset)
{
	PB_CompressionBuffer buffer;
	int bits_in_buffer;
	int i;

	PB_CompressedSequence* sequence_header;
	int stream_offset;

	Varlena* input_slice;
	PB_CompressionBuffer* input_pointer;
	uint8* output_pointer = output;

	PB_DecodingMap* map = get_decoding_map(codeset, PB_NO_SWAP_MAP);

	PB_TRACE(errmsg("->decode_pc_idx()"));

	/*
	 * Calculate the offset of the stream in the compressed data.
	 */
	sequence_header = (PB_CompressedSequence*)
					  PG_DETOAST_DATUM_SLICE(input, 0, sizeof(PB_CompressedSequence) - VARHDRSZ);
	stream_offset = PB_COMPRESSED_SEQUENCE_STREAM_OFFSET(sequence_header) - VARHDRSZ;
	pfree(sequence_header);

	PB_DEBUG3(errmsg("decode_pc_idx():calculated stream offset:%u", stream_offset));

	/*
	 * Detoast the data slice that contains the desired substring.
	 */
	if (codeset->has_equal_length)
	{
		/*
		 * All codes have equal length.
		 *  -> exact slice, that contains the desired sequence can be computed
		 */
		int code_length = codeset->words[0].code_length;
		int64 bits_to_skip = start_position * code_length;
		int64 slice_start = bits_to_skip / PB_COMPRESSION_BUFFER_BIT_SIZE;
		int64 slice_size = (bits_to_skip + output_length * code_length)
								/ PB_COMPRESSION_BUFFER_BIT_SIZE + 1 - slice_start;

		slice_start *= PB_COMPRESSION_BUFFER_BYTE_SIZE;
		slice_start += stream_offset;
		slice_size *= PB_COMPRESSION_BUFFER_BYTE_SIZE;
		input_slice = (Varlena*)
					  PG_DETOAST_DATUM_SLICE(input,
											 slice_start,
											 slice_size);

		PB_DEBUG1(errmsg("decode_pc_idx(): all codes have equal length\n\tskipping %ld bits\n\tslice starts at byte %ld\n\tslice size is %ld bytes", bits_to_skip, slice_start,slice_size));

		i = -1;
		input_pointer = (PB_CompressionBuffer*) VARDATA_ANY(input_slice);
		bits_in_buffer = bits_to_skip % PB_COMPRESSION_BUFFER_BIT_SIZE;
		buffer = *input_pointer << bits_in_buffer;
		bits_in_buffer = PB_COMPRESSION_BUFFER_BIT_SIZE - bits_in_buffer;
		input_pointer++;
	}
	else
	{
		/*
		 * Cannot compute the exact slice.
		 */
		int max_codeword_length = codeset->words[codeset->n_symbols - 1].code_length;
		int raw_size = toast_raw_datum_size((Datum)input);

		if (!start_entry)
		{
			/*
			 * Compute upper bound for slice size.
			 */
			int slice_size = (start_position + output_length) *
							 max_codeword_length /
							 PB_COMPRESSION_BUFFER_BIT_SIZE + 1;
			slice_size *= PB_COMPRESSION_BUFFER_BYTE_SIZE;

			/*
			 * If upper bound exceeds raw size, limit slice size.
			 */
			if (slice_size + stream_offset > raw_size)
				slice_size = raw_size - stream_offset;

			input_slice = (Varlena*)
						  PG_DETOAST_DATUM_SLICE(input,
												 stream_offset,
												 slice_size);

			PB_DEBUG1(errmsg("decode_pc_idx(): skipping through sequence\n\tslice size is %d bytes", slice_size));

			bits_in_buffer = 0;
			buffer = 0;
			i = start_position - 1;
			input_pointer = (PB_CompressionBuffer*) VARDATA_ANY(input_slice);
		}
		else
		{
			/*
			 * An index entry is defined.
			 *  -> get the slice starting at indexed position.
			 */
			int slice_start = stream_offset +
							  start_entry->block * PB_COMPRESSION_BUFFER_BYTE_SIZE;
			int slice_size = (output_length + (start_position % PB_INDEX_PART_SIZE)) *
							 max_codeword_length /
							 PB_COMPRESSION_BUFFER_BIT_SIZE + 1;
			slice_size *= PB_COMPRESSION_BUFFER_BYTE_SIZE;

			/*
			 * If upper bound exceeds raw size, limit slice size.
			 */
			if (slice_size + slice_start > raw_size)
				slice_size = raw_size - slice_start;

			input_slice = (Varlena*)
						  PG_DETOAST_DATUM_SLICE(input,
												 slice_start,
												 slice_size);

			PB_DEBUG1(errmsg("decode_pc_idx(): index entry given\n\tstarting in block %u\n\tslice starts at byte %d\n\tslice size is %d bytes", start_entry->block, slice_start, slice_size));

			input_pointer = (PB_CompressionBuffer*) VARDATA_ANY(input_slice);
			bits_in_buffer = PB_COMPRESSION_BUFFER_BIT_SIZE - start_entry->bit;
			buffer = *(input_pointer) << start_entry->bit;
			input_pointer++;
			i = ((start_position + 1) % PB_INDEX_PART_SIZE) - 1;
		}

		PB_DEBUG1(errmsg("decode_pc_idx(): reading %d chars to skip", i + 1));

		while (i >= 0)
		{
			PB_PrefixCode val;
			int length;

			DECODE(input_pointer, buffer, bits_in_buffer, val, length, map);
			i--;
		}
	}

	i += output_length;

	PB_DEBUG1(errmsg("decode_pc_idx(): decoding %d chars",i));

	while (i >= 0)
	{
		PB_PrefixCode val;
		int length;
		DECODE(input_pointer, buffer, bits_in_buffer, val, length, map);

		PB_DEBUG3(errmsg("decode_pc_idx(): i=%d val=%u len=%u s=%c buf=%08X%08X bib=%d", i, val, length, map[val].symbol, (unsigned int) (buffer >> 32), (unsigned int) buffer, bits_in_buffer));

		*output_pointer = map[val].symbol;
		output_pointer++;
		i--;
	}

	pfree(map);
	pfree(input_slice);

	PB_TRACE(errmsg("<-decode_pc_idx()"));
}

/**
 * decode_pc_rle_idx()
 * 		Decode a rle sequence, possibly indexed
 *
 */
static void decode_pc_rle_idx(Varlena* input,
							  uint8* output,
							  uint32 start_position,
							  uint32 output_length,
							  PB_IndexEntry* start_entry,
							  PB_CodeSet* codeset)
{
	PB_CompressionBuffer buffer;
	int bits_in_buffer;
	int i;

	PB_CompressedSequence* sequence_header;
	int stream_offset;

	Varlena* input_slice;
	PB_CompressionBuffer* input_pointer;
	uint8* output_pointer = output;

	PB_DecodingMap* map = get_decoding_map(codeset, PB_NO_SWAP_MAP);

	const int max_codeword_length = codeset->words[codeset->n_symbols - 1].code_length;
	const int raw_size = toast_raw_datum_size((Datum)input);

	PB_TRACE(errmsg("->decode_pc_rle_idx()"));

	/*
	 * Calculate the offset of the stream in the compressed data.
	 */
	sequence_header = (PB_CompressedSequence*)
					  PG_DETOAST_DATUM_SLICE(input, 0, sizeof(PB_CompressedSequence) - VARHDRSZ);
	stream_offset = PB_COMPRESSED_SEQUENCE_STREAM_OFFSET(sequence_header) - VARHDRSZ;
	pfree(sequence_header);

	PB_DEBUG1(errmsg("decode_pc_rle_idx():calculated stream offset:%u", stream_offset));

	/*
	 * Detoast the data slice that contains the desired substring.
	 */
	if (!start_entry)
	{
		/*
		 * Calculate upper bound for slice size.
		 */
		int slice_size = ((start_position + output_length + 1) *
						 max_codeword_length + 8) /
						 PB_COMPRESSION_BUFFER_BIT_SIZE + 1;
		slice_size *= PB_COMPRESSION_BUFFER_BYTE_SIZE;

		/*
		 * If upper bound exceeds raw size, limit slice size.
		 */
		if (slice_size + stream_offset > raw_size)
			slice_size = raw_size - stream_offset;

		input_slice = (Varlena*)
					  PG_DETOAST_DATUM_SLICE(input,
											 stream_offset,
											 slice_size);

		PB_DEBUG1(errmsg("decode_pc_rle_idx(): skipping through sequence\n\tslice size is %d bytes", slice_size));

		bits_in_buffer = 0;
		buffer = 0;
		i = start_position - 1;
		input_pointer = (PB_CompressionBuffer*) VARDATA_ANY(input_slice);
	}
	else
	{
		/*
		 * An index entry is defined.
		 *  -> get the slice starting at indexed position.
		 */
		int slice_start = stream_offset +
						  start_entry->block * PB_COMPRESSION_BUFFER_BYTE_SIZE;
		int slice_size = ((output_length + (start_position % PB_INDEX_PART_SIZE) + 1) *
						 max_codeword_length  + 8) /
						 PB_COMPRESSION_BUFFER_BIT_SIZE + 1;
		slice_size *= PB_COMPRESSION_BUFFER_BYTE_SIZE;

		/*
		 * If upper bound exceeds raw size, limit slice size.
		 */
		if (slice_size + slice_start > raw_size)
			slice_size = raw_size - slice_start;

		input_slice = (Varlena*)
					  PG_DETOAST_DATUM_SLICE(input,
											 slice_start,
											 slice_size);

		PB_DEBUG1(errmsg("decode_pc_rle_idx(): index entry given\n\tstarting in block %u\n\tslice starts at byte %d\n\tslice size is %d bytes", start_entry->block, slice_start, slice_size));

		input_pointer = (PB_CompressionBuffer*) VARDATA_ANY(input_slice);
		bits_in_buffer = PB_COMPRESSION_BUFFER_BIT_SIZE - start_entry->bit;
		buffer = *(input_pointer) << start_entry->bit;
		input_pointer++;
		i = ((start_position + 1) % PB_INDEX_PART_SIZE) - 1 + start_entry->rle_shift;
	}

	PB_DEBUG1(errmsg("decode_pc_rle_idx(): reading %d chars to skip", i + 1));

	while (i >= 0)
	{
		PB_PrefixCode val;
		int length;
		DECODE(input_pointer, buffer, bits_in_buffer, val, length, map);

		i--;
		if (map[val].symbol == PB_RUN_LENGTH_SYMBOL)
		{
			int repeated_chars = 0;
			READ_N_BITS(input_pointer, buffer, bits_in_buffer, repeated_chars, PB_RUN_LENGTH_BIT_SIZE);

			i -= repeated_chars + PB_MIN_RUN_LENGTH - 2;
			/*
			 * Note: (-2) because 1 has been deducted in this iteration above and
			 * 			1 will be deducted in the next iteration, when the symbol
			 * 			after the RLE is read.
			 */
		}
	}

	/*
	 * Write out last run-length if there was one.
	 */
	if (i < -1)
	{
		int repeated_chars = -i;
		PB_PrefixCode val;
		int length;
		uint8 current;

		if (repeated_chars >= output_length)
		{
			repeated_chars = output_length;
		}

		PB_DEBUG1(errmsg("decode_pc_rle_idx(): %d remaining chars from last run-length", i + 1));

		DECODE(input_pointer, buffer, bits_in_buffer, val, length, map);
		i--;

		current = map[val].symbol;

		memset(output_pointer, current, repeated_chars);
		output_pointer += repeated_chars;
	}

	i += output_length;

	PB_DEBUG1(errmsg("decode_pc_rle_idx(): decoding %d chars", i + 1));

	while (i >= 0)
	{
		PB_PrefixCode val;
		int length;
		DECODE(input_pointer, buffer, bits_in_buffer, val, length, map);

		PB_DEBUG3(errmsg("decode_pc_rle_idx(): i=%d val=%u len=%u s=%c buf=%08X%08X bib=%d", i, val, length, map[val].symbol, (unsigned int) (buffer >> 32), (unsigned int) buffer, bits_in_buffer));

		if (map[val].symbol != PB_RUN_LENGTH_SYMBOL)
		{
			*output_pointer = map[val].symbol;
			output_pointer++;
			i--;
		}
		else
		{
			PB_CompressionBuffer repeated_chars = 0;
			uint8 out;

			READ_N_BITS(input_pointer, buffer, bits_in_buffer, repeated_chars, PB_RUN_LENGTH_BIT_SIZE);

			DECODE(input_pointer, buffer, bits_in_buffer, val, length, map);

			repeated_chars += PB_MIN_RUN_LENGTH;
			out = map[val].symbol;

			PB_DEBUG3(errmsg("decode_pc_rle_idx(): write %lu times %c", repeated_chars, out));

			i -= (int) repeated_chars;
			if (i < (-1))
			{
				repeated_chars -= (PB_CompressionBuffer) (-(i + 1));
			}

			memset(output_pointer, out, repeated_chars);
			output_pointer += repeated_chars;
		}
	}

	pfree(map);
	pfree(input_slice);

	PB_TRACE(errmsg("<-decode_pc_rle_idx()"));
}

/**
 * decode_pc_swp_idx()
 * 		Decode a sequence from encode_pc_swp or encode_pc_swp_idx
 *
 */
static void decode_pc_swp_idx(Varlena* input,
							  uint8* output,
							  uint32 start_position,
							  uint32 output_length,
							  PB_IndexEntry* start_entry,
							  PB_CodeSet* codeset)
{
	PB_CompressionBuffer buffer;
	int bits_in_buffer;
	int i;
	int swap_counter;

	PB_CompressedSequence* sequence_header;
	int stream_offset;

	Varlena* input_slice;
	PB_CompressionBuffer* input_pointer;
	uint8* output_pointer = output;

	PB_DecodingMap* map = get_decoding_map(codeset, PB_NO_SWAP_MAP);
	PB_DecodingMap* swap_map = get_decoding_map(codeset, PB_SWAP_MAP);

	const uint8 master_symbol = codeset->words[codeset->n_symbols - codeset->n_swapped_symbols].symbol;

	const int max_codeword_length = codeset->words[codeset->n_symbols - 1].code_length +
									codeset->words[master_symbol].code_length +
									PB_SWAP_RUN_LENGTH_BIT_SIZE;
	const int raw_size = toast_raw_datum_size((Datum)input);

	PB_TRACE(errmsg("->decode_pc_swp_idx()"));

	/*
	 * Calculate the offset of the stream in the compressed data.
	 */
	sequence_header = (PB_CompressedSequence*)
					  PG_DETOAST_DATUM_SLICE(input, 0, sizeof(PB_CompressedSequence) - VARHDRSZ);
	stream_offset = PB_COMPRESSED_SEQUENCE_STREAM_OFFSET(sequence_header) - VARHDRSZ;
	pfree(sequence_header);

	PB_DEBUG1(errmsg("decode_pc_swp_idx():calculated stream offset:%u", stream_offset));

	/*
	 * Detoast the data slice that contains the desired substring.
	 */
	if (!start_entry)
	{
		/*
		 * Calculate upper bound for slice size.
		 */
		int slice_size = (start_position + output_length) *
						 max_codeword_length /
						 PB_COMPRESSION_BUFFER_BIT_SIZE + 1;
		slice_size *= PB_COMPRESSION_BUFFER_BYTE_SIZE;

		/*
		 * If upper bound exceeds raw size, limit slice size.
		 */
		if (slice_size + stream_offset > raw_size)
			slice_size = raw_size - stream_offset;

		input_slice = (Varlena*)
					  PG_DETOAST_DATUM_SLICE(input,
											 stream_offset,
											 slice_size);

		PB_DEBUG1(errmsg("decode_pc_swp_idx(): skipping through sequence\n\tslice size is %d bytes", slice_size));

		i = start_position - 1;
		input_pointer = (PB_CompressionBuffer*) VARDATA_ANY(input_slice);

		buffer = *input_pointer;
		input_pointer++;
		swap_counter = buffer >> (PB_COMPRESSION_BUFFER_BIT_SIZE - PB_SWAP_RUN_LENGTH_BIT_SIZE);
		bits_in_buffer = PB_COMPRESSION_BUFFER_BIT_SIZE - PB_SWAP_RUN_LENGTH_BIT_SIZE;
		buffer = buffer << PB_SWAP_RUN_LENGTH_BIT_SIZE;
	}
	else
	{
		/*
		 * An index entry is defined.
		 *  -> get the slice starting at indexed position.
		 */
		int slice_start = stream_offset +
						  start_entry->block * PB_COMPRESSION_BUFFER_BYTE_SIZE;
		int slice_size = (output_length + (start_position % PB_INDEX_PART_SIZE)) *
						 max_codeword_length /
						 PB_COMPRESSION_BUFFER_BIT_SIZE + 1;
		slice_size *= PB_COMPRESSION_BUFFER_BYTE_SIZE;

		/*
		 * If upper bound exceeds raw size, limit slice size.
		 */
		if (slice_size + slice_start > raw_size)
			slice_size = raw_size - slice_start;

		input_slice = (Varlena*)
					  PG_DETOAST_DATUM_SLICE(input,
											 slice_start,
											 slice_size);

		PB_DEBUG1(errmsg("decode_pc_swp_idx(): index entry given\n\tstarting in block %u\n\tslice starts at byte %d\n\tslice size is %d bytes\n\tswap shift:%d", start_entry->block, slice_start, slice_size,start_entry->swap_shift));

		input_pointer = (PB_CompressionBuffer*) VARDATA_ANY(input_slice);
		bits_in_buffer = PB_COMPRESSION_BUFFER_BIT_SIZE - start_entry->bit;
		buffer = *(input_pointer) << start_entry->bit;
		input_pointer++;
		i = ((start_position + 1) % PB_INDEX_PART_SIZE) - 1;
		swap_counter = start_entry->swap_shift;
	}

	PB_DEBUG1(errmsg("decode_pc_swp_idx(): reading %d chars to skip, swap_counter = %d, bib=%d", i + 1, swap_counter, bits_in_buffer));

	/*
	 * Skipping start of sequence.
	 */
	while (i >= 0)
	{
		PB_PrefixCode val;
		int length;

		DECODE(input_pointer, buffer, bits_in_buffer, val, length, map);
		i--;
		if (map[val].symbol == master_symbol)
		{
			swap_counter--;
			if (swap_counter < 0)
			{
				DECODE(input_pointer, buffer, bits_in_buffer, val, length, swap_map);
				READ_N_BITS(input_pointer, buffer, bits_in_buffer, swap_counter, PB_SWAP_RUN_LENGTH_BIT_SIZE);
			}
		}
	}

	i += output_length;

	PB_DEBUG1(errmsg("decode_pc_swp_idx(): reading %d chars, swap_counter = %d, bib=%d", i + 1, swap_counter, bits_in_buffer));

	while (i >= 0)
	{
		PB_PrefixCode val;
		int length;

		DECODE(input_pointer, buffer, bits_in_buffer, val, length, map);
		if (map[val].symbol != master_symbol)
		{
			PB_DEBUG3(errmsg("NORMAL: i:%d val:%u len:%u s:%c buf:%08X%08X bib:%d", i, val, length, map[val].symbol, (unsigned int) (buffer >> 32), (unsigned int) buffer, bits_in_buffer));

			*output_pointer = map[val].symbol;
			output_pointer++;
			i--;
		}
		else
		{
			if (swap_counter > 0)
			{
				PB_DEBUG3(errmsg("MASTER: i:%d swap_counter:%d val:%u len:%u s:%c buf:%08X%08X bib:%d", i, swap_counter, val, length, master_symbol, (unsigned int) (buffer >> 32), (unsigned int) buffer, bits_in_buffer));

				*output_pointer = master_symbol;
				output_pointer++;
				i--;
				swap_counter--;
			}
			else
			{
				DECODE(input_pointer, buffer, bits_in_buffer, val, length, swap_map);
				PB_DEBUG3(errmsg("SLAVE i:%d val:%u len:%u s:%c buf:%08X%08X bib:%d", i, val, length, swap_map[val].symbol, (unsigned int) (buffer >> 32), (unsigned int) buffer, bits_in_buffer));

				*output_pointer = swap_map[val].symbol;
				output_pointer++;
				i--;
				READ_N_BITS(input_pointer, buffer, bits_in_buffer, swap_counter, PB_SWAP_RUN_LENGTH_BIT_SIZE);
			}
		}
	}

	pfree(map);
	pfree(swap_map);
	pfree(input_slice);

	PB_TRACE(errmsg("<-decode_pc_swp_idx()"));
}

/**
 * decode_pc_swp_rle_idx()
 * 		Decode a sequence from encode_pc_swp_rle or encode_pc_swp_rle_idx
 *
 */
static void decode_pc_swp_rle_idx(Varlena* input,
								  uint8* output,
								  uint32 start_position,
								  uint32 output_length,
								  PB_IndexEntry* start_entry,
								  PB_CodeSet* codeset)
{
	PB_CompressionBuffer buffer;
	int bits_in_buffer;
	int i;
	int swap_counter;

	PB_CompressedSequence* sequence_header;
	int stream_offset;

	Varlena* input_slice;
	PB_CompressionBuffer* input_pointer;
	uint8* output_pointer = output;

	PB_DecodingMap* map = get_decoding_map(codeset, PB_NO_SWAP_MAP);
	PB_DecodingMap* swap_map = get_decoding_map(codeset, PB_SWAP_MAP);

	const uint8 master_symbol = codeset->words[codeset->n_symbols - codeset->n_swapped_symbols].symbol;
	const int max_codeword_length = codeset->words[codeset->n_symbols - 1].code_length +
									codeset->words[master_symbol].code_length +
									PB_SWAP_RUN_LENGTH_BIT_SIZE;
	const int raw_size = toast_raw_datum_size((Datum)input);

	PB_TRACE(errmsg("->decode_pc_swp_rle_idx(%u,%u)", start_position, output_length));

	/*
	 * Calculate the offset of the stream in the compressed data.
	 */
	sequence_header = (PB_CompressedSequence*)
					  PG_DETOAST_DATUM_SLICE(input, 0, sizeof(PB_CompressedSequence) - VARHDRSZ);
	stream_offset = PB_COMPRESSED_SEQUENCE_STREAM_OFFSET(sequence_header) - VARHDRSZ;
	pfree(sequence_header);

	PB_DEBUG1(errmsg("decode_pc_swp_rle_idx():calculated stream offset:%u", stream_offset));

	/*
	 * Detoast the data slice that contains the desired substring.
	 */
	if (!start_entry)
	{
		/*
		 * Calculate upper bound for slice size.
		 */
		int slice_size = (start_position + output_length) *
						 max_codeword_length /
						 PB_COMPRESSION_BUFFER_BIT_SIZE + 1;
		slice_size *= PB_COMPRESSION_BUFFER_BYTE_SIZE;

		/*
		 * If upper bound exceeds raw size, limit slice size.
		 */
		if (slice_size + stream_offset > raw_size)
			slice_size = raw_size - stream_offset;

		input_slice = (Varlena*)
					  PG_DETOAST_DATUM_SLICE(input,
											 stream_offset,
											 slice_size);

		PB_DEBUG1(errmsg("decode_pc_swp_rle_idx(): skipping through sequence\n\tslice size is %d bytes", slice_size));

		i = start_position - 1;
		input_pointer = (PB_CompressionBuffer*) VARDATA_ANY(input_slice);

		buffer = *input_pointer;
		input_pointer++;
		swap_counter = buffer >> (PB_COMPRESSION_BUFFER_BIT_SIZE - PB_SWAP_RUN_LENGTH_BIT_SIZE);
		bits_in_buffer = PB_COMPRESSION_BUFFER_BIT_SIZE - PB_SWAP_RUN_LENGTH_BIT_SIZE;
		buffer = buffer << PB_SWAP_RUN_LENGTH_BIT_SIZE;

		PB_DEBUG1(errmsg("buf:%08X%08X bib:%d", (unsigned int) (buffer >> 32), (unsigned int) buffer, bits_in_buffer));
	}
	else
	{
		/*
		 * An index entry is defined.
		 *  -> get the slice starting at indexed position.
		 */
		int slice_start = stream_offset +
						  start_entry->block * PB_COMPRESSION_BUFFER_BYTE_SIZE;
		int slice_size = (output_length + (start_position % PB_INDEX_PART_SIZE)) *
						 max_codeword_length /
						 PB_COMPRESSION_BUFFER_BIT_SIZE + 1;
		slice_size *= PB_COMPRESSION_BUFFER_BYTE_SIZE;

		/*
		 * If upper bound exceeds raw size, limit slice size.
		 */
		if (slice_size + slice_start > raw_size)
			slice_size = raw_size - slice_start;

		input_slice = (Varlena*)
					  PG_DETOAST_DATUM_SLICE(input,
											 slice_start,
											 slice_size);

		PB_DEBUG1(errmsg("decode_pc_swp_rle_idx(): index entry given\n\tstarting in block %u\n\tslice starts at byte %d\n\tslice size is %d bytes\nrle_shift:%u\nswap_shift:%u", start_entry->block, slice_start, slice_size,start_entry->rle_shift, start_entry->swap_shift));

		input_pointer = (PB_CompressionBuffer*) VARDATA_ANY(input_slice);
		bits_in_buffer = PB_COMPRESSION_BUFFER_BIT_SIZE - start_entry->bit;
		buffer = *(input_pointer) << start_entry->bit;
		input_pointer++;
		i = ((start_position + 1) % PB_INDEX_PART_SIZE) - 1 + start_entry->rle_shift;
		swap_counter = start_entry->swap_shift;
	}

	PB_DEBUG1(errmsg("decode_pc_swp_rle_idx(): reading %d chars to skip, swap_counter = %d, bib=%d, rle_shift=%u", i + 1, swap_counter, bits_in_buffer, start_entry == NULL ? -1 : start_entry->rle_shift));

	/*
	 * Skipping.
	 */
	while (i >= 0)
	{
		PB_PrefixCode val;
		int length;
		uint8 current;

		DECODE(input_pointer, buffer, bits_in_buffer, val, length, map);
		current = map[val].symbol;
		i--;

		if (current == master_symbol)
		{
			swap_counter--;
			if (swap_counter < 0)
			{
				DECODE(input_pointer, buffer, bits_in_buffer, val, length, swap_map);
				current = swap_map[val].symbol;
				READ_N_BITS(input_pointer, buffer, bits_in_buffer, swap_counter, PB_SWAP_RUN_LENGTH_BIT_SIZE);
			}
		}

		if (current == PB_RUN_LENGTH_SYMBOL)
		{
			int repeated_chars = 0;

			READ_N_BITS(input_pointer, buffer, bits_in_buffer, repeated_chars, PB_RUN_LENGTH_BIT_SIZE);

			i -= repeated_chars + PB_MIN_RUN_LENGTH - 2;
			/*
			 * Note: (-2) because 1 has been deducted in this iteration above and
			 * 			1 will be deducted in the next iteration, when the symbol
			 * 			after the RLE is read.
			 */
		}
	}

	if (i < -1)
	{
		/*
		 * The symbol read was an RLE symbol. Write out remaining chars.
		 */
		int repeated_chars = -i;
		PB_PrefixCode val;
		int length;
		uint8 current;

		if (repeated_chars >= output_length)
		{
			repeated_chars = output_length;
		}

		PB_DEBUG1(errmsg("decode_pc_swp_rle_idx(): %d remaining chars from last rle", repeated_chars));

		DECODE(input_pointer, buffer, bits_in_buffer, val, length, map);
		i--;
		current = map[val].symbol;
		if (current == master_symbol)
		{
			swap_counter--;
			if (swap_counter < 0)
			{
				PB_PrefixCode val;
				int length;

				DECODE(input_pointer, buffer, bits_in_buffer, val, length, swap_map);
				current = swap_map[val].symbol;

				READ_N_BITS(input_pointer, buffer, bits_in_buffer, swap_counter, PB_SWAP_RUN_LENGTH_BIT_SIZE);
			}
		}

		PB_DEBUG1(errmsg("decode_pc_swp_rle_idx(): write %dx%c", repeated_chars,current));

		memset(output_pointer, current, repeated_chars);
		output_pointer += repeated_chars;
	}

	i += output_length;

	PB_DEBUG1(errmsg("decode_pc_swp_rle_idx(): reading %d chars, swap_counter = %d, bib=%d", i + 1, swap_counter, bits_in_buffer));

	while (i >= 0)
	{
		PB_PrefixCode val;
		int length;
		uint8 current;

		DECODE(input_pointer, buffer, bits_in_buffer, val, length, map);
		current = map[val].symbol;

		if (current != master_symbol && current != PB_RUN_LENGTH_SYMBOL)
		{
			*output_pointer = current;
			output_pointer++;
			i--;

			PB_DEBUG3(errmsg("NSNRLE %c (%d) code:%u len:%u i:%d swap_counter:%d bitsInBuffer:%d buffer:%08X%08X", current, current,val, length,  i, swap_counter, bits_in_buffer, (uint32) (buffer >> 32), (uint32) buffer));
		}
		else
		{
			if (current == master_symbol)
			{
				swap_counter--;
				if (swap_counter < 0)
				{
					PB_PrefixCode val;
					int length;

					DECODE(input_pointer, buffer, bits_in_buffer, val, length, swap_map);
					current = swap_map[val].symbol;

					READ_N_BITS(input_pointer, buffer, bits_in_buffer, swap_counter, PB_SWAP_RUN_LENGTH_BIT_SIZE);
				}
			}

			if (current != PB_RUN_LENGTH_SYMBOL)
			{
				*output_pointer = current;
				output_pointer++;
				i--;

				PB_DEBUG3(errmsg("SNRLE %c (%d) code:%u len:%u i:%d swap_counter:%d bitsInBuffer:%d buffer:%08X%08X", current, current,val, length,  i, swap_counter, bits_in_buffer, (uint32) (buffer >> 32), (uint32) buffer));
			}
			else
			{
				PB_CompressionBuffer repeated_chars = 0;

				READ_N_BITS(input_pointer, buffer, bits_in_buffer, repeated_chars, PB_RUN_LENGTH_BIT_SIZE);

				DECODE(input_pointer, buffer, bits_in_buffer, val, length, map);
				repeated_chars += PB_MIN_RUN_LENGTH;
				current = map[val].symbol;

				if (current == master_symbol)
				{
					swap_counter--;
					if (swap_counter < 0)
					{
						PB_PrefixCode val;
						int length;

						DECODE(input_pointer, buffer, bits_in_buffer, val, length, swap_map);
						current = swap_map[val].symbol;

						READ_N_BITS(input_pointer, buffer, bits_in_buffer, swap_counter, PB_SWAP_RUN_LENGTH_BIT_SIZE);
					}
				}

				PB_DEBUG3(errmsg("SRLE %c (%d) code:%u len:%u i:%d swap_counter:%lu bitsInBuffer:%d buffer:%08X%08X", current, current,val, length,  i, repeated_chars, bits_in_buffer, (uint32) (buffer >> 32), (uint32) buffer));

				i -= (int) repeated_chars;
				if (i < (-1))
				{
					repeated_chars -= (PB_CompressionBuffer) (-(i + 1));
				}

				memset(output_pointer, current, repeated_chars);
				output_pointer += repeated_chars;
			}
		}
	}

	pfree((PB_DecodingMap*) map);
	pfree((PB_DecodingMap*) swap_map);

	PB_TRACE(errmsg("<-decode_pc_swp_rle_idx()"))
}

/*
 * public functions
 */

/**
 * get_compressed_size()
 * 		Compute the size of a compressed sequence.
 *
 * 	PB_SequenceInfo* info : stats of input sequence
 * 	PB_CodeSet* codeset : code to test
 */
uint32 get_compressed_size(const PB_SequenceInfo* info,
						   PB_CodeSet* codeset)
{
	uint64 total_stream_size_bits = 0;
	uint32 total_size = 0;

	int i;

	const uint32* frequencies;

	PB_TRACE(errmsg("->get_compressed_stream_size()"));

	if (info->n_symbols > codeset->n_symbols)
		ereport(ERROR,(errmsg("code set can not express given sequence")));

	if (codeset->uses_rle)
	{
		frequencies = info->rle_info->rle_frequencies;
		total_stream_size_bits += frequencies[PB_RUN_LENGTH_SYMBOL] * PB_RUN_LENGTH_BIT_SIZE;
	}
	else
		frequencies = info->frequencies;

	if (codeset->n_swapped_symbols > 0 && codeset->is_fixed == FALSE)
	{
		/*
		 * If code set is swapped
		 */
		uint8 master_symbol = codeset->words[codeset->n_symbols - codeset->n_swapped_symbols].symbol;
		uint8 master_symbol_length = 0;

		for (i = 0; i < codeset->n_symbols - codeset->n_swapped_symbols; i++)
		{
			total_stream_size_bits += frequencies[codeset->words[i].symbol] * codeset->words[i].code_length;
			if (codeset->words[i].symbol == master_symbol)
			{
				master_symbol_length = codeset->words[i].code_length + PB_SWAP_RUN_LENGTH_BIT_SIZE;
			}

			PB_DEBUG1(errmsg("M %c: %u x %u = %u", codeset->words[i].symbol,frequencies[codeset->words[i].symbol], codeset->words[i].code_length, frequencies[codeset->words[i].symbol] * codeset->words[i].code_length));
		}

		for (i = codeset->n_symbols - codeset->n_swapped_symbols + 1; i < codeset->n_symbols; i++)
		{
			total_stream_size_bits += frequencies[codeset->words[i].symbol] *
					(codeset->words[i].code_length + master_symbol_length);

			PB_DEBUG1(errmsg("S %c: %u x %u = %u", codeset->words[i].symbol,frequencies[codeset->words[i].symbol], codeset->words[i].code_length + master_symbol_length, frequencies[codeset->words[i].symbol] * (codeset->words[i].code_length + master_symbol_length)));
		}

		total_stream_size_bits += PB_SWAP_RUN_LENGTH_BIT_SIZE;
		total_stream_size_bits += (frequencies[master_symbol] / PB_MAX_SWAP_RUN_LENGTH) * (PB_SWAP_RUN_LENGTH_BIT_SIZE + 1);
	}
	else
	{
		for (i = 0; i < codeset->n_symbols; i++)
		{
			total_stream_size_bits += frequencies[codeset->words[i].symbol] * codeset->words[i].code_length;

			PB_DEBUG1(errmsg("P %c: %u x %u = %u", codeset->words[i].symbol,frequencies[codeset->words[i].symbol], codeset->words[i].code_length, frequencies[codeset->words[i].symbol] * codeset->words[i].code_length));
		}
	}

	total_stream_size_bits = PB_ALIGN_BIT_SIZE(total_stream_size_bits);

	total_size = sizeof(PB_CompressedSequence);
	total_size += codeset->is_fixed ? 0 : sizeof(PB_Codeword) * codeset->n_symbols;
	total_size += codeset->has_equal_length ? 0 : info->sequence_length / PB_INDEX_PART_SIZE * sizeof(PB_IndexEntry);
	total_size = PB_ALIGN_BYTE_SIZE(total_size);
	total_size +=  total_stream_size_bits / 8;

	PB_DEBUG1(errmsg("get_compressed_stream_size(): totalsize in byte:%u (%lu bits)", total_size, total_stream_size_bits));

	PB_TRACE(errmsg("<-get_compressed_stream_size()"));

	return total_size;
}


/**
 * encode()
 * 		Encode a sequence.
 *
 * 	uint8* input : input sequence, not null-terminated
 * 					length must be in info
 * 	uint32 compressed_size : length of the compressed stream, calculated
 * 									with get_compressed_stream_size()
 * 	PB_CodeSet* codeset : codeset for encoding
 * 	PB_SequenceInfo* info : info about the sequence to compresss
 */
PB_CompressedSequence* encode(uint8* input,
							  uint32 compressed_size,
							  PB_CodeSet* codeset,
							  PB_SequenceInfo* info)
{
	PB_CompressedSequence* result;

	PB_TRACE(errmsg("->encode(): compressed size:%u, uncompressed size:%u", compressed_size, info->sequence_length));

	/*
	 * Initialization of result
	 */
	result = palloc0(compressed_size);
	SET_VARSIZE(result, compressed_size);
	result->sequence_length = info->sequence_length;
	if (codeset->is_fixed)
	{
		result->is_fixed = TRUE;
		result->n_symbols = 0;
		result->n_swapped_symbols = codeset->fixed_id;
		PB_DEBUG1(errmsg("encode(): uses fix code with id %d", codeset->fixed_id));
	}
	else
	{
		PB_Codeword* code = PB_COMPRESSED_SEQUENCE_SYMBOL_POINTER(result);

		result->is_fixed = FALSE;
		result->n_symbols = codeset->n_symbols;
		result->n_swapped_symbols = codeset->n_swapped_symbols;
		memcpy(code,
			   codeset->words,
			   codeset->n_symbols * sizeof(PB_Codeword));
		PB_DEBUG1(errmsg("encode(): copied sequence specific code"));
	}

	result->has_equal_length = codeset->has_equal_length;
	result->uses_rle = codeset->uses_rle;

	if (codeset->has_equal_length || info->sequence_length < PB_INDEX_PART_SIZE)
		result->has_index = FALSE;
	else
		result->has_index = TRUE;

	/*
	 * Choose the encoding function
	 */
	if (result->has_index)
	{
		if (codeset->n_swapped_symbols > 0)
		{
			if (codeset->uses_rle)
				encode_pc_swp_rle_idx(input, result, codeset);
			else
				encode_pc_swp_idx(input, result, codeset);
		}
		else
		{
			if (codeset->uses_rle)
				encode_pc_rle_idx(input, result, codeset);
			else
				encode_pc_idx(input, result, codeset);
		}
	}
	else
	{
		if (codeset->n_swapped_symbols > 0)
		{
			if (codeset->uses_rle)
				encode_pc_swp_rle(input, result, codeset);
			else
				encode_pc_swp(input, result, codeset);
		}
		else
		{
			if (codeset->uses_rle)
				encode_pc_rle(input, result, codeset);
			else
				encode_pc(input, result, codeset);
		}
	}

	PB_TRACE(errmsg("<-encode()"));

	return result;
}

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
		PB_CodeSet** fixed_codesets)
{
	PB_CompressedSequence* input_header;
	PB_CodeSet* codeset;
	PB_IndexEntry* start_entry = NULL;

	PB_TRACE(errmsg("->decode()"));

	/*
	 * detoast only the header of the sequence
	 */
	input_header = (PB_CompressedSequence*) PG_DETOAST_DATUM_SLICE(input,
																   0,
																   sizeof(PB_CompressedSequence) - VARHDRSZ);

	PB_DEBUG1(errmsg("decode(): input header detoasted\n\tsequence_length:%u\n\tn_symbols:%u\n\tn_swapped_symbols:%u\n\thas_equal_length:%d\n\thas_index:%d\n\tis_fixed:%d\n\tuses_rle:%d",
							input_header->sequence_length, input_header->n_symbols, input_header->n_swapped_symbols, input_header->has_equal_length, input_header->has_index, input_header->is_fixed, input_header->uses_rle));


	/*
	 * Restore codeset.
	 */
	if (input_header->is_fixed)
	{
		codeset = fixed_codesets[input_header->n_swapped_symbols];

		PB_DEBUG1(errmsg("decode():uses fixed code with id %u", input_header->n_swapped_symbols));
	}
	else
	{
		int code_size = sizeof(PB_Codeword) * input_header->n_symbols;
		PB_Codeword* code;
		int i;

		codeset = palloc0(sizeof(PB_CodeSet) + code_size);
		codeset->n_symbols = input_header->n_symbols;
		codeset->n_swapped_symbols = input_header->n_swapped_symbols;
		codeset->is_fixed = FALSE;
		codeset->has_equal_length = input_header->has_equal_length;
		codeset->uses_rle = input_header->uses_rle;

		pfree(input_header);
		input_header = (PB_CompressedSequence*)
					   PG_DETOAST_DATUM_SLICE(input,0,sizeof(PB_CompressedSequence) - VARHDRSZ +
													  code_size);

		code = PB_COMPRESSED_SEQUENCE_SYMBOL_POINTER(input_header);
		memcpy(codeset->words, code, code_size);

		for (i = 0; i < codeset->n_symbols - codeset->n_swapped_symbols; i++)
			if (codeset->max_codeword_length < codeset->words[i].code_length)
				codeset->max_codeword_length = codeset->words[i].code_length;

		for (i = codeset->n_symbols - codeset->n_swapped_symbols; i < codeset->n_symbols; i++)
			if (codeset->max_codeword_length < codeset->words[i].code_length)
				codeset->max_codeword_length = codeset->words[i].code_length;

		PB_DEBUG1(errmsg("decode():Sequence specific code copied"));
	}

	if (input_header->has_index)
	{
		int start_entry_no;

		start_entry_no = (start_position + 1) / PB_INDEX_PART_SIZE - 1;
		if (start_entry_no >= 0)
		{
			Varlena* data_slice = (Varlena*) PG_DETOAST_DATUM_SLICE(input,
													  sizeof(PB_CompressedSequence) - VARHDRSZ +
													  sizeof(PB_Codeword) * input_header->n_symbols +
													  sizeof(PB_IndexEntry) * start_entry_no,
													  sizeof(PB_IndexEntry));

			start_entry = palloc0(sizeof(PB_IndexEntry));
			memcpy(start_entry, VARDATA_ANY(data_slice), sizeof(PB_IndexEntry));
			pfree(data_slice);

			PB_DEBUG1(errmsg("decode(): index found, uses entry no %d", start_entry_no));
		}
	}

	if (codeset->n_swapped_symbols > 0)
	{
		if (codeset->uses_rle)
			decode_pc_swp_rle_idx(input,
								  output,
								  start_position,
								  out_length,
								  start_entry,
								  codeset);
		else
			decode_pc_swp_idx(input,
							  output,
							  start_position,
							  out_length,
							  start_entry,
							  codeset);
	}
	else
	{
		if (codeset->uses_rle)
			decode_pc_rle_idx(input,
							  output,
							  start_position,
							  out_length,
							  start_entry,
							  codeset);
		else
			decode_pc_idx(input,
						  output,
						  start_position,
						  out_length,
						  start_entry,
						  codeset);
	}

	if (start_entry)
		pfree(start_entry);

	if (codeset->is_fixed == FALSE)
		pfree(codeset);

	PB_TRACE(errmsg("<-decode()"));

	pfree(input_header);
}

