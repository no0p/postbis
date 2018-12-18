/*-------------------------------------------------------------------------
*
* Copyright (c) 2013, Max Planck Institute for Marine Microbiology
*
* This software is released under the PostgreSQL License
*
* Author: Michael Schneider <mschneid@mpi-bremen.de>
*
* IDENTIFICATION
*   src/sequence/stats.c
*
*-------------------------------------------------------------------------
*/
#include "postgres.h"

#include "sequence/sequence.h"
#include "utils/debug.h"

#include "sequence/stats.h"

/*
 * local function declarations
 */

static PB_SequenceInfo* get_sequence_info_cstring_rle_case_sensitive(const uint8* input);
static PB_SequenceInfo* get_sequence_info_cstring_rle_case_insensitive(const uint8* input);
static PB_SequenceInfo* get_sequence_info_cstring_case_sensitive(const uint8* input);
static PB_SequenceInfo* get_sequence_info_cstring_case_insensitive(const uint8* input);
static PB_SequenceInfo* get_sequence_info_text_rle_case_sensitive(const text* input);
static PB_SequenceInfo* get_sequence_info_text_rle_case_insensitive(const text* input);
static PB_SequenceInfo* get_sequence_info_text_case_sensitive(const text* input);
static PB_SequenceInfo* get_sequence_info_text_case_insensitive(const text* input);

/*
 * macros for heap-sort
 */

/*
 * determines index of left child of a node in a heap
 */
#define HEAP_LEFT(i)	(2 * i + 1)

/*
 * determines index of right child of a node in a heap
 */
#define HEAP_RIGHT(i)	(2 * i + 2)

/*
 * determines index of parent of a node in a heap
 */
#define HEAP_PARENT(i)	((i - 1) / 2)

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
					  uint64* ascii_bitmap_high)
{
	unsigned int i;

	uint8* heap;
	unsigned int heap_size;

	uint64 bitmap_high;
	uint64 bitmap_low;

	PB_TRACE(errmsg("->collect_alphabet"));

	heap = (uint8*) palloc0(PB_ASCII_SIZE * sizeof(uint8));

	/*
	 * Build heap, by inserting element-wise.
	 */
	heap_size = 0;
	for (i = 0; i < PB_ASCII_SIZE; i++)
	{
		if (0 != frequencies[i])
		{
			unsigned int j = heap_size;
			heap[heap_size] = i;
			heap_size++;
			while (j > 0 && frequencies[heap[j]] > frequencies[heap[HEAP_PARENT(j)]])
			{
				unsigned int swap_temp = heap[j];
				heap[j] = heap[HEAP_PARENT(j)];
				heap[HEAP_PARENT(j)] = swap_temp;
				j = HEAP_PARENT(j);
			}
		}
	}

	*symbols = (uint8*) palloc0(heap_size * sizeof(uint8));
	*n_symbols = heap_size;

	bitmap_high = 0;
	bitmap_low = 0;

	/*
	 * Copy into result, using ExtractMax
	 */
	for (i = 0; i < *n_symbols; i++)
	{
		unsigned int j;

		/*
		 * Maximal symbol is a position 0.
		 */
		uint8 c = heap[0];

		(*symbols)[i] = c;

		/*
		 * Update bitmap.
		 */
		if (c >= 64)
		{
			bitmap_high = bitmap_high | ((uint64) 1) << (c - 64);
		}
		else
		{
			bitmap_low = bitmap_low | ((uint64) 1) << c;
		}

		/*
		 * Remove maximal element. Put minimal element on top.
		 */
		heap[0] = heap[heap_size - 1];
		heap_size--;

		/*
		 * Heapify.
		 */
		j = 0;
		do
		{
			unsigned int k = j;

			if (HEAP_LEFT(k) < heap_size && frequencies[heap[HEAP_LEFT(k)]] > frequencies[heap[k]])
			{
				j = HEAP_LEFT(k);
			}
			if (HEAP_RIGHT(k) < heap_size && frequencies[heap[HEAP_RIGHT(k)]] > frequencies[heap[j]])
			{
				j = HEAP_RIGHT(k);
			}
			if (k != j)
			{
				unsigned int swap_temp = heap[k];
				heap[k] = heap[j];
				heap[j] = swap_temp;
			}
			else
			{
				break;
			}
		} while ( 1 );
	}

	pfree(heap);

	if (ascii_bitmap_high)
	{
		*ascii_bitmap_high = bitmap_high;
	}
	if (ascii_bitmap_low)
	{
		*ascii_bitmap_low = bitmap_low;
	}

	PB_TRACE(errmsg("<-collect_alphabet()"));
}

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
void check_ascii(PB_SequenceInfo* sequence_info)
{
	int i;
	int nonAsciiChars = 0;

	PB_TRACE(errmsg("->check_ascii"));

	if (0 < sequence_info->frequencies[0])
	{
		int null_terminators = sequence_info->frequencies[0];

		PB_SEQUENCE_INFO_PFREE(sequence_info);

		ereport(ERROR,(errmsg("input sequence violates alphabet constraints"),
				errdetail("Failing datum contains %d null-terminator(s).", null_terminators)));
	}

	for (i = 194; i <= 244; i++)
	{
		nonAsciiChars += sequence_info->frequencies[i];
	}

	if (0 != nonAsciiChars)
	{
		PB_SEQUENCE_INFO_PFREE(sequence_info);

		ereport(ERROR,(errmsg("input sequence violates alphabet constraints"),
					errdetail("Failing datum contains %d non-ASCII character(s).", nonAsciiChars)));
	}

	PB_TRACE(errmsg("<-check_ascii()"));
}

static PB_SequenceInfo* get_sequence_info_cstring_rle_case_sensitive(const uint8* input)
{
	uint8 current = '\000';
	uint8 recent = '\000';
	unsigned int repeated_chars = -1;
	uint8* input_pointer = (uint8*) input;
	int64 temp_len;

	PB_SequenceInfo* result;

	uint32* frequencies;
	uint32* rle_frequencies;

	result = (PB_SequenceInfo*) palloc0(sizeof(PB_SequenceInfo));
	result->rle_info = (PB_RleInfo*) palloc0(sizeof(PB_RleInfo));
	result->ignore_case = false;

	frequencies = (uint32*) &(result->frequencies);
	rle_frequencies = (uint32*) &(result->rle_info->rle_frequencies);

	while (*input_pointer != '\000')
	{
		current = *input_pointer;
		input_pointer++;
		frequencies[current]++;
		repeated_chars++;

		if (current != recent)
		{
			/*
			 * a series of length 'repeated_chars' equal characters ends
			 */
			if (repeated_chars < PB_MIN_RUN_LENGTH)
			{
				/*
				 * the series is too short for using rle
				 */
				rle_frequencies[recent] += repeated_chars;
				recent = current;
				repeated_chars = 0;
			}
			else
			{
				 /*
				  * number of full PB_MAX_RUN_LENGTH rle 'blocks' needed
				  */
				const int rle_blocks = repeated_chars / (PB_MAX_RUN_LENGTH - 1);
				/*
				 *  remaining characters, that are not covered by full rle blocks
				 */
				const int remainder = repeated_chars % (PB_MAX_RUN_LENGTH - 1);

				rle_frequencies[PB_RUN_LENGTH_SYMBOL] += rle_blocks;
				rle_frequencies[recent] += rle_blocks;

				if (remainder >= PB_MIN_RUN_LENGTH)
				{
					/*
					 * enough remaining characters to put them into an rle block
					 */
					rle_frequencies[PB_RUN_LENGTH_SYMBOL]++;
					rle_frequencies[recent]++;
				}
				else
				{
					/*
					 * too few remaining characters to use RLE
					 */
					rle_frequencies[recent] += remainder;
				}

				recent = current;
				repeated_chars = 0;
			}
		}
	}

	/*
	 * processing last character
	 */
	repeated_chars++;
	if (repeated_chars < PB_MIN_RUN_LENGTH)
	{
		rle_frequencies[recent] += repeated_chars;
	}
	else
	{
		const int rle_blocks = repeated_chars / (PB_MAX_RUN_LENGTH - 1);
		const int remainder = repeated_chars % (PB_MAX_RUN_LENGTH - 1);

		rle_frequencies[PB_RUN_LENGTH_SYMBOL] += rle_blocks;
		rle_frequencies[recent] += rle_blocks;

		if (remainder >= PB_MIN_RUN_LENGTH)
		{
			rle_frequencies[PB_RUN_LENGTH_SYMBOL]++;
			rle_frequencies[recent]++;
		}
		else
		{
			rle_frequencies[recent] += remainder;
		}
	}

	temp_len = (input_pointer - (uint8*) input);

	if (temp_len >= PB_MAX_INPUT_SEQUENCE_LENGTH)
	{
		ereport(ERROR,(errmsg("input sequence violates length constraints"),
				errdetail("Maximum is %ld characters. This sequence has %ld characters,", PB_MAX_INPUT_SEQUENCE_LENGTH, temp_len)));
	}

	check_ascii(result);

	result->sequence_length = (uint32) temp_len;
	result->ignore_case = false;
	collect_alphabet((uint32*) &(result->frequencies), &(result->n_symbols), &(result->symbols), &(result->ascii_bitmap_low), &(result->ascii_bitmap_high));
	collect_alphabet((uint32*) &(result->rle_info->rle_frequencies), &(result->rle_info->n_symbols), &(result->rle_info->symbols), NULL, NULL);

	return result;
}

static PB_SequenceInfo* get_sequence_info_cstring_rle_case_insensitive(const uint8* input)
{
	uint8 current = '\000';
	uint8 recent = '\000';
	unsigned int repeated_chars = -1;
	uint8* input_pointer = (uint8*) input;
	int64 temp_len;

	PB_SequenceInfo* result;

	uint32* frequencies;
	uint32* rle_frequencies;

	result = (PB_SequenceInfo*) palloc0(sizeof(PB_SequenceInfo));
	result->rle_info = (PB_RleInfo*) palloc0(sizeof(PB_RleInfo));
	result->ignore_case = true;

	frequencies = (uint32*) &(result->frequencies);
	rle_frequencies = (uint32*) &(result->rle_info->rle_frequencies);

	while (*input_pointer != '\000')
	{
		current = TO_UPPER(*input_pointer);
		input_pointer++;
		frequencies[current]++;
		repeated_chars++;

		if (current != recent)
		{
			/*
			 * a series of length 'repeated_chars' equal characters ends
			 */
			if (repeated_chars < PB_MIN_RUN_LENGTH)
			{
				/*
				 * the series is too short for using rle
				 */
				rle_frequencies[recent] += repeated_chars;
				recent = current;
				repeated_chars = 0;
			}
			else
			{
				/*
				 * number of full PB_MAX_RUN_LENGTH rle 'blocks' needed
				 */
				const int rle_blocks = repeated_chars / (PB_MAX_RUN_LENGTH - 1);
				/*
				 *  remaining characters, that are not covered by full rle blocks
				 */
				const int remainder = repeated_chars % (PB_MAX_RUN_LENGTH - 1);

				rle_frequencies[PB_RUN_LENGTH_SYMBOL] += rle_blocks;
				rle_frequencies[recent] += rle_blocks;

				if (remainder >= PB_MIN_RUN_LENGTH)
				{
					/*
					 * enough remaining characters to put them into an rle block
					 */
					rle_frequencies[PB_RUN_LENGTH_SYMBOL]++;
					rle_frequencies[recent]++;
				}
				else
				{
					/*
					 * too few remaining characters to use RLE
					 */
					rle_frequencies[recent] += remainder;
				}
				recent = current;
				repeated_chars = 0;
			}
		}
	}

	/*
	 * processing last character
	 */
	repeated_chars++;
	if (repeated_chars < PB_MIN_RUN_LENGTH)
	{
		rle_frequencies[recent] += repeated_chars;
	}
	else
	{
		const int rle_blocks = repeated_chars / (PB_MAX_RUN_LENGTH - 1);
		const int remainder = repeated_chars % (PB_MAX_RUN_LENGTH - 1);
		rle_frequencies[PB_RUN_LENGTH_SYMBOL] += rle_blocks;
		rle_frequencies[recent] += rle_blocks;
		if (remainder >= PB_MIN_RUN_LENGTH)
		{
			rle_frequencies[PB_RUN_LENGTH_SYMBOL]++;
			rle_frequencies[recent]++;
		}
		else
		{
			rle_frequencies[recent] += remainder;
		}
	}

	temp_len = (input_pointer - (uint8*) input);

	if (temp_len >= PB_MAX_INPUT_SEQUENCE_LENGTH)
	{
		ereport(ERROR,(errmsg("input sequence violates length constraints"),
				errdetail("Maximum is %ld characters. This sequence has %ld characters,", PB_MAX_INPUT_SEQUENCE_LENGTH, temp_len)));
	}

	check_ascii(result);

	result->sequence_length = (uint32) temp_len;
	result->ignore_case = true;
	collect_alphabet((uint32*) &(result->frequencies), &(result->n_symbols), &(result->symbols), &(result->ascii_bitmap_low), &(result->ascii_bitmap_high));
	collect_alphabet((uint32*) &(result->rle_info->rle_frequencies), &(result->rle_info->n_symbols), &(result->rle_info->symbols), NULL, NULL);

	return result;
}

static PB_SequenceInfo* get_sequence_info_cstring_case_sensitive(const uint8* input)
{
	uint8* input_pointer = (uint8*) input;

	int64 temp_len;

	PB_SequenceInfo* result = (PB_SequenceInfo*) palloc0(sizeof(PB_SequenceInfo));

	uint32* frequencies = (uint32*) &(result->frequencies);

	while (*input_pointer != '\000')
	{
		frequencies[*input_pointer]++;
		input_pointer++;
	}

	temp_len = (input_pointer - (uint8*) input);

	if (temp_len >= PB_MAX_INPUT_SEQUENCE_LENGTH)
	{
		ereport(ERROR,(errmsg("input sequence violates length constraints"),
				errdetail("Maximum is %ld characters. This sequence has %ld characters,", PB_MAX_INPUT_SEQUENCE_LENGTH, temp_len)));
	}

	check_ascii(result);

	result->sequence_length = (uint32) temp_len;
	result->ignore_case = false;
	collect_alphabet((uint32*) &(result->frequencies),&(result->n_symbols),&(result->symbols), &(result->ascii_bitmap_low), &(result->ascii_bitmap_high));

	return result;
}

/*
 * local functions
 */
static PB_SequenceInfo* get_sequence_info_cstring_case_insensitive(const uint8* input)
{
	uint8* input_pointer = (uint8*) input;

	int64 temp_len;
	int i;

	PB_SequenceInfo* result = (PB_SequenceInfo*) palloc0(sizeof(PB_SequenceInfo));

	uint32* frequencies = (uint32*) &(result->frequencies);

	while (*input_pointer != '\000')
	{
		frequencies[*input_pointer]++;
		input_pointer++;
	}

	temp_len = (input_pointer - (uint8*) input);

	if (temp_len >= PB_MAX_INPUT_SEQUENCE_LENGTH)
	{
		ereport(ERROR,(errmsg("input sequence violates length constraints"),
				errdetail("Maximum is %ld characters. This sequence has %ld characters,", PB_MAX_INPUT_SEQUENCE_LENGTH, temp_len)));
	}

	for (i = 97; i <= 122; i++)
	{
		frequencies[i - 32] += frequencies[i];
		frequencies[i] = 0;
	}

	check_ascii(result);

	result->sequence_length = (uint32) temp_len;
	result->ignore_case = true;
	collect_alphabet((uint32*) &(result->frequencies),&(result->n_symbols),&(result->symbols), &(result->ascii_bitmap_low), &(result->ascii_bitmap_high));

	return result;
}

static PB_SequenceInfo* get_sequence_info_text_rle_case_sensitive(const text* input)
{
	uint8 current = '\000';
	uint8 recent = '\000';
	unsigned int repeated_chars = -1;
	int i = VARSIZE_ANY_EXHDR(input) - 1;
	uint8* input_pointer = (uint8*) VARDATA_ANY(input);

	PB_SequenceInfo* result;

	uint32* frequencies;
	uint32* rle_frequencies;

	result = (PB_SequenceInfo*) palloc0(sizeof(PB_SequenceInfo));
	result->rle_info = (PB_RleInfo*) palloc0(sizeof(PB_RleInfo));
	result->ignore_case = false;

	frequencies = (uint32*) &(result->frequencies);
	rle_frequencies = (uint32*) &(result->rle_info->rle_frequencies);

	while (i >= 0)
	{
		current = *input_pointer;
		input_pointer++;
		i--;
		frequencies[current]++;
		repeated_chars++;

		if (current != recent)
		{
			/*
			 * a series of length 'repeated_chars' equal characters ends
			 */
			if (repeated_chars < PB_MIN_RUN_LENGTH)
			{
				/*
				 * the series is too short for using rle
				 */
				rle_frequencies[recent] += repeated_chars;
				recent = current;
				repeated_chars = 0;
			}
			else
			{
				/*
				 * number of full PB_MAX_RUN_LENGTH rle 'blocks' needed
				 */
				const int rle_blocks = repeated_chars / (PB_MAX_RUN_LENGTH - 1);
				/*
				 *  remaining characters, that are not covered by full rle blocks
				 */
				const int remainder = repeated_chars % (PB_MAX_RUN_LENGTH - 1);

				rle_frequencies[PB_RUN_LENGTH_SYMBOL] += rle_blocks;
				rle_frequencies[recent] += rle_blocks;

				if (remainder >= PB_MIN_RUN_LENGTH)
				{
					/*
					 * enough remaining characters to put them into an rle block
					 */
					rle_frequencies[PB_RUN_LENGTH_SYMBOL]++;
					rle_frequencies[recent]++;
				}
				else
				{
					/*
					 * too few remaining characters to use RLE
					 */
					rle_frequencies[recent] += remainder;
				}
				recent = current;
				repeated_chars = 0;
			}
		}
	}

	/*
	 * processing last character
	 */
	repeated_chars++;
	if (repeated_chars < PB_MIN_RUN_LENGTH)
	{
		rle_frequencies[recent] += repeated_chars;
	}
	else
	{
		const int rle_blocks = repeated_chars / (PB_MAX_RUN_LENGTH - 1);
		const int remainder = repeated_chars % (PB_MAX_RUN_LENGTH - 1);
		rle_frequencies[PB_RUN_LENGTH_SYMBOL] += rle_blocks;
		rle_frequencies[recent] += rle_blocks;
		if (remainder >= PB_MIN_RUN_LENGTH)
		{
			rle_frequencies[PB_RUN_LENGTH_SYMBOL]++;
			rle_frequencies[recent]++;
		}
		else
		{
			rle_frequencies[recent] += remainder;
		}
	}

	check_ascii(result);

	result->sequence_length = VARSIZE_ANY_EXHDR(input);
	result->ignore_case = false;
	collect_alphabet((uint32*) &(result->frequencies), &(result->n_symbols), &(result->symbols), &(result->ascii_bitmap_low), &(result->ascii_bitmap_high));
	collect_alphabet((uint32*) &(result->rle_info->rle_frequencies), &(result->rle_info->n_symbols), &(result->rle_info->symbols), NULL, NULL);

	return result;
}

static PB_SequenceInfo* get_sequence_info_text_rle_case_insensitive(const text* input)
{
	uint8 current = '\000';
	uint8 recent = '\000';
	unsigned int repeated_chars = -1;
	int i = VARSIZE_ANY_EXHDR(input) - 1;
	uint8* input_pointer = (uint8*) VARDATA_ANY(input);

	PB_SequenceInfo* result;

	uint32* frequencies;
	uint32* rle_frequencies;

	result = (PB_SequenceInfo*) palloc0(sizeof(PB_SequenceInfo));
	result->rle_info = (PB_RleInfo*) palloc0(sizeof(PB_RleInfo));
	result->ignore_case = true;

	frequencies = (uint32*) &(result->frequencies);
	rle_frequencies = (uint32*) &(result->rle_info->rle_frequencies);

	while (i >= 0)
	{
		current = TO_UPPER(*input_pointer);
		input_pointer++;
		i--;
		frequencies[current]++;
		repeated_chars++;

		if (current != recent)
		{
			/*
			 * a series of length 'repeated_chars' equal characters ends
			 */
			if (repeated_chars < PB_MIN_RUN_LENGTH)
			{
				/*
				 * the series is too short for using rle
				 */
				rle_frequencies[recent] += repeated_chars;
				recent = current;
				repeated_chars = 0;
			}
			else
			{
				/*
				 * number of full PB_MAX_RUN_LENGTH rle 'blocks' needed
				 */
				const int rle_blocks = repeated_chars / (PB_MAX_RUN_LENGTH - 1);
				/*
				 *  remaining characters, that are not covered by full rle blocks
				 */
				const int remainder = repeated_chars % (PB_MAX_RUN_LENGTH - 1);

				rle_frequencies[PB_RUN_LENGTH_SYMBOL] += rle_blocks;
				rle_frequencies[recent] += rle_blocks;

				if (remainder >= PB_MIN_RUN_LENGTH)
				{
					/*
					 * enough remaining characters to put them into an rle block
					 */
					rle_frequencies[PB_RUN_LENGTH_SYMBOL]++;
					rle_frequencies[recent]++;
				}
				else
				{
					/*
					 * too few remaining characters to use RLE
					 */
					rle_frequencies[recent] += remainder;
				}
				recent = current;
				repeated_chars = 0;
			}
		}
	}

	/*
	 * processing last character
	 */
	repeated_chars++;
	if (repeated_chars < PB_MIN_RUN_LENGTH)
	{
		rle_frequencies[recent] += repeated_chars;
	}
	else
	{
		const int rle_blocks = repeated_chars / (PB_MAX_RUN_LENGTH - 1);
		const int remainder = repeated_chars % (PB_MAX_RUN_LENGTH - 1);
		rle_frequencies[PB_RUN_LENGTH_SYMBOL] += rle_blocks;
		rle_frequencies[recent] += rle_blocks;
		if (remainder >= PB_MIN_RUN_LENGTH)
		{
			rle_frequencies[PB_RUN_LENGTH_SYMBOL]++;
			rle_frequencies[recent]++;
		}
		else
		{
			rle_frequencies[recent] += remainder;
		}
	}

	check_ascii(result);

	result->sequence_length = VARSIZE_ANY_EXHDR(input);
	result->ignore_case = true;
	collect_alphabet((uint32*) &(result->frequencies), &(result->n_symbols), &(result->symbols), &(result->ascii_bitmap_low), &(result->ascii_bitmap_high));
	collect_alphabet((uint32*) &(result->rle_info->rle_frequencies), &(result->rle_info->n_symbols), &(result->rle_info->symbols), NULL, NULL);

	return result;
}

static PB_SequenceInfo* get_sequence_info_text_case_sensitive(const text* input)
{
	uint8* input_pointer = (uint8*) VARDATA_ANY(input);
	int i = VARSIZE_ANY_EXHDR(input) - 1;

	PB_SequenceInfo* result = (PB_SequenceInfo*) palloc0(sizeof(PB_SequenceInfo));

	uint32* frequencies = (uint32*) &(result->frequencies);

	while (i >= 0)
	{
		frequencies[*input_pointer]++;
		input_pointer++;
		i--;
	}

	check_ascii(result);

	result->sequence_length = VARSIZE_ANY_EXHDR(input);
	result->ignore_case = false;
	collect_alphabet((uint32*) &(result->frequencies), &(result->n_symbols), &(result->symbols), &(result->ascii_bitmap_low), &(result->ascii_bitmap_high));

	return result;
}

static PB_SequenceInfo* get_sequence_info_text_case_insensitive(const text* input)
{
	uint8* input_pointer = (uint8*) VARDATA_ANY(input);
	int i = VARSIZE_ANY_EXHDR(input) - 1;

	PB_SequenceInfo* result = (PB_SequenceInfo*) palloc0(sizeof(PB_SequenceInfo));

	uint32* frequencies = (uint32*) &(result->frequencies);

	while (i >= 0)
	{
		frequencies[*input_pointer]++;
		input_pointer++;
		i--;
	}

	for (i = 97; i <= 122; i++)
	{
		frequencies[i - 32] += frequencies[i];
		frequencies[i] = 0;
	}

	check_ascii(result);

	result->sequence_length = VARSIZE_ANY_EXHDR(input);
	result->ignore_case = true;
	collect_alphabet((uint32*) &(result->frequencies), &(result->n_symbols), &(result->symbols), &(result->ascii_bitmap_low), &(result->ascii_bitmap_high));

	return result;
}

/*
 * public functions
 */

/**
 * get_sequence_info_cstring()
 * 		Obtain sequence length, symbol frequencies and alphabet.
 *
 * 		uint8* input : a null-terminated input string
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
										   int mode)
{
	PB_SequenceInfo* result = NULL;

	PB_TRACE(errmsg("->get_sequence_info_cstring(): mode = %d", mode));

	if (mode & PB_SEQUENCE_INFO_WITH_RLE)
	{
		if (mode & PB_SEQUENCE_INFO_CASE_SENSITIVE) {
			result = get_sequence_info_cstring_rle_case_sensitive(input);
		}
		else
		{
			result = get_sequence_info_cstring_rle_case_insensitive(input);
		}
	}
	else
	{
		if (mode & PB_SEQUENCE_INFO_CASE_SENSITIVE) {
			result = get_sequence_info_cstring_case_sensitive(input);
		}
		else
		{
			result = get_sequence_info_cstring_case_insensitive(input);
		}
	}

	PB_TRACE(errmsg("<-get_sequence_info_cstring()"))

	return result;
}

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
										int mode)
{
	PB_SequenceInfo* result = NULL;

	PB_TRACE(errmsg("->get_sequence_info_text(): mode = %d", mode))

	if (mode & PB_SEQUENCE_INFO_WITH_RLE)
	{
		if (mode & PB_SEQUENCE_INFO_CASE_SENSITIVE) {
			result = get_sequence_info_text_rle_case_sensitive(input);
		}
		else
		{
			result = get_sequence_info_text_rle_case_insensitive(input);
		}
	}
	else
	{
		if (mode & PB_SEQUENCE_INFO_CASE_SENSITIVE) {
			result = get_sequence_info_text_case_sensitive(input);
		}
		else
		{
			result = get_sequence_info_text_case_insensitive(input);
		}
	}

	PB_TRACE(errmsg("<-get_sequence_info_text()"))

	return result;
}
