/*-------------------------------------------------------------------------
*
* Copyright (c) 2013, Max Planck Institute for Marine Microbiology
*
* This software is released under the PostgreSQL License
*
* Author: Michael Schneider <mschneid@mpi-bremen.de>
*
* IDENTIFICATION
*   src/sequence/functions.c
*
*-------------------------------------------------------------------------
*/
#include "math.h"

#include "postgres.h"
#include "fmgr.h"
#include "utils/array.h"
#include "catalog/pg_type.h"
#include "access/tuptoaster.h"

#include "sequence/sequence.h"
#include "sequence/stats.h"
#include "sequence/code_set_creation.h"
#include "sequence/sequence.h"
#include "sequence/decompression_iteration.h"
#include "utils/debug.h"

#include "sequence/functions.h"

/**
 * reverse()
 * 		Reverses a detoasted compressed sequence.
 */
PB_CompressedSequence* reverse(PB_CompressedSequence* sequence, PB_CodeSet** fixed_codesets)
{
	PB_CompressedSequence* result;

	uint8* temp;
	uint8* output_pointer;
	PB_CodeSet* codeset;
	PB_SequenceInfo info;

	PB_TRACE(errmsg("->reverse()"));

	temp = palloc0(sequence->sequence_length + 1);

	output_pointer = temp;
	output_pointer += sequence->sequence_length- 1;

	PB_BEGIN_DECODE((Varlena*) sequence, 0, sequence->sequence_length, fixed_codesets, (*output_pointer)) {
		output_pointer--;
	} PB_END_DECODE

	if (sequence->is_fixed)
	{
		codeset = fixed_codesets[sequence->n_swapped_symbols];

		PB_DEBUG1(errmsg("reverse():uses fixed code with id %u", sequence->n_swapped_symbols));
	}
	else
	{
		int code_size = sizeof(PB_Codeword) * sequence->n_symbols;
		PB_Codeword* code;
		int i;

		codeset = palloc0(sizeof(PB_CodeSet) + code_size);
		codeset->n_symbols = sequence->n_symbols;
		codeset->n_swapped_symbols = sequence->n_swapped_symbols;
		codeset->is_fixed = false;
		codeset->has_equal_length = sequence->has_equal_length;
		codeset->uses_rle = sequence->uses_rle;

		code = PB_COMPRESSED_SEQUENCE_SYMBOL_POINTER(sequence);
		memcpy(codeset->words, code, code_size);

		for (i = 0; i < codeset->n_symbols - codeset->n_swapped_symbols; i++)
			if (codeset->max_codeword_length < codeset->words[i].code_length)
				codeset->max_codeword_length = codeset->words[i].code_length;

		for (i = codeset->n_symbols - codeset->n_swapped_symbols; i < codeset->n_symbols; i++)
			if (codeset->max_codeword_length < codeset->words[i].code_length)
				codeset->max_codeword_length = codeset->words[i].code_length;
	}

	info.sequence_length = sequence->sequence_length;

	result = encode(temp, VARSIZE(sequence), codeset, &info);

	pfree(temp);
	if (!codeset->is_fixed)
		pfree(codeset);

	PB_TRACE(errmsg("<-reverse()"));

	return result;
}

/**
 * sequence_equal()
 * 		Compares two compressed sequences. Returns (-1) if equal, 0 if not.
 * 		This one can be faster than sequence_compare, as it terminates if
 * 		sequences have a different length.
 *
 * 	Varlena* seq1 : first possibly toasted sequence
 * 	Varlena* seq2 : second possibly toasted sequence
 */
bool sequence_equal(Varlena* raw_seq1, Varlena* raw_seq2, PB_CodeSet** fixed_codesets)
{
	PB_CompressedSequence* input_header;
	PB_CompressedSequence* seq1;
	PB_CompressedSequence* seq2;
	uint8* decompressed_sequence;
	int len;

	uint8 c;
	uint8* output_pointer;

	PB_TRACE(errmsg("->sequence_equal()"));

	/* Get length of first sequence */
	input_header = (PB_CompressedSequence*)
					PG_DETOAST_DATUM_SLICE(raw_seq1,0,4);

	len = input_header->sequence_length;

	/* Get length of second sequence */
	input_header = (PB_CompressedSequence*)
					PG_DETOAST_DATUM_SLICE(raw_seq2,0,4);

	/* Compare length */
	if (len != input_header->sequence_length)
	{
		PB_TRACE(errmsg("<-sequence_equal(): exits due to different length"));
		return false;
	}

	seq1 = (PB_CompressedSequence*) PG_DETOAST_DATUM(raw_seq1);
	seq2 = (PB_CompressedSequence*) PG_DETOAST_DATUM(raw_seq2);

	decompressed_sequence = palloc0(len);
	decode((Varlena*) seq1, decompressed_sequence, 0, len, fixed_codesets);

	output_pointer = decompressed_sequence;

	PB_BEGIN_DECODE((Varlena*) seq2, 0, len, fixed_codesets, (c)) {
		if (c != *output_pointer)
		{
			PB_TRACE(errmsg("<-sequence_equal(): exits due to different char at position %ld", output_pointer - decompressed_sequence));
			return false;
		}
		output_pointer++;
	} PB_END_DECODE

	pfree(decompressed_sequence);

	PB_TRACE(errmsg("<-sequence_equal()"));

	return true;
}

/**
 * sequence_compare()
 * 		Compares two compressed sequences. Returns (-1) if first is smaller,
 * 		0 if both are equal, 1 if second sequence is smaller (lexicographically).
 *
 * 	Varlena* seq_a : first possibly toasted sequence
 * 	Varlena* seq_b : second possibly toasted sequence
 */
int sequence_compare(Varlena* seq_a, Varlena* seq_b, PB_CodeSet** fixed_codesets)
{
	PB_CompressedSequence* seq1 = (PB_CompressedSequence*) PG_DETOAST_DATUM(seq_a);
	PB_CompressedSequence* seq2 = (PB_CompressedSequence*) PG_DETOAST_DATUM(seq_b);
	int swap = 1;

	uint8* decompressed_sequence;

	uint8 c;
	uint8* output_pointer;

	PB_TRACE(errmsg("->sequence_compare()"));

	/*
	 * seq1 should be the shorter sequence, as it will be fully decompressed
	 */
	if (seq1->sequence_length > seq2->sequence_length)
	{
		PB_CompressedSequence* swapseq = seq1;
		seq1 = seq2;
		seq2 = swapseq;
		swap = (-1);
	}

	decompressed_sequence = palloc0(seq1->sequence_length);
	decode((Varlena*) seq1, decompressed_sequence, 0, seq1->sequence_length, fixed_codesets);

	output_pointer = decompressed_sequence;

	PB_BEGIN_DECODE((Varlena*) seq2, 0, seq1->sequence_length, fixed_codesets, (c)) {
		if (c != *output_pointer)
		{
			int result = 1;

			if (c > *output_pointer)
				result = (-1);

			result *= swap;

			PB_TRACE(errmsg("<-sequence_compare(): exits due to different char at position %ld", output_pointer - decompressed_sequence));

			return result;
		}
		output_pointer++;
	} PB_END_DECODE

	pfree(decompressed_sequence);

	if (seq1->sequence_length < seq2->sequence_length)
	{
		PB_TRACE(errmsg("<-sequence_compare(): exits, seq1 is shorter"));

		return (-swap);
	}

	PB_TRACE(errmsg("<-sequence_compare()"));

	return 0;
}

/* generated using the AUTODIN II polynomial
 * x^32 + x^26 + x^23 + x^22 + x^16 +
 * x^12 + x^11 + x^10 + x^8 + x^7 + x^5 + x^4 + x^2 + x^1 + 1
 */

static const unsigned int crc32tab[256] = {
 0x00000000, 0x77073096, 0xee0e612c, 0x990951ba,
 0x076dc419, 0x706af48f, 0xe963a535, 0x9e6495a3,
 0x0edb8832, 0x79dcb8a4, 0xe0d5e91e, 0x97d2d988,
 0x09b64c2b, 0x7eb17cbd, 0xe7b82d07, 0x90bf1d91,
 0x1db71064, 0x6ab020f2, 0xf3b97148, 0x84be41de,
 0x1adad47d, 0x6ddde4eb, 0xf4d4b551, 0x83d385c7,
 0x136c9856, 0x646ba8c0, 0xfd62f97a, 0x8a65c9ec,
 0x14015c4f, 0x63066cd9, 0xfa0f3d63, 0x8d080df5,
 0x3b6e20c8, 0x4c69105e, 0xd56041e4, 0xa2677172,
 0x3c03e4d1, 0x4b04d447, 0xd20d85fd, 0xa50ab56b,
 0x35b5a8fa, 0x42b2986c, 0xdbbbc9d6, 0xacbcf940,
 0x32d86ce3, 0x45df5c75, 0xdcd60dcf, 0xabd13d59,
 0x26d930ac, 0x51de003a, 0xc8d75180, 0xbfd06116,
 0x21b4f4b5, 0x56b3c423, 0xcfba9599, 0xb8bda50f,
 0x2802b89e, 0x5f058808, 0xc60cd9b2, 0xb10be924,
 0x2f6f7c87, 0x58684c11, 0xc1611dab, 0xb6662d3d,
 0x76dc4190, 0x01db7106, 0x98d220bc, 0xefd5102a,
 0x71b18589, 0x06b6b51f, 0x9fbfe4a5, 0xe8b8d433,
 0x7807c9a2, 0x0f00f934, 0x9609a88e, 0xe10e9818,
 0x7f6a0dbb, 0x086d3d2d, 0x91646c97, 0xe6635c01,
 0x6b6b51f4, 0x1c6c6162, 0x856530d8, 0xf262004e,
 0x6c0695ed, 0x1b01a57b, 0x8208f4c1, 0xf50fc457,
 0x65b0d9c6, 0x12b7e950, 0x8bbeb8ea, 0xfcb9887c,
 0x62dd1ddf, 0x15da2d49, 0x8cd37cf3, 0xfbd44c65,
 0x4db26158, 0x3ab551ce, 0xa3bc0074, 0xd4bb30e2,
 0x4adfa541, 0x3dd895d7, 0xa4d1c46d, 0xd3d6f4fb,
 0x4369e96a, 0x346ed9fc, 0xad678846, 0xda60b8d0,
 0x44042d73, 0x33031de5, 0xaa0a4c5f, 0xdd0d7cc9,
 0x5005713c, 0x270241aa, 0xbe0b1010, 0xc90c2086,
 0x5768b525, 0x206f85b3, 0xb966d409, 0xce61e49f,
 0x5edef90e, 0x29d9c998, 0xb0d09822, 0xc7d7a8b4,
 0x59b33d17, 0x2eb40d81, 0xb7bd5c3b, 0xc0ba6cad,
 0xedb88320, 0x9abfb3b6, 0x03b6e20c, 0x74b1d29a,
 0xead54739, 0x9dd277af, 0x04db2615, 0x73dc1683,
 0xe3630b12, 0x94643b84, 0x0d6d6a3e, 0x7a6a5aa8,
 0xe40ecf0b, 0x9309ff9d, 0x0a00ae27, 0x7d079eb1,
 0xf00f9344, 0x8708a3d2, 0x1e01f268, 0x6906c2fe,
 0xf762575d, 0x806567cb, 0x196c3671, 0x6e6b06e7,
 0xfed41b76, 0x89d32be0, 0x10da7a5a, 0x67dd4acc,
 0xf9b9df6f, 0x8ebeeff9, 0x17b7be43, 0x60b08ed5,
 0xd6d6a3e8, 0xa1d1937e, 0x38d8c2c4, 0x4fdff252,
 0xd1bb67f1, 0xa6bc5767, 0x3fb506dd, 0x48b2364b,
 0xd80d2bda, 0xaf0a1b4c, 0x36034af6, 0x41047a60,
 0xdf60efc3, 0xa867df55, 0x316e8eef, 0x4669be79,
 0xcb61b38c, 0xbc66831a, 0x256fd2a0, 0x5268e236,
 0xcc0c7795, 0xbb0b4703, 0x220216b9, 0x5505262f,
 0xc5ba3bbe, 0xb2bd0b28, 0x2bb45a92, 0x5cb36a04,
 0xc2d7ffa7, 0xb5d0cf31, 0x2cd99e8b, 0x5bdeae1d,
 0x9b64c2b0, 0xec63f226, 0x756aa39c, 0x026d930a,
 0x9c0906a9, 0xeb0e363f, 0x72076785, 0x05005713,
 0x95bf4a82, 0xe2b87a14, 0x7bb12bae, 0x0cb61b38,
 0x92d28e9b, 0xe5d5be0d, 0x7cdcefb7, 0x0bdbdf21,
 0x86d3d2d4, 0xf1d4e242, 0x68ddb3f8, 0x1fda836e,
 0x81be16cd, 0xf6b9265b, 0x6fb077e1, 0x18b74777,
 0x88085ae6, 0xff0f6a70, 0x66063bca, 0x11010b5c,
 0x8f659eff, 0xf862ae69, 0x616bffd3, 0x166ccf45,
 0xa00ae278, 0xd70dd2ee, 0x4e048354, 0x3903b3c2,
 0xa7672661, 0xd06016f7, 0x4969474d, 0x3e6e77db,
 0xaed16a4a, 0xd9d65adc, 0x40df0b66, 0x37d83bf0,
 0xa9bcae53, 0xdebb9ec5, 0x47b2cf7f, 0x30b5ffe9,
 0xbdbdf21c, 0xcabac28a, 0x53b39330, 0x24b4a3a6,
 0xbad03605, 0xcdd70693, 0x54de5729, 0x23d967bf,
 0xb3667a2e, 0xc4614ab8, 0x5d681b02, 0x2a6f2b94,
 0xb40bbe37, 0xc30c8ea1, 0x5a05df1b, 0x2d02ef8d,
};

#define _CRC32_(crc, ch) (crc = (crc >> 8) ^ crc32tab[(crc ^ (ch)) & 0xff])

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
uint32 sequence_crc32(PB_CompressedSequence* seq, PB_CodeSet** fixed_codesets)
{
	uint32 crc32 = ~((uint32) 0);
	uint8 c;

	PB_BEGIN_DECODE((Varlena*) seq, 0, seq->sequence_length, fixed_codesets, (c)) {
		_CRC32_(crc32, c);
	} PB_END_DECODE

	return ~crc32;
}

typedef struct {
	uint32 pos;
	uint32 len;
	void* next;
} PB_PartialMatch;

/*
 * sequence_strpos()
 * 		Find position of given string.
 * 		Baeza-Yates-Gonnet implementation.
 *
 * 	PB_CompressedSequence* seq : detoasted sequence to search in
 * 	Text* search : detoasted text to search for
 * 	PB_CodeSet** fixed_codesets : fixed codes
 */
uint32 sequence_strpos(PB_CompressedSequence* seq, text* search, PB_CodeSet** fixed_codesets)
{
	uint64 search_vector = 0;
	uint64 match;
	uint64* specific_vectors;
	int search_vector_length = 64;

	PB_SequenceInfo* search_str_info;
	uint32 pos = 0;

	PB_PartialMatch* first_match = NULL;
	PB_PartialMatch* last_match = NULL;

	uint8 c;

	PB_TRACE(errmsg("->sequence_strpos()"));

	/* collect info about pattern */
	search_str_info = get_sequence_info_text(search, PB_SEQUENCE_INFO_CASE_SENSITIVE);

	/* terminate if pattern is longer than sequence */
	if (search_str_info->sequence_length > seq->sequence_length || search_str_info->n_symbols > (seq->is_fixed ? fixed_codesets[seq->n_swapped_symbols]->n_symbols : seq->n_symbols))
	{
		PB_TRACE(errmsg("<-sequence_strpos(): too short pl:%u sl:%u ps:%u ss:%u", search_str_info->sequence_length, seq->sequence_length, search_str_info->n_symbols, (seq->is_fixed ? fixed_codesets[seq->n_swapped_symbols]->n_symbols : seq->n_symbols)));

		return 0;
	}

	/* terminate if pattern contains characters the sequence does not */
	{
		int i = 0;
		PB_Codeword* codewords = (seq->is_fixed ? fixed_codesets[seq->n_swapped_symbols]->words : PB_COMPRESSED_SEQUENCE_SYMBOL_POINTER(seq));
		uint64 bitmap_high = 0;
		uint64 bitmap_low = 0;
		int n_symbols =  (seq->is_fixed ? fixed_codesets[seq->n_swapped_symbols]->n_symbols : seq->n_symbols);

		for (i = 0; i < n_symbols; i++) {
			uint8 c = codewords[i].symbol;
			if (c >= 64)
				bitmap_high = bitmap_high | ((uint64) 1) << (c - 64);
			else
				bitmap_low = bitmap_low | ((uint64) 1) << c;
		}
		PB_DEBUG2(errmsg("alphabet matching pl:%lu sl:%lu ph:%lu sh:%lu", search_str_info->ascii_bitmap_high, bitmap_high, search_str_info->ascii_bitmap_low, bitmap_low));

		if (((search_str_info->ascii_bitmap_high ^ bitmap_high) & search_str_info->ascii_bitmap_high) ||
			((search_str_info->ascii_bitmap_low ^ bitmap_low) & search_str_info->ascii_bitmap_low))
		{
			PB_TRACE(errmsg("<-sequence_strpos(): alphabet mismatch pl:%lu sl:%lu ph:%lu sh:%lu", search_str_info->ascii_bitmap_high, bitmap_high, search_str_info->ascii_bitmap_low, bitmap_low));

			return 0;
		}
	}

	/*
	 * Determine search window length
	 */
	if (search_str_info->sequence_length < search_vector_length)
		search_vector_length = search_str_info->sequence_length;

	/*
	 * Build specific vectors
	 */
	{
		uint64 position_bit = 1;
		int i;
		uint8* input_pointer = (uint8*) VARDATA_ANY(search);

		specific_vectors = palloc0(sizeof(uint64) * PB_ASCII_SIZE);

		for (i = 0; i < search_vector_length; i++, input_pointer++)
		{
			specific_vectors[*input_pointer] |= position_bit;
			position_bit <<= 1;
		}
		match = position_bit >> 1;

		for (i = 0; i < PB_ASCII_SIZE; i++)
		{
			if (specific_vectors[i]) {
				PB_DEBUG2(errmsg("%c: %lx", i, specific_vectors[i]));
			}
		}

		PB_DEBUG2(errmsg("match:%lx", match));
	}

	PB_BEGIN_DECODE((Varlena*) seq, 0, seq->sequence_length, fixed_codesets, (c)) {
		/*
		 * Extend partial matches
		 */
		if (first_match) {
			PB_PartialMatch* current = first_match;
			last_match = NULL;

			while (1) {
				uint8* search_str = (uint8*) VARDATA(search) + current->len;

				if (*search_str == c)
				{
					/* extend match */
					current->len++;
					last_match = current;

					/* check if match long enough to terminate */
					if (current->len == search_str_info->sequence_length)
					{
						int result = current->pos;

						/* free all matches */
						while (first_match) {
							PB_PartialMatch* next = first_match->next;
							pfree(first_match);
							first_match = next;
						}

						pfree(specific_vectors);
						PB_SEQUENCE_INFO_PFREE(search_str_info);

						return result;
					}

					/* go to next in list */
					if (current->next)
						current = current->next;
					else
						break;
				}
				else
				{
					/* delete match */
					if (current->next)
					{
						/* current has predecessor in list */
						PB_PartialMatch* next = current->next;

						if (current == first_match)
							first_match = next;
						else
							last_match->next = next;

						pfree(current);

						current = next;
					}
					else
					{
						/* current is last element in list */
						if (current == first_match)
							first_match = NULL;

						pfree(current);

						break;
					}
				}
			}
		}

		/* shift-and */
		search_vector = ((search_vector << 1) | 1) & specific_vectors[c];
		pos++;

		PB_DEBUG3(errmsg("pos: %u, c:%c sv:%lx", pos, c, search_vector));

		if ((search_vector & match) && pos >= search_vector_length)
		{
			PB_DEBUG2(errmsg("found partial match at pos %u", pos));

			/* partial match found */
			if (search_str_info->sequence_length > search_vector_length)
			{
				/* add partial match to end of list */
				PB_PartialMatch* new_match = palloc0(sizeof(PB_PartialMatch));

				new_match->pos = pos - search_vector_length + 1;
				new_match->len = (uint32) search_vector_length;

				if (last_match)
					last_match->next = new_match;
				else
					last_match = new_match;

				if (!first_match)
					first_match = new_match;
			}
			else
			{
				/* partial match is indeed full match */
				pfree(specific_vectors);
				PB_SEQUENCE_INFO_PFREE(search_str_info);
				return pos - search_vector_length + 1;
			}
		}
	} PB_END_DECODE

	pfree(specific_vectors);
	PB_SEQUENCE_INFO_PFREE(search_str_info);

	PB_TRACE(errmsg("<-sequence_strpos()"));

	return 0;
}
