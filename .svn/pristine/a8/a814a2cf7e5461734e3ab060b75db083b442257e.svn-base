/*-------------------------------------------------------------------------
*
* Copyright (c) 2013, Max Planck Institute for Marine Microbiology
*
* This software is released under the PostgreSQL License
*
* Author: Michael Schneider <mschneid@mpi-bremen.de>
*
* IDENTIFICATION
*   src/types/bio_functions.c
*
*-------------------------------------------------------------------------
*/

#include "postgres.h"
#include "fmgr.h"

#include "sequence/sequence.h"
#include "sequence/decompression_iteration.h"
#include "sequence/stats.h"
#include "types/rna_sequence.h"
#include "types/aa_sequence.h"
#include "utils/debug.h"

Datum transcribe_dna(PG_FUNCTION_ARGS);
Datum reverse_transcribe_rna(PG_FUNCTION_ARGS);
Datum translate_rna(PG_FUNCTION_ARGS);

/**
 * transcribe_dna()
 * 		Transcribes DNA to RNA.
 *
 * 	PB_CompressedSequence* input : compressed input sequence
 */
PG_FUNCTION_INFO_V1 (transcribe_dna);
Datum transcribe_dna(PG_FUNCTION_ARGS)
{
	PB_CompressedSequence* input = (PB_CompressedSequence*)
			PG_GETARG_VARLENA_P(0);
	PB_CompressedSequence* result;

	uint32 mem_size = VARSIZE(input);

	PB_TRACE(errmsg("->transcribe_dna()"));

	/*
	 * Make a copy of the original.
	 */
	result = palloc0(mem_size);
	memcpy(result, input, mem_size);

	if (result->is_fixed)
	{
		/*
		 * Transcribed DNA is complement RNA
		 */
		result->n_swapped_symbols = result->n_swapped_symbols ^ 0x4;
	}
	else
	{
		PB_Codeword* codewords = PB_COMPRESSED_SEQUENCE_SYMBOL_POINTER(result);
		int i;

		for (i = 0; i < result->n_symbols; i++)
		{
			uint8 symbol = codewords[i].symbol;
			switch (symbol)
			{
			case 'A':
				symbol = 'U';
				break;
			case 'T':
				symbol = 'A';
				break;
			case 'C':
				symbol = 'G';
				break;
			case 'G':
				symbol = 'C';
				break;
			case 'R':
				symbol = 'Y';
				break;
			case 'Y':
				symbol = 'R';
				break;
			case 'M':
				symbol = 'K';
				break;
			case 'K':
				symbol = 'M';
				break;
			case 'D':
				symbol = 'H';
				break;
			case 'H':
				symbol = 'D';
				break;
			case 'V':
				symbol = 'B';
				break;
			case 'B':
				symbol = 'V';
				break;
			case 'a':
				symbol = 'u';
				break;
			case 't':
				symbol = 'a';
				break;
			case 'c':
				symbol = 'g';
				break;
			case 'g':
				symbol = 'c';
				break;
			case 'r':
				symbol = 'y';
				break;
			case 'y':
				symbol = 'r';
				break;
			case 'm':
				symbol = 'k';
				break;
			case 'k':
				symbol = 'm';
				break;
			case 'd':
				symbol = 'h';
				break;
			case 'h':
				symbol = 'd';
				break;
			case 'v':
				symbol = 'b';
				break;
			case 'b':
				symbol = 'v';
				break;
			}
			codewords[i].symbol = symbol;
		}
	}

	PB_TRACE(errmsg("<-transcribe_dna()"));

	PG_RETURN_POINTER(result);
}

/**
 * reverse_transcribe_rna()
 * 		Reverse transcribes RNA to DNA.
 *
 * 	PB_CompressedSequence* input : compressed input sequence
 */
PG_FUNCTION_INFO_V1 (reverse_transcribe_rna);
Datum reverse_transcribe_rna(PG_FUNCTION_ARGS)
{
	PB_CompressedSequence* input = (PB_CompressedSequence*)
			PG_GETARG_VARLENA_P(0);
	PB_CompressedSequence* result;

	uint32 mem_size = VARSIZE(input);

	PB_TRACE(errmsg("->reverse_transcribe_rna()"));

	/*
	 * Make a copy of the original.
	 */
	result = palloc0(mem_size);
	memcpy(result, input, mem_size);

	if (result->is_fixed)
	{
		/*
		 * Reverse transcribed RNA is complement DNA.
		 */
		result->n_swapped_symbols = result->n_swapped_symbols ^ 0x4;
	}
	else
	{
		PB_Codeword* codewords = PB_COMPRESSED_SEQUENCE_SYMBOL_POINTER(result);
		int i;

		for (i = 0; i < result->n_symbols; i++)
		{
			uint8 symbol = codewords[i].symbol;
			switch (symbol)
			{
			case 'A':
				symbol = 'T';
				break;
			case 'U':
				symbol = 'A';
				break;
			case 'C':
				symbol = 'G';
				break;
			case 'G':
				symbol = 'C';
				break;
			case 'R':
				symbol = 'Y';
				break;
			case 'Y':
				symbol = 'R';
				break;
			case 'M':
				symbol = 'K';
				break;
			case 'K':
				symbol = 'M';
				break;
			case 'D':
				symbol = 'H';
				break;
			case 'H':
				symbol = 'D';
				break;
			case 'V':
				symbol = 'B';
				break;
			case 'B':
				symbol = 'V';
				break;
			case 'a':
				symbol = 't';
				break;
			case 'u':
				symbol = 'a';
				break;
			case 'c':
				symbol = 'g';
				break;
			case 'g':
				symbol = 'c';
				break;
			case 'r':
				symbol = 'y';
				break;
			case 'y':
				symbol = 'r';
				break;
			case 'm':
				symbol = 'k';
				break;
			case 'k':
				symbol = 'm';
				break;
			case 'd':
				symbol = 'h';
				break;
			case 'h':
				symbol = 'd';
				break;
			case 'v':
				symbol = 'b';
				break;
			case 'b':
				symbol = 'v';
				break;
			}
			codewords[i].symbol = symbol;
		}
	}

	PB_TRACE(errmsg("<-reverse_transcribe_rna()"));

	PG_RETURN_POINTER(result);
}

/**
 * translate_rna()
 * 		Translate an RNA sequence to AA sequence.
 *
 * 	The translation table must be a 64 chars long text.
 *
 * 	PB_CompressedSequence* input : RNA sequence
 * 	text* table_input : translation table
 */
PG_FUNCTION_INFO_V1 (translate_rna);
Datum translate_rna(PG_FUNCTION_ARGS)
{
	PB_CompressedSequence* input = (PB_CompressedSequence*) PG_GETARG_VARLENA_P(0);
	text* table_input = (text*) PG_GETARG_VARLENA_P(1);

	uint8* raw_output = palloc0(input->sequence_length / 3 + 1);

	PB_CompressedSequence* result;

	PB_TRACE(errmsg("->translate_rna()"));

	if (VARSIZE(table_input) - VARHDRSZ != 64)
		ereport(ERROR, (errmsg("translation table has invalid size"),
						errdetail("The length of the this translation table is %u.", (VARSIZE(table_input) - VARHDRSZ)),
						errhint("A translation table must be of length 64.")));

	/*
	 * Translate, write output as cstring
	 */
	{
		uint8 c = 0;
		int codon = 0;
		int pos = 0;

		uint8* output_pointer = raw_output;
		uint8* transl_table = (uint8*) VARDATA(table_input);

		PB_BEGIN_DECODE((Varlena*) input, 0, input->sequence_length, get_fixed_rna_codes(), (c)) {
			codon *= 4;

			switch (TO_UPPER(c)) {
			case 'U':
				break;
			case 'C':
				codon += 1;
				break;
			case 'A':
				codon += 2;
				break;
			case 'G':
				codon += 3;
				break;
			default:
				codon = 64;
			}

			pos++;
			if (pos % 3 == 0)
			{
				if (codon < 64)
				{
					uint8 out = transl_table[codon];
					if (out == '*')
						out = '\0';

					*output_pointer = out;
				}
				else
					*output_pointer = 'X';

				output_pointer++;
				codon = 0;
			}

			PB_DEBUG2(errmsg("pos:%d c:%c codon:%d", pos, c, codon));
		} PB_END_DECODE

		*output_pointer = '\0';
	}

	/*
	 * Compress raw output
	 */
	{
		PB_SequenceInfo* info = get_sequence_info_cstring(raw_output, PB_SEQUENCE_INFO_CASE_INSENSITIVE | PB_SEQUENCE_INFO_WITHOUT_RLE);
		PB_AaSequenceTypMod typmod;

		typmod.case_sensitive = PB_AA_TYPMOD_CASE_INSENSITIVE;
		typmod.restricting_alphabet = PB_AA_TYPMOD_IUPAC;

		result = compress_aa_sequence(raw_output, typmod, info);
	}

	PB_TRACE(errmsg("<-translate_rna()"));

	PG_RETURN_POINTER(result);
}
