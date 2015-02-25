/*-------------------------------------------------------------------------
*
* Copyright (c) 2013, Max Planck Institute for Marine Microbiology
*
* This software is released under the PostgreSQL License
*
* Author: Michael Schneider <mschneid@mpi-bremen.de>
*
* IDENTIFICATION
*   src/sequence/generation.c
*
*-------------------------------------------------------------------------
*/
#include <stdlib.h>

#include "postgres.h"
#include "fmgr.h"

#include "sequence/sequence.h"
#include "types/alphabet.h"
#include "utils/debug.h"

Datum generate_sequence (PG_FUNCTION_ARGS);

/**
 * generate_sequence()
 * 		Generate a random sequence.
 *
 * 	PB_Alphabet* input : Alphabet
 * 	int32 : sequence length
 */
PG_FUNCTION_INFO_V1 (generate_sequence);
Datum generate_sequence (PG_FUNCTION_ARGS) {
	PB_Alphabet* input =
			(PB_Alphabet*) PG_DETOAST_DATUM(PG_GETARG_POINTER(0));
	int32 sequence_length = PG_GETARG_INT32(1);

	text* result;
	char* output_pointer;

	PB_Symbol* symbols;
	PB_SymbolProbability* probabilities;
	PB_SymbolProbability* cumulative_probabilities;
	PB_SymbolProbability sum;

	int i,j;
	int mem_size;

	PB_TRACE(errmsg("->generate_sequence() of length %d", sequence_length));

	symbols = PB_ALPHABET_SYMBOL_POINTER(input);

	if (PB_ALPHABET_TYPE(input) == PB_ALPHABET_WITH_PROBABILITIES)
		probabilities = PB_ALPHABET_SYMBOL_PROBABILITY_POINTER(input);
	else
	{
		/*
		 * If there are no probabilities in the alphabet, assume equal distribution
		 */
		PB_SymbolProbability equal_probability = 1.0f / PB_ALPHABET_SIZE(input);
		probabilities = palloc0(PB_ALPHABET_SIZE(input) * sizeof(PB_SymbolProbability));
		for (i = 0; i < PB_ALPHABET_SIZE(input); i++)
		{
			probabilities[i] = equal_probability;
		}
	}

	/*
	 * Compute cumulative probabilities.
	 */
	cumulative_probabilities = palloc0(PB_ALPHABET_SIZE(input) * sizeof(PB_SymbolProbability));

	sum = 0.0f;
	for (i = 0; i < PB_ALPHABET_SIZE(input); i++)
	{
		cumulative_probabilities[i] = sum;
		sum += probabilities[i];
		PB_DEBUG2(errmsg("generate_sequence(): symbol %c cumprob %f", symbols[i],sum));
	}

	mem_size = VARHDRSZ + sequence_length * sizeof(PB_Symbol);
	result = palloc0(mem_size);
	SET_VARSIZE(result,mem_size);
	output_pointer = VARDATA(result);

	/*
	 * Generate random number between 0 and 1, use cumulative
	 * probabilites to find the respective symbol.
	 */
	for (j = sequence_length - 1; j >= 0; output_pointer++, j--)
	{
		PB_SymbolProbability random_number;

		random_number = (PB_SymbolProbability) rand() / (PB_SymbolProbability) RAND_MAX;

		for (i = PB_ALPHABET_SIZE(input) - 1; i >= 0 ; i--)
		{
			if (random_number >= cumulative_probabilities[i])
			{
				*output_pointer = symbols[i];
				break;
			}
		}
	}

	pfree(cumulative_probabilities);

	if (PB_ALPHABET_TYPE(input) == PB_ALPHABET_WITHOUT_PROBABILITIES)
		pfree(probabilities);

	PB_TRACE(errmsg("<-generate_sequence()"));

	PG_RETURN_TEXT_P(result);
}

