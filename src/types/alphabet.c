/*-------------------------------------------------------------------------
*
* Copyright (c) 2013, Max Planck Institute for Marine Microbiology
*
* This software is released under the PostgreSQL License
*
* Author: Michael Schneider <mschneid@mpi-bremen.de>
*
* IDENTIFICATION
*   src/types/alphabet.c
*
*-------------------------------------------------------------------------
*/
#include "types/alphabet.h"

#include "postgres.h"
#include "fmgr.h"
#include "utils/builtins.h"

#include "utils/bitarray.h"
#include "sequence/sequence.h"
#include "sequence/stats.h"
#include "sequence/decompression_iteration.h"

/**
 * parse_alphabet_from_text()
 * 		Parses an alphabet from a text. The following formats
 * 		are allowed.
 * 			{[symbol],[...],[symbol]}
 * 			{{[symbol],[...],[symbol]},{[probability],[...],[probability]}}
 */
PB_Alphabet* parse_alphabet_from_text(text* input)
{
	PB_Alphabet* result = NULL;
	PB_Symbol* symbols = NULL;
	PB_SymbolProbability* probabilities = NULL;

	PB_BitArray characters_found;

	int i;

	int n_elements = 0;

	int n_curlies_open = 0;
	int n_curlies_close = 0;
	int n_commas = 0;

	int parse_state;

	char* input_pointer = (char*) VARDATA_ANY(input);
	int text_length = (int) VARSIZE_ANY_EXHDR(input);

	/*
	 * Bitarray for ASCII characters, to check for duplicates
	 */
	characters_found = palloc0(PB_BITARRAY_CALC_MEMSIZE(128));

	/*
	 * Count all curly brackets and commas.
	 */
	for (i = 0; i < text_length; input_pointer++, i++)
	{
		if (*input_pointer == '{')
			n_curlies_open++;
		else if (*input_pointer == '}')
			n_curlies_close++;
		else if (*input_pointer == ',')
			n_commas++;
	}

	if (n_curlies_open == 1 && n_curlies_close == 1)
	{
		/*
		 * One pair of curly brackets means only symbols are
		 * given. The number of elements must be the number
		 * of commas + 1.
		 */
		n_elements = n_commas + 1;

		PB_ALPHABET_CREATE_WITHOUT_PROBABILITIES(result, n_elements);
	}
	else if	(n_curlies_open == 3 && n_curlies_close == 3)
	{
		/*
		 * Three pairs of curly brackets means symbols and
		 * probabilities are given. The number of elements must be the
		 * number of commas - 1 (the comma separating symbols and
		 * probabilities) divided by 2 (for symbols and probablities)
		 * + 1.
		 */
		if (n_commas % 2 == 0)
			ereport(ERROR,(errmsg("multidimensional arrays must have array expressions with matching dimensions")));

		n_elements = (n_commas - 1) / 2 + 1;

		PB_ALPHABET_CREATE_WITH_PROBABILITIES(result, n_elements);

		probabilities = PB_ALPHABET_SYMBOL_PROBABILITY_POINTER(result);
	}
	else
	{
		/*
		 * Any other number of curly brackets implies wrong dimensions
		 * of the input array(s).
		 */
		ereport(ERROR,(errmsg("input array has wrong dimensions"),
					errdetail("Only one- and two-dimensional arrays are acceptable.")));
	}

	symbols = PB_ALPHABET_SYMBOL_POINTER(result);
	input_pointer = (char*) VARDATA_ANY(input);
	n_curlies_open = 0;

	PB_DEBUG1(errmsg("%d elements found", n_elements));

	/*
	 * Parse first array.
	 *
	 * parse_state -4: start state with probabilities
	 * parse_state -3: start state without
	 * parse_state -2: expecting to read a symbol
	 * parse_state -1: symbol read, expecting a comma
	 * parse_state 0: comma read, or end of array
	 */

	if (PB_ALPHABET_TYPE(result) == PB_ALPHABET_WITH_PROBABILITIES)
		parse_state = -4;
	else
		parse_state = -3;

	for (i = 0; i < n_elements; i++)
	{
		while (parse_state < 0)
		{
			char c = *input_pointer;
			input_pointer++;

			PB_DEBUG2(errmsg("Element #%d Character:%c parse_state:%d", i, c, parse_state));

			/*
			 * Ignore spaces
			 */
			if (c == ' ' || c == '\t')
				continue;

			/*
			 * Quotes, control-characters and non-ASCII are not allowed symbols
			 */
			else if (c < 32 || c >= 127 || c == '\'' || c == '"')
				ereport(ERROR,(errmsg("input sequence violates alphabet constraints"),
						errdetail("Failing datum contains symbol \"%c\" (%d).", c, c)));

			/*
			 * Open curlies expected
			 */
			if (parse_state == -4 || parse_state == -3)
			{
				if (c == '{')
				{
					n_curlies_open++;
					parse_state++;
				}
				else
					ereport(ERROR,(errmsg("input array has wrong dimensions"),
								errdetail("Only one- and two-dimensional arrays are acceptable.")));
			}
			/*
			 * Expected to read a symbol.
			 */
			else if (parse_state == -2)
			{
				if (c == ',' || c == '}')
					ereport(ERROR,(errmsg("input array contains null-value")));
				else if (c == '{')
					ereport(ERROR,(errmsg("input array has wrong dimensions"),
								errdetail("Only one- and two-dimensional arrays are acceptable.")));
				else
				{
					if (PB_BITARRAY_GET(characters_found, (unsigned int)c) == 1)
						ereport(ERROR,(errmsg("input array contains duplicate symbol"),
								errdetail("Input array contains duplicate symbol \"%c\".", c)));

					PB_BITARRAY_SET(characters_found, (unsigned int) c, 1);
					parse_state++;
					symbols[i] = c;
				}
			}
			/*
			 * Expected to read a delimiter or the end of the array
			 */
			else if (parse_state == -1)
			{
				if (c == ',')
				{
					if (i == n_elements - 1) {
						ereport(ERROR,(errmsg("multidimensional arrays must have array expressions with matching dimensions")));					}
					parse_state++;
				}
				else if (c == '}')
				{
					if (i != n_elements - 1) {
						ereport(ERROR,(errmsg("multidimensional arrays must have array expressions with matching dimensions")));
					}
					n_curlies_open--;
					parse_state++;
				}
				else
					ereport(ERROR,(errmsg("input array contains multi-character symbol")));
			}
		}
		parse_state -= 2;

		PB_DEBUG2(errmsg("Element #%d = %c", i, symbols[i]));
	}

	/*
	 * Parse second array if applicable.
	 *
	 * parse_state -4: expect separating comma
	 * parse_state -3: start state
	 * parse_state -2: expecting to read a symbol
	 * parse_state -1: symbol read, expecting a comma
	 * parse_state 0: comma read, or end of array
	 */
	if (PB_ALPHABET_TYPE(result) == PB_ALPHABET_WITH_PROBABILITIES)
	{
		float probability_sum = 0.0f;

		parse_state = -4;

		for (i = 0; i < n_elements; i++)
		{
			while (parse_state < 0)
			{
				char c = *input_pointer;
				input_pointer++;

				PB_DEBUG2(errmsg("Element #%d Character:%c parse_state:%d", i, c, parse_state));

				if (c == ' ' || c == '\t')
					continue;
				else if (c < 32 || c >= 127 || c == '\'' || c == '"')
					ereport(ERROR,(errmsg("input sequence violates alphabet constraints"),
							errdetail("Failing datum contains symbol \"%c\" (%d).", c, c)));

				if (parse_state == -4)
				{
					if (c != ',')
						ereport(ERROR,(errmsg("input array has wrong dimensions"),
									errdetail("Only one- and two-dimensional arrays are acceptable.")));
					else
						parse_state++;
				}
				else if (parse_state == -3)
				{
					if (c == '{')
					{
						n_curlies_open++;
						parse_state++;
					}
					else
						ereport(ERROR,(errmsg("input array has wrong dimensions"),
									errdetail("Only one- and two-dimensional arrays are acceptable.")));
				}
				else if (parse_state == -2)
				{
					if (c == ',' || c == '}')
						ereport(ERROR,(errmsg("input array contains null-value")));
					else if (c == '{')
						ereport(ERROR,(errmsg("multidimensional arrays must have array expressions with matching dimensions")));
					else
					{
						char* end;
						parse_state++;
						errno = 0;
						probabilities[i] = (PB_SymbolProbability) strtod(input_pointer - 1, &end);
						if (errno != 0)
							ereport(ERROR,(errmsg("probability array contains non-numeric value"),
									errdetail("Failing element contains non-numeric value at \"%s\"", (input_pointer - 1) )));
						input_pointer = end;
					}
				}
				else if (parse_state == -1)
				{
					if (c == ',')
					{
						if (i == n_elements - 1)
							ereport(ERROR,(errmsg("multidimensional arrays must have array expressions with matching dimensions")));
						parse_state++;
					}
					else if (c == '}')
					{
						if (i != n_elements - 1)
							ereport(ERROR,(errmsg("multidimensional arrays must have array expressions with matching dimensions")));
						parse_state++;
					}
					else
						ereport(ERROR,(errmsg("probabilities must be separated by commas")));
				}
			}
			parse_state -= 2;
		}

		for (i = 0; i < n_elements; i++)
		{
			probability_sum += probabilities[i];
			PB_DEBUG2(errmsg("Element #%d Symbol:%c Probability:%f", i, symbols[i], probabilities[i]));
		}

		if (probability_sum < 0.999999f || probability_sum > 1.000001f )
			ereport(ERROR,(errmsg("probabilities do not add up to 1 (%f)", probability_sum)));
	}

	pfree(characters_found);

	return result;
}

/**
 * alphabet_to_text()
 * 		Returns an array formatted text for a given alphabet.
 *
 * 	PB_Alphabet* input : alphabet to output
 */
text* alphabet_to_text(PB_Alphabet* input)
{
	text* result;

	char* float_buffer;
	char* output_pointer;
	uint8* probability_print_lengths = NULL;
	int i;

	int alphabet_size = PB_ALPHABET_SIZE(input);
	PB_SymbolProbability* probabilities = PB_ALPHABET_SYMBOL_PROBABILITY_POINTER(input);
	PB_Symbol* symbols = PB_ALPHABET_SYMBOL_POINTER(input);
	int size = 0;

	/*
	 * +2 curly brackets
	 * +1 symbol + 1 comma per symbol
	 * -1 comma for the last symbol
	 * +4 varlena header
	 */
	size = 2 + 2 * alphabet_size - 1 + 4;

	float_buffer = (char*) palloc0(32);

	if (PB_ALPHABET_TYPE(input) == PB_ALPHABET_WITH_PROBABILITIES)
	{
		probability_print_lengths = palloc(sizeof(uint8) * alphabet_size);
		/*
		 * +2 outer curly brackets
		 * +1 comma between symbols and probabilities
		 * +2 inner curly brackets
		 * +1 comma per probability
		 * -1 comma for the last probability
		 */
		size += 2 + 1 + 2 + alphabet_size -1;

		/*
		 * Count lengths of printed versions of the probabilities.
		 */
		for (i = 0; i < alphabet_size; i++)
		{
			uint8 print_length;
			sprintf((char*) float_buffer, "%f", probabilities[i]);
			print_length = (uint8) strlen(float_buffer);
			size += print_length;
			probability_print_lengths[i] = print_length;
		}
	}

	ereport(DEBUG1,(errmsg("Alphabet size as text:%d", size)));

	result = palloc(size);
	SET_VARSIZE(result,size);
	output_pointer = (char*) VARDATA(result);

	/*
	 * String construction
	 */
	*output_pointer = '{';
	output_pointer++;

	if (PB_ALPHABET_TYPE(input) == PB_ALPHABET_WITH_PROBABILITIES)
	{
		*output_pointer = '{';
		output_pointer++;
	}

	for (i = 0; i < alphabet_size - 1; i++)
	{
		sprintf(output_pointer, "%c,", symbols[i]);
		output_pointer += 2;
	}

	sprintf(output_pointer, "%c}", symbols[alphabet_size - 1]);
	output_pointer += 2;

	if (PB_ALPHABET_TYPE(input) == PB_ALPHABET_WITH_PROBABILITIES)
	{
		sprintf((char*)output_pointer, ",{");
		output_pointer += 2;

		for (i = 0; i < alphabet_size - 1; i++)
		{
			sprintf((char*)output_pointer, "%f,", probabilities[i]);
			output_pointer += probability_print_lengths[i] + 1;
		}

		sprintf((char*)output_pointer, "%f}}", probabilities[alphabet_size - 1]);
		output_pointer += probability_print_lengths[alphabet_size - 1] + 2;
	}

	pfree(float_buffer);

	if (PB_ALPHABET_TYPE(input) == PB_ALPHABET_WITH_PROBABILITIES)
		pfree(probability_print_lengths);

	return result;
}

/*
 *	pgsql interface functions
 */

/**
 * alphabet_in()
 * 		Creates an alphabet from a cstring.
 *
 * 	char* input : null-terminated input sequence (cstring)
 *
 */
PG_FUNCTION_INFO_V1 (alphabet_in);
Datum alphabet_in (PG_FUNCTION_ARGS)
{
	text* input_text = (text*) CStringGetTextDatum(PG_GETARG_CSTRING(0));
	PB_Alphabet* alphabet = parse_alphabet_from_text(input_text);

	PG_RETURN_POINTER(alphabet);
}

/**
 * alphabet_in_text()
 * 		Creates an alphabet from a text.
 *
 * 	text* input : text input sequence
 *
 */
PG_FUNCTION_INFO_V1 (alphabet_in_text);
Datum alphabet_in_text (PG_FUNCTION_ARGS)
{
	text* input_text = (text*) PG_GETARG_VARLENA_PP(0);
	PB_Alphabet* alphabet = parse_alphabet_from_text(input_text);

	PG_RETURN_POINTER(alphabet);
}

/**
 * alphabet_out()
 * 		Alphabet to cstring.
 *
 * 	PB_Alphabet* input : Alphabet
 */
PG_FUNCTION_INFO_V1 (alphabet_out);
Datum alphabet_out (PG_FUNCTION_ARGS) {
	PB_Alphabet* input =
			(PB_Alphabet*) PG_DETOAST_DATUM(PG_GETARG_POINTER(0));

	ereport(DEBUG1,(errmsg("as:%d", PB_ALPHABET_SIZE(input))));

	PG_RETURN_CSTRING(TextDatumGetCString(alphabet_to_text(input)));
}

/**
 * alphabet_out_text()
 * 		Alphabet to text.
 *
 * 	PB_Alphabet* input : Alphabet
 */
PG_FUNCTION_INFO_V1 (alphabet_out_text);
Datum alphabet_out_text (PG_FUNCTION_ARGS) {
	PB_Alphabet* input =
			(PB_Alphabet*) PG_DETOAST_DATUM(PG_GETARG_POINTER(0));

	ereport(DEBUG1,(errmsg("as:%d", PB_ALPHABET_SIZE(input))));

	PG_RETURN_POINTER(alphabet_to_text(input));
}

/**
 * Creates an alphabet from a text sequence.
 *
 * text* input : the input sequence
 */
PG_FUNCTION_INFO_V1(get_alphabet_text_sequence);
Datum get_alphabet_text_sequence (PG_FUNCTION_ARGS)
{
	text* input = (text*) PG_GETARG_VARLENA_PP(0);
	PB_Alphabet* result;
	PB_SequenceInfo* info;
	PB_Symbol* symbols;
	PB_SymbolProbability* probabilities;
	PB_SymbolProbability sequence_length;
	int i;

	info = get_sequence_info_text(input, PB_SEQUENCE_INFO_WITHOUT_RLE | PB_SEQUENCE_INFO_CASE_SENSITIVE);

	PB_ALPHABET_CREATE_WITH_PROBABILITIES(result, info->n_symbols);

	symbols = PB_ALPHABET_SYMBOL_POINTER(result);
	probabilities = PB_ALPHABET_SYMBOL_PROBABILITY_POINTER(result);
	sequence_length = (PB_SymbolProbability) info->sequence_length;

	for (i = 0; i < info->n_symbols; i++)
	{
		uint8 c = info->symbols[i];
		symbols[i] = c;
		probabilities[i] = ((PB_SymbolProbability) info->frequencies[c] / sequence_length);
		ereport(DEBUG2,(errmsg("%d %c:%d %f", i, c, info->frequencies[c], probabilities[i])));
	}

	PB_SEQUENCE_INFO_PFREE(info);

	PG_RETURN_POINTER(result);
}

/**
 * Creates an alphabet from compressed sequence.
 *
 * PB_CompressedSequence* input : the input sequence
 */
PB_Alphabet* get_alphabet_compressed_sequence (PB_CompressedSequence* input, PB_CodeSet** fixed_codesets)
{
	PB_Alphabet* result;
	PB_SequenceInfo* info;
	PB_Symbol* symbols;
	PB_SymbolProbability* probabilities;
	PB_SymbolProbability sequence_length;
	int i;
	uint8 c;

	info = palloc0(sizeof(PB_SequenceInfo));
	PB_BEGIN_DECODE((Varlena*) input, 0, input->sequence_length, fixed_codesets, c) {
		info->frequencies[c]++;
	} PB_END_DECODE

	collect_alphabet((uint32*) &(info->frequencies),&(info->n_symbols),&(info->symbols), NULL, NULL);

	PB_ALPHABET_CREATE_WITH_PROBABILITIES(result, info->n_symbols);

	symbols = PB_ALPHABET_SYMBOL_POINTER(result);
	probabilities = PB_ALPHABET_SYMBOL_PROBABILITY_POINTER(result);
	sequence_length = (PB_SymbolProbability) input->sequence_length;

	for (i = 0; i < info->n_symbols; i++)
	{
		uint8 c = info->symbols[i];
		symbols[i] = c;
		probabilities[i] = ((PB_SymbolProbability) info->frequencies[c] / sequence_length);
		ereport(DEBUG2,(errmsg("%d %c:%d %f", i, c, info->frequencies[c], probabilities[i])));
	}

	PB_SEQUENCE_INFO_PFREE(info);

	return result;
}


