/*-------------------------------------------------------------------------
*
* Copyright (c) 2013, Max Planck Institute for Marine Microbiology
*
* This software is released under the PostgreSQL License
*
* Author: Michael Schneider <mschneid@mpi-bremen.de>
*
* IDENTIFICATION
*   src/types/aa_sequence.c
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
#include "sequence/compression.h"
#include "sequence/functions.h"
#include "utils/debug.h"
#include "types/alphabet.h"

#include "types/aa_sequence.h"

/**
 * Section 1 - static fields
 */

static PB_CodeSet aa_iupac = {
		.n_symbols = (uint8) 23,
		.max_codeword_length = (uint8) 6,
		.n_swapped_symbols = (uint8) 0,
		.max_swapped_codeword_length = (uint8) 0,
		.has_equal_length = (bool) false,
		.is_fixed = (bool) true,
		.uses_rle = (bool) false,
		.ignore_case = (bool) true,
		.fixed_id = (bool) 0,
		.swap_savings = (uint64) 0,
		.ascii_bitmap_low = (uint64) 0,
		.ascii_bitmap_high = (uint64) 132086782,
		.words = {  {0,3,'L'}, {32,4,'A'}, {48,4,'G'}, {64,4,'S'},
				   {80,4,'V'}, {96,4,'E'},{112,4,'T'},{128,4,'K'},
				  {144,5,'X'},{152,5,'I'},{160,5,'P'},{168,5,'R'},
				  {176,5,'N'},{184,5,'Q'},{192,5,'F'},{200,5,'Y'},
				  {208,5,'M'},{216,5,'H'},{224,5,'D'},{232,5,'Z'},
				  {240,5,'B'},{248,6,'C'},{252,6,'W'}}
/*		.words = {  {0,4,'L'}, {16,4,'A'}, {32,4,'G'}, {48,4,'S'},
				   {64,4,'V'}, {80,4,'E'}, {96,4,'T'},{112,4,'K'},
				  {128,4,'X'},{144,5,'I'},{152,5,'P'},{160,5,'R'},
				  {168,5,'N'},{176,5,'Q'},{184,5,'F'},{192,5,'Y'},
				  {200,5,'M'},{208,5,'H'},{216,5,'C'},{224,5,'W'},
				  {232,5,'B'},{240,5,'D'},{248,5,'Z'}}*/
};

static PB_CodeSet aa_iupac_cs = {
		.n_symbols = (uint8) 46,
		.max_codeword_length = (uint8) 6,
		.n_swapped_symbols = (uint8) 0,
		.max_swapped_codeword_length = (uint8) 0,
		.has_equal_length = (bool) false,
		.is_fixed = (bool) true,
		.uses_rle = (bool) false,
		.ignore_case = (bool) false,
		.fixed_id = (bool) 1,
		.swap_savings = (uint64) 0,
		.ascii_bitmap_low = (uint64) 0,
		.ascii_bitmap_high = (uint64) 567308409055968254,
		.words = {  {0,5,'L'}, {16,5,'A'}, {32,5,'G'}, {48,5,'S'},
				   {64,5,'V'}, {80,5,'E'}, {96,5,'T'},{112,5,'K'},
				  {128,5,'X'},{144,6,'I'},{152,6,'P'},{160,6,'R'},
				  {168,6,'N'},{176,6,'Q'},{184,6,'F'},{192,6,'Y'},
				  {200,6,'M'},{208,6,'H'},{216,6,'C'},{224,6,'W'},
				  {232,6,'B'},{240,6,'D'},{248,6,'Z'},
				    {8,5,'l'}, {24,5,'a'}, {40,5,'g'}, {56,5,'s'},
				   {72,5,'v'}, {88,5,'e'},{104,5,'t'},{120,5,'k'},
				  {136,5,'x'},{148,6,'i'},{156,6,'p'},{164,6,'r'},
				  {172,6,'n'},{180,6,'q'},{188,6,'f'},{196,6,'y'},
				  {204,6,'m'},{212,6,'h'},{220,6,'c'},{228,6,'w'},
				  {236,6,'b'},{244,6,'d'},{252,6,'z'}
		}
};

static unsigned int n_fixed_aa_codes = 2;
static PB_CodeSet* fixed_aa_codes[] = {
		&aa_iupac,
		&aa_iupac_cs
};

/**
 * get_fixed_aa_codes()
 * 		Returns pointer to fixed AA codes.
 */
PB_CodeSet** get_fixed_aa_codes(void)
{
	return fixed_aa_codes;
}

PB_AaSequenceTypMod non_restricting_aa_typmod = {
	.case_sensitive = PB_AA_TYPMOD_CASE_SENSITIVE,
	.restricting_alphabet = PB_AA_TYPMOD_ASCII
};

/*
 * Section 2 - other public functions
 */

/**
 * aa_sequence_typmod_to_int()
 * 		Convert from PB_AaSequenceTypMod to int
 *
 * 	PB_AaSequenceTypMod typmod : type modifier
 */
int aa_sequence_typmod_to_int(PB_AaSequenceTypMod typmod)
{
	return *((int*) &typmod);
}

/**
 * int_to_aa_sequence_typmod()
 * 		Convert from int to PB_AaSequenceTypMod
 *
 * 	int typmod : type modifier
 */
PB_AaSequenceTypMod int_to_aa_sequence_typmod(int typmod)
{
	if (-1 == typmod) {
		typmod = 0;
	}
	return *((PB_AaSequenceTypMod*) &typmod);
}

/**
 * get_fixed_aa_code()
 * 		Returns a fixed code for the specified id.
 *
 * 	Id	|	Code Description
 * -----------------------------------------------------------
 * 	0	|	AA IUPAC code
 * 	1	|	AA IUPAC code, case sensitive
 *
 */
PB_CodeSet* get_fixed_aa_code(unsigned int fixed_code_id)
{
	PB_CodeSet* result = NULL;

	if (fixed_code_id < n_fixed_aa_codes)
		result = fixed_aa_codes[fixed_code_id];

	return result;
}

/**
 * compress_aa_sequence()
 * 		Compress an AA sequence.
 *
 *	First it will be checked whether the given sequence matches
 *	the restricting alphabet specified with the type modifiers.
 *	Then one of three types of compression will be performed,
 *	depending on which one is most efficient.
 *
 *	A) Simple Huffman-Codes
 *	B) Equal lengths codes
 *	C) Fixed pre-built codes
 *
 * 	uint8* input : unterminated or null-terminated input sequence
 * 	PB_AaSequenceTypMod typmod : target type modifier
 * 	PB_SequenceInfo info : given by input function
 */
PB_CompressedSequence* compress_aa_sequence(uint8* input,
											PB_AaSequenceTypMod typmod,
											PB_SequenceInfo* info)
{
	PB_CompressedSequence* result;
	PB_CodeSet* codeset = NULL;
	uint32 compressed_size = 0xFFFFFFFF;

	PB_TRACE(errmsg("->compress_aa_sequence()"));

	PB_DEBUG1(errmsg("compress_aa_sequence(): low bitmap:%ld high bitmap:%ld",  info->ascii_bitmap_low, info->ascii_bitmap_high));

	/*
	 * Check alphabet constraints.
	 */
	if (typmod.restricting_alphabet == PB_AA_TYPMOD_IUPAC  && (!PB_CHECK_CODESET((&aa_iupac_cs),info)))
	{
		ereport(ERROR,(errmsg("input sequence violates alphabet restrictions")));
	}

	/*
	 * Build sequence specific code set if:
	 * 		-sequence is longer than 512 chars or
	 * 		-sequence can not be expressed with a pre-built codeset
	 */
	if (info->sequence_length > 512 || !PB_CHECK_CODESET((&aa_iupac_cs),info)) {
		codeset = get_huffman_code(info);

		if (!codeset)
			codeset = get_equal_lengths_code(info);

		compressed_size = get_compressed_size(info, codeset);
	}

	/*
	 * Replace sequence specific code set with pre-built code set
	 * if advantageous.
	 */
	if (PB_CHECK_CODESET((&aa_iupac_cs),info))
	{
		PB_CodeSet* fixed_codeset = &aa_iupac_cs;
		uint32 compressed_size_fixed;

		if (PB_CHECK_CODESET((&aa_iupac),info))
			fixed_codeset = &aa_iupac;

		compressed_size_fixed = get_compressed_size(info, fixed_codeset);

		if (compressed_size > compressed_size_fixed)
		{
			if (codeset)
				pfree(codeset);

			codeset = fixed_codeset;
			compressed_size = compressed_size_fixed;
		}
	}

	/*
	 * Compress.
	 */
	result = encode(input, compressed_size, codeset, info);

	if (!codeset->is_fixed)
		pfree(codeset);

	return result;
}

/**
 * decompress_aa_sequence()
 * 		Decompress an AA sequence
 *
 * 	PB_CompressedSequence* input : compressed input sequence (takes unTOASTed data)
 * 	uint8* output : already allocated appropriately sized target memory area
 * 	uint32 from_position : position to start decompressing, first = 0
 * 	uint32 length : length of sequence to decompress
 */
void decompress_aa_sequence(PB_CompressedSequence* input,
							uint8* output,
							uint32 from_position,
							uint32 length)
{
	decode((Varlena*) input, output, from_position, length, fixed_aa_codes);
}

/*
 *	Section 3 - pgsql interface functions
 */


/**
 *	aa_sequence_typmod_in()
 *		Condense type modifier keywords into single integer value.
 *
 *	cstring[] input : lower-case keywords separated into array
 */
PG_FUNCTION_INFO_V1 (aa_sequence_typmod_in);
Datum aa_sequence_typmod_in(PG_FUNCTION_ARGS) {
	ArrayType* input = PG_GETARG_ARRAYTYPE_P(0);
	char* read_pointer = ((char*) input) + ARR_DATA_OFFSET(input);

	PB_AaSequenceTypMod result;

	bool typeModCaseInsensitive = false;
	bool typeModCaseSensitive = false;
	bool typeModIupac = false;
	bool typeModAscii = false;

	int i;

	PB_TRACE(errmsg("->aa_sequence_typmod_in()"));

	/*
	 * Parsing
	 */
	for (i = ARR_DIMS(input)[0] - 1; i >= 0; i--) {
		PB_DEBUG2(errmsg("aa_sequence_typmod_in(): %d= \"%s\".", i, read_pointer));
		if (!strcmp(read_pointer, "case_insensitive")) {
			typeModCaseInsensitive = true;
		} else if (!strcmp(read_pointer, "case_sensitive")) {
			typeModCaseSensitive = true;
		} else if (!strcmp(read_pointer, "iupac")) {
			typeModIupac = true;
		} else if (!strcmp(read_pointer, "ascii")) {
			typeModAscii = true;
		} else {
			ereport(ERROR,(errmsg("type modifier invalid"),
					errdetail("Can not recognize type modifier \"%s\".", read_pointer)));
		}
		read_pointer += strlen(read_pointer) + 1;
	}

	/*
	 * Syntax checking.
	 */
	if (typeModCaseInsensitive && typeModCaseSensitive) {
		ereport(ERROR,(errmsg("CASE_INSENSITIVE and CASE_SENSITIVE are mutually exclusive type modifiers")));
	}

	if ((int) typeModIupac + (int) typeModAscii > 1) {
		ereport(ERROR,(errmsg("IUPAC, and ASCII are mutually exclusive type modifiers")));
	}

	/*
	 * Build integer value from parsed type modifiers.
	 */
	if (typeModCaseSensitive) {
		result.case_sensitive = PB_AA_TYPMOD_CASE_SENSITIVE;
	} else {
		result.case_sensitive = PB_AA_TYPMOD_CASE_INSENSITIVE;
	}

	if (typeModAscii) {
		result.restricting_alphabet = PB_AA_TYPMOD_ASCII;
	} else {
		result.restricting_alphabet = PB_AA_TYPMOD_IUPAC;
	}

	PB_TRACE(errmsg("<-aa_sequence_typmod_in() returning %d", aa_sequence_typmod_to_int(result)));

	PG_RETURN_INT32(aa_sequence_typmod_to_int(result));
}

/**
 * 	aa_sequence_typmod_out()
 * 		Restore type modifier keywords from single integer value.
 *
 * 	int typmod : single value representing the type modifiers
 */
PG_FUNCTION_INFO_V1 (aa_sequence_typmod_out);
Datum aa_sequence_typmod_out(PG_FUNCTION_ARGS) {
	int input = PG_GETARG_INT32(0);
	PB_AaSequenceTypMod typmod = int_to_aa_sequence_typmod(input);
	int len = 3; /* 2 parentheses + null-terminator */
	char* result;
	char* out;

	/*
	 * Compute length of result.
	 */
	if (typmod.case_sensitive == PB_AA_TYPMOD_CASE_SENSITIVE) {
		len += 15; /* strlen('CASE_SENSITIVE,') = 15 */
	} else {
		len += 17; /* strlen('CASE_INSENSITIVE,') = 17 */
	}
	len += 5; /* strlen('ASCII') = 5, strlen('IUPAC') = 5 */

	result = palloc0(len);
	out = result;

	*out = '(';
	out++;

	if (typmod.case_sensitive == PB_AA_TYPMOD_CASE_SENSITIVE) {
		strcpy(out, "CASE_SENSITIVE,");
		out+=15;
	} else {
		strcpy(out, "CASE_INSENSITIVE,");
		out+=17;
	}

	if (typmod.restricting_alphabet == PB_AA_TYPMOD_ASCII) {
		strcpy(out, "ASCII");
		out+=5;
	} else {
		strcpy(out, "IUPAC");
		out+=5;
	}

	*out = ')';
	out++;
	*out = 0;

	PB_TRACE(errmsg("<-aa_sequence_typmod_out() returning '%s'", out));

	PG_RETURN_INT32(result);
}

/**
 * aa_sequence_in()
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
PG_FUNCTION_INFO_V1 (aa_sequence_in);
Datum aa_sequence_in (PG_FUNCTION_ARGS)
{
	uint8* input = (uint8*) PG_GETARG_CSTRING(0);
	//Oid oid = (Oid) PG_GETARG_OID(1);
	int32 typmod_int = PG_GETARG_INT32(2);
	PB_AaSequenceTypMod typmod;

	PB_SequenceInfo* info;
	PB_CompressedSequence* result;

	int mode;

	PB_TRACE(errmsg("->aa_sequence_in()"));

	if ((-1) == typmod_int)
		typmod = non_restricting_aa_typmod;
	else
		 typmod = int_to_aa_sequence_typmod(typmod_int);

	/*
	 * Determine sequence info collection mode.
	 */
	mode = 0;
	if (typmod.case_sensitive == PB_AA_TYPMOD_CASE_SENSITIVE)
		mode |= PB_SEQUENCE_INFO_CASE_SENSITIVE;

	info = get_sequence_info_cstring(input, mode);

	result = compress_aa_sequence(input, typmod, info);

	PB_SEQUENCE_INFO_PFREE(info);

	PB_TRACE(errmsg("<-aa_sequence_in()"));

	PG_RETURN_POINTER(result);
}

/**
 * aa_sequence_in_varlena()
 * 		Compress a given input sequence.
 *
 * 	This function expects a varlena input sequence, that is
 * 	text, varchar or char. It will be called by the respective
 * 	cast functions.
 *
 * 	varlena* input : input sequence
 * 	int typmod : single value representing target type modifier
 */
PG_FUNCTION_INFO_V1 (aa_sequence_in_varlena);
Datum aa_sequence_in_varlena (PG_FUNCTION_ARGS)
{
	Varlena* input = (Varlena*) PG_GETARG_VARLENA_P(0);
	PB_AaSequenceTypMod typmod = int_to_aa_sequence_typmod(PG_GETARG_INT32(1));

	PB_SequenceInfo* info;
	PB_CompressedSequence* result;

	int mode;

	PB_TRACE(errmsg("->aa_sequence_in_varlena()"));

	/*
	 * Determine sequence info collection mode.
	 */
	mode = 0;
	if (typmod.case_sensitive == PB_AA_TYPMOD_CASE_SENSITIVE)
		mode |= PB_SEQUENCE_INFO_CASE_SENSITIVE;

	info = get_sequence_info_text(input, mode);

	result = compress_aa_sequence((uint8*) VARDATA(input), typmod, info);

	PB_SEQUENCE_INFO_PFREE(info);

	PB_TRACE(errmsg("<-aa_sequence_in_varlena()"));

	PG_RETURN_POINTER(result);
}

/**
 * aa_sequence_cast()
 * 		Decompress a given sequence and compress it again
 * 		using a different compression
 *
 * 	varlena* input : input sequence
 * 	int typmod : single value representing target type modifier
 */
PG_FUNCTION_INFO_V1 (aa_sequence_cast);
Datum aa_sequence_cast (PG_FUNCTION_ARGS)
{
	PB_CompressedSequence* input = (PB_CompressedSequence*) PG_GETARG_VARLENA_P(0);
	PB_AaSequenceTypMod typmod = int_to_aa_sequence_typmod(PG_GETARG_INT32(1));

	uint8* plain;
	PB_SequenceInfo* info;
	PB_CompressedSequence* result;

	int mode ;

	PB_TRACE(errmsg("->aa_sequence_cast()"));

	/*
	 * Determine sequence info collection mode.
	 */
	mode = 0;
	if (typmod.case_sensitive == PB_AA_TYPMOD_CASE_SENSITIVE)
		mode |= PB_SEQUENCE_INFO_CASE_SENSITIVE;

	/*
	 * Decompress sequence.
	 */
	plain = palloc0(input->sequence_length + 1);
	decompress_aa_sequence(input,plain,0,input->sequence_length);
	plain[input->sequence_length] = '\0';

	/*
	 * Compress again.
	 */
	info = get_sequence_info_cstring(plain, mode);
	result = compress_aa_sequence(plain, typmod, info);

	PB_SEQUENCE_INFO_PFREE(info);

	pfree(plain);

	PB_TRACE(errmsg("<-aa_sequence_cast()"));

	PG_RETURN_POINTER(result);
}

/**
 * aa_sequence_out()
 * 		Decompress a sequence.
 *
 * 	CompressedSequence* input : compressed input sequence
 */
PG_FUNCTION_INFO_V1 (aa_sequence_out);
Datum aa_sequence_out (PG_FUNCTION_ARGS)
{
	PB_CompressedSequence* input = (PB_CompressedSequence*) PG_GETARG_VARLENA_P(0);
	uint8* result;

	PB_TRACE(errmsg("->aa_sequence_out()"));

	result = palloc0(input->sequence_length + 1);
	decompress_aa_sequence(input, result, 0, input->sequence_length);
	result[input->sequence_length] = '\000';

	PB_TRACE(errmsg("<-aa_sequence_out()"));

	PG_RETURN_CSTRING ((char*)result);
}

/**
 * aa_sequence_out_varlena()
 * 		Decompress a sequence into varlena.
 *
 * 	CompressedSequence* input : compressed input sequence
 */
PG_FUNCTION_INFO_V1 (aa_sequence_out_varlena);
Datum aa_sequence_out_varlena (PG_FUNCTION_ARGS)
{
	PB_CompressedSequence* input = (PB_CompressedSequence*) PG_GETARG_VARLENA_P(0);
	text* result;

	PB_TRACE(errmsg("->aa_sequence_out_varlena()"));

	result = palloc0(input->sequence_length + VARHDRSZ);
	SET_VARSIZE (result, input->sequence_length + VARHDRSZ);
	decompress_aa_sequence(input, (uint8*) VARDATA_ANY(result), 0, input->sequence_length);

	PB_TRACE(errmsg("<-aa_sequence_out_varlena()"));

	PG_RETURN_POINTER(result);
}

/**
 * aa_sequence_substring()
 * 		Decompress a substring of a sequence.
 *
 * 	This function mimics the originals substr function's
 * 	behaviour. The first position is 1.
 *
 * 	PB_CompressedSequence* input : compressed input sequence
 * 	int start : position to start from
 * 	int len : length of substring
 */
PG_FUNCTION_INFO_V1 (aa_sequence_substring);
Datum aa_sequence_substring (PG_FUNCTION_ARGS) {
	Varlena* input = (Varlena*) PG_GETARG_RAW_VARLENA_P(0);
	int start = PG_GETARG_DATUM(1);
	int len = PG_GETARG_DATUM(2);

	PB_CompressedSequence* input_header;

	text* result;

	PB_TRACE(errmsg("->aa_sequence_substring()"));

	input_header = (PB_CompressedSequence*)
			PG_DETOAST_DATUM_SLICE(input, 0, sizeof(PB_CompressedSequence) - VARHDRSZ);

	if (len < 0) {
		ereport(ERROR,(errmsg("negative substring length not allowed")));
	}

	/*
	 * SQL's first position is 1, our first position is 0
	 */
	start--;
	if (start < 0) {
		len += start;
		start = 0;
	}
	if (start >= input_header->sequence_length || len < 1) {
		result = palloc0(4);
		SET_VARSIZE (result, 4);
		PG_RETURN_POINTER(result);
	}
	if (start + len > input_header->sequence_length) {
		len = input_header->sequence_length - start;
	}

	result = palloc0(len + VARHDRSZ);
	SET_VARSIZE (result, len + VARHDRSZ);
	decompress_aa_sequence((PB_CompressedSequence*) input,(uint8*) VARDATA(result),start,len);

	pfree(input_header);

	PB_TRACE(errmsg("<-aa_sequence_substring()"));

	PG_RETURN_POINTER(result);
}

/**
 * aa_sequence_char_length()
 * 		Get length of sequence.
 *
 * 	PB_CompressedSequence* input : compressed input sequence
 */
PG_FUNCTION_INFO_V1 (aa_sequence_char_length);
Datum aa_sequence_char_length (PG_FUNCTION_ARGS)
{
	PB_CompressedSequence* seq = (PB_CompressedSequence*)
			PG_DETOAST_DATUM_SLICE(PG_GETARG_RAW_VARLENA_P(0),0,4);

	PG_RETURN_INT32(seq->sequence_length);
}

/**
 * aa_sequence_compression_ratio()
 * 		Get compression ratio.
 *
 * 	The ratio between the size of the sequence as pgsql text
 * 	type and the size of the compressed sequence including all
 * 	required meta-data, such as the substring-index.
 *
 * 	CompressedSequence* input : compressed input sequence
 */
PG_FUNCTION_INFO_V1 (aa_sequence_compression_ratio);
Datum aa_sequence_compression_ratio (PG_FUNCTION_ARGS)
{
	Varlena* input = (Varlena*) PG_GETARG_RAW_VARLENA_P(0);
	PB_CompressedSequence* input_header = (PB_CompressedSequence*)
		PG_DETOAST_DATUM_SLICE(input, 0, sizeof(PB_CompressedSequence) - VARHDRSZ);

	double cr = (double) toast_raw_datum_size((Datum)input) / (double) (input_header->sequence_length + VARHDRSZ);

	PB_DEBUG1(errmsg("aa_sequence_compression_ratio(): cr=%f memsize=%u len=%u", cr, (uint32) toast_raw_datum_size((Datum)input), input_header->sequence_length + VARHDRSZ));

	PG_RETURN_FLOAT8(cr);
}

/**
 * aa_sequence_reverse()
 * 		Returns the reverse of an aa sequence.
 *
 * 	PB_CompressedSequence* input : compressed input sequence
 */
PG_FUNCTION_INFO_V1 (aa_sequence_reverse);
Datum aa_sequence_reverse(PG_FUNCTION_ARGS)
{
	PB_CompressedSequence* input =
			(PB_CompressedSequence*) PG_GETARG_VARLENA_P(0);

	PB_CompressedSequence* result;

	result = reverse(input, fixed_aa_codes);

	PG_RETURN_POINTER(result);
}

/**
 * get_alphabet_aa_sequence()
 * 		Calculates alphabet from aa sequence.
 *
 * 	PB_CompressedSequence* input : compressed input sequence
 */
PG_FUNCTION_INFO_V1 (get_alphabet_aa_sequence);
Datum get_alphabet_aa_sequence(PG_FUNCTION_ARGS)
{
	PB_CompressedSequence* input =
			(PB_CompressedSequence*) PG_GETARG_VARLENA_P(0);
	PB_Alphabet* result = get_alphabet_compressed_sequence(input, fixed_aa_codes);

	PG_RETURN_POINTER(result);
}

/**
 * equal_aa()
 * 		Compares two AA sequence for equality
 *
 * 	Varlena* seq1 : first sequence
 * 	Varlena* seq2 : second sequence
 */
PG_FUNCTION_INFO_V1 (equal_aa);
Datum equal_aa(PG_FUNCTION_ARGS)
{
	Varlena* seq1 = (Varlena*) PG_GETARG_RAW_VARLENA_P(0);
	Varlena* seq2 = (Varlena*) PG_GETARG_RAW_VARLENA_P(1);
	bool result;

	PB_TRACE(errmsg("->equal_aa()"));

	result = sequence_equal(seq1, seq2, fixed_aa_codes);

	PB_TRACE(errmsg("<-equal_aa() exists with %d", result));

	PG_RETURN_BOOL(result);
}

/**
 * compare_aa_lt()
 * 		Compares two AA sequence for equality Less-than.
 *
 * 	Varlena* seq1 : first sequence
 * 	Varlena* seq2 : second sequence
 */
PG_FUNCTION_INFO_V1 (compare_aa_lt);
Datum compare_aa_lt(PG_FUNCTION_ARGS)
{
	Varlena* seq1 = (Varlena*) PG_GETARG_RAW_VARLENA_P(0);
	Varlena* seq2 = (Varlena*) PG_GETARG_RAW_VARLENA_P(1);
	bool result = false;

	PB_TRACE(errmsg("->compare_aa_lt()"));

	if (sequence_compare(seq1, seq2, fixed_aa_codes) < 0)
		result = true;

	PB_TRACE(errmsg("<-compare_aa_lt() exists with %d", result));

	PG_RETURN_BOOL(result);
}

/**
 * compare_aa_le()
 * 		Compares two AA sequence for equality Less or equal.
 *
 * 	Varlena* seq1 : first sequence
 * 	Varlena* seq2 : second sequence
 */
PG_FUNCTION_INFO_V1 (compare_aa_le);
Datum compare_aa_le(PG_FUNCTION_ARGS)
{
	Varlena* seq1 = (Varlena*) PG_GETARG_RAW_VARLENA_P(0);
	Varlena* seq2 = (Varlena*) PG_GETARG_RAW_VARLENA_P(1);
	bool result = false;

	PB_TRACE(errmsg("->compare_aa_le()"));

	if (sequence_compare(seq1, seq2, fixed_aa_codes) <= 0)
		result = true;

	PB_TRACE(errmsg("<-compare_aa_le() exists with %d", result));

	PG_RETURN_BOOL(result);
}

/**
 * compare_aa_gt()
 * 		Compares two AA sequence for equality. greater than
 *
 * 	Varlena* seq1 : first sequence
 * 	Varlena* seq2 : second sequence
 */
PG_FUNCTION_INFO_V1 (compare_aa_gt);
Datum compare_aa_gt(PG_FUNCTION_ARGS)
{
	Varlena* seq1 = (Varlena*) PG_GETARG_RAW_VARLENA_P(0);
	Varlena* seq2 = (Varlena*) PG_GETARG_RAW_VARLENA_P(1);
	bool result = false;

	PB_TRACE(errmsg("->compare_aa_gt()"));

	if (sequence_compare(seq1, seq2, fixed_aa_codes) > 0)
		result = true;

	PB_TRACE(errmsg("<-compare_aa_gt() exists with %d", result));

	PG_RETURN_BOOL(result);
}

/**
 * compare_aa_ge()
 * 		Compares two AA sequence for equality. greater or equal
 *
 * 	Varlena* seq1 : first sequence
 * 	Varlena* seq2 : second sequence
 */
PG_FUNCTION_INFO_V1 (compare_aa_ge);
Datum compare_aa_ge(PG_FUNCTION_ARGS)
{
	Varlena* seq1 = (Varlena*) PG_GETARG_RAW_VARLENA_P(0);
	Varlena* seq2 = (Varlena*) PG_GETARG_RAW_VARLENA_P(1);
	bool result = false;

	PB_TRACE(errmsg("->compare_aa_ge()"));

	if (sequence_compare(seq1, seq2, fixed_aa_codes) >= 0)
		result = true;

	PB_TRACE(errmsg("<-compare_aa_ge() exists with %d", result));

	PG_RETURN_BOOL(result);
}

/**
 * compare_aa()
 * 		Compares two AA sequence for equality.
 *
 * 	Varlena* seq1 : first possibly toasted sequence
 * 	Varlena* seq2 : second possibly toasted sequence
 */
PG_FUNCTION_INFO_V1 (compare_aa);
Datum compare_aa(PG_FUNCTION_ARGS)
{
	Varlena* seq1 = (Varlena*) PG_GETARG_RAW_VARLENA_P(0);
	Varlena* seq2 = (Varlena*) PG_GETARG_RAW_VARLENA_P(1);
	int result;

	PB_TRACE(errmsg("->compare_aa()"));

	result = sequence_compare(seq1, seq2, fixed_aa_codes);

	PB_TRACE(errmsg("<-compare_aa() exists with %d", result));

	PG_RETURN_INT32(result);
}

/**
 * hash_aa()
 * 		Returns a CRC32 for a AA sequence.
 *
 * 	PB_CompressedSequence* seq1 : input sequence
 */
PG_FUNCTION_INFO_V1 (hash_aa);
Datum hash_aa(PG_FUNCTION_ARGS)
{
	PB_CompressedSequence* seq1 = (PB_CompressedSequence*) PG_GETARG_VARLENA_P(0);
	uint32 result;

	PB_TRACE(errmsg("->hash_aa()"));

	result = sequence_crc32(seq1, fixed_aa_codes);

	PB_TRACE(errmsg("<-hash_aa() exits with %u", result));

	PG_RETURN_UINT32(result);
}

/**
 * strpos_aa()
 * 		Finds the first occurrence of a pattern in a sequence.
 */
PG_FUNCTION_INFO_V1 (strpos_aa);
Datum strpos_aa(PG_FUNCTION_ARGS)
{
	PB_CompressedSequence* seq = (PB_CompressedSequence*) PG_GETARG_VARLENA_P(0);
	text* search = (text*) PG_GETARG_VARLENA_P(1);
	uint32 result;

	PB_TRACE(errmsg("->strpos_aa()"));

	result = sequence_strpos(seq, search, fixed_aa_codes);

	PB_TRACE(errmsg("<-strpos_aa()"));

	PG_RETURN_UINT32(result);
}

/**
 * octet_length_aa()
 * 		Returns byte size of datum.
 */
PG_FUNCTION_INFO_V1 (octet_length_aa);
Datum octet_length_aa(PG_FUNCTION_ARGS)
{
	PG_RETURN_UINT32((uint32) toast_raw_datum_size((Datum) PG_GETARG_RAW_VARLENA_P(0)));
}
