/*-------------------------------------------------------------------------
*
* Copyright (c) 2013, Max Planck Institute for Marine Microbiology
*
* This software is released under the PostgreSQL License
*
* Author: Michael Schneider <mschneid@mpi-bremen.de>
*
* IDENTIFICATION
*   src/types/aligned_aa_sequence.c
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

#include "types/aligned_aa_sequence.h"

/**
 * Section 1 - static fields
 */

static PB_CodeSet aligned_aa_iupac = {
		.n_symbols = (uint8) 25,
		.max_codeword_length = (uint8) 5,
		.n_swapped_symbols = (uint8) 0,
		.max_swapped_codeword_length = (uint8) 0,
		.has_equal_length = (bool) FALSE,
		.is_fixed = (bool) TRUE,
		.uses_rle = (bool) FALSE,
		.ignore_case = (bool) TRUE,
		.fixed_id = (bool) 0,
		.swap_savings = (uint64) 0,
		.ascii_bitmap_low = (uint64) 105553116266496,
		.ascii_bitmap_high = (uint64) 132086782,
		.words = {  {0,2,'-'},
				   {64,5,'L'}, {72,5,'A'}, {80,5,'G'}, {88,5,'S'},
				   {96,5,'V'},{104,5,'E'},{112,5,'T'},{120,5,'K'},
				  {128,5,'D'},{136,5,'I'},{144,5,'P'},{152,5,'R'},
				  {160,5,'N'},{168,5,'Q'},{176,5,'F'},{184,5,'Y'},
				  {192,5,'M'},{200,5,'H'},{208,5,'C'},{216,5,'W'},
				  {224,5,'B'},{232,5,'Z'},{240,5,'X'},{248,5,'.'}}
};

static PB_CodeSet aligned_aa_iupac_cs = {
		.n_symbols = (uint8) 48,
		.max_codeword_length = (uint8) 6,
		.n_swapped_symbols = (uint8) 0,
		.max_swapped_codeword_length = (uint8) 0,
		.has_equal_length = (bool) FALSE,
		.is_fixed = (bool) TRUE,
		.uses_rle = (bool) FALSE,
		.ignore_case = (bool) FALSE,
		.fixed_id = (bool) 1,
		.swap_savings = (uint64) 0,
		.ascii_bitmap_low = (uint64) 105553116266496,
		.ascii_bitmap_high = (uint64) 567308409055968254,
		.words = {  {0,2,'-'},
				   {64,6,'L'}, {72,6,'A'}, {80,6,'G'}, {88,6,'S'},
				   {96,6,'V'},{104,6,'E'},{112,6,'T'},{120,6,'K'},
				  {128,6,'D'},{136,6,'I'},{144,6,'P'},{152,6,'R'},
				  {160,6,'N'},{168,6,'Q'},{176,6,'F'},{184,6,'Y'},
				  {192,6,'M'},{200,6,'H'},{208,6,'C'},{216,6,'W'},
				  {224,6,'B'},{232,6,'Z'},{240,6,'X'},{248,5,'.'},
				   {68,6,'l'}, {76,6,'a'}, {84,6,'g'}, {92,6,'s'},
				  {100,6,'v'},{108,6,'e'},{116,6,'t'},{124,6,'k'},
				  {132,6,'d'},{140,6,'i'},{148,6,'p'},{155,6,'r'},
				  {164,6,'n'},{172,6,'q'},{180,6,'f'},{188,6,'y'},
				  {196,6,'m'},{204,6,'h'},{212,6,'c'},{220,6,'w'},
				  {228,6,'b'},{236,6,'z'},{244,6,'x'}}
};

static unsigned int n_fixed_aligned_aa_codes = 2;
static PB_CodeSet* fixed_aligned_aa_codes[] = {
		&aligned_aa_iupac,
		&aligned_aa_iupac_cs
};

PB_AlignedAaSequenceTypMod non_restricting_aligned_aa_typmod = {
	.case_sensitive = PB_ALIGNED_AA_TYPMOD_CASE_SENSITIVE,
	.restricting_alphabet = PB_ALIGNED_AA_TYPMOD_ASCII
};

/*
 * Section 2 - other public functions
 */

/**
 * aligned_aa_sequence_typmod_to_int()
 * 		Convert from PB_AlignedAaSequenceTypMod to int
 *
 * 	PB_AlignedAaSequenceTypMod typmod : type modifier
 */
int aligned_aa_sequence_typmod_to_int(PB_AlignedAaSequenceTypMod typmod)
{
	return *((int*) &typmod);
}

/**
 * int_to_aligned_aa_sequence_typmod()
 * 		Convert from int to PB_AlignedAaSequenceTypMod
 *
 * 	int typmod : type modifier
 */
PB_AlignedAaSequenceTypMod int_to_aligned_aa_sequence_typmod(int typmod)
{
	if (-1 == typmod) {
		typmod = 0;
	}
	return *((PB_AlignedAaSequenceTypMod*) &typmod);
}

/**
 * get_fixed_aligned_aa_code()
 * 		Returns a fixed code for the specified id.
 *
 * 	Id	|	Code Description
 * -----------------------------------------------------------
 * 	0	|	AA IUPAC code
 * 	1	|	AA IUPAC code, case sensitive
 *
 */
PB_CodeSet* get_fixed_aligned_aa_code(unsigned int fixed_code_id)
{
	PB_CodeSet* result = NULL;

	if (fixed_code_id < n_fixed_aligned_aa_codes)
		result = fixed_aligned_aa_codes[fixed_code_id];

	return result;
}

/**
 * get_fixed_aligned_aa_codes()
 * 		Returns pointer to fixed aligned AA codes.
 */
PB_CodeSet** get_fixed_aligned_aa_codes(void);

/**
 * compress_aligned_aa_sequence()
 * 		Compress an aligned AA sequence.
 *
 *	First it will be checked whether the given sequence matches
 *	the restricting alphabet specified with the type modifiers.
 *	Then the sequence will be optimally compressed.
 *
 * 	uint8* input : unterminated or null-terminated input sequence
 * 	PB_RnaSequenceTypMod typmod : target type modifier
 * 	PB_SequenceInfo info : given by input function
 */
PB_CompressedSequence* compress_aligned_aa_sequence(uint8* input,
													PB_AlignedAaSequenceTypMod typmod,
													PB_SequenceInfo* info)
{
	PB_CompressedSequence* result;
	PB_CodeSet* codeset;
	uint32 compressed_size;

	PB_TRACE(errmsg("->compress_aligned_aa_sequence()"));

	PB_DEBUG1(errmsg("compress_aligned_aa_sequence(): low bitmap:%ld high bitmap:%ld",  info->ascii_bitmap_low, info->ascii_bitmap_high));

	/*
	 * Check alphabet constraints.
	 */
	if (typmod.restricting_alphabet == PB_ALIGNED_AA_TYPMOD_IUPAC  && (!PB_CHECK_CODESET((&aligned_aa_iupac_cs),info)))
	{
		ereport(ERROR,(errmsg("input sequence violates alphabet restrictions")));
	}

	/*
	 * Build sequence specific code set.
	 */
	codeset = get_optimal_code(info);

	compressed_size = get_compressed_size(info, codeset);

	/*
	 * Replace sequence specific code set with pre-built code set
	 * if advantageous.
	 */
	if (PB_CHECK_CODESET((&aligned_aa_iupac_cs),info))
	{
		PB_CodeSet* fixed_codeset = &aligned_aa_iupac_cs;
		uint32 compressed_size_fixed;

		if (PB_CHECK_CODESET((&aligned_aa_iupac),info))
			fixed_codeset = &aligned_aa_iupac;

		compressed_size_fixed = get_compressed_size(info, fixed_codeset);

		if (compressed_size > compressed_size_fixed)
		{
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
 * decompress_aligned_aa_sequence()
 * 		Decompress an aligned AA sequence
 *
 * 	PB_CompressedSequence* input : compressed input sequence (takes unTOASTed data)
 * 	uint8* output : already allocated appropriately sized target memory area
 * 	uint32 from_position : position to start decompressing, first = 0
 * 	uint32 length : length of sequence to decompress
 */
void decompress_aligned_aa_sequence(PB_CompressedSequence* input,
									 uint8* output,
									 uint32 from_position,
									 uint32 length)
{
	decode((Varlena*) input, output, from_position, length, fixed_aligned_aa_codes);
}

/*
 *	Section 3 - pgsql interface functions
 */

/**
 *	aligned_aa_sequence_typmod_in()
 *		Condense type modifier keywords into single integer value.
 *
 *	cstring[] input : lower-case keywords separated into array
 */
PG_FUNCTION_INFO_V1 (aligned_aa_sequence_typmod_in);
Datum aligned_aa_sequence_typmod_in(PG_FUNCTION_ARGS)
{
	ArrayType* input = PG_GETARG_ARRAYTYPE_P(0);
	char* read_pointer = ((char*) input) + ARR_DATA_OFFSET(input);

	PB_AlignedAaSequenceTypMod result;
	bool typeModCaseInsensitive = false;
	bool typeModCaseSensitive = false;
	bool typeModIupac = false;
	bool typeModFlc = false;
	bool typeModAscii = false;

	int i;

	PB_TRACE(errmsg("->aligned_aa_sequence_typmod_in()"));

	/*
	 * Parsing
	 */
	for (i = ARR_DIMS(input)[0] - 1; i >= 0; i--) {
		PB_DEBUG2(errmsg("aligned_rna_sequence_typmod_in(): %d= \"%s\".", i, read_pointer));
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
		ereport(ERROR,(errmsg("CASE_INSENSITIVE and CASE_SENSITIVE are mutually exclusive type modifiers.")));
	}

	if ((int) typeModIupac + (int) typeModFlc + (int) typeModAscii > 1) {
		ereport(ERROR,(errmsg("IUPAC and ASCII are mutually exclusive type modifiers.")));
	}

	/*
	 * Build integer value from parsed type modifiers.
	 */
	if (typeModCaseSensitive) {
		result.case_sensitive = PB_ALIGNED_AA_TYPMOD_CASE_SENSITIVE;
	} else {
		result.case_sensitive = PB_ALIGNED_AA_TYPMOD_CASE_INSENSITIVE;
	}

	if (typeModAscii) {
		result.restricting_alphabet = PB_ALIGNED_AA_TYPMOD_ASCII;
	} else {
		result.restricting_alphabet = PB_ALIGNED_AA_TYPMOD_IUPAC;
	}

	PB_TRACE(errmsg("<-aa_sequence_typmod_in() returning %d", aligned_aa_sequence_typmod_to_int(result)));

	PG_RETURN_INT32(aligned_aa_sequence_typmod_to_int(result));
}

/**
 * 	aligned_aa_sequence_typmod_out()
 * 		Restore type modifier keywords from single integer value.
 *
 * 	int typmod : single value representing the type modifiers
 */
PG_FUNCTION_INFO_V1 (aligned_aa_sequence_typmod_out);
Datum aligned_aa_sequence_typmod_out(PG_FUNCTION_ARGS)
{
	int input = PG_GETARG_INT32(0);
	PB_AlignedAaSequenceTypMod typmod = int_to_aligned_aa_sequence_typmod(input);
	int len = 3; /* 2 parentheses + null-terminator */
	char* result;
	char* out;

	/*
	 * Compute length of result.
	 */
	if (typmod.case_sensitive == PB_ALIGNED_AA_TYPMOD_CASE_SENSITIVE) {
		len += 15; /* strlen('CASE_SENSITIVE,') = 15 */
	} else {
		len += 17; /* strlen('CASE_INSENSITIVE,') = 17 */
	}
	len += 5; /* strlen('ASCII') = 5, strlen('IUPAC') = 5  */

	result = palloc0(len);
	out = result;

	*out = '(';
	out++;

	if (typmod.case_sensitive == PB_ALIGNED_AA_TYPMOD_CASE_SENSITIVE) {
		strcpy(out, "CASE_SENSITIVE,");
		out+=15;
	} else {
		strcpy(out, "CASE_INSENSITIVE,");
		out+=17;
	}

	if (typmod.restricting_alphabet == PB_ALIGNED_AA_TYPMOD_ASCII) {
		strcpy(out, "ASCII");
		out+=5;
	} else {
		strcpy(out, "IUPAC");
		out+=5;
	}

	*out = ')';
	out++;
	*out = 0;

	PG_RETURN_INT32(result);
}

/**
 * aligned_aa_sequence_in()
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
PG_FUNCTION_INFO_V1 (aligned_aa_sequence_in);
Datum aligned_aa_sequence_in (PG_FUNCTION_ARGS)
{
	uint8* input = (uint8*) PG_GETARG_CSTRING(0);
	//Oid oid = (Oid) PG_GETARG_OID(1);
	int32 typmod_int = PG_GETARG_INT32(2);
	PB_AlignedAaSequenceTypMod typmod;

	PB_SequenceInfo* info;
	PB_CompressedSequence* result;

	int mode;

	PB_TRACE(errmsg("->aligned_aa_sequence_in()"));

	if ((-1) == typmod_int)
		typmod = non_restricting_aligned_aa_typmod;
	else
		 typmod = int_to_aligned_aa_sequence_typmod(typmod_int);

	/*
	 * Determine sequence info collection mode.
	 */
	mode = 0;
	if (typmod.case_sensitive == PB_ALIGNED_AA_TYPMOD_CASE_SENSITIVE)
		mode |= PB_SEQUENCE_INFO_CASE_SENSITIVE;

	info = get_sequence_info_cstring(input, mode);

	result = compress_aligned_aa_sequence(input, typmod, info);

	PB_SEQUENCE_INFO_PFREE(info);

	PB_TRACE(errmsg("<-aligned_aa_sequence_in()"));

	PG_RETURN_POINTER(result);
}


/**
 * aligned_aa_sequence_in_varlena()
 * 		Compress a given input sequence.
 *
 * 	This function expects a varlena input sequence, that is
 * 	text, varchar or char. It will be called by the respective
 * 	cast functions.
 *
 * 	varlena* input : input sequence
 * 	int typmod : single value representing target type modifier
 */
PG_FUNCTION_INFO_V1 (aligned_aa_sequence_in_varlena);
Datum aligned_aa_sequence_in_varlena (PG_FUNCTION_ARGS)
{
	Varlena* input = (Varlena*) PG_GETARG_VARLENA_P(0);
	PB_AlignedAaSequenceTypMod typmod = int_to_aligned_aa_sequence_typmod(PG_GETARG_INT32(1));

	PB_SequenceInfo* info;
	PB_CompressedSequence* result;

	int mode;

	PB_TRACE(errmsg("->aligned_aa_sequence_in_varlena()"));

	/*
	 * Determine sequence info collection mode.
	 */
	mode = 0;
	if (typmod.case_sensitive == PB_ALIGNED_AA_TYPMOD_CASE_SENSITIVE)
		mode |= PB_SEQUENCE_INFO_CASE_SENSITIVE;

	info = get_sequence_info_text(input, mode);

	result = compress_aligned_aa_sequence((uint8*) VARDATA(input), typmod, info);

	PB_SEQUENCE_INFO_PFREE(info);

	PB_TRACE(errmsg("<-aligned_aa_sequence_in_varlena()"));

	PG_RETURN_POINTER(result);
}

/**
 * aligned_aa_sequence_cast()
 * 		Decompress a given sequence and compress it again
 * 		using a different compression
 *
 * 	varlena* input : input sequence
 * 	int typmod : single value representing target type modifier
 */
PG_FUNCTION_INFO_V1 (aligned_aa_sequence_cast);
Datum aligned_aa_sequence_cast (PG_FUNCTION_ARGS)
{
	PB_CompressedSequence* input = (PB_CompressedSequence*) PG_GETARG_VARLENA_P(0);
	PB_AlignedAaSequenceTypMod typmod = int_to_aligned_aa_sequence_typmod(PG_GETARG_INT32(1));

	uint8* plain;
	PB_SequenceInfo* info;
	PB_CompressedSequence* result;

	int mode ;

	PB_TRACE(errmsg("->aligned_aa_sequence_cast()"));

	/*
	 * Determine sequence info collection mode.
	 */
	mode = 0;
	if (typmod.case_sensitive == PB_ALIGNED_AA_TYPMOD_CASE_SENSITIVE)
		mode |= PB_SEQUENCE_INFO_CASE_SENSITIVE;

	/*
	 * Decompress sequence.
	 */
	plain = palloc0(input->sequence_length + 1);
	decompress_aligned_aa_sequence(input,plain,0,input->sequence_length);
	plain[input->sequence_length] = '\0';

	/*
	 * Compress again.
	 */
	info = get_sequence_info_cstring(plain, mode);
	result = compress_aligned_aa_sequence(plain, typmod, info);

	PB_SEQUENCE_INFO_PFREE(info);

	pfree(plain);

	PB_TRACE(errmsg("<-aligned_aa_sequence_cast()"));

	PG_RETURN_POINTER(result);
}

/**
 * aligned_aa_sequence_out()
 * 		Decompress a sequence.
 *
 * 	CompressedSequence* input : compressed input sequence
 */
PG_FUNCTION_INFO_V1 (aligned_aa_sequence_out);
Datum aligned_aa_sequence_out (PG_FUNCTION_ARGS)
{
	PB_CompressedSequence* input = (PB_CompressedSequence*) PG_GETARG_VARLENA_P(0);
	uint8* result;

	PB_TRACE(errmsg("->aligned_aa_sequence_out()"));

	result = palloc0(input->sequence_length + 1);
	decompress_aligned_aa_sequence(input, result, 0, input->sequence_length);
	result[input->sequence_length] = '\000';

	PB_TRACE(errmsg("<-aligned_aa_sequence_out()"));

	PG_RETURN_CSTRING ((char*)result);
}

/**
 * aligned_aa_sequence_out_varlena()
 * 		Decompress a sequence into varlena.
 *
 * 	CompressedSequence* input : compressed input sequence
 */
PG_FUNCTION_INFO_V1 (aligned_aa_sequence_out_varlena);
Datum aligned_aa_sequence_out_varlena (PG_FUNCTION_ARGS)
{
	PB_CompressedSequence* input = (PB_CompressedSequence*) PG_GETARG_VARLENA_P(0);
	text* result;

	PB_TRACE(errmsg("->aligned_aa_sequence_out_varlena()"));

	result = palloc0(input->sequence_length + VARHDRSZ);
	SET_VARSIZE (result, input->sequence_length + VARHDRSZ);
	decompress_aligned_aa_sequence(input, (uint8*) VARDATA_ANY(result), 0, input->sequence_length);

	PB_TRACE(errmsg("<-aligned_aa_sequence_out_varlena()"));

	PG_RETURN_POINTER(result);
}

/**
 * aligned_aa_sequence_substring()
 * 		Decompress a substring of a sequence.
 *
 * 	This function mimics the originals substr function's
 * 	behaviour. The first position is 1.
 *
 * 	PB_CompressedSequence* input : compressed input sequence
 * 	int start : position to start from
 * 	int len : length of substring
 */
PG_FUNCTION_INFO_V1 (aligned_aa_sequence_substring);
Datum aligned_aa_sequence_substring (PG_FUNCTION_ARGS) {
	Varlena* input = (Varlena*) PG_GETARG_RAW_VARLENA_P(0);
	int start = PG_GETARG_DATUM(1);
	int len = PG_GETARG_DATUM(2);

	PB_CompressedSequence* input_header;

	text* result;

	PB_TRACE(errmsg("->aligned_aa_sequence_substring()"));

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
	decompress_aligned_aa_sequence((PB_CompressedSequence*) input,(uint8*) VARDATA(result),start,len);

	pfree(input_header);

	PB_TRACE(errmsg("<-aligned_aa_sequence_substring()"));

	PG_RETURN_POINTER(result);
}

/**
 * aligned_aa_sequence_char_length()
 * 		Get length of sequence.
 *
 * 	PB_CompressedSequence* input : compressed input sequence
 */
PG_FUNCTION_INFO_V1 (aligned_aa_sequence_char_length);
Datum aligned_aa_sequence_char_length (PG_FUNCTION_ARGS)
{
	PB_CompressedSequence* seq = (PB_CompressedSequence*)
			PG_DETOAST_DATUM_SLICE(PG_GETARG_RAW_VARLENA_P(0),0,4);

	PG_RETURN_INT32(seq->sequence_length);
}

/**
 * aligned_aa_sequence_compression_ratio()
 * 		Get compression ratio.
 *
 * 	The ratio between the size of the sequence as pgsql text
 * 	type and the size of the compressed sequence including all
 * 	required meta-data, such as the substring-index.
 *
 * 	PB_CompressedSequence* input : compressed input sequence
 */
PG_FUNCTION_INFO_V1 (aligned_aa_sequence_compression_ratio);
Datum aligned_aa_sequence_compression_ratio (PG_FUNCTION_ARGS)
{
	Varlena* input = (Varlena*) PG_GETARG_RAW_VARLENA_P(0);
	PB_CompressedSequence* input_header = (PB_CompressedSequence*)
		PG_DETOAST_DATUM_SLICE(input, 0, sizeof(PB_CompressedSequence) - VARHDRSZ);

	double cr = (double) toast_raw_datum_size((Datum)input) / (double) (input_header->sequence_length + VARHDRSZ);

	PB_DEBUG1(errmsg("aligned_aa_sequence_compression_ratio(): cr=%f memsize=%u len=%u", cr, (uint32) toast_raw_datum_size((Datum)input), input_header->sequence_length + VARHDRSZ));

	PG_RETURN_FLOAT8(cr);
}

/**
 * aligned_aa_sequence_reverse()
 * 		Returns the reverse of an aligned aa sequence.
 *
 * 	PB_CompressedSequence* input : compressed input sequence
 */
PG_FUNCTION_INFO_V1 (aligned_aa_sequence_reverse);
Datum aligned_aa_sequence_reverse(PG_FUNCTION_ARGS)
{
	PB_CompressedSequence* input =
			(PB_CompressedSequence*) PG_GETARG_VARLENA_P(0);

	PB_CompressedSequence* result;

	result = reverse(input, fixed_aligned_aa_codes);

	PG_RETURN_POINTER(result);
}

/**
 * get_alphabet_aligned_aa_sequence()
 * 		Calculates alphabet from aligned aa sequence.
 *
 * 	PB_CompressedSequence* input : compressed input sequence
 */
PG_FUNCTION_INFO_V1 (get_alphabet_aligned_aa_sequence);
Datum get_alphabet_aligned_aa_sequence(PG_FUNCTION_ARGS)
{
	PB_CompressedSequence* input =
			(PB_CompressedSequence*) PG_GETARG_VARLENA_P(0);
	PB_Alphabet* result = get_alphabet_compressed_sequence(input, fixed_aligned_aa_codes);

	PG_RETURN_POINTER(result);
}

/**
 * equal_aligned_aa()
 * 		Compares two aligned AA sequence for equality
 *
 * 	Varlena* seq1 : first sequence
 * 	Varlena* seq2 : second sequence
 */
PG_FUNCTION_INFO_V1 (equal_aligned_aa);
Datum equal_aligned_aa(PG_FUNCTION_ARGS)
{
	Varlena* seq1 = (Varlena*) PG_GETARG_RAW_VARLENA_P(0);
	Varlena* seq2 = (Varlena*) PG_GETARG_RAW_VARLENA_P(1);
	bool result;

	PB_TRACE(errmsg("->equal_aligned_aa()"));

	result = sequence_equal(seq1, seq2, fixed_aligned_aa_codes);

	PB_TRACE(errmsg("<-equal_aligned_aa() exists with %d", result));

	PG_RETURN_BOOL(result);
}

/**
 * compare_aligned_aa_lt()
 * 		Compares two aligned AA sequence for equality Less-than.
 *
 * 	Varlena* seq1 : first sequence
 * 	Varlena* seq2 : second sequence
 */
PG_FUNCTION_INFO_V1 (compare_aligned_aa_lt);
Datum compare_aligned_aa_lt(PG_FUNCTION_ARGS)
{
	Varlena* seq1 = (Varlena*) PG_GETARG_RAW_VARLENA_P(0);
	Varlena* seq2 = (Varlena*) PG_GETARG_RAW_VARLENA_P(1);
	bool result = FALSE;

	PB_TRACE(errmsg("->compare_aligned_aa_lt()"));

	if (sequence_compare(seq1, seq2, fixed_aligned_aa_codes) < 0)
		result = TRUE;

	PB_TRACE(errmsg("<-compare_aligned_aa_lt() exists with %d", result));

	PG_RETURN_BOOL(result);
}

/**
 * compare_aligned_aa_le()
 * 		Compares two aligned AA sequence for equality Less or equal.
 *
 * 	Varlena* seq1 : first sequence
 * 	Varlena* seq2 : second sequence
 */
PG_FUNCTION_INFO_V1 (compare_aligned_aa_le);
Datum compare_aligned_aa_le(PG_FUNCTION_ARGS)
{
	Varlena* seq1 = (Varlena*) PG_GETARG_RAW_VARLENA_P(0);
	Varlena* seq2 = (Varlena*) PG_GETARG_RAW_VARLENA_P(1);
	bool result = FALSE;

	PB_TRACE(errmsg("->compare_aligned_aa_le()"));

	if (sequence_compare(seq1, seq2, fixed_aligned_aa_codes) <= 0)
		result = TRUE;

	PB_TRACE(errmsg("<-compare_aligned_aa_le() exists with %d", result));

	PG_RETURN_BOOL(result);
}

/**
 * compare_aligned_aa_gt()
 * 		Compares two aligned AA sequence for equality. greater than
 *
 * 	Varlena* seq1 : first sequence
 * 	Varlena* seq2 : second sequence
 */
PG_FUNCTION_INFO_V1 (compare_aligned_aa_gt);
Datum compare_aligned_aa_gt(PG_FUNCTION_ARGS)
{
	Varlena* seq1 = (Varlena*) PG_GETARG_RAW_VARLENA_P(0);
	Varlena* seq2 = (Varlena*) PG_GETARG_RAW_VARLENA_P(1);
	bool result = FALSE;

	PB_TRACE(errmsg("->compare_aligned_aa_gt()"));

	if (sequence_compare(seq1, seq2, fixed_aligned_aa_codes) > 0)
		result = TRUE;

	PB_TRACE(errmsg("<-compare_aligned_aa_gt() exists with %d", result));

	PG_RETURN_BOOL(result);
}

/**
 * compare_aligned_aa_ge()
 * 		Compares two aligned AA sequence for equality. greater or equal
 *
 * 	Varlena* seq1 : first sequence
 * 	Varlena* seq2 : second sequence
 */
PG_FUNCTION_INFO_V1 (compare_aligned_aa_ge);
Datum compare_aligned_aa_ge(PG_FUNCTION_ARGS)
{
	Varlena* seq1 = (Varlena*) PG_GETARG_RAW_VARLENA_P(0);
	Varlena* seq2 = (Varlena*) PG_GETARG_RAW_VARLENA_P(1);
	bool result = FALSE;

	PB_TRACE(errmsg("->compare_aligned_aa_ge()"));

	if (sequence_compare(seq1, seq2, fixed_aligned_aa_codes) >= 0)
		result = TRUE;

	PB_TRACE(errmsg("<-compare_aligned_aa_ge() exists with %d", result));

	PG_RETURN_BOOL(result);
}

/**
 * compare_aligned_aa()
 * 		Compares two aligned AA sequence for equality.
 *
 * 	Varlena* seq1 : first possibly toasted sequence
 * 	Varlena* seq2 : second possibly toasted sequence
 */
PG_FUNCTION_INFO_V1 (compare_aligned_aa);
Datum compare_aligned_aa(PG_FUNCTION_ARGS)
{
	Varlena* seq1 = (Varlena*) PG_GETARG_RAW_VARLENA_P(0);
	Varlena* seq2 = (Varlena*) PG_GETARG_RAW_VARLENA_P(1);
	int result;

	PB_TRACE(errmsg("->compare_aligned_aa()"));

	result = sequence_compare(seq1, seq2, fixed_aligned_aa_codes);

	PB_TRACE(errmsg("<-compare_aligned_aa() exists with %d", result));

	PG_RETURN_INT32(result);
}

/**
 * hash_aligned_aa()
 * 		Returns a CRC32 for a aligned AA sequence.
 *
 * 	PB_CompressedSequence* seq1 : input sequence
 */
PG_FUNCTION_INFO_V1 (hash_aligned_aa);
Datum hash_aligned_aa(PG_FUNCTION_ARGS)
{
	PB_CompressedSequence* seq1 = (PB_CompressedSequence*) PG_GETARG_VARLENA_P(0);
	uint32 result;

	PB_TRACE(errmsg("->hash_aligned_aa()"));

	result = sequence_crc32(seq1, fixed_aligned_aa_codes);

	PB_TRACE(errmsg("<-hash_aligned_aa() exits with %u", result));

	PG_RETURN_UINT32(result);
}

/**
 * strpos_aligned_aa()
 * 		Finds the first occurrence of a pattern in a sequence.
 */
PG_FUNCTION_INFO_V1 (strpos_aligned_aa);
Datum strpos_aligned_aa(PG_FUNCTION_ARGS)
{
	PB_CompressedSequence* seq = (PB_CompressedSequence*) PG_GETARG_VARLENA_P(0);
	text* search = (text*) PG_GETARG_VARLENA_P(1);
	uint32 result;

	PB_TRACE(errmsg("->strpos_aligned_aa()"));

	result = sequence_strpos(seq, search, fixed_aligned_aa_codes);

	PB_TRACE(errmsg("<-strpos_aligned_aa()"));

	PG_RETURN_UINT32(result);
}

/**
 * octet_length_aligned_aa()
 * 		Returns byte size of datum.
 */
PG_FUNCTION_INFO_V1 (octet_length_aligned_aa);
Datum octet_length_aligned_aa(PG_FUNCTION_ARGS)
{
	PG_RETURN_UINT32((uint32) toast_raw_datum_size((Datum) PG_GETARG_RAW_VARLENA_P(0)));
}
