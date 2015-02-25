/*-------------------------------------------------------------------------
*
* Copyright (c) 2013, Max Planck Institute for Marine Microbiology
*
* This software is released under the PostgreSQL License
*
* Author: Michael Schneider <mschneid@mpi-bremen.de>
*
* IDENTIFICATION
*   src/types/dna_sequence.c
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

#include "types/dna_sequence.h"

/*
 * Section 1 - static fields
 */
static PB_CodeSet dna_flc = {
		.n_symbols = (uint8) 4,
		.max_codeword_length = (uint8) 2,
		.n_swapped_symbols = (uint8) 0,
		.max_swapped_codeword_length = (uint8) 0,
		.has_equal_length = (bool) TRUE,
		.is_fixed = (bool) TRUE,
		.uses_rle = (bool) FALSE,
		.ignore_case = (bool) TRUE,
		.fixed_id = (bool) 0,
		.swap_savings = (uint64) 0,
		.ascii_bitmap_low = (uint64) 0,
		.ascii_bitmap_high = (uint64) 1048714,
		.words = {{0,2,'A'},{64,2,'C'},{128,2,'G'},{192,2,'T'}}
};

static PB_CodeSet dna_flc_cs = {
		.n_symbols = (uint8) 8,
		.max_codeword_length = (uint8) 3,
		.n_swapped_symbols = (uint8) 0,
		.max_swapped_codeword_length = (uint8) 0,
		.has_equal_length = (bool) TRUE,
		.is_fixed = (bool) TRUE,
		.uses_rle = (bool) FALSE,
		.ignore_case = (bool) FALSE,
		.fixed_id = (bool) 1,
		.swap_savings = (uint64) 0,
		.ascii_bitmap_low = (uint64) 0,
		.ascii_bitmap_high = (uint64) 4504192333906058,
		.words = {{0,3,'A'},  {32,3,'C'},  {64,3,'G'}, {96,3,'T'},
				  {128,3,'a'},{160,3,'c'},{192,3,'g'},{224,3,'t'}}
};

static PB_CodeSet dna_iupac = {
		.n_symbols = (uint8) 15,
		.max_codeword_length = (uint8) 8,
		.n_swapped_symbols = (uint8) 0,
		.max_swapped_codeword_length = (uint8) 0,
		.has_equal_length = (bool) FALSE,
		.is_fixed = (bool) TRUE,
		.uses_rle = (bool) FALSE,
		.ignore_case = (bool) TRUE,
		.fixed_id = (bool) 2,
		.swap_savings = (uint64) 0,
		.ascii_bitmap_low = (uint64) 0,
		.ascii_bitmap_high = (uint64) 47999390,
		.words = {{  0,2,'A'},{ 64,2,'C'},{128,2,'G'},{192,3,'T'},
				  {224,4,'N'},{240,7,'M'},{242,7,'R'},{244,7,'Y'},
				  {246,7,'W'},{248,7,'B'},{250,7,'V'},{252,8,'S'},
				  {253,8,'K'},{254,8,'D'},{255,8,'H'}}
};

static PB_CodeSet dna_iupac_cs = {
		.n_symbols = (uint8) 30,
		.max_codeword_length = (uint8) 8,
		.n_swapped_symbols = (uint8) 0,
		.max_swapped_codeword_length = (uint8) 0,
		.has_equal_length = (bool) FALSE,
		.is_fixed = (bool) TRUE,
		.uses_rle = (bool) FALSE,
		.ignore_case = (bool) FALSE,
		.fixed_id = (bool) 3,
		.swap_savings = (uint64) 0,
		.ascii_bitmap_low = (uint64) 0,
		.ascii_bitmap_high = (uint64) 206155810325948830,
		.words = {{  0,3,'A'},{ 64,3,'C'},{128,3,'G'},{192,4,'T'},
				  {224,6,'N'},{232,7,'Y'},{236,7,'R'},{240,8,'M'},
				  {242,8,'W'},{244,8,'B'},{246,8,'V'},{248,8,'S'},
				  {250,8,'K'},{252,8,'D'},{254,8,'H'},
				  { 32,3,'a'},{ 96,3,'c'},{160,3,'g'},{208,4,'t'},
				  {228,6,'n'},{234,7,'y'},{238,7,'r'},{241,8,'m'},
				  {243,8,'w'},{245,8,'b'},{247,8,'v'},{249,8,'s'},
				  {251,8,'k'},{253,8,'d'},{255,8,'h'}}
};

static PB_CodeSet dna_flc_complement = {
		.n_symbols = (uint8) 4,
		.max_codeword_length = (uint8) 2,
		.n_swapped_symbols = (uint8) 0,
		.max_swapped_codeword_length = (uint8) 0,
		.has_equal_length = (bool) TRUE,
		.is_fixed = (bool) TRUE,
		.uses_rle = (bool) FALSE,
		.ignore_case = (bool) TRUE,
		.fixed_id = (bool) 4,
		.swap_savings = (uint64) 0,
		.ascii_bitmap_low = (uint64) 0,
		.ascii_bitmap_high = (uint64) 1048714,
		.words = {{0,2,'T'},{64,2,'G'},{128,2,'C'},{192,2,'A'}}
};

static PB_CodeSet dna_flc_cs_complement = {
		.n_symbols = (uint8) 8,
		.max_codeword_length = (uint8) 3,
		.n_swapped_symbols = (uint8) 0,
		.max_swapped_codeword_length = (uint8) 0,
		.has_equal_length = (bool) TRUE,
		.is_fixed = (bool) TRUE,
		.uses_rle = (bool) FALSE,
		.ignore_case = (bool) FALSE,
		.fixed_id = (bool) 5,
		.swap_savings = (uint64) 0,
		.ascii_bitmap_low = (uint64) 0,
		.ascii_bitmap_high = (uint64) 4504192333906058,
		.words = {{0,3,'T'},  {32,3,'G'},  {64,3,'C'}, {96,3,'A'},
				  {128,3,'t'},{160,3,'g'},{192,3,'c'},{224,3,'a'}}
};

static PB_CodeSet dna_iupac_complement = {
		.n_symbols = (uint8) 15,
		.max_codeword_length = (uint8) 8,
		.n_swapped_symbols = (uint8) 0,
		.max_swapped_codeword_length = (uint8) 0,
		.has_equal_length = (bool) FALSE,
		.is_fixed = (bool) TRUE,
		.uses_rle = (bool) FALSE,
		.ignore_case = (bool) TRUE,
		.fixed_id = (bool) 6,
		.swap_savings = (uint64) 0,
		.ascii_bitmap_low = (uint64) 0,
		.ascii_bitmap_high = (uint64) 47999390,
		.words = {{  0,2,'T'},{ 64,2,'G'},{128,2,'C'},{192,3,'A'},
				  {224,4,'N'},{240,7,'K'},{242,7,'Y'},{244,7,'R'},
				  {246,7,'W'},{248,7,'V'},{250,7,'B'},{252,8,'S'},
				  {253,8,'M'},{254,8,'H'},{255,8,'D'}}
};

static PB_CodeSet dna_iupac_cs_complement = {
		.n_symbols = (uint8) 30,
		.max_codeword_length = (uint8) 8,
		.n_swapped_symbols = (uint8) 0,
		.max_swapped_codeword_length = (uint8) 0,
		.has_equal_length = (bool) FALSE,
		.is_fixed = (bool) TRUE,
		.uses_rle = (bool) FALSE,
		.ignore_case = (bool) FALSE,
		.fixed_id = (bool) 7,
		.swap_savings = (uint64) 0,
		.ascii_bitmap_low = (uint64) 0,
		.ascii_bitmap_high = (uint64) 206155810325948830,
		.words = {{  0,3,'T'},{ 64,3,'G'},{128,3,'C'},{192,4,'A'},
				  {224,6,'N'},{232,7,'R'},{236,7,'Y'},{240,8,'K'},
				  {242,8,'W'},{244,8,'V'},{246,8,'B'},{248,8,'S'},
				  {250,8,'M'},{252,8,'H'},{254,8,'D'},
				  { 32,3,'t'},{ 96,3,'g'},{160,3,'c'},{208,4,'a'},
				  {228,6,'n'},{234,7,'r'},{238,7,'y'},{241,8,'k'},
				  {243,8,'w'},{245,8,'v'},{247,8,'b'},{249,8,'s'},
				  {251,8,'m'},{253,8,'h'},{255,8,'d'}}
};

static unsigned int n_fixed_dna_codes = 8;
static PB_CodeSet* fixed_dna_codes[] = {
		&dna_flc,
		&dna_flc_cs,
		&dna_iupac,
		&dna_iupac_cs,
		&dna_flc_complement,
		&dna_flc_cs_complement,
		&dna_iupac_complement,
		&dna_iupac_cs_complement
};

PB_DnaSequenceTypMod non_restricting_dna_typmod = {
	.case_sensitive = PB_DNA_TYPMOD_CASE_SENSITIVE,
	.restricting_alphabet = PB_DNA_TYPMOD_ASCII,
	.compression_strategy = PB_DNA_TYPMOD_SHORT
};

/*
 * Section 2 - other public functions
 */

/**
 * dna_sequence_typmod_to_int()
 * 		Convert from PB_DnaSequenceTypMod to int
 *
 * 	PB_DnaSequenceTypMod typmod : type modifier
 */
int dna_sequence_typmod_to_int(PB_DnaSequenceTypMod typmod)
{
	return *((int*) &typmod);
}

/**
 * int_to_dna_sequence_typmod()
 * 		Convert from int to PB_DnaSequenceTypMod
 *
 * 	int typmod : type modifier
 */
PB_DnaSequenceTypMod int_to_dna_sequence_typmod(int typmod)
{
	if (-1 == typmod) {
		typmod = 0;
	}
	return *((PB_DnaSequenceTypMod*) &typmod);
}

/**
 * get_fixed_dna_code()
 * 		Returns a fixed code for the specified id.
 *
 * 	Id	|	Code Description
 * -----------------------------------------------------------
 * 	0	|	DNA four-letter code
 * 	1	|	DNA four-letter code, case sensitive
 * 	2	|	DNA IUPAC code
 * 	3	|	DNA IUPAC code, case sensitive
 * 	4	|	DNA four-letter code complement
 * 	5	|	DNA four-letter code complement, case sensitive
 * 	6	|	DNA IUPAC code complement
 * 	7	|	DNA IUPAC code complement, case sensitive
 *
 */
PB_CodeSet* get_fixed_dna_code(unsigned int fixed_code_id)
{
	PB_CodeSet* result = NULL;

	if (fixed_code_id < n_fixed_dna_codes)
		result = fixed_dna_codes[fixed_code_id];

	return result;
}

/**
 * get_fixed_dna_codes()
 * 		Returns pointer to fixed DNA codes.
 */
PB_CodeSet** get_fixed_dna_codes(void)
{
	return fixed_dna_codes;
}

/**
 * compress_dna_sequence()
 * 		Compress a DNA sequence.
 *
 *	First it will be checked whether the given sequence matches
 *	the restricting alphabet specified by the type modifiers.
 *	Then one of three types of compression will be performed.
 *
 *	A) Compression with an applicable built-in code-set will be
 *		performed on sequences:
 *		* with type modifier FLC
 *		* shorter than 128 chars and no ASCII type modifier
 *		* with type modifier SHORT and no ASCII type modifier
 *	B) Huffman-Coding, RLE and Rare-symbol-swapping will be performed
 *		on sequences:
 *		* with type modifier REFERENCE
 *	C) Huffman-Coding and Rare-symbol-swapping will be user for
 *		sequences:
 *		* with type modifier DEFAULT
 *		* all other
 *
 *	Run-length encoding and rare-symbol-swapping will, of course,
 *	only be employed if it actually reduces total size.
 *
 * 	uint8* input : unterminated or null-terminated input sequence
 * 	PB_DnaSequenceTypMod typmod : target type modifier
 * 	PB_SequenceInfo info : given by input function
 */
PB_CompressedSequence* compress_dna_sequence(uint8* input,
											 PB_DnaSequenceTypMod typmod,
											 PB_SequenceInfo* info)
{
	PB_CompressedSequence* result;
	PB_CodeSet* code_set = NULL;

	PB_TRACE(errmsg("->compress_dna_sequence(), typmod=%d",dna_sequence_typmod_to_int(typmod)));

	PB_DEBUG1(errmsg("compress_dna_sequence(): low bitmap:%ld high bitmap:%ld",  info->ascii_bitmap_low, info->ascii_bitmap_high));

	/*
	 * Check alphabet constraints.
	 */
	if ((typmod.restricting_alphabet == PB_DNA_TYPMOD_FLC && (!PB_CHECK_CODESET((&dna_flc_cs),info)) ) ||
		(typmod.restricting_alphabet == PB_DNA_TYPMOD_IUPAC  && (!PB_CHECK_CODESET((&dna_iupac_cs),info)) ))
	{
		ereport(ERROR,(errmsg("input sequence violates alphabet restrictions")));
	}

	/*
	 * Choose code set.
	 */
	if (typmod.restricting_alphabet == PB_DNA_TYPMOD_FLC ||
		((info->sequence_length < 128  || typmod.compression_strategy == PB_DNA_TYPMOD_SHORT)
		   && typmod.restricting_alphabet != PB_DNA_TYPMOD_ASCII))
	{
		uint8 fixed_code_id;

		/*
		 * Choose applicable pre-built code set.
		 */
		if (typmod.restricting_alphabet == PB_DNA_TYPMOD_IUPAC) {
			if (typmod.case_sensitive == PB_DNA_TYPMOD_CASE_SENSITIVE) {
				fixed_code_id = 3;
			} else {
				fixed_code_id = 2;
			}
		} else {
			if (typmod.case_sensitive == PB_DNA_TYPMOD_CASE_SENSITIVE) {
				fixed_code_id = 1;
			} else {
				fixed_code_id = 0;
			}
		}

		/*
		 * Check, if cheaper code is qualified to encode sequence.
		 */
		if (fixed_code_id > 0 && PB_CHECK_CODESET((&dna_flc), info)) {
			fixed_code_id = 0;
		} else if (fixed_code_id > 1 && PB_CHECK_CODESET((&dna_flc_cs), info)) {
			fixed_code_id = 1;
		} else if (fixed_code_id > 2 && PB_CHECK_CODESET((&dna_iupac), info)) {
			fixed_code_id = 2;
		}

		code_set = fixed_dna_codes[fixed_code_id];

		PB_DEBUG1(errmsg("compress_dna_sequence(): using fixed_code %d", fixed_code_id));

	}
	else
	{
		/*
		 * Build optimal code.
		 */
		code_set = get_optimal_code(info);
	}

	/*
	 * Compress.
	 */
	result = encode(input, get_compressed_size(info, code_set), code_set, info);

	if (!code_set->is_fixed)
		pfree(code_set);

	PB_TRACE(errmsg("<-compress_dna_sequence()"));

	return result;
}

/**
 * decompress_dna_sequence()
 * 		Decompress a DNA sequence
 *
 * 	PB_CompressedSequence* input : compressed input sequence (takes unTOASTed data)
 * 	uint8* output : already allocated appropriately sized target memory area
 * 	uint32 from_position : position to start decompressing, first = 0
 * 	uint32 length : length of sequence to decompress
 */
void decompress_dna_sequence(PB_CompressedSequence* input,
							 uint8* output,
							 uint32 from_position,
							 uint32 length)
{
	decode((Varlena*) input, output, from_position, length, fixed_dna_codes);
}

/*
 *	Section 3 - pgsql interface functions
 */

/**
 *	dna_sequence_typmod_in()
 *		Condense type modifier keywords into single integer value.
 *
 *	cstring[] input : lower-case keywords separated into array
 */
PG_FUNCTION_INFO_V1 (dna_sequence_typmod_in);
Datum dna_sequence_typmod_in(PG_FUNCTION_ARGS)
{
	ArrayType* input = PG_GETARG_ARRAYTYPE_P(0);
	char* read_pointer = ((char*) input) + ARR_DATA_OFFSET(input);

	PB_DnaSequenceTypMod result;

	bool typeModCaseInsensitive = false;
	bool typeModCaseSensitive = false;
	bool typeModIupac = false;
	bool typeModFlc = false;
	bool typeModAscii = false;
	bool typeModDefault = false;
	bool typeModShortRead = false;
	bool typeModRef = false;

	int i;

	PB_TRACE(errmsg("->dna_sequence_typmod_in()"));

	/*
	 * Parsing
	 */
	for (i = ARR_DIMS(input)[0] - 1; i >= 0; i--)
	{
		PB_DEBUG2(errmsg("dna_sequence_typmod_in(): %d= \"%s\".", i, read_pointer));
		if (!strcmp(read_pointer, "case_insensitive")) {
			typeModCaseInsensitive = true;
		} else if (!strcmp(read_pointer, "case_sensitive")) {
			typeModCaseSensitive = true;
		} else if (!strcmp(read_pointer, "iupac")) {
			typeModIupac = true;
		} else if (!strcmp(read_pointer, "flc")) {
			typeModFlc = true;
		} else if (!strcmp(read_pointer, "ascii")) {
			typeModAscii = true;
		} else if (!strcmp(read_pointer, "default")) {
			typeModDefault = true;
		} else if (!strcmp(read_pointer, "short")) {
			typeModShortRead = true;
		} else if (!strcmp(read_pointer, "reference")) {
			typeModRef = true;
		} else {
			ereport(ERROR,(errmsg("type modifier invalid"),
					errdetail("Can not recognize type modifier \"%s\".", read_pointer)));
		}
		read_pointer += strlen(read_pointer) + 1;
	}

	/*
	 * Syntax checking.
	 */
	if (typeModCaseInsensitive && typeModCaseSensitive)
	{
		ereport(ERROR,(errmsg("CASE_INSENSITIVE and CASE_SENSITIVE are mutually exclusive type modifiers")));
	}

	if ((int) typeModIupac + (int) typeModFlc + (int) typeModAscii > 1)
	{
		ereport(ERROR,(errmsg("IUPAC, FLC and ASCII are mutually exclusive type modifiers")));
	}

	if ((int) typeModDefault + (int) typeModShortRead + (int) typeModRef > 1)
	{
		ereport(ERROR,(errmsg("DEFAULT, SHORT and REFERENCE are mutually exclusive type modifiers")));
	}

	/*
	 * Build integer value from parsed type modifiers.
	 */
	if (typeModCaseSensitive) {
		result.case_sensitive = PB_DNA_TYPMOD_CASE_SENSITIVE;
	} else {
		result.case_sensitive = PB_DNA_TYPMOD_CASE_INSENSITIVE;
	}

	if (typeModShortRead) {
		result.compression_strategy = PB_DNA_TYPMOD_SHORT;
	} else if (typeModRef) {
		result.compression_strategy = PB_DNA_TYPMOD_REFERENCE;
	} else {
		result.compression_strategy = PB_DNA_TYPMOD_DEFAULT;
	}

	if (typeModFlc) {
		result.restricting_alphabet = PB_DNA_TYPMOD_FLC;
	} else if (typeModAscii) {
		result.restricting_alphabet = PB_DNA_TYPMOD_ASCII;
	} else {
		result.restricting_alphabet = PB_DNA_TYPMOD_IUPAC;
	}

	PB_TRACE(errmsg("<-dna_sequence_typmod_in() returning %d", dna_sequence_typmod_to_int(result)));

	PG_RETURN_INT32(dna_sequence_typmod_to_int(result));
}

/**
 * 	dna_sequence_typmod_out()
 * 		Restore type modifier keywords from single integer value.
 *
 * 	int typmod : single value representing the type modifiers
 */
PG_FUNCTION_INFO_V1 (dna_sequence_typmod_out);
Datum dna_sequence_typmod_out(PG_FUNCTION_ARGS)
{
	int input = PG_GETARG_INT32(0);
	PB_DnaSequenceTypMod typmod = int_to_dna_sequence_typmod(input);
	int len = 3; /* 2 parentheses + null-terminator */
	char* result;
	char* out;

	PB_TRACE(errmsg("->dna_sequence_typmod_out(%d)", input));

	/*
	 * Compute length of result.
	 */
	if (typmod.case_sensitive == PB_DNA_TYPMOD_CASE_SENSITIVE) {
		len += 15; /* strlen('CASE_SENSITIVE,') = 15 */
	} else {
		len += 17; /* strlen('CASE_INSENSITIVE,') = 17 */
	}

	if (typmod.compression_strategy == PB_DNA_TYPMOD_SHORT) {
		len += 11; /* strlen('SHORT_READ,') = 11 */
	} else if (typmod.compression_strategy == PB_DNA_TYPMOD_REFERENCE) {
		len += 9; /* strlen('REFERENCE,') = 9 */
	} else {
		len += 8; /* strlen('DEFAULT,') = 8 */
	}

	if (typmod.restricting_alphabet == PB_DNA_TYPMOD_FLC) {
		len += 3; /* strlen('FLC') = 3 */
	} else {
		len += 5; /* strlen('IUPAC') = 5, strlen('ASCII') = 5  */
	}

	result = palloc0(len);
	out = result;

	*out = '(';
	out++;

	if (typmod.case_sensitive == PB_DNA_TYPMOD_CASE_SENSITIVE) {
		strcpy(out, "CASE_SENSITIVE,");
		out+=15;
	} else {
		strcpy(out, "CASE_INSENSITIVE,");
		out+=17;
	}

	if (typmod.compression_strategy == PB_DNA_TYPMOD_SHORT) {
		strcpy(out, "SHORT_READ,");
		out+=11;
	} else if (typmod.compression_strategy == PB_DNA_TYPMOD_REFERENCE) {
		strcpy(out, "REFERENCE,");
		out+=19;
	} else {
		strcpy(out, "DEFAULT,");
		out+=8;
	}

	if (typmod.restricting_alphabet == PB_DNA_TYPMOD_FLC) {
		strcpy(out, "FLC");
		out+=3;
	} else if (typmod.restricting_alphabet == PB_DNA_TYPMOD_ASCII) {
		strcpy(out, "ASCII");
		out+=5;
	} else {
		strcpy(out, "IUPAC");
		out+=5;
	}

	*out = ')';
	out++;
	*out = 0;

	PB_TRACE(errmsg("<-dna_sequence_typmod_out() returning '%s'", out));

	PG_RETURN_CSTRING(result);
}

/**
 * dna_sequence_in()
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
PG_FUNCTION_INFO_V1 (dna_sequence_in);
Datum dna_sequence_in (PG_FUNCTION_ARGS)
{
	uint8* input = (uint8*) PG_GETARG_CSTRING(0);
	//Oid oid = (Oid) PG_GETARG_OID(1);
	int32 typmod_int = PG_GETARG_INT32(2);
	PB_DnaSequenceTypMod typmod;

	PB_SequenceInfo* info;
	PB_CompressedSequence* result;

	int mode;

	PB_TRACE(errmsg("->dna_sequence_in()"));

	if ((-1) == typmod_int)
		typmod = non_restricting_dna_typmod;
	else
		 typmod = int_to_dna_sequence_typmod(typmod_int);
	/*
	 * Determine sequence info collection mode.
	 */
	mode = 0;
	if (typmod.compression_strategy == PB_DNA_TYPMOD_REFERENCE)
		mode = PB_SEQUENCE_INFO_WITH_RLE;
	if (typmod.case_sensitive == PB_DNA_TYPMOD_CASE_SENSITIVE)
		mode |= PB_SEQUENCE_INFO_CASE_SENSITIVE;

	info = get_sequence_info_cstring(input, mode);

	result = compress_dna_sequence(input, typmod, info);

	PB_SEQUENCE_INFO_PFREE(info);

	PB_TRACE(errmsg("<-dna_sequence_in()"));

	PG_RETURN_POINTER(result);
}

/**
 * dna_sequence_in_varlena()
 * 		Compress a given input sequence.
 *
 * 	This function expects a varlena input sequence, that is
 * 	text, varchar or char. It will be called by the respective
 * 	cast functions.
 *
 * 	varlena* input : input sequence
 * 	int typmod : single value representing target type modifier
 */
PG_FUNCTION_INFO_V1 (dna_sequence_in_varlena);
Datum dna_sequence_in_varlena (PG_FUNCTION_ARGS)
{
	Varlena* input = (Varlena*) PG_GETARG_VARLENA_P(0);
	PB_DnaSequenceTypMod typmod = int_to_dna_sequence_typmod(PG_GETARG_INT32(1));

	PB_SequenceInfo* info;
	PB_CompressedSequence* result;

	int mode;

	PB_TRACE(errmsg("->dna_sequence_in_varlena()"));

	/*
	 * Determine sequence info collection mode.
	 */
	mode = 0;
	if (typmod.compression_strategy == PB_DNA_TYPMOD_REFERENCE)
		mode = PB_SEQUENCE_INFO_WITH_RLE;
	if (typmod.case_sensitive == PB_DNA_TYPMOD_CASE_SENSITIVE)
		mode |= PB_SEQUENCE_INFO_CASE_SENSITIVE;

	info = get_sequence_info_text(input, mode);

	result = compress_dna_sequence((uint8*) VARDATA(input), typmod, info);

	PB_SEQUENCE_INFO_PFREE(info);

	PB_TRACE(errmsg("<-dna_sequence_in_varlena()"));

	PG_RETURN_POINTER(result);
}

/**
 * dna_sequence_cast()
 * 		Decompress a given sequence and compress it again
 * 		using a different compression.
 *
 * 	PB_CompressedSequence* input : input sequence
 * 	int typmod : single value representing target type modifier
 */
PG_FUNCTION_INFO_V1 (dna_sequence_cast);
Datum dna_sequence_cast (PG_FUNCTION_ARGS)
{
	PB_CompressedSequence* input =
		(PB_CompressedSequence*) PG_GETARG_VARLENA_P(0);
	PB_DnaSequenceTypMod typmod = int_to_dna_sequence_typmod(PG_GETARG_INT32(1));

	uint8* plain;
	PB_SequenceInfo* info;
	PB_CompressedSequence* result;

	int mode ;

	PB_TRACE(errmsg("->dna_sequence_cast()"));

	/*
	 * Determine sequence info collection mode.
	 */
	mode = 0;
	if (typmod.compression_strategy == PB_DNA_TYPMOD_REFERENCE)
		mode = PB_SEQUENCE_INFO_WITH_RLE;

	if (typmod.case_sensitive == PB_DNA_TYPMOD_CASE_SENSITIVE)
		mode |= PB_SEQUENCE_INFO_CASE_SENSITIVE;

	/*
	 * Decompress sequence.
	 */
	plain = palloc0(input->sequence_length + 1);
	decompress_dna_sequence(input,plain,0,input->sequence_length);
	plain[input->sequence_length] = '\0';

	/*
	 * Compress again.
	 */
	info = get_sequence_info_cstring(plain, mode);
	result = compress_dna_sequence(plain, typmod, info);

	PB_SEQUENCE_INFO_PFREE(info);

	pfree(plain);

	PB_TRACE(errmsg("<-dna_sequence_cast()"));

	PG_RETURN_POINTER(result);
}

/**
 * dna_sequence_out()
 * 		Decompress a sequence.
 *
 * 	PB_CompressedSequence* input : compressed input sequence
 */
PG_FUNCTION_INFO_V1 (dna_sequence_out);
Datum dna_sequence_out (PG_FUNCTION_ARGS)
{
	PB_CompressedSequence* input = (PB_CompressedSequence*) PG_GETARG_VARLENA_P(0);
	uint8* result;

	PB_TRACE(errmsg("->dna_sequence_out()"));

	result = palloc0(input->sequence_length + 1);
	decompress_dna_sequence(input,result,0,input->sequence_length);
	result[input->sequence_length] = '\000';

	PB_TRACE(errmsg("<-dna_sequence_out()"));

	PG_RETURN_CSTRING ((char*)result);
}

/**
 * dna_sequence_out_varlena()
 * 		Decompress a sequence into varlena.
 *
 * 	PB_CompressedSequence* input : compressed input sequence
 */
PG_FUNCTION_INFO_V1 (dna_sequence_out_varlena);
Datum dna_sequence_out_varlena (PG_FUNCTION_ARGS)
{
	PB_CompressedSequence* input = (PB_CompressedSequence*) PG_GETARG_VARLENA_P(0);
	text* result;

	PB_TRACE(errmsg("->dna_sequence_out_varlena()"));

	result = palloc0(input->sequence_length + VARHDRSZ);
	SET_VARSIZE (result, input->sequence_length + VARHDRSZ);
	decompress_dna_sequence(input,(uint8*)VARDATA_ANY(result),0,input->sequence_length);

	PB_TRACE(errmsg("<-dna_sequence_out_varlena()"));

	PG_RETURN_POINTER(result);
}

/**
 * dna_sequence_substring()
 * 		Decompress a substring of a sequence.
 *
 * 	This function mimics the behaviour of the originals substr function.
 * 	The first position is 1.
 *
 * 	Varlena* input : compressed input sequence
 * 	int start : position to start from
 * 	int len : length of substring
 */
PG_FUNCTION_INFO_V1 (dna_sequence_substring);
Datum dna_sequence_substring (PG_FUNCTION_ARGS)
{
	Varlena* input = (Varlena*) PG_GETARG_RAW_VARLENA_P(0);
	int start = PG_GETARG_DATUM(1);
	int len = PG_GETARG_DATUM(2);

	PB_CompressedSequence* input_header;

	text* result;

	PB_TRACE(errmsg("->dna_sequence_substring()"));

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
	decompress_dna_sequence((PB_CompressedSequence*) input,(uint8*) VARDATA(result),start,len);

	pfree(input_header);

	PB_TRACE(errmsg("<-dna_sequence_substring()"));

	PG_RETURN_POINTER(result);
}

/**
 * dna_sequence_char_length()
 * 		Get length of sequence.
 *
 * 	PB_CompressedSequence* input : compressed input sequence
 */
PG_FUNCTION_INFO_V1 (dna_sequence_char_length);
Datum dna_sequence_char_length (PG_FUNCTION_ARGS)
{
	PB_CompressedSequence* dna_seq = (PB_CompressedSequence*)
			PG_DETOAST_DATUM_SLICE(PG_GETARG_RAW_VARLENA_P(0),0,4);

	PG_RETURN_INT32(dna_seq->sequence_length);
}

/**
 * dna_sequence_compression_ratio()
 * 		Get compression ratio.
 *
 * 	The ratio between the size of the sequence as pgsql text
 * 	type and the size of the compressed sequence including all
 * 	required meta-data, such as the substring-index.
 *
 * 	PB_CompressedSequence* input : compressed input sequence
 */
PG_FUNCTION_INFO_V1 (dna_sequence_compression_ratio);
Datum dna_sequence_compression_ratio (PG_FUNCTION_ARGS)
{
	Varlena* input = (Varlena*) PG_GETARG_RAW_VARLENA_P(0);
	PB_CompressedSequence* input_header = (PB_CompressedSequence*)
		PG_DETOAST_DATUM_SLICE(input, 0, sizeof(PB_CompressedSequence) - VARHDRSZ);

	double cr = (double) toast_raw_datum_size((Datum)input) / (double) (input_header->sequence_length + VARHDRSZ);

	PB_DEBUG1(errmsg("dna_sequence_compression_ratio(): cr=%f memsize=%u len=%u", cr, (uint32) toast_raw_datum_size((Datum)input), input_header->sequence_length + VARHDRSZ));

	PG_RETURN_FLOAT8(cr);
}

static void complement_dna(PB_CompressedSequence* sequence)
{
	PB_Codeword* codewords = PB_COMPRESSED_SEQUENCE_SYMBOL_POINTER(sequence);

	if (sequence->is_fixed)
	{
		sequence->n_swapped_symbols = sequence->n_swapped_symbols ^ 0x4;
	}
	else
	{
		int i;

		for (i = 0; i < sequence->n_symbols; i++)
		{
			uint8 symbol = codewords[i].symbol;
			switch (symbol)
			{
			case 'A':
				symbol = 'T';
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
				symbol = 't';
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
}

/**
 * dna_sequence_complement()
 * 		Returns the complement of a DNA sequence.
 *
 * 	PB_CompressedSequence* input : compressed input sequence
 */
PG_FUNCTION_INFO_V1 (dna_sequence_complement);
Datum dna_sequence_complement(PG_FUNCTION_ARGS)
{
	PB_CompressedSequence* input = (PB_CompressedSequence*)
			PG_GETARG_VARLENA_P(0);
	PB_CompressedSequence* result;

	uint32 mem_size = VARSIZE(input);

	PB_TRACE(errmsg("->dna_sequence_complement()"));

	/*
	 * Make a copy of the original.
	 */
	result = palloc0(mem_size);
	memcpy(result, input, mem_size);

	complement_dna(result);

	PB_TRACE(errmsg("<-dna_sequence_complement()"));

	PG_RETURN_POINTER(result);
}

/**
 * dna_sequence_reverse()
 * 		Returns the reverse of a DNA sequence.
 *
 * 	PB_CompressedSequence* input : compressed input sequence
 */
PG_FUNCTION_INFO_V1 (dna_sequence_reverse);
Datum dna_sequence_reverse(PG_FUNCTION_ARGS)
{
	PB_CompressedSequence* input =
			(PB_CompressedSequence*) PG_GETARG_VARLENA_P(0);

	PB_CompressedSequence* result;

	result = reverse(input, fixed_dna_codes);

	PG_RETURN_POINTER(result);
}

/**
 * dna_sequence_reverse_complement()
 * 		Returns the reverse of a DNA sequence.
 *
 * 	PB_CompressedSequence* input : compressed input sequence
 */
PG_FUNCTION_INFO_V1 (dna_sequence_reverse_complement);
Datum dna_sequence_reverse_complement(PG_FUNCTION_ARGS)
{
	PB_CompressedSequence* input =
			(PB_CompressedSequence*) PG_GETARG_VARLENA_P(0);

	PB_CompressedSequence* result;

	result = reverse(input, fixed_dna_codes);

	complement_dna(result);

	PG_RETURN_POINTER(result);
}

/**
 * get_alphabet_dna_sequence()
 * 		Calculates alphabet from DNA sequence.
 *
 * 	PB_CompressedSequence* input : compressed input sequence
 */
PG_FUNCTION_INFO_V1 (get_alphabet_dna_sequence);
Datum get_alphabet_dna_sequence(PG_FUNCTION_ARGS)
{
	PB_CompressedSequence* input =
			(PB_CompressedSequence*) PG_GETARG_VARLENA_P(0);
	PB_Alphabet* result = get_alphabet_compressed_sequence(input, fixed_dna_codes);

	PG_RETURN_POINTER(result);
}

/**
 * equal_dna()
 * 		Compares two DNA sequence for equality
 *
 * 	Varlena* seq1 : first sequence
 * 	Varlena* seq2 : second sequence
 */
PG_FUNCTION_INFO_V1 (equal_dna);
Datum equal_dna(PG_FUNCTION_ARGS)
{
	Varlena* seq1 = (Varlena*) PG_GETARG_RAW_VARLENA_P(0);
	Varlena* seq2 = (Varlena*) PG_GETARG_RAW_VARLENA_P(1);
	bool result;

	PB_TRACE(errmsg("->equal_dna()"));

	result = sequence_equal(seq1, seq2, fixed_dna_codes);

	PB_TRACE(errmsg("<-equal_dna() exists with %d", result));

	PG_RETURN_BOOL(result);
}

/**
 * compare_dna_lt()
 * 		Compares two DNA sequence for equality Less-than.
 *
 * 	Varlena* seq1 : first sequence
 * 	Varlena* seq2 : second sequence
 */
PG_FUNCTION_INFO_V1 (compare_dna_lt);
Datum compare_dna_lt(PG_FUNCTION_ARGS)
{
	Varlena* seq1 = (Varlena*) PG_GETARG_RAW_VARLENA_P(0);
	Varlena* seq2 = (Varlena*) PG_GETARG_RAW_VARLENA_P(1);
	bool result = FALSE;

	PB_TRACE(errmsg("->compare_dna_lt()"));

	if (sequence_compare(seq1, seq2, fixed_dna_codes) < 0)
		result = TRUE;

	PB_TRACE(errmsg("<-compare_dna_lt() exists with %d", result));

	PG_RETURN_BOOL(result);
}

/**
 * compare_dna_le()
 * 		Compares two DNA sequence for equality Less or equal.
 *
 * 	Varlena* seq1 : first sequence
 * 	Varlena* seq2 : second sequence
 */
PG_FUNCTION_INFO_V1 (compare_dna_le);
Datum compare_dna_le(PG_FUNCTION_ARGS)
{
	Varlena* seq1 = (Varlena*) PG_GETARG_RAW_VARLENA_P(0);
	Varlena* seq2 = (Varlena*) PG_GETARG_RAW_VARLENA_P(1);
	bool result = FALSE;

	PB_TRACE(errmsg("->compare_dna_le()"));

	if (sequence_compare(seq1, seq2, fixed_dna_codes) <= 0)
		result = TRUE;

	PB_TRACE(errmsg("<-compare_dna_le() exists with %d", result));

	PG_RETURN_BOOL(result);
}

/**
 * compare_dna_gt()
 * 		Compares two DNA sequence for equality. greater than
 *
 * 	Varlena* seq1 : first sequence
 * 	Varlena* seq2 : second sequence
 */
PG_FUNCTION_INFO_V1 (compare_dna_gt);
Datum compare_dna_gt(PG_FUNCTION_ARGS)
{
	Varlena* seq1 = (Varlena*) PG_GETARG_RAW_VARLENA_P(0);
	Varlena* seq2 = (Varlena*) PG_GETARG_RAW_VARLENA_P(1);
	bool result = FALSE;

	PB_TRACE(errmsg("->compare_dna_gt()"));

	if (sequence_compare(seq1, seq2, fixed_dna_codes) > 0)
		result = TRUE;

	PB_TRACE(errmsg("<-compare_dna_gt() exists with %d", result));

	PG_RETURN_BOOL(result);
}

/**
 * compare_dna_ge()
 * 		Compares two DNA sequence for equality. greater or equal
 *
 * 	Varlena* seq1 : first sequence
 * 	Varlena* seq2 : second sequence
 */
PG_FUNCTION_INFO_V1 (compare_dna_ge);
Datum compare_dna_ge(PG_FUNCTION_ARGS)
{
	Varlena* seq1 = (Varlena*) PG_GETARG_RAW_VARLENA_P(0);
	Varlena* seq2 = (Varlena*) PG_GETARG_RAW_VARLENA_P(1);
	bool result = FALSE;

	PB_TRACE(errmsg("->compare_dna_ge()"));

	if (sequence_compare(seq1, seq2, fixed_dna_codes) >= 0)
		result = TRUE;

	PB_TRACE(errmsg("<-compare_dna_ge() exists with %d", result));

	PG_RETURN_BOOL(result);
}

/**
 * compare_dna()
 * 		Compares two DNA sequence for equality.
 *
 * 	Varlena* seq1 : first possibly toasted sequence
 * 	Varlena* seq2 : second possibly toasted sequence
 */
PG_FUNCTION_INFO_V1 (compare_dna);
Datum compare_dna(PG_FUNCTION_ARGS)
{
	Varlena* seq1 = (Varlena*) PG_GETARG_RAW_VARLENA_P(0);
	Varlena* seq2 = (Varlena*) PG_GETARG_RAW_VARLENA_P(1);
	int result;

	PB_TRACE(errmsg("->compare_dna()"));

	result = sequence_compare(seq1, seq2, fixed_dna_codes);

	PB_TRACE(errmsg("<-compare_dna() exists with %d", result));

	PG_RETURN_INT32(result);
}

/**
 * hash_dna()
 * 		Returns a CRC32 for a DNA sequence.
 *
 * 	PB_CompressedSequence* seq1 : input sequence
 */
PG_FUNCTION_INFO_V1 (hash_dna);
Datum hash_dna(PG_FUNCTION_ARGS)
{
	PB_CompressedSequence* seq1 = (PB_CompressedSequence*) PG_GETARG_VARLENA_P(0);
	uint32 result;

	PB_TRACE(errmsg("->hash_dna()"));

	result = sequence_crc32(seq1, fixed_dna_codes);

	PB_TRACE(errmsg("<-hash_dna() exits with %u", result));

	PG_RETURN_UINT32(result);
}

/**
 * strpos_dna()
 * 		Finds the first occurrence of a pattern in a sequence.
 */
PG_FUNCTION_INFO_V1 (strpos_dna);
Datum strpos_dna(PG_FUNCTION_ARGS)
{
	PB_CompressedSequence* seq = (PB_CompressedSequence*) PG_GETARG_VARLENA_P(0);
	text* search = (text*) PG_GETARG_VARLENA_P(1);
	uint32 result;

	PB_TRACE(errmsg("->strpos_dna()"));

	result = sequence_strpos(seq, search, fixed_dna_codes);

	PB_TRACE(errmsg("<-strpos_dna()"));

	PG_RETURN_UINT32(result);
}

/**
 * octet_length_dna()
 * 		Returns byte size of datum.
 */
PG_FUNCTION_INFO_V1 (octet_length_dna);
Datum octet_length_dna(PG_FUNCTION_ARGS)
{
	PG_RETURN_UINT32((uint32) toast_raw_datum_size((Datum) PG_GETARG_RAW_VARLENA_P(0)));
}
