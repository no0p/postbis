/*-------------------------------------------------------------------------
*
* Copyright (c) 2013, Max Planck Institute for Marine Microbiology
*
* This software is released under the PostgreSQL License
*
* Author: Michael Schneider <mschneid@mpi-bremen.de>
*
* IDENTIFICATION
*   include/utils/bitarray.h
*
*-------------------------------------------------------------------------
*/
#ifndef UTILS_BITARRAY_H_
#define UTILS_BITARRAY_H_

/**
 * One chunk of a PB_BitArray
 */
typedef unsigned long int PB_BitArrayChunk;

/**
 * A simple bit array implementation
 *
 * Example:
 * 	PB_BitArray ba = malloc(PB_BITARRAY_CALC_MEMSIZE(5000));
 * 	PB_BITARRAY_SET(ba, 1234, 1);
 * 	if (PB_BITARRAY_GET(ba, 1234))
 * 		do something
 * 	free(ba);
 */
typedef PB_BitArrayChunk* PB_BitArray;

/**
 * Size of a PB_BitArrayChunk in bytes
 */
#define PB_BITARRAY_CHUNK_SIZE	((PB_BitArrayChunk) sizeof(PB_BitArrayChunk))

/**
 * Size of a PB_BitArrayChunk in bits
 */
#define PB_BITARRAY_BITS_PER_CHUNK	((PB_BitArrayChunk)(8 * PB_BITARRAY_CHUNK_SIZE))

/**
 * Creates a PB_BitArrayChunk that is all 0 except bit x.
 */
#define PB_BITARRAY_MASK(x)		(((PB_BitArrayChunk)1) << x)

/**
 * Creates a PB_BitArrayChunk that is all 1 except bit x.
 */
#define PB_BITARRAY_NMASK(x)	(~(((PB_BitArrayChunk)1) << x))

/**
 * Calculates the size in memory of a PB_BitArray that takes 'bits' bits.
 */
#define PB_BITARRAY_CALC_MEMSIZE(bits)	(((bits) / PB_BITARRAY_BITS_PER_CHUNK + 1) * PB_BITARRAY_CHUNK_SIZE)

/**
 * Gets the PB_BitArrayChunk that contains bit 'bit'.
 */
#define PB_BITARRAY_CHUNK(bitarray,bit) \
	(((PB_BitArray) (bitarray))[((bit) / PB_BITARRAY_BITS_PER_CHUNK)])

/**
 * Sets a bit in a PB_BitArray.
 */
#define PB_BITARRAY_SET(bitarray, bit, value) \
	PB_BITARRAY_CHUNK(bitarray,bit) = \
	(PB_BITARRAY_CHUNK(bitarray,bit) & \
	PB_BITARRAY_NMASK(bit % PB_BITARRAY_BITS_PER_CHUNK)) | \
	((PB_BitArrayChunk)(value & 1)) << (bit % PB_BITARRAY_BITS_PER_CHUNK)

/**
 * Gets a bit from a PB_BitArray.
 */
#define PB_BITARRAY_GET(bitarray,bit) \
	(PB_BITARRAY_CHUNK(bitarray,bit) & \
	PB_BITARRAY_MASK(bit % PB_BITARRAY_BITS_PER_CHUNK)) >> \
	(bit % PB_BITARRAY_BITS_PER_CHUNK)

#endif /* UTILS_BITARRAY_H_ */
