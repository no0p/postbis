/*-------------------------------------------------------------------------
*
* Copyright (c) 2013, Max Planck Institute for Marine Microbiology
*
* This software is released under the PostgreSQL License
*
* Author: Michael Schneider <mschneid@mpi-bremen.de>
*
* IDENTIFICATION
*   src/sequence/code_set_creation
*-------------------------------------------------------------------------
*/
#include <math.h>

#include "postgres.h"

#include "sequence/sequence.h"
#include "sequence/compression.h"
#include "utils/debug.h"

#include "sequence/code_set_creation.h"

/*
 * Local types.
 */

/**
 * Type to store one node of a huffman tree.
 */
typedef struct {
	uint8 zero;		/* child node connected with "zero"-labeled edge */
	uint8 one;		/* child node connected with "one"-labeled edge */
	uint8 symbol;	/* the symbol represented by the leaf node or
						0 if internal node */
	int16 next;		/* linked list to save order of frequency */
	uint32 frequency;	/* frequency */
} PB_HuffmanTreeNode;

/**
 * Stores a huffman tree with up to 256 nodes.
 */
typedef struct {
	uint8 nNodes;
	PB_HuffmanTreeNode nodes[];
} PB_HuffmanTree;

/*
 * Local function declarations
 */
static PB_HuffmanTree* get_huffman_tree(const uint8 n_symbols,
										const uint8* symbols,
										const uint32* frequencies);

static PB_CodeSet* get_huffman_code_dfs(const PB_HuffmanTree* tree,
										const uint8 n_symbols);

/*
 * Local functions
 */

/**
 * get_huffman_tree()
 * 		Build a huffman tree from symbol frequencies.
 *
 * 	uint8 n_symbols : number of symbols
 * 	uint8* symbols : array of symbols, sorted according to frequencies
 * 	uint32* frequencies : symbol frequencies
 */
static PB_HuffmanTree* get_huffman_tree(const uint8 n_symbols,
										const uint8* symbols,
										const uint32* frequencies)
{
	PB_HuffmanTree* tree;

	int next_free_node;
	int min_unvisited_node;

	int i;

	/*
	 * upper boundary for number of nodes in this tree is 2n - 1
	 */
	tree = palloc0(sizeof(PB_HuffmanTree) + (2 * n_symbols - 1) * sizeof(PB_HuffmanTreeNode));

	/*
	 * initially every symbol builds its own one-noded tree
	 */
	for (i = 0; i < n_symbols; i++)
	{
		tree->nodes[i].zero = 0;
		tree->nodes[i].one = 0;
		tree->nodes[i].next = i - 1;
		tree->nodes[i].frequency = frequencies[symbols[i]];
		tree->nodes[i].symbol = symbols[i];

		PB_DEBUG2(errmsg("get_huffman_tree(): i:%d c:%d f:%u",i,symbols[i],frequencies[symbols[i]]));
	}

	/*
	 * no need to start combining nodes if only one node is present
	 */
	if (1 == n_symbols)
	{
		tree->nNodes = 1;
		return tree;
	}

	min_unvisited_node = n_symbols - 1;
	next_free_node = n_symbols;

	PB_DEBUG1(errmsg("get_huffman_tree(): n_symbols=%u min_unvisited_node=%u", n_symbols, min_unvisited_node));

	do
	{
		/*
		 * second minimal node is the next bigger node of the min_unvisited_node
		 */
		int min_2nd_unvisited_node  = tree->nodes[min_unvisited_node].next;

		/*
		 * create parent node with two min_unvisited_node nodes as childs
		 */
		tree->nodes[next_free_node].zero = min_unvisited_node;
		tree->nodes[next_free_node].one = min_2nd_unvisited_node;
		tree->nodes[next_free_node].frequency = tree->nodes[min_unvisited_node].frequency + tree->nodes[min_2nd_unvisited_node].frequency;
		tree->nodes[next_free_node].symbol = 0;

		PB_DEBUG2(errmsg("get_huffman_tree(): combined node %u (%u) and %u (%u) to %d", min_unvisited_node, tree->nodes[min_unvisited_node].frequency, min_2nd_unvisited_node, tree->nodes[min_2nd_unvisited_node].frequency, next_free_node));

		/*
		 * now the minimal node is the one after the 2nd minimal in the sorted list
		 */
		min_unvisited_node = tree->nodes[min_2nd_unvisited_node].next;

		/*
		 * update sorted linked list and insert new inner node
		 */
		if (tree->nodes[next_free_node].frequency < tree->nodes[min_unvisited_node].frequency && min_unvisited_node >= 0)
		{
			/*
			 * if new node is minimal
			 */
			PB_DEBUG2(errmsg("get_huffman_tree(): is minimal, %d %u %d", min_unvisited_node, tree->nodes[min_unvisited_node].frequency, next_free_node));

			tree->nodes[next_free_node].next = (uint8) min_unvisited_node;
			min_unvisited_node = next_free_node;
		}
		else
		{
			/*
			 * new node is not minimal
			 */
			int i = min_unvisited_node;
			int predecessor = min_2nd_unvisited_node;

			/*
			 * go through list (bottom-up) until current node has greater frequency than new node
			 */
			while (tree->nodes[next_free_node].frequency >= tree->nodes[i].frequency && i >= 0)
			{
				PB_DEBUG2(errmsg("get_huffman_tree(): compared node %u (%u) and %d (%u)", next_free_node, tree->nodes[next_free_node].frequency, i, tree->nodes[i].frequency));

				predecessor = i;
				i = tree->nodes[i].next;
			}

			/*
			 * insert new node between current node and its predecessor
			 */
			tree->nodes[next_free_node].next = i;
			tree->nodes[predecessor].next = next_free_node;
		}
		next_free_node++;
	} while (min_unvisited_node >= 0);

	tree->nNodes = next_free_node;

	return tree;
}

/**
 * get_huffman_code_dfs()
 * 		Perform depth-first search through huffman tree to obtain code.
 *
 * 	Only sets n_symbols, symbols and max_codeword_length
 *
 * 	PB_HuffmanTree* tree : huffman tree
 * 	PB_SequenceInfo* info : symbol frequencies
 */
static PB_CodeSet* get_huffman_code_dfs(const PB_HuffmanTree* tree,
										  const uint8 n_symbols)
{
	/*
	 * Stack
	 */
	uint8* stack;				/* points at nodes */
	PB_PrefixCode* stack_code;	/* stores path to nodes */
	uint8* stack_tree_depth;	/* stores length of paths */
	int top;					/* number of elements in stack, also pointer to top */

	int max_tree_depth;			/* maximal depth of tree during DFS */
	int equal_word_length = -1;
	PB_CodeSet* result;

	result = (PB_CodeSet*) palloc0(sizeof(PB_CodeSet) + n_symbols * sizeof(PB_Codeword));
	result->n_symbols = n_symbols;

	stack = palloc0((n_symbols + 2) * sizeof(uint8));
	stack_code = palloc0((n_symbols + 2) * sizeof(PB_PrefixCode));
	stack_tree_depth = palloc0((n_symbols + 2) * sizeof(uint8));

	/*
	 * init stack with root node
	 */
	stack[0] = tree->nNodes - 1;
	stack_code[0] = 0;
	stack_tree_depth[0] = 0;
	top = 0;
	max_tree_depth = 0;

	/*
	 * DFS-Visit
	 */
	do
	{
		/*
		 * pull from stack
		 */
		int cNode = stack[top];
		PB_PrefixCode code = stack_code[top];
		int height = stack_tree_depth[top];
		top--;

		if (tree->nodes[cNode].symbol > 0)
		{
			/*
			 * if current node is leaf node, output code
			 */

			if (height > max_tree_depth)
				max_tree_depth = height;

			result->words[cNode].code = code << (PB_PREFIX_CODE_BIT_SIZE - height);
			result->words[cNode].code_length = height;
			result->words[cNode].symbol = tree->nodes[cNode].symbol;

			/*
			 * Check if length of this code is equal to length of earlier
			 * constructed codes
			 */
			if (equal_word_length == -1 || equal_word_length == height)
				equal_word_length = height;
			else
				equal_word_length = 0;

			PB_DEBUG2(errmsg("get_huffman_code_dfs(): leaf node %d %c top:%d val:%u (%u) %d %u", cNode, tree->nodes[cNode].symbol,top, code, result->words[cNode].code, height, tree->nodes[cNode].frequency));
		}
		else
		{
			/*
			 * if current node is internal node, push childs on stack
			 */
			PB_DEBUG2(errmsg("get_huffman_code_dfs(): internal node %d top:%d val:%u %d", cNode, top, code, height));

			/*
			 * left node
			 */
			top++;
			stack[top] = tree->nodes[cNode].zero;
			stack_code[top] = code << 1;
			stack_tree_depth[top] = height + 1;

			/*
			 * right node
			 */
			top++;
			stack[top] = tree->nodes[cNode].one;
			stack_code[top] = (code << 1) + 1;
			stack_tree_depth[top] = height + 1;
		}
	} while (top >= 0);

	result->max_codeword_length = max_tree_depth;
	result->has_equal_length = equal_word_length == false ? false : true;

	PB_DEBUG1(errmsg("get_huffman_code_dfs(): maxworldlen %u maxlen %lu", result->max_codeword_length, PB_PREFIX_CODE_BIT_SIZE));

	/*
	 * if there are codes that are too long for the PrefixCode type
	 *  return null
	 */
	if (max_tree_depth > PB_PREFIX_CODE_BIT_SIZE)
	{
		pfree(result);
		result = NULL;
	}

	pfree(stack);
	pfree(stack_code);
	pfree(stack_tree_depth);

	return result;
}

/*
 * Public functions.
 */

/**
 * get_equal_lengths_code()
 * 		Creates a code set where all codewords have equal lengths.
 *
 * 	RLE statistics will be ignored.
 *
 * 	PB_SequenceInfo* info : info about sequence at hand
 */
PB_CodeSet* get_equal_lengths_code(const PB_SequenceInfo* info)
{
	PB_CodeSet* result;

	int i;

	PB_TRACE(errmsg("->get_equal_lengths_code()"));

	result = (PB_CodeSet*) palloc0(sizeof(PB_CodeSet) + info->n_symbols * sizeof(PB_Codeword));

	result->n_symbols = info->n_symbols;
	result->max_codeword_length = (int) ceil(log2((double)info->n_symbols));
	result->has_equal_length = true;
	result->ascii_bitmap_high = info->ascii_bitmap_high;
	result->ascii_bitmap_low = info->ascii_bitmap_low;
	result->ignore_case = info->ignore_case;

	PB_DEBUG1(errmsg("get_equal_lengths_code(): equal length = %d", result->max_codeword_length));

	for (i = 0; i < result->n_symbols; i++)
	{
		result->words[i].code = i << (PB_PREFIX_CODE_BIT_SIZE - result->max_codeword_length);
		result->words[i].code_length = result->max_codeword_length;
		result->words[i].symbol = info->symbols[i];

		PB_DEBUG2(errmsg("get_equal_lengths_code(): i:%d word:%u len:%u symbol:%c freq:%u", i, result->words[i].code, result->words[i].code_length, result->words[i].symbol, info->frequencies[(uint8)result->words[i].symbol]));
	}

	PB_TRACE(errmsg("<-get_equal_lengths_code()"));

	return result;
}

/**
 * get_huffman_code()
 * 		Build huffman code from given sequence info.
 *
 * If the huffman tree, that is created during this process, is deeper
 * than the size of PrefixCode allows the function returns NULL.
 *
 * 	PB_SequenceInfo* info : stats of the sequence
 */
PB_CodeSet* get_huffman_code(const PB_SequenceInfo* info)
{
	PB_HuffmanTree* tree;
	PB_CodeSet* result = NULL;

	PB_TRACE(errmsg("->get_huffman_code()"));

	if (info->n_symbols > 0)
	{
		tree = get_huffman_tree(info->n_symbols, info->symbols, info->frequencies);
		result = get_huffman_code_dfs(tree, info->n_symbols);
		pfree(tree);
	}

	if (result)
	{
		result->ascii_bitmap_high = info->ascii_bitmap_high;
		result->ascii_bitmap_low = info->ascii_bitmap_low;
		result->ignore_case = info->ignore_case;
	}

	PB_TRACE(errmsg("<-get_huffman_code()"));

	return result;
}

/**
 * get_huffman_code_rle()
 * 		Build huffman code from given rle sequence info.
 *
 * If the huffman tree, that is created during this process, is deeper
 * than the size of PrefixCode allows the function returns NULL.
 *
 * 	PB_SequenceInfo* info : stats of the sequence
 */
PB_CodeSet* get_huffman_code_rle(const PB_SequenceInfo* info) {
	PB_CodeSet* result = NULL;

	PB_TRACE(errmsg("->get_huffman_code_rle()"));

	if (info->rle_info)
	{
		if (info->n_symbols > 0)
		{
			PB_HuffmanTree* tree;
			tree = get_huffman_tree(info->rle_info->n_symbols,
											info->rle_info->symbols,
											info->rle_info->rle_frequencies);
			result = get_huffman_code_dfs(tree, info->rle_info->n_symbols);
			pfree(tree);
		}

		if (result)
		{
			result->uses_rle = true;
			result->ascii_bitmap_high = info->ascii_bitmap_high;
			result->ascii_bitmap_low = info->ascii_bitmap_low;
			result->ignore_case = info->ignore_case;
			result->has_equal_length = false;
		}
	}

	PB_TRACE(errmsg("<-get_huffman_code_rle()"));

	return result;
}

/**
 * truncate_huffman_code()
 * 		Truncate huffman code if useful.
 *
 * 	If it is advantageous the huffman tree will be truncated. If
 * 	not the result is NULL.
 *
 * 	PB_CodeSet* codeset : original huffman code created from stats
 * 	PB_SequenceInfo* info : stats of input sequence
 */
PB_CodeSet* truncate_huffman_code(const PB_CodeSet* codeset,
								  const PB_SequenceInfo* info)
{
	PB_CodeSet* result = NULL;

	int64 bits_saved;

	const uint32* frequencies;

	int i,j;

	uint8* subtree = palloc(sizeof(uint8) * codeset->n_symbols);
	memset(subtree, 0xFF, codeset->n_symbols);

	PB_TRACE(errmsg("->truncate_huffman_code()"));

	if (codeset->uses_rle)
		frequencies = info->rle_info->rle_frequencies;
	else
		frequencies = info->frequencies;

	/*
	 * try all codes (or leave nodes of the huffman tree) except the one with
	 *   the least symbol frequency
	 */
	for (i = 0; i < codeset->n_symbols - 1; i++)
	{
		/*
		 * length of path from parent node to root
		 *   =  part of code that is equal to "slave" symbol's code
		 */
		const int master_length = (PB_PREFIX_CODE_BIT_SIZE - codeset->words[i].code_length + 1);
		int n_symbols_in_subtree = 0;

		PB_DEBUG2(errmsg("truncate_huffman_code(): checking code %d",i));

		/*
		 * One bit per master symbol can be saved.
		 */
		bits_saved = frequencies[codeset->words[i].symbol];
		/* Minus upper limit for master symbol swap run-length */
		bits_saved -= (frequencies[codeset->words[i].symbol] / PB_MAX_SWAP_RUN_LENGTH) * (PB_SWAP_RUN_LENGTH_BIT_SIZE + 1);
		/* Minus counter at the beginning */
		bits_saved -= PB_SWAP_RUN_LENGTH_BIT_SIZE;

		/*
		 * try all nodes that are in the same sub-tree as the master-symbol node
		 */
		for (j = i + 1; j < codeset->n_symbols; j++)
		{
			if (codeset->words[i].code >> master_length == codeset->words[j].code >> master_length)
			{
				PB_DEBUG2(errmsg("truncate_huffman_code(): %c::%c freq:%u saved:%lu", codeset->words[i].symbol, codeset->words[j].symbol, frequencies[codeset->words[j].symbol],
								frequencies[codeset->words[j].symbol] * PB_SWAP_RUN_LENGTH_BIT_SIZE));

				/*
				 * Minus PB_SWAP_RUN_LENGTH_BIT_SIZE extra bits for swapped symbol characters
				 */
				bits_saved -= frequencies[codeset->words[j].symbol] * PB_SWAP_RUN_LENGTH_BIT_SIZE;

				subtree[n_symbols_in_subtree] = j;
				n_symbols_in_subtree++;
			}
		}

		/*
		 * if we can actually save bits, build new codes and terminate loop
		 */
		if (n_symbols_in_subtree > 0 && bits_saved > 0)
		{
			int max_word_length;
			int symbol_output_pointer;

			/*
			 * one additional symbol, the "swap" symbol or master symbol
			 */
			result = palloc0(sizeof(PB_CodeSet) + (codeset->n_symbols + 1) * sizeof(PB_Codeword));

			result->n_symbols = codeset->n_symbols + 1;
			result->n_swapped_symbols = n_symbols_in_subtree + 1;
			result->ascii_bitmap_high = codeset->ascii_bitmap_high;
			result->ascii_bitmap_low = codeset->ascii_bitmap_low;
			result->ignore_case = codeset->ignore_case;
			result->swap_savings = bits_saved;
			result->uses_rle = codeset->uses_rle;
			result->has_equal_length = false;

			/*
			 * Create master code.
			 */
			max_word_length = 0;
			symbol_output_pointer = 0;
			for (j = 0; j < codeset->n_symbols; j++)
			{
				if (i == j)
				{
					/*
					 * Insert truncated master symbol.
					 */
					PB_PrefixCode truncation_map;
					int truncated_length;

					truncated_length = codeset->words[i].code_length - 1;
					truncation_map = ((PB_PrefixCode) 0xFFFFFFFF) << (PB_PREFIX_CODE_BIT_SIZE - truncated_length);

					result->words[symbol_output_pointer].symbol = codeset->words[i].symbol;
					result->words[symbol_output_pointer].code = codeset->words[i].code & truncation_map;
					result->words[symbol_output_pointer].code_length = truncated_length;

					if (truncated_length > max_word_length)
						max_word_length = result->words[symbol_output_pointer].code_length;

					symbol_output_pointer++;

					PB_DEBUG2(errmsg("truncate_huffman_code(): codeold: %u new:%u lenold:%u",codeset->words[i].code, result->words[symbol_output_pointer].code,truncated_length));
				} else {
					/*
					 * Insert other symbols.
					 */
					if (codeset->words[i].code >> master_length != codeset->words[j].code >> master_length)
					{
						result->words[symbol_output_pointer].symbol = codeset->words[j].symbol;
						result->words[symbol_output_pointer].code = codeset->words[j].code;
						result->words[symbol_output_pointer].code_length = codeset->words[j].code_length;
						if (result->words[symbol_output_pointer].code_length > max_word_length)
							max_word_length = result->words[symbol_output_pointer].code_length;

						symbol_output_pointer++;
					}
				}
			}
			result->max_codeword_length = max_word_length;

			max_word_length = 0;

			/*
			 * First symbol is master symbol.
			 */
			result->words[symbol_output_pointer].symbol = codeset->words[i].symbol;
			result->words[symbol_output_pointer].code_length = 1;
			result->words[symbol_output_pointer].code = codeset->words[i].code << (codeset->words[i].code_length - 1);
			symbol_output_pointer++;
			/*
			 * Insert swapped symbols.
			 */
			for (j = 0; j < n_symbols_in_subtree; j++)
			{
				result->words[symbol_output_pointer].symbol = codeset->words[subtree[j]].symbol;
				result->words[symbol_output_pointer].code_length = codeset->words[subtree[j]].code_length - codeset->words[i].code_length + 1;
				result->words[symbol_output_pointer].code = codeset->words[subtree[j]].code << (codeset->words[i].code_length - 1);
				if (result->words[symbol_output_pointer].code_length > max_word_length)
					max_word_length = result->words[symbol_output_pointer].code_length;

				symbol_output_pointer++;
			}
			result->max_swapped_codeword_length = max_word_length;

			break;
		}
	}

#ifdef DEBUG
	if (result)
	{
		int i;
		int master = result->n_symbols - result->n_swapped_symbols;

		PB_DEBUG1(errmsg("swap master %c (%d)", result->words[master].symbol, result->words[master].symbol));
		PB_DEBUG1(errmsg("swap save %lu", result->swap_savings));
		PB_DEBUG1(errmsg("mastercodelen %u maxwordlen %u", result->n_symbols - result->n_swapped_symbols, result->max_codeword_length));
		for (i = 0; i < result->n_symbols - result->n_swapped_symbols; i++)
		{
			PB_DEBUG2(errmsg("MC %d: %c (%d) - %d len %d", i, result->words[i].symbol, result->words[i].symbol, result->words[i].code, result->words[i].code_length));
		}
		PB_DEBUG1(errmsg("swapcodelen %u maxwordlen %u", result->n_swapped_symbols, result->max_swapped_codeword_length));
		for (i = result->n_symbols - result->n_swapped_symbols; i < result->n_symbols; i++)
		{
			PB_DEBUG2(errmsg("MC %d: %c (%d) - %d len %d", i, result->words[i].symbol, result->words[i].symbol, result->words[i].code, result->words[i].code_length));
		}
	}
#endif

	pfree(subtree);

	PB_TRACE(errmsg("<-truncate_huffman_code()"));

	return result;
}

/**
 * get_optimal_code()
 * 		Creates an optimal code for a given sequence.
 *
 * 	PB_SequenceInfo* info : info about sequence at hand
 */
PB_CodeSet* get_optimal_code(const PB_SequenceInfo* info)
{
	PB_CodeSet* result = NULL;

	PB_TRACE(errmsg("->get_optimal_code()"));

	result = get_huffman_code(info);

	if (!result)
		result = get_equal_lengths_code(info);
	else
	{
		if (info->sequence_length >= PB_MIN_LENGTH_FOR_SWAPPING)
		{
			PB_CodeSet* truncated_code = truncate_huffman_code(result, info);

			if (truncated_code)
			{
				pfree(result);
				result = truncated_code;
			}
		}
	}

	if (info->rle_info)
	{
		PB_CodeSet* rle_code;
		rle_code = get_huffman_code_rle(info);

		if (rle_code && info->sequence_length >= PB_MIN_LENGTH_FOR_SWAPPING)
		{
			PB_CodeSet* truncated_code = truncate_huffman_code(rle_code, info);

			if (truncated_code)
			{
				pfree(rle_code);
				rle_code = truncated_code;
			}
		}

		if (get_compressed_size(info, result) < get_compressed_size(info, rle_code))
			pfree(rle_code);
		else
		{
			pfree(result);
			result = rle_code;
		}
	}

	PB_TRACE(errmsg("<-get_optimal_code()"));

	return result;
}
