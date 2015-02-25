#-------------------------------------------------------------------------
#
# Copyright (c) 2013, Max Planck Institute for Marine Microbiology
#
# This software is released under the PostgreSQL License
#
# Author: Michael Schneider <mschneid@mpi-bremen.de>
#
# IDENTIFICATION
#   Makefile
#
#-------------------------------------------------------------------------
OBJS = src/pg_magic.o \
		src/sequence/stats.o \
		src/sequence/code_set_creation.o \
		src/sequence/compression.o \
		src/sequence/generation.o \
		src/sequence/functions.o \
		src/types/dna_sequence.o \
		src/types/rna_sequence.o \
		src/types/aa_sequence.o \
		src/types/aligned_rna_sequence.o \
		src/types/aligned_dna_sequence.o \
		src/types/aligned_aa_sequence.o \
		src/types/alphabet.o \
		src/types/bio_functions.o
MODULE_big = postbis
DATA = sql/postbis--1.0.sql
REGRESS = dna_sequence.test \
		rna_sequence.test \
		aa_sequence.test \
		btree_hash_index.test
REGRESS_OPTS = --inputdir=test
REGRESS_OPTS += --output=test
EXTENSION = postbis
PG_CPPFLAGS=-I./include
#PG_CPPFLAGS+=-D DEBUG
#PG_CPPFLAGS+=-g
#PG_CPPFLAGS+=-O0
PG_CONFIG = pg_config
PGXS := $(shell $(PG_CONFIG) --pgxs)
include $(PGXS)

