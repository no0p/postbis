/*-------------------------------------------------------------------------
*
* Copyright (c) 2013, Max Planck Institute for Marine Microbiology
*
* This software is released under the PostgreSQL License
*
* Author: Michael Schneider <mschneid@mpi-bremen.de>
*
* IDENTIFICATION
*   sql/btree_hash_index.test.sql
*
*-------------------------------------------------------------------------
*/
DROP EXTENSION IF EXISTS postbis CASCADE;
CREATE EXTENSION postbis;

SELECT 'ACGT'::dna_sequence = 'ACGT'::dna_sequence; /* t */
SELECT 'ACGT'::dna_sequence = 'TGCA'::dna_sequence; /* f */
SELECT 'ACGT'::dna_sequence != 'ACGT'::dna_sequence; /* f */
SELECT 'ACGT'::dna_sequence != 'TGCA'::dna_sequence; /* t */
SELECT 'ACGT'::dna_sequence < 'ACGT'::dna_sequence; /* f */
SELECT 'ACGG'::dna_sequence < 'ACGT'::dna_sequence; /* t */
SELECT 'ACGT'::dna_sequence < 'ACGTA'::dna_sequence; /* t */
SELECT 'ACGT'::dna_sequence < 'ACG'::dna_sequence; /* f */
SELECT 'ACGT'::dna_sequence <= 'ACGG'::dna_sequence; /* f */
SELECT 'ACGT'::dna_sequence <= 'ACGT'::dna_sequence; /* t */
SELECT 'ACGT'::dna_sequence <= 'ACTT'::dna_sequence; /* t */
SELECT 'ACGT'::dna_sequence <= 'ACGTA'::dna_sequence; /* t */
SELECT 'ACGT'::dna_sequence <= 'ACG'::dna_sequence; /* f */
SELECT 'ACGT'::dna_sequence > 'ACGT'::dna_sequence; /* f */
SELECT 'ACGG'::dna_sequence > 'ACGT'::dna_sequence; /* f */
SELECT 'ACGT'::dna_sequence > 'ACGTA'::dna_sequence; /* f */
SELECT 'ACGT'::dna_sequence > 'ACG'::dna_sequence; /* t */
SELECT 'ACGT'::dna_sequence >= 'ACGG'::dna_sequence; /* t */
SELECT 'ACGT'::dna_sequence >= 'ACGT'::dna_sequence; /* t */
SELECT 'ACGT'::dna_sequence >= 'ACTT'::dna_sequence; /* f */
SELECT 'ACGT'::dna_sequence >= 'ACGTA'::dna_sequence; /* f */
SELECT 'ACGT'::dna_sequence >= 'ACG'::dna_sequence; /* t */

DROP TABLE IF EXISTS idx_test;
CREATE TABLE idx_test (
  id serial primary key,
  raw_sequence text,
  len int,
  compressed_sequence dna_sequence
);

/* 10,000 DNA sequences between 5,000 and 10,000nts */
INSERT INTO idx_test (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::dna_sequence FROM (
    SELECT generate_sequence(dna_flc(), random_length) AS seq, random_length AS len FROM (
      SELECT (random() * 5000 + 5000)::int AS random_length, generate_series(1, 10000)
    ) AS b
  ) AS a;
INSERT INTO idx_test (raw_sequence, len, compressed_sequence)
  SELECT seq, char_length(seq), seq::dna_sequence FROM (
    SELECT 'AACGCTGACTAGAGCATCAGACT'::text AS seq
  ) AS a;

\timing on

EXPLAIN SELECT id, len FROM idx_test WHERE compressed_sequence = 'AACGCTGACTAGAGCATCAGACT'::dna_sequence;
SELECT id, len FROM idx_test WHERE compressed_sequence = 'AACGCTGACTAGAGCATCAGACT'::dna_sequence;

EXPLAIN SELECT id, len FROM idx_test WHERE compressed_sequence = 'ACGT'::dna_sequence;
SELECT id, len FROM idx_test WHERE compressed_sequence = 'ACGT'::dna_sequence;

CREATE INDEX seq_btree_idx ON idx_test USING btree (compressed_sequence);

EXPLAIN SELECT id, len FROM idx_test WHERE compressed_sequence = 'AACGCTGACTAGAGCATCAGACT'::dna_sequence;
SELECT id, len FROM idx_test WHERE compressed_sequence = 'AACGCTGACTAGAGCATCAGACT'::dna_sequence;

EXPLAIN SELECT id, len FROM idx_test WHERE compressed_sequence = 'ACGT'::dna_sequence;
SELECT id, len FROM idx_test WHERE compressed_sequence = 'ACGT'::dna_sequence;

DROP INDEX seq_btree_idx;

CREATE INDEX seq_hash_idx ON idx_test USING hash (compressed_sequence);

EXPLAIN SELECT id, len FROM idx_test WHERE compressed_sequence = 'AACGCTGACTAGAGCATCAGACT'::dna_sequence;
SELECT id, len FROM idx_test WHERE compressed_sequence = 'AACGCTGACTAGAGCATCAGACT'::dna_sequence;

EXPLAIN SELECT id, len FROM idx_test WHERE compressed_sequence = 'ACGT'::dna_sequence;
SELECT id, len FROM idx_test WHERE compressed_sequence = 'ACGT'::dna_sequence;

DROP INDEX seq_hash_idx;

DROP TABLE idx_test;


