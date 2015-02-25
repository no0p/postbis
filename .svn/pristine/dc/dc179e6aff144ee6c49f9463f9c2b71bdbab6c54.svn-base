/*-------------------------------------------------------------------------
*
* Copyright (c) 2013, Max Planck Institute for Marine Microbiology
*
* This software is released under the PostgreSQL License
*
* Author: Michael Schneider <mschneid@mpi-bremen.de>
*
* IDENTIFICATION
*   sql/rna_sequence.test.sql
*
*-------------------------------------------------------------------------
*/
DROP EXTENSION IF EXISTS postbis CASCADE;
CREATE EXTENSION postbis;

DROP TABLE IF EXISTS rna_sequence_errors;
CREATE TABLE rna_sequence_errors (
  id serial primary key,
  test_set text,
  test_type text,
  raw_sequence text,
  details text
);

/*
* Type modifier combination 1: FLC, CASE_INSENSITIVE
*/
DROP TABLE IF EXISTS rna_sequence_test_flc_ic;

CREATE TABLE rna_sequence_test_flc_ic (
  id serial primary key,
  raw_sequence text,
  len int,
  compressed_sequence rna_sequence(FLC, CASE_INSENSITIVE)
);

/* 1000 RNA sequences with four-letter code and GC-content = 0.5, 10 <= len < 2010 */
INSERT INTO rna_sequence_test_flc_ic (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::rna_sequence(FLC,CASE_INSENSITIVE) FROM (
    SELECT generate_sequence(rna_flc(), random_length) AS seq, random_length AS len FROM (
      SELECT (random() * 2000 + 10)::int AS random_length, generate_series(1, 1000)
    ) AS b
  ) AS a;

/* 100 RNA sequences with four-letter code and GC-content = 0.5, 1000000 <= len < 2000000 */
INSERT INTO rna_sequence_test_flc_ic (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::rna_sequence(FLC,CASE_INSENSITIVE) FROM (
    SELECT generate_sequence(rna_flc(), random_length) AS seq, random_length AS len FROM (
      SELECT (random() * 1000000 + 1000000)::int AS random_length, generate_series(1, 100)
    ) AS b
  ) AS a;

/* 1000 RNA sequences with four-letter code and GC-content = 0.2, 10 <= len < 2010 */
INSERT INTO rna_sequence_test_flc_ic (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::rna_sequence(FLC,CASE_INSENSITIVE) FROM (
    SELECT generate_sequence(rna_flc('{0.1,0.1,0.4,0.4}'), random_length) AS seq, random_length AS len FROM (
      SELECT (random() * 2000 + 10)::int AS random_length, generate_series(1, 1000)
    ) AS b
  ) AS a;

/* 100 RNA sequences with four-letter code and GC-content = 0.2, 1000000 <= len < 2000000 */
INSERT INTO rna_sequence_test_flc_ic (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::rna_sequence(FLC,CASE_INSENSITIVE) FROM (
    SELECT generate_sequence(rna_flc('{0.1,0.1,0.4,0.4}'), random_length) AS seq, random_length AS len FROM (
      SELECT (random() * 1000000 + 1000000)::int AS random_length, generate_series(1, 100)
    ) AS b
  ) AS a;

/* should output error */
INSERT INTO rna_sequence_test_flc_ic (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::rna_sequence(FLC,CASE_INSENSITIVE) FROM (
    SELECT generate_sequence('{B,D}'::alphabet, 100) AS seq, 100 AS len
  ) AS a;

/* full sequence decompression */
INSERT INTO rna_sequence_errors (test_set, test_type, raw_sequence)
  SELECT 'rna_sequence_test_flc_ic' AS test_set,
         'full_sequence_decode' AS test_type,
         seq AS raw_sequence
  FROM (
    SELECT seq FROM (
      SELECT raw_sequence AS seq, compressed_sequence::text = raw_sequence AS result
      FROM rna_sequence_test_flc_ic
    ) AS b
    WHERE result = FALSE
  ) AS a;

/* substr function*/
INSERT INTO rna_sequence_errors (test_set, test_type, raw_sequence, details)
  SELECT 'rna_sequence_test_flc_ic' AS test_set,
         'random_access_sequence_decode' AS test_type,
         seq AS raw_sequence,
         det AS details
  FROM (
    SELECT seq, det FROM (
      SELECT raw_sequence AS seq,
             substr(compressed_sequence, start, substr_len) = substr(raw_sequence, start, substr_len) AS result,
             ('start: ' || start || ' len: ' || substr_len) AS det FROM (
        SELECT compressed_sequence,
               raw_sequence,
               (random() * len * 1.2 - len * 0.1)::int AS start,
               (random() * len * 0.1)::int AS substr_len,
               generate_series(1,100)
        FROM rna_sequence_test_flc_ic
      ) AS c
    ) AS b
    WHERE result = FALSE
  ) AS a;

/* char_length function */
INSERT INTO rna_sequence_errors (test_set, test_type, raw_sequence, details)
  SELECT 'rna_sequence_test_flc_ic' AS test_set,
         'sequence_length' AS test_type,
         seq AS raw_sequence,
         det AS details
  FROM (
    SELECT seq, det FROM (
      SELECT raw_sequence AS seq,
             char_length(compressed_sequence) = len AS result,
             (char_length(compressed_sequence)::text || ' vs ' || len) as det
      FROM rna_sequence_test_flc_ic
    ) AS b
    WHERE result = FALSE
  ) AS a;

/* reverse, complement and reverse_complement functions */
INSERT INTO rna_sequence_errors (test_set, test_type, raw_sequence)
  SELECT 'rna_sequence_test_flc_ic' AS test_set,
         'complement' AS test_type,
         seq AS raw_sequence
  FROM (
    SELECT seq FROM (
      SELECT raw_sequence AS seq, reverse_complement(reverse(complement(compressed_sequence)))::text = raw_sequence AS result
      FROM rna_sequence_test_flc_ic
    ) AS b
    WHERE result = FALSE
  ) AS a;

/* get_alphabet function */
INSERT INTO rna_sequence_errors (test_set, test_type, raw_sequence, details)
  SELECT 'rna_sequence_test_flc_ic' AS test_set,
         'get_alphabet' AS test_type,
         raw_sequence,
         det AS details
  FROM (
    SELECT raw_sequence,
          (get_alphabet(raw_sequence)::text  || get_alphabet(compressed_sequence)::text) AS det,
          get_alphabet(raw_sequence)::text = get_alphabet(compressed_sequence)::text AS result
    FROM rna_sequence_test_flc_ic
  ) AS a
  WHERE result = FALSE;

DROP TABLE rna_sequence_test_flc_ic;

/*
* Type modifier combination 2: FLC, CASE_SENSITIVE
*/
DROP TABLE IF EXISTS rna_sequence_test_flc_cs;
CREATE TABLE rna_sequence_test_flc_cs (
  id serial primary key,
  raw_sequence text,
  len int,
  compressed_sequence rna_sequence(FLC, CASE_SENSITIVE)
);

/* 1000 RNA sequences with four-letter code and GC-content = 0.5, 10 <= len < 2010 */
INSERT INTO rna_sequence_test_flc_cs (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::rna_sequence(FLC,CASE_SENSITIVE) FROM (
    SELECT generate_sequence('{A,C,G,U,a,c,g,u}'::alphabet, random_length) AS seq, random_length AS len FROM (
      SELECT (random() * 2000 + 10)::int AS random_length, generate_series(1, 1000)
    ) AS b
  ) AS a;

/* 100 RNA sequences with four-letter code and GC-content = 0.5, 1000000 <= len < 2000000 */
INSERT INTO rna_sequence_test_flc_cs (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::rna_sequence(FLC,CASE_SENSITIVE) FROM (
    SELECT generate_sequence('{A,C,G,U,a,c,g,u}'::alphabet, random_length) AS seq, random_length AS len FROM (
      SELECT (random() * 1000000 + 1000000)::int AS random_length, generate_series(1, 100)
    ) AS b
  ) AS a;

/* 1000 RNA sequences with four-letter code and GC-content = 0.2, 10 <= len < 2010 */
INSERT INTO rna_sequence_test_flc_cs (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::rna_sequence(FLC,CASE_SENSITIVE) FROM (
    SELECT generate_sequence('{{A,C,G,U,a,c,g,u},{0.05,0.05,0.2,0.2,0.05,0.05,0.2,0.2}}'::alphabet, random_length) AS seq, random_length AS len FROM (
      SELECT (random() * 2000 + 10)::int AS random_length, generate_series(1, 1000)
    ) AS b
  ) AS a;

/* 100 RNA sequences with four-letter code and GC-content = 0.2, 1000000 <= len < 2000000 */
INSERT INTO rna_sequence_test_flc_cs (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::rna_sequence(FLC,CASE_SENSITIVE) FROM (
    SELECT generate_sequence('{{A,C,G,U,a,c,g,u},{0.05,0.05,0.2,0.2,0.05,0.05,0.2,0.2}}'::alphabet, random_length) AS seq, random_length AS len FROM (
      SELECT (random() * 1000000 + 1000000)::int AS random_length, generate_series(1, 100)
    ) AS b
  ) AS a;

/* should output error */
INSERT INTO rna_sequence_test_flc_cs (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::rna_sequence(FLC,CASE_SENSITIVE) FROM (
    SELECT generate_sequence('{B,D}'::alphabet, 100) AS seq, 100 AS len
  ) AS a;

/* full sequence decompression */
INSERT INTO rna_sequence_errors (test_set, test_type, raw_sequence)
  SELECT 'rna_sequence_test_flc_cs' AS test_set,
         'full_sequence_decode' AS test_type,
         seq AS raw_sequence
  FROM (
    SELECT seq FROM (
      SELECT raw_sequence AS seq, compressed_sequence::text = raw_sequence AS result
      FROM rna_sequence_test_flc_cs
    ) AS b
    WHERE result = FALSE
  ) AS a;

/* substr function */
INSERT INTO rna_sequence_errors (test_set, test_type, raw_sequence, details)
  SELECT 'rna_sequence_test_flc_cs' AS test_set,
         'random_access_sequence_decode' AS test_type,
         seq AS raw_sequence,
         det AS details
  FROM (
    SELECT seq, det FROM (
      SELECT raw_sequence AS seq,
             substr(compressed_sequence, start, substr_len) = substr(raw_sequence, start, substr_len) AS result,
             ('start: ' || start || ' len: ' || substr_len) AS det FROM (
        SELECT compressed_sequence,
               raw_sequence,
               (random() * len * 1.2 - len * 0.1)::int AS start,
               (random() * len * 0.1)::int AS substr_len,
               generate_series(1,100)
        FROM rna_sequence_test_flc_cs
      ) AS c
    ) AS b
    WHERE result = FALSE
  ) AS a;

/* char_length function */
INSERT INTO rna_sequence_errors (test_set, test_type, raw_sequence, details)
  SELECT 'rna_sequence_test_flc_cs' AS test_set,
         'sequence_length' AS test_type,
         seq AS raw_sequence,
         det AS details
  FROM (
    SELECT seq, det FROM (
      SELECT raw_sequence AS seq,
             char_length(compressed_sequence) = len AS result,
             (char_length(compressed_sequence)::text || ' vs ' || len) as det
      FROM rna_sequence_test_flc_cs
    ) AS b
    WHERE result = FALSE
  ) AS a;

/* reverse, complement and reverse_complement functions */
INSERT INTO rna_sequence_errors (test_set, test_type, raw_sequence)
  SELECT 'rna_sequence_test_flc_cs' AS test_set,
         'complement' AS test_type,
         seq AS raw_sequence
  FROM (
    SELECT seq FROM (
      SELECT raw_sequence AS seq, reverse_complement(reverse(complement(compressed_sequence)))::text = raw_sequence AS result
      FROM rna_sequence_test_flc_cs
    ) AS b
    WHERE result = FALSE
  ) AS a;

/* get_alphabet function */
INSERT INTO rna_sequence_errors (test_set, test_type, raw_sequence, details)
  SELECT 'rna_sequence_test_flc_cs' AS test_set,
         'get_alphabet' AS test_type,
         raw_sequence,
         det AS details
  FROM (
    SELECT raw_sequence,
          (get_alphabet(raw_sequence)::text  || get_alphabet(compressed_sequence)::text) AS det,
          get_alphabet(raw_sequence)::text = get_alphabet(compressed_sequence)::text AS result
    FROM rna_sequence_test_flc_cs
  ) AS a
  WHERE result = FALSE;

DROP TABLE rna_sequence_test_flc_cs;

/*
* Type modifier combination 3: IUPAC, CASE_INSENSITIVE
*/
DROP TABLE IF EXISTS rna_sequence_test_iupac_ic;
CREATE TABLE rna_sequence_test_iupac_ic (
  id serial primary key,
  raw_sequence text,
  len int,
  compressed_sequence rna_sequence(iupac, CASE_INSENSITIVE)
);

/* 1000 RNA sequences with iupac code and GC-content = 0.5, 10 <= len < 2010 */
INSERT INTO rna_sequence_test_iupac_ic (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::rna_sequence(iupac,CASE_INSENSITIVE) FROM (
    SELECT generate_sequence(rna_iupac(), random_length) AS seq, random_length AS len FROM (
      SELECT (random() * 2000 + 10)::int AS random_length, generate_series(1, 1000)
    ) AS b
  ) AS a;

/* 100 RNA sequences with iupac code and GC-content = 0.5, 1000000 <= len < 2000000 */
INSERT INTO rna_sequence_test_iupac_ic (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::rna_sequence(iupac,CASE_INSENSITIVE) FROM (
    SELECT generate_sequence(rna_iupac(), random_length) AS seq, random_length AS len FROM (
      SELECT (random() * 1000000 + 1000000)::int AS random_length, generate_series(1, 100)
    ) AS b
  ) AS a;

/* 1000 RNA sequences with iupac code and GC-content = 0.2, 10 <= len < 2010 */
INSERT INTO rna_sequence_test_iupac_ic (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::rna_sequence(iupac,CASE_INSENSITIVE) FROM (
    SELECT generate_sequence(rna_iupac('{0.3,0.3,0.1,0.1,0.1,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01}'), random_length) AS seq, random_length AS len FROM (
      SELECT (random() * 2000 + 10)::int AS random_length, generate_series(1, 1000)
    ) AS b
  ) AS a;

/* 100 RNA sequences with iupac code and GC-content = 0.2, 1000000 <= len < 2000000 */
INSERT INTO rna_sequence_test_iupac_ic (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::rna_sequence(iupac,CASE_INSENSITIVE) FROM (
    SELECT generate_sequence(rna_iupac('{0.3,0.3,0.1,0.1,0.1,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01}'), random_length) AS seq, random_length AS len FROM (
      SELECT (random() * 1000000 + 1000000)::int AS random_length, generate_series(1, 100)
    ) AS b
  ) AS a;

/* should output error */
INSERT INTO rna_sequence_test_iupac_ic (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::rna_sequence(iupac,CASE_INSENSITIVE) FROM (
    SELECT generate_sequence('{I,J}'::alphabet, 100) AS seq, 100 AS len
  ) AS a;

/* full sequence decompression */
INSERT INTO rna_sequence_errors (test_set, test_type, raw_sequence)
  SELECT 'rna_sequence_test_iupac_ic' AS test_set,
         'full_sequence_decode' AS test_type,
         seq AS raw_sequence
  FROM (
    SELECT seq FROM (
      SELECT raw_sequence AS seq, compressed_sequence::text = raw_sequence AS result
      FROM rna_sequence_test_iupac_ic
    ) AS b
    WHERE result = FALSE
  ) AS a;

/* substr function */
INSERT INTO rna_sequence_errors (test_set, test_type, raw_sequence, details)
  SELECT 'rna_sequence_test_iupac_ic' AS test_set,
         'random_access_sequence_decode' AS test_type,
         seq AS raw_sequence,
         det AS details
  FROM (
    SELECT seq, det FROM (
      SELECT raw_sequence AS seq,
             substr(compressed_sequence, start, substr_len) = substr(raw_sequence, start, substr_len) AS result,
             ('start: ' || start || ' len: ' || substr_len) AS det FROM (
        SELECT compressed_sequence,
               raw_sequence,
               (random() * len * 1.2 - len * 0.1)::int AS start,
               (random() * len * 0.1)::int AS substr_len,
               generate_series(1,100)
        FROM rna_sequence_test_iupac_ic
      ) AS c
    ) AS b
    WHERE result = FALSE
  ) AS a;

/* char_length function */
INSERT INTO rna_sequence_errors (test_set, test_type, raw_sequence, details)
  SELECT 'rna_sequence_test_iupac_ic' AS test_set,
         'sequence_length' AS test_type,
         seq AS raw_sequence,
         det AS details
  FROM (
    SELECT seq, det FROM (
      SELECT raw_sequence AS seq,
             char_length(compressed_sequence) = len AS result,
             (char_length(compressed_sequence)::text || ' vs ' || len) as det
      FROM rna_sequence_test_iupac_ic
    ) AS b
    WHERE result = FALSE
  ) AS a;

/* complement, reverse and reverse_complement functions */
INSERT INTO rna_sequence_errors (test_set, test_type, raw_sequence)
  SELECT 'rna_sequence_test_iupac_ic' AS test_set,
         'complement' AS test_type,
         seq AS raw_sequence
  FROM (
    SELECT seq FROM (
      SELECT raw_sequence AS seq, reverse_complement(reverse(complement(compressed_sequence)))::text = raw_sequence AS result
      FROM rna_sequence_test_iupac_ic
    ) AS b
    WHERE result = FALSE
  ) AS a;

/* get_alphabet function */
INSERT INTO rna_sequence_errors (test_set, test_type, raw_sequence, details)
  SELECT 'rna_sequence_test_iupac_ic' AS test_set,
         'get_alphabet' AS test_type,
         raw_sequence,
         det AS details
  FROM (
    SELECT raw_sequence,
          (get_alphabet(raw_sequence)::text  || get_alphabet(compressed_sequence)::text) AS det,
          get_alphabet(raw_sequence)::text = get_alphabet(compressed_sequence)::text AS result
    FROM rna_sequence_test_iupac_ic
  ) AS a
  WHERE result = FALSE;

DROP TABLE rna_sequence_test_iupac_ic;

/*
* Type modifier combination 4: IUPAC, CASE_SENSITIVE
*/
DROP TABLE IF EXISTS rna_sequence_test_iupac_cs;
CREATE TABLE rna_sequence_test_iupac_cs (
  id serial primary key,
  raw_sequence text,
  len int,
  compressed_sequence rna_sequence(iupac, CASE_SENSITIVE)
);

/* 1000 RNA sequences with iupac code and GC-content = 0.5, 10 <= len < 2010 */
INSERT INTO rna_sequence_test_iupac_cs (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::rna_sequence(iupac,CASE_SENSITIVE) FROM (
    SELECT generate_sequence('{{A,C,G,U,N,Y,R,M,W,B,V,S,K,D,H,a,c,g,u,n,y,r,m,w,b,v,s,k,d,h},{0.1,0.1,0.1,0.1,0.05,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.1,0.1,0.1,0.1,0.05,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005}}'::alphabet, random_length) AS seq, random_length AS len FROM (
      SELECT (random() * 2000 + 10)::int AS random_length, generate_series(1, 1000)
    ) AS b
  ) AS a;

/* 100 RNA sequences with iupac code and GC-content = 0.5, 1000000 <= len < 2000000 */
INSERT INTO rna_sequence_test_iupac_cs (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::rna_sequence(iupac,CASE_SENSITIVE) FROM (
    SELECT generate_sequence('{{A,C,G,U,N,Y,R,M,W,B,V,S,K,D,H,a,c,g,u,n,y,r,m,w,b,v,s,k,d,h},{0.1,0.1,0.1,0.1,0.05,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.1,0.1,0.1,0.1,0.05,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005}}'::alphabet, random_length) AS seq, random_length AS len FROM (
      SELECT (random() * 1000000 + 1000000)::int AS random_length, generate_series(1, 100)
    ) AS b
  ) AS a;

/* 1000 RNA sequences with iupac code and GC-content = 0.2, 10 <= len < 2010 */
INSERT INTO rna_sequence_test_iupac_cs (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::rna_sequence(iupac,CASE_SENSITIVE) FROM (
    SELECT generate_sequence('{{A,C,G,U,N,Y,R,M,W,B,V,S,K,D,H,a,c,g,u,n,y,r,m,w,b,v,s,k,d,h},{0.05,0.05,0.15,0.15,0.05,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.05,0.05,0.15,0.15,0.05,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005}}'::alphabet, random_length) AS seq, random_length AS len FROM (
      SELECT (random() * 2000 + 10)::int AS random_length, generate_series(1, 1000)
    ) AS b
  ) AS a;

/* 100 RNA sequences with iupac code and GC-content = 0.2, 1000000 <= len < 2000000 */
INSERT INTO rna_sequence_test_iupac_cs (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::rna_sequence(iupac,CASE_SENSITIVE) FROM (
    SELECT generate_sequence('{{A,C,G,U,N,Y,R,M,W,B,V,S,K,D,H,a,c,g,u,n,y,r,m,w,b,v,s,k,d,h},{0.05,0.05,0.15,0.15,0.05,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.05,0.05,0.15,0.15,0.05,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005}}'::alphabet, random_length) AS seq, random_length AS len FROM (
      SELECT (random() * 1000000 + 1000000)::int AS random_length, generate_series(1, 100)
    ) AS b
  ) AS a;

/* should output error */
INSERT INTO rna_sequence_test_iupac_cs (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::rna_sequence(iupac,CASE_SENSITIVE) FROM (
    SELECT generate_sequence('{I,J}'::alphabet, 100) AS seq, 100 AS len
  ) AS a;

/* full sequence decompression */
INSERT INTO rna_sequence_errors (test_set, test_type, raw_sequence)
  SELECT 'rna_sequence_test_iupac_cs' AS test_set,
         'full_sequence_decode' AS test_type,
         seq AS raw_sequence
  FROM (
    SELECT seq FROM (
      SELECT raw_sequence AS seq, compressed_sequence::text = raw_sequence AS result
      FROM rna_sequence_test_iupac_cs
    ) AS b
    WHERE result = FALSE
  ) AS a;

/* substr function */
INSERT INTO rna_sequence_errors (test_set, test_type, raw_sequence, details)
  SELECT 'rna_sequence_test_iupac_cs' AS test_set,
         'random_access_sequence_decode' AS test_type,
         seq AS raw_sequence,
         det AS details
  FROM (
    SELECT seq, det FROM (
      SELECT raw_sequence AS seq,
             substr(compressed_sequence, start, substr_len) = substr(raw_sequence, start, substr_len) AS result,
             ('start: ' || start || ' len: ' || substr_len) AS det FROM (
        SELECT compressed_sequence,
               raw_sequence,
               (random() * len * 1.2 - len * 0.1)::int AS start,
               (random() * len * 0.1)::int AS substr_len,
               generate_series(1,100)
        FROM rna_sequence_test_iupac_cs
      ) AS c
    ) AS b
    WHERE result = FALSE
  ) AS a;

/* char_length function */
INSERT INTO rna_sequence_errors (test_set, test_type, raw_sequence, details)
  SELECT 'rna_sequence_test_iupac_cs' AS test_set,
         'sequence_length' AS test_type,
         seq AS raw_sequence,
         det AS details
  FROM (
    SELECT seq, det FROM (
      SELECT raw_sequence AS seq,
             char_length(compressed_sequence) = len AS result,
             (char_length(compressed_sequence)::text || ' vs ' || len) as det
      FROM rna_sequence_test_iupac_cs
    ) AS b
    WHERE result = FALSE
  ) AS a;

/* reverse, complement and reverse_complement functions */
INSERT INTO rna_sequence_errors (test_set, test_type, raw_sequence)
  SELECT 'rna_sequence_test_iupac_cs' AS test_set,
         'complement' AS test_type,
         seq AS raw_sequence
  FROM (
    SELECT seq FROM (
      SELECT raw_sequence AS seq, reverse_complement(reverse(complement(compressed_sequence)))::text = raw_sequence AS result
      FROM rna_sequence_test_iupac_cs
    ) AS b
    WHERE result = FALSE
  ) AS a;

/* get_alphabet function */
INSERT INTO rna_sequence_errors (test_set, test_type, raw_sequence, details)
  SELECT 'rna_sequence_test_iupac_cs' AS test_set,
         'get_alphabet' AS test_type,
         raw_sequence,
         det AS details
  FROM (
    SELECT raw_sequence,
          (get_alphabet(raw_sequence)::text  || get_alphabet(compressed_sequence)::text) AS det,
          get_alphabet(raw_sequence)::text = get_alphabet(compressed_sequence)::text AS result
    FROM rna_sequence_test_iupac_cs
  ) AS a
  WHERE result = FALSE;

DROP TABLE rna_sequence_test_iupac_cs;

/*
* Type modifier combination 5: ASCII, CASE_INSENSITIVE
*/
DROP TABLE IF EXISTS rna_sequence_test_ascii_ic;
CREATE TABLE rna_sequence_test_ascii_ic (
  id serial primary key,
  raw_sequence text,
  len int,
  compressed_sequence rna_sequence(ascii, CASE_INSENSITIVE)
);

/* 1000 RNA sequences with ascii code and GC-content = 0.5, 10 <= len < 2010 */
INSERT INTO rna_sequence_test_ascii_ic (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::rna_sequence(ascii,CASE_INSENSITIVE) FROM (
    SELECT generate_sequence('{A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z}'::alphabet, random_length) AS seq, random_length AS len FROM (
      SELECT (random() * 2000 + 10)::int AS random_length, generate_series(1, 1000)
    ) AS b
  ) AS a;

/* 100 RNA sequences with ascii code and GC-content = 0.5, 1000000 <= len < 2000000 */
INSERT INTO rna_sequence_test_ascii_ic (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::rna_sequence(ascii,CASE_INSENSITIVE) FROM (
    SELECT generate_sequence('{A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z}'::alphabet, random_length) AS seq, random_length AS len FROM (
      SELECT (random() * 1000000 + 1000000)::int AS random_length, generate_series(1, 100)
    ) AS b
  ) AS a;

/* full sequence decompression */
INSERT INTO rna_sequence_errors (test_set, test_type, raw_sequence)
  SELECT 'rna_sequence_test_ascii_ic' AS test_set,
         'full_sequence_decode' AS test_type,
         seq AS raw_sequence
  FROM (
    SELECT seq FROM (
      SELECT raw_sequence AS seq, compressed_sequence::text = raw_sequence AS result
      FROM rna_sequence_test_ascii_ic
    ) AS b
    WHERE result = FALSE
  ) AS a;

/* substr function */
INSERT INTO rna_sequence_errors (test_set, test_type, raw_sequence, details)
  SELECT 'rna_sequence_test_ascii_ic' AS test_set,
         'random_access_sequence_decode' AS test_type,
         seq AS raw_sequence,
         det AS details
  FROM (
    SELECT seq, det FROM (
      SELECT raw_sequence AS seq,
             substr(compressed_sequence, start, substr_len) = substr(raw_sequence, start, substr_len) AS result,
             ('start: ' || start || ' len: ' || substr_len) AS det FROM (
        SELECT compressed_sequence,
               raw_sequence,
               (random() * len * 1.2 - len * 0.1)::int AS start,
               (random() * len * 0.1)::int AS substr_len,
               generate_series(1,100)
        FROM rna_sequence_test_ascii_ic
      ) AS c
    ) AS b
    WHERE result = FALSE
  ) AS a;

/* char_length function */
INSERT INTO rna_sequence_errors (test_set, test_type, raw_sequence, details)
  SELECT 'rna_sequence_test_ascii_ic' AS test_set,
         'sequence_length' AS test_type,
         seq AS raw_sequence,
         det AS details
  FROM (
    SELECT seq, det FROM (
      SELECT raw_sequence AS seq,
             char_length(compressed_sequence) = len AS result,
             (char_length(compressed_sequence)::text || ' vs ' || len) as det
      FROM rna_sequence_test_ascii_ic
    ) AS b
    WHERE result = FALSE
  ) AS a;

/* complement, reverse and reverse_complement functions */
INSERT INTO rna_sequence_errors (test_set, test_type, raw_sequence)
  SELECT 'rna_sequence_test_ascii_ic' AS test_set,
         'complement' AS test_type,
         seq AS raw_sequence
  FROM (
    SELECT seq FROM (
      SELECT raw_sequence AS seq, reverse_complement(reverse(complement(compressed_sequence)))::text = raw_sequence AS result
      FROM rna_sequence_test_ascii_ic
    ) AS b
    WHERE result = FALSE
  ) AS a;

/* get_alphabet function */
INSERT INTO rna_sequence_errors (test_set, test_type, raw_sequence, details)
  SELECT 'rna_sequence_test_ascii_ic' AS test_set,
         'get_alphabet' AS test_type,
         raw_sequence,
         det AS details
  FROM (
    SELECT raw_sequence,
          (get_alphabet(raw_sequence)::text  || get_alphabet(compressed_sequence)::text) AS det,
          get_alphabet(raw_sequence)::text = get_alphabet(compressed_sequence)::text AS result
    FROM rna_sequence_test_ascii_ic
  ) AS a
  WHERE result = FALSE;

DROP TABLE rna_sequence_test_ascii_ic;

/*
* Type modifier combination 6: ASCII, CASE_SENSITIVE
*/
DROP TABLE IF EXISTS rna_sequence_test_ascii_cs;
CREATE TABLE rna_sequence_test_ascii_cs (
  id serial primary key,
  raw_sequence text,
  len int,
  compressed_sequence rna_sequence(ascii, CASE_SENSITIVE)
);

/* 1000 RNA sequences with ascii code, 10 <= len < 2010 */
INSERT INTO rna_sequence_test_ascii_cs (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::rna_sequence(ascii,CASE_SENSITIVE) FROM (
    SELECT generate_sequence('{A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z}'::alphabet, random_length) AS seq, random_length AS len FROM (
      SELECT (random() * 2000 + 10)::int AS random_length, generate_series(1, 1000)
    ) AS b
  ) AS a;

/* 100 RNA sequences with ascii code, 1000000 <= len < 2000000 */
INSERT INTO rna_sequence_test_ascii_cs (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::rna_sequence(ascii,CASE_SENSITIVE) FROM (
    SELECT generate_sequence('{A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z}'::alphabet, random_length) AS seq, random_length AS len FROM (
      SELECT (random() * 1000000 + 1000000)::int AS random_length, generate_series(1, 100)
    ) AS b
  ) AS a;

/* full sequence decompression */
INSERT INTO rna_sequence_errors (test_set, test_type, raw_sequence)
  SELECT 'rna_sequence_test_ascii_cs' AS test_set,
         'full_sequence_decode' AS test_type,
         seq AS raw_sequence
  FROM (
    SELECT seq FROM (
      SELECT raw_sequence AS seq, compressed_sequence::text = raw_sequence AS result
      FROM rna_sequence_test_ascii_cs
    ) AS b
    WHERE result = FALSE
  ) AS a;

/* substr function */
INSERT INTO rna_sequence_errors (test_set, test_type, raw_sequence, details)
  SELECT 'rna_sequence_test_ascii_cs' AS test_set,
         'random_access_sequence_decode' AS test_type,
         seq AS raw_sequence,
         det AS details
  FROM (
    SELECT seq, det FROM (
      SELECT raw_sequence AS seq,
             substr(compressed_sequence, start, substr_len) = substr(raw_sequence, start, substr_len) AS result,
             ('start: ' || start || ' len: ' || substr_len) AS det FROM (
        SELECT compressed_sequence,
               raw_sequence,
               (random() * len * 1.2 - len * 0.1)::int AS start,
               (random() * len * 0.1)::int AS substr_len,
               generate_series(1,100)
        FROM rna_sequence_test_ascii_cs
      ) AS c
    ) AS b
    WHERE result = FALSE
  ) AS a;

/* char_length function */
INSERT INTO rna_sequence_errors (test_set, test_type, raw_sequence, details)
  SELECT 'rna_sequence_test_ascii_cs' AS test_set,
         'sequence_length' AS test_type,
         seq AS raw_sequence,
         det AS details
  FROM (
    SELECT seq, det FROM (
      SELECT raw_sequence AS seq,
             char_length(compressed_sequence) = len AS result,
             (char_length(compressed_sequence)::text || ' vs ' || len) as det
      FROM rna_sequence_test_ascii_cs
    ) AS b
    WHERE result = FALSE
  ) AS a;

/* complement, reverse and reverse_complement functions */
INSERT INTO rna_sequence_errors (test_set, test_type, raw_sequence)
  SELECT 'rna_sequence_test_ascii_cs' AS test_set,
         'complement' AS test_type,
         seq AS raw_sequence
  FROM (
    SELECT seq FROM (
      SELECT raw_sequence AS seq, reverse_complement(reverse(complement(compressed_sequence)))::text = raw_sequence AS result
      FROM rna_sequence_test_ascii_cs
    ) AS b
    WHERE result = FALSE
  ) AS a;

/* get_alphabet function */
INSERT INTO rna_sequence_errors (test_set, test_type, raw_sequence, details)
  SELECT 'rna_sequence_test_ascii_cs' AS test_set,
         'get_alphabet' AS test_type,
         raw_sequence,
         det AS details
  FROM (
    SELECT raw_sequence,
          (get_alphabet(raw_sequence)::text  || get_alphabet(compressed_sequence)::text) AS det,
          get_alphabet(raw_sequence)::text = get_alphabet(compressed_sequence)::text AS result
    FROM rna_sequence_test_ascii_cs
  ) AS a
  WHERE result = FALSE;

DROP TABLE rna_sequence_test_ascii_cs;

SELECT test_set, test_type, count(*) FROM rna_sequence_errors GROUP BY test_set, test_type;

/*
*SET client_min_messages=DEBUG1;
*SELECT show_counts();
*/

DROP EXTENSION postbis;