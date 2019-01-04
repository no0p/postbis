/*-------------------------------------------------------------------------
*
* Copyright (c) 2013, Max Planck Institute for Marine Microbiology
*
* This software is released under the PostgreSQL License
*
* Author: Michael Schneider <mschneid@mpi-bremen.de>
*
* IDENTIFICATION
*   sql/dna_sequence.test.sql
*
*-------------------------------------------------------------------------
*/
DROP EXTENSION IF EXISTS postbis CASCADE;
CREATE EXTENSION postbis;

DROP TABLE IF EXISTS dna_sequence_errors;
CREATE TABLE dna_sequence_errors (
  id serial primary key,
  test_set text,
  test_type text,
  raw_sequence text,
  details text
);

/*
* Type modifier combination 1: SHORT, FLC, CASE_INSENSITIVE
*/
DROP TABLE IF EXISTS dna_sequence_test_short_flc_ic;

CREATE TABLE dna_sequence_test_short_flc_ic (
  id serial primary key,
  raw_sequence text,
  len int,
  compressed_sequence dna_sequence(SHORT, FLC, CASE_INSENSITIVE)
);

/* 1000 DNA sequences with four-letter code and GC-content = 0.5, 10 <= len < 2010 */
INSERT INTO dna_sequence_test_short_flc_ic (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::dna_sequence(SHORT,FLC,CASE_INSENSITIVE) FROM (
    SELECT generate_sequence(dna_flc(), random_length) AS seq, random_length AS len FROM (
      SELECT (random() * 2000 + 10)::int AS random_length, generate_series(1, 1000)
    ) AS b
  ) AS a;

/* 100 DNA sequences with four-letter code and GC-content = 0.5, 1000000 <= len < 2000000 */
INSERT INTO dna_sequence_test_short_flc_ic (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::dna_sequence(SHORT,FLC,CASE_INSENSITIVE) FROM (
    SELECT generate_sequence(dna_flc(), random_length) AS seq, random_length AS len FROM (
      SELECT (random() * 1000000 + 1000000)::int AS random_length, generate_series(1, 100)
    ) AS b
  ) AS a;

/* 1000 DNA sequences with four-letter code and GC-content = 0.2, 10 <= len < 2010 */
INSERT INTO dna_sequence_test_short_flc_ic (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::dna_sequence(SHORT,FLC,CASE_INSENSITIVE) FROM (
    SELECT generate_sequence(dna_flc('{0.1,0.1,0.4,0.4}'), random_length) AS seq, random_length AS len FROM (
      SELECT (random() * 2000 + 10)::int AS random_length, generate_series(1, 1000)
    ) AS b
  ) AS a;

/* 100 DNA sequences with four-letter code and GC-content = 0.2, 1000000 <= len < 2000000 */
INSERT INTO dna_sequence_test_short_flc_ic (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::dna_sequence(SHORT,FLC,CASE_INSENSITIVE) FROM (
    SELECT generate_sequence(dna_flc('{0.1,0.1,0.4,0.4}'), random_length) AS seq, random_length AS len FROM (
      SELECT (random() * 1000000 + 1000000)::int AS random_length, generate_series(1, 100)
    ) AS b
  ) AS a;

/* should output error */
INSERT INTO dna_sequence_test_short_flc_ic (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::dna_sequence(SHORT,FLC,CASE_INSENSITIVE) FROM (
    SELECT generate_sequence('{B,D}'::alphabet, 100) AS seq, 100 AS len
  ) AS a;

/* full sequence decompression */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence)
  SELECT 'dna_sequence_test_short_flc_ic' AS test_set,
         'full_sequence_decode' AS test_type,
         seq AS raw_sequence
  FROM (
    SELECT seq FROM (
      SELECT raw_sequence AS seq, compressed_sequence::text = raw_sequence AS result
      FROM dna_sequence_test_short_flc_ic
    ) AS b
    WHERE result = false
  ) AS a;

/* substr function*/
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence, details)
  SELECT 'dna_sequence_test_short_flc_ic' AS test_set,
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
        FROM dna_sequence_test_short_flc_ic
      ) AS c
    ) AS b
    WHERE result = false
  ) AS a;

/* char_length function */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence, details)
  SELECT 'dna_sequence_test_short_flc_ic' AS test_set,
         'sequence_length' AS test_type,
         seq AS raw_sequence,
         det AS details
  FROM (
    SELECT seq, det FROM (
      SELECT raw_sequence AS seq,
             char_length(compressed_sequence) = len AS result,
             (char_length(compressed_sequence)::text || ' vs ' || len) as det
      FROM dna_sequence_test_short_flc_ic
    ) AS b
    WHERE result = false
  ) AS a;

/* reverse, complement and reverse_complement functions */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence)
  SELECT 'dna_sequence_test_short_flc_ic' AS test_set,
         'complement' AS test_type,
         seq AS raw_sequence
  FROM (
    SELECT seq FROM (
      SELECT raw_sequence AS seq, reverse_complement(reverse(complement(compressed_sequence)))::text = raw_sequence AS result
      FROM dna_sequence_test_short_flc_ic
    ) AS b
    WHERE result = false
  ) AS a;

/* get_alphabet function */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence, details)
  SELECT 'dna_sequence_test_short_flc_ic' AS test_set,
         'get_alphabet' AS test_type,
         raw_sequence,
         det AS details
  FROM (
    SELECT raw_sequence,
          (get_alphabet(raw_sequence)::text  || get_alphabet(compressed_sequence)::text) AS det,
          get_alphabet(raw_sequence)::text = get_alphabet(compressed_sequence)::text AS result
    FROM dna_sequence_test_short_flc_ic
  ) AS a
  WHERE result = false;

/* transcribe, reverse_transcribe functions */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence)
  SELECT 'dna_sequence_test_short_flc_ic' AS test_set,
         'transcription' AS test_type,
         seq AS raw_sequence
  FROM (
    SELECT seq FROM (
      SELECT raw_sequence AS seq, reverse_transcribe(transcribe(compressed_sequence))::text = raw_sequence AS result
      FROM dna_sequence_test_short_flc_ic
    ) AS b
    WHERE result = false
  ) AS a;

/* strpos function */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence)
  SELECT 'dna_sequence_test_short_flc_ic' AS test_set,
         'string matching' AS test_type,
         seq AS raw_sequence
  FROM (
    SELECT seq FROM (
      SELECT raw_sequence AS seq, (strpos(compressed_sequence, 'ACGTA') = strpos(raw_sequence, 'ACGTA')) AS result
      FROM dna_sequence_test_short_flc_ic
    ) AS b
    WHERE result = false
  ) AS a;

DROP TABLE dna_sequence_test_short_flc_ic;

/*
* Type modifier combination 2: SHORT, FLC, CASE_SENSITIVE
*/
DROP TABLE IF EXISTS dna_sequence_test_short_flc_cs;
CREATE TABLE dna_sequence_test_short_flc_cs (
  id serial primary key,
  raw_sequence text,
  len int,
  compressed_sequence dna_sequence(SHORT, FLC, CASE_SENSITIVE)
);

/* 1000 DNA sequences with four-letter code and GC-content = 0.5, 10 <= len < 2010 */
INSERT INTO dna_sequence_test_short_flc_cs (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::dna_sequence(SHORT,FLC,CASE_SENSITIVE) FROM (
    SELECT generate_sequence('{A,C,G,T,a,c,g,t}'::alphabet, random_length) AS seq, random_length AS len FROM (
      SELECT (random() * 2000 + 10)::int AS random_length, generate_series(1, 1000)
    ) AS b
  ) AS a;

/* 100 DNA sequences with four-letter code and GC-content = 0.5, 1000000 <= len < 2000000 */
INSERT INTO dna_sequence_test_short_flc_cs (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::dna_sequence(SHORT,FLC,CASE_SENSITIVE) FROM (
    SELECT generate_sequence('{A,C,G,T,a,c,g,t}'::alphabet, random_length) AS seq, random_length AS len FROM (
      SELECT (random() * 1000000 + 1000000)::int AS random_length, generate_series(1, 100)
    ) AS b
  ) AS a;

/* 1000 DNA sequences with four-letter code and GC-content = 0.2, 10 <= len < 2010 */
INSERT INTO dna_sequence_test_short_flc_cs (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::dna_sequence(SHORT,FLC,CASE_SENSITIVE) FROM (
    SELECT generate_sequence('{{A,C,G,T,a,c,g,t},{0.05,0.05,0.2,0.2,0.05,0.05,0.2,0.2}}'::alphabet, random_length) AS seq, random_length AS len FROM (
      SELECT (random() * 2000 + 10)::int AS random_length, generate_series(1, 1000)
    ) AS b
  ) AS a;

/* 100 DNA sequences with four-letter code and GC-content = 0.2, 1000000 <= len < 2000000 */
INSERT INTO dna_sequence_test_short_flc_cs (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::dna_sequence(SHORT,FLC,CASE_SENSITIVE) FROM (
    SELECT generate_sequence('{{A,C,G,T,a,c,g,t},{0.05,0.05,0.2,0.2,0.05,0.05,0.2,0.2}}'::alphabet, random_length) AS seq, random_length AS len FROM (
      SELECT (random() * 1000000 + 1000000)::int AS random_length, generate_series(1, 100)
    ) AS b
  ) AS a;

/* should output error */
INSERT INTO dna_sequence_test_short_flc_cs (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::dna_sequence(SHORT,FLC,CASE_SENSITIVE) FROM (
    SELECT generate_sequence('{B,D}'::alphabet, 100) AS seq, 100 AS len
  ) AS a;

/* full sequence decompression */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence)
  SELECT 'dna_sequence_test_short_flc_cs' AS test_set,
         'full_sequence_decode' AS test_type,
         seq AS raw_sequence
  FROM (
    SELECT seq FROM (
      SELECT raw_sequence AS seq, compressed_sequence::text = raw_sequence AS result
      FROM dna_sequence_test_short_flc_cs
    ) AS b
    WHERE result = false
  ) AS a;

/* substr function */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence, details)
  SELECT 'dna_sequence_test_short_flc_cs' AS test_set,
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
        FROM dna_sequence_test_short_flc_cs
      ) AS c
    ) AS b
    WHERE result = false
  ) AS a;

/* char_length function */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence, details)
  SELECT 'dna_sequence_test_short_flc_cs' AS test_set,
         'sequence_length' AS test_type,
         seq AS raw_sequence,
         det AS details
  FROM (
    SELECT seq, det FROM (
      SELECT raw_sequence AS seq,
             char_length(compressed_sequence) = len AS result,
             (char_length(compressed_sequence)::text || ' vs ' || len) as det
      FROM dna_sequence_test_short_flc_cs
    ) AS b
    WHERE result = false
  ) AS a;

/* reverse, complement and reverse_complement functions */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence)
  SELECT 'dna_sequence_test_short_flc_cs' AS test_set,
         'complement' AS test_type,
         seq AS raw_sequence
  FROM (
    SELECT seq FROM (
      SELECT raw_sequence AS seq, reverse_complement(reverse(complement(compressed_sequence)))::text = raw_sequence AS result
      FROM dna_sequence_test_short_flc_cs
    ) AS b
    WHERE result = false
  ) AS a;

/* get_alphabet function */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence, details)
  SELECT 'dna_sequence_test_short_flc_cs' AS test_set,
         'get_alphabet' AS test_type,
         raw_sequence,
         det AS details
  FROM (
    SELECT raw_sequence,
          (get_alphabet(raw_sequence)::text  || get_alphabet(compressed_sequence)::text) AS det,
          get_alphabet(raw_sequence)::text = get_alphabet(compressed_sequence)::text AS result
    FROM dna_sequence_test_short_flc_cs
  ) AS a
  WHERE result = false;

/* transcribe, reverse_transcribe functions */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence)
  SELECT 'dna_sequence_test_short_flc_cs' AS test_set,
         'transcription' AS test_type,
         seq AS raw_sequence
  FROM (
    SELECT seq FROM (
      SELECT raw_sequence AS seq, reverse_transcribe(transcribe(compressed_sequence))::text = raw_sequence AS result
      FROM dna_sequence_test_short_flc_cs
    ) AS b
    WHERE result = false
  ) AS a;

/* strpos function */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence)
  SELECT 'dna_sequence_test_short_flc_cs' AS test_set,
         'string matching' AS test_type,
         seq AS raw_sequence
  FROM (
    SELECT seq FROM (
      SELECT raw_sequence AS seq, (strpos(compressed_sequence, 'ACGTA') = strpos(raw_sequence, 'ACGTA')) AS result
      FROM dna_sequence_test_short_flc_cs
    ) AS b
    WHERE result = false
  ) AS a;

DROP TABLE dna_sequence_test_short_flc_cs;

/*
* Type modifier combination 3: SHORT, IUPAC, CASE_INSENSITIVE
*/
DROP TABLE IF EXISTS dna_sequence_test_short_iupac_ic;
CREATE TABLE dna_sequence_test_short_iupac_ic (
  id serial primary key,
  raw_sequence text,
  len int,
  compressed_sequence dna_sequence(SHORT, iupac, CASE_INSENSITIVE)
);

/* 1000 DNA sequences with iupac code and GC-content = 0.5, 10 <= len < 2010 */
INSERT INTO dna_sequence_test_short_iupac_ic (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::dna_sequence(SHORT,iupac,CASE_INSENSITIVE) FROM (
    SELECT generate_sequence(dna_iupac(), random_length) AS seq, random_length AS len FROM (
      SELECT (random() * 2000 + 10)::int AS random_length, generate_series(1, 1000)
    ) AS b
  ) AS a;

/* 100 DNA sequences with iupac code and GC-content = 0.5, 1000000 <= len < 2000000 */
INSERT INTO dna_sequence_test_short_iupac_ic (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::dna_sequence(SHORT,iupac,CASE_INSENSITIVE) FROM (
    SELECT generate_sequence(dna_iupac(), random_length) AS seq, random_length AS len FROM (
      SELECT (random() * 1000000 + 1000000)::int AS random_length, generate_series(1, 100)
    ) AS b
  ) AS a;

/* 1000 DNA sequences with iupac code and GC-content = 0.2, 10 <= len < 2010 */
INSERT INTO dna_sequence_test_short_iupac_ic (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::dna_sequence(SHORT,iupac,CASE_INSENSITIVE) FROM (
    SELECT generate_sequence(dna_iupac('{0.3,0.3,0.1,0.1,0.1,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01}'), random_length) AS seq, random_length AS len FROM (
      SELECT (random() * 2000 + 10)::int AS random_length, generate_series(1, 1000)
    ) AS b
  ) AS a;

/* 100 DNA sequences with iupac code and GC-content = 0.2, 1000000 <= len < 2000000 */
INSERT INTO dna_sequence_test_short_iupac_ic (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::dna_sequence(SHORT,iupac,CASE_INSENSITIVE) FROM (
    SELECT generate_sequence(dna_iupac('{0.3,0.3,0.1,0.1,0.1,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01}'), random_length) AS seq, random_length AS len FROM (
      SELECT (random() * 1000000 + 1000000)::int AS random_length, generate_series(1, 100)
    ) AS b
  ) AS a;

/* should output error */
INSERT INTO dna_sequence_test_short_iupac_ic (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::dna_sequence(SHORT,iupac,CASE_INSENSITIVE) FROM (
    SELECT generate_sequence('{I,J}'::alphabet, 100) AS seq, 100 AS len
  ) AS a;

/* full sequence decompression */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence)
  SELECT 'dna_sequence_test_short_iupac_ic' AS test_set,
         'full_sequence_decode' AS test_type,
         seq AS raw_sequence
  FROM (
    SELECT seq FROM (
      SELECT raw_sequence AS seq, compressed_sequence::text = raw_sequence AS result
      FROM dna_sequence_test_short_iupac_ic
    ) AS b
    WHERE result = false
  ) AS a;

/* substr function */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence, details)
  SELECT 'dna_sequence_test_short_iupac_ic' AS test_set,
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
        FROM dna_sequence_test_short_iupac_ic
      ) AS c
    ) AS b
    WHERE result = false
  ) AS a;

/* char_length function */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence, details)
  SELECT 'dna_sequence_test_short_iupac_ic' AS test_set,
         'sequence_length' AS test_type,
         seq AS raw_sequence,
         det AS details
  FROM (
    SELECT seq, det FROM (
      SELECT raw_sequence AS seq,
             char_length(compressed_sequence) = len AS result,
             (char_length(compressed_sequence)::text || ' vs ' || len) as det
      FROM dna_sequence_test_short_iupac_ic
    ) AS b
    WHERE result = false
  ) AS a;

/* complement, reverse and reverse_complement functions */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence)
  SELECT 'dna_sequence_test_short_iupac_ic' AS test_set,
         'complement' AS test_type,
         seq AS raw_sequence
  FROM (
    SELECT seq FROM (
      SELECT raw_sequence AS seq, reverse_complement(reverse(complement(compressed_sequence)))::text = raw_sequence AS result
      FROM dna_sequence_test_short_iupac_ic
    ) AS b
    WHERE result = false
  ) AS a;

/* get_alphabet function */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence, details)
  SELECT 'dna_sequence_test_short_iupac_ic' AS test_set,
         'get_alphabet' AS test_type,
         raw_sequence,
         det AS details
  FROM (
    SELECT raw_sequence,
          (get_alphabet(raw_sequence)::text  || get_alphabet(compressed_sequence)::text) AS det,
          get_alphabet(raw_sequence)::text = get_alphabet(compressed_sequence)::text AS result
    FROM dna_sequence_test_short_iupac_ic
  ) AS a
  WHERE result = false;

/* transcribe, reverse_transcribe functions */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence)
  SELECT 'dna_sequence_test_short_iupac_ic' AS test_set,
         'transcription' AS test_type,
         seq AS raw_sequence
  FROM (
    SELECT seq FROM (
      SELECT raw_sequence AS seq, reverse_transcribe(transcribe(compressed_sequence))::text = raw_sequence AS result
      FROM dna_sequence_test_short_iupac_ic
    ) AS b
    WHERE result = false
  ) AS a;

/* strpos function */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence)
  SELECT 'dna_sequence_test_short_iupac_ic' AS test_set,
         'string matching' AS test_type,
         seq AS raw_sequence
  FROM (
    SELECT seq FROM (
      SELECT raw_sequence AS seq, (strpos(compressed_sequence, 'ACGTA') = strpos(raw_sequence, 'ACGTA')) AS result
      FROM dna_sequence_test_short_iupac_ic
    ) AS b
    WHERE result = false
  ) AS a;

DROP TABLE dna_sequence_test_short_iupac_ic;

/*
* Type modifier combination 4: SHORT, IUPAC, CASE_SENSITIVE
*/
DROP TABLE IF EXISTS dna_sequence_test_short_iupac_cs;
CREATE TABLE dna_sequence_test_short_iupac_cs (
  id serial primary key,
  raw_sequence text,
  len int,
  compressed_sequence dna_sequence(SHORT, iupac, CASE_SENSITIVE)
);

/* 1000 DNA sequences with iupac code and GC-content = 0.5, 10 <= len < 2010 */
INSERT INTO dna_sequence_test_short_iupac_cs (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::dna_sequence(SHORT,iupac,CASE_SENSITIVE) FROM (
    SELECT generate_sequence('{{A,C,G,T,N,Y,R,M,W,B,V,S,K,D,H,a,c,g,t,n,y,r,m,w,b,v,s,k,d,h},{0.1,0.1,0.1,0.1,0.05,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.1,0.1,0.1,0.1,0.05,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005}}'::alphabet, random_length) AS seq, random_length AS len FROM (
      SELECT (random() * 2000 + 10)::int AS random_length, generate_series(1, 1000)
    ) AS b
  ) AS a;

/* 100 DNA sequences with iupac code and GC-content = 0.5, 1000000 <= len < 2000000 */
INSERT INTO dna_sequence_test_short_iupac_cs (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::dna_sequence(SHORT,iupac,CASE_SENSITIVE) FROM (
    SELECT generate_sequence('{{A,C,G,T,N,Y,R,M,W,B,V,S,K,D,H,a,c,g,t,n,y,r,m,w,b,v,s,k,d,h},{0.1,0.1,0.1,0.1,0.05,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.1,0.1,0.1,0.1,0.05,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005}}'::alphabet, random_length) AS seq, random_length AS len FROM (
      SELECT (random() * 1000000 + 1000000)::int AS random_length, generate_series(1, 100)
    ) AS b
  ) AS a;

/* 1000 DNA sequences with iupac code and GC-content = 0.2, 10 <= len < 2010 */
INSERT INTO dna_sequence_test_short_iupac_cs (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::dna_sequence(SHORT,iupac,CASE_SENSITIVE) FROM (
    SELECT generate_sequence('{{A,C,G,T,N,Y,R,M,W,B,V,S,K,D,H,a,c,g,t,n,y,r,m,w,b,v,s,k,d,h},{0.05,0.05,0.15,0.15,0.05,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.05,0.05,0.15,0.15,0.05,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005}}'::alphabet, random_length) AS seq, random_length AS len FROM (
      SELECT (random() * 2000 + 10)::int AS random_length, generate_series(1, 1000)
    ) AS b
  ) AS a;

/* 100 DNA sequences with iupac code and GC-content = 0.2, 1000000 <= len < 2000000 */
INSERT INTO dna_sequence_test_short_iupac_cs (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::dna_sequence(SHORT,iupac,CASE_SENSITIVE) FROM (
    SELECT generate_sequence('{{A,C,G,T,N,Y,R,M,W,B,V,S,K,D,H,a,c,g,t,n,y,r,m,w,b,v,s,k,d,h},{0.05,0.05,0.15,0.15,0.05,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.05,0.05,0.15,0.15,0.05,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005}}'::alphabet, random_length) AS seq, random_length AS len FROM (
      SELECT (random() * 1000000 + 1000000)::int AS random_length, generate_series(1, 100)
    ) AS b
  ) AS a;

/* should output error */
INSERT INTO dna_sequence_test_short_iupac_cs (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::dna_sequence(SHORT,iupac,CASE_SENSITIVE) FROM (
    SELECT generate_sequence('{I,J}'::alphabet, 100) AS seq, 100 AS len
  ) AS a;

/* full sequence decompression */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence)
  SELECT 'dna_sequence_test_short_iupac_cs' AS test_set,
         'full_sequence_decode' AS test_type,
         seq AS raw_sequence
  FROM (
    SELECT seq FROM (
      SELECT raw_sequence AS seq, compressed_sequence::text = raw_sequence AS result
      FROM dna_sequence_test_short_iupac_cs
    ) AS b
    WHERE result = false
  ) AS a;

/* substr function */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence, details)
  SELECT 'dna_sequence_test_short_iupac_cs' AS test_set,
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
        FROM dna_sequence_test_short_iupac_cs
      ) AS c
    ) AS b
    WHERE result = false
  ) AS a;

/* char_length function */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence, details)
  SELECT 'dna_sequence_test_short_iupac_cs' AS test_set,
         'sequence_length' AS test_type,
         seq AS raw_sequence,
         det AS details
  FROM (
    SELECT seq, det FROM (
      SELECT raw_sequence AS seq,
             char_length(compressed_sequence) = len AS result,
             (char_length(compressed_sequence)::text || ' vs ' || len) as det
      FROM dna_sequence_test_short_iupac_cs
    ) AS b
    WHERE result = false
  ) AS a;

/* reverse, complement and reverse_complement functions */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence)
  SELECT 'dna_sequence_test_short_iupac_cs' AS test_set,
         'complement' AS test_type,
         seq AS raw_sequence
  FROM (
    SELECT seq FROM (
      SELECT raw_sequence AS seq, reverse_complement(reverse(complement(compressed_sequence)))::text = raw_sequence AS result
      FROM dna_sequence_test_short_iupac_cs
    ) AS b
    WHERE result = false
  ) AS a;

/* get_alphabet function */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence, details)
  SELECT 'dna_sequence_test_short_iupac_cs' AS test_set,
         'get_alphabet' AS test_type,
         raw_sequence,
         det AS details
  FROM (
    SELECT raw_sequence,
          (get_alphabet(raw_sequence)::text  || get_alphabet(compressed_sequence)::text) AS det,
          get_alphabet(raw_sequence)::text = get_alphabet(compressed_sequence)::text AS result
    FROM dna_sequence_test_short_iupac_cs
  ) AS a
  WHERE result = false;

/* transcribe, reverse_transcribe functions */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence)
  SELECT 'dna_sequence_test_short_iupac_cs' AS test_set,
         'transcription' AS test_type,
         seq AS raw_sequence
  FROM (
    SELECT seq FROM (
      SELECT raw_sequence AS seq, reverse_transcribe(transcribe(compressed_sequence))::text = raw_sequence AS result
      FROM dna_sequence_test_short_iupac_cs
    ) AS b
    WHERE result = false
  ) AS a;

/* strpos function */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence)
  SELECT 'dna_sequence_test_short_iupac_cs' AS test_set,
         'string matching' AS test_type,
         seq AS raw_sequence
  FROM (
    SELECT seq FROM (
      SELECT raw_sequence AS seq, (strpos(compressed_sequence, 'ACGTA') = strpos(raw_sequence, 'ACGTA')) AS result
      FROM dna_sequence_test_short_iupac_cs
    ) AS b
    WHERE result = false
  ) AS a;

DROP TABLE dna_sequence_test_short_iupac_cs;

/*
* Type modifier combination 5: SHORT, ASCII, CASE_INSENSITIVE
*/
DROP TABLE IF EXISTS dna_sequence_test_short_ascii_ic;
CREATE TABLE dna_sequence_test_short_ascii_ic (
  id serial primary key,
  raw_sequence text,
  len int,
  compressed_sequence dna_sequence(SHORT, ascii, CASE_INSENSITIVE)
);

/* 1000 DNA sequences with ascii code and GC-content = 0.5, 10 <= len < 2010 */
INSERT INTO dna_sequence_test_short_ascii_ic (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::dna_sequence(SHORT,ascii,CASE_INSENSITIVE) FROM (
    SELECT generate_sequence('{A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z}'::alphabet, random_length) AS seq, random_length AS len FROM (
      SELECT (random() * 2000 + 10)::int AS random_length, generate_series(1, 1000)
    ) AS b
  ) AS a;

/* 100 DNA sequences with ascii code and GC-content = 0.5, 1000000 <= len < 2000000 */
INSERT INTO dna_sequence_test_short_ascii_ic (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::dna_sequence(SHORT,ascii,CASE_INSENSITIVE) FROM (
    SELECT generate_sequence('{A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z}'::alphabet, random_length) AS seq, random_length AS len FROM (
      SELECT (random() * 1000000 + 1000000)::int AS random_length, generate_series(1, 100)
    ) AS b
  ) AS a;

/* full sequence decompression */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence)
  SELECT 'dna_sequence_test_short_ascii_ic' AS test_set,
         'full_sequence_decode' AS test_type,
         seq AS raw_sequence
  FROM (
    SELECT seq FROM (
      SELECT raw_sequence AS seq, compressed_sequence::text = raw_sequence AS result
      FROM dna_sequence_test_short_ascii_ic
    ) AS b
    WHERE result = false
  ) AS a;

/* substr function */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence, details)
  SELECT 'dna_sequence_test_short_ascii_ic' AS test_set,
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
        FROM dna_sequence_test_short_ascii_ic
      ) AS c
    ) AS b
    WHERE result = false
  ) AS a;

/* char_length function */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence, details)
  SELECT 'dna_sequence_test_short_ascii_ic' AS test_set,
         'sequence_length' AS test_type,
         seq AS raw_sequence,
         det AS details
  FROM (
    SELECT seq, det FROM (
      SELECT raw_sequence AS seq,
             char_length(compressed_sequence) = len AS result,
             (char_length(compressed_sequence)::text || ' vs ' || len) as det
      FROM dna_sequence_test_short_ascii_ic
    ) AS b
    WHERE result = false
  ) AS a;

/* complement, reverse and reverse_complement functions */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence)
  SELECT 'dna_sequence_test_short_ascii_ic' AS test_set,
         'complement' AS test_type,
         seq AS raw_sequence
  FROM (
    SELECT seq FROM (
      SELECT raw_sequence AS seq, reverse_complement(reverse(complement(compressed_sequence)))::text = raw_sequence AS result
      FROM dna_sequence_test_short_ascii_ic
    ) AS b
    WHERE result = false
  ) AS a;

/* get_alphabet function */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence, details)
  SELECT 'dna_sequence_test_short_ascii_ic' AS test_set,
         'get_alphabet' AS test_type,
         raw_sequence,
         det AS details
  FROM (
    SELECT raw_sequence,
          (get_alphabet(raw_sequence)::text  || get_alphabet(compressed_sequence)::text) AS det,
          get_alphabet(raw_sequence)::text = get_alphabet(compressed_sequence)::text AS result
    FROM dna_sequence_test_short_ascii_ic
  ) AS a
  WHERE result = false;

/* strpos function */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence)
  SELECT 'dna_sequence_test_short_ascii_ic' AS test_set,
         'string matching' AS test_type,
         seq AS raw_sequence
  FROM (
    SELECT seq FROM (
      SELECT raw_sequence AS seq, (strpos(compressed_sequence, 'ACGTA') = strpos(raw_sequence, 'ACGTA')) AS result
      FROM dna_sequence_test_short_ascii_ic
    ) AS b
    WHERE result = false
  ) AS a;

DROP TABLE dna_sequence_test_short_ascii_ic;

/*
* Type modifier combination 6: SHORT, ASCII, CASE_SENSITIVE
*/
DROP TABLE IF EXISTS dna_sequence_test_short_ascii_cs;
CREATE TABLE dna_sequence_test_short_ascii_cs (
  id serial primary key,
  raw_sequence text,
  len int,
  compressed_sequence dna_sequence(SHORT, ascii, CASE_SENSITIVE)
);

/* 1000 DNA sequences with ascii code, 10 <= len < 2010 */
INSERT INTO dna_sequence_test_short_ascii_cs (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::dna_sequence(SHORT,ascii,CASE_SENSITIVE) FROM (
    SELECT generate_sequence('{A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z}'::alphabet, random_length) AS seq, random_length AS len FROM (
      SELECT (random() * 2000 + 10)::int AS random_length, generate_series(1, 1000)
    ) AS b
  ) AS a;

/* 100 DNA sequences with ascii code, 1000000 <= len < 2000000 */
INSERT INTO dna_sequence_test_short_ascii_cs (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::dna_sequence(SHORT,ascii,CASE_SENSITIVE) FROM (
    SELECT generate_sequence('{A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z}'::alphabet, random_length) AS seq, random_length AS len FROM (
      SELECT (random() * 1000000 + 1000000)::int AS random_length, generate_series(1, 100)
    ) AS b
  ) AS a;

/* full sequence decompression */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence)
  SELECT 'dna_sequence_test_short_ascii_cs' AS test_set,
         'full_sequence_decode' AS test_type,
         seq AS raw_sequence
  FROM (
    SELECT seq FROM (
      SELECT raw_sequence AS seq, compressed_sequence::text = raw_sequence AS result
      FROM dna_sequence_test_short_ascii_cs
    ) AS b
    WHERE result = false
  ) AS a;

/* substr function */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence, details)
  SELECT 'dna_sequence_test_short_ascii_cs' AS test_set,
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
        FROM dna_sequence_test_short_ascii_cs
      ) AS c
    ) AS b
    WHERE result = false
  ) AS a;

/* char_length function */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence, details)
  SELECT 'dna_sequence_test_short_ascii_cs' AS test_set,
         'sequence_length' AS test_type,
         seq AS raw_sequence,
         det AS details
  FROM (
    SELECT seq, det FROM (
      SELECT raw_sequence AS seq,
             char_length(compressed_sequence) = len AS result,
             (char_length(compressed_sequence)::text || ' vs ' || len) as det
      FROM dna_sequence_test_short_ascii_cs
    ) AS b
    WHERE result = false
  ) AS a;

/* complement, reverse and reverse_complement functions */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence)
  SELECT 'dna_sequence_test_short_ascii_cs' AS test_set,
         'complement' AS test_type,
         seq AS raw_sequence
  FROM (
    SELECT seq FROM (
      SELECT raw_sequence AS seq, reverse_complement(reverse(complement(compressed_sequence)))::text = raw_sequence AS result
      FROM dna_sequence_test_short_ascii_cs
    ) AS b
    WHERE result = false
  ) AS a;

/* get_alphabet function */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence, details)
  SELECT 'dna_sequence_test_short_ascii_cs' AS test_set,
         'get_alphabet' AS test_type,
         raw_sequence,
         det AS details
  FROM (
    SELECT raw_sequence,
          (get_alphabet(raw_sequence)::text  || get_alphabet(compressed_sequence)::text) AS det,
          get_alphabet(raw_sequence)::text = get_alphabet(compressed_sequence)::text AS result
    FROM dna_sequence_test_short_ascii_cs
  ) AS a
  WHERE result = false;

/* strpos function */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence)
  SELECT 'dna_sequence_test_short_ascii_cs' AS test_set,
         'string matching' AS test_type,
         seq AS raw_sequence
  FROM (
    SELECT seq FROM (
      SELECT raw_sequence AS seq, (strpos(compressed_sequence, 'ACGTA') = strpos(raw_sequence, 'ACGTA')) AS result
      FROM dna_sequence_test_short_ascii_cs
    ) AS b
    WHERE result = false
  ) AS a;

DROP TABLE dna_sequence_test_short_ascii_cs;
/*
* Type modifier combination 7: DEFAULT
*/
DROP TABLE IF EXISTS dna_sequence_test_default;
CREATE TABLE dna_sequence_test_default (
  id serial primary key,
  raw_sequence text,
  len int,
  compressed_sequence dna_sequence
);

/* 1000 DNA sequences with iupac code and GC-content = 0.5, 10 <= len < 2010 */
INSERT INTO dna_sequence_test_default (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::dna_sequence FROM (
    SELECT generate_sequence(dna_iupac(), random_length) AS seq, random_length AS len FROM (
      SELECT (random() * 2000 + 10)::int AS random_length, generate_series(1, 1000)
    ) AS b
  ) AS a;

/* 100 DNA sequences with iupac code and GC-content = 0.5, 1000000 <= len < 2000000 */
INSERT INTO dna_sequence_test_default (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::dna_sequence FROM (
    SELECT generate_sequence(dna_flc(), random_length) AS seq, random_length AS len FROM (
      SELECT (random() * 1000000 + 1000000)::int AS random_length, generate_series(1, 100)
    ) AS b
  ) AS a;

/* 1000 DNA sequences with four-letter code and GC-content = 0.2, 10 <= len < 2010 */
INSERT INTO dna_sequence_test_default (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::dna_sequence FROM (
    SELECT generate_sequence(dna_flc('{0.1,0.1,0.4,0.4}'), random_length) AS seq, random_length AS len FROM (
      SELECT (random() * 2000 + 10)::int AS random_length, generate_series(1, 1000)
    ) AS b
  ) AS a;

/* 100 DNA sequences with four-letter code and GC-content = 0.2, 1000000 <= len < 2000000 */
INSERT INTO dna_sequence_test_default (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::dna_sequence FROM (
    SELECT generate_sequence(dna_flc('{0.1,0.1,0.4,0.4}'), random_length) AS seq, random_length AS len FROM (
      SELECT (random() * 1000000 + 1000000)::int AS random_length, generate_series(1, 100)
    ) AS b
  ) AS a;

/* 100 long DNA sequences that should result in swapping */
INSERT INTO dna_sequence_test_default (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::dna_sequence FROM (
    SELECT generate_sequence(dna_iupac('{0.25,0.25,0.25,0.249, 0.001, 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}'), random_length) AS seq, random_length AS len FROM (
      SELECT (random() * 1000000 + 1000000)::int AS random_length, generate_series(1, 100)
    ) AS b
  ) AS a;
  
/* 1000 short DNA sequences that should result in swapping */
INSERT INTO dna_sequence_test_default (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::dna_sequence FROM (
    SELECT generate_sequence(dna_iupac('{0.25,0.25,0.25,0.249, 0.001, 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}'), random_length) AS seq, random_length AS len FROM (
      SELECT (random() * 40000 + 20000)::int AS random_length, generate_series(1, 1000)
    ) AS b
  ) AS a;

/* full sequence decompression */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence)
  SELECT 'dna_sequence_test_default' AS test_set,
         'full_sequence_decode' AS test_type,
         seq AS raw_sequence
  FROM (
    SELECT seq FROM (
      SELECT raw_sequence AS seq, compressed_sequence::text = raw_sequence AS result
      FROM dna_sequence_test_default
    ) AS b
    WHERE result = false
  ) AS a;

/* substr function */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence, details)
  SELECT 'dna_sequence_test_default' AS test_set,
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
        FROM dna_sequence_test_default
      ) AS c
    ) AS b
    WHERE result = false
  ) AS a;

/* char_length function */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence, details)
  SELECT 'dna_sequence_test_default' AS test_set,
         'sequence_length' AS test_type,
         seq AS raw_sequence,
         det AS details
  FROM (
    SELECT seq, det FROM (
      SELECT raw_sequence AS seq,
             char_length(compressed_sequence) = len AS result,
             (char_length(compressed_sequence)::text || ' vs ' || len) as det
      FROM dna_sequence_test_default
    ) AS b
    WHERE result = false
  ) AS a;

/* complement, reverse and reverse_complement functions */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence)
  SELECT 'dna_sequence_test_default' AS test_set,
         'complement' AS test_type,
         seq AS raw_sequence
  FROM (
    SELECT seq FROM (
      SELECT raw_sequence AS seq, reverse_complement(reverse(complement(compressed_sequence)))::text = raw_sequence AS result
      FROM dna_sequence_test_default
    ) AS b
    WHERE result = false
  ) AS a;

/* get_alphabet function */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence, details)
  SELECT 'dna_sequence_test_default' AS test_set,
         'get_alphabet' AS test_type,
         raw_sequence,
         det AS details
  FROM (
    SELECT raw_sequence,
          (get_alphabet(raw_sequence)::text  || get_alphabet(compressed_sequence)::text) AS det,
          get_alphabet(raw_sequence)::text = get_alphabet(compressed_sequence)::text AS result
    FROM dna_sequence_test_default
  ) AS a
  WHERE result = false;

/* transcribe, reverse_transcribe functions */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence)
  SELECT 'dna_sequence_test_default' AS test_set,
         'transcription' AS test_type,
         seq AS raw_sequence
  FROM (
    SELECT seq FROM (
      SELECT raw_sequence AS seq, reverse_transcribe(transcribe(compressed_sequence))::text = raw_sequence AS result
      FROM dna_sequence_test_default
    ) AS b
    WHERE result = false
  ) AS a;

/* strpos function */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence)
         'string matching' AS test_type,
         seq AS raw_sequence
  FROM (
    SELECT seq FROM (
      SELECT raw_sequence AS seq, (strpos(compressed_sequence, 'ACGTA') = strpos(raw_sequence, 'ACGTA')) AS result
      FROM dna_sequence_test_default
    ) AS b
    WHERE result = false
  ) AS a;

DROP TABLE dna_sequence_test_default;

/*
* Type modifier combination 8: DEFAULT, CASE_SENSITIVE
*/
DROP TABLE IF EXISTS dna_sequence_test_default_cs;
CREATE TABLE dna_sequence_test_default_cs (
  id serial primary key,
  raw_sequence text,
  len int,
  compressed_sequence dna_sequence(CASE_SENSITIVE)
);

/* 100 long DNA sequences that should result in swapping */
INSERT INTO dna_sequence_test_default_cs (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::dna_sequence(CASE_SENSITIVE) FROM (
    SELECT generate_sequence('{{A,C,G,T,a,c,g,t},{0.25,0.25,0.25,0.246,0.001,0.001,0.001,0.001}}'::alphabet, random_length) AS seq, random_length AS len FROM (
      SELECT (random() * 1000000 + 1000000)::int AS random_length, generate_series(1, 100)
    ) AS b
  ) AS a;
  
/* 1000 short DNA sequences that should result in swapping */
INSERT INTO dna_sequence_test_default_cs (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::dna_sequence(CASE_SENSITIVE) FROM (
    SELECT generate_sequence('{{A,C,G,T,a,c,g,t},{0.25,0.25,0.25,0.246,0.001,0.001,0.001,0.001}}'::alphabet, random_length) AS seq, random_length AS len FROM (
      SELECT (random() * 40000 + 20000)::int AS random_length, generate_series(1, 1000)
    ) AS b
  ) AS a;

/* full sequence decompression */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence)
  SELECT 'dna_sequence_test_default_cs' AS test_set,
         'full_sequence_decode' AS test_type,
         seq AS raw_sequence
  FROM (
    SELECT seq FROM (
      SELECT raw_sequence AS seq, compressed_sequence::text = raw_sequence AS result
      FROM dna_sequence_test_default_cs
    ) AS b
    WHERE result = false
  ) AS a;

/* substr function */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence, details)
  SELECT 'dna_sequence_test_default_cs' AS test_set,
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
        FROM dna_sequence_test_default_cs
      ) AS c
    ) AS b
    WHERE result = false
  ) AS a;

/* char_length function */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence, details)
  SELECT 'dna_sequence_test_default_cs' AS test_set,
         'sequence_length' AS test_type,
         seq AS raw_sequence,
         det AS details
  FROM (
    SELECT seq, det FROM (
      SELECT raw_sequence AS seq,
             char_length(compressed_sequence) = len AS result,
             (char_length(compressed_sequence)::text || ' vs ' || len) as det
      FROM dna_sequence_test_default_cs
    ) AS b
    WHERE result = false
  ) AS a;

/* complement, reverse and reverse_complement functions */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence)
  SELECT 'dna_sequence_test_default_cs' AS test_set,
         'complement' AS test_type,
         seq AS raw_sequence
  FROM (
    SELECT seq FROM (
      SELECT raw_sequence AS seq, reverse_complement(reverse(complement(compressed_sequence)))::text = raw_sequence AS result
      FROM dna_sequence_test_default_cs
    ) AS b
    WHERE result = false
  ) AS a;

/* get_alphabet function */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence, details)
  SELECT 'dna_sequence_test_default_cs' AS test_set,
         'get_alphabet' AS test_type,
         raw_sequence,
         det AS details
  FROM (
    SELECT raw_sequence,
          (get_alphabet(raw_sequence)::text  || get_alphabet(compressed_sequence)::text) AS det,
          get_alphabet(raw_sequence)::text = get_alphabet(compressed_sequence)::text AS result
    FROM dna_sequence_test_default_cs
  ) AS a
  WHERE result = false;

/* transcribe, reverse_transcribe functions */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence)
  SELECT 'dna_sequence_test_default_cs' AS test_set,
         'transcription' AS test_type,
         seq AS raw_sequence
  FROM (
    SELECT seq FROM (
      SELECT raw_sequence AS seq, reverse_transcribe(transcribe(compressed_sequence))::text = raw_sequence AS result
      FROM dna_sequence_test_default_cs
    ) AS b
    WHERE result = false
  ) AS a;

/* strpos function */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence)
  SELECT 'dna_sequence_test_default_cs' AS test_set,
         'string matching' AS test_type,
         seq AS raw_sequence
  FROM (
    SELECT seq FROM (
      SELECT raw_sequence AS seq, (strpos(compressed_sequence, 'ACGTA') = strpos(raw_sequence, 'ACGTA')) AS result
      FROM dna_sequence_test_default_cs
    ) AS b
    WHERE result = false
  ) AS a;

DROP TABLE dna_sequence_test_default_cs;

/*
* Type modifier combination 9: REFERENCE
*/
DROP TABLE IF EXISTS dna_sequence_test_reference;
CREATE TABLE dna_sequence_test_reference (
  id serial primary key,
  raw_sequence text,
  len int,
  compressed_sequence dna_sequence(REFERENCE)
);

/* 100 long DNA sequences that should result in run-length encoding */
INSERT INTO dna_sequence_test_reference (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::dna_sequence(REFERENCE) FROM (
    SELECT (generate_sequence(dna_flc(), random_length * 4) || repeat('NNNNNNNNNNACGT', random_length) || generate_sequence(dna_flc(), random_length * 4)) AS seq,
           random_length * 22 AS len
    FROM (
      SELECT (random() * 10000 + 10000)::int AS random_length, generate_series(1, 100)
    ) AS b
  ) AS a;

/* 10 long DNA sequences that should result in run-length encoding, with RLE at the end */
INSERT INTO dna_sequence_test_reference (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::dna_sequence(REFERENCE) FROM (
    SELECT (generate_sequence(dna_flc(), random_length * 4) || repeat('NNNNNNNNNNACGT', random_length) || generate_sequence(dna_flc(), random_length * 4) || 'NNNNNNNNNN') AS seq,
           random_length * 22 + 10 AS len
    FROM (
      SELECT (random() * 10000 + 10000)::int AS random_length, generate_series(1, 10)
    ) AS b
  ) AS a;

/* 100 short DNA sequences that should result in run-length encoding */
INSERT INTO dna_sequence_test_reference (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::dna_sequence(REFERENCE) FROM (
    SELECT (generate_sequence(dna_flc(), random_length) || repeat('N', random_length) || generate_sequence(dna_flc(), random_length)) AS seq,
           random_length * 3 AS len
    FROM (
      SELECT (random() * 1000 + 1000)::int AS random_length, generate_series(1, 100)
    ) AS b
  ) AS a;

/* 10 short DNA sequences that should result in run-length encoding, with RLE at the end */
INSERT INTO dna_sequence_test_reference (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::dna_sequence(REFERENCE) FROM (
    SELECT (generate_sequence(dna_flc(), random_length) || repeat('N', random_length) || generate_sequence(dna_flc(), random_length) || 'NNNNNNNNNN') AS seq,
           random_length * 3 + 10 AS len
    FROM (
      SELECT (random() * 1000 + 1000)::int AS random_length, generate_series(1, 10)
    ) AS b
  ) AS a;

/* 100 long DNA sequences that should result in run-length encoding and swapping */
INSERT INTO dna_sequence_test_reference (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::dna_sequence(REFERENCE) FROM (
    SELECT (generate_sequence('{{A,C,G,T,N},{0.1,0.4,0.4,0.099,0.001}}'::alphabet, random_length) || repeat('N', 250) || generate_sequence('{{A,C,G,T,N},{0.1,0.4,0.4,0.099,0.001}}'::alphabet, random_length)) AS seq,
           random_length * 2 + 250 AS len
    FROM (
      SELECT (random() * 100000 + 100000)::int AS random_length, generate_series(1, 100)
    ) AS b
  ) AS a;
  
/* 10 long DNA sequences that should result in run-length encoding and swapping, with RLE at the end */
INSERT INTO dna_sequence_test_reference (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::dna_sequence(REFERENCE) FROM (
    SELECT (generate_sequence(dna_flc(), random_length) || repeat('N', 250) || generate_sequence(dna_flc(), random_length) || 'NNNNNNNNNN') AS seq,
           random_length * 2 + 260 AS len
    FROM (
      SELECT (random() * 100000 + 100000)::int AS random_length, generate_series(1, 10)
    ) AS b
  ) AS a;

/* 100 short DNA sequences that should result in run-length encoding and swapping */
INSERT INTO dna_sequence_test_reference (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::dna_sequence(REFERENCE) FROM (
    SELECT (generate_sequence('{{A,C,G,T,N},{0.1,0.4,0.4,0.099,0.001}}'::alphabet, random_length) || repeat('N', 100) || generate_sequence('{{A,C,G,T,N},{0.1,0.4,0.4,0.099,0.001}}'::alphabet, random_length)) AS seq,
           random_length * 2 + 100 AS len
    FROM (
      SELECT (random() * 10000 + 10000)::int AS random_length, generate_series(1, 100)
    ) AS b
  ) AS a;

/* 10 short DNA sequences that should result in run-length encoding and swapping, with RLE at the end */
INSERT INTO dna_sequence_test_reference (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::dna_sequence(REFERENCE) FROM (
    SELECT (generate_sequence(dna_flc(), random_length) || repeat('N', 100) || generate_sequence(dna_flc(), random_length) || 'NNNNNNNNNN') AS seq,
           random_length * 2 + 110 AS len
    FROM (
      SELECT (random() * 10000 + 10000)::int AS random_length, generate_series(1, 10)
    ) AS b
  ) AS a;

/* 100 short DNA sequences that should result in run-length encoding and swapping, with RLE as master symbol */
INSERT INTO dna_sequence_test_reference (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::dna_sequence(REFERENCE) FROM (
    SELECT (generate_sequence(dna_flc(), random_length) || 
             'N' || repeat('A', random_length) ||
             'N' || repeat('C', random_length) ||
             'N' || repeat('G', random_length) ||
             'N' || repeat('T', random_length))
             AS seq,
           random_length * 5 + 4 AS len
    FROM (
      SELECT (random() * 5000 + 5000)::int AS random_length, generate_series(1, 100)
    ) AS b
  ) AS a;
  
/* 100 long DNA sequences that should result in run-length encoding and swapping, with RLE as master symbol */
INSERT INTO dna_sequence_test_reference (raw_sequence, len, compressed_sequence)
  SELECT seq, len, seq::dna_sequence(REFERENCE) FROM (
    SELECT (generate_sequence(dna_flc(), random_length) || 
             'N' || repeat('A', random_length) ||
             'N' || repeat('C', random_length) ||
             'N' || repeat('G', random_length) ||
             'N' || repeat('T', random_length))
             AS seq,
           random_length * 5 + 4 AS len
    FROM (
      SELECT (random() * 100000 + 100000)::int AS random_length, generate_series(1, 100)
    ) AS b
  ) AS a;

/* full sequence decompression */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence)
  SELECT 'dna_sequence_test_reference' AS test_set,
         'full_sequence_decode' AS test_type,
         seq AS raw_sequence
  FROM (
    SELECT seq FROM (
      SELECT raw_sequence AS seq, compressed_sequence::text = raw_sequence AS result
      FROM dna_sequence_test_reference
    ) AS b
    WHERE result = false
  ) AS a;

/* substr function */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence, details)
  SELECT 'dna_sequence_test_reference' AS test_set,
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
        FROM dna_sequence_test_reference
      ) AS c
    ) AS b
    WHERE result = false
  ) AS a;

/* char_length function */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence, details)
  SELECT 'dna_sequence_test_reference' AS test_set,
         'sequence_length' AS test_type,
         seq AS raw_sequence,
         det AS details
  FROM (
    SELECT seq, det FROM (
      SELECT raw_sequence AS seq,
             char_length(compressed_sequence) = len AS result,
             (char_length(compressed_sequence)::text || ' vs ' || len) as det
      FROM dna_sequence_test_reference
    ) AS b
    WHERE result = false
  ) AS a;

/* reverse, complement and reverse_complement functions */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence)
  SELECT 'dna_sequence_test_reference' AS test_set,
         'complement' AS test_type,
         seq AS raw_sequence
  FROM (
    SELECT seq FROM (
      SELECT raw_sequence AS seq, reverse_complement(reverse(complement(compressed_sequence)))::text = raw_sequence AS result
      FROM dna_sequence_test_reference
    ) AS b
    WHERE result = false
  ) AS a;

/* get_alphabet function */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence, details)
  SELECT 'dna_sequence_test_reference' AS test_set,
         'get_alphabet' AS test_type,
         raw_sequence,
         det AS details
  FROM (
    SELECT raw_sequence,
          (get_alphabet(raw_sequence)::text  || get_alphabet(compressed_sequence)::text) AS det,
          get_alphabet(raw_sequence)::text = get_alphabet(compressed_sequence)::text AS result
    FROM dna_sequence_test_reference
  ) AS a
  WHERE result = false;

/* transcribe, reverse_transcribe functions */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence)
  SELECT 'dna_sequence_test_reference' AS test_set,
         'transcription' AS test_type,
         seq AS raw_sequence
  FROM (
    SELECT seq FROM (
      SELECT raw_sequence AS seq, reverse_transcribe(transcribe(compressed_sequence))::text = raw_sequence AS result
      FROM dna_sequence_test_reference
    ) AS b
    WHERE result = false
  ) AS a;

/* strpos function */
INSERT INTO dna_sequence_errors (test_set, test_type, raw_sequence)
  SELECT 'dna_sequence_test_reference' AS test_set,
         'string matching' AS test_type,
         seq AS raw_sequence
  FROM (
    SELECT seq FROM (
      SELECT raw_sequence AS seq, (strpos(compressed_sequence, 'ACGTA') = strpos(raw_sequence, 'ACGTA')) AS result
      FROM dna_sequence_test_reference
    ) AS b
    WHERE result = false
  ) AS a;

DROP TABLE dna_sequence_test_reference;

SELECT test_set, test_type, count(*) FROM dna_sequence_errors GROUP BY test_set, test_type;

/*
*SET client_min_messages=DEBUG1;
*SELECT show_counts();
*/

DROP EXTENSION postbis;