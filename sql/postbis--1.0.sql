/*-------------------------------------------------------------------------
*
* Copyright (c) 2013, Max Planck Institute for Marine Microbiology
*
* This software is released under the PostgreSQL License
*
* Author: Michael Schneider <mschneid@mpi-bremen.de>
*
* IDENTIFICATION
*   sql/postbis--1.0.sql
*
*-------------------------------------------------------------------------
*/
\echo Use "CREATE EXTENSION postbis" to load this file. \quit

/*
*	Type: dna_sequence
*/
CREATE TYPE dna_sequence;

CREATE FUNCTION dna_sequence_typmod_in(cstring[]) RETURNS int4 AS
  '$libdir/postbis','dna_sequence_typmod_in'
  LANGUAGE c IMMUTABLE STRICT;
	
CREATE FUNCTION dna_sequence_typmod_out(int4)
  RETURNS cstring AS 
  '$libdir/postbis','dna_sequence_typmod_out'
  LANGUAGE c IMMUTABLE STRICT;

CREATE FUNCTION dna_sequence_in(cstring, oid, int4)
  RETURNS dna_sequence
  AS '$libdir/postbis','dna_sequence_in'
  LANGUAGE c IMMUTABLE STRICT;

CREATE FUNCTION dna_sequence_out(dna_sequence)
  RETURNS cstring
  AS '$libdir/postbis', 'dna_sequence_out'
  LANGUAGE c IMMUTABLE STRICT;

CREATE TYPE dna_sequence (
  input = dna_sequence_in,
  output = dna_sequence_out,
  typmod_in = dna_sequence_typmod_in,
  typmod_out = dna_sequence_typmod_out,
  internallength = VARIABLE,
  storage = EXTERNAL
);

CREATE FUNCTION dna_sequence_cast(dna_sequence, int4)
  RETURNS dna_sequence AS
  '$libdir/postbis', 'dna_sequence_cast'
  LANGUAGE c IMMUTABLE STRICT;

CREATE CAST (dna_sequence AS dna_sequence)
  WITH FUNCTION dna_sequence_cast(dna_sequence, int4) AS ASSIGNMENT;

CREATE FUNCTION dna_sequence_in(text, int4)
  RETURNS dna_sequence AS
  '$libdir/postbis', 'dna_sequence_in_varlena'
  LANGUAGE c IMMUTABLE STRICT;

CREATE CAST (text AS dna_sequence)
  WITH FUNCTION dna_sequence_in(text, int4) AS ASSIGNMENT;

CREATE FUNCTION dna_sequence_in(varchar, int4)
  RETURNS dna_sequence AS
  '$libdir/postbis', 'dna_sequence_in_varlena'
  LANGUAGE c IMMUTABLE STRICT;

CREATE CAST (varchar AS dna_sequence)
  WITH FUNCTION dna_sequence_in(varchar, int4) AS ASSIGNMENT;

CREATE FUNCTION dna_sequence_in(char, int4)
  RETURNS dna_sequence AS
  '$libdir/postbis', 'dna_sequence_in_varlena'
  LANGUAGE c IMMUTABLE STRICT;

CREATE CAST (char AS dna_sequence)
  WITH FUNCTION dna_sequence_in(char, int4) AS ASSIGNMENT;

CREATE FUNCTION dna_sequence_out_text(dna_sequence)
  RETURNS text AS
  '$libdir/postbis', 'dna_sequence_out_varlena'
  LANGUAGE c IMMUTABLE STRICT;

CREATE CAST (dna_sequence AS text)
  WITH FUNCTION dna_sequence_out_text(dna_sequence) AS ASSIGNMENT;

CREATE FUNCTION dna_sequence_out_varchar(dna_sequence)
  RETURNS varchar AS
  '$libdir/postbis', 'dna_sequence_out_varlena'
  LANGUAGE c IMMUTABLE STRICT;

CREATE CAST (dna_sequence AS varchar)
  WITH FUNCTION dna_sequence_out_varchar(dna_sequence) AS ASSIGNMENT;

CREATE FUNCTION dna_sequence_out_char(dna_sequence)
  RETURNS char AS
  '$libdir/postbis', 'dna_sequence_out_varlena'
  LANGUAGE c IMMUTABLE STRICT;

CREATE CAST (dna_sequence as char)
  WITH FUNCTION dna_sequence_out_char(dna_sequence) AS ASSIGNMENT;

CREATE FUNCTION substr(dna_sequence, int4, int4)
  RETURNS text AS
  '$libdir/postbis', 'dna_sequence_substring'
  LANGUAGE c VOLATILE STRICT;

CREATE FUNCTION char_length(dna_sequence)
  RETURNS int4 AS
  '$libdir/postbis', 'dna_sequence_char_length'
  LANGUAGE c IMMUTABLE STRICT;

CREATE FUNCTION compression_ratio(dna_sequence)
  RETURNS float8 AS
  '$libdir/postbis', 'dna_sequence_compression_ratio'
  LANGUAGE c IMMUTABLE STRICT;

CREATE FUNCTION complement(dna_sequence)
  RETURNS dna_sequence AS
  '$libdir/postbis', 'dna_sequence_complement'
  LANGUAGE c IMMUTABLE STRICT;

CREATE FUNCTION reverse(dna_sequence)
  RETURNS dna_sequence AS
  '$libdir/postbis', 'dna_sequence_reverse'
  LANGUAGE c IMMUTABLE STRICT;

CREATE FUNCTION reverse_complement(dna_sequence)
  RETURNS dna_sequence AS
  '$libdir/postbis', 'dna_sequence_reverse_complement'
  LANGUAGE c IMMUTABLE STRICT;

CREATE FUNCTION equal_dna(dna_sequence, dna_sequence)
  RETURNS bool AS
  '$libdir/postbis', 'equal_dna'
  LANGUAGE c IMMUTABLE STRICT;

CREATE OPERATOR = (
  leftarg = dna_sequence,
  rightarg = dna_sequence,
  procedure = equal_dna,
  commutator = =,
  negator = !=,
  restrict = eqsel,
  join = eqjoinsel,
  hashes,
  merges
);

CREATE FUNCTION not_equal_dna(dna_sequence, dna_sequence)
  RETURNS bool AS $$
    SELECT NOT equal_dna($1,$2);
$$ LANGUAGE SQL IMMUTABLE STRICT;

CREATE OPERATOR != (
  leftarg = dna_sequence,
  rightarg = dna_sequence,
  procedure = not_equal_dna,
  commutator = !=,
  negator = =,
  restrict = neqsel,
  join = neqjoinsel
);

CREATE FUNCTION compare_dna_lt(dna_sequence, dna_sequence)
  RETURNS bool AS
    '$libdir/postbis', 'compare_dna_lt'
  LANGUAGE c IMMUTABLE STRICT;

CREATE OPERATOR < (
  leftarg = dna_sequence,
  rightarg = dna_sequence,
  procedure = compare_dna_lt,
  commutator = >,
  negator = >=,
  restrict = scalarltsel,
  join = scalarltjoinsel
);

CREATE FUNCTION compare_dna_le(dna_sequence, dna_sequence)
  RETURNS bool AS
    '$libdir/postbis', 'compare_dna_le'
  LANGUAGE c IMMUTABLE STRICT;

CREATE OPERATOR <= (
  leftarg = dna_sequence,
  rightarg = dna_sequence,
  procedure = compare_dna_le,
  commutator = >=,
  negator = >,
  restrict = scalarltsel,
  join = scalarltjoinsel
);

CREATE FUNCTION compare_dna_gt(dna_sequence, dna_sequence)
  RETURNS bool AS
    '$libdir/postbis', 'compare_dna_gt'
  LANGUAGE c IMMUTABLE STRICT;

CREATE OPERATOR > (
  leftarg = dna_sequence,
  rightarg = dna_sequence,
  procedure = compare_dna_gt,
  commutator = <,
  negator = <=,
  restrict = scalargtsel,
  join = scalargtjoinsel
);

CREATE FUNCTION compare_dna_ge(dna_sequence, dna_sequence)
  RETURNS bool AS
    '$libdir/postbis', 'compare_dna_ge'
  LANGUAGE c IMMUTABLE STRICT;

CREATE OPERATOR >= (
  leftarg = dna_sequence,
  rightarg = dna_sequence,
  procedure = compare_dna_ge,
  commutator = <=,
  negator = <,
  restrict = scalargtsel,
  join = scalargtjoinsel
);

CREATE FUNCTION compare_dna(dna_sequence, dna_sequence)
  RETURNS integer AS
  '$libdir/postbis', 'compare_dna'
  LANGUAGE c IMMUTABLE STRICT;

CREATE OPERATOR CLASS dna_sequence_btree_ops
  DEFAULT FOR TYPE dna_sequence USING btree AS
    OPERATOR 1 < (dna_sequence, dna_sequence),
    OPERATOR 2 <= (dna_sequence, dna_sequence),
    OPERATOR 3 = (dna_sequence, dna_sequence),
    OPERATOR 4 >= (dna_sequence, dna_sequence),
    OPERATOR 5 > (dna_sequence, dna_sequence),
    FUNCTION 1 compare_dna(dna_sequence, dna_sequence);

CREATE FUNCTION hash_dna(dna_sequence)
  RETURNS integer AS
  '$libdir/postbis', 'hash_dna'
  LANGUAGE c IMMUTABLE STRICT;

CREATE OPERATOR CLASS dna_sequence_hash_ops
  DEFAULT FOR TYPE dna_sequence USING hash AS
    OPERATOR 1 = (dna_sequence, dna_sequence),
    FUNCTION 1 hash_dna(dna_sequence);

CREATE FUNCTION concat_dna(dna_sequence, dna_sequence)
  RETURNS dna_sequence AS $$
    SELECT ($1::text || $2::text)::dna_sequence;
  $$ LANGUAGE sql IMMUTABLE STRICT;

CREATE OPERATOR || (
  leftarg = dna_sequence,
  rightarg = dna_sequence,
  procedure = concat_dna,
  commutator = ||
);

CREATE FUNCTION strpos(dna_sequence, text)
  RETURNS int4 AS
  '$libdir/postbis', 'strpos_dna'
  LANGUAGE c IMMUTABLE STRICT;

CREATE FUNCTION strpos(dna_sequence, dna_sequence)
  RETURNS int4 AS $$
    SELECT strpos($1, $2::text);
  $$ LANGUAGE sql IMMUTABLE STRICT;

CREATE FUNCTION octet_length(dna_sequence)
  RETURNS int4 AS
  '$libdir/postbis', 'octet_length_dna'
  LANGUAGE c IMMUTABLE STRICT;

/*
*	Type: rna_sequence
*/
CREATE TYPE rna_sequence;

CREATE FUNCTION rna_sequence_typmod_in(cstring[])
  RETURNS int4 AS
  '$libdir/postbis', 'rna_sequence_typmod_in'
  LANGUAGE c IMMUTABLE STRICT;
	
CREATE FUNCTION rna_sequence_typmod_out(int4)
  RETURNS cstring AS
  '$libdir/postbis', 'rna_sequence_typmod_out'
  LANGUAGE c IMMUTABLE STRICT;

CREATE FUNCTION rna_sequence_in(cstring, oid, int4)
  RETURNS rna_sequence AS
  '$libdir/postbis', 'rna_sequence_in'
  LANGUAGE c IMMUTABLE STRICT;

CREATE FUNCTION rna_sequence_out(rna_sequence)
  RETURNS cstring AS
  '$libdir/postbis', 'rna_sequence_out'
  LANGUAGE c IMMUTABLE STRICT;

CREATE TYPE rna_sequence (
  input = rna_sequence_in,
  output = rna_sequence_out,
  typmod_in = rna_sequence_typmod_in,
  typmod_out = rna_sequence_typmod_out,
  internallength = VARIABLE,
  storage = EXTERNAL
);

CREATE FUNCTION rna_sequence_cast(rna_sequence,int4)
  RETURNS rna_sequence AS
  '$libdir/postbis', 'rna_sequence_cast'
  LANGUAGE c IMMUTABLE STRICT;
  
CREATE CAST (rna_sequence AS rna_sequence)
  WITH FUNCTION rna_sequence_cast(rna_sequence, int4) AS ASSIGNMENT;

CREATE FUNCTION rna_sequence_in(text, int4)
  RETURNS rna_sequence AS
  '$libdir/postbis', 'rna_sequence_in_varlena'
  LANGUAGE c IMMUTABLE STRICT;

CREATE CAST (text AS rna_sequence)
  WITH FUNCTION rna_sequence_in(text, int4) AS ASSIGNMENT;

CREATE FUNCTION rna_sequence_in(varchar, int4)
  RETURNS rna_sequence AS
  '$libdir/postbis', 'rna_sequence_in_varlena'
  LANGUAGE c IMMUTABLE STRICT;

CREATE CAST (varchar AS rna_sequence)
  WITH FUNCTION rna_sequence_in(varchar, int4) AS ASSIGNMENT;

CREATE FUNCTION rna_sequence_in(char, int4)
  RETURNS rna_sequence AS
  '$libdir/postbis', 'rna_sequence_in_varlena'
  LANGUAGE c IMMUTABLE STRICT;

CREATE CAST (char AS rna_sequence)
  WITH FUNCTION rna_sequence_in(char,int4) AS ASSIGNMENT;

CREATE FUNCTION rna_sequence_out_text(rna_sequence)
  RETURNS text AS
  '$libdir/postbis', 'rna_sequence_out_varlena'
  LANGUAGE c IMMUTABLE STRICT;

CREATE CAST (rna_sequence AS text)
  WITH FUNCTION rna_sequence_out_text(rna_sequence) AS ASSIGNMENT;

CREATE FUNCTION rna_sequence_out_varchar(rna_sequence)
  RETURNS varchar AS
  '$libdir/postbis', 'rna_sequence_out_varlena'
  LANGUAGE c IMMUTABLE STRICT;

CREATE CAST (rna_sequence AS varchar)
  WITH FUNCTION rna_sequence_out_varchar(rna_sequence) AS ASSIGNMENT;

CREATE FUNCTION rna_sequence_out_char(rna_sequence)
  RETURNS char AS
  '$libdir/postbis', 'rna_sequence_out_varlena'
  LANGUAGE c IMMUTABLE STRICT;

CREATE CAST (rna_sequence as char)
  WITH FUNCTION rna_sequence_out_char(rna_sequence) AS ASSIGNMENT;

CREATE FUNCTION substr(rna_sequence, int4, int4)
  RETURNS text AS
  '$libdir/postbis', 'rna_sequence_substring'
  LANGUAGE c IMMUTABLE STRICT;

CREATE FUNCTION char_length(rna_sequence)
  RETURNS int4 AS
  '$libdir/postbis', 'rna_sequence_char_length'
  LANGUAGE c IMMUTABLE STRICT;

CREATE FUNCTION compression_ratio(rna_sequence)
  RETURNS float8 AS
  '$libdir/postbis', 'rna_sequence_compression_ratio'
  LANGUAGE c IMMUTABLE STRICT;

CREATE FUNCTION complement(rna_sequence)
  RETURNS rna_sequence AS
  '$libdir/postbis', 'rna_sequence_complement'
  LANGUAGE c IMMUTABLE STRICT;

CREATE FUNCTION reverse(rna_sequence)
  RETURNS rna_sequence AS
  '$libdir/postbis', 'rna_sequence_reverse'
  LANGUAGE c IMMUTABLE STRICT;

CREATE FUNCTION reverse_complement(rna_sequence)
  RETURNS rna_sequence AS
  '$libdir/postbis', 'rna_sequence_reverse_complement'
  LANGUAGE c IMMUTABLE STRICT;

CREATE FUNCTION transcribe(dna_sequence)
  RETURNS rna_sequence AS
  '$libdir/postbis', 'transcribe_dna'
  LANGUAGE c IMMUTABLE STRICT;

CREATE FUNCTION reverse_transcribe(rna_sequence)
  RETURNS dna_sequence AS
  '$libdir/postbis', 'reverse_transcribe_rna'
  LANGUAGE c IMMUTABLE STRICT;

CREATE FUNCTION equal_rna(rna_sequence, rna_sequence)
  RETURNS bool AS
  '$libdir/postbis', 'equal_rna'
  LANGUAGE c IMMUTABLE STRICT;

CREATE OPERATOR = (
  leftarg = rna_sequence,
  rightarg = rna_sequence,
  procedure = equal_rna,
  commutator = =,
  negator = !=,
  restrict = eqsel,
  join = eqjoinsel,
  hashes,
  merges
);

CREATE FUNCTION not_equal_rna(rna_sequence, rna_sequence)
  RETURNS bool AS $$
    SELECT NOT equal_rna($1,$2);
$$ LANGUAGE SQL IMMUTABLE STRICT;

CREATE OPERATOR != (
  leftarg = rna_sequence,
  rightarg = rna_sequence,
  procedure = not_equal_rna,
  commutator = !=,
  negator = =,
  restrict = neqsel,
  join = neqjoinsel
);

CREATE FUNCTION compare_rna_lt(rna_sequence, rna_sequence)
  RETURNS bool AS
    '$libdir/postbis', 'compare_rna_lt'
  LANGUAGE c IMMUTABLE STRICT;

CREATE OPERATOR < (
  leftarg = rna_sequence,
  rightarg = rna_sequence,
  procedure = compare_rna_lt,
  commutator = >,
  negator = >=,
  restrict = scalarltsel,
  join = scalarltjoinsel
);

CREATE FUNCTION compare_rna_le(rna_sequence, rna_sequence)
  RETURNS bool AS
    '$libdir/postbis', 'compare_rna_le'
  LANGUAGE c IMMUTABLE STRICT;

CREATE OPERATOR <= (
  leftarg = rna_sequence,
  rightarg = rna_sequence,
  procedure = compare_rna_le,
  commutator = >=,
  negator = >,
  restrict = scalarltsel,
  join = scalarltjoinsel
);

CREATE FUNCTION compare_rna_gt(rna_sequence, rna_sequence)
  RETURNS bool AS
    '$libdir/postbis', 'compare_rna_gt'
  LANGUAGE c IMMUTABLE STRICT;

CREATE OPERATOR > (
  leftarg = rna_sequence,
  rightarg = rna_sequence,
  procedure = compare_rna_gt,
  commutator = <,
  negator = <=,
  restrict = scalargtsel,
  join = scalargtjoinsel
);

CREATE FUNCTION compare_rna_ge(rna_sequence, rna_sequence)
  RETURNS bool AS
    '$libdir/postbis', 'compare_rna_ge'
  LANGUAGE c IMMUTABLE STRICT;

CREATE OPERATOR >= (
  leftarg = rna_sequence,
  rightarg = rna_sequence,
  procedure = compare_rna_ge,
  commutator = <=,
  negator = <,
  restrict = scalargtsel,
  join = scalargtjoinsel
);

CREATE FUNCTION compare_rna(rna_sequence, rna_sequence)
  RETURNS integer AS
  '$libdir/postbis', 'compare_rna'
  LANGUAGE c IMMUTABLE STRICT;

CREATE OPERATOR CLASS rna_sequence_btree_ops
  DEFAULT FOR TYPE rna_sequence USING btree AS
    OPERATOR 1 < (rna_sequence, rna_sequence),
    OPERATOR 2 <= (rna_sequence, rna_sequence),
    OPERATOR 3 = (rna_sequence, rna_sequence),
    OPERATOR 4 >= (rna_sequence, rna_sequence),
    OPERATOR 5 > (rna_sequence, rna_sequence),
    FUNCTION 1 compare_rna(rna_sequence, rna_sequence);

CREATE FUNCTION hash_rna(rna_sequence)
  RETURNS integer AS
  '$libdir/postbis', 'hash_rna'
  LANGUAGE c IMMUTABLE STRICT;

CREATE OPERATOR CLASS rna_sequence_hash_ops
  DEFAULT FOR TYPE rna_sequence USING hash AS
    OPERATOR 1 = (rna_sequence, rna_sequence),
    FUNCTION 1 hash_rna(rna_sequence);

CREATE FUNCTION concat_rna(rna_sequence, rna_sequence)
  RETURNS rna_sequence AS $$
    SELECT ($1::text || $2::text)::rna_sequence;
  $$ LANGUAGE sql IMMUTABLE STRICT;

CREATE OPERATOR || (
  leftarg = rna_sequence,
  rightarg = rna_sequence,
  procedure = concat_rna,
  commutator = ||
);

CREATE FUNCTION strpos(rna_sequence, text)
  RETURNS int4 AS
  '$libdir/postbis', 'strpos_rna'
  LANGUAGE c IMMUTABLE STRICT;

CREATE FUNCTION strpos(rna_sequence, rna_sequence)
  RETURNS int4 AS $$
    SELECT strpos($1, $2::text);
  $$ LANGUAGE sql IMMUTABLE STRICT;

CREATE FUNCTION octet_length(rna_sequence)
  RETURNS int4 AS
  '$libdir/postbis', 'octet_length_rna'
  LANGUAGE c IMMUTABLE STRICT;

/*
*	Type: aa_sequence
*/
CREATE TYPE aa_sequence;

CREATE FUNCTION aa_sequence_typmod_in(cstring[])
  RETURNS int4 AS 
  '$libdir/postbis','aa_sequence_typmod_in'
  LANGUAGE c IMMUTABLE STRICT;
	
CREATE FUNCTION aa_sequence_typmod_out(int4)
  RETURNS cstring AS 
  '$libdir/postbis','aa_sequence_typmod_out'
  LANGUAGE c IMMUTABLE STRICT;

CREATE FUNCTION aa_sequence_in(cstring, oid, int4)
  RETURNS aa_sequence AS
  '$libdir/postbis','aa_sequence_in'
  LANGUAGE c IMMUTABLE STRICT;

CREATE FUNCTION aa_sequence_out(aa_sequence)
  RETURNS cstring AS
  '$libdir/postbis', 'aa_sequence_out'
  LANGUAGE c IMMUTABLE STRICT;

CREATE TYPE aa_sequence (
  input = aa_sequence_in,
  output = aa_sequence_out,
  typmod_in = aa_sequence_typmod_in,
  typmod_out = aa_sequence_typmod_out,
  internallength = VARIABLE,
  storage = EXTERNAL
);

CREATE FUNCTION aa_sequence_cast(aa_sequence,int4)
  RETURNS aa_sequence AS
  '$libdir/postbis', 'aa_sequence_cast'
  LANGUAGE c IMMUTABLE STRICT;

CREATE CAST (aa_sequence AS aa_sequence)
  WITH FUNCTION aa_sequence_cast(aa_sequence, int4) AS ASSIGNMENT;

CREATE FUNCTION aa_sequence_in(text, int4)
  RETURNS aa_sequence AS
  '$libdir/postbis', 'aa_sequence_in_varlena'
  LANGUAGE c IMMUTABLE STRICT;

CREATE CAST (text AS aa_sequence)
  WITH FUNCTION aa_sequence_in(text,int4) AS ASSIGNMENT;

CREATE FUNCTION aa_sequence_in(varchar, int4)
  RETURNS aa_sequence AS
  '$libdir/postbis', 'aa_sequence_in_varlena'
  LANGUAGE c IMMUTABLE STRICT;

CREATE CAST (varchar AS aa_sequence)
  WITH FUNCTION aa_sequence_in(varchar, int4) AS ASSIGNMENT;

CREATE FUNCTION aa_sequence_in(char, int4)
  RETURNS aa_sequence AS
  '$libdir/postbis', 'aa_sequence_in_varlena'
  LANGUAGE c IMMUTABLE STRICT;

CREATE CAST (char AS aa_sequence)
  WITH FUNCTION aa_sequence_in(char, int4) AS ASSIGNMENT;

CREATE FUNCTION aa_sequence_out_text(aa_sequence)
  RETURNS text AS
  '$libdir/postbis', 'aa_sequence_out_varlena'
  LANGUAGE c IMMUTABLE STRICT;

CREATE CAST (aa_sequence AS text)
  WITH FUNCTION aa_sequence_out_text(aa_sequence) AS ASSIGNMENT;

CREATE FUNCTION aa_sequence_out_varchar(aa_sequence)
  RETURNS varchar AS
  '$libdir/postbis', 'aa_sequence_out_varlena'
  LANGUAGE c IMMUTABLE STRICT;

CREATE CAST (aa_sequence AS varchar)
  WITH FUNCTION aa_sequence_out_varchar(aa_sequence) AS ASSIGNMENT;

CREATE FUNCTION aa_sequence_out_char(aa_sequence)
  RETURNS char AS
  '$libdir/postbis', 'aa_sequence_out_varlena'
  LANGUAGE c IMMUTABLE STRICT;

CREATE CAST (aa_sequence AS char)
  WITH FUNCTION aa_sequence_out_char(aa_sequence) AS ASSIGNMENT;

CREATE FUNCTION substr(aa_sequence, int4, int4)
  RETURNS text AS
  '$libdir/postbis', 'aa_sequence_substring'
  LANGUAGE c IMMUTABLE STRICT;

CREATE FUNCTION char_length(aa_sequence)
  RETURNS int4 AS
  '$libdir/postbis', 'aa_sequence_char_length'
  LANGUAGE c IMMUTABLE STRICT;

CREATE FUNCTION compression_ratio(aa_sequence)
  RETURNS float8 AS
  '$libdir/postbis', 'aa_sequence_compression_ratio'
  LANGUAGE c IMMUTABLE STRICT;

CREATE FUNCTION reverse(aa_sequence)
  RETURNS aa_sequence AS
  '$libdir/postbis', 'aa_sequence_reverse'
  LANGUAGE c IMMUTABLE STRICT;

CREATE FUNCTION equal_aa(aa_sequence, aa_sequence)
  RETURNS bool AS
  '$libdir/postbis', 'equal_aa'
  LANGUAGE c IMMUTABLE STRICT;

CREATE OPERATOR = (
  leftarg = aa_sequence,
  rightarg = aa_sequence,
  procedure = equal_aa,
  commutator = =,
  negator = !=,
  restrict = eqsel,
  join = eqjoinsel,
  hashes,
  merges
);

CREATE FUNCTION not_equal_aa(aa_sequence, aa_sequence)
  RETURNS bool AS $$
    SELECT NOT equal_aa($1,$2);
$$ LANGUAGE SQL IMMUTABLE STRICT;

CREATE OPERATOR != (
  leftarg = aa_sequence,
  rightarg = aa_sequence,
  procedure = not_equal_aa,
  commutator = !=,
  negator = =,
  restrict = neqsel,
  join = neqjoinsel
);

CREATE FUNCTION compare_aa_lt(aa_sequence, aa_sequence)
  RETURNS bool AS
    '$libdir/postbis', 'compare_aa_lt'
  LANGUAGE c IMMUTABLE STRICT;

CREATE OPERATOR < (
  leftarg = aa_sequence,
  rightarg = aa_sequence,
  procedure = compare_aa_lt,
  commutator = >,
  negator = >=,
  restrict = scalarltsel,
  join = scalarltjoinsel
);

CREATE FUNCTION compare_aa_le(aa_sequence, aa_sequence)
  RETURNS bool AS
    '$libdir/postbis', 'compare_aa_le'
  LANGUAGE c IMMUTABLE STRICT;

CREATE OPERATOR <= (
  leftarg = aa_sequence,
  rightarg = aa_sequence,
  procedure = compare_aa_le,
  commutator = >=,
  negator = >,
  restrict = scalarltsel,
  join = scalarltjoinsel
);

CREATE FUNCTION compare_aa_gt(aa_sequence, aa_sequence)
  RETURNS bool AS
    '$libdir/postbis', 'compare_aa_gt'
  LANGUAGE c IMMUTABLE STRICT;

CREATE OPERATOR > (
  leftarg = aa_sequence,
  rightarg = aa_sequence,
  procedure = compare_aa_gt,
  commutator = <,
  negator = <=,
  restrict = scalargtsel,
  join = scalargtjoinsel
);

CREATE FUNCTION compare_aa_ge(aa_sequence, aa_sequence)
  RETURNS bool AS
    '$libdir/postbis', 'compare_aa_ge'
  LANGUAGE c IMMUTABLE STRICT;

CREATE OPERATOR >= (
  leftarg = aa_sequence,
  rightarg = aa_sequence,
  procedure = compare_aa_ge,
  commutator = <=,
  negator = <,
  restrict = scalargtsel,
  join = scalargtjoinsel
);

CREATE FUNCTION compare_aa(aa_sequence, aa_sequence)
  RETURNS integer AS
  '$libdir/postbis', 'compare_aa'
  LANGUAGE c IMMUTABLE STRICT;

CREATE OPERATOR CLASS aa_sequence_btree_ops
  DEFAULT FOR TYPE aa_sequence USING btree AS
    OPERATOR 1 < (aa_sequence, aa_sequence),
    OPERATOR 2 <= (aa_sequence, aa_sequence),
    OPERATOR 3 = (aa_sequence, aa_sequence),
    OPERATOR 4 >= (aa_sequence, aa_sequence),
    OPERATOR 5 > (aa_sequence, aa_sequence),
    FUNCTION 1 compare_aa(aa_sequence, aa_sequence);

CREATE FUNCTION hash_aa(aa_sequence)
  RETURNS integer AS
  '$libdir/postbis', 'hash_aa'
  LANGUAGE c IMMUTABLE STRICT;

CREATE OPERATOR CLASS aa_sequence_hash_ops
  DEFAULT FOR TYPE aa_sequence USING hash AS
    OPERATOR 1 = (aa_sequence, aa_sequence),
    FUNCTION 1 hash_aa(aa_sequence);

CREATE FUNCTION concat_aa(aa_sequence, aa_sequence)
  RETURNS aa_sequence AS $$
    SELECT ($1::text || $2::text)::aa_sequence;
  $$ LANGUAGE sql IMMUTABLE STRICT;

CREATE OPERATOR || (
  leftarg = aa_sequence,
  rightarg = aa_sequence,
  procedure = concat_aa,
  commutator = ||
);

CREATE FUNCTION strpos(aa_sequence, text)
  RETURNS int4 AS
  '$libdir/postbis', 'strpos_aa'
  LANGUAGE c IMMUTABLE STRICT;

CREATE FUNCTION strpos(aa_sequence, aa_sequence)
  RETURNS int4 AS $$
    SELECT strpos($1, $2::text);
  $$ LANGUAGE sql IMMUTABLE STRICT;

CREATE FUNCTION octet_length(aa_sequence)
  RETURNS int4 AS
  '$libdir/postbis', 'octet_length_aa'
  LANGUAGE c IMMUTABLE STRICT;

CREATE FUNCTION translate(rna_sequence, text)
  RETURNS aa_sequence AS
  '$libdir/postbis', 'translate_rna'
  LANGUAGE c IMMUTABLE STRICT;

CREATE FUNCTION standard_code()
  RETURNS text AS $$
    SELECT 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'::TEXT;
  $$ LANGUAGE sql IMMUTABLE STRICT;

CREATE FUNCTION vertebrate_mitochondrial_code()
  RETURNS text AS $$
    SELECT 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG'::TEXT;
  $$ LANGUAGE sql IMMUTABLE STRICT;

CREATE FUNCTION yeast_mitochondrial_code()
  RETURNS text AS $$
    SELECT 'FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG'::TEXT;
  $$ LANGUAGE sql IMMUTABLE STRICT;

CREATE FUNCTION mold_protozoan_coelenterate_mitochondrial_and_mycoplasma_spiroplasma_code()
  RETURNS text AS $$
    SELECT 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'::TEXT;
  $$ LANGUAGE sql IMMUTABLE STRICT;

CREATE FUNCTION invertebrate_mitochondrial_code()
  RETURNS text AS $$
    SELECT 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG'::TEXT;
  $$ LANGUAGE sql IMMUTABLE STRICT;

CREATE FUNCTION ciliate_dasycladacean_hexamita_nuclear_code()
  RETURNS text AS $$
    SELECT 'FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'::TEXT;
  $$ LANGUAGE sql IMMUTABLE STRICT;

CREATE FUNCTION echinodem_flatworm_mitochondrial_code()
  RETURNS text AS $$
    SELECT 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG'::TEXT;
  $$ LANGUAGE sql IMMUTABLE STRICT;

CREATE FUNCTION euplotid_nuclear_code()
  RETURNS text AS $$
    SELECT 'FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'::TEXT;
  $$ LANGUAGE sql IMMUTABLE STRICT;

CREATE FUNCTION bacterial_archaeal_plant_plastid_code()
  RETURNS text AS $$
    SELECT 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'::TEXT;
  $$ LANGUAGE sql IMMUTABLE STRICT;

CREATE FUNCTION alternative_yeast_nuclear_code()
  RETURNS text AS $$
    SELECT 'FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'::TEXT;
  $$ LANGUAGE sql IMMUTABLE STRICT;

CREATE FUNCTION ascidian_mitochondrial_code()
  RETURNS text AS $$
    SELECT 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG'::TEXT;
  $$ LANGUAGE sql IMMUTABLE STRICT;

CREATE FUNCTION alternative_flatworm_mitochondrial_code()
  RETURNS text AS $$
    SELECT 'FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG'::TEXT;
  $$ LANGUAGE sql IMMUTABLE STRICT;

CREATE FUNCTION blepharisma_nuclear_code()
  RETURNS text AS $$
    SELECT 'FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'::TEXT;
  $$ LANGUAGE sql IMMUTABLE STRICT;

CREATE FUNCTION chlorophycean_mitochondrial_code()
  RETURNS text AS $$
    SELECT 'FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'::TEXT;
  $$ LANGUAGE sql IMMUTABLE STRICT;

CREATE FUNCTION trematode_mitochondrial_code()
  RETURNS text AS $$
    SELECT 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG'::TEXT;
  $$ LANGUAGE sql IMMUTABLE STRICT;

CREATE FUNCTION scenedesmus_obliquus_mitochondrial_code()
  RETURNS text AS $$
    SELECT 'FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'::TEXT;
  $$ LANGUAGE sql IMMUTABLE STRICT;

CREATE FUNCTION thraustochytrium_mitochondrial_code()
  RETURNS text AS $$
    SELECT 'FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'::TEXT;
  $$ LANGUAGE sql IMMUTABLE STRICT;

CREATE FUNCTION pterobranchia_mitochondrial_code()
  RETURNS text AS $$
    SELECT 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG'::TEXT;
  $$ LANGUAGE sql IMMUTABLE STRICT;

CREATE FUNCTION candidate_division_sr1_gracilibacteria_code()
  RETURNS text AS $$
    SELECT 'FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'::TEXT;
  $$ LANGUAGE sql IMMUTABLE STRICT;

CREATE FUNCTION get_transl_table(int)
  RETURNS text AS $$
  DECLARE
    result text;
  BEGIN
    CASE $1
      WHEN 1 THEN
        result := standard_code();
      WHEN 2 THEN
        result := vertebrate_mitochondrial_code();
      WHEN 3 THEN
        result := yeast_mitochondrial_code();
      WHEN 4 THEN
        result := mold_protozoan_coelenterate_mitochondrial_and_mycoplasma_spiroplasma_code();
      WHEN 5 THEN
        result := invertebrate_mitochondrial_code();
      WHEN 6 THEN
        result := ciliate_dasycladacean_hexamita_nuclear_code();
      WHEN 9 THEN
        result := echinodem_flatworm_mitochondrial_code();
      WHEN 10 THEN
        result := euplotid_nuclear_code();
      WHEN 11 THEN
        result := bacterial_archaeal_plant_plastid_code();
      WHEN 12 THEN
        result := alternative_yeast_nuclear_code();
      WHEN 13 THEN
        result := ascidian_mitochondrial_code();
      WHEN 14 THEN
        result := alternative_flatworm_mitochondrial_code();
      WHEN 15 THEN
        result := blepharisma_nuclear_code();
      WHEN 16 THEN
        result := chlorophycean_mitochondrial_code();
      WHEN 21 THEN
        result := trematode_mitochondrial_code();
      WHEN 22 THEN
        result := scenedesmus_obliquus_mitochondrial_code();
      WHEN 23 THEN
        result := thraustochytrium_mitochondrial_code();
      WHEN 24 THEN
        result := pterobranchia_mitochondrial_code();
      WHEN 25 THEN
        result := candidate_division_sr1_gracilibacteria_code();
      ELSE
        RAISE EXCEPTION 'translation table % does not exists', $1;
      END CASE;
    RETURN result;
  END;
  $$ LANGUAGE plpgsql IMMUTABLE STRICT;

CREATE FUNCTION translate(rna_sequence)
  RETURNS aa_sequence AS $$
    SELECT translate($1, standard_code());
  $$ LANGUAGE sql IMMUTABLE STRICT;

CREATE FUNCTION six_frame(rna_sequence)
  RETURNS SETOF rna_sequence AS $$
    SELECT ($1::rna_sequence)
    UNION ALL
    SELECT (reverse_complement($1)::rna_sequence)
    UNION ALL
    SELECT (substr($1, 2, char_length($1) - 1)::rna_sequence )
    UNION ALL
    SELECT (reverse_complement(substr($1, 2, char_length($1) - 1)::rna_sequence ))
    UNION ALL
    SELECT (substr($1, 3, char_length($1) - 2)::rna_sequence )
    UNION ALL
    SELECT (reverse_complement(substr($1, 3, char_length($1) - 2)::rna_sequence ))
    ;
  $$ LANGUAGE sql IMMUTABLE STRICT;

/*
*	Type: aligned_dna_sequence
*/
CREATE TYPE aligned_dna_sequence;

CREATE FUNCTION aligned_dna_sequence_typmod_in(cstring[])
  RETURNS int4 AS 
  '$libdir/postbis','aligned_dna_sequence_typmod_in'
  LANGUAGE c IMMUTABLE STRICT;
	
CREATE FUNCTION aligned_dna_sequence_typmod_out(int4)
  RETURNS cstring AS 
  '$libdir/postbis','aligned_dna_sequence_typmod_out'
  LANGUAGE c IMMUTABLE STRICT;

CREATE FUNCTION aligned_dna_sequence_in(cstring, oid, int4)
  RETURNS aligned_dna_sequence
  AS '$libdir/postbis','aligned_dna_sequence_in'
  LANGUAGE c IMMUTABLE STRICT;

CREATE FUNCTION aligned_dna_sequence_out(aligned_dna_sequence)
  RETURNS cstring
  AS '$libdir/postbis', 'aligned_dna_sequence_out'
  LANGUAGE c IMMUTABLE STRICT;

CREATE TYPE aligned_dna_sequence (
  input = aligned_dna_sequence_in,
  output = aligned_dna_sequence_out,
  typmod_in = aligned_dna_sequence_typmod_in,
  typmod_out = aligned_dna_sequence_typmod_out,
  internallength = VARIABLE,
  storage = EXTERNAL
);

CREATE FUNCTION aligned_dna_sequence_cast(aligned_dna_sequence,int4)
  RETURNS aligned_dna_sequence AS
  '$libdir/postbis', 'aligned_dna_sequence_cast'
  LANGUAGE c IMMUTABLE STRICT;

CREATE CAST (aligned_dna_sequence AS aligned_dna_sequence)
  WITH FUNCTION aligned_dna_sequence_cast(aligned_dna_sequence, int4) AS ASSIGNMENT;

CREATE FUNCTION aligned_dna_sequence_in(text, int4)
  RETURNS aligned_dna_sequence AS
  '$libdir/postbis', 'aligned_dna_sequence_in_varlena'
  LANGUAGE c IMMUTABLE STRICT;

CREATE CAST (text AS aligned_dna_sequence)
  WITH FUNCTION aligned_dna_sequence_in(text, int4) AS ASSIGNMENT;

CREATE FUNCTION aligned_dna_sequence_in(varchar, int4)
  RETURNS aligned_dna_sequence AS
  '$libdir/postbis', 'aligned_dna_sequence_in_varlena'
  LANGUAGE c IMMUTABLE STRICT;

CREATE CAST (varchar AS aligned_dna_sequence)
  WITH FUNCTION aligned_dna_sequence_in(varchar, int4) AS ASSIGNMENT;

CREATE FUNCTION aligned_dna_sequence_in(char, int4)
  RETURNS aligned_dna_sequence AS
  '$libdir/postbis', 'aligned_dna_sequence_in_varlena'
  LANGUAGE c IMMUTABLE STRICT;

CREATE CAST (char AS aligned_dna_sequence)
  WITH FUNCTION aligned_dna_sequence_in(char, int4) AS ASSIGNMENT;

CREATE FUNCTION aligned_dna_sequence_out_text(aligned_dna_sequence)
  RETURNS text AS
  '$libdir/postbis', 'aligned_dna_sequence_out_varlena'
  LANGUAGE c IMMUTABLE STRICT;

CREATE CAST (aligned_dna_sequence AS text)
  WITH FUNCTION aligned_dna_sequence_out_text(aligned_dna_sequence) AS ASSIGNMENT;

CREATE FUNCTION aligned_dna_sequence_out_varchar(aligned_dna_sequence)
  RETURNS varchar AS
  '$libdir/postbis', 'aligned_dna_sequence_out_varlena'
  LANGUAGE c IMMUTABLE STRICT;

CREATE CAST (aligned_dna_sequence as varchar)
  WITH FUNCTION aligned_dna_sequence_out_varchar(aligned_dna_sequence) AS ASSIGNMENT;

CREATE FUNCTION aligned_dna_sequence_out_char(aligned_dna_sequence)
  RETURNS char AS
  '$libdir/postbis', 'aligned_dna_sequence_out_varlena'
  LANGUAGE c IMMUTABLE STRICT;

CREATE CAST (aligned_dna_sequence as char)
  WITH FUNCTION aligned_dna_sequence_out_char(aligned_dna_sequence) AS ASSIGNMENT;

CREATE FUNCTION substr(aligned_dna_sequence, int4, int4)
  RETURNS text AS
  '$libdir/postbis', 'aligned_dna_sequence_substring'
  LANGUAGE c IMMUTABLE STRICT;

CREATE FUNCTION char_length(aligned_dna_sequence)
  RETURNS int4 AS
  '$libdir/postbis', 'aligned_dna_sequence_char_length'
  LANGUAGE c IMMUTABLE STRICT;

CREATE FUNCTION compression_ratio(aligned_dna_sequence)
  RETURNS float8 AS
  '$libdir/postbis', 'aligned_dna_sequence_compression_ratio'
  LANGUAGE c IMMUTABLE STRICT;

CREATE TYPE pairwise_dna_alignment AS (
  sequence1 aligned_dna_sequence,
  sequence2 aligned_dna_sequence,
  score int4
);

CREATE FUNCTION complement(aligned_dna_sequence)
  RETURNS aligned_dna_sequence AS
  '$libdir/postbis', 'aligned_dna_sequence_complement'
  LANGUAGE c IMMUTABLE STRICT;

CREATE FUNCTION reverse(aligned_dna_sequence)
  RETURNS aligned_dna_sequence AS
  '$libdir/postbis', 'aligned_dna_sequence_reverse'
  LANGUAGE c IMMUTABLE STRICT;

CREATE FUNCTION reverse_complement(aligned_dna_sequence)
  RETURNS aligned_dna_sequence AS
  '$libdir/postbis', 'aligned_dna_sequence_reverse_complement'
  LANGUAGE c IMMUTABLE STRICT;

CREATE FUNCTION equal_aligned_dna(aligned_dna_sequence, aligned_dna_sequence)
  RETURNS bool AS
  '$libdir/postbis', 'equal_aligned_dna'
  LANGUAGE c IMMUTABLE STRICT;

CREATE OPERATOR = (
  leftarg = aligned_dna_sequence,
  rightarg = aligned_dna_sequence,
  procedure = equal_aligned_dna,
  commutator = =,
  negator = !=,
  restrict = eqsel,
  join = eqjoinsel,
  hashes,
  merges
);

CREATE FUNCTION not_equal_aligned_dna(aligned_dna_sequence, aligned_dna_sequence)
  RETURNS bool AS $$
    SELECT NOT equal_aligned_dna($1,$2);
$$ LANGUAGE SQL IMMUTABLE STRICT;

CREATE OPERATOR != (
  leftarg = aligned_dna_sequence,
  rightarg = aligned_dna_sequence,
  procedure = not_equal_aligned_dna,
  commutator = !=,
  negator = =,
  restrict = neqsel,
  join = neqjoinsel
);

CREATE FUNCTION compare_aligned_dna_lt(aligned_dna_sequence, aligned_dna_sequence)
  RETURNS bool AS
    '$libdir/postbis', 'compare_aligned_dna_lt'
  LANGUAGE c IMMUTABLE STRICT;

CREATE OPERATOR < (
  leftarg = aligned_dna_sequence,
  rightarg = aligned_dna_sequence,
  procedure = compare_aligned_dna_lt,
  commutator = >,
  negator = >=,
  restrict = scalarltsel,
  join = scalarltjoinsel
);

CREATE FUNCTION compare_aligned_dna_le(aligned_dna_sequence, aligned_dna_sequence)
  RETURNS bool AS
    '$libdir/postbis', 'compare_aligned_dna_le'
  LANGUAGE c IMMUTABLE STRICT;

CREATE OPERATOR <= (
  leftarg = aligned_dna_sequence,
  rightarg = aligned_dna_sequence,
  procedure = compare_aligned_dna_le,
  commutator = >=,
  negator = >,
  restrict = scalarltsel,
  join = scalarltjoinsel
);

CREATE FUNCTION compare_aligned_dna_gt(aligned_dna_sequence, aligned_dna_sequence)
  RETURNS bool AS
    '$libdir/postbis', 'compare_aligned_dna_gt'
  LANGUAGE c IMMUTABLE STRICT;

CREATE OPERATOR > (
  leftarg = aligned_dna_sequence,
  rightarg = aligned_dna_sequence,
  procedure = compare_aligned_dna_gt,
  commutator = <,
  negator = <=,
  restrict = scalargtsel,
  join = scalargtjoinsel
);

CREATE FUNCTION compare_aligned_dna_ge(aligned_dna_sequence, aligned_dna_sequence)
  RETURNS bool AS
    '$libdir/postbis', 'compare_aligned_dna_ge'
  LANGUAGE c IMMUTABLE STRICT;

CREATE OPERATOR >= (
  leftarg = aligned_dna_sequence,
  rightarg = aligned_dna_sequence,
  procedure = compare_aligned_dna_ge,
  commutator = <=,
  negator = <,
  restrict = scalargtsel,
  join = scalargtjoinsel
);

CREATE FUNCTION compare_aligned_dna(aligned_dna_sequence, aligned_dna_sequence)
  RETURNS integer AS
  '$libdir/postbis', 'compare_aligned_dna'
  LANGUAGE c IMMUTABLE STRICT;

CREATE OPERATOR CLASS aligned_dna_sequence_btree_ops
  DEFAULT FOR TYPE aligned_dna_sequence USING btree AS
    OPERATOR 1 < (aligned_dna_sequence, aligned_dna_sequence),
    OPERATOR 2 <= (aligned_dna_sequence, aligned_dna_sequence),
    OPERATOR 3 = (aligned_dna_sequence, aligned_dna_sequence),
    OPERATOR 4 >= (aligned_dna_sequence, aligned_dna_sequence),
    OPERATOR 5 > (aligned_dna_sequence, aligned_dna_sequence),
    FUNCTION 1 compare_aligned_dna(aligned_dna_sequence, aligned_dna_sequence);

CREATE FUNCTION hash_aligned_dna(aligned_dna_sequence)
  RETURNS integer AS
  '$libdir/postbis', 'hash_aligned_dna'
  LANGUAGE c IMMUTABLE STRICT;

CREATE OPERATOR CLASS aligned_dna_sequence_hash_ops
  DEFAULT FOR TYPE aligned_dna_sequence USING hash AS
    OPERATOR 1 = (aligned_dna_sequence, aligned_dna_sequence),
    FUNCTION 1 hash_aligned_dna(aligned_dna_sequence);

CREATE FUNCTION concat_aligned_dna(aligned_dna_sequence, aligned_dna_sequence)
  RETURNS aligned_dna_sequence AS $$
    SELECT ($1::text || $2::text)::aligned_dna_sequence;
  $$ LANGUAGE sql IMMUTABLE STRICT;

CREATE OPERATOR || (
  leftarg = aligned_dna_sequence,
  rightarg = aligned_dna_sequence,
  procedure = concat_aligned_dna,
  commutator = ||
);

CREATE FUNCTION strpos(aligned_dna_sequence, text)
  RETURNS int4 AS
  '$libdir/postbis', 'strpos_aligned_dna'
  LANGUAGE c IMMUTABLE STRICT;

CREATE FUNCTION strpos(aligned_dna_sequence, aligned_dna_sequence)
  RETURNS int4 AS $$
    SELECT strpos($1, $2::text);
  $$ LANGUAGE sql IMMUTABLE STRICT;

CREATE FUNCTION octet_length(aligned_dna_sequence)
  RETURNS int4 AS
  '$libdir/postbis', 'octet_length_aligned_dna'
  LANGUAGE c IMMUTABLE STRICT;

/*
*	Type: aligned_rna_sequence
*/
CREATE TYPE aligned_rna_sequence;

CREATE FUNCTION aligned_rna_sequence_typmod_in(cstring[])
  RETURNS int4 AS 
  '$libdir/postbis','aligned_rna_sequence_typmod_in'
  LANGUAGE c IMMUTABLE STRICT;

CREATE FUNCTION aligned_rna_sequence_typmod_out(int4) 
  RETURNS cstring AS 
  '$libdir/postbis','aligned_rna_sequence_typmod_out'
  LANGUAGE c IMMUTABLE STRICT;

CREATE FUNCTION aligned_rna_sequence_in(cstring, oid, int4)
  RETURNS aligned_rna_sequence AS
  '$libdir/postbis','aligned_rna_sequence_in'
  LANGUAGE c IMMUTABLE STRICT;

CREATE FUNCTION aligned_rna_sequence_out(aligned_rna_sequence)
  RETURNS cstring AS
  '$libdir/postbis', 'aligned_rna_sequence_out'
  LANGUAGE c IMMUTABLE STRICT;

CREATE TYPE aligned_rna_sequence (
  input = aligned_rna_sequence_in,
  output = aligned_rna_sequence_out,
  typmod_in = aligned_rna_sequence_typmod_in,
  typmod_out = aligned_rna_sequence_typmod_out,
  internallength = VARIABLE,
  storage = EXTERNAL
);

CREATE FUNCTION aligned_rna_sequence_cast(aligned_rna_sequence, int4)
  RETURNS aligned_rna_sequence AS
  '$libdir/postbis', 'aligned_rna_sequence_cast'
  LANGUAGE c IMMUTABLE STRICT;

CREATE CAST (aligned_rna_sequence AS aligned_rna_sequence)
  WITH FUNCTION aligned_rna_sequence_cast(aligned_rna_sequence, int4) AS ASSIGNMENT;

CREATE FUNCTION aligned_rna_sequence_in(text, int4)
  RETURNS aligned_rna_sequence AS
  '$libdir/postbis', 'aligned_rna_sequence_in_varlena'
  LANGUAGE c IMMUTABLE STRICT;

CREATE CAST (text AS aligned_rna_sequence)
  WITH FUNCTION aligned_rna_sequence_in(text, int4) AS ASSIGNMENT;

CREATE FUNCTION aligned_rna_sequence_in(varchar, int4)
  RETURNS aligned_rna_sequence AS
  '$libdir/postbis', 'aligned_rna_sequence_in_varlena'
  LANGUAGE c IMMUTABLE STRICT;

CREATE CAST (varchar AS aligned_rna_sequence)
  WITH FUNCTION aligned_rna_sequence_in(varchar,int4) AS ASSIGNMENT;

CREATE FUNCTION aligned_rna_sequence_in(char, int4)
  RETURNS aligned_rna_sequence AS
  '$libdir/postbis', 'aligned_rna_sequence_in_varlena'
  LANGUAGE c IMMUTABLE STRICT;

CREATE CAST (char AS aligned_rna_sequence)
  WITH FUNCTION aligned_rna_sequence_in(char, int4) AS ASSIGNMENT;

CREATE FUNCTION aligned_rna_sequence_out_text(aligned_rna_sequence)
  RETURNS text AS
  '$libdir/postbis', 'aligned_rna_sequence_out_varlena'
  LANGUAGE c IMMUTABLE STRICT;

CREATE CAST (aligned_rna_sequence AS text)
  WITH FUNCTION aligned_rna_sequence_out_text(aligned_rna_sequence) AS ASSIGNMENT;

CREATE FUNCTION aligned_rna_sequence_out_varchar(aligned_rna_sequence)
  RETURNS varchar AS
  '$libdir/postbis', 'aligned_rna_sequence_out_varlena'
  LANGUAGE c IMMUTABLE STRICT;

CREATE CAST (aligned_rna_sequence AS varchar)
  WITH FUNCTION aligned_rna_sequence_out_varchar(aligned_rna_sequence) AS ASSIGNMENT;

CREATE FUNCTION aligned_rna_sequence_out_char(aligned_rna_sequence)
  RETURNS char AS
  '$libdir/postbis', 'aligned_rna_sequence_out_varlena'
  LANGUAGE c IMMUTABLE STRICT;

CREATE CAST (aligned_rna_sequence AS char)
  WITH FUNCTION aligned_rna_sequence_out_char(aligned_rna_sequence) AS ASSIGNMENT;

CREATE FUNCTION substr(aligned_rna_sequence, int4, int4)
  RETURNS text AS
  '$libdir/postbis', 'aligned_rna_sequence_substring'
  LANGUAGE c VOLATILE STRICT;

CREATE FUNCTION char_length(aligned_rna_sequence)
  RETURNS int4 AS
  '$libdir/postbis', 'aligned_rna_sequence_char_length'
  LANGUAGE c IMMUTABLE STRICT;

CREATE FUNCTION compression_ratio(aligned_rna_sequence)
  RETURNS float8 AS
  '$libdir/postbis', 'aligned_rna_sequence_compression_ratio'
  LANGUAGE c IMMUTABLE STRICT;

CREATE TYPE pairwise_rna_alignment AS (
  sequence1 aligned_rna_sequence,
  sequence2 aligned_rna_sequence,
  score int4
);

CREATE FUNCTION complement(aligned_rna_sequence)
  RETURNS aligned_rna_sequence AS
  '$libdir/postbis', 'aligned_rna_sequence_complement'
  LANGUAGE c IMMUTABLE STRICT;

CREATE FUNCTION reverse(aligned_rna_sequence)
  RETURNS aligned_rna_sequence AS
  '$libdir/postbis', 'aligned_rna_sequence_reverse'
  LANGUAGE c IMMUTABLE STRICT;

CREATE FUNCTION reverse_complement(aligned_rna_sequence)
  RETURNS aligned_rna_sequence AS
  '$libdir/postbis', 'aligned_rna_sequence_reverse_complement'
  LANGUAGE c IMMUTABLE STRICT;

CREATE FUNCTION equal_aligned_rna(aligned_rna_sequence, aligned_rna_sequence)
  RETURNS bool AS
  '$libdir/postbis', 'equal_aligned_rna'
  LANGUAGE c IMMUTABLE STRICT;

CREATE OPERATOR = (
  leftarg = aligned_rna_sequence,
  rightarg = aligned_rna_sequence,
  procedure = equal_aligned_rna,
  commutator = =,
  negator = !=,
  restrict = eqsel,
  join = eqjoinsel,
  hashes,
  merges
);

CREATE FUNCTION not_equal_aligned_rna(aligned_rna_sequence, aligned_rna_sequence)
  RETURNS bool AS $$
    SELECT NOT equal_aligned_rna($1,$2);
$$ LANGUAGE SQL IMMUTABLE STRICT;

CREATE OPERATOR != (
  leftarg = aligned_rna_sequence,
  rightarg = aligned_rna_sequence,
  procedure = not_equal_aligned_rna,
  commutator = !=,
  negator = =,
  restrict = neqsel,
  join = neqjoinsel
);

CREATE FUNCTION compare_aligned_rna_lt(aligned_rna_sequence, aligned_rna_sequence)
  RETURNS bool AS
    '$libdir/postbis', 'compare_aligned_rna_lt'
  LANGUAGE c IMMUTABLE STRICT;

CREATE OPERATOR < (
  leftarg = aligned_rna_sequence,
  rightarg = aligned_rna_sequence,
  procedure = compare_aligned_rna_lt,
  commutator = >,
  negator = >=,
  restrict = scalarltsel,
  join = scalarltjoinsel
);

CREATE FUNCTION compare_aligned_rna_le(aligned_rna_sequence, aligned_rna_sequence)
  RETURNS bool AS
    '$libdir/postbis', 'compare_aligned_rna_le'
  LANGUAGE c IMMUTABLE STRICT;

CREATE OPERATOR <= (
  leftarg = aligned_rna_sequence,
  rightarg = aligned_rna_sequence,
  procedure = compare_aligned_rna_le,
  commutator = >=,
  negator = >,
  restrict = scalarltsel,
  join = scalarltjoinsel
);

CREATE FUNCTION compare_aligned_rna_gt(aligned_rna_sequence, aligned_rna_sequence)
  RETURNS bool AS
    '$libdir/postbis', 'compare_aligned_rna_gt'
  LANGUAGE c IMMUTABLE STRICT;

CREATE OPERATOR > (
  leftarg = aligned_rna_sequence,
  rightarg = aligned_rna_sequence,
  procedure = compare_aligned_rna_gt,
  commutator = <,
  negator = <=,
  restrict = scalargtsel,
  join = scalargtjoinsel
);

CREATE FUNCTION compare_aligned_rna_ge(aligned_rna_sequence, aligned_rna_sequence)
  RETURNS bool AS
    '$libdir/postbis', 'compare_aligned_rna_ge'
  LANGUAGE c IMMUTABLE STRICT;

CREATE OPERATOR >= (
  leftarg = aligned_rna_sequence,
  rightarg = aligned_rna_sequence,
  procedure = compare_aligned_rna_ge,
  commutator = <=,
  negator = <,
  restrict = scalargtsel,
  join = scalargtjoinsel
);

CREATE FUNCTION compare_aligned_rna(aligned_rna_sequence, aligned_rna_sequence)
  RETURNS integer AS
  '$libdir/postbis', 'compare_aligned_rna'
  LANGUAGE c IMMUTABLE STRICT;

CREATE OPERATOR CLASS aligned_rna_sequence_btree_ops
  DEFAULT FOR TYPE aligned_rna_sequence USING btree AS
    OPERATOR 1 < (aligned_rna_sequence, aligned_rna_sequence),
    OPERATOR 2 <= (aligned_rna_sequence, aligned_rna_sequence),
    OPERATOR 3 = (aligned_rna_sequence, aligned_rna_sequence),
    OPERATOR 4 >= (aligned_rna_sequence, aligned_rna_sequence),
    OPERATOR 5 > (aligned_rna_sequence, aligned_rna_sequence),
    FUNCTION 1 compare_aligned_rna(aligned_rna_sequence, aligned_rna_sequence);

CREATE FUNCTION hash_aligned_rna(aligned_rna_sequence)
  RETURNS integer AS
  '$libdir/postbis', 'hash_aligned_rna'
  LANGUAGE c IMMUTABLE STRICT;

CREATE OPERATOR CLASS aligned_rna_sequence_hash_ops
  DEFAULT FOR TYPE aligned_rna_sequence USING hash AS
    OPERATOR 1 = (aligned_rna_sequence, aligned_rna_sequence),
    FUNCTION 1 hash_aligned_rna(aligned_rna_sequence);

CREATE FUNCTION concat_aligned_rna(aligned_rna_sequence, aligned_rna_sequence)
  RETURNS aligned_rna_sequence AS $$
    SELECT ($1::text || $2::text)::aligned_rna_sequence;
  $$ LANGUAGE sql IMMUTABLE STRICT;

CREATE OPERATOR || (
  leftarg = aligned_rna_sequence,
  rightarg = aligned_rna_sequence,
  procedure = concat_aligned_rna,
  commutator = ||
);

CREATE FUNCTION strpos(aligned_rna_sequence, text)
  RETURNS int4 AS
  '$libdir/postbis', 'strpos_aligned_rna'
  LANGUAGE c IMMUTABLE STRICT;

CREATE FUNCTION strpos(aligned_rna_sequence, aligned_rna_sequence)
  RETURNS int4 AS $$
    SELECT strpos($1, $2::text);
  $$ LANGUAGE sql IMMUTABLE STRICT;

CREATE FUNCTION octet_length(aligned_rna_sequence)
  RETURNS int4 AS
  '$libdir/postbis', 'octet_length_aligned_rna'
  LANGUAGE c IMMUTABLE STRICT;

/*
*	Type: aligned_aa_sequence
*/
CREATE TYPE aligned_aa_sequence;

CREATE FUNCTION aligned_aa_sequence_typmod_in(cstring[]) 
  RETURNS int4 AS 
  '$libdir/postbis','aligned_aa_sequence_typmod_in'
  LANGUAGE c IMMUTABLE STRICT;

CREATE FUNCTION aligned_aa_sequence_typmod_out(int4) 
  RETURNS cstring AS 
  '$libdir/postbis','aligned_aa_sequence_typmod_out'
  LANGUAGE c IMMUTABLE STRICT;

CREATE FUNCTION aligned_aa_sequence_in(cstring, oid, int4)
  RETURNS aligned_aa_sequence AS
  '$libdir/postbis','aligned_aa_sequence_in'
  LANGUAGE c IMMUTABLE STRICT;

CREATE FUNCTION aligned_aa_sequence_out(aligned_aa_sequence)
  RETURNS cstring
  AS '$libdir/postbis', 'aligned_aa_sequence_out'
  LANGUAGE c IMMUTABLE STRICT;

CREATE TYPE aligned_aa_sequence (
  input = aligned_aa_sequence_in,
  output = aligned_aa_sequence_out,
  typmod_in = aligned_aa_sequence_typmod_in,
  typmod_out = aligned_aa_sequence_typmod_out,
  internallength = VARIABLE,
  storage = EXTERNAL
);

CREATE FUNCTION aligned_aa_sequence_cast(aligned_aa_sequence,int4)
  RETURNS aligned_aa_sequence as
  '$libdir/postbis', 'aligned_aa_sequence_cast'
  LANGUAGE c IMMUTABLE STRICT;

CREATE CAST (aligned_aa_sequence as aligned_aa_sequence)
  WITH FUNCTION aligned_aa_sequence_cast(aligned_aa_sequence, int4) AS ASSIGNMENT ;

CREATE FUNCTION aligned_aa_sequence_in(text, int4)
  RETURNS aligned_aa_sequence AS
  '$libdir/postbis', 'aligned_aa_sequence_in_varlena'
  LANGUAGE c IMMUTABLE STRICT;

CREATE CAST (text AS aligned_aa_sequence)
  WITH FUNCTION aligned_aa_sequence_in(text,int4) AS ASSIGNMENT;

CREATE FUNCTION aligned_aa_sequence_in(varchar, int4)
  RETURNS aligned_aa_sequence AS
  '$libdir/postbis', 'aligned_aa_sequence_in_varlena'
  LANGUAGE c IMMUTABLE STRICT;

CREATE CAST (varchar AS aligned_aa_sequence)
  WITH FUNCTION aligned_aa_sequence_in(varchar,int4) AS ASSIGNMENT;

CREATE FUNCTION aligned_aa_sequence_in(char, int4)
  RETURNS aligned_aa_sequence AS
  '$libdir/postbis', 'aligned_aa_sequence_in_varlena'
  LANGUAGE c IMMUTABLE STRICT;

CREATE CAST (char AS aligned_aa_sequence)
  WITH FUNCTION aligned_aa_sequence_in(char, int4) AS ASSIGNMENT;

CREATE FUNCTION aligned_aa_sequence_out_text(aligned_aa_sequence)
  RETURNS text AS
  '$libdir/postbis', 'aligned_aa_sequence_out_varlena'
  LANGUAGE c IMMUTABLE STRICT;

CREATE CAST (aligned_aa_sequence as text)
  WITH FUNCTION aligned_aa_sequence_out_text(aligned_aa_sequence) AS ASSIGNMENT;

CREATE FUNCTION aligned_aa_sequence_out_varchar(aligned_aa_sequence)
  RETURNS varchar AS
  '$libdir/postbis', 'aligned_aa_sequence_out_varlena'
  LANGUAGE c IMMUTABLE STRICT;

CREATE CAST (aligned_aa_sequence as varchar)
  WITH FUNCTION aligned_aa_sequence_out_varchar(aligned_aa_sequence) AS ASSIGNMENT;

CREATE FUNCTION aligned_aa_sequence_out_char(aligned_aa_sequence)
  RETURNS char AS
  '$libdir/postbis', 'aligned_aa_sequence_out_varlena'
  LANGUAGE c IMMUTABLE STRICT;

CREATE CAST (aligned_aa_sequence as char)
  WITH FUNCTION aligned_aa_sequence_out_char(aligned_aa_sequence) AS ASSIGNMENT;

CREATE FUNCTION substr(aligned_aa_sequence, int4, int4)
  RETURNS text AS
  '$libdir/postbis', 'aligned_aa_sequence_substring'
  LANGUAGE c VOLATILE STRICT;

CREATE FUNCTION char_length(aligned_aa_sequence)
  RETURNS int4 AS
  '$libdir/postbis', 'aligned_aa_sequence_char_length'
  LANGUAGE c IMMUTABLE STRICT;

CREATE FUNCTION compression_ratio(aligned_aa_sequence)
  RETURNS float8 AS
  '$libdir/postbis', 'aligned_aa_sequence_compression_ratio'
  LANGUAGE c IMMUTABLE STRICT;

CREATE FUNCTION reverse(aligned_aa_sequence)
  RETURNS aligned_aa_sequence AS
  '$libdir/postbis', 'aligned_aa_sequence_reverse'
  LANGUAGE c IMMUTABLE STRICT;

CREATE TYPE pairwise_aa_alignment AS (
  sequence1 aligned_aa_sequence,
  sequence2 aligned_aa_sequence,
  score int4
);

CREATE FUNCTION equal_aligned_aa(aligned_aa_sequence, aligned_aa_sequence)
  RETURNS bool AS
  '$libdir/postbis', 'equal_aligned_aa'
  LANGUAGE c IMMUTABLE STRICT;

CREATE OPERATOR = (
  leftarg = aligned_aa_sequence,
  rightarg = aligned_aa_sequence,
  procedure = equal_aligned_aa,
  commutator = =,
  negator = !=,
  restrict = eqsel,
  join = eqjoinsel,
  hashes,
  merges
);

CREATE FUNCTION not_equal_aligned_aa(aligned_aa_sequence, aligned_aa_sequence)
  RETURNS bool AS $$
    SELECT NOT equal_aligned_aa($1,$2);
$$ LANGUAGE SQL IMMUTABLE STRICT;

CREATE OPERATOR != (
  leftarg = aligned_aa_sequence,
  rightarg = aligned_aa_sequence,
  procedure = not_equal_aligned_aa,
  commutator = !=,
  negator = =,
  restrict = neqsel,
  join = neqjoinsel
);

CREATE FUNCTION compare_aligned_aa_lt(aligned_aa_sequence, aligned_aa_sequence)
  RETURNS bool AS
    '$libdir/postbis', 'compare_aligned_aa_lt'
  LANGUAGE c IMMUTABLE STRICT;

CREATE OPERATOR < (
  leftarg = aligned_aa_sequence,
  rightarg = aligned_aa_sequence,
  procedure = compare_aligned_aa_lt,
  commutator = >,
  negator = >=,
  restrict = scalarltsel,
  join = scalarltjoinsel
);

CREATE FUNCTION compare_aligned_aa_le(aligned_aa_sequence, aligned_aa_sequence)
  RETURNS bool AS
    '$libdir/postbis', 'compare_aligned_aa_le'
  LANGUAGE c IMMUTABLE STRICT;

CREATE OPERATOR <= (
  leftarg = aligned_aa_sequence,
  rightarg = aligned_aa_sequence,
  procedure = compare_aligned_aa_le,
  commutator = >=,
  negator = >,
  restrict = scalarltsel,
  join = scalarltjoinsel
);

CREATE FUNCTION compare_aligned_aa_gt(aligned_aa_sequence, aligned_aa_sequence)
  RETURNS bool AS
    '$libdir/postbis', 'compare_aligned_aa_gt'
  LANGUAGE c IMMUTABLE STRICT;

CREATE OPERATOR > (
  leftarg = aligned_aa_sequence,
  rightarg = aligned_aa_sequence,
  procedure = compare_aligned_aa_gt,
  commutator = <,
  negator = <=,
  restrict = scalargtsel,
  join = scalargtjoinsel
);

CREATE FUNCTION compare_aligned_aa_ge(aligned_aa_sequence, aligned_aa_sequence)
  RETURNS bool AS
    '$libdir/postbis', 'compare_aligned_aa_ge'
  LANGUAGE c IMMUTABLE STRICT;

CREATE OPERATOR >= (
  leftarg = aligned_aa_sequence,
  rightarg = aligned_aa_sequence,
  procedure = compare_aligned_aa_ge,
  commutator = <=,
  negator = <,
  restrict = scalargtsel,
  join = scalargtjoinsel
);

CREATE FUNCTION compare_aligned_aa(aligned_aa_sequence, aligned_aa_sequence)
  RETURNS integer AS
  '$libdir/postbis', 'compare_aligned_aa'
  LANGUAGE c IMMUTABLE STRICT;

CREATE OPERATOR CLASS aligned_aa_sequence_btree_ops
  DEFAULT FOR TYPE aligned_aa_sequence USING btree AS
    OPERATOR 1 < (aligned_aa_sequence, aligned_aa_sequence),
    OPERATOR 2 <= (aligned_aa_sequence, aligned_aa_sequence),
    OPERATOR 3 = (aligned_aa_sequence, aligned_aa_sequence),
    OPERATOR 4 >= (aligned_aa_sequence, aligned_aa_sequence),
    OPERATOR 5 > (aligned_aa_sequence, aligned_aa_sequence),
    FUNCTION 1 compare_aligned_aa(aligned_aa_sequence, aligned_aa_sequence);

CREATE FUNCTION hash_aligned_aa(aligned_aa_sequence)
  RETURNS integer AS
  '$libdir/postbis', 'hash_aligned_aa'
  LANGUAGE c IMMUTABLE STRICT;

CREATE OPERATOR CLASS aligned_aa_sequence_hash_ops
  DEFAULT FOR TYPE aligned_aa_sequence USING hash AS
    OPERATOR 1 = (aligned_aa_sequence, aligned_aa_sequence),
    FUNCTION 1 hash_aligned_aa(aligned_aa_sequence);

CREATE FUNCTION concat_aligned_aa(aligned_aa_sequence, aligned_aa_sequence)
  RETURNS aligned_aa_sequence AS $$
    SELECT ($1::text || $2::text)::aligned_aa_sequence;
  $$ LANGUAGE sql IMMUTABLE STRICT;

CREATE OPERATOR || (
  leftarg = aligned_aa_sequence,
  rightarg = aligned_aa_sequence,
  procedure = concat_aligned_aa,
  commutator = ||
);

CREATE FUNCTION strpos(aligned_aa_sequence, text)
  RETURNS int4 AS
  '$libdir/postbis', 'strpos_aligned_aa'
  LANGUAGE c IMMUTABLE STRICT;

CREATE FUNCTION strpos(aligned_aa_sequence, aligned_aa_sequence)
  RETURNS int4 AS $$
    SELECT strpos($1, $2::text);
  $$ LANGUAGE sql IMMUTABLE STRICT;

CREATE FUNCTION octet_length(aligned_aa_sequence)
  RETURNS int4 AS
  '$libdir/postbis', 'octet_length_aligned_aa'
  LANGUAGE c IMMUTABLE STRICT;

/*
*	Type: alphabet
*/
CREATE TYPE alphabet;

CREATE FUNCTION alphabet_in(cstring)
  RETURNS alphabet AS
  '$libdir/postbis', 'alphabet_in'
  LANGUAGE c IMMUTABLE STRICT;
  
CREATE FUNCTION alphabet_out(alphabet)
  RETURNS cstring AS
  '$libdir/postbis', 'alphabet_out'
  LANGUAGE c IMMUTABLE STRICT;
  
CREATE TYPE alphabet (
   input = alphabet_in,
   output = alphabet_out,
   internallength = VARIABLE
);

CREATE FUNCTION alphabet_in_text(text)
  RETURNS alphabet AS
  '$libdir/postbis', 'alphabet_in_text'
  LANGUAGE c IMMUTABLE STRICT;

CREATE CAST (text as alphabet)
  WITH FUNCTION alphabet_in_text(text) AS ASSIGNMENT;
  
CREATE FUNCTION alphabet_out_text(alphabet)
  RETURNS text AS
  '$libdir/postbis', 'alphabet_out_text'
  LANGUAGE c IMMUTABLE STRICT;

CREATE CAST (alphabet as text)
  WITH FUNCTION alphabet_out_text(alphabet) AS ASSIGNMENT;

CREATE FUNCTION alphabet_in_textarray(text[]) RETURNS alphabet AS $$
  SELECT $1::text::alphabet;
  $$ LANGUAGE sql IMMUTABLE STRICT;
  
CREATE CAST (text[] AS alphabet)
  WITH FUNCTION alphabet_in_textarray(text[]) AS ASSIGNMENT;
  
CREATE FUNCTION alphabet_out_textarray(alphabet) RETURNS text[] AS $$
  SELECT $1::text::text[];
  $$ LANGUAGE sql IMMUTABLE STRICT;

CREATE CAST (alphabet AS text[])
  WITH FUNCTION alphabet_out_textarray(alphabet) AS ASSIGNMENT;

CREATE FUNCTION get_alphabet(text)
  RETURNS alphabet AS
  '$libdir/postbis', 'get_alphabet_text_sequence'
  LANGUAGE c IMMUTABLE STRICT;

CREATE FUNCTION get_alphabet(dna_sequence)
  RETURNS alphabet AS
  '$libdir/postbis', 'get_alphabet_dna_sequence'
  LANGUAGE c IMMUTABLE STRICT;

CREATE FUNCTION get_alphabet(rna_sequence)
  RETURNS alphabet AS
  '$libdir/postbis', 'get_alphabet_rna_sequence'
  LANGUAGE c IMMUTABLE STRICT;

CREATE FUNCTION get_alphabet(aa_sequence)
  RETURNS alphabet AS
  '$libdir/postbis', 'get_alphabet_aa_sequence'
  LANGUAGE c IMMUTABLE STRICT;

CREATE FUNCTION get_alphabet(aligned_dna_sequence)
  RETURNS alphabet AS
  '$libdir/postbis', 'get_alphabet_aligned_dna_sequence'
  LANGUAGE c IMMUTABLE STRICT;

CREATE FUNCTION get_alphabet(aligned_rna_sequence)
  RETURNS alphabet AS
  '$libdir/postbis', 'get_alphabet_aligned_rna_sequence'
  LANGUAGE c IMMUTABLE STRICT;

CREATE FUNCTION get_alphabet(aligned_aa_sequence)
  RETURNS alphabet AS
  '$libdir/postbis', 'get_alphabet_aligned_aa_sequence'
  LANGUAGE c IMMUTABLE STRICT;

/*
* DNA alphabet generator functions
*/
CREATE FUNCTION dna_flc() RETURNS alphabet AS $$
  SELECT '{A,C,G,T}'::alphabet
  $$ LANGUAGE SQL IMMUTABLE STRICT;

CREATE FUNCTION dna_flc(float4[]) RETURNS alphabet AS $$
  SELECT ('{{A,C,G,T},' || $1::text || '}')::alphabet
  $$ LANGUAGE SQL IMMUTABLE STRICT;

CREATE FUNCTION dna_flc_random() RETURNS alphabet AS $$
  SELECT dna_flc(random_probabilities) FROM (
    SELECT ('{' || p1 || ',' || p2 || ',' || p4 || ',' || p3 || '}')::float4[] AS random_probabilities FROM (
      SELECT (a*b) AS p1, (a*(1-b)) AS p2, ((1-a)*c) AS p3, (1 - a - (1-a)*c) AS p4 FROM (
        SELECT round((random() / 2 + 0.25)::numeric,2) AS a, round((random() / 5 + 0.4)::numeric,2) AS b, round((random() / 5 + 0.4)::numeric,2) AS c
      ) AS q1
    ) AS q2
  ) AS q3
  $$ LANGUAGE SQL VOLATILE STRICT;
  
CREATE FUNCTION dna_iupac() RETURNS alphabet AS $$
  SELECT '{A,C,G,T,N,R,M,K,Y,W,B,V,S,D,H}'::alphabet
  $$ LANGUAGE SQL IMMUTABLE STRICT;

CREATE FUNCTION dna_iupac(float4[]) RETURNS alphabet AS $$
  SELECT ('{{A,C,G,T,N,R,M,K,Y,W,B,V,S,D,H},' || $1::text || '}')::alphabet
  $$ LANGUAGE SQL IMMUTABLE STRICT;

CREATE FUNCTION dna_iupac_random() RETURNS alphabet AS $$
  SELECT dna_iupac(random_probabilities) FROM (
    SELECT ('{' || pA || ',' || pc || ',' || pG || ',' || pT || ',' || pN || ',' || pAmb || ',' || pAmb || ',' || pAmb || ',' || pAmb || ',' || pAmb || ',' || pAmb || ',' || pAmb || ',' || pAmb || ',' || pAmb || ',' || pAmb || '}')::float4[] AS random_probabilities FROM (
      SELECT (mainvsamb*(1-gccont)*atratio) AS pA,
             (mainvsamb*(1-gccont)*(1-atratio)) AS pT,
             (mainvsamb*gccont*gcratio) AS pG,
             (mainvsamb*gccont*(1-gcratio)) AS pC,
             ((1-mainvsamb)*nshare) AS pN,
             ((1-mainvsamb)*(1-nshare)/10) AS pAmb
      FROM (
        SELECT round((random() / 5 + 0.8)::numeric,2) AS mainvsamb,
               round((random() / 2 + 0.25)::numeric,2) AS gccont,
               round((random() / 5 + 0.45)::numeric,2) AS atratio,
               round((random() / 5 + 0.45)::numeric,2) AS gcratio,
               round((random() / 5 + 0.6)::numeric,2) AS nshare 
      ) AS q1
    ) AS q2
  ) AS q3
  $$ LANGUAGE SQL VOLATILE STRICT;

/*
* RNA alphabet generator functions
*/
CREATE FUNCTION rna_flc() RETURNS alphabet AS $$
  SELECT '{A,C,G,U}'::alphabet
  $$ LANGUAGE SQL IMMUTABLE STRICT;

CREATE FUNCTION rna_flc(float4[]) RETURNS alphabet AS $$
  SELECT ('{{A,C,G,U},' || $1::text || '}')::alphabet
  $$ LANGUAGE SQL IMMUTABLE STRICT;

CREATE FUNCTION rna_flc_random() RETURNS alphabet AS $$
  SELECT rna_flc(random_probabilities) FROM (
    SELECT ('{' || p1 || ',' || p2 || ',' || p4 || ',' || p3 || '}')::float4[] AS random_probabilities FROM (
      SELECT (a*b) AS p1, (a*(1-b)) AS p2, ((1-a)*c) AS p3, (1 - a - (1-a)*c) AS p4 FROM (
        SELECT round((random() / 2 + 0.25)::numeric,2) AS a, round((random() / 5 + 0.4)::numeric,2) AS b, round((random() / 5 + 0.4)::numeric,2) AS c
      ) AS q1
    ) AS q2
  ) AS q3
  $$ LANGUAGE SQL VOLATILE STRICT;
  
CREATE FUNCTION rna_iupac() RETURNS alphabet AS $$
  SELECT '{A,C,G,U,N,R,M,K,Y,W,B,V,S,D,H}'::alphabet
  $$ LANGUAGE SQL IMMUTABLE STRICT;

CREATE FUNCTION rna_iupac(float4[]) RETURNS alphabet AS $$
  SELECT ('{{A,C,G,U,N,R,M,K,Y,W,B,V,S,D,H},' || $1::text || '}')::alphabet
  $$ LANGUAGE SQL IMMUTABLE STRICT;

CREATE FUNCTION rna_iupac_random() RETURNS alphabet AS $$
  SELECT rna_iupac(random_probabilities) FROM (
    SELECT ('{' || pA || ',' || pc || ',' || pG || ',' || pU || ',' || pN || ',' || pAmb || ',' || pAmb || ',' || pAmb || ',' || pAmb || ',' || pAmb || ',' || pAmb || ',' || pAmb || ',' || pAmb || ',' || pAmb || ',' || pAmb || '}')::float4[] AS random_probabilities FROM (
      SELECT (mainvsamb*(1-gccont)*atratio) AS pA,
             (mainvsamb*(1-gccont)*(1-atratio)) AS pU,
             (mainvsamb*gccont*gcratio) AS pG,
             (mainvsamb*gccont*(1-gcratio)) AS pC,
             ((1-mainvsamb)*nshare) AS pN,
             ((1-mainvsamb)*(1-nshare)/10) AS pAmb
      FROM (
        SELECT round((random() / 5 + 0.8)::numeric,2) AS mainvsamb,
               round((random() / 2 + 0.25)::numeric,2) AS gccont,
               round((random() / 5 + 0.45)::numeric,2) AS atratio,
               round((random() / 5 + 0.45)::numeric,2) AS gcratio,
               round((random() / 5 + 0.6)::numeric,2) AS nshare 
      ) AS q1
    ) AS q2
  ) AS q3
  $$ LANGUAGE SQL VOLATILE STRICT;
  
/*
*	Amino acid alphabet generator functions
*/
CREATE FUNCTION aa_iupac() RETURNS alphabet AS $$
  SELECT '{A,V,M,I,L,P,W,F,Y,T,Q,G,S,C,N,K,R,H,E,D,B,Z,X}'::alphabet
  $$ LANGUAGE SQL IMMUTABLE STRICT;

CREATE FUNCTION aa_iupac(float4[]) RETURNS alphabet AS $$
  SELECT ('{{A,V,M,I,L,P,W,F,Y,T,Q,G,S,C,N,K,R,H,E,D,B,Z,X},' || $1::text || '}')::alphabet
  $$ LANGUAGE SQL IMMUTABLE STRICT;

CREATE FUNCTION aa_iupac_random() RETURNS alphabet AS $$
  SELECT aa_iupac(random_probabilities) FROM (
    SELECT ('{' ||
          p1 || ',' ||
          p2 || ',' ||
          p3 || ',' ||
          p4 || ',' ||
          p5 || ',' ||
          p6 || ',' ||
          p7 || ',' ||
          p8 || ',' ||
          p4 || ',' ||
          p3 || ',' ||
          p2 || ',' ||
          p1 || ',' ||
          p8 || ',' ||
          p7 || ',' ||
          p6 || ',' ||
          p5 || ',' ||
          p3 || ',' ||
          p2 || ',' ||
          p4 || ',' ||
          p5 || ',' ||
          p1 || ',' ||
          p7 || ',' ||
          p6 ||
          '}')::float4[] AS random_probabilities FROM (
      SELECT (a*b*c/3) as p1,
             (a*b*(1-c)/3) as p2,
             (a*(1-b)*c/3) as p3,
             (a*(1-b)*(1-c)/3) as p4,
             ((1-a)*b*c/3) as p5,
             ((1-a)*b*(1-c)/3) as p6,
             ((1-a)*(1-b)*c/3) as p7,
             ((1-a)*(1-b)*(1-c)/2) as p8
      FROM (
        SELECT round((random() / 2 + 0.25)::numeric,2) AS a,
               round((random() / 5 + 0.4)::numeric,2) AS b,
               round((random() )::numeric,2) AS c
      ) AS q1
    ) AS q2
  ) AS q3
  $$ LANGUAGE SQL VOLATILE STRICT;

/*
* Aligned DNA alphabet generator functions
*/
CREATE FUNCTION aligned_dna_flc() RETURNS alphabet AS $$
  SELECT '{A,C,G,T,-,.}'::alphabet
  $$ LANGUAGE SQL IMMUTABLE STRICT;

CREATE FUNCTION aligned_dna_flc(float4[]) RETURNS alphabet AS $$
  SELECT ('{{A,C,G,T,-,.},' || $1::text || '}')::alphabet
  $$ LANGUAGE SQL IMMUTABLE STRICT;

CREATE FUNCTION aligned_dna_flc_random() RETURNS alphabet AS $$
  SELECT aligned_dna_flc(random_probabilities) FROM (
    SELECT ('{' || pA || ',' || pC || ',' || pG || ',' || pT || ',' || pMinus || ',' || pDot || '}')::float4[] AS random_probabilities FROM (
      SELECT (mainvsal*gccontent*ratio) AS pA,
             (mainvsal*gccontent*(1-ratio)) AS pT,
             (mainvsal*(1-gccontent)*ratio) AS pG,
             (mainvsal*(1-gccontent)*(1-ratio)) AS pC,
             ((1-mainvsal)*ratio) AS pMinus,
             ((1-mainvsal)*(1-ratio)) AS pDot
      FROM (
        SELECT round((random())::numeric,2) AS mainvsal,
               round((random() / 2 + 0.25)::numeric,2) AS gccontent,
               round((random() / 5 + 0.4)::numeric,2) AS ratio
      ) AS q1
    ) AS q2
  ) AS q3
  $$ LANGUAGE SQL VOLATILE STRICT;
  
CREATE FUNCTION aligned_dna_iupac() RETURNS alphabet AS $$
  SELECT '{A,C,G,T,N,R,M,K,Y,W,B,V,S,D,H,-,.}'::alphabet
  $$ LANGUAGE SQL IMMUTABLE STRICT;

CREATE FUNCTION aligned_dna_iupac(float4[]) RETURNS alphabet AS $$
  SELECT ('{{A,C,G,T,N,R,M,K,Y,W,B,V,S,D,H,-,.},' || $1::text || '}')::alphabet
  $$ LANGUAGE SQL IMMUTABLE STRICT;

CREATE FUNCTION aligned_dna_iupac_random() RETURNS alphabet AS $$
  SELECT aligned_dna_iupac(random_probabilities) FROM (
    SELECT ('{' || pAT || ',' || pGC || ',' || pGC || ',' || pAT || ',' || pN || ',' || pAmb || ',' || pAmb || ',' || pAmb || ',' || pAmb || ',' || pAmb || ',' || pAmb || ',' || pAmb || ',' || pAmb || ',' || pAmb || ',' || pAmb || ',' || pMinus || ',' || pDot || '}')::float4[] AS random_probabilities FROM (
      SELECT (mainvsal*mainvsamb*(1-gccont)/2) AS pAT,
             (mainvsal*mainvsamb*gccont/2) AS pGC,
             (mainvsal*(1-mainvsamb)*nshare) AS pN,
             (mainvsal*(1-mainvsamb)*(1-nshare)/10) AS pAmb,
             ((1-mainvsal)*nshare) AS pMinus,
             ((1-mainvsal)*(1-nshare)) AS pDot
      FROM (
        SELECT round((random())::numeric,2) AS mainvsal,
               round((random() / 5 + 0.8)::numeric,2) AS mainvsamb,
               round((random() / 2 + 0.25)::numeric,2) AS gccont,
               round((random() / 5 + 0.6)::numeric,2) AS nshare 
      ) AS q1
    ) AS q2
  ) AS q3
  $$ LANGUAGE SQL VOLATILE STRICT;

CREATE FUNCTION entropy(alphabet) RETURNS float4 AS $$
  DECLARE
    alphabet_array TEXT[];
    dimensions int;
    elements int;
    entropy float4;
  BEGIN
    SELECT $1::TEXT::TEXT[] INTO alphabet_array;
    SELECT array_ndims(alphabet_array) INTO dimensions;
    IF dimensions = 1
    THEN
      SELECT array_length(alphabet_array,1) INTO elements;
      RAISE NOTICE '%d',elements;
      SELECT (1.0 / elements::numeric) * log(2, (1.0 / elements::numeric)) * elements::numeric * (-1) INTO entropy;
    ELSIF dimensions = 2
    THEN 
      SELECT array_length(alphabet_array,2) INTO elements;
      SELECT SUM(prob * log(2,prob)) * (-1)
      FROM (
        SELECT unnest(alphabet_array[2:2]::numeric[]) AS prob
      ) AS q1 INTO entropy;
    END IF;
    RETURN entropy;
  END;
  $$ LANGUAGE plpgsql IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION gc_content(alphabet)
  RETURNS real AS $$
    DECLARE
      alpha text[][];
      n int;
      gc real = 0.0;
    BEGIN
      SELECT $1::text[][] INTO alpha;
      SELECT array_length(alpha, 2) INTO n;
      FOR i IN 1..n LOOP
        IF alpha[1][i] = 'G' OR alpha[1][i] = 'C' THEN
          gc := gc + alpha[2][i]::real;
        END IF;
      END LOOP;
      RETURN gc;
    END;
  $$ LANGUAGE plpgsql IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION gc_content(dna_sequence)
  RETURNS real AS $$
    SELECT gc_content(get_alphabet($1));
  $$ LANGUAGE sql IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION gc_content(rna_sequence)
  RETURNS real AS $$
    SELECT gc_content(get_alphabet($1));
  $$ LANGUAGE sql IMMUTABLE STRICT;

/*
*	Test functions
*/

CREATE FUNCTION generate_sequence(alphabet, int)
  RETURNS text AS
  '$libdir/postbis', 'generate_sequence'
  LANGUAGE c VOLATILE STRICT;

