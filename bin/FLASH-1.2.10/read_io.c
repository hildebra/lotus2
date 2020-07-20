/*
 * read_io.c: Code for input and output of reads, e.g. from/to FASTQ files.
 */

/*
 * Copyright (C) 2012 Tanja Magoc
 * Copyright (C) 2012, 2013, 2014 Eric Biggers
 *
 * This file is part of FLASH, a fast tool to merge overlapping paired-end
 * reads.
 *
 * FLASH is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 3 of the License, or (at your option)
 * any later version.
 *
 * FLASH is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with FLASH; if not, see http://www.gnu.org/licenses/.
 */

#include "iostream.h"
#include "read.h"
#include "read_io.h"
#include "util.h"

#include <assert.h>
#include <limits.h>

/******************************************
 * Read input                             *
 ******************************************/

static bool
load_fastq_read(struct input_stream *in, struct read *r,
		uint64_t *line_no_p)
{
	ssize_t ret;

	/* tag <newline>
	 * sequence <newline>
	 * + <newline>
	 * quality <newline>  */

	ret = input_stream_getline(in, &r->tag, &r->tag_bufsz);
	if (ret <= 0)
		return false;
	if (ret > INT_MAX)
		goto too_long;
	++*line_no_p;
	r->tag_len = ret;

	ret = input_stream_getline(in, &r->seq, &r->seq_bufsz);
	if (ret <= 0)
		goto unexpected_eof;
	if (ret > INT_MAX)
		goto too_long;
	++*line_no_p;
	r->seq_len = ret;

	ret = input_stream_getline(in, &r->qual, &r->qual_bufsz);
	if (ret <= 0)
		goto unexpected_eof;
	if (r->qual[0] != '+')
		goto expected_plus;
	++*line_no_p;

	ret = input_stream_getline(in, &r->qual, &r->qual_bufsz);
	if (ret <= 0)
		goto unexpected_eof;
	if (ret > INT_MAX)
		goto too_long;
	++*line_no_p;
	r->qual_len = ret;

	return true;


unexpected_eof:
	fatal_error("Unexpected EOF reading \"%s\" (line %"PRIu64")",
		    input_stream_get_name(in), *line_no_p);

expected_plus:
	fatal_error("Expected '+' character in FASTQ separator in \"%s\" "
		    "(line %"PRIu64")", input_stream_get_name(in),
		    *line_no_p);

too_long:
	fatal_error("Line %"PRIu64" in \"%s\" is too long",
		    *line_no_p, input_stream_get_name(in));
}

static bool
load_tab_delimited_read(struct input_stream *in, struct read *r,
			uint64_t *line_no_p)
{
	ssize_t ret;
	const char *delims = "\t\n";

	/* tag <tab> sequence <tab> quality <newline>  */

	ret = input_stream_getdelims(in, &r->tag, &r->tag_bufsz, delims);
	if (ret <= 0)
		return false;
	if (ret > INT_MAX)
		goto too_long;
	if (r->tag[ret - 1] != '\t')
		goto expected_tab;
	r->tag_len = ret;

	ret = input_stream_getdelims(in, &r->seq, &r->seq_bufsz, delims);
	if (ret <= 0)
		goto unexpected_eof;
	if (ret > INT_MAX)
		goto too_long;
	if (r->seq[ret - 1] != '\t')
		goto expected_tab;
	r->seq_len = ret;

	ret = input_stream_getdelims(in, &r->qual, &r->qual_bufsz, delims);
	if (ret <= 0)
		goto unexpected_eof;
	if (ret > INT_MAX)
		goto too_long;
	if (r->qual[ret - 1] == '\t')
		goto expected_newline;
	r->qual_len = ret;
	++*line_no_p;

	return true;

too_long:
	fatal_error("Field in line %"PRIu64" of \"%s\" is too long",
		    *line_no_p, input_stream_get_name(in));

unexpected_eof:
	fatal_error("Unexpected EOF reading \"%s\" (line %"PRIu64")",
		    input_stream_get_name(in), *line_no_p);

expected_tab:
	fatal_error("Invalid data in \"%s\": "
		    "expected tab character (line %"PRIu64")",
		    input_stream_get_name(in), *line_no_p);

expected_newline:
	fatal_error("Invalid data in \"%s\": "
		    "expected newline character (line %"PRIu64")",
		    input_stream_get_name(in), *line_no_p);
}

static bool
load_tab_delimited_pair(struct input_stream *in,
			struct read *r1, struct read *r2, uint64_t *line_no_p)
{
	ssize_t ret;
	const char *delims = "\t\n";

	/* tag <tab> seq_1 <tab> qual_1 <tab> seq_2 <tab> qual_2 <newline>  */

	ret = input_stream_getdelims(in, &r1->tag, &r1->tag_bufsz, delims);
	if (ret <= 0)
		return false;
	if (ret > INT_MAX)
		goto too_long;
	if (r1->tag[ret - 1] != '\t')
		goto expected_tab;
	r1->tag_len = ret;

	ret = input_stream_getdelims(in, &r1->seq, &r1->seq_bufsz, delims);
	if (ret <= 0)
		goto unexpected_eof;
	if (ret > INT_MAX)
		goto too_long;
	if (r1->seq[ret - 1] != '\t')
		goto expected_tab;
	r1->seq_len = ret;

	ret = input_stream_getdelims(in, &r1->qual, &r1->qual_bufsz, delims);
	if (ret <= 0)
		goto unexpected_eof;
	if (ret > INT_MAX)
		goto too_long;
	r1->qual_len = ret;

	if (r1->qual[ret - 1] == '\n') {
		/* Actually just a single read; use a void second read.  */
		r2->tag_len = 0;
		r2->seq_len = 0;
		r2->qual_len = 0;
		++*line_no_p;
		return true;
	}

	/* Set tag of read 2 to be the same as the tag of read 1  */
	copy_tag(r2, r1);

	ret = input_stream_getdelims(in, &r2->seq, &r2->seq_bufsz, delims);
	if (ret <= 0)
		goto unexpected_eof;
	if (ret > INT_MAX)
		goto too_long;
	if (r2->seq[ret - 1] != '\t')
		goto expected_tab;
	r2->seq_len = ret;

	ret = input_stream_getdelims(in, &r2->qual, &r2->qual_bufsz, delims);
	if (ret <= 0)
		goto unexpected_eof;
	if (ret > INT_MAX)
		goto too_long;
	if (r2->qual[ret - 1] == '\t')
		goto expected_newline;
	r2->qual_len = ret;
	++*line_no_p;

	return true;

too_long:
	fatal_error("Field in line %"PRIu64" of \"%s\" is too long",
		    *line_no_p, input_stream_get_name(in));

unexpected_eof:
	fatal_error("Unexpected EOF reading \"%s\" (line %"PRIu64")",
		    input_stream_get_name(in), *line_no_p);

expected_tab:
	fatal_error("Invalid data in \"%s\": "
		    "expected tab character (line %"PRIu64")",
		    input_stream_get_name(in), *line_no_p);

expected_newline:
	fatal_error("Invalid data in \"%s\": "
		    "expected newline character (line %"PRIu64")",
		    input_stream_get_name(in), *line_no_p);
}

/*
 * Loads the next read from the stream @in.
 *
 * @iparams specifies the format being used; e.g. FASTQ with a certain phred
 * offset.
 *
 * In each resulting read, whitespace is stripped from the end of the sequence,
 * tag, and quality scores.  The sequence is translated into only the characters
 * A, C, G, T, and N, and the quality values are re-scaled to start at 0.
 *
 * Returns true on success, false on end-of-file.  Aborts on read error or if
 * the data is invalid.
 */
bool
load_read(struct input_stream *in, const struct read_format_params *iparams,
	  struct read *r, uint64_t *line_no_p)
{
	bool ret;

	switch (iparams->fmt) {
	case READ_FORMAT_FASTQ:
		ret = load_fastq_read(in, r, line_no_p);
		break;
	case READ_FORMAT_TAB_DELIMITED:
		ret = load_tab_delimited_read(in, r, line_no_p);
		break;
	default:
		assert(0);
		ret = false;
	}

	if (ret)
		clean_read(r, iparams->phred_offset, in, *line_no_p);
	return ret;
}

/*
 * Similar to load_read(), but loads a pair of reads from the file instead.
 * This is only relevant (and must only be called) for file formats that store
 * both reads of the pair in the same sequential file.
 *
 * As a special case, this function may only fill in @r1, and set @r2->seq_len
 * to 0, to indicate that the next record in the file was actually an unpaired
 * read, not a read pair.  This is possible in formats for which
 * read_format_supports_mixed_reads() returns true (e.g. tab-delimited).
 */
bool
load_read_pair(struct input_stream *in, const struct read_format_params *iparams,
	       struct read *r1, struct read *r2,
	       uint64_t *line_no_p)
{
	bool ret;

	switch (iparams->fmt) {
	case READ_FORMAT_FASTQ:
		ret = load_fastq_read(in, r1, line_no_p);
		if (ret && !load_fastq_read(in, r2, line_no_p))
			fatal_error("Interleaved FASTQ file \"%s\" has an "
				    "odd number of reads",
				    input_stream_get_name(in));
		break;
	case READ_FORMAT_TAB_DELIMITED:
		ret = load_tab_delimited_pair(in, r1, r2, line_no_p);
		break;
	default:
		assert(0);
		ret = false;
	}
	if (ret) {
		clean_read(r1, iparams->phred_offset, in, *line_no_p);
		clean_read(r2, iparams->phred_offset, in, *line_no_p);
	}
	return ret;
}

/******************************************
 * Read output                            *
 ******************************************/

static void
write_fastq_read(struct output_stream *out, const struct read *r)
{
	/* Add '@' to tag if missing  */
	if (r->tag_len == 0 || r->tag[0] != '@')
		output_stream_fputc(out, '@');

	output_stream_write(out, r->tag, r->tag_len);
	output_stream_fputc(out, '\n');
	output_stream_write(out, r->seq, r->seq_len);
	output_stream_fputc(out, '\n');
	output_stream_fputc(out, '+');
	output_stream_fputc(out, '\n');
	output_stream_write(out, r->qual, r->qual_len);
	output_stream_fputc(out, '\n');
}

static void
write_tab_delimited_read(struct output_stream *out, const struct read *r)
{
	const char *tag = r->tag;
	int tag_len = r->tag_len;

	/* Strip '@' from tag  */
	if (tag_len > 0 && tag[0] == '@')
		tag++, tag_len--;

	output_stream_write(out, tag, tag_len);
	output_stream_fputc(out, '\t');
	output_stream_write(out, r->seq, r->seq_len);
	output_stream_fputc(out, '\t');
	output_stream_write(out, r->qual, r->qual_len);
	output_stream_fputc(out, '\n');
}

static void
write_tab_delimited_pair(struct output_stream *out,
			 const struct read *r1, const struct read *r2)
{
	const char *tag = r1->tag;
	int tag_len = r1->tag_len;

	/* Strip '@' and /1 or /2 from tag  */
	if (tag_len > 0 && tag[0] == '@')
		tag++, tag_len--;
	if (tag_len >= 2 && tag[tag_len - 2] == '/' &&
	    (tag[tag_len - 1] == '1' || tag[tag_len - 1] == '2'))
		tag_len -= 2;

	output_stream_write(out, tag, tag_len);
	output_stream_fputc(out, '\t');

	output_stream_write(out, r1->seq, r1->seq_len);
	output_stream_fputc(out, '\t');
	output_stream_write(out, r1->qual, r1->qual_len);
	output_stream_fputc(out, '\t');
	output_stream_write(out, r2->seq, r2->seq_len);
	output_stream_fputc(out, '\t');
	output_stream_write(out, r2->qual, r2->qual_len);
	output_stream_fputc(out, '\n');
}

/* Writes a read to the specified output stream in the format specified by
 * @oparams.
 *
 * Modifies the qual string of @r!  */
void
write_read(struct output_stream *out, const struct read_format_params *oparams,
	   struct read *r)
{
	clean_read_for_write(r, oparams->phred_offset);

	switch (oparams->fmt) {
	case READ_FORMAT_FASTQ:
		write_fastq_read(out, r);
		break;
	case READ_FORMAT_TAB_DELIMITED:
		write_tab_delimited_read(out, r);
		break;
	default:
		assert(0);
	}
}

/* Writes a read pair to the specified output stream in the format specified by
 * @oparams.
 *
 * Modifies the qual string of @r1 and @r2!  */
void
write_read_pair(struct output_stream *out,
		const struct read_format_params *oparams,
		struct read *r1, struct read *r2)
{
	clean_read_for_write(r1, oparams->phred_offset);
	clean_read_for_write(r2, oparams->phred_offset);

	switch (oparams->fmt) {
	case READ_FORMAT_FASTQ:
		/* Interleaved FASTQ format  */
		write_fastq_read(out, r1);
		write_fastq_read(out, r2);
		break;
	case READ_FORMAT_TAB_DELIMITED:
		/* Tab-delimited format, with two reads in a pair on one line  */
		write_tab_delimited_pair(out, r1, r2);
		break;
	default:
		assert(0);
	}
}
