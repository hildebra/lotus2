/*
 * read_util.c: Utility functions for processing reads
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

#include "read.h"
#include "iostream.h"
#include "util.h"

#include <assert.h>
#include <ctype.h>
#include <string.h>

/*
 * A table mapping ASCII characters (actually, 8-bit bytes) to "canonical" form
 * for base processing:
 *
 * a, A => A
 * c, C => C
 * g, G => G
 * t, T => T
 * everything else => N
 */
static const char canonical_ascii_tab[256] = {
        'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
        'N', /* A */ 'A', 'N', /* C */ 'C', 'N', 'N', 'N', /* G */ 'G',
                'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N', /* T */ 'T', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
        'N', /* a */ 'A', 'N', /* c */ 'C', 'N', 'N', 'N', /* g */ 'G',
                'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N', /* t */ 'T', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
};

/* Turns lowercase a, c, g, t into uppercase;
 * uppercase A, C, G, T stay the same;
 * everything else turns into 'N'.  */
static inline char
canonical_ascii_char(char c)
{
	return canonical_ascii_tab[(unsigned char)c];
}

/* A table mapping an ASCII base character in canonical form to its complement.
 * Unknown bases (N's) stay unknown.  */
static const char complement_tab[] = {
	['A'] = 'T',
	['C'] = 'G',
	['G'] = 'C',
	['T'] = 'A',
	['N'] = 'N',
};

/* Complements a canonical ASCII base (A, C, G, T, N). */
static inline char
complement(char c)
{
	return complement_tab[(unsigned char)c];
}

static inline char
identity_mapping(char c)
{
	return c;
}

/* Reverse a sequence of @len chars, applying the mapping function @map_char to
 * each, including the middle char if @len is odd.  */
static inline void
reverse_with_mapping(char *p, size_t len, char (*map_char)(char))
{
	char tmp;
	char *pp = p + len;
	while (pp > p) {
		--pp;
		tmp = *p;
		*p = (*map_char)(*pp);
		*pp = (*map_char)(tmp);
		++p;
	}
}

/* Reverse-complement a read in place.  */
void
reverse_complement(struct read *r)
{
	reverse_with_mapping(r->seq, r->seq_len, complement);
	reverse_with_mapping(r->qual, r->seq_len, identity_mapping);
}

/* Remove all whitespace from the end of the line/string.  Return the length of
 * the trimmed string. */
static inline int
trim(char *s, int len)
{
	while (len > 0 && isspace((unsigned char)s[len - 1]))
		s[--len] = '\0';
	return len;
}

void
clean_read(struct read *r, int phred_offset, struct input_stream *in,
	   uint64_t line_no)
{
	int seq_len;
	char *seq;
	char *qual;

	r->seq_len = trim(r->seq, r->seq_len);
	r->tag_len = trim(r->tag, r->tag_len);
	r->qual_len = trim(r->qual, r->qual_len);

	seq_len = r->seq_len;

	if (r->qual_len != seq_len) {
		fatal_error("Qual string length (%d) not the same as sequence "
			    "length (%d) (file \"%s\", near line %"PRIu64")",
			    r->qual_len, seq_len,
			    input_stream_get_name(in), line_no);
	}

	seq = r->seq;
	for (int i = 0; i < seq_len; i++) {
		if (isspace((unsigned char)seq[i])) {
			fatal_error("Invalid sequence string: "
				    "contains whitespace "
				    "(file \"%s\", near line %"PRIu64")",
				    input_stream_get_name(in), line_no);
		}
		seq[i] = canonical_ascii_char(seq[i]);
	}

	qual = r->qual;

	if (phred_offset > 0) {
		for (int i = 0; i < seq_len; i++) {
			if (qual[i] < phred_offset) {
				fatal_error("Qual string contains character "
					    "under phred_offset = %d "
					    "(file \"%s\", near line %"PRIu64")",
					    phred_offset,
					    input_stream_get_name(in), line_no);
			}
			qual[i] -= phred_offset;
		}
	}
}

void
clean_read_for_write(struct read *r, int phred_offset)
{

	assert(r->seq_len == r->qual_len);

	if (phred_offset > 0) {
		char *qual = r->qual;
		int qual_len = r->qual_len;
		for (int i = 0; i < qual_len; i++)
			qual[i] += phred_offset;
	}
}

void
copy_tag(struct read *to, const struct read *from)
{
	if (to->tag_bufsz < from->tag_len + 1) {
		to->tag = xrealloc(to->tag, from->tag_len + 1);
		to->tag_bufsz = from->tag_len + 1;
	}
	to->tag_len = from->tag_len;
	memcpy(to->tag, from->tag, from->tag_len + 1);
}

/*
 * Given the FASTQ tags of two paired-end reads, find the FASTQ tag to give to
 * the combined read.
 *
 * This is done by stripping off the characters trailing the '/' (e.g. "/1" and
 * "/2"), unless there is a "barcode" beginning with the '#' character, which is
 * kept.
 */
void
get_combined_tag(const struct read *read_1,
		 const struct read *read_2,
		 struct read *combined_read)
{
	char *p;
	copy_tag(combined_read, read_1);
	for (p = &combined_read->tag[combined_read->tag_len - 1];
	     p >= combined_read->tag;
	     p--)
	{
		if (*p == '/') {
			/* Tags are different, and there's a forward slash in
			 * the first tag.  Remove everything after the forward
			 * slash, unless there's a barcode, which we keep. */
			if (*(p + 1) != '\0' && *(p + 2) == '#') {
				/* read ID has a barcode. */
				do {
					*p = *(p + 2);
				} while (*(++p + 2) != '\0');
			}
			*p = '\0';
			combined_read->tag_len = p - combined_read->tag;
			break;
		}
	}
}
