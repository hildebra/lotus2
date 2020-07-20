#ifndef _FLASH_READ_IO_H_
#define _FLASH_READ_IO_H_

#include <stdbool.h>
#include <inttypes.h>

struct output_stream;
struct input_stream;
struct read;

struct read_format_params {
	enum {
		READ_FORMAT_FASTQ,
		READ_FORMAT_TAB_DELIMITED,
	} fmt;
	int phred_offset;
};

/* Returns true iff the specified read format supports both unpaired and paired
 * reads in the same file.  */
static inline bool
read_format_supports_mixed_reads(const struct read_format_params *params)
{
	return params->fmt == READ_FORMAT_TAB_DELIMITED;
}

extern void
write_read(struct output_stream *out,
	   const struct read_format_params *oparams,
	   struct read *read);

extern void
write_read_pair(struct output_stream *out,
		const struct read_format_params *oparams,
		struct read *read_1, struct read *read_2);

extern bool
load_read(struct input_stream *in,
	  const struct read_format_params *iparams,
	  struct read *read, uint64_t *line_no_p);

extern bool
load_read_pair(struct input_stream *in,
	       const struct read_format_params *iparams,
	       struct read *read_1, struct read *read_2,
	       uint64_t *line_no_p);

#endif /* _FLASH_READ_IO_H_  */
