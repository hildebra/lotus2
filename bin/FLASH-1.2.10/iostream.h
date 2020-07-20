#ifndef _FLASH_IOSTREAM_H_
#define _FLASH_IOSTREAM_H_

#include <stdio.h>
#include <stddef.h>
#include <sys/types.h>

/* Input stream functions  */

struct input_stream;

extern struct input_stream *
new_input_stream(const char *filename);

extern const char *
input_stream_get_name(struct input_stream *in);

extern ssize_t
input_stream_getdelims(struct input_stream *in, char **lineptr,
		       size_t *n, const char *delims);
extern ssize_t
input_stream_getline(struct input_stream *in, char **lineptr, size_t *n);

extern void
free_input_stream(struct input_stream *in);

/* Output stream functions  */

struct output_stream;

enum out_compression_type {
	OUT_COMPRESSION_NONE,
	OUT_COMPRESSION_GZIP,
};

extern struct output_stream *
new_output_stream(enum out_compression_type ctype,
		  const char *path,
		  const char *filter_prog,
		  const char *filter_prog_args);
extern const char *
output_stream_get_name(struct output_stream *out);

extern void
output_stream_write(struct output_stream *out,
		    const void *buf, size_t count);

extern void
output_stream_fputs(struct output_stream *out, const char *s);

extern void
output_stream_fputc(struct output_stream *out, char c);

extern void
free_output_stream(struct output_stream *out);

/* fopen() and fclose() wrappers  */

extern void *
xfopen(const char *filename, const char *mode);

extern void
xfclose(FILE *fp, const char *name);

#endif /* _FLASH_IOSTREAM_H_  */
