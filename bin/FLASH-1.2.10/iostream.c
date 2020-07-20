/*
 * iostream.c:  Buffered file input and output streams with automatic error
 * checking and extra features such as gzip compression.
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
#include "util.h"

#include <assert.h>
#include <errno.h>
#include <fcntl.h>
#include <limits.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>

/* Return true iff the specified string is a single hyphen, which may represent
 * standard input or standard output.  */
static bool
string_is_hyphen(const char *str)
{
	return str[0] == '-' && str[1] == '\0';
}

static bool
mode_is_writing(const char *mode)
{
	return strchr(mode, 'w') != NULL;
}

static bool
flags_is_writing(int flags)
{
	int accmode = (flags & O_ACCMODE);
	return (accmode == O_WRONLY || accmode == O_RDWR);
}

static const char *
access_mode_string(const char *mode)
{
	return mode_is_writing(mode) ? "writing" : "reading";
}

static const char *
access_flags_string(int flags)
{
	return flags_is_writing(flags) ? "writing" : "reading";
}

static int
mode_to_flags(const char *mode)
{
	return mode_is_writing(mode) ? O_WRONLY | O_TRUNC | O_CREAT : O_RDONLY;
}

/* Returns a standard file descriptor depending on the flags:
 *
 *	reading => standard input
 *	writing => standard output
 */
static int
standard_fd_from_flags(int flags)
{
	return flags_is_writing(flags) ? STDOUT_FILENO : STDIN_FILENO;
}

/* Like fopen(), but aborts on error.  */
void *
xfopen(const char *path, const char *mode)
{
	FILE *fp = fopen(path, mode);
	if (!fp)
		fatal_error_with_errno("Failed to open \"%s\" for %s",
				       path, access_mode_string(mode));
	return fp;
}

/* Like fclose(), but aborts on error.  */
void
xfclose(FILE *fp, const char *name)
{
	if (fclose(fp))
		fatal_error_with_errno("Error closing \"%s\"", name);
}

#ifndef O_BINARY
#  define O_BINARY 0
#endif

/* Like open(), but aborts on error, and also interprets "-" as standard output
 * or standard input rather than a path, depending on the requested mode.  Also
 * automatically provides O_BINARY on Windows.  */
static void *
xopen(const char *path, int flags, mode_t mode)
{
	int fd;

	if (string_is_hyphen(path)) {
		fd = standard_fd_from_flags(flags);
	#ifdef __WIN32__
		_setmode(fd, O_BINARY);
	#endif
	} else {
		fd = open(path, flags | O_BINARY, mode);
	}

	if (fd < 0)
		fatal_error_with_errno("Failed to open \"%s\" for %s",
				       path, access_flags_string(flags));

	/* XXX: autodetect whether posix_fadvise() is available or not.  */
#ifdef __linux__
	/* Advise the operating system that the file will be read or written
	 * sequentially.  */
	posix_fadvise(fd, 0, 0, POSIX_FADV_SEQUENTIAL);
#endif

	return (void*)(intptr_t)fd;
}

/* Like gzopen(), but aborts on error, and also interprets "-" as standard
 * output or standard input rather than a path, depending on the requested mode.
 */
static void *
xgzopen(const char *path, const char *mode)
{
	gzFile gzf;
	int fd;

	fd = (int)(intptr_t)xopen(path, mode_to_flags(mode), 0644);
	errno = 0;
	gzf = gzdopen(fd, mode);
	if (!gzf)
		fatal_error_with_errno("Failed to open \"%s\" for %s",
				       path, access_mode_string(mode));
	return gzf;
}

/* Like gzclose(), but aborts on error.  */
static void
xgzclose(void *_fp, const char *name)
{
	gzFile gzf = (gzFile)_fp;
	errno = 0;
	if (gzclose(gzf) != Z_OK)
		fatal_error_with_errno("Error closing \"%s\"", name);
}

/* Checks the error status on a gzFile.  Returns only if no error occurred or if
 * the error was caused by an interruption (EINTR).  */
static void
check_gzerror(gzFile gzf, bool was_writing, const char *name)
{
	int errnum;
	const char *error_str;

	error_str = gzerror(gzf, &errnum);
	switch (errnum) {
	case Z_OK:
		return;
	case Z_ERRNO:
		if (errno == EINTR)
			return;
		fatal_error_with_errno("Error %s \"%s\"",
				       (was_writing ? "writing" : "reading"),
				       name);
	default:
		fatal_error("zlib error %s \"%s\": %s",
			    (was_writing ? "writing" : "reading"),
			    name, error_str);
	}
}

/* Reads data from gzFile, with error checking.  Returns the number of bytes
 * successfully read, which will be less than or equal to @count; 0 implies the
 * stream is at end-of-file.  Aborts on error.  */
static size_t
xgzread(void *_fp, void *buf, size_t count, const char *name)
{
	gzFile gzf = (gzFile)_fp;
	int trycount = min(count, INT_MAX);

	for (;;) {
		int ret = gzread(gzf, buf, trycount);
		if (ret >= 0)
			return ret;
		if (gzeof(gzf))
			return 0;
		check_gzerror(gzf, false, name);
	}
}

/* Writes data to a gzFile, with error checking.  */
static void
xgzwrite(void *_fp, const void *_buf, size_t count, const char *name)
{
	gzFile gzf = (gzFile)_fp;
	const char *ptr = _buf;

	while (count) {
		int trycount = min(count, INT_MAX);
		int ret = gzwrite(gzf, ptr, trycount);
		if (!ret) {
			check_gzerror(gzf, true, name);
			continue;
		}
		ptr += ret;
		count -= ret;
	}
}

/* Reads data from a file descriptor, with error checking.  Returns the number
 * of bytes successfully read, which will be less than or equal to @count; 0
 * implies end-of-file has been reached.  Aborts on error.  */
static size_t
xread(void *_fp, void *buf, size_t count, const char *name)
{
	int fd = (int)(intptr_t)_fp;
	ssize_t trycount = min(count, SSIZE_MAX);

	for (;;) {
		ssize_t ret = read(fd, buf, trycount);
		if (ret >= 0)
			return ret;
		if (errno != EINTR)
			fatal_error_with_errno("Error reading \"%s\"", name);
	}
}

/* Similar to xread(), but retries on short reads.  */
static size_t
full_xread(void *fp, void *_buf, size_t count, const char *name)
{
	char *ptr = _buf;
	while (count) {
		size_t ret = xread(fp, ptr, count, name);
		if (ret == 0)
			break;
		ptr += ret;
		count -= ret;
	}
	return ptr - (char*)_buf;
}

/* Writes data to a file descriptor, with error checking.  */
static void
xwrite(void *_fp, const void *_buf, size_t count, const char *name)
{
	int fd = (int)(intptr_t)_fp;
	const char *ptr = _buf;

	while (count) {
		ssize_t trycount = min(count, SSIZE_MAX);
		ssize_t ret = write(fd, ptr, trycount);
		if (ret < 0) {
			if (errno == EINTR)
				continue;
			fatal_error_with_errno("Error writing \"%s\"", name);
		}
		ptr += ret;
		count -= ret;
	}
}

/* Closes a file descriptor, with error checking.  */
static void
xclose(void *_fp, const char *name)
{
	int fd = (int)(intptr_t)_fp;
	if (close(fd))
		fatal_error_with_errno("Error closing \"%s\"", name);
}

/* Writes data to a FILE *, with error checking.  */
static void
xfwrite(void *_fp, const void *buf, size_t count, const char *name)
{
	FILE *fp = (FILE*)_fp;

	if (fwrite(buf, 1, count, fp) != count)
		fatal_error_with_errno("Error writing \"%s\"", name);
}

/* Closes a FILE * opened with popen(), with error checking.  */
static void
xpclose(void *_fp, const char *name)
{
	FILE *fp = (FILE*)_fp;
	int status = pclose(fp);

	if (status == -1)
		fatal_error_with_errno("Error closing pipe to \"%s\"", name);
	if (status)
		fatal_error("Program writing to \"%s\" "
			    "exited with failure status", name);
}

#define IOSTREAM_BUFSIZE 32768

/*****************************
 * Input stream functions    *
 *****************************/

static void *
input_gzfile_open(const char *path)
{
	return xgzopen(path, "rb");
}

static void *
input_fd_open(const char *path)
{
	return xopen(path, O_RDONLY, 0);
}

struct input_stream {
	const struct input_stream_operations *ops;
	void *fp;
	char *name;
	char *buf_begin;
	char *buf_end;
	char *buf_cur_begin;
	char *buf_cur_end;
};

struct input_stream_operations {
	/* Open the specified file for reading.  */
	void * (*open)(const char *path);

	/* Read data from the specified file (up to @count bytes;
	 * return 0 if at end of file).  */
	size_t (*read)(void *fp, void *buf, size_t count, const char *name);

	/* Close the specified file.  */
	void (*close)(void *fp, const char *name);
};

/* Input stream operations for reading raw data from a file descriptor  */
static const struct input_stream_operations fd_input_stream_ops = {
	.open  = input_fd_open,
	.read  = xread,
	.close = xclose,
};

/* Input stream operations for reading a gzip compressed file using zlib  */
static const struct input_stream_operations gzfile_input_stream_ops = {
	.open  = input_gzfile_open,
	.read  = xgzread,
	.close = xgzclose,
};

/* Auto-detects the correct input_stream_operations to use for the specified
 * file.  */
static const struct input_stream_operations *
select_input_stream_ops(const char *path)
{
	/* XXX:  We can't rewind standard input after checking for magic bytes,
	 * so for that case we rely on the fact that gzread returns the literal
	 * data if the stream does not, in fact, contain gzipped data.  */
	if (string_is_hyphen(path))
		return &gzfile_input_stream_ops;

	/* Test for gzip magic bytes { 0x1f, 0x8b}  */

	unsigned char magic[2] = {0, 0};
	void *tmp_fp;

	tmp_fp = input_fd_open(path);
	full_xread(tmp_fp, magic, sizeof(magic), path);
	xclose(tmp_fp, path);

	if (magic[0] == 0x1f && magic[1] == 0x8b)
		return &gzfile_input_stream_ops;

	/* Default to reading the raw data  */

	return &fd_input_stream_ops;
}

/* Creates an input stream to read lines from the file specified by @path.
 *
 * Gzip files are auto-detected.  */
struct input_stream *
new_input_stream(const char *path)
{
	struct input_stream *in = xmalloc(sizeof(*in));

	assert(path != NULL);

	/* Select input_stream_operations and open stream  */
	in->ops = select_input_stream_ops(path);
	in->fp = (*in->ops->open)(path);
	in->name = xstrdup(path);

	/* Allocate internal buffer  */
	in->buf_begin     = xmalloc(IOSTREAM_BUFSIZE);
	in->buf_end       = in->buf_begin + IOSTREAM_BUFSIZE;
	in->buf_cur_begin = in->buf_begin;
	in->buf_cur_end   = in->buf_begin;

	return in;
}

/* Returns the name of the file being read.  */
const char *
input_stream_get_name(struct input_stream *in)
{
	return in->name;
}

/* Returns a pointer to the next instance in any of the @delims in @buf of
 * length @size.
 *
 * There must be at least one delimiter.  */
static char *
find_delim(char *buf, size_t size, const char *delims)
{
	/* Fast case: just one delimiter.  */
	if (delims[1] == '\0')
		return memchr(buf, delims[0], size);

	/* Multiple delimiters.  */
	for (size_t i = 0; i < size; i++) {
		const char *p = delims;
		do {
			if (buf[i] == *p)
				return &buf[i];
		} while (*++p);
	}
	return NULL;
}

/* Reads delimited data from an input stream.  Semantics are like getdelim(),
 * but aborts on read error and also allows multiple delimiters.  */
ssize_t
input_stream_getdelims(struct input_stream *in, char **lineptr, size_t *n,
		       const char *delims)
{
	assert(*delims != '\0');

	/* offset = number of bytes copied to *lineptr buffer, excluding
	 *          terminating null byte  */
	size_t offset = 0;

	for (;;) {
		size_t navail;
		char *delim_ptr;
		size_t copysize;

		navail = in->buf_cur_end - in->buf_cur_begin;

		if (navail == 0) {
			/* No more data in internal buffer; try to fill it  */

			in->buf_cur_begin = in->buf_begin;
			navail = (*in->ops->read)(in->fp,
						  in->buf_cur_begin,
						  in->buf_end - in->buf_cur_begin,
						  in->name);
			in->buf_cur_end = in->buf_cur_begin + navail;

			if (navail == 0)  /* At end-of-file  */
				break;
		}

		/* Find the first delimiter in the internal buffer.  If found,
		 * copy up to and including the delimiter, then return (break
		 * loop).  If not found, copy all the data and try to read more
		 * (continue loop).  */
		delim_ptr = find_delim(in->buf_cur_begin, navail, delims);
		if (delim_ptr)
			copysize = delim_ptr - in->buf_cur_begin + 1;
		else
			copysize = navail;

		if (offset + copysize + 1 < offset) {
			/* Very unlikely: size would overflow.  */
			fatal_error("Line or field in \"%s\" is too long!",
				    in->name);
		}

		if (*n < offset + copysize + 1) {
			*n = max(*n * 3 / 2, offset + copysize + 1);
			*n = max(*n, 128);
			*lineptr = xrealloc(*lineptr, *n);
		}

		memcpy(*lineptr + offset, in->buf_cur_begin, copysize);
		offset += copysize;
		in->buf_cur_begin += copysize;
		if (delim_ptr)
			break;
	}

	if (offset == 0)
		return -1;

	(*lineptr)[offset] = '\0';
	return offset;
}

/* Reads a line from an input stream.  Semantics are like getline(), but aborts
 * on read error.  */
ssize_t
input_stream_getline(struct input_stream *in, char **lineptr, size_t *n)
{
	return input_stream_getdelims(in, lineptr, n, "\n");
}

/* Closes and frees an input stream.  */
void
free_input_stream(struct input_stream *in)
{
	if (in) {
		(*in->ops->close)(in->fp, in->name);
		xfree(in->name, strlen(in->name));
		xfree(in->buf_begin, IOSTREAM_BUFSIZE);
		xfree(in, sizeof(*in));
	}
}

/*****************************
 * Output stream functions   *
 *****************************/

static void *
output_fd_open(const char *path, const char *filter_prog,
	       const char *filter_prog_args)
{
	assert(filter_prog == NULL);
	return xopen(path, O_WRONLY | O_CREAT | O_TRUNC, 0644);
}

static void *
output_filter_open(const char *path, const char *filter_prog,
		   const char *filter_prog_args)
{
	assert(filter_prog != NULL);

	if (filter_prog_args == NULL)
		filter_prog_args = "";

	size_t len = strlen(filter_prog) + 32 +
		     strlen(path) + strlen(filter_prog_args);
	char command[len + 1];
	FILE *f;
	char *p = command;

	p += sprintf(p, "%s %s", filter_prog, filter_prog_args);
	if (!string_is_hyphen(path))
		sprintf(p, " > '%s'", path);

	f = popen(command, "w");
	if (!f)
		fatal_error_with_errno("Failed to launch the command \"%s\"",
				       command);
	return f;
}

static void *
output_gzfile_open(const char *path, const char *filter_prog,
		   const char *filter_prog_args)
{
	assert(filter_prog == NULL);
	return xgzopen(path, "wb");
}

struct output_stream_operations {
	/* Open the specified file for writing, possibly filtering the data
	 * through a filter program.  Must return the open file handle,
	 * descriptor, or pointer cast to a void *.  */
	void *(*open)(const char *path,
		      const char *filter_prog, const char *filter_prog_args);

	/* Write a buffer of data to the stream.  */
	void (*write)(void *fp, const void *buf, size_t count,
		      const char *name);

	/* Flush and close the stream.  */
	void (*close)(void *fp, const char *name);
};

/* Operations to write raw data to a file descriptor.  */
static const struct output_stream_operations fd_output_stream_ops = {
	.open  = output_fd_open,
	.write = xwrite,
	.close = xclose,
};

/* Operations to write data through a filter program to the file.  */
static const struct output_stream_operations filter_output_stream_ops = {
	.open  = output_filter_open,
	.write = xfwrite,
	.close = xpclose,
};

/* Operations to write gzip-compressed data.  */
static const struct output_stream_operations gzfile_output_stream_ops = {
	.open  = output_gzfile_open,
	.write = xgzwrite,
	.close = xgzclose,
};

struct output_stream {
	const struct output_stream_operations *ops;
	void *fp;
	char *name;
	char *buf_begin;
	char *buf_end;
	char *buf_cur_end;
};

/* Select the appropriate output_stream_operations based on the requested
 * compression type and whether a filter program was specified.  */
static const struct output_stream_operations *
select_output_stream_ops(enum out_compression_type ctype, bool have_filter_prog)
{
	switch (ctype) {
	case OUT_COMPRESSION_NONE:
		if (have_filter_prog)
			return &filter_output_stream_ops;
		else
			return &fd_output_stream_ops;
	case OUT_COMPRESSION_GZIP:
		assert(!have_filter_prog);
		return &gzfile_output_stream_ops;
	default:
		assert(0);
		return NULL;
	}
}

/*
 * Creates a new output stream.
 *
 * @ctype
 *	The compression type to use.
 * @path
 *	Path to the file to write, or "-" for standard output.
 * @filter_prog
 *	If non-NULL, the name of a program through which to filter the output,
 *	and @ctype must be OUT_COMPRESSION_NONE.  The filter program must read
 *	from standard input and write to standard output.
 * @filter_prog_args
 *	Additional text (arguments) to place on the command line for
 *	@filter_prog.
 *
 * Returns the new, opened output stream; aborts on error.
 */
struct output_stream *
new_output_stream(enum out_compression_type ctype,
		  const char *path,
		  const char *filter_prog,
		  const char *filter_prog_args)
{
	struct output_stream *out = xmalloc(sizeof(*out));

	assert(path != NULL);

	/* Select output_stream_operations and open stream.  */
	out->ops = select_output_stream_ops(ctype, filter_prog != NULL);
	out->fp = (*out->ops->open)(path, filter_prog, filter_prog_args);
	out->name = xstrdup(path);

	/* Allocate internal buffer.  */
	out->buf_begin = xmalloc(IOSTREAM_BUFSIZE);
	out->buf_end = out->buf_begin + IOSTREAM_BUFSIZE;
	out->buf_cur_end = out->buf_begin;

	return out;
}

/* Returns the name of the file to which the output stream is writing.  */
const char *
output_stream_get_name(struct output_stream *out)
{
	return out->name;
}

static void
flush_output_stream(struct output_stream *out)
{
	(*out->ops->write)(out->fp, out->buf_begin,
			   out->buf_cur_end - out->buf_begin, out->name);
	out->buf_cur_end = out->buf_begin;
}

/* Writes a buffer of data to an output stream.  */
void
output_stream_write(struct output_stream *out, const void *_buf, size_t count)
{
	const char *ptr = _buf;

	while (count) {
		if (out->buf_cur_end == out->buf_end) {
			/* Output buffer full; flush it.  */
			flush_output_stream(out);
		}

		/* Buffer as much data as possible.  */

		size_t tocopy = min(count, out->buf_end - out->buf_cur_end);

		memcpy(out->buf_cur_end, ptr, tocopy);
		out->buf_cur_end += tocopy;
		ptr += tocopy;
		count -= tocopy;
	}
}

/* Writes a null-terminated string to an output stream.  */
void
output_stream_fputs(struct output_stream *out, const char *s)
{
	output_stream_write(out, s, strlen(s));
}

/* Writes a byte to an output stream.  */
void
output_stream_fputc(struct output_stream *out, char c)
{
	if (out->buf_cur_end == out->buf_end)
		flush_output_stream(out);
	*out->buf_cur_end++ = c;
}

/* Flushes, closes, and frees an output stream.  */
void
free_output_stream(struct output_stream *out)
{
	if (out) {
		if (out->buf_cur_end > out->buf_begin)
			flush_output_stream(out);
		(*out->ops->close)(out->fp, out->name);
		xfree(out->name, strlen(out->name));
		xfree(out->buf_begin, IOSTREAM_BUFSIZE);
		xfree(out, sizeof(*out));
	}
}
