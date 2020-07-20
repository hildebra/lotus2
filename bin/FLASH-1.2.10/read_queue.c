/*
 * read_queue.c: Code to set up reader/writer threads and shared queues to pass
 * reads between threads in memory.
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
#include "read_queue.h"
#include "util.h"

#include <assert.h>
#include <errno.h>
#include <inttypes.h>
#include <pthread.h>
#include <stdlib.h>
#include <string.h>

static struct read *
new_read(void)
{
	return xzalloc(sizeof(struct read));
}

static void
free_read(struct read *r)
{
	if (r) {
		xfree(r->tag, r->tag_bufsz);
		xfree(r->seq, r->seq_bufsz);
		xfree(r->qual, r->qual_bufsz);
		xfree(r, sizeof(*r));
	}
}

static struct read_set *
new_read_set(size_t num_reads, bool full)
{
	struct read_set *s = xmalloc(sizeof(*s) + num_reads * sizeof(s->reads[0]));
	if (full) {
		for (size_t i = 0; i < num_reads; i++)
			s->reads[i] = new_read();
	} else {
		for (size_t i = 0; i < num_reads; i++)
			s->reads[i] = NULL;
	}
	s->filled = 0;
	s->num_reads = num_reads;
	return s;
}

void
free_read_set(struct read_set *s)
{
	if (s) {
		for (size_t i = 0; i < s->num_reads; i++)
			free_read(s->reads[i]);
		xfree(s, sizeof(*s));
	}
}

static void
init_mutex(pthread_mutex_t *mutex)
{
	if (pthread_mutex_init(mutex, NULL))
		fatal_error_with_errno("Failed to initialize mutex");
}

static void
init_cond(pthread_cond_t *cond)
{
	if (pthread_cond_init(cond, NULL))
		fatal_error_with_errno("Failed to initialize condition variable");
}

/*
 * Producer-consumer queue; it holds pointers to `struct read_sets', which can
 * be added or removed from the queue in a thread-safe manner using
 * read_queue_put() and read_queue_get(), respectively.
 */
struct read_queue {
	size_t size;
	size_t front;
	size_t filled;
	bool terminated;
	struct read_set **read_sets;
	pthread_mutex_t lock;
	pthread_cond_t read_set_avail_cond;
	pthread_cond_t space_avail_cond;
};

static struct read_queue *
new_read_queue(size_t size, size_t reads_per_set, bool full)
{
	struct read_queue *q = xmalloc(sizeof(*q));

	q->read_sets = xmalloc(size * sizeof(q->read_sets[0]));
	q->size = size;
	q->front = 0;
	if (full) {
		for (size_t i = 0; i < size; i++)
			q->read_sets[i] = new_read_set(reads_per_set, true);
		q->filled = size;
	} else {
		for (size_t i = 0; i < size; i++)
			q->read_sets[i] = NULL;
		q->filled = 0;
	}
	q->terminated = false;
	init_mutex(&q->lock);
	init_cond(&q->read_set_avail_cond);
	init_cond(&q->space_avail_cond);
	return q;
}

static void
free_read_queue(struct read_queue *q)
{
	if (q) {
		size_t filled = q->filled;
		size_t i = q->front;

		while (filled--) {
			free_read_set(q->read_sets[i]);
			i = (i + 1) % q->size;
		}

		xfree(q->read_sets, q->size * sizeof(q->read_sets[0]));

		pthread_mutex_destroy(&q->lock);
		pthread_cond_destroy(&q->read_set_avail_cond);
		pthread_cond_destroy(&q->space_avail_cond);

		xfree(q, sizeof(*q));
	}
}

/* Retrieves the next available read set from the queue, blocking until one is
 * available.  Or, returns NULL if the queue has terminated and no more read
 * sets are available.  */
static struct read_set *
read_queue_get(struct read_queue *q)
{
	struct read_set *s;

	pthread_mutex_lock(&q->lock);
	while (q->filled == 0 && !q->terminated)
		pthread_cond_wait(&q->read_set_avail_cond, &q->lock);

	if (q->filled != 0) {
		s = q->read_sets[q->front];
		q->front = (q->front + 1) % q->size;
		q->filled--;
		pthread_cond_signal(&q->space_avail_cond);
	} else
		s = NULL;

	pthread_mutex_unlock(&q->lock);
	return s;
}

/* Put a read set into the queue, blocking until there is an empty space
 * available. */
static void
read_queue_put(struct read_queue *q, struct read_set *s)
{
	pthread_mutex_lock(&q->lock);
	while (q->filled == q->size)
		pthread_cond_wait(&q->space_avail_cond, &q->lock);

	q->read_sets[(q->front + q->filled) % q->size] = s;
	q->filled++;

	pthread_cond_signal(&q->read_set_avail_cond);
	pthread_mutex_unlock(&q->lock);
}

/* "Terminate" the specified queue.  This will cause read_queue_get() to return
 * NULL once the queue is empty.  */
static void
read_queue_terminate(struct read_queue *q)
{
	pthread_mutex_lock(&q->lock);
	q->terminated = true;
	pthread_cond_broadcast(&q->read_set_avail_cond);
	pthread_mutex_unlock(&q->lock);
}

struct reader_params {
	struct input_stream *in;
	const struct read_format_params *iparams;
	bool verbose;
	struct read_queue *avail_read_q;
	struct read_queue *unprocessed_read_1_q;
	struct read_queue *unprocessed_read_2_q;
	struct read_queue *unpaired_read_q;
};

struct writer_params {
	struct output_stream *out;
	const struct read_format_params *oparams;
	struct read_queue *to_write_queue_1;
	struct read_queue *to_write_queue_2;
	struct read_queue *avail_queue;
};

static void
processed(uint64_t pair_no)
{
	info("Processed %"PRIu64" read pairs", pair_no);
}

static void *
reader1_proc(void *_params)
{
	struct reader_params *params = _params;
	uint64_t pair_no = 0;
	uint64_t line_no = 1;
	struct read_set *s;

	for (;;) {

		s = read_queue_get(params->avail_read_q);

		for (s->filled = 0;
		     s->filled < s->num_reads;
		     s->filled++)
		{
			if (!load_read(params->in, params->iparams,
				       s->reads[s->filled], &line_no))
				goto eof_reached;

			if (params->verbose && ++pair_no % 25000 == 0)
				processed(pair_no);
		}

		/* Note: although we're placing the set in
		 * 'unprocessed_read_1_q', the set may in fact be read 2, not
		 * read 1.  This procedure works the same way in both cases.  */

		read_queue_put(params->unprocessed_read_1_q, s);
	}

eof_reached:
	if (params->verbose && pair_no % 25000 != 0)
		processed(pair_no);

	if (s->filled)
		read_queue_put(params->unprocessed_read_1_q, s);
	else
		free_read_set(s);

	read_queue_terminate(params->unprocessed_read_1_q);

	free_input_stream(params->in);
	xfree(params, sizeof(*params));
	return NULL;
}

static void *
reader2_proc(void *_params)
{
	struct reader_params *params = _params;
	struct read_set *s_read1, *s_read2, *s_unpaired = NULL;
	uint64_t pair_no = 0;
	uint64_t line_no = 1;

	s_read1 = read_queue_get(params->avail_read_q);
	s_read1->filled = 0;
	s_read2 = read_queue_get(params->avail_read_q);
	s_read2->filled = 0;
	if (params->unpaired_read_q) {
		s_unpaired = read_queue_get(params->avail_read_q);
		s_unpaired->filled = 0;
	}

	while (load_read_pair(params->in, params->iparams,
			      s_read1->reads[s_read1->filled],
			      s_read2->reads[s_read1->filled],
			      &line_no))
	{
		if (s_read2->reads[s_read1->filled]->seq_len) {
			/* Read pair.  */
			++s_read1->filled;
			++s_read2->filled;
			if (s_read1->filled == s_read1->num_reads) {
				read_queue_put(params->unprocessed_read_1_q, s_read1);
				read_queue_put(params->unprocessed_read_2_q, s_read2);
				s_read1 = read_queue_get(params->avail_read_q);
				s_read1->filled = 0;
				s_read2 = read_queue_get(params->avail_read_q);
				s_read2->filled = 0;
			}
			if (params->verbose && ++pair_no % 25000 == 0)
				processed(pair_no);
		} else if (params->unpaired_read_q) {
			/* Actually an unpaired read.  */
			struct read *r = s_read1->reads[s_read1->filled];
			s_read1->reads[s_read1->filled] =
				s_unpaired->reads[s_unpaired->filled];
			s_unpaired->reads[s_unpaired->filled] = r;
			++s_unpaired->filled;
			if (s_unpaired->filled == s_unpaired->num_reads) {
				s_unpaired->type = READS_UNPAIRED;
				read_queue_put(params->unpaired_read_q, s_unpaired);
				s_unpaired = read_queue_get(params->avail_read_q);
				s_unpaired->filled = 0;
			}
		}
	}

	if (params->verbose && pair_no % 25000 != 0)
		processed(pair_no);

	if (s_read1->filled)
		read_queue_put(params->unprocessed_read_1_q, s_read1);
	else
		free_read_set(s_read1);

	if (s_read2->filled)
		read_queue_put(params->unprocessed_read_2_q, s_read2);
	else
		free_read_set(s_read2);

	if (s_unpaired) {
		if (s_unpaired->filled) {
			s_unpaired->type = READS_UNPAIRED;
			read_queue_put(params->unpaired_read_q, s_unpaired);
		} else {
			free_read_set(s_unpaired);
		}
	}

	read_queue_terminate(params->unprocessed_read_1_q);
	read_queue_terminate(params->unprocessed_read_2_q);

	free_input_stream(params->in);
	xfree(params, sizeof(*params));
	return NULL;
}

static void *
writer_proc(void *_params)
{
	struct writer_params *params = _params;
	struct read_set *s1, *s2;

	for (;;) {
		s1 = read_queue_get(params->to_write_queue_1);
		if (!s1)
			break;

		if (params->to_write_queue_2 && s1->type == READS_UNCOMBINED) {
			/* Get other read in uncombined pair  */
			s2 = read_queue_get(params->to_write_queue_2);
			assert(s2);
			assert(s1->filled == s2->filled);
		} else {
			s2 = NULL;
		}

		for (size_t i = 0; i < s1->filled; i++) {
			if (s2)
				write_read_pair(params->out, params->oparams,
						s1->reads[i], s2->reads[i]);
			else
				write_read(params->out, params->oparams,
					   s1->reads[i]);
		}
		read_queue_put(params->avail_queue, s1);
		if (s2)
			read_queue_put(params->avail_queue, s2);
	}
	free_output_stream(params->out);
	xfree(params, sizeof(*params));
	return NULL;
}

static pthread_t
start_reader2(struct input_stream *in,
	      const struct read_format_params *iparams,
	      bool verbose,
	      struct read_queue *avail_read_q,
	      struct read_queue *unprocessed_read_1_q,
	      struct read_queue *unprocessed_read_2_q,
	      struct read_queue *unpaired_read_q)
{
	struct reader_params *params = xmalloc(sizeof(*params));

	params->in = in;
	params->iparams = iparams;
	params->verbose = verbose;
	params->avail_read_q = avail_read_q;
	params->unprocessed_read_1_q = unprocessed_read_1_q;
	params->unprocessed_read_2_q = unprocessed_read_2_q;
	params->unpaired_read_q = unpaired_read_q;

	return create_thread(reader2_proc, params);
}

static pthread_t
start_reader1(struct input_stream *in,
	      const struct read_format_params *iparams,
	      bool verbose,
	      struct read_queue *avail_read_q,
	      struct read_queue *unprocessed_read_q)
{
	struct reader_params *params = xmalloc(sizeof(*params));

	params->in = in;
	params->iparams = iparams;
	params->verbose = verbose;
	params->avail_read_q = avail_read_q;
	params->unprocessed_read_1_q = unprocessed_read_q;
	params->unprocessed_read_2_q = NULL;
	params->unpaired_read_q = NULL;

	return create_thread(reader1_proc, params);
}

static pthread_t
start_writer2(struct output_stream *out,
	      const struct read_format_params *oparams,
	      struct read_queue *to_write_queue_1,
	      struct read_queue *to_write_queue_2,
	      struct read_queue *avail_queue)
{
	struct writer_params *params = xmalloc(sizeof(*params));

	params->out = out;
	params->oparams = oparams;
	params->to_write_queue_1 = to_write_queue_1;
	params->to_write_queue_2 = to_write_queue_2;
	params->avail_queue = avail_queue;

	return create_thread(writer_proc, params);
}

static pthread_t
start_writer1(struct output_stream *out,
	      const struct read_format_params *oparams,
	      struct read_queue *to_write_queue,
	      struct read_queue *avail_queue)
{
	return start_writer2(out, oparams, to_write_queue, NULL, avail_queue);
}


struct read_io_handle {

	pthread_t reader_1;
	pthread_t reader_2;
	pthread_t writer_1;
	pthread_t writer_2;
	pthread_t writer_3;
	bool reader_1_started;
	bool reader_2_started;
	bool writer_1_started;
	bool writer_2_started;
	bool writer_3_started;

	unsigned combiner_threads_remaining;
	pthread_mutex_t combiner_threads_remaining_mutex;

	struct read_queue *avail_read_q;
	struct read_queue *unprocessed_read_1_q;
	struct read_queue *unprocessed_read_2_q;

	struct read_queue *combined_read_q;
	struct read_queue *uncombined_read_1_q;
	struct read_queue *uncombined_read_2_q;

	pthread_mutex_t get_unprocessed_pair_mutex;
	pthread_mutex_t put_uncombined_pair_mutex;
};

/* Retrieves some unprocessed read pairs from the I/O layer.  Returns %true iff
 * more reads were available; returns false if end of file was reached.  */
bool
get_unprocessed_read_pairs(struct read_io_handle *h, struct read_set **s1_p,
			   struct read_set **s2_p)
{
	/* get_unprocessed_pair_mutex ensures the reads are paired up correctly.
	 */
	struct read_set *s1, *s2;

	pthread_mutex_lock(&h->get_unprocessed_pair_mutex);

	s1 = read_queue_get(h->unprocessed_read_1_q);
	s2 = read_queue_get(h->unprocessed_read_2_q);

	pthread_mutex_unlock(&h->get_unprocessed_pair_mutex);

	if (s1 && s2) {
		if (s1->filled != s2->filled)
			goto mismatch;
		*s1_p = s1;
		*s2_p = s2;
		return true;
	}

	if (s1 || s2)
		goto mismatch;
	return false;

mismatch:
	fatal_error("Input files do not contain the same number of reads");
}

/* Submits a set of combined reads to the I/O layer to be written.  */
void
put_combined_reads(struct read_io_handle *h, struct read_set *s)
{
	s->type = READS_COMBINED;

	read_queue_put(h->combined_read_q, s);
}

/* Submits a set of uncombined read pairs to the I/O layer to be written.  */
void
put_uncombined_read_pairs(struct read_io_handle *h,
			  struct read_set *s1, struct read_set *s2)
{
	s1->type = READS_UNCOMBINED;
	s2->type = READS_UNCOMBINED;

	/* put_unprocessed_pair_mutex ensures the reads are paired up correctly.
	 */

	pthread_mutex_lock(&h->put_uncombined_pair_mutex);

	read_queue_put(h->uncombined_read_1_q, s1);
	read_queue_put(h->uncombined_read_2_q, s2);

	pthread_mutex_unlock(&h->put_uncombined_pair_mutex);
}


/* Retrieve a read set (full of read structures) that is ready to be reused.  */
struct read_set *
get_avail_read_set(struct read_io_handle *h)
{
	struct read_set *s;

	s = read_queue_get(h->avail_read_q);
	s->filled = 0;
	return s;
}

/* Return a set of read pairs to the pool for reuse.  */
void
put_avail_read_pairs(struct read_io_handle *h,
		     struct read_set *s1, struct read_set *s2)
{
	read_queue_put(h->avail_read_q, s1);
	read_queue_put(h->avail_read_q, s2);
}

/* Notify the I/O layer that a combiner thread has terminated.
 * When all the combiner threads have been terminated, the writers will shut
 * down.  */
void
notify_combiner_terminated(struct read_io_handle *h)
{
	pthread_mutex_lock(&h->combiner_threads_remaining_mutex);

	if (--h->combiner_threads_remaining == 0) {

		/* Terminate the writer queues.  */

		read_queue_terminate(h->combined_read_q);

		if (h->uncombined_read_1_q != h->avail_read_q &&
		    h->uncombined_read_1_q != h->combined_read_q)
			read_queue_terminate(h->uncombined_read_1_q);

		if (h->uncombined_read_2_q != h->avail_read_q)
			read_queue_terminate(h->uncombined_read_2_q);
	}

	pthread_mutex_unlock(&h->combiner_threads_remaining_mutex);
}

struct read_set *
new_empty_read_set(struct read_io_handle *h)
{
	return new_read_set(BASE_READS_PER_READ_SET +
			    (h->combiner_threads_remaining * PERTHREAD_READS_PER_READ_SET),
			    false);
}

/* Starts the FLASH I/O layer, which is responsible for input/output of reads.
 *
 * If @in_2 is not NULL, then @in_1 and @in_2 are the input files for read 1 and
 * read 2 of the pairs, respectively.  Otherwise @in_1 contains both read 1 and
 * read 2 of the pairs interleaved.
 *
 * Either 1, 2, or 3 output files may be specified --- see below for more
 * details.  */
struct read_io_handle *
start_readers_and_writers(struct input_stream *in_1,
			  struct input_stream *in_2,
			  struct output_stream *out_combined,
			  struct output_stream *out_uncombined_1,
			  struct output_stream *out_uncombined_2,
			  const struct read_format_params *iparams,
			  const struct read_format_params *oparams,
			  unsigned num_combiner_threads,
			  bool verbose)
{
	assert(in_1 != NULL);
	assert(out_combined != NULL &&
	       (out_uncombined_1 != NULL || out_uncombined_2 == NULL));
	assert(iparams != NULL);
	assert(oparams != NULL);
	assert(num_combiner_threads > 0);

	if (verbose)
		info("Starting reader and writer threads");

	struct read_io_handle *h = xzalloc(sizeof(*h));

	size_t reads_per_set = BASE_READS_PER_READ_SET +
				(num_combiner_threads * PERTHREAD_READS_PER_READ_SET);
	size_t queue_size = num_combiner_threads * QUEUE_SIZE_PER_THREAD;

	h->avail_read_q = new_read_queue(queue_size * 3, reads_per_set, true);
	h->unprocessed_read_1_q = new_read_queue(queue_size, reads_per_set, false);
	h->unprocessed_read_2_q = new_read_queue(queue_size, reads_per_set, false);
	h->combined_read_q = new_read_queue(queue_size, reads_per_set, false);

	init_mutex(&h->get_unprocessed_pair_mutex);
	init_mutex(&h->put_uncombined_pair_mutex);

	h->combiner_threads_remaining = num_combiner_threads;
	init_mutex(&h->combiner_threads_remaining_mutex);

	/* Start writers.  */

	if (out_uncombined_2) {
		/* All 3 output files specified: one for combined reads, one for
		 * read 1 of uncombined pairs, and one for read 2 of uncombined
		 * pairs.  */

		h->uncombined_read_1_q = new_read_queue(queue_size, reads_per_set, false);
		h->uncombined_read_2_q = new_read_queue(queue_size, reads_per_set, false);

		h->writer_1 = start_writer1(out_combined, oparams,
					    h->combined_read_q,
					    h->avail_read_q);
		h->writer_1_started = true;

		h->writer_2 = start_writer1(out_uncombined_1, oparams,
					    h->uncombined_read_1_q,
					    h->avail_read_q);
		h->writer_2_started = true;

		h->writer_3 = start_writer1(out_uncombined_2, oparams,
					    h->uncombined_read_2_q,
					    h->avail_read_q);
		h->writer_3_started = true;
	} else if (out_uncombined_1) {
		/* 2 output files specified: one for combined reads and one for
		 * uncombined pairs.  */

		h->uncombined_read_1_q = new_read_queue(queue_size, reads_per_set, false);
		h->uncombined_read_2_q = new_read_queue(queue_size, reads_per_set, false);

		h->writer_1 = start_writer1(out_combined, oparams,
					    h->combined_read_q,
					    h->avail_read_q);
		h->writer_1_started = true;

		h->writer_2 = start_writer2(out_uncombined_1, oparams,
					    h->uncombined_read_1_q,
					    h->uncombined_read_2_q,
					    h->avail_read_q);
		h->writer_2_started = true;
	} else {
		/* 1 output file specified: combined reads, plus optionally
		 * uncombined pairs if supported by the format.  */

		if (read_format_supports_mixed_reads(oparams)) {
			h->uncombined_read_1_q = h->combined_read_q;
			h->uncombined_read_2_q = new_read_queue(queue_size, reads_per_set, false);

			h->writer_1 = start_writer2(out_combined, oparams,
						    h->combined_read_q,
						    h->uncombined_read_2_q,
						    h->avail_read_q);
			h->writer_1_started = true;
		} else {
			/* Can only output combined reads.
			 * Reroute uncombined reads back to the queue of
			 * available (for reuse) reads.  */
			h->uncombined_read_1_q = h->avail_read_q;
			h->uncombined_read_2_q = h->avail_read_q;

			h->writer_1 = start_writer1(out_combined, oparams,
						    h->combined_read_q,
						    h->avail_read_q);
			h->writer_1_started = true;
		}
	}

	/* Start readers.  */

	if (in_2) {
		/* Two input files:  read 1 in each pair comes from the first
		 * file, and read 2 in each pair comes from the second file.
		 *
		 * Only set @verbose for one.  */
		h->reader_1 = start_reader1(in_1,
					    iparams,
					    verbose,
					    h->avail_read_q,
					    h->unprocessed_read_1_q);
		h->reader_1_started = true;

		h->reader_2 = start_reader1(in_2,
					    iparams,
					    false,
					    h->avail_read_q,
					    h->unprocessed_read_2_q);
		h->reader_2_started = true;
	} else {
		/* One input file:  both reads in each pair come from the same
		 * file.  */
		struct read_queue *unpaired_read_q = NULL;

		if (read_format_supports_mixed_reads(iparams)) {
			if (!out_uncombined_2 &&
			    read_format_supports_mixed_reads(oparams))
				unpaired_read_q = h->uncombined_read_1_q;
			else
				warning("Any unpaired reads in the input file "
					"will be ignored!\n\t"
					"Use tab-delimited output to "
					"preserve them.");
		}

		h->reader_1 = start_reader2(in_1,
					    iparams,
					    verbose,
					    h->avail_read_q,
					    h->unprocessed_read_1_q,
					    h->unprocessed_read_2_q,
					    unpaired_read_q);
		h->reader_1_started = true;
	}


	return h;
}

/* Terminates the FLASH I/O layer, which is responsible for input/output of
 * reads.
 */
void
stop_readers_and_writers(struct read_io_handle *h)
{
	if (h->reader_1_started)
		join_thread(h->reader_1);
	if (h->reader_2_started)
		join_thread(h->reader_2);
	if (h->writer_1_started)
		join_thread(h->writer_1);
	if (h->writer_2_started)
		join_thread(h->writer_2);
	if (h->writer_3_started)
		join_thread(h->writer_3);

	free_read_queue(h->avail_read_q);
	free_read_queue(h->unprocessed_read_1_q);
	free_read_queue(h->unprocessed_read_2_q);
	free_read_queue(h->combined_read_q);

	if (h->uncombined_read_1_q != h->avail_read_q &&
	    h->uncombined_read_1_q != h->combined_read_q)
		free_read_queue(h->uncombined_read_1_q);

	if (h->uncombined_read_2_q != h->avail_read_q)
		free_read_queue(h->uncombined_read_2_q);

	pthread_mutex_destroy(&h->put_uncombined_pair_mutex);
	pthread_mutex_destroy(&h->get_unprocessed_pair_mutex);
	pthread_mutex_destroy(&h->combiner_threads_remaining_mutex);

	xfree(h, sizeof(*h));
}
