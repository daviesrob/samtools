/*
    tmp_file.h - write to and read from a temporary binary file
    for fast storage plus added compression.
    
    Andrew Whitwham, August 2017
*/

#ifndef _TMP_SAM_FILE_H_
#define _TMP_SAM_FILE_H_

#include <lz4.h>
#include "htslib/sam.h"

#ifdef _cplusplus
extern "C" {
#endif

#define TMP_SAM_GROUP_SIZE 100
#define TMP_SAM_MAX_DATA 1024
#define TMP_SAM_RING_SIZE 1048576
#define TMP_SAM_OK 0
#define TMP_SAM_MEM_ERROR -1
#define TMP_SAM_FILE_ERROR -2
#define TMP_SAM_LZ4_ERROR -3
#define TMP_SAM_INPUT_ERROR -4 

typedef struct {
    FILE *fp;
    LZ4_stream_t *stream;
    LZ4_streamDecode_t *dstream;
    size_t data_size;
    size_t max_data_size;
    size_t ring_buffer_size;
    size_t comp_buffer_size;
    size_t offset;
    uint8_t *data;
    uint8_t *ring_buffer;
    uint8_t *ring_index;
    char *comp_buffer;
    char *name;
    size_t group_size;
    size_t input_size;
    size_t read_size;
    size_t output_size;
    size_t entry_number;
    int verbose;
    char *dict;
    size_t groups_written;
} tmp_file_t;


/*
 * Opens the temp file and initialises memory.
 * Verbose mode prints out error messages to stderr.
 * Returns 0 on success, a negative number on failure.
 */
int tmp_file_open_write(tmp_file_t *tmp, char *tmp_name, int verbose);


/* 
 * Stores an in memory bam structure for writing and if enough are gathered together writes
 * it to disk.  Mulitiple alignments compress better that single ones though after a certain number
 * there is a law of diminishing returns.
 * Returns 0 on success, a negative number on failure.
 */
int tmp_file_write(tmp_file_t *tmp, bam1_t *inbam);


/*
 * Closes the file after writing out any remaining alignments.  Adds a size_t 0 to
 * mark the end of the file.  Companion function to tmp_file_open_read below.
 * Returns 0 on success, a negative number on failure.
 */
int tmp_file_close_write(tmp_file_t *tmp);

    
/*
 * Opens the file for reading.  Optionally, if given a pointer to an existing
 * bam1_t structure, it will free the data entry to prevent memory leaks.
 * Companion function to tmp_file_close_write above.
 * Returns 0 on success, a negative number on failure.
 */
int tmp_file_open_read(tmp_file_t *tmp, bam1_t *inbam);


/*
 * An alternative to tmp_file_close_write that does the same job without actually
 * closing the file. Companion function to tmp_file_begin_read below.
 * Returns 0 on success, a negative number on failure.
 */
int tmp_file_end_write(tmp_file_t *tmp);
    

/*
 * An alternative to tmp_file_open_read but works on an open file.
 * Companion function to tmp_file_end_write above.
 * Returns 0 on success, a negative number on failure.
 */
int tmp_file_begin_read(tmp_file_t *tmp, bam1_t *inbam);
    

/*
 * Read the next alignment, either from memory or from disk.
 * Returns size of entry on success, 0 on end of file or a negative on error.
 */
int tmp_file_read(tmp_file_t *tmp, bam1_t *inbam);


/*
 * Frees up memory, closes the file and optionally deletes it.  Giving this function
 * pointer to the bam1_t structure used for reading will set its data value to null,
 * preventing bam_destroy1() from trying to free already freed memory.
 * Returns 0 on success, a negative number or EOF on failure.
 */
int tmp_file_destroy(tmp_file_t *tmp, bam1_t *inbam, int delete);

#ifdef __cplusplus
}
#endif

#endif /* _TMP_SAM_FILE_H_ */
