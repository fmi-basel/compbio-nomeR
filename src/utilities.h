#ifndef _utilities_h_
#define _utilities_h_

#include <stdio.h>
#include <math.h>
#include <Rinternals.h>
#include <samtools-1.7-compat.h>
#include <stdbool.h>

/*! @function
 @abstract  Get whether the query is fist mate
 @param  b  pointer to an alignment
 @return    boolean true if query is on the reverse strand
 */
#define bam_is_first(b) (((b)->core.flag&BAM_FREAD1) != 0)

/*! @function
 @abstract  Get whether the query is second mate
 @param  b  pointer to an alignment
 @return    boolean true if query is on the reverse strand
 */
#define bam_is_second(b) (((b)->core.flag&BAM_FREAD2) != 0)


/*! @function
 * @abstrat Get whether the query is not primary alignment
 * @param b pointer to an alignment
 * @return boolean true if query is not primary
 */

#define bam_is_not_primary(b) (((b)->core.flag&BAM_FSECONDARY) != 0)

/*! @function
 * @abstrat Get whether the query failed QC
 * @param b pointer to an alignment
 * @return boolean true if query failed QC
 */

#define bam_is_qc_fail(b) (((b)->core.flag&BAM_FQCFAIL) != 0)


/*! @function
 * @abstrat Get whether the query is optical or PCR duplicate
 * @param b pointer to an alignment
 * @return boolean true if query is optical or PCR duplicate
 */

#define bam_is_pcr_dupl(b) (((b)->core.flag&BAM_FDUP) != 0)

/*! @function
 * @abstrat Get whether the query is supplementary alignment
 * @param b pointer to an alignment
 * @return boolean true if query is supplementary alignment
 */

#define bam_is_suppl(b) (((b)->core.flag&BAM_FSUPPLEMENTARY) != 0)



samfile_t * _bam_tryopen(const char *filename, const char *filemode, void *aux);

int _fromQueryToRefPos(int aln_start, // most left position of alignment, 0-based
                       int query_pos, // position in a query sequence, 0-based
                       int n_cigar, // length of cigar array
                       const uint32_t *cigar // pointer to cigar array
);


int _fromRefToQueryPos(int aln_start, // most left position of alignment, 0-based
                      int ref_pos, // position in reference, 0-based
                      int n_cigar, // length of cigar array
                      const uint32_t *cigar // pointer to cigar array
                      
                      );


#endif