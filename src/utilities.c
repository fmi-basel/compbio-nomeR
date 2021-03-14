#include "utilities.h"

// open bamfile. modified code from QuasR

samfile_t * _bam_tryopen(const char *filename, const char *filemode, void *aux)
{
	samfile_t *sfile = samopen(filename, filemode, aux);
	if (sfile == 0)
		Rf_error("failed to open SAM/BAM file\n  file: '%s'", 
           filename);
	if (sfile->header == 0 || sfile->header->n_targets == 0) {
		samclose(sfile);
		Rf_error("SAM/BAM header missing or empty\n  file: '%s'", 
           filename);
	}
	return sfile;
}


/* convert position in query into position in reference
 * return position in reference or -1 if not found
 */
int _fromQueryToRefPos(int aln_start, // most left position of alignment, 0-based
                       int query_pos, // position in a query sequence, 0-based
                       int n_cigar, // length of cigar array
                       const uint32_t *cigar // pointer to cigar array
){
	int qpos = 0, rpos = aln_start;
	bool found = 0;
	for (int k = 0; k < n_cigar; ++k) {
		int type = bam_cigar_type(bam_cigar_op(cigar[k]));
		int len = bam_cigar_oplen(cigar[k]);
		// cigar is MATCH
		if(bam_cigar_op(cigar[k]) == BAM_CMATCH || bam_cigar_op(cigar[k]) == BAM_CEQUAL || bam_cigar_op(cigar[k]) == BAM_CDIFF){
			// go across cigar len and search for query_pos
			for(int p = 0; p < len; ++p){
				if (type & 1) qpos++;
				if (type & 2) rpos++;
				if(qpos == query_pos){
					found = 1;
					break;
				}
			}
		} else{
			if (type & 1) qpos += len;
			if (type & 2) rpos += len;
		}
		
		if(found)
			break;
		
	}
	
	if(found)
		return(rpos);
	else
		return(-1);
}


/* convert position in reference into position in a query
 * return position in query or -1 if not found
 */
int _fromRefToQueryPos(int aln_start, // most left position of alignment, 0-based
                       int ref_pos, // position in reference, 0-based
                       int n_cigar, // length of cigar array
                       const uint32_t *cigar // pointer to cigar array
                         
){
	if(ref_pos < aln_start)
		return(-1);
	int qpos = 0, rpos = aln_start;
	bool found = 0;
	
	for (int k = 0; k < n_cigar; ++k) {
		int type = bam_cigar_type(bam_cigar_op(cigar[k]));
		int len = bam_cigar_oplen(cigar[k]);
		// cigar is MATCH
		if(bam_cigar_op(cigar[k]) == BAM_CMATCH || bam_cigar_op(cigar[k]) == BAM_CEQUAL || bam_cigar_op(cigar[k]) == BAM_CDIFF){
			// go across cigar len and search for ref_pos
			for(int p = 0; p < len; ++p){
				if (type & 1) qpos++;
				if (type & 2) rpos++;
				if(rpos == ref_pos){
					found = 1;
					break;
				}
			}
		} else{
			if (type & 1) qpos += len;
			if (type & 2) rpos += len;
		}
		
		if(found)
			break;
		
	}
	
	if(found)
		return(qpos);
	else
		return(-1);
	
}


