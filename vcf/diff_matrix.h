#ifndef DIFF_MATRIX_H
#define DIFF_MATRIX_H

#include <stdint.h>
#include <stdbool.h>
#include <htslib/synced_bcf_reader.h>

typedef enum {
    GT_MATCH,
    GT_MISMATCH,
    GT_UNKNOWN,
    GT_DIFF_PLOIDY 
} GtComparisonResult;

typedef struct {
    int *allele_counts;
    int ploidy;
    int missing_count;
} GenotypeInfo;

typedef struct {
    int n_samples;
    char **samples;
} PairwiseMatrixInfo;

typedef struct {
    int **matrix;
    PairwiseMatrixInfo pm_info;
} PairwiseMatrix;

// Main function
PairwiseMatrix *get_bcf_pairwise_matrix(const char *filename, char *regions);

// Genotype functions
GtComparisonResult gt_cmp(const GenotypeInfo *gt_i, const GenotypeInfo *gt_j, int n_allele);
GenotypeInfo *get_genotype_info(const int32_t *gt, int max_ploidy, int n_allele);
void free_genotype_info(GenotypeInfo *gt_info);

// Matrix functions
PairwiseMatrix *init_pairwise_matrix(const bcf_hdr_t *hdr);
void free_pairwise_matrix(PairwiseMatrix *pairwise_matrix);
bool update_pairwise_matrix(PairwiseMatrix *pairwise_matrix, int i, int j);
int save_pairwise_matrix_bin(const char *filename, PairwiseMatrix *pairwise_matrix);
int save_pairwise_matrix_txt(const char *filename, PairwiseMatrix *pairwise_matrix);

// Others
char *format_regions(char *list[], int start_idx, int n);
bcf_srs_t *set_reader_for_pairwise(const char *filename, const char *regions);

#endif 
