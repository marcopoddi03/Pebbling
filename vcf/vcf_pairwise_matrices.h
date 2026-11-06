#ifndef VCF_PAIRWISE_MATRICES_H
#define VCF_PAIRWISE_MATRICES_H

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
    int **diss_matrix;
    int **sim_matrix;
    PairwiseMatrixInfo pm_info;
} PairwiseMatrices;

// Main function
PairwiseMatrices *get_vcf_pairwise_matrices(const char *filename, char *regions);

// Genotype functions
GtComparisonResult gt_cmp(const GenotypeInfo *gt_i, const GenotypeInfo *gt_j, int n_allele);
GenotypeInfo *get_genotype_info(const int32_t *gt, int max_ploidy, int n_allele);
void free_genotype_info(GenotypeInfo *gt_info);

// Matrix functions
PairwiseMatrices *init_pairwise_matrices(const bcf_hdr_t *hdr);
void free_pairwise_matrices(PairwiseMatrices *pairwise_matrices);
bool update_pairwise_matrices(PairwiseMatrices *pairwise_matrices, int i, int j);
int save_pairwise_matrices_bin(const char *filename, PairwiseMatrices *pairwise_matrices);
int save_pairwise_matrices_txt(const char *filename, PairwiseMatrices *pairwise_matrices);

// Others
char *format_regions(char *list[], int start_idx, int n);
bcf_srs_t *set_reader_for_pairwise(const char *filename, const char *regions);

#endif 
