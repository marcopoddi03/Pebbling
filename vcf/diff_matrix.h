#ifndef DIFF_MATRIX_H
#define DIFF_MATRIX_H

#include <stdint.h>
#include <stdbool.h>
#include <htslib/synced_bcf_reader.h>

typedef enum {
    GT_MATCH,
    GT_MISMATCH,
    GT_UNKNOWN 
} GtComparisonResult;

typedef struct {
    int *allele_counts;
    int ploidy;
    int missing_count;
} GenotypeInfo;

typedef struct {
    int **matrix;
    int n_samples;
    char **samples;
} DiffMatrix;

GtComparisonResult gt_cmp(const GenotypeInfo *gt_i, const GenotypeInfo *gt_j, int n_allele);
bool save_diff_matrix_bin(const char *filename, DiffMatrix *diff_matrix);
bool save_diff_matrix_txt(const char *filename, DiffMatrix *diff_matrix);
DiffMatrix *init_diff_matrix(const bcf_hdr_t *hdr);
char *format_regions(char *list[], int start_idx, int n);
void free_diff_matrix(DiffMatrix *diff_matrix);
bcf_srs_t *set_reader_for_diff(const char *filename, const char *regions);
void update_diff_matrix(DiffMatrix *diff_matrix, int i, int j);
DiffMatrix *get_bcf_diff_matrix(const char *filename, char *regions);
void free_genotype_info(GenotypeInfo *gt_info);
GenotypeInfo *get_genotype_info(const int32_t *gt, int max_ploidy, int n_allele);

#endif 
