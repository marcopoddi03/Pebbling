#ifndef DIFF_MATRIX_H
#define DIFF_MATRIX_H

bool gt_cmp(const int32_t *gt_i, const int32_t *gt_j, int max_ploidy);
int int_cmp(const void *a, const void *b);
bool save_diff_matrix(const char *filename, int **matrix, int n);
int **init_diff_matrix(int n_samples);
char *format_regions(char *list[], int start_idx, int n);
void free_diff_matrix(int **diff_matrix);
bcf_srs_t *set_reader_for_diff(const char *filename, const char *regions);
void update_diff_matrix(int **diff_matrix, int i, int j);

#endif 
