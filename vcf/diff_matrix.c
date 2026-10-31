#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <htslib/hts.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/vcf.h>


int int_cmp(const void *a, const void *b) {
    int int_a = *((const int *)a);
    int int_b = *((const int *)b);
    return (int_a > int_b) - (int_a < int_b);
}

bool save_diff_matrix(const char *filename, int **matrix, int n) {
    FILE *file = fopen(filename, "wb");
    if (!file) {
        fprintf(stderr, "Error: could not open file\n");
        return false;
    }

    size_t dim = (size_t)n * n;
    size_t writed = fwrite(matrix[0], sizeof(int), dim, file);

    fclose(file);

    return (writed == dim);
}

bool gt_cmp(const int32_t *gt_i, const int32_t *gt_j, int max_ploidy) {
    if (!gt_i || !gt_j || max_ploidy <= 0)
        return false;
    int32_t tmp_i[max_ploidy];
    int32_t tmp_j[max_ploidy];
    for (int k = 0; k < max_ploidy; k++) {

        if (gt_i[k] == bcf_int32_vector_end)
            tmp_i[k] = -1;
        else if (bcf_gt_is_missing(gt_i[k]))
            return false;
        else
            tmp_i[k] = bcf_gt_allele(gt_i[k]);

        if (gt_j[k] == bcf_int32_vector_end)
            tmp_j[k] = -1;
        else if (bcf_gt_is_missing(gt_j[k]))
            return false;
        else
            tmp_j[k] = bcf_gt_allele(gt_j[k]);

    }

    qsort(tmp_i, max_ploidy, sizeof(int32_t), int_cmp);
    qsort(tmp_j, max_ploidy, sizeof(int32_t), int_cmp); 

    for (int k = 0; k < max_ploidy; k++) {
        if (tmp_i[k] != tmp_j[k])
            return false;
    }
    
    return true; 
}


char *format_regions(char *list[], int start_idx, int n) {
    if (start_idx >= n) {
        char *empty = malloc(1);
        if (empty) 
            empty[0] = '\0';
        return empty;
    }
    size_t len = 1; 
    for (int i = start_idx; i < n; i++)
        len += strlen(list[i]) + (i > start_idx ? 1 : 0);
    char *regions = malloc(len);
    if (!regions) {
        fprintf(stderr, "Error: memory allocation failed\n");
        return NULL;
    }

    regions[0] = '\0';

    for (int i = start_idx; i < n; i++) {
        if (i > start_idx)
            strcat(regions, ",");

        strcat(regions, list[i]);
    }

    return regions;  
}




int **init_diff_matrix(int n_samples) {
    int **diff_matrix = malloc(n_samples * sizeof(int *));
    if (!diff_matrix) {
        fprintf(stderr, "Error: memory allocation failed\n");
        return NULL;
    }

    int *data = calloc(n_samples * n_samples, sizeof(int));
    if (!data) {
        fprintf(stderr, "Error: memory allocation failed\n");
        free(diff_matrix);
        return NULL;
    }

    for (int i = 0; i < n_samples; i++) {
        diff_matrix[i] = data + (i * n_samples);
    }

    return diff_matrix;
}

void free_diff_matrix(int **diff_matrix) {
    if (diff_matrix) {
        free(diff_matrix[0]);
        free(diff_matrix);
    }
}


bcf_srs_t *set_reader_for_diff(const char *filename, const char *regions) {
    bcf_srs_t *sr = bcf_sr_init();
    if (!sr) {
        fprintf(stderr, "Error: could not inizialize reader\n");
        return NULL;
    }

    if (bcf_sr_set_opt(sr, BCF_SR_REQUIRE_IDX)) {
        fprintf(stderr, "Error: could not recognise option\n");
        bcf_sr_destroy(sr);
        return NULL;
    }

    if (!bcf_sr_add_reader(sr, filename)) {
        fprintf(stderr, "Error: %s\n", bcf_sr_strerror(sr->errnum));
        bcf_sr_destroy(sr);
        return NULL;
    }

    if (bcf_sr_set_regions(sr, regions, 0) < 0) {
        fprintf(stderr, "Error: could not set regions\n");
        bcf_sr_destroy(sr);
        return NULL;
    }

    return sr;
}


void update_diff_matrix(int **diff_matrix, int i, int j) {
    //da implementare
}