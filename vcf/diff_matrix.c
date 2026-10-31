#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <htslib/hts.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/vcf.h>
#include "diff_matrix.h"

int int32_t_cmp(const void *a, const void *b) {
    int32_t x = *((const int32_t *)a);
    int32_t y = *((const int32_t *)b);
    return (x > y) - (x < y);
}

bool save_diff_matrix(const char *filename, DiffMatrix *diff_matrix) {

    FILE *file = fopen(filename, "wb");
    if (!file) {
        fprintf(stderr, "Error: could not open file\n");
        return false;
    }

    size_t dim = (size_t)diff_matrix->n_samples * diff_matrix->n_samples;
    size_t writed = fwrite(diff_matrix->matrix[0], sizeof(int), dim, file);

    fclose(file);

    return (writed == dim);
}

GtComparisonResult gt_cmp(const int32_t *gt_i, const int32_t *gt_j, int max_ploidy) {

    if (!gt_i || !gt_j || max_ploidy <= 0)
        return GT_UNKNOWN;
    int32_t tmp_i[max_ploidy];
    int32_t tmp_j[max_ploidy];

    for (int k = 0; k < max_ploidy; k++) {

        if (gt_i[k] == bcf_int32_vector_end)
            tmp_i[k] = -1;
        else if (bcf_gt_is_missing(gt_i[k]))
            return GT_UNKNOWN;
        else
            tmp_i[k] = bcf_gt_allele(gt_i[k]);

        if (gt_j[k] == bcf_int32_vector_end)
            tmp_j[k] = -1;
        else if (bcf_gt_is_missing(gt_j[k]))
            return GT_UNKNOWN;
        else
            tmp_j[k] = bcf_gt_allele(gt_j[k]);

    }

    qsort(tmp_i, max_ploidy, sizeof(int32_t), int32_t_cmp);
    qsort(tmp_j, max_ploidy, sizeof(int32_t), int32_t_cmp); 

    for (int k = 0; k < max_ploidy; k++) {
        if (tmp_i[k] != tmp_j[k])
            return GT_MISMATCH;
    }
    
    return GT_MATCH; 
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

DiffMatrix *init_diff_matrix(int n_samples) {
    int **matrix = malloc(n_samples * sizeof(int *));
    if (!matrix) {
        fprintf(stderr, "Error: memory allocation failed\n");
        return NULL;
    }

    if (n_samples <= 0)
        return NULL;


    int n = (size_t)n_samples;

    int *data = calloc(n * n, sizeof(int));
    if (!data) {
        fprintf(stderr, "Error: memory allocation failed\n");
        free(matrix);
        return NULL;
    }

    DiffMatrix *diff_matrix = malloc(sizeof(DiffMatrix));
    if (!diff_matrix) {
        fprintf(stderr, "Error: memory allocation failed\n");
        free(data);
        free(matrix);
        return NULL;
    }

    for (int i = 0; i < n_samples; i++) {
        matrix[i] = data + (i * n_samples);
    }

    diff_matrix->matrix = matrix;
    diff_matrix->n_samples = n_samples;

    return diff_matrix;
}

void free_diff_matrix(DiffMatrix *diff_matrix) {
    if (diff_matrix) {
        free(diff_matrix->matrix[0]);
        free(diff_matrix->matrix);
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


void update_diff_matrix(DiffMatrix *diff_matrix, int i, int j) {
    diff_matrix->matrix[i][j]++;
    diff_matrix->matrix[j][i]++;
}

DiffMatrix *get_bcf_diff_matrix(const char *filename, char *regions) {

    bcf_srs_t *sr = set_reader_for_diff(filename, regions);
    if (!sr)
        return NULL;

    bcf_hdr_t *hdr = bcf_sr_get_header(sr, 0);
    if (!hdr) {
        bcf_sr_destroy(sr);
        return NULL;
    }
        
    int n_samples = bcf_hdr_nsamples(hdr);
    if (n_samples <= 0)
        return NULL;    

    DiffMatrix *diff_matrix = init_diff_matrix(n_samples);
    if(!diff_matrix) {
        bcf_sr_destroy(sr);
        return NULL;
    }

    int32_t *gt_arr = NULL;
    int ngt_arr = 0;

    while (bcf_sr_next_line(sr)) {
        bcf1_t *rec = bcf_sr_get_line(sr, 0);
        if (!rec)
            continue;


        int ngt = bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr);
        if (ngt <= 0) {
            free(gt_arr);
            continue;
        }

        if (ngt % n_samples != 0) {
            return NULL;
        }

        int max_ploidy = ngt / n_samples;

        for (int i = 0; i < n_samples; i++) {
            int32_t *gt_i = gt_arr + i * max_ploidy;

            for (int j = i + 1; j < n_samples; j++) {

                int32_t *gt_j = gt_arr + j * max_ploidy;

                if(gt_cmp(gt_i, gt_j, max_ploidy) == GT_MISMATCH) 
                    update_diff_matrix(diff_matrix, i, j);
            }
        }  
    }

    if (sr->errnum != 0) {
        fprintf(stderr, "Error reading file: %s\n", bcf_sr_strerror(sr->errnum));

        free(gt_arr);
        free_diff_matrix(diff_matrix);
        bcf_sr_destroy(sr);

        return NULL;
    }

    bcf_sr_destroy(sr);
    free(gt_arr); 

    return diff_matrix;
}