#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <htslib/synced_bcf_reader.h>
#include <string.h>
#include "diff_matrix.h"

// Main function
DiffMatrix *get_bcf_diff_matrix(const char *filename, char *regions) {

    bcf_srs_t *sr = set_reader_for_diff(filename, regions);
    if (!sr)
        return NULL;

    bcf_hdr_t *hdr = bcf_sr_get_header(sr, 0);
    if (!hdr) {
        bcf_sr_destroy(sr);
        return NULL;
    }
        

    DiffMatrix *diff_matrix = init_diff_matrix(hdr);
    if(!diff_matrix) {
        bcf_sr_destroy(sr);
        return NULL;
    }

    int n_samples = bcf_hdr_nsamples(hdr);
    if (n_samples <= 0)
    {
        bcf_sr_destroy(sr);
        return NULL;
    } 

    int32_t *gt_arr = NULL;
    int ngt_arr = 0;

    GenotypeInfo **gt_info = malloc((size_t)n_samples * sizeof(*gt_info));
    if (!gt_info) {
        free_diff_matrix(diff_matrix);
        bcf_sr_destroy(sr);
        return NULL;
    }

    while (bcf_sr_next_line(sr)) {
        bcf1_t *rec = bcf_sr_get_line(sr, 0);
        if (!rec)
            continue;

        int ngt = bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr);
        if (ngt <= 0) {
            continue;
        }

        if (ngt % n_samples != 0) {
            free(gt_arr);
            free_diff_matrix(diff_matrix);
            bcf_sr_destroy(sr);
            free(gt_info);
            return NULL;
        }

        int max_ploidy = ngt / n_samples;
        int n_allele = rec->n_allele;

        for (int k = 0; k < n_samples; k++) {
            int32_t *gt_k = gt_arr + k * max_ploidy;
            gt_info[k] = get_genotype_info(gt_k, max_ploidy, n_allele);

            if (!gt_info[k]) {
                for (int j = 0; j < k; j++) {
                    free_genotype_info(gt_info[j]);
                }
                free_diff_matrix(diff_matrix);
                bcf_sr_destroy(sr);
                free(gt_info);
                free(gt_arr);
                return NULL;
            }
        }

        for (int i = 0; i < n_samples; i++) {
            for (int j = i + 1; j < n_samples; j++) {
                if(gt_cmp(gt_info[i], gt_info[j], n_allele) == GT_MISMATCH) 
                    update_diff_matrix(diff_matrix, i, j);
            }

        }  

        for (int k=0; k < n_samples; k++) {
            free_genotype_info(gt_info[k]);
            gt_info[k] = NULL;
        }

    }

    if (sr->errnum != 0) {
        fprintf(stderr, "Error reading file: %s\n", bcf_sr_strerror(sr->errnum));

        free(gt_arr);
        free_diff_matrix(diff_matrix);
        bcf_sr_destroy(sr);
        free(gt_info);
        return NULL;
    }

    bcf_sr_destroy(sr);
    free(gt_arr); 
    free(gt_info);

    return diff_matrix;
}



// Genotype functions

GtComparisonResult gt_cmp(const GenotypeInfo *gt_i, const GenotypeInfo *gt_j, int n_allele) {

    if (!gt_i || !gt_j || n_allele <= 0)
        return GT_UNKNOWN;

    if (gt_i->ploidy != gt_j->ploidy)
        return GT_MISMATCH;

    int required = 0;

    for (int k = 0; k < n_allele; k++) {
        int a1 = gt_i->allele_counts[k];
        int a2 = gt_j->allele_counts[k];

        if (a1 < 0 || a2 < 0)
            return GT_UNKNOWN;

        required += (a1 > a2) ? a1 : a2;
    }

    if (required > gt_i->ploidy)
        return GT_MISMATCH;

    if (gt_i->missing_count > 0 || gt_j->missing_count > 0) {
        return GT_UNKNOWN;
    }

    return GT_MATCH; 
}

GenotypeInfo *get_genotype_info(const int32_t *gt, int max_ploidy, int n_allele) {
    if (!gt || max_ploidy <= 0 || n_allele <= 0)
        return NULL;

    GenotypeInfo *gt_info = malloc(sizeof(*gt_info));
    if (!gt_info)
        return NULL;

    int *allele_counts = calloc((size_t)n_allele, sizeof(int));

    if (!allele_counts) {
        free(gt_info);
        return NULL;
    }

    int missing_count = 0;
    int ploidy = 0;

    for (int k = 0; k < max_ploidy; k++) {

        if (gt[k] == bcf_int32_vector_end)
            break;

        ploidy++;

        if (bcf_gt_is_missing(gt[k])) {
            missing_count++;
            continue;
        }

        int allele = bcf_gt_allele(gt[k]);

        if (allele < 0 || allele >= n_allele) {
            free(allele_counts);
            free(gt_info);
            return NULL;
        }

        allele_counts[allele]++;
    }

    gt_info->allele_counts = allele_counts;
    gt_info->missing_count = missing_count;
    gt_info->ploidy = ploidy;

    return gt_info;
}

void free_genotype_info(GenotypeInfo *gt_info) {
    if (!gt_info)
        return;
    free(gt_info->allele_counts);
    free(gt_info);
}



// Matrix functions

DiffMatrix *init_diff_matrix(const bcf_hdr_t *hdr) {
    if (!hdr || !hdr->samples)
        return NULL;

    int n_samples = bcf_hdr_nsamples(hdr);
    if (n_samples <=0)
        return NULL;

    int **matrix = malloc(n_samples * sizeof(int *));
    if (!matrix) {
        fprintf(stderr, "Error: memory allocation failed\n");
        return NULL;
    }

    size_t n = (size_t)n_samples;

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

    char **samples = malloc(n_samples * sizeof(char *));
    if(!samples) {
        fprintf(stderr, "Error: memory allocation failed\n");
        free(data);
        free(matrix);
        free(diff_matrix);
        return NULL;
    }

    for (int i = 0; i < n_samples; i++) {
        matrix[i] = data + (i * n_samples);

        samples[i] = strdup(hdr->samples[i]);
        if (!samples[i]) {
            fprintf(stderr, "Error: memory allocation failed\n");
            for (int j = 0; j < i; j++)
                free(samples[j]);
            free(data);
            free(matrix);
            free(diff_matrix);
            free(samples);
            return NULL;
        }
    }

    diff_matrix->matrix = matrix;
    diff_matrix->n_samples = n_samples;
    diff_matrix->samples = samples;

    return diff_matrix;
}

void free_diff_matrix(DiffMatrix *diff_matrix) {
    
    if (!diff_matrix)
     return;

    if (diff_matrix->matrix) {
        free(diff_matrix->matrix[0]);
        free(diff_matrix->matrix);
    }

    if (diff_matrix->samples) {
        for (int i = 0; i < diff_matrix->n_samples; i++)
            free(diff_matrix->samples[i]);
        free(diff_matrix->samples);
    }

    free(diff_matrix);
}

void update_diff_matrix(DiffMatrix *diff_matrix, int i, int j) {
    diff_matrix->matrix[i][j]++;
    diff_matrix->matrix[j][i]++;
}

bool save_diff_matrix_bin(const char *filename, DiffMatrix *diff_matrix) {

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

bool save_diff_matrix_txt(const char *filename, DiffMatrix *diff_matrix) {

    FILE *file = fopen(filename, "w");
    if (!file) {
        fprintf(stderr, "Error: could not open file\n");
        return false;
    }

    int n = diff_matrix->n_samples;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            fprintf(file, "%d", diff_matrix->matrix[i][j]);

            if (j < n - 1) {
                fprintf(file, " ");
            }
        }
        fprintf(file, "\n");
    }

    fclose(file);

    return true;
}



// Others

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