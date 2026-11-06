#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <htslib/synced_bcf_reader.h>
#include <string.h>
#include "diff_matrix.h"

// Main function
PairwiseMatrix *get_bcf_pairwise_matrix(const char *filename, char *regions) {

    bcf_srs_t *sr = set_reader_for_pairwise(filename, regions);
    if (!sr)
        return NULL;

    bcf_hdr_t *hdr = bcf_sr_get_header(sr, 0);
    if (!hdr) {
        bcf_sr_destroy(sr);
        return NULL;
    }
       
    int n_samples = bcf_hdr_nsamples(hdr);
    if (n_samples <= 0)
    {
        bcf_sr_destroy(sr);
        return NULL;
    } 

    PairwiseMatrix *pairwise_matrix = init_pairwise_matrix(hdr);
    if(!pairwise_matrix) {
        bcf_sr_destroy(sr);
        return NULL;
    }

    int32_t *gt_arr = NULL;
    int ngt_arr = 0;
    size_t n = (size_t)n_samples;

    GenotypeInfo **gt_info = malloc(n * sizeof(*gt_info));
    if (!gt_info) {
        free_pairwise_matrix(pairwise_matrix);
        bcf_sr_destroy(sr);
        return NULL;
    }

    while (bcf_sr_next_line(sr)) {
        bcf1_t *rec = bcf_sr_get_line(sr, 0);
        if (!rec) {
            free(gt_arr);
            free_pairwise_matrix(pairwise_matrix);
            bcf_sr_destroy(sr);
            free(gt_info);
            return NULL;
        }

        int ngt = bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr);
        if (ngt <= 0) {
            continue;
        }

        if (ngt % n != 0) {
            free(gt_arr);
            free_pairwise_matrix(pairwise_matrix);
            bcf_sr_destroy(sr);
            free(gt_info);
            return NULL;
        }

        int max_ploidy = ngt / n;
        int n_allele = rec->n_allele;

        for (int k = 0; k < n; k++) {
            int32_t *gt_k = gt_arr + k * max_ploidy;
            gt_info[k] = get_genotype_info(gt_k, max_ploidy, n_allele);

            if (!gt_info[k]) {
                for (int j = 0; j < k; j++) {
                    free_genotype_info(gt_info[j]);
                }
                free_pairwise_matrix(pairwise_matrix);
                bcf_sr_destroy(sr);
                free(gt_info);
                free(gt_arr);
                return NULL;
            }
        }

        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                if(gt_cmp(gt_info[i], gt_info[j], n_allele) == GT_MISMATCH) 
                    update_pairwise_matrix(pairwise_matrix, i, j);
            }

        }  

        for (int k=0; k < n; k++) {
            free_genotype_info(gt_info[k]);
            gt_info[k] = NULL;
        }

    }

    if (sr->errnum != 0) {
        free(gt_arr);
        free_pairwise_matrix(pairwise_matrix);
        bcf_sr_destroy(sr);
        free(gt_info);
        return NULL;
    }

    bcf_sr_destroy(sr);
    free(gt_arr); 
    free(gt_info);

    return pairwise_matrix;
}



// Genotype functions

GtComparisonResult gt_cmp(const GenotypeInfo *gt_i, const GenotypeInfo *gt_j, int n_allele) {

    if (!gt_i || !gt_j || n_allele <= 0 || !gt_i->allele_counts || !gt_j->allele_counts || gt_i->ploidy <= 0 || gt_j->ploidy <= 0)
        return GT_UNKNOWN;

    if (gt_i->ploidy != gt_j->ploidy)
        return GT_DIFF_PLOIDY;

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

PairwiseMatrix *init_pairwise_matrix(const bcf_hdr_t *hdr) {
    if (!hdr || !hdr->samples)
        return NULL;

    int n_samples = bcf_hdr_nsamples(hdr);
    if (n_samples <=0)
        return NULL;

    size_t n = (size_t)n_samples;

    int **matrix = malloc(n * sizeof(int *));
    if (!matrix) {
        return NULL;
    }

    int *data = calloc(n * n, sizeof(int));
    if (!data) {
        free(matrix);
        return NULL;
    }

    PairwiseMatrix *pairwise_matrix = malloc(sizeof(PairwiseMatrix));
    if (!pairwise_matrix) {
        free(data);
        free(matrix);
        return NULL;
    }

    char **samples = malloc(n * sizeof(char *));
    if(!samples) {
        free(data);
        free(matrix);
        free(pairwise_matrix);
        return NULL;
    }

    for (int i = 0; i < n; i++) {
        matrix[i] = data + (i * n);

        samples[i] = strdup(hdr->samples[i]);
        if (!samples[i]) {
            for (int j = 0; j < i; j++)
                free(samples[j]);
            free(data);
            free(matrix);
            free(pairwise_matrix);
            free(samples);
            return NULL;
        }
    }

    pairwise_matrix->matrix = matrix;
    pairwise_matrix->pm_info.n_samples = n_samples;
    pairwise_matrix->pm_info.samples = samples;

    return pairwise_matrix;
}

void free_pairwise_matrix(PairwiseMatrix *pairwise_matrix) {
    
    if (!pairwise_matrix)
        return;

    if (pairwise_matrix->matrix) {
        free(pairwise_matrix->matrix[0]);
        free(pairwise_matrix->matrix);
    }

    if (pairwise_matrix->pm_info.samples) {
        for (size_t i = 0; i < pairwise_matrix->pm_info.n_samples; i++)
            free(pairwise_matrix->pm_info.samples[i]);
        free(pairwise_matrix->pm_info.samples);
    }

    free(pairwise_matrix);
}

bool update_pairwise_matrix(PairwiseMatrix *pairwise_matrix, int i, int j) {

    if (!pairwise_matrix || i < 0 || j < 0 || !pairwise_matrix->matrix || !pairwise_matrix->matrix[0])
        return false;

    size_t n = (size_t)pairwise_matrix->pm_info.n_samples;
    if (i >= n || j >= n)
        return false;

    if (i == j)
        return true;

    pairwise_matrix->matrix[i][j]++;
    pairwise_matrix->matrix[j][i]++;

    return true;
}

int save_pairwise_matrix_bin(const char *filename, PairwiseMatrix *pairwise_matrix) {

    if (!pairwise_matrix || !pairwise_matrix->matrix || !pairwise_matrix->matrix[0] || pairwise_matrix->pm_info.n_samples <= 0){
        return -1;
    }

    size_t n = (size_t)pairwise_matrix->pm_info.n_samples;

    FILE *file = fopen(filename, "wb");
    if (!file) {
        return -2;
    }

    size_t dim = n * n;
    size_t writed = fwrite(pairwise_matrix->matrix[0], sizeof(int), dim, file);

    fclose(file);

    if (writed != dim)
        return -2;

    return 0;  // -1 matrice non valida, -2 problema file, 0 successo
}

int save_pairwise_matrix_txt(const char *filename, PairwiseMatrix *pairwise_matrix) {

    if (!pairwise_matrix || !pairwise_matrix->matrix || !pairwise_matrix->matrix[0] || pairwise_matrix->pm_info.n_samples <= 0){
        return -1;
    }

    FILE *file = fopen(filename, "w");
    if (!file) {
        return -2;
    }

    int n = pairwise_matrix->pm_info.n_samples;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            fprintf(file, "%d", pairwise_matrix->matrix[i][j]);

            if (j < n - 1)
                fprintf(file, " ");
        }
        fprintf(file, "\n");
    }

    fclose(file);

    return 0; // -1 matrice non valida, -2 problema file, 0 successo
}



// Others

char *format_regions(char *list[], int start_idx, int n) { 
    if (!list || start_idx < 0 || n < 0 || start_idx >= n)
        return NULL;

    size_t len = 1; 
    for (int i = start_idx; i < n; i++) {
        if (!list[i])
            return NULL;
        len += strlen(list[i]) + (i > start_idx ? 1 : 0);
    }
        
    char *regions = malloc(len);
    if (!regions)
        return NULL;

    regions[0] = '\0';

    for (int i = start_idx; i < n; i++) {
        if (i > start_idx)
            strcat(regions, ",");

        strcat(regions, list[i]);
    }

    return regions;  
}

bcf_srs_t *set_reader_for_pairwise(const char *filename, const char *regions) {

    if(!regions)
        return NULL;

    bcf_srs_t *sr = bcf_sr_init();
    if (!sr) {
        return NULL;
    }

    if (bcf_sr_set_opt(sr, BCF_SR_REQUIRE_IDX)) {
        bcf_sr_destroy(sr);
        return NULL;
    }

    if (!bcf_sr_add_reader(sr, filename)) {
        bcf_sr_destroy(sr);
        return NULL;
    }

    if (bcf_sr_set_regions(sr, regions, 0) < 0) {
        bcf_sr_destroy(sr);
        return NULL;
    }

    return sr;
}