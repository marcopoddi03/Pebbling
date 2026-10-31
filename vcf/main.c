#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <htslib/hts.h>
#include <htslib/synced_bcf_reader.h>


bool gt_cmp(const int32_t *gt_i, const int32_t *gt_j, int max_ploidy);
int int_cmp(const void *a, const void *b);

int int_cmp(const void *a, const void *b) {
    int int_a = *((const int *)a);
    int int_b = *((const int *)b);
    return (int_a > int_b) - (int_a < int_b);
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



int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <file.vcf.gz> <chr1> [chr2 chr3 ...]\n", argv[0]);
        fprintf(stderr, "Example: %s data.vcf.gz chr3 chr8 chr15 chr19 chr21\n", argv[0]);
        return 1;
    }

    char *filename = argv[1];

    //List of chroms
    size_t len = 1;

    for (int i = 2; i < argc; i++)
        len += strlen(argv[i]) + (i > 2 ? 1 : 0);

    char *regions = malloc(len);
    if (!regions)
        return 1;

    regions[0] = '\0';

    for (int i = 2; i < argc; i++) {
        if (i > 2)
            strcat(regions, ",");

        strcat(regions, argv[i]);
    }

    bcf_srs_t *sr = bcf_sr_init();
    if (!sr) {
        fprintf(stderr, "Error: could not inizialize reader\n");
        return 1;
    }

    if (bcf_sr_set_opt(sr, BCF_SR_REQUIRE_IDX)) {
        fprintf(stderr, "Error: could not recognise option\n");
        bcf_sr_destroy(sr);
        return 1;
    }

    if (!bcf_sr_add_reader(sr, filename)) {
        fprintf(stderr, "Error: %s\n", bcf_sr_strerror(sr->errnum));
        bcf_sr_destroy(sr);
        return 1;
    }

    if (bcf_sr_set_regions(sr, regions, 0) < 0) {
        fprintf(stderr, "Error: could not set regions\n");
        bcf_sr_destroy(sr);
        return 1;
    }

    bcf_hdr_t *hdr = bcf_sr_get_header(sr, 0);
    int n_samples = bcf_hdr_nsamples(hdr);

    int **diff_matrix = malloc(n_samples * sizeof(int *));
    if (!diff_matrix) {
        fprintf(stderr, "Error: memory allocation failed\n");
        return 1;
    }

    int *data = calloc(n_samples * n_samples, sizeof(int));
    if (!data) {
        fprintf(stderr, "Error: memory allocation failed\n");
        return 1;
    }

    for (int i = 0; i < n_samples; i++) {
        diff_matrix[i] = data + (i * n_samples);
    }

    while (bcf_sr_next_line(sr)) {
        bcf1_t *rec = bcf_sr_get_line(sr, 0);
        if (!rec)
            continue;
        bcf_unpack(rec, BCF_UN_ALL);

        int32_t *gt_arr = NULL;
        int ngt_arr = 0;

        int ngt = bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr);
        if (ngt <= 0)
            continue;


        int max_ploidy = ngt / n_samples;

        for (int i = 0; i < n_samples; i++) {
            int32_t *gt_i = gt_arr + i * max_ploidy;

            for (int j = i + 1; j < n_samples; j++) {

                int32_t *gt_j = gt_arr + j * max_ploidy;

                gt_cmp(gt_i, gt_j, max_ploidy);
            }
        } 
        free(gt_arr);  
    }


    bcf_sr_destroy(sr);
    free(diff_matrix);
    free(data);
    
    return 0;
}