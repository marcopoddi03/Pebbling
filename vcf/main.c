#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <htslib/hts.h>
#include <stdint.h>
#include <htslib/synced_bcf_reader.h>
#include "diff_matrix.h"


int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <file.vcf.gz> <chr1> [chr2 chr3 ...]\n", argv[0]);
        fprintf(stderr, "Example: %s data.vcf.gz chr3 chr8 chr15 chr19 chr21\n", argv[0]);
        return 1;
    }

    char *filename = argv[1];

    char *regions = format_regions(argv, 2, argc);
    if (!regions)
        return 1;

    bcf_srs_t *sr = set_reader_for_diff(filename, regions);
    if (!sr) {
        free(regions);
        return 1;
    }

    bcf_hdr_t *hdr = bcf_sr_get_header(sr, 0);
    int n_samples = bcf_hdr_nsamples(hdr);

    int **diff_matrix = init_diff_matrix(n_samples);
    
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

                if(gt_cmp(gt_i, gt_j, max_ploidy) == BCF_CMP_MISMATCH) 
                    update_diff_matrix(diff_matrix, i, j);
            }
        } 
        free(gt_arr);  
    }

    bool saved = save_diff_matrix("matrix.bin", diff_matrix, n_samples);
    bcf_sr_destroy(sr);
    free_diff_matrix(diff_matrix);
    free(regions);
    
    if(!saved)
        return 1;

    return 0;
}