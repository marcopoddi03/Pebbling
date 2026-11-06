#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "diff_matrix.h"


int main(int argc, char **argv) {
    if (argc < 4) {
        fprintf(stderr, "Usage: %s <file_in.vcf.gz> <matrix_out.txt> <chr1> [chr2 chr3 ...]\n", argv[0]);
        fprintf(stderr, "Example: %s data.vcf.gz chr3 chr8 chr15 chr19 chr21\n", argv[0]);
        return 1;
    }

    char *filename_in = argv[1];
    char *filename_out= argv[2];

    char *regions = format_regions(argv, 3, argc);
    if (!regions)
        return 1;

    PairwiseMatrix *pairwise_matrix = get_bcf_pairwise_matrix(filename_in, regions);

    if (!pairwise_matrix) {
        free(regions);
        return 1;
    }

    int saved = save_pairwise_matrix_txt(filename_out, pairwise_matrix);

    free_pairwise_matrix(pairwise_matrix);
    free(regions);
    
    if(saved)
        return 1;

    return 0;
}