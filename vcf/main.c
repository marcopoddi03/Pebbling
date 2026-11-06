#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "vcf_pairwise_matrices.h"


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

    PairwiseMatrices *pairwise_matrices = get_vcf_pairwise_matrices(filename_in, regions);

    if (!pairwise_matrices) {
        free(regions);
        return 1;
    }

    int saved = save_pairwise_matrices_txt(filename_out, pairwise_matrices);

    free_pairwise_matrices(pairwise_matrices);
    free(regions);
    
    if(saved)
        return 1;

    return 0;
}