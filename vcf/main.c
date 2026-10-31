#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
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

    DiffMatrix *diff_matrix = get_bcf_diff_matrix(filename, regions);

    if (!diff_matrix) {
        free(regions);
        return 1;
    }

    bool saved = save_diff_matrix("matrix.bin", diff_matrix);

    free_diff_matrix(diff_matrix);
    free(regions);
    
    if(!saved)
        return 1;

    return 0;
}