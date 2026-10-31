#include <stdio.h>
#include <stdlib.h>
#include <htslib/hts.h>
#include <htslib/synced_bcf_reader.h>

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
    if (!matrix) {
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
    }


    bcf_sr_destroy(sr);
    free(diff_matrix);
    free(data);
    return 0;
}