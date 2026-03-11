#include <errno.h>
#include <getopt.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#define CMD_MAX 8192

typedef struct {
    const char *reference;
    const char *read1;
    const char *read2;
    const char *output_prefix;
    int threads;
    int min_mapping_quality;
    int min_base_quality;
} Options;

static void print_usage(const char *prog) {
    fprintf(stderr,
            "Fast SNP calling wrapper (bwa + samtools + bcftools)\\n\\n"
            "Usage:\\n"
            "  %s -r ref.fa -1 reads_1.fq.gz -2 reads_2.fq.gz -o sample [-t 8] [--min-mq 20] [--min-bq 20]\\n\\n"
            "Required arguments:\\n"
            "  -r, --reference      Reference FASTA\\n"
            "  -1, --read1          Read1 FASTQ(.gz)\\n"
            "  -2, --read2          Read2 FASTQ(.gz)\\n"
            "  -o, --output-prefix  Output prefix (usually sample name)\n"
            "  -s, --sample         Sample name alias of --output-prefix\n\n"
            "Optional arguments:\\n"
            "  -t, --threads        Thread number (default: 4)\\n"
            "      --min-mq         Minimum mapping quality for mpileup (default: 20)\\n"
            "      --min-bq         Minimum base quality for mpileup (default: 20)\\n"
            "  -h, --help           Show this help message\\n\\n"
            "Outputs:\\n"
            "  <prefix>.sorted.bam\\n"
            "  <prefix>.sorted.bam.bai\\n"
            "  <prefix>.snp.vcf.gz\\n"
            "  <prefix>.snp.vcf.gz.tbi\\n",
            prog);
}

static bool file_exists(const char *path) {
    struct stat st;
    return stat(path, &st) == 0;
}

static int run_command(const char *label, const char *cmd) {
    int rc;

    fprintf(stderr, "[INFO] %s\\n", label);
    fprintf(stderr, "[CMD ] %s\\n", cmd);

    rc = system(cmd);
    if (rc != 0) {
        fprintf(stderr, "[ERROR] Command failed with exit code %d: %s\\n", rc, label);
        return 1;
    }
    return 0;
}

static int validate_options(const Options *opt) {
    if (!opt->reference || !opt->read1 || !opt->read2 || !opt->output_prefix) {
        fprintf(stderr, "[ERROR] Missing required arguments.\\n");
        return 1;
    }
    if (opt->threads < 1) {
        fprintf(stderr, "[ERROR] Thread number must be >= 1.\\n");
        return 1;
    }
    if (opt->min_mapping_quality < 0 || opt->min_base_quality < 0) {
        fprintf(stderr, "[ERROR] Quality thresholds must be >= 0.\\n");
        return 1;
    }

    if (!file_exists(opt->reference)) {
        fprintf(stderr, "[ERROR] Reference not found: %s\\n", opt->reference);
        return 1;
    }
    if (!file_exists(opt->read1)) {
        fprintf(stderr, "[ERROR] Read1 not found: %s\\n", opt->read1);
        return 1;
    }
    if (!file_exists(opt->read2)) {
        fprintf(stderr, "[ERROR] Read2 not found: %s\\n", opt->read2);
        return 1;
    }

    return 0;
}

int main(int argc, char **argv) {
    Options opt = {0};
    int c;
    char bwt_index[4096];
    char bam[4096];
    char vcf[4096];
    char cmd[CMD_MAX];

    static struct option long_options[] = {
        {"reference", required_argument, 0, 'r'},
        {"read1", required_argument, 0, '1'},
        {"read2", required_argument, 0, '2'},
        {"output-prefix", required_argument, 0, 'o'},
        {"sample", required_argument, 0, 's'},
        {"threads", required_argument, 0, 't'},
        {"min-mq", required_argument, 0, 1000},
        {"min-bq", required_argument, 0, 1001},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    opt.threads = 4;
    opt.min_mapping_quality = 20;
    opt.min_base_quality = 20;

    while ((c = getopt_long(argc, argv, "r:1:2:o:s:t:h", long_options, NULL)) != -1) {
        switch (c) {
            case 'r':
                opt.reference = optarg;
                break;
            case '1':
                opt.read1 = optarg;
                break;
            case '2':
                opt.read2 = optarg;
                break;
            case 'o':
                opt.output_prefix = optarg;
                break;
            case 's':
                opt.output_prefix = optarg;
                break;
            case 't':
                opt.threads = atoi(optarg);
                break;
            case 1000:
                opt.min_mapping_quality = atoi(optarg);
                break;
            case 1001:
                opt.min_base_quality = atoi(optarg);
                break;
            case 'h':
                print_usage(argv[0]);
                return 0;
            default:
                print_usage(argv[0]);
                return 1;
        }
    }

    if (validate_options(&opt) != 0) {
        print_usage(argv[0]);
        return 1;
    }

    if (snprintf(bwt_index, sizeof(bwt_index), "%s.bwt", opt.reference) >= (int)sizeof(bwt_index)) {
        fprintf(stderr, "[ERROR] Reference path too long.\\n");
        return 1;
    }
    if (snprintf(bam, sizeof(bam), "%s.sorted.bam", opt.output_prefix) >= (int)sizeof(bam)) {
        fprintf(stderr, "[ERROR] Output prefix too long.\\n");
        return 1;
    }
    if (snprintf(vcf, sizeof(vcf), "%s.snp.vcf.gz", opt.output_prefix) >= (int)sizeof(vcf)) {
        fprintf(stderr, "[ERROR] Output prefix too long.\\n");
        return 1;
    }

    if (!file_exists(bwt_index)) {
        if (snprintf(cmd, sizeof(cmd), "bwa index '%s'", opt.reference) >= (int)sizeof(cmd)) {
            fprintf(stderr, "[ERROR] Command too long when building bwa index command.\\n");
            return 1;
        }
        if (run_command("Build bwa index", cmd) != 0) {
            return 1;
        }
    } else {
        fprintf(stderr, "[INFO] BWA index detected, skip indexing.\\n");
    }

    if (snprintf(cmd,
                 sizeof(cmd),
                 "bash -lc \"set -euo pipefail; bwa mem -t %d '%s' '%s' '%s' | "
                 "samtools view -@ %d -bS - | samtools sort -@ %d -o '%s' -\"",
                 opt.threads,
                 opt.reference,
                 opt.read1,
                 opt.read2,
                 opt.threads,
                 opt.threads,
                 bam) >= (int)sizeof(cmd)) {
        fprintf(stderr, "[ERROR] Command too long when building align/sort command.\\n");
        return 1;
    }
    if (run_command("Align reads and sort BAM", cmd) != 0) {
        return 1;
    }

    if (snprintf(cmd, sizeof(cmd), "samtools index -@ %d '%s'", opt.threads, bam) >= (int)sizeof(cmd)) {
        fprintf(stderr, "[ERROR] Command too long when building bam index command.\\n");
        return 1;
    }
    if (run_command("Index BAM", cmd) != 0) {
        return 1;
    }

    if (snprintf(cmd,
                 sizeof(cmd),
                 "bash -lc \"set -euo pipefail; "
                 "bcftools mpileup --threads %d -q %d -Q %d -Ou -f '%s' '%s' | "
                 "bcftools call --threads %d -mv -Oz -o '%s'\"",
                 opt.threads,
                 opt.min_mapping_quality,
                 opt.min_base_quality,
                 opt.reference,
                 bam,
                 opt.threads,
                 vcf) >= (int)sizeof(cmd)) {
        fprintf(stderr, "[ERROR] Command too long when building SNP calling command.\\n");
        return 1;
    }
    if (run_command("Call SNPs", cmd) != 0) {
        return 1;
    }

    if (snprintf(cmd, sizeof(cmd), "bcftools index -t '%s'", vcf) >= (int)sizeof(cmd)) {
        fprintf(stderr, "[ERROR] Command too long when building VCF index command.\\n");
        return 1;
    }
    if (run_command("Index VCF", cmd) != 0) {
        return 1;
    }

    fprintf(stderr, "[DONE] SNP calling finished. Output: %s\\n", vcf);
    return 0;
}
