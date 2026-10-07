#!/usr/bin/env bash
###
#  Extract FORMAT/LRR and FORMAT/THETA from a gtc2vcf or affy2vcf VCF/BCF into two
#  headerless, whitespace-delimited, gzipped matrices (e.g. as input for HI-CNV)
#
#  Outputs (PREFIX given with -o):
#    PREFIX.LRR.txt.gz     LRR matrix
#    PREFIX.theta.txt.gz   THETA matrix
#    PREFIX.variants.txt   one line per variant: CHROM POS ID REF ALT ALLELE_A ALLELE_B
#    PREFIX.samples.txt    one line per sample
#
#  By default each matrix row is a variant and each column a sample (the layout used
#  by HI-CNV), in the order of PREFIX.variants.txt and PREFIX.samples.txt. With -t the
#  matrices are transposed (rows are samples), which requires holding them in memory.
#
#  THETA is written as stored in the VCF, i.e. oriented towards allele B (see the
#  ALLELE_A/ALLELE_B columns of PREFIX.variants.txt). Missing values, NaN, and
#  infinite values are written as NA.
###

set -euo pipefail

usage() {
    cat <<EOF
Usage: $(basename "$0") [options] -o <prefix> <in.vcf|in.bcf>

Options:
    -o PREFIX    output prefix (required)
    -r REGION    restrict to region(s), e.g. 1 or 1:1000000-2000000 (requires an indexed input)
    -i EXPR      only include variants for which EXPR is true (bcftools filtering expression)
    -S FILE      file with the samples to include, one per line, in the desired column order
    -d INT       number of decimal places [3]
    -t           transpose: rows are samples and columns are variants
    -h           print this help message
EOF
    exit "${1:-0}"
}

prefix=""
region=""
include=""
samples_file=""
decimals=3
transpose=0
while getopts "o:r:i:S:d:th" opt; do
    case $opt in
    o) prefix=$OPTARG ;;
    r) region=$OPTARG ;;
    i) include=$OPTARG ;;
    S) samples_file=$OPTARG ;;
    d) decimals=$OPTARG ;;
    t) transpose=1 ;;
    h) usage 0 ;;
    *) usage 1 >&2 ;;
    esac
done
shift $((OPTIND - 1))
[ $# -eq 1 ] || usage 1 >&2
[ -n "$prefix" ] || { echo "Error: missing -o <prefix>" >&2; exit 1; }
[[ $decimals =~ ^[0-9]+$ ]] || { echo "Error: -d expects a non-negative integer" >&2; exit 1; }
vcf=$1

command -v bcftools >/dev/null || { echo "Error: bcftools not found in PATH" >&2; exit 1; }
for tag in FORMAT/LRR FORMAT/THETA INFO/ALLELE_A INFO/ALLELE_B; do
    if ! bcftools view -h "$vcf" | grep -q "^##${tag%%/*}=<ID=${tag#*/},"; then
        echo "Error: $tag is not defined in the header of $vcf" >&2
        [ "$tag" = FORMAT/THETA ] && echo "       (for affy2vcf, request it with --tags ...,THETA)" >&2
        exit 1
    fi
done

sample_opts=()
[ -n "$samples_file" ] && sample_opts+=(--samples-file "$samples_file" --force-samples)
query_opts=("${sample_opts[@]+"${sample_opts[@]}"}")
[ -n "$region" ] && query_opts+=(--regions "$region")
[ -n "$include" ] && query_opts+=(--include "$include")

# take sample names from the header after subsetting, as bcftools query --list-samples does not follow the
# order of the --samples-file while the queried values do
bcftools view --header-only "${sample_opts[@]+"${sample_opts[@]}"}" "$vcf" | tail -n 1 | cut -f 10- | tr '\t' '\n' |
    { grep -v '^$' || true; } >"$prefix.samples.txt"
nsmpl=$(wc -l <"$prefix.samples.txt" | tr -d ' ')
[ "$nsmpl" -gt 0 ] || { echo "Error: no samples selected" >&2; exit 1; }

# one pass over the VCF: variant info, then nsmpl LRR values, then nsmpl THETA values
bcftools query "${query_opts[@]+"${query_opts[@]}"}" \
    -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/ALLELE_A\t%INFO/ALLELE_B[\t%LRR][\t%THETA]\n' "$vcf" |
    awk -F '\t' -v n="$nsmpl" -v d="$decimals" -v t="$transpose" \
        -v var="$prefix.variants.txt" -v lrr="$prefix.LRR.txt.gz" -v theta="$prefix.theta.txt.gz" '
    function fmt(x) {
        # missing (.), NaN, and +/-Inf are written as NA
        if (x == "." || x == "" || tolower(x) ~ /nan|inf/) return "NA"
        y = sprintf(numfmt, x)
        return y ~ /^-0\.?0*$/ ? substr(y, 2) : y  # avoid "-0.000"
    }
    BEGIN {
        numfmt = "%." d "f"
        lrr_cmd = "gzip -c > \"" lrr "\""
        theta_cmd = "gzip -c > \"" theta "\""
    }
    {
        if (NF != 7 + 2 * n) {
            printf "Error: expected %d columns but found %d at %s:%s\n", 7 + 2 * n, NF, $1, $2 > "/dev/stderr"
            err = 1
            exit 1
        }
        print $1, $2, $3, $4, $5, $6, $7 > var
        if (t) {
            # store values to write the transposed matrices at the end
            for (i = 1; i <= n; i++) {
                L[i] = L[i] (NR > 1 ? " " : "") fmt($(7 + i))
                T[i] = T[i] (NR > 1 ? " " : "") fmt($(7 + n + i))
            }
        } else {
            row = fmt($8)
            for (i = 2; i <= n; i++) row = row " " fmt($(7 + i))
            print row | lrr_cmd
            row = fmt($(8 + n))
            for (i = 2; i <= n; i++) row = row " " fmt($(7 + n + i))
            print row | theta_cmd
        }
    }
    END {
        if (err) exit 1
        if (NR == 0) {
            print "Error: no variants selected" > "/dev/stderr"
            exit 1
        }
        if (t) for (i = 1; i <= n; i++) {
            print L[i] | lrr_cmd
            print T[i] | theta_cmd
        }
        close(lrr_cmd)
        close(theta_cmd)
        printf "Wrote %d variants x %d samples to %s and %s\n", NR, n, lrr, theta > "/dev/stderr"
    }'
