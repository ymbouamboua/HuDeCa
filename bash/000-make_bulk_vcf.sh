#!/usr/bin/env bash
# ============================================================
# make_bulk_vcf.sh — FINAL HARDENED VERSION
# ------------------------------------------------------------
# Merge multiple per-donor VCFs into a single bulk VCF,
# filter by MAC and F_MISSING, index the VCF,
# and output a summary.
#
# USAGE:
#   bash make_bulk_vcf.sh <vcf_dir> "<pattern>" <output.vcf.gz> [MAC] [F_MISSING]
#
# EXAMPLE:
#   bash make_bulk_vcf.sh data "*.vcf.gz" bulk.vcf.gz 1 0.2
# ============================================================

set -euo pipefail

# ------------------------
# Arguments
# ------------------------
vcf_dir="$1"
pattern="$2"
output_vcf="$3"
mac_thr="${4:-1}"
fmiss_thr="${5:-0.2}"

summary_file="${output_vcf%.vcf.gz}_summary.txt"

# ------------------------
# Sanity checks
# ------------------------
if [[ ! -d "$vcf_dir" ]]; then
    echo "ERROR: VCF directory not found: $vcf_dir"
    exit 1
fi

if [[ -z "$pattern" || -z "$output_vcf" ]]; then
    echo "Usage: $0 <vcf_dir> \"<pattern>\" <output.vcf.gz> [MAC] [F_MISSING]"
    exit 1
fi

# ------------------------
# Resolve VCF list SAFELY
# ------------------------
shopt -s nullglob
vcf_files=( "$vcf_dir"/$pattern )
shopt -u nullglob

if [[ ${#vcf_files[@]} -eq 0 ]]; then
    echo "ERROR: No VCF files found matching pattern '$pattern' in $vcf_dir"
    exit 1
fi

echo "=================================================="
echo "VCF files detected (${#vcf_files[@]}):"
printf "  - %s\n" "${vcf_files[@]}"
echo "=================================================="

# ------------------------
# Temp files (safe)
# ------------------------
tmp_joint="$(mktemp bulk_joint.vcf.gz)"
tmp_filt="$(mktemp bulk_joint_filtered.vcf.gz)"

cleanup() {
    rm -f "$tmp_joint" "$tmp_joint.tbi" "$tmp_joint.csi" \
          "$tmp_filt" "$tmp_filt.tbi" "$tmp_filt.csi"
}
trap cleanup EXIT

# ==================================================
# STEP 1 — Index individual VCFs
# ==================================================
echo "Step 1. Index individual VCFs if needed"
echo "--------------------------------------------------"

for f in "${vcf_files[@]}"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: Missing file $f"
        exit 1
    fi

    if [[ ! -f "${f}.tbi" && ! -f "${f}.csi" ]]; then
        echo "  Indexing $f"
        bcftools index --force "$f"
    else
        echo "  Already indexed: $f"
    fi
done

# ==================================================
# STEP 2 — Merge + sort (CRITICAL)
# ==================================================
echo "=================================================="
echo "Step 2. Merge + sort VCFs"
echo "=================================================="

bcftools merge --threads 8 -m all "${vcf_files[@]}" \
  | bcftools sort -Oz -o "$tmp_joint"

bcftools index --force "$tmp_joint"
bcftools index -n "$tmp_joint" >/dev/null

echo "Merged VCF: $tmp_joint"

# ==================================================
# STEP 3 — Filter variants
# ==================================================
echo "=================================================="
echo "Step 3. Filter (MAC >= $mac_thr, F_MISSING < $fmiss_thr)"
echo "=================================================="

bcftools view \
  -i "MAC>=${mac_thr} && F_MISSING<${fmiss_thr}" \
  -Oz -o "$tmp_filt" "$tmp_joint"

bcftools index --force "$tmp_filt"
bcftools index -n "$tmp_filt" >/dev/null

echo "Filtered VCF: $tmp_filt"

# ==================================================
# STEP 4 — Move final output
# ==================================================
echo "=================================================="
echo "Step 4. Final output"
echo "=================================================="

mv "$tmp_filt" "$output_vcf"
mv "${tmp_filt}.tbi" "${output_vcf}.tbi" 2>/dev/null || true
mv "${tmp_filt}.csi" "${output_vcf}.csi" 2>/dev/null || true

echo "Final VCF written to: $output_vcf"

# ==================================================
# STEP 5 — Summary stats
# ==================================================
echo "=================================================="
echo "Step 5. Summary statistics"
echo "=================================================="

nsamples=$(bcftools query -l "$output_vcf" | wc -l)
nsnps=$(bcftools view -H "$output_vcf" | wc -l)

mean_missing=$(bcftools query -f '[%GT]\n' "$output_vcf" \
  | awk '
    {
      miss=0
      for(i=1;i<=NF;i++) if($i=="./.") miss++
      print miss/NF
    }' | awk '{s+=$1} END {if(NR>0) print s/NR; else print 0}'
)

cat << EOF | tee "$summary_file"
Bulk VCF Summary
====================
Samples: $nsamples
SNPs: $nsnps
Mean missing per SNP: $mean_missing
EOF

# Example:
#bash ./make_bulk_vcf.sh \
#  ~/projects/singlecell/000-hudeca-nose/data/processed/FreeBayes_bulk_OUTDIR \
#  "*.vcf.gz" \
#  ~/projects/singlecell/000-hudeca/data/processed/FreeBayes_bulk_OUTDIR/bulk_rna.vcf.gz
