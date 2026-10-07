#!/usr/bin/env bash
set -euo pipefail

# Public/portfolio version of a genotype-imputation workflow.
# Replace placeholder paths before use.

BCFTOOLS="${BCFTOOLS:-/path/to/bcftools}"
PLINK="${PLINK:-/path/to/plink}"
GLIMPSE_PHASE="${GLIMPSE_PHASE:-GLIMPSE_phase}"
GLIMPSE_LIGATE="${GLIMPSE_LIGATE:-GLIMPSE_ligate}"
REFERENCE_FASTA="${REFERENCE_FASTA:-/path/to/hs37d5.fa}"
REF_ROOT="${REF_ROOT:-/path/to/1000G/reference}"
BAM_LIST="${BAM_LIST:-/path/to/strict_bamlist}"
THREADS_GL="${THREADS_GL:-10}"
THREADS_INDEX="${THREADS_INDEX:-30}"
THREADS_GLIMPSE="${THREADS_GLIMPSE:-30}"
THREADS_LIGATE="${THREADS_LIGATE:-22}"
THREADS_MERGE="${THREADS_MERGE:-30}"
WORKDIR="${WORKDIR:-imputation}"

mkdir -p "${WORKDIR}"/{gl_calling,glimpse/glimpse_output}
cd "${WORKDIR}"

sample_name_from_bam() {
    local bam="$1"
    basename "${bam}" | cut -d'.' -f1 | cut -d'-' -f1
}

genotype_likelihoods_autosomes() {
    : > calling.sh
    for chr in $(seq 1 22); do
        while IFS= read -r bam; do
            [[ -z "${bam}" ]] && continue
            sample="$(sample_name_from_bam "${bam}")"
            sites_vcf="${REF_ROOT}/calling_bcftools/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.biallelic_snps.sites.vcf.gz"
            sites_tsv="${REF_ROOT}/calling_bcftools/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.biallelic_snps.sites.tsv.gz"
            out_bcf="gl_calling/${sample}.1kg.phase3.v5a.20130502.biallelic_snps.chr${chr}.bcf"
            printf '%q mpileup -f %q -Q 30 -I -E -a FORMAT/DP -T %q -r %q %q | %q call -Aim -C alleles -T %q -Ob -o %q\n' \
                "${BCFTOOLS}" "${REFERENCE_FASTA}" "${sites_vcf}" "${chr}" "${bam}" \
                "${BCFTOOLS}" "${sites_tsv}" "${out_bcf}" >> calling.sh
        done < "${BAM_LIST}"
    done
    parallel -j "${THREADS_GL}" < calling.sh
}

index_gl_bcf() {
    find gl_calling -maxdepth 1 -name '*.bcf' -print0 | \
        xargs -0 -I{} -n1 echo "${BCFTOOLS} index {}" | \
        parallel -j "${THREADS_INDEX}"
}

glimpse_autosomes() {
    : > glimpse.commands
    shopt -s nullglob
    for bcf in gl_calling/*.chr1.bcf; do
        sample="$(basename "${bcf}" | cut -d'.' -f1)"
        mkdir -p "glimpse/glimpse_output/${sample}"
        for chr in $(seq 1 22); do
            input_bcf="gl_calling/${sample}.1kg.phase3.v5a.20130502.biallelic_snps.chr${chr}.bcf"
            ref_bcf="${REF_ROOT}/glimpse_format/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.biallelic_snps.bcf"
            map_file="${REF_ROOT}/glimpse_format/maps/genetic_maps.b37/chr${chr}.b37.gmap.gz"
            chunks="${REF_ROOT}/glimpse_format/chunks.chr${chr}.2MB_200KB.txt"
            while IFS= read -r line || [[ -n "${line}" ]]; do
                [[ -z "${line}" ]] && continue
                id="$(awk '{printf "%02d", $1}' <<< "${line}")"
                input_region="$(awk '{print $3}' <<< "${line}")"
                output_region="$(awk '{print $4}' <<< "${line}")"
                out_vcf="glimpse/glimpse_output/${sample}/${sample}.1kg.phase3.v5a.20130502.biallelic_snps.glimpse.chr${chr}.${id}.vcf.gz"
                out_log="glimpse/glimpse_output/${sample}/${sample}.1kg.phase3.v5a.20130502.biallelic_snps.glimpse.chr${chr}.${id}.log"
                printf '%q --input %q --reference %q --map %q --input-region %q --output-region %q --output %q > %q\n' \
                    "${GLIMPSE_PHASE}" "${input_bcf}" "${ref_bcf}" "${map_file}" \
                    "${input_region}" "${output_region}" "${out_vcf}" "${out_log}" >> glimpse.commands
            done < "${chunks}"
        done
    done
    parallel -j "${THREADS_GLIMPSE}" < glimpse.commands
}

postprocess_autosomes() {
    mkdir -p imputed_chromosomes final_merge
    : > ligate.commands
    for sample_dir in glimpse/glimpse_output/*; do
        [[ -d "${sample_dir}" ]] || continue
        sample="$(basename "${sample_dir}")"
        for chr in $(seq 1 22); do
            list_file="${sample}.${chr}.list"
            find "${sample_dir}" -maxdepth 1 -name "*chr${chr}.*.vcf.gz" | sort > "${list_file}"
            out_bcf="imputed_chromosomes/${sample}.1kg.phase3.v5a.20130502.biallelic_snps.glimpse.chr${chr}.imputed.bcf"
            log_file="imputed_chromosomes/${sample}.${chr}.ligate.log"
            printf '%q --input %q --output %q > %q\n' \
                "${GLIMPSE_LIGATE}" "${list_file}" "${out_bcf}" "${log_file}" >> ligate.commands
        done
    done
    parallel -j "${THREADS_LIGATE}" < ligate.commands

    : > merging.commands
    : > concat.list
    for chr in $(seq 1 22); do
        merge_list="merge.${chr}.list"
        find imputed_chromosomes -maxdepth 1 -name "*chr${chr}.imputed.bcf" | sort > "${merge_list}"
        merged_bcf="imputed_chromosomes/merged.1kg.phase3.v5a.20130502.biallelic_snps.glimpse.chr${chr}.imputed.bcf"
        printf '%q merge -l %q -Ob -o %q\n' "${BCFTOOLS}" "${merge_list}" "${merged_bcf}" >> merging.commands
        echo "${merged_bcf}" >> concat.list
    done
    parallel -j "${THREADS_MERGE}" < merging.commands

    final_vcf="final_merge/merged.1kg.phase3.v5a.20130502.biallelic_snps.glimpse.imputed.vcf.gz"
    "${BCFTOOLS}" concat -f concat.list -Oz -o "${final_vcf}"
    "${PLINK}" --vcf "${final_vcf}" --vcf-min-gp 0.99 --make-bed --double-id \
        --threads "${THREADS_LIGATE}" \
        --out "final_merge/merged.1kg.phase3.v5a.20130502.biallelic_snps.glimpse.imputed.GP99"
}

update_metadata() {
    local prefix="${1:?Usage: update-metadata <plink_prefix>}"
    if [[ -f update-ids-glimpse.txt ]]; then
        "${PLINK}" --bfile "${prefix}" --update-ids update-ids-glimpse.txt --make-bed --out "${prefix}.ids"
        prefix="${prefix}.ids"
    fi
    if [[ -f update-sex.txt ]]; then
        "${PLINK}" --bfile "${prefix}" --update-sex update-sex.txt --make-bed --out "${prefix}.sex"
    fi
}

convert_snp_ids_to_position() {
    local prefix="${1:?Usage: convert-ids <plink_prefix>}"
    awk '{print $1"\t"$1":"$4"\t"$3"\t"$4"\t"$5"\t"$6}' "${prefix}.bim" > "${prefix}.bim.tmp"
    mv "${prefix}.bim.tmp" "${prefix}.bim"
}

chrX_workflow() {
    : > callingX.sh
    while IFS= read -r bam; do
        [[ -z "${bam}" ]] && continue
        sample="$(sample_name_from_bam "${bam}")"
        sites_vcf="${REF_ROOT}/calling_bcftools/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.biallelic_snps.sites.vcf.gz"
        sites_tsv="${REF_ROOT}/calling_bcftools/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.biallelic_snps.sites.tsv.gz"
        out_bcf="gl_calling/${sample}.1kg.phase3.v1c.20130502.biallelic_snps.chrX.bcf"
        printf '%q mpileup -f %q -Q 30 -I -E -a FORMAT/DP -T %q -r X %q | %q call -Aim -C alleles -T %q -Ob -o %q\n' \
            "${BCFTOOLS}" "${REFERENCE_FASTA}" "${sites_vcf}" "${bam}" "${BCFTOOLS}" "${sites_tsv}" "${out_bcf}" >> callingX.sh
    done < "${BAM_LIST}"
    parallel -j "${THREADS_GL}" < callingX.sh

    find gl_calling -maxdepth 1 -name '*chrX.bcf' -print0 | \
        xargs -0 -I{} -n1 echo "${BCFTOOLS} index {}" | parallel -j "${THREADS_INDEX}"

    : > glimpseX.commands
    shopt -s nullglob
    for bcf in gl_calling/*chrX.bcf; do
        sample="$(basename "${bcf}" | cut -d'.' -f1)"
        mkdir -p "glimpse/glimpse_output/${sample}"
        ref_bcf="${REF_ROOT}/glimpse_format/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.biallelic_snps.bcf"
        map_file="${REF_ROOT}/glimpse_format/maps/genetic_maps.b37/chrX.b37.gmap.gz"
        chunks="${REF_ROOT}/glimpse_format/chunks.chrX.2MB_200KB.txt"
        while IFS= read -r line || [[ -n "${line}" ]]; do
            [[ -z "${line}" ]] && continue
            id="$(awk '{printf "%02d", $1}' <<< "${line}")"
            input_region="$(awk '{print $3}' <<< "${line}")"
            output_region="$(awk '{print $4}' <<< "${line}")"
            out_vcf="glimpse/glimpse_output/${sample}/${sample}.1kg.phase3.v1c.20130502.biallelic_snps.glimpse.chrX.${id}.vcf.gz"
            out_log="glimpse/glimpse_output/${sample}/${sample}.1kg.phase3.v1c.20130502.biallelic_snps.glimpse.chrX.${id}.log"
            printf '%q --input %q --reference %q --map %q --input-region %q --output-region %q --output %q > %q\n' \
                "${GLIMPSE_PHASE}" "${bcf}" "${ref_bcf}" "${map_file}" "${input_region}" "${output_region}" "${out_vcf}" "${out_log}" >> glimpseX.commands
        done < "${chunks}"
    done
    parallel -j "${THREADS_LIGATE}" < glimpseX.commands

    mkdir -p imputed_chromosomes final_merge
    : > ligate.commandsX
    for sample_dir in glimpse/glimpse_output/*; do
        [[ -d "${sample_dir}" ]] || continue
        sample="$(basename "${sample_dir}")"
        list_file="${sample}.X.list"
        find "${sample_dir}" -maxdepth 1 -name '*chrX.*.vcf.gz' | sort > "${list_file}"
        out_bcf="imputed_chromosomes/${sample}.1kg.phase3.v1c.20130502.biallelic_snps.glimpse.chrX.imputed.bcf"
        log_file="imputed_chromosomes/${sample}.X.ligate.log"
        printf '%q --input %q --output %q > %q\n' "${GLIMPSE_LIGATE}" "${list_file}" "${out_bcf}" "${log_file}" >> ligate.commandsX
    done
    parallel -j "${THREADS_LIGATE}" < ligate.commandsX

    find imputed_chromosomes -maxdepth 1 -name '*chrX.imputed.bcf' | sort > merge.X.list
    merged_x="imputed_chromosomes/merged.1kg.phase3.v1c.20130502.biallelic_snps.glimpse.chrX.imputed.bcf"
    "${BCFTOOLS}" merge -l merge.X.list -Ob -o "${merged_x}"
    echo "${merged_x}" > concatX.list
    x_vcf="final_merge/merged.1kg.phase3.20130502.biallelic_snps.glimpse.imputed.justX.vcf.gz"
    "${BCFTOOLS}" concat -f concatX.list -Oz -o "${x_vcf}"
    "${PLINK}" --vcf "${x_vcf}" --vcf-min-gp 0.99 --make-bed --double-id \
        --out "final_merge/merged.1kg.phase3.20130502.biallelic_snps.glimpse.imputed.justX.GP99"
}

usage() {
    cat <<'USAGE'
Usage:
  bash imputation_pipeline.sh genotype-likelihoods
  bash imputation_pipeline.sh index
  bash imputation_pipeline.sh glimpse
  bash imputation_pipeline.sh postprocess
  bash imputation_pipeline.sh chrX
  bash imputation_pipeline.sh update-metadata <plink_prefix>
  bash imputation_pipeline.sh convert-ids <plink_prefix>
USAGE
}

case "${1:-}" in
    genotype-likelihoods) genotype_likelihoods_autosomes ;;
    index) index_gl_bcf ;;
    glimpse) glimpse_autosomes ;;
    postprocess) postprocess_autosomes ;;
    chrX) chrX_workflow ;;
    update-metadata) shift; update_metadata "$@" ;;
    convert-ids) shift; convert_snp_ids_to_position "$@" ;;
    *) usage; exit 1 ;;
esac
