#!/usr/bin/env bash

echo Download the SILVA database assets

export RDIR="SILVA_v138"
mkdir -p "${RDIR}"
chmod -R 777 "${RDIR}"
cd "${RDIR}"

for URL in \
  "https://www.arb-silva.de/fileadmin/silva_databases/release_138/Exports/taxonomy/tax_slv_ssu_138.txt.gz" \
  "https://www.arb-silva.de/fileadmin/silva_databases/release_138/Exports/taxonomy/tax_slv_ssu_138.tre.gz" \
  "https://www.arb-silva.de/fileadmin/silva_databases/release_138/Exports/taxonomy/taxmap_slv_ssu_ref_nr_138.txt.gz" \
  "https://www.arb-silva.de/fileadmin/silva_databases/release_138/Exports/SILVA_138_SSURef_NR99_tax_silva_trunc.fasta.gz" \
  "https://www.arb-silva.de/fileadmin/silva_databases/release_138/Exports/SILVA_138_SSURef_NR99_tax_silva_full_align_trunc.fasta.gz"
  do
    FEXT=$(basename ${URL})
    echo "Fetch ${FEXT}"
    curl -fsSLO ${URL}
    echo "Unpack ${FEXT}"
    gunzip "${FEXT}"
  done

echo The SILVA database assets download is complete
