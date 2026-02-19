"""
VCF Parser for SurakshaRx v4.1 — FIXED
Key fix: Now reads the GT (genotype) field and skips 0/0 (homozygous reference)
variants. A 0/0 genotype means the patient does NOT carry the alternate allele,
so those variants must not be counted toward diplotype calling.

Supported genotypes:
  0/0 or 0|0  → homozygous reference → SKIP (patient does not carry variant)
  0/1 or 0|1  → heterozygous        → INCLUDE (patient carries one copy)
  1/0 or 1|0  → heterozygous        → INCLUDE
  1/1 or 1|1  → homozygous alt      → INCLUDE (patient carries two copies)
  ./.         → missing / no call   → SKIP
"""

import re
from typing import Dict, List, Optional

TARGET_GENES = {"CYP2D6", "CYP2C19", "CYP2C9", "SLCO1B1", "TPMT", "DPYD"}


def parse_info_field(info_str: str) -> Dict[str, str]:
    info_dict = {}
    for item in info_str.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            info_dict[key.strip()] = value.strip()
        else:
            info_dict[item.strip()] = True
    return info_dict


def parse_genotype(format_str: str, sample_str: str) -> Optional[str]:
    """
    Extract the GT (genotype) value from FORMAT and SAMPLE columns.
    Returns the raw GT string e.g. '0/0', '0/1', '1/1', or None if not found.
    """
    if not format_str or not sample_str:
        return None
    fields = format_str.split(":")
    values = sample_str.split(":")
    if "GT" not in fields:
        return None
    gt_idx = fields.index("GT")
    if gt_idx >= len(values):
        return None
    return values[gt_idx]


def patient_carries_variant(gt: Optional[str]) -> bool:
    """
    Returns True if the genotype indicates the patient carries the ALT allele.

    Rules:
      0/0, 0|0  → False  (homozygous reference — patient does NOT carry variant)
      ./.       → False  (missing / no call)
      0/1, 1/0  → True   (heterozygous — one copy of variant)
      1/1, 1|1  → True   (homozygous alt — two copies)
      None      → True   (no FORMAT column present; include by default for
                          older VCF files that lack sample columns)
    """
    if gt is None:
        # No genotype info at all — older/simple VCFs without sample columns.
        # Include the variant (legacy behaviour).
        return True

    # Normalise separator
    gt_norm = gt.replace("|", "/")

    # Missing / no-call
    if "./." in gt_norm or gt_norm == ".":
        return False

    # Parse allele values
    alleles = gt_norm.split("/")
    try:
        allele_ints = [int(a) for a in alleles if a != "."]
    except ValueError:
        # Unexpected format — include conservatively
        return True

    # Patient carries the variant if ANY allele is non-reference (> 0)
    return any(a > 0 for a in allele_ints)


def is_homozygous_alt(gt: Optional[str]) -> bool:
    """Returns True when the patient carries TWO copies of the ALT allele (1/1)."""
    if gt is None:
        return False
    gt_norm = gt.replace("|", "/")
    alleles = gt_norm.split("/")
    try:
        allele_ints = [int(a) for a in alleles if a != "."]
    except ValueError:
        return False
    return len(allele_ints) >= 2 and all(a > 0 for a in allele_ints)


def parse_vcf(file_content: str) -> Dict:
    lines = file_content.strip().split("\n")
    metadata: Dict[str, str] = {}
    variants: List[Dict] = []
    parse_errors: List[str] = []
    header_cols: List[str] = []

    for line_num, line in enumerate(lines, 1):
        line = line.strip()
        if not line:
            continue

        # ── Meta-information lines ──────────────────────────────────────────
        if line.startswith("##"):
            if "=" in line:
                key = line[2:line.index("=")]
                val = line[line.index("=") + 1:]
                metadata[key] = val
            continue

        # ── Header line ─────────────────────────────────────────────────────
        if line.startswith("#CHROM"):
            header_cols = line[1:].split("\t")
            continue

        # ── Data lines ──────────────────────────────────────────────────────
        try:
            parts = line.split("\t")
            if len(parts) < 8:
                parts = line.split()
            if len(parts) < 8:
                parse_errors.append(f"Line {line_num}: insufficient columns ({len(parts)})")
                continue

            chrom      = parts[0]
            pos        = parts[1]
            rsid_raw   = parts[2]
            ref        = parts[3]
            alt        = parts[4]
            qual       = parts[5]
            filter_val = parts[6]
            info_str   = parts[7]

            # FORMAT and SAMPLE columns (columns 8 and 9)
            format_str  = parts[8]  if len(parts) > 8 else None
            sample_str  = parts[9]  if len(parts) > 9 else None

            rsid = rsid_raw if rsid_raw != "." else None

            # ── KEY FIX: parse and check genotype ──────────────────────────
            gt = parse_genotype(format_str, sample_str)

            if not patient_carries_variant(gt):
                # Patient is homozygous reference (0/0) or no-call — skip
                continue

            # ── Parse INFO ─────────────────────────────────────────────────
            info = parse_info_field(info_str)

            gene        = info.get("GENE", info.get("gene", None))
            star_allele = info.get("STAR", info.get("star", None))
            if star_allele and not star_allele.startswith("*"):
                star_allele = f"*{star_allele}"

            # Fall back to rsID-based gene inference if INFO lacks GENE=
            if not gene and rsid:
                gene = infer_gene_from_rsid(rsid)

            if gene and gene.upper() in TARGET_GENES:
                # Determine zygosity for diplotype weighting
                zygosity = "homozygous" if is_homozygous_alt(gt) else "heterozygous"

                variants.append({
                    "chrom":             chrom,
                    "pos":               pos,
                    "rsid":              rsid or f"chr{chrom}:{pos}",
                    "ref":               ref,
                    "alt":               alt,
                    "gene":              gene.upper(),
                    "star_allele":       star_allele,
                    "quality":           qual,
                    "filter":            filter_val,
                    "info":              info,
                    "functional_status": info.get("FUNC", info.get("FUNCTION",
                                                  info.get("function", "Unknown"))),
                    "genotype":          gt,
                    "zygosity":          zygosity,
                })

        except Exception as e:
            parse_errors.append(f"Line {line_num}: {str(e)}")

    detected_genes   = list(set(v["gene"] for v in variants))
    variants_by_gene: Dict[str, List[Dict]] = {}
    for v in variants:
        variants_by_gene.setdefault(v["gene"], []).append(v)

    return {
        "variants":            variants,
        "variants_by_gene":    variants_by_gene,
        "metadata":            metadata,
        "detected_genes":      detected_genes,
        "parse_errors":        parse_errors,
        "total_variants":      len(variants),
        "vcf_parsing_success": len(parse_errors) == 0 or len(variants) > 0,
    }


def infer_gene_from_rsid(rsid: str) -> Optional[str]:
    rsid_gene_map = {
        # CYP2D6
        "rs3892097":  "CYP2D6", "rs5030655":  "CYP2D6", "rs35742686": "CYP2D6",
        "rs1065852":  "CYP2D6", "rs28371725":  "CYP2D6", "rs16947":    "CYP2D6",
        "rs28371706": "CYP2D6", "rs59421388":  "CYP2D6", "rs1135840":  "CYP2D6",
        # CYP2C19
        "rs4244285":  "CYP2C19", "rs4986893":  "CYP2C19", "rs28399504": "CYP2C19",
        "rs56337013": "CYP2C19", "rs12248560":  "CYP2C19", "rs12769205": "CYP2C19",
        "rs17884712": "CYP2C19",
        # CYP2C9
        "rs1799853":  "CYP2C9", "rs1057910":  "CYP2C9", "rs28371686": "CYP2C9",
        "rs9332131":  "CYP2C9", "rs72558187": "CYP2C9",
        # SLCO1B1
        "rs4149056":  "SLCO1B1", "rs2306283":  "SLCO1B1", "rs11045819": "SLCO1B1",
        # TPMT
        "rs1800460":  "TPMT", "rs1142345":  "TPMT", "rs1800462":  "TPMT",
        "rs1800584":  "TPMT",
        # DPYD
        "rs3918290":  "DPYD", "rs55886062": "DPYD", "rs67376798": "DPYD",
        "rs75017182": "DPYD", "rs1801265":  "DPYD", "rs1801159":  "DPYD",
    }
    return rsid_gene_map.get(rsid)


def determine_diplotype(variants_for_gene: List[Dict]) -> str:
    """
    Build a diplotype string from the variants the patient actually carries.
    Homozygous alt (1/1) variants contribute the same star allele twice.
    """
    star_alleles: List[str] = []
    for v in variants_for_gene:
        sa = v.get("star_allele")
        if not sa:
            continue
        star_alleles.append(sa)
        # If homozygous, add the allele a second time
        if v.get("zygosity") == "homozygous":
            star_alleles.append(sa)

    if len(star_alleles) >= 2:
        return f"{star_alleles[0]}/{star_alleles[1]}"
    elif len(star_alleles) == 1:
        return f"{star_alleles[0]}/*1"
    return "*1/*1"


def get_sample_vcf() -> str:
    """Return a minimal valid VCF for UI demo when no file is uploaded."""
    return """##fileformat=VCFv4.2
##fileDate=20240101
##source=SurakshaRxTest
##reference=GRCh38
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">
##INFO=<ID=STAR,Number=1,Type=String,Description="Star allele">
##INFO=<ID=FUNCTION,Number=1,Type=String,Description="Functional status">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
chr22\t42522613\trs4244285\tG\tA\t100\tPASS\tGENE=CYP2C19;STAR=*2;FUNCTION=no_function\tGT\t1/1
chr22\t42523943\trs4986893\tG\tA\t100\tPASS\tGENE=CYP2C19;STAR=*3;FUNCTION=no_function\tGT\t0/1
chr22\t42526694\trs1065852\tG\tA\t95\tPASS\tGENE=CYP2D6;STAR=*4;FUNCTION=no_function\tGT\t0/1
chr12\t21331549\trs4149056\tT\tC\t98\tPASS\tGENE=SLCO1B1;STAR=*5;FUNCTION=decreased_function\tGT\t0/1
chr6\t18143955\trs1800460\tC\tT\t99\tPASS\tGENE=TPMT;STAR=*3B;FUNCTION=no_function\tGT\t0/1
chr1\t97915614\trs3918290\tC\tT\t97\tPASS\tGENE=DPYD;STAR=*2A;FUNCTION=no_function\tGT\t0/1
"""