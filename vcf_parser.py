"""
VCF Parser for PharmaGuard v1.0
Parses Variant Call Format files and extracts pharmacogenomic variants
"""

import re
from typing import Dict, List, Optional

# Target genes for pharmacogenomics analysis
TARGET_GENES = {"CYP2D6", "CYP2C19", "CYP2C9", "SLCO1B1", "TPMT", "DPYD"}


def parse_info_field(info_str: str) -> Dict[str, str]:
    """Parse the INFO field of a VCF record into a dictionary."""
    info_dict = {}
    for item in info_str.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            info_dict[key.strip()] = value.strip()
        else:
            info_dict[item.strip()] = True
    return info_dict


def parse_vcf(file_content: str) -> Dict:
    """
    Parse VCF file content and extract pharmacogenomic variants.
    
    Returns:
        dict with keys: variants, metadata, detected_genes, parse_errors
    """
    lines = file_content.strip().split("\n")
    
    metadata = {}
    variants = []
    parse_errors = []
    header_cols = []
    
    for line_num, line in enumerate(lines, 1):
        line = line.strip()
        
        if not line:
            continue
            
        # Parse metadata lines
        if line.startswith("##"):
            if "=" in line:
                key = line[2:line.index("=")]
                val = line[line.index("=") + 1:]
                metadata[key] = val
            continue
        
        # Parse header line
        if line.startswith("#CHROM"):
            header_cols = line[1:].split("\t")
            continue
        
        # Parse variant lines
        try:
            parts = line.split("\t")
            if len(parts) < 8:
                # Try space-separated
                parts = line.split()
            
            if len(parts) < 8:
                parse_errors.append(f"Line {line_num}: insufficient columns")
                continue
            
            chrom = parts[0]
            pos = parts[1]
            rsid = parts[2] if parts[2] != "." else None
            ref = parts[3]
            alt = parts[4]
            qual = parts[5]
            filter_val = parts[6]
            info_str = parts[7]
            
            # Parse INFO field
            info = parse_info_field(info_str)
            
            # Extract gene info
            gene = info.get("GENE", info.get("gene", None))
            star_allele = info.get("STAR", info.get("star", None))
            
            # Try to infer gene from rsID if not in INFO
            if not gene and rsid:
                gene = infer_gene_from_rsid(rsid)
            
            # Only keep target gene variants
            if gene and gene.upper() in TARGET_GENES:
                variant = {
                    "chrom": chrom,
                    "pos": pos,
                    "rsid": rsid or f"chr{chrom}:{pos}",
                    "ref": ref,
                    "alt": alt,
                    "gene": gene.upper(),
                    "star_allele": star_allele,
                    "quality": qual,
                    "filter": filter_val,
                    "info": info,
                    "functional_status": info.get("FUNCTION", info.get("function", "Unknown")),
                }
                variants.append(variant)
        
        except Exception as e:
            parse_errors.append(f"Line {line_num}: {str(e)}")
    
    detected_genes = list(set(v["gene"] for v in variants))
    
    # Group variants by gene
    variants_by_gene = {}
    for v in variants:
        gene = v["gene"]
        if gene not in variants_by_gene:
            variants_by_gene[gene] = []
        variants_by_gene[gene].append(v)
    
    return {
        "variants": variants,
        "variants_by_gene": variants_by_gene,
        "metadata": metadata,
        "detected_genes": detected_genes,
        "parse_errors": parse_errors,
        "total_variants": len(variants),
        "vcf_parsing_success": len(parse_errors) == 0 or len(variants) > 0,
    }


def infer_gene_from_rsid(rsid: str) -> Optional[str]:
    """Infer gene from known rsIDs for pharmacogenomic variants."""
    # Common pharmacogenomic rsID to gene mappings
    rsid_gene_map = {
        # CYP2D6
        "rs3892097": "CYP2D6",
        "rs5030655": "CYP2D6",
        "rs35742686": "CYP2D6",
        "rs1065852": "CYP2D6",
        "rs28371725": "CYP2D6",
        "rs16947": "CYP2D6",
        # CYP2C19
        "rs4244285": "CYP2C19",
        "rs4986893": "CYP2C19",
        "rs28399504": "CYP2C19",
        "rs56337013": "CYP2C19",
        "rs12248560": "CYP2C19",
        # CYP2C9
        "rs1799853": "CYP2C9",
        "rs1057910": "CYP2C9",
        "rs28371686": "CYP2C9",
        "rs9332131": "CYP2C9",
        # SLCO1B1
        "rs4149056": "SLCO1B1",
        "rs2306283": "SLCO1B1",
        "rs11045819": "SLCO1B1",
        # TPMT
        "rs1800460": "TPMT",
        "rs1142345": "TPMT",
        "rs1800462": "TPMT",
        "rs1800584": "TPMT",
        # DPYD
        "rs3918290": "DPYD",
        "rs55886062": "DPYD",
        "rs67376798": "DPYD",
        "rs75017182": "DPYD",
    }
    return rsid_gene_map.get(rsid)


def determine_diplotype(variants_for_gene: List[Dict]) -> str:
    """Determine diplotype from variants for a given gene."""
    star_alleles = [v["star_allele"] for v in variants_for_gene if v.get("star_allele")]
    
    if len(star_alleles) >= 2:
        return f"{star_alleles[0]}/{star_alleles[1]}"
    elif len(star_alleles) == 1:
        return f"{star_alleles[0]}/*1"
    else:
        return "*1/*1"  # Assume normal if no variants detected


def get_sample_vcf() -> str:
    """Return a sample VCF file for testing."""
    return """##fileformat=VCFv4.2
##fileDate=20240101
##source=PharmaGuardTest
##reference=GRCh38
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">
##INFO=<ID=STAR,Number=1,Type=String,Description="Star allele">
##INFO=<ID=FUNCTION,Number=1,Type=String,Description="Functional status">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr22\t42522613\trs4244285\tG\tA\t100\tPASS\tGENE=CYP2C19;STAR=*2;FUNCTION=no_function
chr22\t42523943\trs4986893\tG\tA\t100\tPASS\tGENE=CYP2C19;STAR=*3;FUNCTION=no_function
chr22\t42526694\trs1065852\tG\tA\t95\tPASS\tGENE=CYP2D6;STAR=*4;FUNCTION=no_function
chr12\t21331549\trs4149056\tT\tC\t98\tPASS\tGENE=SLCO1B1;STAR=*5;FUNCTION=decreased_function
chr6\t18143955\trs1800460\tC\tT\t99\tPASS\tGENE=TPMT;STAR=*3B;FUNCTION=no_function
chr1\t97915614\trs3918290\tC\tT\t97\tPASS\tGENE=DPYD;STAR=*2A;FUNCTION=no_function
"""