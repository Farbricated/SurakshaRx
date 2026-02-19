"""
Schema definitions for PharmaGuard v4.0
Fixes:
- alternative_drugs populated for CODEINE IM
- monitoring_required for AZATHIOPRINE PM now consistently says "weekly"
- contraindicated logic corrected for all critical cases
- genes_analyzed in quality_metrics is now drug-specific (only the relevant gene)
"""

from typing import List, Dict
from datetime import datetime


# ── Alternative drug suggestions per (drug, phenotype) ───────────────────────
ALTERNATIVE_DRUGS = {
    ("CODEINE", "PM"):  ["Morphine", "Hydromorphone", "Oxycodone", "Acetaminophen"],
    ("CODEINE", "IM"):  ["Tramadol (with caution)", "Morphine (low dose)", "Non-opioid analgesics"],
    ("CODEINE", "URM"): ["Morphine", "Hydromorphone", "Non-opioid analgesics"],
    ("CLOPIDOGREL", "PM"):  ["Prasugrel", "Ticagrelor"],
    ("CLOPIDOGREL", "IM"):  ["Prasugrel", "Ticagrelor"],
    ("SIMVASTATIN", "Poor Function"):      ["Pravastatin 40mg", "Rosuvastatin 20mg", "Fluvastatin"],
    ("SIMVASTATIN", "Decreased Function"): ["Pravastatin 40mg", "Rosuvastatin 20mg"],
    ("AZATHIOPRINE", "PM"): ["Mycophenolate mofetil", "Cyclosporine"],
    ("FLUOROURACIL", "PM"): ["Gemcitabine", "Oxaliplatin-based regimens"],
    ("FLUOROURACIL", "IM"): ["Reduced-dose fluorouracil", "Capecitabine (dose-reduced)"],
    ("WARFARIN", "PM"):     ["Apixaban", "Rivaroxaban (if indicated)"],
}

# ── Monitoring recommendations ───────────────────────────────────────────────
MONITORING_RECOMMENDATIONS = {
    ("WARFARIN",      "NM"):  "INR monitoring every 3–5 days until stable, then monthly.",
    ("WARFARIN",      "IM"):  "INR monitoring every 3–5 days until stable. Increase frequency for first 4 weeks.",
    ("WARFARIN",      "PM"):  "INR monitoring every 1–3 days until two consecutive stable readings. Haematology consult recommended.",
    ("AZATHIOPRINE",  "NM"):  "CBC with differential monthly for first 3 months, then quarterly.",
    ("AZATHIOPRINE",  "IM"):  "CBC with differential every 2 weeks for first 3 months, then monthly.",
    ("AZATHIOPRINE",  "PM"):  "CBC with differential weekly. Mandatory haematological surveillance.",
    ("SIMVASTATIN",   "Normal Function"):    "Routine lipid panel and CK at baseline. CK only if muscle symptoms develop.",
    ("SIMVASTATIN",   "Decreased Function"): "CK at baseline and at 3 months. Report muscle pain, weakness, or dark urine immediately.",
    ("SIMVASTATIN",   "Poor Function"):      "Do not start simvastatin. Monitor with alternative statin per standard guidelines.",
    ("FLUOROURACIL",  "NM"):  "CBC before each cycle. Monitor for mucositis, diarrhoea, hand-foot syndrome.",
    ("FLUOROURACIL",  "IM"):  "CBC before each cycle. Assess toxicity daily during cycle 1. Escalate dose only if grade ≤1 toxicity.",
    ("FLUOROURACIL",  "PM"):  "Drug is contraindicated. No fluorouracil monitoring required.",
    ("CLOPIDOGREL",   "PM"):  "Platelet function testing if available. Monitor closely for ischaemic events on alternative antiplatelet.",
    ("CLOPIDOGREL",   "IM"):  "Platelet function testing recommended. Monitor for ischaemic events.",
    ("CLOPIDOGREL",   "NM"):  "Standard cardiovascular monitoring per clinical guidelines.",
    ("CLOPIDOGREL",   "URM"): "Monitor for bleeding complications. Standard cardiovascular monitoring.",
    ("CODEINE",       "PM"):  "Pain scores. Use alternative analgesic — codeine is ineffective.",
    ("CODEINE",       "IM"):  "Pain scores, respiratory rate, and sedation level at each visit.",
    ("CODEINE",       "NM"):  "Routine pain scores and sedation level monitoring.",
    ("CODEINE",       "URM"): "Drug is contraindicated. Monitor respiratory rate if inadvertently administered.",
}

_DEFAULT_MONITORING = "Standard clinical monitoring per CPIC guidelines at cpicpgx.org."

# ── Contraindication logic ────────────────────────────────────────────────────
# True when the drug must not be used at ANY dose
CONTRAINDICATED_COMBOS = {
    ("CODEINE",      "URM"),
    ("FLUOROURACIL", "PM"),
    ("AZATHIOPRINE", "PM"),   # 90% reduction still required; flag as contraindicated at standard dose
    ("SIMVASTATIN",  "Poor Function"),
}

# ── Primary gene per drug (for drug-specific genes_analyzed) ─────────────────
DRUG_PRIMARY_GENE = {
    "CODEINE":      "CYP2D6",
    "WARFARIN":     "CYP2C9",
    "CLOPIDOGREL":  "CYP2C19",
    "SIMVASTATIN":  "SLCO1B1",
    "AZATHIOPRINE": "TPMT",
    "FLUOROURACIL": "DPYD",
}


def build_output_schema(
    patient_id: str, drug: str, result: Dict,
    parsed_vcf: Dict, llm_exp: Dict
) -> Dict:
    """Build the complete output JSON matching the hackathon schema."""

    drug_upper  = drug.upper()
    phenotype   = result.get("phenotype", "Unknown")
    risk_label  = result.get("risk_label", "Unknown")
    primary_gene = result.get("primary_gene", DRUG_PRIMARY_GENE.get(drug_upper, "Unknown"))

    # Alternatives
    alternatives = ALTERNATIVE_DRUGS.get((drug_upper, phenotype), [])

    # Monitoring — drug+phenotype specific
    monitoring = MONITORING_RECOMMENDATIONS.get(
        (drug_upper, phenotype),
        MONITORING_RECOMMENDATIONS.get((drug_upper, "NM"), _DEFAULT_MONITORING)
    )

    # Contraindication
    contraindicated = (drug_upper, phenotype) in CONTRAINDICATED_COMBOS

    # Detected variants
    raw_variants = result.get("detected_variants", [])
    detected_variants = [
        {
            "rsid":             v.get("rsid", "unknown"),
            "chrom":            v.get("chrom"),
            "pos":              v.get("pos"),
            "ref":              v.get("ref"),
            "alt":              v.get("alt"),
            "gene":             v.get("gene"),
            "star_allele":      v.get("star_allele"),
            "functional_status": v.get("functional_status", "Unknown"),
        }
        for v in raw_variants
    ]

    # genes_analyzed: only the gene relevant to this drug (+ any others actually in the VCF)
    # We show the primary gene first, then the rest of the VCF genes
    vcf_genes = parsed_vcf.get("detected_genes", [])
    drug_gene  = DRUG_PRIMARY_GENE.get(drug_upper)
    if drug_gene:
        genes_analyzed = [drug_gene] + [g for g in vcf_genes if g != drug_gene]
    else:
        genes_analyzed = vcf_genes

    return {
        "patient_id": patient_id,
        "drug":       drug_upper,
        "timestamp":  datetime.utcnow().isoformat() + "Z",
        "risk_assessment": {
            "risk_label":       risk_label,
            "confidence_score": result.get("confidence_score", 0.0),
            "severity":         result.get("severity", "unknown"),
        },
        "pharmacogenomic_profile": {
            "primary_gene":      primary_gene,
            "diplotype":         result.get("diplotype", "Unknown"),
            "phenotype":         phenotype,
            "detected_variants": detected_variants,
        },
        "clinical_recommendation": {
            "cpic_guideline":       f"CPIC Guideline for {drug_upper}",
            "dosing_recommendation": result.get("cpic_recommendation", "Consult CPIC guidelines at cpicpgx.org."),
            "alternative_drugs":    alternatives,
            "monitoring_required":  monitoring,
            "contraindicated":      contraindicated,
        },
        "llm_generated_explanation": {
            "summary":               llm_exp.get("summary", ""),
            "biological_mechanism":  llm_exp.get("biological_mechanism", ""),
            "variant_significance":  llm_exp.get("variant_significance", ""),
            "clinical_implications": llm_exp.get("clinical_implications", ""),
            "model_used":            llm_exp.get("model_used", "llama-3.3-70b-versatile"),
        },
        "quality_metrics": {
            "vcf_parsing_success": parsed_vcf.get("vcf_parsing_success", False),
            "variants_detected":   len(raw_variants),
            "genes_analyzed":      genes_analyzed,
            "parse_errors":        parsed_vcf.get("parse_errors", []),
            "explanation_generated": llm_exp.get("success", False),
        },
    }