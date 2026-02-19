"""
Risk Engine for PharmaGuard v1.0
Maps gene diplotypes → phenotypes → drug-specific risks
Based on CPIC (Clinical Pharmacogenomics Implementation Consortium) guidelines
"""

from typing import Dict, List, Optional, Tuple

# ─────────────────────────────────────────────────────────────
# DIPLOTYPE → PHENOTYPE MAPPINGS (CPIC-based)
# PM = Poor Metabolizer, IM = Intermediate, NM = Normal
# RM = Rapid, URM = Ultrarapid, Unknown = insufficient data
# ─────────────────────────────────────────────────────────────

DIPLOTYPE_PHENOTYPE = {
    "CYP2D6": {
        "*1/*1": "NM",   "*1/*2": "NM",   "*2/*2": "NM",
        "*1/*4": "IM",   "*1/*5": "IM",   "*1/*6": "IM",
        "*4/*4": "PM",   "*4/*5": "PM",   "*5/*5": "PM",
        "*4/*6": "PM",   "*3/*4": "PM",   "*3/*5": "PM",
        "*1/*1xN": "URM", "*2xN/*2xN": "URM",
        "*10/*10": "IM", "*41/*41": "IM",
        "*1/*41": "IM",  "*2/*41": "IM",
        "*4/*41": "PM",
    },
    "CYP2C19": {
        "*1/*1": "NM",   "*1/*2": "IM",   "*1/*3": "IM",
        "*2/*2": "PM",   "*2/*3": "PM",   "*3/*3": "PM",
        "*1/*17": "RM",  "*17/*17": "URM",
        "*2/*17": "IM",  "*1/*4": "IM",
        "*4/*4": "PM",
    },
    "CYP2C9": {
        "*1/*1": "NM",   "*1/*2": "IM",   "*1/*3": "IM",
        "*2/*2": "IM",   "*2/*3": "PM",   "*3/*3": "PM",
        "*1/*5": "IM",   "*1/*6": "IM",
        "*5/*5": "PM",   "*6/*6": "PM",
    },
    "SLCO1B1": {
        "*1/*1": "Normal Function",
        "*1/*5": "Decreased Function",
        "*5/*5": "Poor Function",
        "*1/*15": "Decreased Function",
        "*15/*15": "Poor Function",
        "*1/*1b": "Normal Function",
    },
    "TPMT": {
        "*1/*1": "NM",   "*1/*2": "IM",   "*1/*3A": "IM",
        "*1/*3B": "IM",  "*1/*3C": "IM",
        "*2/*3A": "PM",  "*3A/*3A": "PM", "*3C/*3C": "PM",
        "*2/*2": "PM",
    },
    "DPYD": {
        "*1/*1": "NM",
        "*1/*2A": "IM",
        "*2A/*2A": "PM",
        "*1/*13": "IM",
        "*13/*13": "PM",
        "*1/HapB3": "IM",
        "HapB3/HapB3": "PM",
    },
}

# ─────────────────────────────────────────────────────────────
# DRUG × GENE × PHENOTYPE → RISK ASSESSMENT
# ─────────────────────────────────────────────────────────────

DRUG_RISK_TABLE = {
    "CODEINE": {
        "gene": "CYP2D6",
        "risks": {
            "PM":  ("Ineffective",   "moderate", 0.92, "CYP2D6 Poor Metabolizers cannot convert codeine to morphine. Pain relief will be absent. Consider alternative opioid."),
            "IM":  ("Adjust Dosage", "low",       0.85, "Reduced conversion to morphine. Lower efficacy expected. Monitor pain control closely."),
            "NM":  ("Safe",          "none",      0.95, "Normal codeine metabolism expected. Standard dosing appropriate."),
            "RM":  ("Adjust Dosage", "moderate",  0.88, "Higher-than-normal morphine levels possible. Monitor for side effects."),
            "URM": ("Toxic",         "critical",  0.97, "Ultra-rapid metabolizers convert codeine to morphine dangerously fast. Risk of respiratory depression and death. CONTRAINDICATED."),
            "Unknown": ("Unknown",   "low",       0.50, "Insufficient variant data. Genetic testing recommended before prescribing."),
        }
    },
    "WARFARIN": {
        "gene": "CYP2C9",
        "risks": {
            "NM":  ("Safe",          "none",      0.93, "Normal warfarin metabolism. Standard dosing per INR monitoring."),
            "IM":  ("Adjust Dosage", "moderate",  0.91, "Reduced warfarin metabolism. Start with 25-50% lower dose. Frequent INR monitoring required."),
            "PM":  ("Adjust Dosage", "high",      0.95, "Severely reduced metabolism. Very high bleeding risk at standard doses. Start with 50-75% dose reduction."),
            "Unknown": ("Unknown",   "low",       0.50, "Phenotype unknown. Use standard clinical monitoring protocols."),
        }
    },
    "CLOPIDOGREL": {
        "gene": "CYP2C19",
        "risks": {
            "PM":  ("Ineffective",   "high",      0.94, "Poor Metabolizers cannot activate clopidogrel. Platelet inhibition severely reduced. High risk of cardiovascular events. Use prasugrel or ticagrelor."),
            "IM":  ("Adjust Dosage", "moderate",  0.87, "Reduced activation. Suboptimal platelet inhibition. Consider alternative antiplatelet therapy."),
            "NM":  ("Safe",          "none",      0.95, "Normal activation of clopidogrel. Standard dosing appropriate."),
            "RM":  ("Safe",          "none",      0.90, "Slightly enhanced activation. Standard dosing appropriate."),
            "URM": ("Adjust Dosage", "low",       0.82, "Possibly enhanced activation. Monitor for bleeding risk."),
            "Unknown": ("Unknown",   "low",       0.50, "Phenotype unknown. Genetic testing recommended for high-risk patients."),
        }
    },
    "SIMVASTATIN": {
        "gene": "SLCO1B1",
        "risks": {
            "Normal Function":    ("Safe",          "none",     0.93, "Normal hepatic uptake of simvastatin. Standard dosing appropriate."),
            "Decreased Function": ("Adjust Dosage", "moderate", 0.90, "Reduced hepatic uptake leads to higher plasma simvastatin. Increased myopathy risk. Limit dose to 20mg or consider alternative statin."),
            "Poor Function":      ("Toxic",         "high",     0.95, "Severely impaired hepatic uptake. Very high risk of myopathy and rhabdomyolysis. Avoid simvastatin. Use pravastatin or rosuvastatin."),
            "Unknown":            ("Unknown",       "low",      0.50, "Phenotype unknown. Monitor for muscle pain and weakness."),
        }
    },
    "AZATHIOPRINE": {
        "gene": "TPMT",
        "risks": {
            "NM":  ("Safe",          "none",      0.94, "Normal TPMT activity. Standard azathioprine dosing appropriate. Routine monitoring."),
            "IM":  ("Adjust Dosage", "moderate",  0.91, "Reduced TPMT activity. Start at 30-70% of standard dose. Monitor for myelosuppression."),
            "PM":  ("Toxic",         "critical",  0.97, "Absent TPMT activity. Standard doses cause life-threatening myelosuppression. Reduce dose by 90% or use alternative immunosuppressant."),
            "Unknown": ("Unknown",   "low",       0.50, "Phenotype unknown. TPMT testing strongly recommended before initiating therapy."),
        }
    },
    "FLUOROURACIL": {
        "gene": "DPYD",
        "risks": {
            "NM":  ("Safe",          "none",      0.93, "Normal DPD activity. Standard fluorouracil dosing appropriate."),
            "IM":  ("Adjust Dosage", "high",      0.92, "Reduced DPD activity. Start at 50% dose reduction. Titrate based on toxicity. Life-threatening toxicity risk at standard doses."),
            "PM":  ("Toxic",         "critical",  0.98, "Absent DPD activity. Standard doses cause severe, potentially fatal toxicity. Fluorouracil is CONTRAINDICATED. Use alternative chemotherapy."),
            "Unknown": ("Unknown",   "low",       0.50, "Phenotype unknown. DPYD genotyping strongly recommended before fluorouracil-based chemotherapy."),
        }
    },
}

# CPIC dosing recommendations
CPIC_RECOMMENDATIONS = {
    ("CODEINE", "PM"):  "Avoid codeine. Use non-opioid analgesics or opioids not metabolized by CYP2D6 (e.g., morphine, hydromorphone).",
    ("CODEINE", "IM"):  "Use codeine with caution. Consider lower starting dose. Monitor for inadequate analgesia.",
    ("CODEINE", "NM"):  "Use label-recommended dosing.",
    ("CODEINE", "URM"): "Avoid codeine. Life-threatening respiratory depression risk. Use non-opioid or alternative opioid.",
    ("CODEINE", "RM"):  "Use label-recommended dosing. Monitor for opioid side effects.",
    
    ("WARFARIN", "NM"):  "Use standard ACCP/CPIC dosing algorithm. Target INR 2.0-3.0.",
    ("WARFARIN", "IM"):  "Initiate at 25-50% lower dose. Increase INR monitoring frequency. Use clinical decision support tools.",
    ("WARFARIN", "PM"):  "Initiate at 50-75% lower dose. Very frequent INR monitoring. Consider hematology consult.",
    
    ("CLOPIDOGREL", "PM"):  "Use alternative antiplatelet: prasugrel (if no contraindication) or ticagrelor.",
    ("CLOPIDOGREL", "IM"):  "Consider alternative antiplatelet therapy, especially for high-risk ACS/PCI patients.",
    ("CLOPIDOGREL", "NM"):  "Use label-recommended dosing (75mg/day maintenance).",
    ("CLOPIDOGREL", "URM"): "Use label-recommended dosing. Monitor for bleeding.",
    
    ("SIMVASTATIN", "Normal Function"):    "Use label-recommended dosing. Max 40mg/day.",
    ("SIMVASTATIN", "Decreased Function"): "Limit dose to 20mg/day. Consider switching to pravastatin, rosuvastatin, or fluvastatin.",
    ("SIMVASTATIN", "Poor Function"):      "Avoid simvastatin. Use pravastatin 40mg or rosuvastatin 20mg.",
    
    ("AZATHIOPRINE", "NM"):  "Use standard dosing (2-3 mg/kg/day). Monitor CBC monthly.",
    ("AZATHIOPRINE", "IM"):  "Reduce dose by 30-70%. Monitor CBC every 2 weeks for first 3 months.",
    ("AZATHIOPRINE", "PM"):  "Reduce dose by 90% or use alternative. Weekly CBC monitoring mandatory.",
    
    ("FLUOROURACIL", "NM"):  "Use label-recommended dosing.",
    ("FLUOROURACIL", "IM"):  "Reduce starting dose by 50%. Escalate based on tolerance. Monitor for severe toxicity.",
    ("FLUOROURACIL", "PM"):  "Avoid fluorouracil and capecitabine. Use alternative chemotherapy regimen.",
}


def get_phenotype(gene: str, diplotype: str) -> str:
    """Look up phenotype for a gene/diplotype combination."""
    gene_map = DIPLOTYPE_PHENOTYPE.get(gene, {})
    
    # Direct lookup
    if diplotype in gene_map:
        return gene_map[diplotype]
    
    # Try reversed diplotype
    parts = diplotype.split("/")
    if len(parts) == 2:
        reversed_dip = f"{parts[1]}/{parts[0]}"
        if reversed_dip in gene_map:
            return gene_map[reversed_dip]
    
    return "Unknown"


def assess_drug_risk(drug: str, gene: str, phenotype: str) -> Dict:
    """Get risk assessment for a drug/gene/phenotype combination."""
    drug_upper = drug.upper()
    
    if drug_upper not in DRUG_RISK_TABLE:
        return {
            "risk_label": "Unknown",
            "severity": "low",
            "confidence_score": 0.0,
            "clinical_note": f"Drug '{drug}' is not in the supported drug list.",
            "supported": False,
        }
    
    drug_info = DRUG_RISK_TABLE[drug_upper]
    risks = drug_info["risks"]
    
    risk_data = risks.get(phenotype, risks.get("Unknown", ("Unknown", "low", 0.0, "No data available.")))
    
    return {
        "risk_label": risk_data[0],
        "severity": risk_data[1],
        "confidence_score": risk_data[2],
        "clinical_note": risk_data[3],
        "primary_gene": drug_info["gene"],
        "supported": True,
    }


def get_cpic_recommendation(drug: str, phenotype: str) -> str:
    """Get CPIC dosing recommendation."""
    return CPIC_RECOMMENDATIONS.get(
        (drug.upper(), phenotype),
        "Consult CPIC guidelines at cpicpgx.org for specific dosing recommendations."
    )


def run_risk_assessment(parsed_vcf: Dict, drugs: List[str]) -> List[Dict]:
    """
    Run full risk assessment for all drugs against parsed VCF data.
    
    Returns list of complete risk assessment results per drug.
    """
    results = []
    variants_by_gene = parsed_vcf.get("variants_by_gene", {})
    
    for drug in drugs:
        drug_upper = drug.strip().upper()
        
        if drug_upper not in DRUG_RISK_TABLE:
            results.append({
                "drug": drug_upper,
                "error": f"Drug '{drug_upper}' not supported. Supported: {', '.join(DRUG_RISK_TABLE.keys())}",
                "risk_label": "Unknown",
                "severity": "low",
                "confidence_score": 0.0,
            })
            continue
        
        primary_gene = DRUG_RISK_TABLE[drug_upper]["gene"]
        gene_variants = variants_by_gene.get(primary_gene, [])
        
        # Determine diplotype
        if gene_variants:
            star_alleles = [v["star_allele"] for v in gene_variants if v.get("star_allele")]
            if len(star_alleles) >= 2:
                diplotype = f"{star_alleles[0]}/{star_alleles[1]}"
            elif len(star_alleles) == 1:
                diplotype = f"{star_alleles[0]}/*1"
            else:
                diplotype = "*1/*1"
        else:
            diplotype = "*1/*1"
        
        phenotype = get_phenotype(primary_gene, diplotype)
        risk = assess_drug_risk(drug_upper, primary_gene, phenotype)
        cpic_rec = get_cpic_recommendation(drug_upper, phenotype)
        
        results.append({
            "drug": drug_upper,
            "primary_gene": primary_gene,
            "diplotype": diplotype,
            "phenotype": phenotype,
            "risk_label": risk["risk_label"],
            "severity": risk["severity"],
            "confidence_score": risk["confidence_score"],
            "clinical_note": risk.get("clinical_note", ""),
            "cpic_recommendation": cpic_rec,
            "detected_variants": gene_variants,
        })
    
    return results


SEVERITY_ORDER = {"none": 0, "low": 1, "moderate": 2, "high": 3, "critical": 4}

def get_overall_severity(results: List[Dict]) -> str:
    """Get the highest severity across all drug assessments."""
    max_severity = "none"
    for r in results:
        s = r.get("severity", "none")
        if SEVERITY_ORDER.get(s, 0) > SEVERITY_ORDER.get(max_severity, 0):
            max_severity = s
    return max_severity