"""
Schema definitions for PharmaGuard v1.0
Pydantic models matching the exact required hackathon JSON output schema
"""

from pydantic import BaseModel, Field
from typing import List, Optional, Dict, Any
from datetime import datetime
import uuid


class DetectedVariant(BaseModel):
    rsid: str
    chrom: Optional[str] = None
    pos: Optional[str] = None
    ref: Optional[str] = None
    alt: Optional[str] = None
    gene: Optional[str] = None
    star_allele: Optional[str] = None
    functional_status: Optional[str] = None


class RiskAssessment(BaseModel):
    risk_label: str  # Safe | Adjust Dosage | Toxic | Ineffective | Unknown
    confidence_score: float
    severity: str  # none | low | moderate | high | critical


class PharmacogenomicProfile(BaseModel):
    primary_gene: str
    diplotype: str
    phenotype: str  # PM | IM | NM | RM | URM | Unknown
    detected_variants: List[DetectedVariant]


class ClinicalRecommendation(BaseModel):
    cpic_guideline: str
    dosing_recommendation: str
    alternative_drugs: Optional[List[str]] = []
    monitoring_required: Optional[str] = None
    contraindicated: bool = False


class LLMExplanation(BaseModel):
    summary: str
    biological_mechanism: str
    variant_significance: str
    clinical_implications: str
    model_used: Optional[str] = "llama-3.3-70b-versatile"


class QualityMetrics(BaseModel):
    vcf_parsing_success: bool
    variants_detected: int
    genes_analyzed: List[str]
    parse_errors: List[str] = []
    explanation_generated: bool = True


class PharmaGuardResult(BaseModel):
    patient_id: str
    drug: str
    timestamp: str
    risk_assessment: RiskAssessment
    pharmacogenomic_profile: PharmacogenomicProfile
    clinical_recommendation: ClinicalRecommendation
    llm_generated_explanation: LLMExplanation
    quality_metrics: QualityMetrics


# Alternative drug suggestions per risk scenario
ALTERNATIVE_DRUGS = {
    ("CODEINE", "PM"):  ["Morphine", "Hydromorphone", "Oxycodone", "Acetaminophen"],
    ("CODEINE", "URM"): ["Morphine", "Hydromorphone", "Non-opioid analgesics"],
    ("CLOPIDOGREL", "PM"):  ["Prasugrel", "Ticagrelor"],
    ("CLOPIDOGREL", "IM"):  ["Prasugrel", "Ticagrelor"],
    ("SIMVASTATIN", "Poor Function"): ["Pravastatin", "Rosuvastatin", "Fluvastatin"],
    ("SIMVASTATIN", "Decreased Function"): ["Pravastatin 40mg", "Rosuvastatin 20mg"],
    ("AZATHIOPRINE", "PM"): ["Mycophenolate mofetil", "Cyclosporine"],
    ("FLUOROURACIL", "PM"): ["Gemcitabine", "Oxaliplatin-based regimens"],
    ("FLUOROURACIL", "IM"): ["Reduced dose fluorouracil", "Capecitabine with dose reduction"],
}

MONITORING_RECOMMENDATIONS = {
    "WARFARIN": "INR monitoring every 3-5 days until stable, then monthly",
    "AZATHIOPRINE": "CBC with differential every 1-2 weeks for first 3 months",
    "SIMVASTATIN": "CK levels at baseline and if muscle symptoms develop",
    "FLUOROURACIL": "CBC before each cycle, monitor for mucositis, diarrhea, hand-foot syndrome",
    "CLOPIDOGREL": "Platelet function testing if available; monitor for ischemic events",
    "CODEINE": "Pain scores, respiratory rate, sedation levels",
}


def build_output_schema(patient_id: str, drug: str, result: Dict,
                         parsed_vcf: Dict, llm_exp: Dict) -> Dict:
    """Build the complete output JSON matching hackathon schema."""
    
    drug_upper = drug.upper()
    phenotype = result.get("phenotype", "Unknown")
    risk_label = result.get("risk_label", "Unknown")
    
    alternatives = ALTERNATIVE_DRUGS.get((drug_upper, phenotype), [])
    monitoring = MONITORING_RECOMMENDATIONS.get(drug_upper, "Standard clinical monitoring")
    contraindicated = risk_label in ("Toxic",) and result.get("severity") == "critical"
    
    # Build detected variants list
    raw_variants = result.get("detected_variants", [])
    detected_variants = []
    for v in raw_variants:
        detected_variants.append({
            "rsid": v.get("rsid", "unknown"),
            "chrom": v.get("chrom"),
            "pos": v.get("pos"),
            "ref": v.get("ref"),
            "alt": v.get("alt"),
            "gene": v.get("gene"),
            "star_allele": v.get("star_allele"),
            "functional_status": v.get("functional_status", "Unknown"),
        })
    
    output = {
        "patient_id": patient_id,
        "drug": drug_upper,
        "timestamp": datetime.utcnow().isoformat() + "Z",
        "risk_assessment": {
            "risk_label": risk_label,
            "confidence_score": result.get("confidence_score", 0.0),
            "severity": result.get("severity", "unknown"),
        },
        "pharmacogenomic_profile": {
            "primary_gene": result.get("primary_gene", "Unknown"),
            "diplotype": result.get("diplotype", "Unknown"),
            "phenotype": phenotype,
            "detected_variants": detected_variants,
        },
        "clinical_recommendation": {
            "cpic_guideline": f"CPIC Guideline for {drug_upper}",
            "dosing_recommendation": result.get("cpic_recommendation", "Consult CPIC guidelines."),
            "alternative_drugs": alternatives,
            "monitoring_required": monitoring,
            "contraindicated": contraindicated,
        },
        "llm_generated_explanation": {
            "summary": llm_exp.get("summary", ""),
            "biological_mechanism": llm_exp.get("biological_mechanism", ""),
            "variant_significance": llm_exp.get("variant_significance", ""),
            "clinical_implications": llm_exp.get("clinical_implications", ""),
            "model_used": llm_exp.get("model_used", "llama-3.3-70b-versatile"),
        },
        "quality_metrics": {
            "vcf_parsing_success": parsed_vcf.get("vcf_parsing_success", False),
            "variants_detected": len(raw_variants),
            "genes_analyzed": parsed_vcf.get("detected_genes", []),
            "parse_errors": parsed_vcf.get("parse_errors", []),
            "explanation_generated": llm_exp.get("success", False),
        },
    }
    
    return output