"""
LLM Explainer for PharmaGuard v1.0
Uses Groq API to generate clinical explanations for pharmacogenomic risk assessments
"""

import os
from typing import Dict, List, Optional
from groq import Groq

def get_groq_client(api_key: str) -> Groq:
    return Groq(api_key=api_key)


def build_clinical_prompt(drug: str, gene: str, diplotype: str, phenotype: str,
                           risk_label: str, severity: str, variants: List[Dict]) -> str:
    """Build a detailed clinical prompt for LLM explanation generation."""
    
    variant_details = ""
    if variants:
        for v in variants[:5]:  # limit to 5 variants
            variant_details += f"  - {v.get('rsid','N/A')} | {v.get('ref','?')}>{v.get('alt','?')} | Star: {v.get('star_allele','N/A')} | Function: {v.get('functional_status','Unknown')}\n"
    else:
        variant_details = "  - No variants detected (wild-type assumed)\n"
    
    prompt = f"""You are a clinical pharmacogenomics expert. Generate a concise, accurate clinical explanation for the following patient pharmacogenomic risk assessment.

PATIENT DATA:
- Drug: {drug}
- Gene: {gene}
- Diplotype: {diplotype}
- Phenotype: {phenotype}
- Risk Assessment: {risk_label} (Severity: {severity})
- Detected Variants:
{variant_details}

Generate a clinical explanation with these EXACT sections:

SUMMARY:
Write 2-3 sentences summarizing the overall risk and key clinical implication for this patient.

BIOLOGICAL_MECHANISM:
Explain in 2-3 sentences the biological mechanism: how the genetic variants affect the enzyme/protein, and how that impacts drug metabolism or response.

VARIANT_SIGNIFICANCE:
Explain the significance of the specific variants detected. Reference the rsIDs and star alleles. Keep to 2-3 sentences.

CLINICAL_IMPLICATIONS:
State the specific clinical implications for prescribing {drug} to this patient. What should the clinician do? 2-3 sentences.

Be precise, cite the specific variants (rsIDs), use correct pharmacological terminology. Do not add disclaimers or preambles."""
    
    return prompt


def parse_llm_response(response_text: str) -> Dict[str, str]:
    """Parse the structured LLM response into sections."""
    sections = {
        "summary": "",
        "biological_mechanism": "",
        "variant_significance": "",
        "clinical_implications": "",
    }
    
    section_map = {
        "SUMMARY:": "summary",
        "BIOLOGICAL_MECHANISM:": "biological_mechanism",
        "VARIANT_SIGNIFICANCE:": "variant_significance",
        "CLINICAL_IMPLICATIONS:": "clinical_implications",
    }
    
    current_section = None
    current_text = []
    
    for line in response_text.split("\n"):
        line = line.strip()
        matched = False
        for header, key in section_map.items():
            if line.startswith(header):
                if current_section:
                    sections[current_section] = " ".join(current_text).strip()
                current_section = key
                remainder = line[len(header):].strip()
                current_text = [remainder] if remainder else []
                matched = True
                break
        if not matched and current_section and line:
            current_text.append(line)
    
    if current_section:
        sections[current_section] = " ".join(current_text).strip()
    
    # Fallback: if parsing fails, put everything in summary
    if not any(sections.values()):
        sections["summary"] = response_text.strip()
    
    return sections


def generate_explanation(api_key: str, drug: str, gene: str, diplotype: str,
                          phenotype: str, risk_label: str, severity: str,
                          variants: List[Dict]) -> Dict[str, str]:
    """
    Generate clinical explanation using Groq LLM.
    Returns dict with summary, biological_mechanism, variant_significance, clinical_implications.
    """
    try:
        client = get_groq_client(api_key)
        prompt = build_clinical_prompt(drug, gene, diplotype, phenotype, risk_label, severity, variants)
        
        response = client.chat.completions.create(
            model="llama-3.3-70b-versatile",
            messages=[
                {
                    "role": "system",
                    "content": "You are a board-certified clinical pharmacologist and pharmacogenomics specialist. Provide accurate, evidence-based clinical explanations. Be concise and specific."
                },
                {
                    "role": "user",
                    "content": prompt
                }
            ],
            max_tokens=600,
            temperature=0.2,
        )
        
        raw_text = response.choices[0].message.content
        parsed = parse_llm_response(raw_text)
        parsed["model_used"] = "llama-3.3-70b-versatile"
        parsed["raw_response"] = raw_text
        parsed["success"] = True
        return parsed
    
    except Exception as e:
        return {
            "summary": f"LLM explanation unavailable: {str(e)}",
            "biological_mechanism": "",
            "variant_significance": "",
            "clinical_implications": "",
            "model_used": "llama-3.3-70b-versatile",
            "success": False,
            "error": str(e),
        }


def generate_all_explanations(api_key: str, risk_results: List[Dict]) -> List[Dict]:
    """Generate LLM explanations for all drug risk assessments."""
    enriched = []
    for result in risk_results:
        if result.get("error"):
            result["llm_explanation"] = {
                "summary": result["error"],
                "biological_mechanism": "",
                "variant_significance": "",
                "clinical_implications": "",
                "success": False,
            }
            enriched.append(result)
            continue
        
        explanation = generate_explanation(
            api_key=api_key,
            drug=result["drug"],
            gene=result["primary_gene"],
            diplotype=result["diplotype"],
            phenotype=result["phenotype"],
            risk_label=result["risk_label"],
            severity=result["severity"],
            variants=result.get("detected_variants", []),
        )
        result["llm_explanation"] = explanation
        enriched.append(result)
    
    return enriched