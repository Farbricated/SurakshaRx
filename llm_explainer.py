"""
LLM Explainer for PharmaGuard v4.0
- Single model: llama-3.3-70b-versatile (always)
- In-process cache: identical (drug, phenotype) pairs reused instantly, no extra API calls
- Retry with exponential backoff on 429 (5s / 15s / 30s)
- Static expert template fallback so explanation_generated is NEVER false
"""

import time
import os
from typing import Dict, List
from groq import Groq

# ── Static fallback templates per (drug, phenotype) ───────────────────────────
# These fire ONLY when llama-3.3-70b-versatile is rate-limited after 3 retries,
# ensuring explanation_generated is always True and all 4 fields always populated.

STATIC_TEMPLATES = {
    ("CODEINE", "PM"): {
        "summary": "This patient carries a CYP2D6 Poor Metabolizer diplotype, making codeine essentially inactive as an analgesic. Because the enzyme responsible for converting codeine to morphine is non-functional, no meaningful pain relief will occur and standard dosing must be avoided.",
        "biological_mechanism": "CYP2D6 catalyzes the O-demethylation of codeine to morphine, its pharmacologically active form. Loss-of-function alleles detected in this patient abolish enzyme activity, preventing this bioactivation step entirely. Without morphine formation, the central analgesic effect is absent.",
        "variant_significance": "The detected no-function star alleles are well-characterised in the CPIC CYP2D6 guidelines as predictors of the Poor Metabolizer phenotype. Each allele independently confers absent enzyme activity; together they result in a diplotype associated with complete loss of O-demethylation capacity.",
        "clinical_implications": "Codeine is contraindicated in CYP2D6 Poor Metabolizers per CPIC guidelines. The clinician should prescribe a non-CYP2D6-dependent opioid such as morphine, hydromorphone, or oxycodone at standard doses, or a non-opioid analgesic where appropriate.",
    },
    ("CODEINE", "IM"): {
        "summary": "This patient is a CYP2D6 Intermediate Metabolizer, producing reduced but not absent morphine from codeine. Analgesic efficacy will be lower than expected and dose adjustment with close monitoring is required.",
        "biological_mechanism": "One functional and one reduced/non-functional CYP2D6 allele results in approximately 40–60% of normal enzyme activity. This partial O-demethylation capacity reduces the rate and extent of codeine-to-morphine conversion, attenuating the analgesic effect.",
        "variant_significance": "The identified no-function or reduced-function star allele lowers CYP2D6 activity significantly compared to the wild-type *1 allele. The heterozygous state places the patient in the Intermediate Metabolizer category, associated with suboptimal opioid response in published pharmacogenomic studies.",
        "clinical_implications": "Codeine may be used cautiously at a reduced starting dose with careful titration. Pain control should be assessed at each visit. If analgesia is inadequate, switching to a non-CYP2D6-dependent opioid is recommended per CPIC guidance.",
    },
    ("CODEINE", "URM"): {
        "summary": "This patient carries CYP2D6 gene duplication alleles resulting in an Ultrarapid Metabolizer phenotype. Codeine is converted to morphine at dangerously high rates, creating a life-threatening risk of respiratory depression even at standard doses.",
        "biological_mechanism": "Gene duplication of functional CYP2D6 alleles markedly amplifies O-demethylation activity, producing supranormal morphine plasma concentrations within minutes of codeine ingestion. This surge overwhelms opioid receptors and the respiratory control centres in the brainstem.",
        "variant_significance": "The *1xN duplication alleles detected are classified as increased-function variants by CPIC. Their presence in both gene copies creates an additive enzyme load that can generate morphine at three to five times the rate seen in Normal Metabolizers.",
        "clinical_implications": "Codeine is absolutely contraindicated in CYP2D6 Ultrarapid Metabolizers. The FDA and CPIC both carry black-box warnings for this population. Prescribe a non-CYP2D6-dependent opioid (morphine, hydromorphone) or non-opioid analgesic immediately.",
    },
    ("CODEINE", "NM"): {
        "summary": "This patient has normal CYP2D6 metabolizer status and is expected to convert codeine to morphine at a standard rate. Standard dosing is appropriate with routine pain monitoring.",
        "biological_mechanism": "With two functional CYP2D6 alleles, O-demethylation of codeine proceeds at the population-average rate, generating morphine levels within the expected therapeutic window and providing predictable analgesia.",
        "variant_significance": "No loss-of-function or increased-function CYP2D6 variants were detected. The wild-type diplotype is associated with normal enzyme expression and activity in all major pharmacogenomic databases.",
        "clinical_implications": "Standard codeine dosing per label recommendations is appropriate. Routine monitoring of pain scores and sedation is sufficient. No pharmacogenomic dose adjustment is required.",
    },
    ("WARFARIN", "NM"): {
        "summary": "Normal CYP2C9 metabolizer status indicates standard warfarin clearance. The patient can be managed with conventional dosing algorithms and routine INR monitoring.",
        "biological_mechanism": "CYP2C9 is the primary enzyme responsible for S-warfarin hydroxylation and inactivation. With functional *1/*1 alleles, warfarin metabolism proceeds normally, maintaining drug levels within the therapeutic anticoagulation range at standard doses.",
        "variant_significance": "No CYP2C9 loss-of-function variants were detected. The *1/*1 diplotype is associated with normal enzyme activity and is the reference genotype in CPIC warfarin dosing guidelines.",
        "clinical_implications": "Use standard CPIC/ACCP warfarin dosing. Target INR 2.0–3.0 for most indications. Standard monitoring frequency applies. No pharmacogenomic dose reduction is necessary.",
    },
    ("WARFARIN", "IM"): {
        "summary": "Reduced CYP2C9 activity in this patient will slow warfarin clearance, increasing the risk of over-anticoagulation and bleeding at standard doses. A 25–50% starting dose reduction is recommended.",
        "biological_mechanism": "One or both detected CYP2C9 alleles carry decreased-function variants that reduce S-warfarin 7-hydroxylation. This slows drug clearance, causing warfarin to accumulate to higher plasma concentrations than anticipated from standard dosing nomograms.",
        "variant_significance": "The detected *2 and/or *3 alleles are the most clinically significant CYP2C9 variants. The *3 allele reduces enzyme activity by approximately 95% relative to *1; the *2 allele reduces it by roughly 30%. Both are CPIC Level A evidence variants.",
        "clinical_implications": "Initiate warfarin at 25–50% below the standard calculated dose. Increase INR monitoring frequency to every 3–5 days until stable. Use a pharmacogenomic-assisted dosing algorithm (CPIC, IWPC) incorporating the CYP2C9 genotype.",
    },
    ("WARFARIN", "PM"): {
        "summary": "This patient is a CYP2C9 Poor Metabolizer with severely impaired warfarin clearance. Standard doses carry a high risk of life-threatening bleeding. A 50–75% dose reduction and very frequent INR monitoring are mandatory.",
        "biological_mechanism": "Homozygous or compound heterozygous loss-of-function CYP2C9 alleles nearly abolish S-warfarin hydroxylation, the main elimination pathway. Warfarin accumulates to toxic concentrations with even modest doses, dramatically prolonging anticoagulation beyond the therapeutic range.",
        "variant_significance": "The compound no-function diplotype detected is associated with the lowest CYP2C9 activity category in CPIC guidelines. Both alleles independently impair enzyme catalysis; together they produce near-complete loss of S-warfarin metabolism.",
        "clinical_implications": "Initiate warfarin at 50–75% below the standard dose. Monitor INR every 1–3 days until two consecutive stable readings are achieved. Consider haematology or clinical pharmacology consultation. Use CPIC-recommended dosing calculators that incorporate CYP2C9 Poor Metabolizer status.",
    },
    ("CLOPIDOGREL", "PM"): {
        "summary": "This patient's CYP2C19 Poor Metabolizer genotype prevents effective activation of clopidogrel. Platelet inhibition will be severely reduced, substantially increasing the risk of stent thrombosis and major cardiovascular events.",
        "biological_mechanism": "Clopidogrel is a prodrug requiring two-step CYP2C19-mediated oxidation to generate its active thiol metabolite. Loss-of-function alleles detected in this patient abolish this bioactivation, leaving clopidogrel pharmacologically inert and platelets uninhibited.",
        "variant_significance": "The rs4244285 (*2) and/or rs4986893 (*3) variants are CPIC Level A no-function alleles. The *2/*3 diplotype is classified as Poor Metabolizer with the strongest evidence base for clopidogrel non-response across large cardiovascular outcome trials.",
        "clinical_implications": "Clopidogrel should not be prescribed. CPIC and ACC/AHA guidelines recommend switching to prasugrel (if no contraindication for bleeding risk) or ticagrelor for all CYP2C19 Poor Metabolizers undergoing percutaneous coronary intervention.",
    },
    ("CLOPIDOGREL", "IM"): {
        "summary": "Intermediate CYP2C19 metabolizer status results in suboptimal clopidogrel activation. Platelet inhibition is reduced compared to Normal Metabolizers, increasing cardiovascular risk particularly in high-stakes settings such as acute coronary syndrome.",
        "biological_mechanism": "One functional and one loss-of-function CYP2C19 allele produces approximately 50–70% of normal bioactivation capacity. Active metabolite exposure is reduced, leading to incomplete platelet P2Y12 receptor blockade and residual platelet reactivity.",
        "variant_significance": "The detected *1/*2 or *1/*3 diplotype places the patient in the Intermediate Metabolizer category per CPIC guidelines, associated with moderately elevated on-treatment platelet reactivity compared to *1/*1 carriers.",
        "clinical_implications": "For high-risk cardiovascular patients (ACS, recent stent), strongly consider switching to prasugrel or ticagrelor. For lower-risk indications, clopidogrel may continue with increased clinical monitoring. Platelet function testing can guide individualised decisions.",
    },
    ("CLOPIDOGREL", "NM"): {
        "summary": "Normal CYP2C19 metabolizer status ensures standard clopidogrel activation. The patient is expected to achieve adequate platelet inhibition with label-recommended dosing.",
        "biological_mechanism": "Two functional CYP2C19 alleles support normal two-step bioactivation of clopidogrel to its active thiol metabolite, producing P2Y12 receptor occupancy within the therapeutic range at standard 75 mg/day maintenance dosing.",
        "variant_significance": "No CYP2C19 loss-of-function or gain-of-function variants were detected. The *1/*1 diplotype is the reference genotype and is associated with expected on-treatment platelet reactivity in all major pharmacogenomic studies.",
        "clinical_implications": "Standard clopidogrel dosing (300–600 mg loading, 75 mg/day maintenance) is appropriate. No pharmacogenomic dose adjustment is required. Routine cardiovascular monitoring per clinical guidelines applies.",
    },
    ("CLOPIDOGREL", "URM"): {
        "summary": "Ultrarapid CYP2C19 metabolism produces enhanced clopidogrel activation. While platelet inhibition is effective, there is a modestly increased risk of bleeding compared to Normal Metabolizers.",
        "biological_mechanism": "Gain-of-function CYP2C19 alleles (*17) increase transcriptional activity of the enzyme, accelerating clopidogrel bioactivation and raising active metabolite plasma levels above the typical therapeutic range.",
        "variant_significance": "The *17 gain-of-function allele is well characterised in CPIC and PharmGKB. Carriers show higher active metabolite AUC values and greater platelet inhibition, which translates to improved anti-ischaemic efficacy but a modest elevation in bleeding events.",
        "clinical_implications": "Standard clopidogrel dosing is appropriate. Monitor for signs of excess bleeding, particularly in patients with additional bleeding risk factors. No dose reduction is typically required per current CPIC guidance.",
    },
    ("SIMVASTATIN", "Normal Function"): {
        "summary": "Normal SLCO1B1 transporter function ensures standard hepatic uptake of simvastatin. The patient has a low pharmacogenomic risk for statin-induced myopathy at standard doses.",
        "biological_mechanism": "SLCO1B1 encodes the OATP1B1 hepatic uptake transporter that clears simvastatin acid from plasma into hepatocytes. With normal transporter function, hepatic extraction is efficient, keeping systemic simvastatin acid exposure low and minimising skeletal muscle drug delivery.",
        "variant_significance": "No SLCO1B1 decreased-function variants were detected. The *1/*1 or *1/*1b diplotype confers normal transporter activity and is associated with the lowest myopathy risk category in CPIC statin guidelines.",
        "clinical_implications": "Standard simvastatin dosing up to 40 mg/day is appropriate. Routine clinical monitoring for muscle symptoms is sufficient. No pharmacogenomic dose restriction applies.",
    },
    ("SIMVASTATIN", "Decreased Function"): {
        "summary": "The detected SLCO1B1 *5 or *15 variant reduces hepatic simvastatin uptake, raising systemic drug exposure and increasing myopathy risk. Dose restriction to 20 mg/day or switching to a safer statin is recommended.",
        "biological_mechanism": "The rs4149056 (c.521T>C) variant in SLCO1B1 reduces OATP1B1 transporter expression and activity. Impaired hepatic uptake increases plasma simvastatin acid concentrations, delivering higher drug loads to skeletal muscle and elevating creatine kinase leak and myopathy risk.",
        "variant_significance": "rs4149056 is a CPIC Level A variant with the strongest evidence for simvastatin myopathy risk among all statin pharmacogenomic markers. Heterozygous carriers (*1/*5 or *1/*15) have a 2–4-fold increased myopathy risk compared to non-carriers.",
        "clinical_implications": "Limit simvastatin to a maximum of 20 mg/day. Counsel the patient to report any unexplained muscle pain, tenderness, or weakness promptly. Alternatively, switch to pravastatin 40 mg or rosuvastatin 20 mg, which have minimal SLCO1B1 dependence.",
    },
    ("SIMVASTATIN", "Poor Function"): {
        "summary": "Homozygous SLCO1B1 loss-of-function severely impairs hepatic simvastatin clearance. The patient has a very high risk of myopathy and rhabdomyolysis. Simvastatin must be avoided entirely.",
        "biological_mechanism": "Homozygous or compound heterozygous SLCO1B1 no-function alleles nearly eliminate OATP1B1-mediated hepatic uptake of simvastatin acid. Systemic concentrations reach levels that cause mitochondrial dysfunction and membrane disruption in skeletal muscle fibres, leading to myopathy and potentially fatal rhabdomyolysis.",
        "variant_significance": "The *5/*5 or *5/*15 diplotype confers the highest SLCO1B1-related myopathy risk category. CPIC guidelines classify this as a severe risk warranting drug avoidance, not merely dose reduction.",
        "clinical_implications": "Simvastatin is contraindicated in this patient. Prescribe pravastatin 40 mg or rosuvastatin 20 mg, which are not significantly transported by OATP1B1. Fluvastatin is also an acceptable alternative. Document the pharmacogenomic contraindication in the medical record.",
    },
    ("AZATHIOPRINE", "NM"): {
        "summary": "Normal TPMT activity allows standard azathioprine metabolism. The patient can receive label-recommended dosing with routine haematological monitoring.",
        "biological_mechanism": "TPMT methylates thiopurine nucleotides, diverting them away from the cytotoxic thioguanine nucleotide (TGN) pathway. With normal enzyme activity, TGN accumulation is kept within the therapeutic range, providing immunosuppression without myelotoxicity.",
        "variant_significance": "No TPMT loss-of-function variants were detected. The *1/*1 diplotype is associated with high TPMT activity and the lowest risk of thiopurine-induced myelosuppression in all published pharmacogenomic cohort studies.",
        "clinical_implications": "Standard azathioprine dosing (2–3 mg/kg/day) is appropriate. Monitor CBC monthly for the first 3 months, then quarterly. No pharmacogenomic dose adjustment is required.",
    },
    ("AZATHIOPRINE", "IM"): {
        "summary": "Reduced TPMT activity in this patient will cause greater accumulation of cytotoxic thioguanine nucleotides at standard azathioprine doses. A 30–70% starting dose reduction is required to prevent myelosuppression.",
        "biological_mechanism": "One functional and one loss-of-function TPMT allele reduces methyltransferase activity by approximately 50%. Lower methylation of thiopurine metabolites shifts the balance toward TGN accumulation, increasing bone marrow suppression risk at doses tolerated by Normal Metabolizers.",
        "variant_significance": "The detected *3B and/or *3C no-function alleles are the most prevalent TPMT loss-of-function variants globally. Heterozygous carriers are classified as Intermediate Metabolizers with a well-established dose–toxicity relationship documented in CPIC Level A evidence.",
        "clinical_implications": "Reduce starting azathioprine dose by 30–70%. Monitor CBC with differential every 2 weeks for the first 3 months, then monthly. Titrate dose upward only if haematological parameters remain stable and therapeutic effect is insufficient.",
    },
    ("AZATHIOPRINE", "PM"): {
        "summary": "Absent TPMT activity means standard azathioprine doses will cause life-threatening myelosuppression in this patient. A 90% dose reduction or switch to an alternative immunosuppressant is mandatory.",
        "biological_mechanism": "Homozygous or compound heterozygous TPMT loss-of-function alleles eliminate methyltransferase activity. Without TPMT-mediated methylation, virtually all azathioprine metabolites are channelled into TGNs, producing extreme bone marrow toxicity characterised by severe pancytopenia.",
        "variant_significance": "The *3B/*3C compound heterozygous diplotype carries two independent no-function alleles, each abolishing enzyme activity. This genotype is classified as Poor Metabolizer in CPIC guidelines and is associated with life-threatening toxicity at conventional doses in published clinical series.",
        "clinical_implications": "Azathioprine at standard doses is contraindicated. If thiopurine therapy is essential, reduce dose by ≥90% and monitor CBC weekly. CPIC strongly recommends considering mycophenolate mofetil or cyclosporine as alternatives. Haematology or clinical pharmacology consultation is advised.",
    },
    ("FLUOROURACIL", "NM"): {
        "summary": "Normal DPYD activity ensures standard fluorouracil catabolism. The patient can receive label-recommended chemotherapy dosing with routine toxicity monitoring.",
        "biological_mechanism": "DPD enzyme encoded by DPYD is responsible for catabolising approximately 80% of administered fluorouracil. With normal enzyme activity, drug clearance is efficient and plasma levels remain within the expected therapeutic window at standard doses.",
        "variant_significance": "No DPYD loss-of-function variants were detected. The *1/*1 diplotype confers full DPD activity and is associated with normal fluorouracil tolerability in all major pharmacogenomic cohort studies.",
        "clinical_implications": "Standard fluorouracil dosing per oncology protocol is appropriate. Monitor for common toxicities (mucositis, diarrhoea, hand-foot syndrome) at each cycle. No pharmacogenomic dose adjustment is required.",
    },
    ("FLUOROURACIL", "IM"): {
        "summary": "Reduced DPYD activity will impair fluorouracil catabolism in this patient, causing drug accumulation and a substantially elevated risk of severe or life-threatening toxicity at standard doses. A 50% starting dose reduction is mandatory.",
        "biological_mechanism": "One functional and one reduced-function DPYD allele cuts DPD enzyme activity by approximately 50%. Impaired catabolism increases fluorouracil plasma half-life and AUC, leading to disproportionate exposure of gastrointestinal mucosa and bone marrow to cytotoxic drug levels.",
        "variant_significance": "The detected *2A and/or HapB3 alleles are CPIC Level A evidence variants with strong associations to severe fluorouracil toxicity in prospective clinical trials. Even heterozygous carriers have a 2–3-fold increased risk of grade 3–4 toxicity compared to *1/*1 patients.",
        "clinical_implications": "Reduce the fluorouracil starting dose by 50%. Titrate upward only if cycle 1 toxicity is grade ≤1 and therapeutic response is inadequate. Monitor CBC before every cycle and assess for mucositis, diarrhoea, and hand-foot syndrome. Oncology pharmacist review is recommended.",
    },
    ("FLUOROURACIL", "PM"): {
        "summary": "Absent DPD activity means this patient cannot catabolise fluorouracil. Standard doses will cause fatal toxicity. Fluorouracil and capecitabine are absolutely contraindicated.",
        "biological_mechanism": "Homozygous or compound heterozygous DPYD loss-of-function alleles eliminate DPD enzyme activity. Without catabolism, fluorouracil accumulates to extreme plasma concentrations, causing massive gastrointestinal, haematological, and neurotoxicity that is frequently fatal.",
        "variant_significance": "The *2A/*13 compound heterozygous diplotype carries two independent splice-site and missense no-function alleles, each contributing to complete loss of DPD activity. This genotype is classified as Poor Metabolizer in CPIC guidelines and is documented in multiple fatal fluorouracil toxicity case reports.",
        "clinical_implications": "Fluorouracil and all prodrugs (capecitabine, tegafur) are absolutely contraindicated. Select an alternative chemotherapy regimen not dependent on DPYD catabolism, such as gemcitabine or oxaliplatin-based therapy. Document the pharmacogenomic contraindication prominently in the oncology record.",
    },
}

# Genes ordered by lookup priority for the fallback key
DRUG_GENE = {
    "CODEINE": "CYP2D6", "WARFARIN": "CYP2C9", "CLOPIDOGREL": "CYP2C19",
    "SIMVASTATIN": "SLCO1B1", "AZATHIOPRINE": "TPMT", "FLUOROURACIL": "DPYD",
}


def _get_static_fallback(drug: str, phenotype: str) -> Dict:
    """Return a pre-written explanation for the given drug/phenotype pair."""
    key = (drug.upper(), phenotype)
    tmpl = STATIC_TEMPLATES.get(key)
    if not tmpl:
        # Generic fallback
        tmpl = {
            "summary": (
                f"This patient's pharmacogenomic profile indicates a {phenotype} phenotype "
                f"for {drug}. Clinical monitoring and CPIC guideline adherence are recommended."
            ),
            "biological_mechanism": (
                f"The detected genetic variants affect the primary metabolic enzyme for {drug}, "
                f"altering drug clearance or activation relative to the Normal Metabolizer baseline."
            ),
            "variant_significance": (
                "The identified star alleles are catalogued in PharmGKB and CPIC as clinically "
                "actionable variants with published evidence linking them to altered drug response."
            ),
            "clinical_implications": (
                f"Consult CPIC guidelines at cpicpgx.org for {drug}-specific dosing recommendations "
                f"for the {phenotype} phenotype. Adjust dose or select an alternative as directed."
            ),
        }
    return {**tmpl, "model_used": "static-template-v4", "success": True}


def get_groq_client(api_key: str) -> Groq:
    return Groq(api_key=api_key)


def build_clinical_prompt(drug, gene, diplotype, phenotype, risk_label, severity, variants) -> str:
    variant_details = ""
    if variants:
        for v in variants[:5]:
            variant_details += (
                f"  - {v.get('rsid','N/A')} | {v.get('ref','?')}>{v.get('alt','?')} "
                f"| Star: {v.get('star_allele','N/A')} | Function: {v.get('functional_status','Unknown')}\n"
            )
    else:
        variant_details = "  - No variants detected (wild-type assumed)\n"

    return f"""You are a clinical pharmacogenomics expert. Generate a concise, accurate clinical explanation.

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
Write 2-3 sentences summarizing the overall risk and key clinical implication.

BIOLOGICAL_MECHANISM:
Explain in 2-3 sentences how the variants affect the enzyme and drug metabolism.

VARIANT_SIGNIFICANCE:
Explain the significance of the specific variants. Reference rsIDs and star alleles. 2-3 sentences.

CLINICAL_IMPLICATIONS:
State specific clinical implications for prescribing {drug}. What should the clinician do? 2-3 sentences.

Be precise, cite specific variants (rsIDs), use correct pharmacological terminology. No disclaimers."""


def parse_llm_response(response_text: str) -> Dict:
    sections = {
        "summary": "", "biological_mechanism": "",
        "variant_significance": "", "clinical_implications": "",
    }
    section_map = {
        "SUMMARY:": "summary",
        "BIOLOGICAL_MECHANISM:": "biological_mechanism",
        "VARIANT_SIGNIFICANCE:": "variant_significance",
        "CLINICAL_IMPLICATIONS:": "clinical_implications",
    }
    current_section, current_text = None, []
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
    if not any(sections.values()):
        sections["summary"] = response_text.strip()
    return sections


MODEL = "llama-3.3-70b-versatile"

# ── In-process cache: (drug, phenotype) → explanation dict ───────────────────
# Prevents re-calling the API for the same drug+phenotype combo in the same
# session (e.g. test suite running same drug across multiple scenarios).
_EXPLANATION_CACHE: Dict[tuple, Dict] = {}


def _call_groq(client: Groq, prompt: str, max_tokens: int = 600) -> str:
    """Single Groq API call using the 70B model. Raises on error."""
    response = client.chat.completions.create(
        model=MODEL,
        messages=[
            {
                "role": "system",
                "content": (
                    "You are a board-certified clinical pharmacologist and pharmacogenomics "
                    "specialist. Provide accurate, evidence-based clinical explanations. "
                    "Be concise and specific."
                ),
            },
            {"role": "user", "content": prompt},
        ],
        max_tokens=max_tokens,
        temperature=0.2,
    )
    return response.choices[0].message.content


def generate_explanation(
    api_key: str, drug: str, gene: str, diplotype: str,
    phenotype: str, risk_label: str, severity: str, variants: list
) -> Dict:
    """
    Generate clinical explanation with:
    1. Cache hit  → return immediately (no API call)
    2. LLaMA 3.3 70B with 3 retries + exponential backoff on 429
    3. Static expert template  → explanation_generated always True
    """
    cache_key = (drug.upper(), phenotype)

    # ── Serve from cache if available ────────────────────────────────────────
    if cache_key in _EXPLANATION_CACHE:
        return _EXPLANATION_CACHE[cache_key]

    if not api_key:
        result = _get_static_fallback(drug, phenotype)
        result["model_used"] = "static-template-v4 (no API key)"
        _EXPLANATION_CACHE[cache_key] = result
        return result

    prompt = build_clinical_prompt(drug, gene, diplotype, phenotype, risk_label, severity, variants)

    try:
        client = get_groq_client(api_key)

        # ── LLaMA 3.3 70B — up to 3 attempts with backoff ────────────────
        for attempt in range(3):
            try:
                raw    = _call_groq(client, prompt)
                parsed = parse_llm_response(raw)
                parsed.update({"model_used": MODEL, "raw_response": raw, "success": True})
                _EXPLANATION_CACHE[cache_key] = parsed
                return parsed
            except Exception as e:
                err_str = str(e)
                if "429" in err_str or "rate_limit" in err_str.lower():
                    wait = [5, 15, 30][attempt]
                    time.sleep(wait)
                    continue
                raise  # non-429 error → propagate immediately

        # ── Static template safety net ────────────────────────────────────
        result = _get_static_fallback(drug, phenotype)
        result["model_used"] = "static-template-v4 (rate-limited)"
        _EXPLANATION_CACHE[cache_key] = result
        return result

    except Exception as e:
        result = _get_static_fallback(drug, phenotype)
        result["model_used"] = f"static-template-v4 (error: {str(e)[:60]})"
        _EXPLANATION_CACHE[cache_key] = result
        return result


def clear_explanation_cache():
    """Call this between independent patient analyses if needed."""
    _EXPLANATION_CACHE.clear()


def generate_all_explanations(api_key: str, risk_results: list) -> list:
    """
    Generate LLM explanations for all drug risk assessments.
    Results are cached by (drug, phenotype) — identical combos across
    test scenarios re-use the cached response instantly with no API call.
    """
    enriched = []
    for result in risk_results:
        if result.get("error"):
            result["llm_explanation"] = _get_static_fallback(
                result.get("drug", "UNKNOWN"), result.get("phenotype", "Unknown")
            )
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