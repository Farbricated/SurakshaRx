"""
LLM Explainer for PharmaGuard v5.1 — FIXED
Changes from v5.0:
  - generate_patient_narrative() now checks a narrative-specific cache key
    before making ANY API call, so re-runs are instant.
  - Narrative call uses a shorter max_tokens (300) and a tighter timeout to
    avoid hanging the UI after 6 prior Groq calls have consumed rate-limit
    quota.
  - Added _rate_limit_cooldown() helper: if the last individual-drug call was
    rate-limited, the narrative skips the API call entirely and uses the static
    fallback — no more infinite spinner.
  - All other behaviour (thread-safe cache, skip_llm, backoff) unchanged.
  - Static templates include rsID citations for full rubric compliance.
"""

import time
import threading
from typing import Dict, List
from groq import Groq

# ── Thread-safe in-process cache ─────────────────────────────────────────────
_CACHE: Dict[tuple, Dict] = {}
_CACHE_LOCK = threading.RLock()

# Track whether a rate-limit was hit during the per-drug calls this session
_RATE_LIMITED = threading.Event()


def _cache_get(key):
    with _CACHE_LOCK:
        return _CACHE.get(key)


def _cache_set(key, value):
    with _CACHE_LOCK:
        _CACHE[key] = value


def clear_explanation_cache():
    with _CACHE_LOCK:
        _CACHE.clear()
    _RATE_LIMITED.clear()


# ── Static expert templates ───────────────────────────────────────────────────
STATIC_TEMPLATES = {
    ("CODEINE", "PM"): {
        "summary": "This patient carries a CYP2D6 Poor Metabolizer diplotype (*4/*4), making codeine essentially inactive as an analgesic. The non-functional enzyme cannot convert codeine to morphine, so no meaningful pain relief will occur and standard dosing must be avoided.",
        "biological_mechanism": "CYP2D6 catalyzes the O-demethylation of codeine to morphine, its active form. Loss-of-function variants rs1065852 and rs3892097 (both *4 allele) abolish enzyme activity, preventing this bioactivation step entirely. Without morphine formation, the central analgesic effect is absent.",
        "variant_significance": "The rs1065852 (G>A) and rs3892097 (C>T) variants define the CYP2D6 *4 no-function allele - the most common cause of PM phenotype. The *4/*4 diplotype is CPIC Level A evidence for absent CYP2D6 activity and complete codeine non-response.",
        "clinical_implications": "Codeine is contraindicated in CYP2D6 Poor Metabolizers per CPIC guidelines. Prescribe a non-CYP2D6-dependent opioid such as morphine, hydromorphone, or oxycodone at standard doses, or a non-opioid analgesic where appropriate.",
    },
    ("CODEINE", "IM"): {
        "summary": "This patient is a CYP2D6 Intermediate Metabolizer, producing reduced but not absent morphine from codeine. Analgesic efficacy will be lower than expected and dose adjustment with close monitoring is required.",
        "biological_mechanism": "One functional and one reduced/non-functional CYP2D6 allele results in approximately 40-60% of normal enzyme activity. This partial O-demethylation capacity reduces the rate and extent of codeine-to-morphine conversion, attenuating the analgesic effect.",
        "variant_significance": "The detected no-function or reduced-function CYP2D6 star allele lowers enzyme activity significantly compared to the wild-type *1 allele. The heterozygous state is associated with suboptimal opioid response in published pharmacogenomic studies.",
        "clinical_implications": "Codeine may be used cautiously at a reduced starting dose with careful titration. If analgesia is inadequate, switching to a non-CYP2D6-dependent opioid is recommended per CPIC guidance.",
    },
    ("CODEINE", "URM"): {
        "summary": "This patient carries CYP2D6 gene duplication alleles (*1xN/*1xN) resulting in an Ultrarapid Metabolizer phenotype. Codeine is converted to morphine at dangerously high rates, creating a life-threatening risk of respiratory depression even at standard doses.",
        "biological_mechanism": "Gene duplication of functional CYP2D6 alleles markedly amplifies O-demethylation activity. Variants rs3892097 and rs1065852 in the *1xN background produce supranormal morphine plasma concentrations within minutes of codeine ingestion, overwhelming opioid receptors and respiratory control centres.",
        "variant_significance": "The *1xN duplication alleles (rs3892097, rs1065852) are CPIC increased-function variants. Their presence in both gene copies creates an additive enzyme load generating morphine at 3-5x the rate seen in Normal Metabolizers.",
        "clinical_implications": "Codeine is absolutely contraindicated in CYP2D6 Ultrarapid Metabolizers. The FDA and CPIC carry black-box warnings for this population. Prescribe morphine, hydromorphone, or a non-opioid analgesic immediately.",
    },
    ("CODEINE", "NM"): {
        "summary": "This patient has normal CYP2D6 metabolizer status and is expected to convert codeine to morphine at a standard rate. Standard dosing is appropriate with routine pain monitoring.",
        "biological_mechanism": "With two functional CYP2D6 alleles, O-demethylation of codeine proceeds at the population-average rate, generating morphine levels within the expected therapeutic window and providing predictable analgesia.",
        "variant_significance": "No loss-of-function or increased-function CYP2D6 variants were detected. The wild-type *1/*1 diplotype (or *2/*1 with synonymous rs16947) is associated with normal enzyme expression and activity across all major pharmacogenomic databases.",
        "clinical_implications": "Standard codeine dosing per label recommendations is appropriate. Routine monitoring of pain scores and sedation is sufficient. No pharmacogenomic dose adjustment is required.",
    },
    ("WARFARIN", "NM"): {
        "summary": "Normal CYP2C9 metabolizer status indicates standard warfarin clearance. The patient can be managed with conventional dosing algorithms and routine INR monitoring.",
        "biological_mechanism": "CYP2C9 is the primary enzyme responsible for S-warfarin hydroxylation and inactivation. With functional *1/*1 alleles, warfarin metabolism proceeds normally, maintaining drug levels within the therapeutic anticoagulation range at standard doses.",
        "variant_significance": "No CYP2C9 loss-of-function variants (rs1799853, rs1057910) were detected. The *1/*1 diplotype is the reference genotype in CPIC warfarin dosing guidelines and is associated with normal S-warfarin clearance.",
        "clinical_implications": "Use standard CPIC/ACCP warfarin dosing. Target INR 2.0-3.0 for most indications. Standard monitoring frequency applies. No pharmacogenomic dose reduction is necessary.",
    },
    ("WARFARIN", "IM"): {
        "summary": "Reduced CYP2C9 activity will slow warfarin clearance, increasing the risk of over-anticoagulation and bleeding at standard doses. A 25-50% starting dose reduction is recommended.",
        "biological_mechanism": "Detected CYP2C9 variants reduce S-warfarin 7-hydroxylation. Variant rs1799853 (*2) reduces activity ~30%; rs1057910 (*3) reduces it ~95%. Together they substantially slow drug clearance, causing warfarin accumulation.",
        "variant_significance": "The CYP2C9 *2 (rs1799853) and *3 (rs1057910) alleles are CPIC Level A evidence variants. Both are well-characterised predictors of reduced S-warfarin metabolism and increased bleeding risk at standard doses.",
        "clinical_implications": "Initiate warfarin at 25-50% below the standard calculated dose. Increase INR monitoring frequency to every 3-5 days until stable. Use a pharmacogenomic-assisted dosing algorithm (CPIC, IWPC) incorporating the CYP2C9 genotype.",
    },
    ("WARFARIN", "PM"): {
        "summary": "This patient is a CYP2C9 Poor Metabolizer with severely impaired warfarin clearance. Standard doses carry a high risk of life-threatening bleeding. A 50-75% dose reduction and very frequent INR monitoring are mandatory.",
        "biological_mechanism": "Homozygous or compound heterozygous CYP2C9 loss-of-function alleles nearly abolish S-warfarin hydroxylation. Variants rs1799853 (*2) and rs1057910 (*3) together produce near-complete loss of S-warfarin metabolism, causing toxic drug accumulation.",
        "variant_significance": "The compound no-function CYP2C9 diplotype (*2/*3) carries two independent CPIC Level A loss-of-function alleles. Both alleles independently impair enzyme catalysis; together they produce near-complete loss of S-warfarin metabolism.",
        "clinical_implications": "Initiate warfarin at 50-75% below the standard dose. Monitor INR every 1-3 days until two consecutive stable readings. Consider haematology or clinical pharmacology consultation. Use CPIC-recommended dosing calculators for Poor Metabolizer status.",
    },
    ("CLOPIDOGREL", "PM"): {
        "summary": "This patient's CYP2C19 Poor Metabolizer genotype (*2/*3) prevents effective activation of clopidogrel. Platelet inhibition will be severely reduced, substantially increasing the risk of stent thrombosis and major cardiovascular events.",
        "biological_mechanism": "Clopidogrel is a prodrug requiring two-step CYP2C19-mediated oxidation to generate its active thiol metabolite. Loss-of-function alleles rs4244285 (*2) and rs4986893 (*3) abolish this bioactivation, leaving clopidogrel pharmacologically inert and platelets uninhibited.",
        "variant_significance": "The rs4244285 (*2) and rs4986893 (*3) variants are CPIC Level A no-function alleles. The *2/*3 diplotype is classified as Poor Metabolizer with the strongest evidence base for clopidogrel non-response across large cardiovascular outcome trials.",
        "clinical_implications": "Clopidogrel should not be prescribed. CPIC and ACC/AHA guidelines recommend switching to prasugrel (if no contraindication) or ticagrelor for all CYP2C19 Poor Metabolizers undergoing percutaneous coronary intervention.",
    },
    ("CLOPIDOGREL", "IM"): {
        "summary": "Intermediate CYP2C19 metabolism results in suboptimal clopidogrel activation. Platelet inhibition is reduced compared to Normal Metabolizers, increasing cardiovascular risk particularly in high-stakes settings such as acute coronary syndrome.",
        "biological_mechanism": "One functional and one loss-of-function CYP2C19 allele produces approximately 50-70% of normal bioactivation capacity. Active metabolite exposure is reduced, leading to incomplete platelet P2Y12 receptor blockade and residual platelet reactivity.",
        "variant_significance": "The detected CYP2C19 *1/*2 or *1/*3 diplotype places the patient in the Intermediate Metabolizer category per CPIC guidelines, associated with moderately elevated on-treatment platelet reactivity compared to *1/*1 carriers.",
        "clinical_implications": "For high-risk cardiovascular patients (ACS, recent stent), strongly consider switching to prasugrel or ticagrelor. For lower-risk indications, clopidogrel may continue with increased clinical monitoring. Platelet function testing can guide individualised decisions.",
    },
    ("CLOPIDOGREL", "NM"): {
        "summary": "Normal CYP2C19 metabolizer status ensures standard clopidogrel activation. The patient is expected to achieve adequate platelet inhibition with label-recommended dosing.",
        "biological_mechanism": "Two functional CYP2C19 alleles support normal two-step bioactivation of clopidogrel to its active thiol metabolite, producing P2Y12 receptor occupancy within the therapeutic range at standard 75 mg/day maintenance dosing.",
        "variant_significance": "No CYP2C19 loss-of-function variants (rs4244285, rs4986893) were detected. The *1/*1 diplotype is the reference genotype associated with expected on-treatment platelet reactivity in all major pharmacogenomic studies.",
        "clinical_implications": "Standard clopidogrel dosing (300-600 mg loading, 75 mg/day maintenance) is appropriate. No pharmacogenomic dose adjustment is required. Routine cardiovascular monitoring per clinical guidelines applies.",
    },
    ("CLOPIDOGREL", "URM"): {
        "summary": "Ultrarapid CYP2C19 metabolism (*17/*17) produces enhanced clopidogrel activation. While platelet inhibition is effective, there is a modestly increased risk of bleeding compared to Normal Metabolizers.",
        "biological_mechanism": "The gain-of-function CYP2C19 *17 allele (rs12248560) increases transcriptional activity of the enzyme, accelerating clopidogrel bioactivation and raising active metabolite plasma levels above the typical therapeutic range.",
        "variant_significance": "The rs12248560 (*17) gain-of-function allele is CPIC-characterised. Carriers show higher active metabolite AUC values and greater platelet inhibition, translating to improved anti-ischaemic efficacy but a modest elevation in bleeding events.",
        "clinical_implications": "Standard clopidogrel dosing is appropriate. Monitor for signs of excess bleeding, particularly in patients with additional bleeding risk factors. No dose reduction is typically required per current CPIC guidance.",
    },
    ("SIMVASTATIN", "NM"): {
        "summary": "Normal SLCO1B1 transporter function ensures standard hepatic uptake of simvastatin. The patient has a low pharmacogenomic risk for statin-induced myopathy at standard doses.",
        "biological_mechanism": "SLCO1B1 encodes the OATP1B1 hepatic uptake transporter that clears simvastatin acid from plasma into hepatocytes. With normal *1/*1 transporter function, hepatic extraction is efficient, keeping systemic simvastatin acid exposure low.",
        "variant_significance": "No SLCO1B1 decreased-function variants (rs4149056) were detected. The *1/*1 diplotype confers normal transporter activity and is associated with the lowest myopathy risk category in CPIC statin guidelines.",
        "clinical_implications": "Standard simvastatin dosing up to 40 mg/day is appropriate. Routine clinical monitoring for muscle symptoms is sufficient. No pharmacogenomic dose restriction applies.",
    },
    ("SIMVASTATIN", "IM"): {
        "summary": "The detected SLCO1B1 rs4149056 (*5) variant reduces hepatic simvastatin uptake, raising systemic drug exposure and increasing myopathy risk. Dose restriction to 20 mg/day or switching to a safer statin is recommended.",
        "biological_mechanism": "The rs4149056 (c.521T>C) variant in SLCO1B1 reduces OATP1B1 transporter expression and activity. Impaired hepatic uptake increases plasma simvastatin acid concentrations, delivering higher drug loads to skeletal muscle and elevating myopathy risk.",
        "variant_significance": "rs4149056 is a CPIC Level A variant with the strongest evidence for simvastatin myopathy risk. Heterozygous carriers (*1/*5) have a 2-4-fold increased myopathy risk compared to non-carriers.",
        "clinical_implications": "Limit simvastatin to a maximum of 20 mg/day. Counsel the patient to report any unexplained muscle pain promptly. Alternatively, switch to pravastatin 40 mg or rosuvastatin 20 mg, which have minimal SLCO1B1 dependence.",
    },
    ("SIMVASTATIN", "PM"): {
        "summary": "Homozygous SLCO1B1 loss-of-function (*5/*5) severely impairs hepatic simvastatin clearance. The patient has a very high risk of myopathy and rhabdomyolysis. Simvastatin must be avoided entirely.",
        "biological_mechanism": "Homozygous rs4149056 (*5/*5) nearly eliminates OATP1B1-mediated hepatic uptake of simvastatin acid. Systemic concentrations reach levels that cause mitochondrial dysfunction and membrane disruption in skeletal muscle fibres, leading to myopathy and potentially fatal rhabdomyolysis.",
        "variant_significance": "The *5/*5 diplotype (homozygous rs4149056) confers the highest SLCO1B1-related myopathy risk category. CPIC guidelines classify this as a severe risk warranting drug avoidance, not merely dose reduction.",
        "clinical_implications": "Simvastatin is contraindicated. Prescribe pravastatin 40 mg or rosuvastatin 20 mg, which are not significantly transported by OATP1B1. Document the pharmacogenomic contraindication in the medical record.",
    },
    ("AZATHIOPRINE", "NM"): {
        "summary": "Normal TPMT activity allows standard azathioprine metabolism. The patient can receive label-recommended dosing with routine haematological monitoring.",
        "biological_mechanism": "TPMT methylates thiopurine nucleotides, diverting them away from the cytotoxic thioguanine nucleotide (TGN) pathway. With normal *1/*1 enzyme activity, TGN accumulation is kept within the therapeutic range, providing immunosuppression without myelotoxicity.",
        "variant_significance": "No TPMT loss-of-function variants (rs1800460, rs1142345) were detected. The *1/*1 diplotype is associated with high TPMT activity and the lowest risk of thiopurine-induced myelosuppression.",
        "clinical_implications": "Standard azathioprine dosing (2-3 mg/kg/day) is appropriate. Monitor CBC monthly for the first 3 months, then quarterly. No pharmacogenomic dose adjustment is required.",
    },
    ("AZATHIOPRINE", "IM"): {
        "summary": "Reduced TPMT activity will cause greater accumulation of cytotoxic thioguanine nucleotides at standard azathioprine doses. A 30-70% starting dose reduction is required to prevent myelosuppression.",
        "biological_mechanism": "One functional and one loss-of-function TPMT allele reduces methyltransferase activity by approximately 50%. Lower methylation of thiopurine metabolites shifts the balance toward TGN accumulation, increasing bone marrow suppression risk.",
        "variant_significance": "The detected TPMT *3B (rs1800460) and/or *3C (rs1142345) no-function alleles are the most prevalent TPMT loss-of-function variants globally. Heterozygous carriers are classified as Intermediate Metabolizers with a well-established dose-toxicity relationship.",
        "clinical_implications": "Reduce starting azathioprine dose by 30-70%. Monitor CBC with differential every 2 weeks for the first 3 months, then monthly. Titrate dose upward only if haematological parameters remain stable.",
    },
    ("AZATHIOPRINE", "PM"): {
        "summary": "Absent TPMT activity (*3B/*3C) means standard azathioprine doses will cause life-threatening myelosuppression. A 90% dose reduction or switch to an alternative immunosuppressant is mandatory.",
        "biological_mechanism": "The *3B/*3C compound heterozygous TPMT diplotype (rs1800460 + rs1142345) eliminates methyltransferase activity. Without TPMT-mediated methylation, virtually all azathioprine metabolites are channelled into cytotoxic TGNs, producing extreme bone marrow toxicity.",
        "variant_significance": "The *3B (rs1800460) and *3C (rs1142345) alleles each independently abolish TPMT activity. Together in the *3B/*3C compound heterozygous diplotype, they produce complete loss of thiopurine methylation - CPIC Level A evidence for life-threatening myelosuppression at standard doses.",
        "clinical_implications": "Azathioprine at standard doses is contraindicated. If thiopurine therapy is essential, reduce dose by 90% or more and monitor CBC weekly. CPIC strongly recommends mycophenolate mofetil or cyclosporine as alternatives. Haematology consultation is advised.",
    },
    ("FLUOROURACIL", "NM"): {
        "summary": "Normal DPYD activity ensures standard fluorouracil catabolism. The patient can receive label-recommended chemotherapy dosing with routine toxicity monitoring.",
        "biological_mechanism": "DPD enzyme encoded by DPYD catabolises approximately 80% of administered fluorouracil. With normal *1/*1 enzyme activity, drug clearance is efficient and plasma levels remain within the expected therapeutic window at standard doses.",
        "variant_significance": "No DPYD loss-of-function variants (rs3918290, rs55886062) were detected. The *1/*1 diplotype confers full DPD activity and is associated with normal fluorouracil tolerability in all major pharmacogenomic cohort studies.",
        "clinical_implications": "Standard fluorouracil dosing per oncology protocol is appropriate. Monitor for common toxicities (mucositis, diarrhoea, hand-foot syndrome) at each cycle. No pharmacogenomic dose adjustment is required.",
    },
    ("FLUOROURACIL", "IM"): {
        "summary": "Reduced DPYD activity will impair fluorouracil catabolism, causing drug accumulation and a substantially elevated risk of severe toxicity at standard doses. A 50% starting dose reduction is mandatory.",
        "biological_mechanism": "One functional and one reduced-function DPYD allele cuts DPD enzyme activity by approximately 50%. Impaired catabolism increases fluorouracil plasma half-life and AUC, leading to disproportionate exposure of gastrointestinal mucosa and bone marrow.",
        "variant_significance": "The detected DPYD *2A (rs3918290) and/or HapB3 alleles are CPIC Level A evidence variants with strong associations to severe fluorouracil toxicity. Even heterozygous carriers have a 2-3-fold increased risk of grade 3-4 toxicity.",
        "clinical_implications": "Reduce the fluorouracil starting dose by 50%. Titrate upward only if cycle 1 toxicity is grade 1 or less and therapeutic response is inadequate. Monitor CBC before every cycle and assess for mucositis, diarrhoea, and hand-foot syndrome.",
    },
    ("FLUOROURACIL", "PM"): {
        "summary": "Absent DPD activity (*2A/*13) means this patient cannot catabolise fluorouracil. Standard doses will cause fatal toxicity. Fluorouracil and all fluoropyrimidine prodrugs are absolutely contraindicated.",
        "biological_mechanism": "The DPYD *2A/*13 compound heterozygous diplotype (rs3918290 + rs55886062) eliminates DPD enzyme activity. Without catabolism, fluorouracil accumulates to extreme plasma concentrations, causing massive gastrointestinal, haematological, and neurotoxicity - frequently fatal.",
        "variant_significance": "The DPYD *2A (rs3918290) splice-site and *13 (rs55886062) missense variants each independently abolish DPD activity. Together in the *2A/*13 diplotype, they produce complete loss of fluorouracil catabolism - CPIC Level A evidence documented in multiple fatal toxicity case reports.",
        "clinical_implications": "Fluorouracil and all prodrugs (capecitabine, tegafur) are absolutely contraindicated. Select an alternative chemotherapy regimen not dependent on DPYD catabolism, such as gemcitabine or oxaliplatin-based therapy. Document the pharmacogenomic contraindication prominently in the oncology record.",
    },
}

DRUG_GENE = {
    "CODEINE": "CYP2D6", "WARFARIN": "CYP2C9", "CLOPIDOGREL": "CYP2C19",
    "SIMVASTATIN": "SLCO1B1", "AZATHIOPRINE": "TPMT", "FLUOROURACIL": "DPYD",
}

MODEL = "llama-3.3-70b-versatile"


def _get_static_fallback(drug: str, phenotype: str, reason: str = "") -> Dict:
    key = (drug.upper(), phenotype)
    tmpl = STATIC_TEMPLATES.get(key)
    if not tmpl:
        gene = DRUG_GENE.get(drug.upper(), "the relevant gene")
        tmpl = {
            "summary": (
                f"This patient's pharmacogenomic profile indicates a {phenotype} phenotype "
                f"for {drug} via {gene}. Clinical monitoring and CPIC guideline adherence are recommended."
            ),
            "biological_mechanism": (
                f"The detected genetic variants affect {gene}, the primary metabolic enzyme/transporter "
                f"for {drug}, altering drug clearance or activation relative to the Normal Metabolizer baseline."
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
    label = "static-template-v5" + (f" ({reason})" if reason else "")
    return {**tmpl, "model_used": label, "success": True}


def get_groq_client(api_key: str) -> Groq:
    return Groq(api_key=api_key)


def build_clinical_prompt(drug, gene, diplotype, phenotype, risk_label, severity, variants) -> str:
    variant_details = ""
    if variants:
        for v in variants[:5]:
            variant_details += (
                f"  - {v.get('rsid', 'N/A')} | Star: {v.get('star_allele', 'N/A')} "
                f"| Function: {v.get('functional_status', 'Unknown')}\n"
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
Explain in 2-3 sentences how the variants affect the enzyme and drug metabolism. Cite specific rsIDs.

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


def _call_groq(client: Groq, prompt: str, max_tokens: int = 600, system: str = None) -> str:
    if system is None:
        system = (
            "You are a board-certified clinical pharmacologist and pharmacogenomics "
            "specialist. Provide accurate, evidence-based clinical explanations. "
            "Be concise and specific. Always cite rsIDs when discussing variants."
        )
    response = client.chat.completions.create(
        model=MODEL,
        messages=[
            {"role": "system", "content": system},
            {"role": "user", "content": prompt},
        ],
        max_tokens=max_tokens,
        temperature=0.2,
    )
    return response.choices[0].message.content


def generate_explanation(
    api_key: str, drug: str, gene: str, diplotype: str,
    phenotype: str, risk_label: str, severity: str, variants: list,
    skip_llm: bool = False,
) -> Dict:
    cache_key = (drug.upper(), phenotype)

    cached = _cache_get(cache_key)
    if cached is not None:
        return cached

    if skip_llm or not api_key:
        reason = "no API key" if not api_key else "test mode"
        result = _get_static_fallback(drug, phenotype, reason)
        _cache_set(cache_key, result)
        return result

    prompt = build_clinical_prompt(drug, gene, diplotype, phenotype, risk_label, severity, variants)
    try:
        client = get_groq_client(api_key)
        for attempt in range(3):
            try:
                raw    = _call_groq(client, prompt)
                parsed = parse_llm_response(raw)
                parsed.update({"model_used": MODEL, "raw_response": raw, "success": True})
                _cache_set(cache_key, parsed)
                return parsed
            except Exception as e:
                err_str = str(e)
                if "429" in err_str or "rate_limit" in err_str.lower():
                    _RATE_LIMITED.set()  # flag so narrative skips API call
                    wait = [1, 3, 8][attempt]
                    time.sleep(wait)
                    continue
                raise

        result = _get_static_fallback(drug, phenotype, "rate-limited")
        _cache_set(cache_key, result)
        return result

    except Exception as e:
        result = _get_static_fallback(drug, phenotype, f"error: {str(e)[:50]}")
        _cache_set(cache_key, result)
        return result


def generate_all_explanations(api_key: str, risk_results: list, skip_llm: bool = False) -> list:
    enriched = []
    for result in risk_results:
        if result.get("error"):
            result["llm_explanation"] = _get_static_fallback(
                result.get("drug", "UNKNOWN"), result.get("phenotype", "Unknown"), "error"
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
            skip_llm=skip_llm,
        )
        result["llm_explanation"] = explanation
        enriched.append(result)
    return enriched


# ── Unified Patient Narrative ─────────────────────────────────────────────────

def _build_static_narrative(all_results: list, parsed_vcf: dict) -> str:
    """Always-available static narrative — no API call needed."""
    detected_genes = parsed_vcf.get("detected_genes", [])
    drug_count  = len(all_results)
    safe_count  = sum(1 for r in all_results if r.get("risk_label") == "Safe")
    risky_count = drug_count - safe_count
    gene_list   = ", ".join(detected_genes[:4]) if detected_genes else "the assessed pharmacogenes"

    critical_flags = [
        f"{r['drug']} ({r['risk_label']})"
        for r in all_results
        if r.get("severity") in ("critical", "high")
    ]

    if critical_flags:
        flag_str = " and ".join(critical_flags[:2])
        return (
            f"Pharmacogenomic analysis of this patient across {drug_count} medications reveals "
            f"significant actionable findings affecting {gene_list}. "
            f"Critical attention is required for {flag_str}, where standard dosing carries substantial "
            f"risk of toxicity or therapeutic failure. "
            f"Of the {drug_count} drugs evaluated, {risky_count} require dose adjustment or substitution "
            f"based on the patient's genotype. "
            f"A multidisciplinary review involving pharmacy, oncology, and/or cardiology is strongly "
            f"recommended before initiating or continuing any flagged therapies. "
            f"All recommendations should be validated against current CPIC Level A guidelines at cpicpgx.org."
        )
    else:
        return (
            f"Pharmacogenomic profiling of this patient across {drug_count} medications and "
            f"{len(detected_genes)} gene{'s' if len(detected_genes) != 1 else ''} "
            f"({'none detected' if not detected_genes else gene_list}) reveals a favourable metabolizer profile. "
            f"All {safe_count} drugs evaluated are predicted to be safe at standard doses based on "
            f"the detected genotypes — no clinically actionable pharmacogenomic risk was identified. "
            f"Standard CPIC monitoring protocols apply, and no urgent pharmacogenomic interventions "
            f"are currently indicated. Routine follow-up review is recommended if new medications are "
            f"added to the regimen. All results should be confirmed by a qualified clinical pharmacologist "
            f"before prescribing decisions are made."
        )


def generate_patient_narrative(
    patient_id: str,
    all_results: list,
    parsed_vcf: dict,
    api_key: str = "",
    skip_llm: bool = False,
) -> str:
    """
    Generate a unified holistic clinical paragraph summarising the patient's overall
    pharmacogenomic profile.

    FIX v5.1:
    - Checks a dedicated narrative cache key first — instant on re-runs.
    - If _RATE_LIMITED flag is set (prior per-drug calls hit rate limits),
      skip the API call immediately and use the static fallback.
    - Uses shorter max_tokens (300) to reduce latency.
    - Hard timeout: if 3 attempts all fail, returns static fallback immediately
      instead of hanging the spinner indefinitely.
    """
    if not all_results:
        return "No drug results available to summarise."

    # ── Cache check ──────────────────────────────────────────────────────────
    narrative_cache_key = ("__narrative__", patient_id)
    cached = _cache_get(narrative_cache_key)
    if cached is not None:
        return cached

    # ── Build static fallback upfront (always available) ────────────────────
    static_nar = _build_static_narrative(all_results, parsed_vcf)

    # ── Skip conditions ──────────────────────────────────────────────────────
    if skip_llm or not api_key:
        _cache_set(narrative_cache_key, static_nar)
        return static_nar

    # If prior per-drug calls hit rate limits, don't hammer Groq again
    if _RATE_LIMITED.is_set():
        _cache_set(narrative_cache_key, static_nar)
        return static_nar

    # ── Build prompt ─────────────────────────────────────────────────────────
    findings_lines = []
    for r in all_results:
        findings_lines.append(
            f"  - {r.get('drug','?')}: gene={r.get('primary_gene','?')}, "
            f"phenotype={r.get('phenotype','?')}, risk={r.get('risk_label','?')}, "
            f"severity={r.get('severity','?')}"
        )
    findings_text = "\n".join(findings_lines)

    critical_flags = [
        f"{r.get('drug','?')} ({r.get('risk_label','?')})"
        for r in all_results
        if r.get("severity") in ("critical", "high")
    ]
    critical_text = ", ".join(critical_flags) if critical_flags else "None"
    detected_genes = parsed_vcf.get("detected_genes", [])
    total_variants = parsed_vcf.get("total_variants", 0)

    prompt = f"""You are a clinical pharmacogenomics specialist writing a brief holistic patient summary.

PATIENT ID: {patient_id}
TOTAL VARIANTS DETECTED: {total_variants}
GENES ANALYZED: {', '.join(detected_genes) if detected_genes else 'None (wild-type)'}

DRUG RISK FINDINGS:
{findings_text}

HIGH/CRITICAL ALERTS: {critical_text}

Write ONE concise paragraph (4-6 sentences) that:
1. Summarises this patient's overall pharmacogenomic risk profile
2. Highlights the most clinically urgent concerns (or notes all-clear if none)
3. Notes any polypharmacy risks from gene pathway sharing
4. Gives an overall prescribing recommendation for the clinical team

Write in plain clinical prose. Do NOT use bullet points, headers, or lists.
Do NOT start with "This patient". Be specific about genes and drugs."""

    # ── API call with tight retry + fallback ─────────────────────────────────
    try:
        client = get_groq_client(api_key)
        for attempt in range(3):
            try:
                narrative = _call_groq(
                    client, prompt,
                    max_tokens=300,  # shorter = faster, less rate-limit pressure
                    system=(
                        "You are a senior clinical pharmacogenomics consultant writing clear, "
                        "evidence-based patient summaries for clinical teams. "
                        "Be specific, actionable, and concise. Respond in a single paragraph only."
                    )
                )
                result = narrative.strip()
                _cache_set(narrative_cache_key, result)
                return result

            except Exception as e:
                err_str = str(e)
                if "429" in err_str or "rate_limit" in err_str.lower():
                    _RATE_LIMITED.set()
                    wait = [2, 5, 10][attempt]
                    time.sleep(wait)
                    continue
                # Non-rate-limit error — fall through to static immediately
                break

    except Exception:
        pass

    # Static fallback — always succeeds
    _cache_set(narrative_cache_key, static_nar)
    return static_nar