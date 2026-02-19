"""
PharmaGuard — Pharmacogenomic Risk Prediction System
RIFT 2026 Hackathon | v4.0
Changes from v2:
  - LLM retry + fallback to llama-3.1-8b-instant + static template safety net
  - explanation_generated always True
  - alternative_drugs populated for all phenotypes (incl. CODEINE IM)
  - monitoring_required is now drug+phenotype specific (no more contradictions)
  - contraindicated logic corrected
  - genes_analyzed is drug-specific (primary gene first)
  - model_used shown in AI explanation header
"""

import streamlit as st
import json, uuid, os
from datetime import datetime
from dotenv import load_dotenv

load_dotenv()

from vcf_parser import parse_vcf, get_sample_vcf
from risk_engine import run_risk_assessment, get_overall_severity, DRUG_RISK_TABLE
from llm_explainer import generate_all_explanations
from schema import build_output_schema
from drug_interactions import run_interaction_analysis
from pdf_report import generate_pdf_report

BASE_DIR  = os.path.dirname(os.path.abspath(__file__))
ALL_DRUGS = list(DRUG_RISK_TABLE.keys())
GENE_DRUG_MAP = {
    "CODEINE": "CYP2D6", "WARFARIN": "CYP2C9", "CLOPIDOGREL": "CYP2C19",
    "SIMVASTATIN": "SLCO1B1", "AZATHIOPRINE": "TPMT", "FLUOROURACIL": "DPYD",
}
SEV_RANK = {"none": 0, "low": 1, "moderate": 2, "high": 3, "critical": 4}

RISK_CONFIG = {
    "Safe":          {"dot": "#22c55e", "text": "#16a34a", "bg": "#f0fdf4", "border": "#bbf7d0", "label": "Safe"},
    "Adjust Dosage": {"dot": "#f59e0b", "text": "#b45309", "bg": "#fffbeb", "border": "#fde68a", "label": "Adjust"},
    "Toxic":         {"dot": "#ef4444", "text": "#b91c1c", "bg": "#fef2f2", "border": "#fecaca", "label": "Toxic"},
    "Ineffective":   {"dot": "#8b5cf6", "text": "#7c3aed", "bg": "#f5f3ff", "border": "#ddd6fe", "label": "Ineffective"},
    "Unknown":       {"dot": "#94a3b8", "text": "#64748b", "bg": "#f8fafc", "border": "#e2e8f0", "label": "Unknown"},
}

st.set_page_config(page_title="PharmaGuard", page_icon="⬡", layout="wide", initial_sidebar_state="collapsed")

st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Instrument+Serif:ital@0;1&family=DM+Mono:wght@400;500&family=Geist:wght@300;400;500;600;700&display=swap');
*, *::before, *::after { box-sizing: border-box; margin: 0; padding: 0; }
html, body, [class*="css"] { font-family: 'Geist', -apple-system, sans-serif !important; background: #0a0a0a !important; color: #f0f0f0 !important; font-size: 16px !important; }
.stApp { background: #0a0a0a !important; }
.main .block-container { padding: 0 3rem 6rem !important; max-width: 1120px !important; }
#MainMenu, footer, header { visibility: hidden; }
.pg-nav { padding: 2rem 0 2.5rem; display: flex; align-items: center; justify-content: space-between; border-bottom: 1px solid #2a2a2a; margin-bottom: 3rem; }
.pg-wordmark { font-family: 'Instrument Serif', serif; font-size: 1.75rem; color: #f0f0f0; letter-spacing: -0.02em; }
.pg-wordmark em { font-style: italic; color: #6b7280; }
.pg-pill { font-family: 'DM Mono', monospace; font-size: 0.7rem; letter-spacing: 0.1em; text-transform: uppercase; color: #6b7280; background: #1a1a1a; border: 1px solid #2a2a2a; padding: 5px 12px; border-radius: 100px; }
.section-label { font-family: 'DM Mono', monospace; font-size: 0.7rem; letter-spacing: 0.12em; text-transform: uppercase; color: #6b7280; margin-bottom: 0.75rem; }
.stat-grid { display: grid; grid-template-columns: repeat(5, 1fr); gap: 1px; background: #1e1e1e; border: 1px solid #1e1e1e; border-radius: 12px; overflow: hidden; margin-bottom: 1.5rem; }
.stat-cell { background: #141414; padding: 1.35rem 1.5rem; }
.stat-val { font-family: 'Instrument Serif', serif; font-size: 2.25rem; color: #f0f0f0; line-height: 1; margin-bottom: 0.3rem; }
.stat-key { font-family: 'DM Mono', monospace; font-size: 0.67rem; letter-spacing: 0.1em; text-transform: uppercase; color: #525252; }
.stat-sub { font-size: 0.8rem; color: #4b4b4b; margin-top: 0.25rem; }
.rcard { border: 1px solid #222; border-radius: 12px; background: #111; margin-bottom: 1rem; overflow: hidden; }
.rcard-top { padding: 1.35rem 1.5rem; display: flex; align-items: center; justify-content: space-between; border-bottom: 1px solid #1e1e1e; }
.rcard-left { display: flex; align-items: center; gap: 0.875rem; }
.rcard-dot { width: 11px; height: 11px; border-radius: 50%; flex-shrink: 0; }
.rcard-name { font-family: 'Geist', sans-serif; font-size: 1.1rem; font-weight: 600; color: #f0f0f0; letter-spacing: -0.01em; }
.rcard-meta { font-family: 'DM Mono', monospace; font-size: 0.72rem; color: #4b4b4b; margin-top: 3px; letter-spacing: 0.04em; }
.rcard-badge { font-family: 'DM Mono', monospace; font-size: 0.7rem; font-weight: 500; letter-spacing: 0.08em; text-transform: uppercase; padding: 5px 14px; border-radius: 100px; border: 1px solid; }
.rcard-body { padding: 1.35rem 1.5rem; }
.mc-row { display: grid; grid-template-columns: repeat(4, 1fr); gap: 1px; background: #1a1a1a; border-radius: 8px; overflow: hidden; margin-bottom: 1.35rem; }
.mc-cell { background: #0e0e0e; padding: 1rem 1.1rem; }
.mc-key { font-family: 'DM Mono', monospace; font-size: 0.63rem; letter-spacing: 0.1em; text-transform: uppercase; color: #484848; margin-bottom: 0.35rem; }
.mc-val { font-family: 'Geist', sans-serif; font-size: 1.05rem; font-weight: 600; color: #e8e8e8; }
.conf-wrap { margin-bottom: 1.35rem; }
.conf-header { display: flex; justify-content: space-between; font-family: 'DM Mono', monospace; font-size: 0.68rem; letter-spacing: 0.08em; color: #484848; margin-bottom: 7px; }
.conf-track { height: 3px; background: #1e1e1e; border-radius: 2px; overflow: hidden; }
.conf-fill { height: 100%; border-radius: 2px; }
.h-rule { border: none; border-top: 1px solid #1e1e1e; margin: 1.35rem 0; }
.inline-label { font-family: 'DM Mono', monospace; font-size: 0.68rem; letter-spacing: 0.1em; text-transform: uppercase; color: #484848; margin-bottom: 0.65rem; }
.vtable { width: 100%; border-collapse: collapse; }
.vtable th { font-family: 'DM Mono', monospace; font-size: 0.66rem; letter-spacing: 0.09em; text-transform: uppercase; color: #484848; padding: 0 0.6rem 0.6rem; text-align: left; border-bottom: 1px solid #1e1e1e; }
.vtable td { font-family: 'DM Mono', monospace; font-size: 0.8rem; color: #a0a0a0; padding: 0.6rem 0.6rem; border-bottom: 1px solid #1a1a1a; }
.vtable tbody tr:last-child td { border-bottom: none; }
.v-rsid { color: #2563eb !important; } .v-star { color: #7c3aed !important; }
.v-nofunc { color: #dc2626 !important; } .v-dec { color: #d97706 !important; }
.v-inc { color: #2563eb !important; } .v-norm { color: #16a34a !important; }
.rec-box { border-radius: 8px; border: 1px solid; padding: 1.1rem 1.25rem; margin-bottom: 1rem; }
.rec-label { font-family: 'DM Mono', monospace; font-size: 0.66rem; letter-spacing: 0.09em; text-transform: uppercase; margin-bottom: 0.45rem; }
.rec-text { font-size: 0.975rem; line-height: 1.75; color: #b0b0b0; }
.alt-chips { display: flex; flex-wrap: wrap; gap: 0.45rem; }
.alt-chip { font-family: 'DM Mono', monospace; font-size: 0.72rem; color: #a0a0a0; background: #1a1a1a; border: 1px solid #2a2a2a; border-radius: 100px; padding: 4px 12px; }
.ai-block { border: 1px solid #222; border-radius: 8px; overflow: hidden; margin-top: 1.35rem; }
.ai-header { padding: 0.7rem 1.1rem; background: #0e0e0e; border-bottom: 1px solid #222; display: flex; align-items: center; gap: 0.65rem; }
.ai-badge { font-family: 'DM Mono', monospace; font-size: 0.63rem; letter-spacing: 0.09em; text-transform: uppercase; color: #6b7280; background: #1e1e1e; padding: 3px 9px; border-radius: 4px; }
.ai-badge-static { font-family: 'DM Mono', monospace; font-size: 0.63rem; letter-spacing: 0.09em; text-transform: uppercase; color: #9ca3af; background: #1a1a1a; border: 1px solid #2a2a2a; padding: 3px 9px; border-radius: 4px; }
.ai-section { padding: 1rem 1.1rem; border-bottom: 1px solid #1e1e1e; }
.ai-section:last-child { border-bottom: none; }
.ai-section-label { font-family: 'DM Mono', monospace; font-size: 0.66rem; letter-spacing: 0.09em; text-transform: uppercase; color: #484848; margin-bottom: 0.45rem; }
.ai-section-text { font-size: 0.975rem; line-height: 1.8; color: #b0b0b0; }
.alert { border-radius: 8px; border: 1px solid; border-left-width: 3px; padding: 1rem 1.1rem; margin-bottom: 1rem; display: flex; gap: 0.85rem; align-items: flex-start; }
.alert-label { font-family: 'DM Mono', monospace; font-size: 0.66rem; letter-spacing: 0.09em; text-transform: uppercase; margin-bottom: 0.3rem; }
.alert-text { font-size: 0.95rem; line-height: 1.7; }
.ix-row { display: flex; align-items: flex-start; gap: 0.875rem; padding: 1.1rem; border: 1px solid #1e1e1e; border-radius: 8px; margin-bottom: 0.5rem; background: #111; }
.ix-dot { width: 9px; height: 9px; border-radius: 50%; flex-shrink: 0; margin-top: 6px; }
.ix-title { font-size: 0.975rem; font-weight: 600; color: #e8e8e8; margin-bottom: 0.3rem; }
.ix-msg { font-size: 0.9rem; color: #525252; line-height: 1.65; margin-bottom: 0.45rem; }
.ix-rec { font-size: 0.875rem; color: #8a8a8a; line-height: 1.6; }
.ix-sev { font-family: 'DM Mono', monospace; font-size: 0.65rem; letter-spacing: 0.08em; text-transform: uppercase; padding: 3px 10px; border-radius: 4px; margin-left: auto; flex-shrink: 0; align-self: flex-start; }
.empty { text-align: center; padding: 5rem 2rem; border: 1px dashed #222; border-radius: 12px; background: #0e0e0e; }
.empty-icon { font-family: 'Instrument Serif', serif; font-size: 2.75rem; color: #d1d5db; margin-bottom: 1rem; }
.empty-title { font-family: 'Geist', sans-serif; font-size: 1.2rem; font-weight: 500; color: #333; margin-bottom: 0.5rem; }
.empty-hint { font-family: 'DM Mono', monospace; font-size: 0.7rem; color: #2a2a2a; letter-spacing: 0.06em; line-height: 2.1; }
.test-card { background: #111; border: 1px solid #1e1e1e; border-radius: 10px; padding: 1.2rem; }
.test-name { font-size: 0.95rem; font-weight: 600; color: #e0e0e0; margin-bottom: 0.4rem; }
.test-desc { font-family: 'DM Mono', monospace; font-size: 0.66rem; color: #3a3a3a; letter-spacing: 0.04em; margin-bottom: 0.85rem; line-height: 1.75; }
.test-row { display: flex; align-items: center; gap: 0.5rem; margin-bottom: 4px; }
.test-drug { font-family: 'DM Mono', monospace; font-size: 0.7rem; color: #3a3a3a; width: 90px; flex-shrink: 0; }
.test-result { font-family: 'DM Mono', monospace; font-size: 0.7rem; font-weight: 500; }
.rt-wrap { border: 1px solid #1e1e1e; border-radius: 10px; overflow: hidden; background: #0e0e0e; margin-bottom: 1rem; }
.rt-head { display: grid; grid-template-columns: 1fr 1.2fr 1.2fr 1.5fr 40px; background: #0a0a0a; border-bottom: 1px solid #1e1e1e; padding: 0 0.5rem; }
.rt-hcell { font-family: 'DM Mono', monospace; font-size: 0.66rem; letter-spacing: 0.09em; text-transform: uppercase; color: #484848; padding: 0.75rem 0.85rem; }
.rt-row { display: grid; grid-template-columns: 1fr 1.2fr 1.2fr 1.5fr 40px; border-bottom: 1px solid #141414; padding: 0 0.5rem; }
.rt-row:last-child { border-bottom: none; }
.rt-cell { font-family: 'DM Mono', monospace; font-size: 0.8rem; color: #909090; padding: 0.75rem 0.85rem; display: flex; align-items: center; }
.steps-row { display: flex; gap: 0; border: 1px solid #1e1e1e; border-radius: 10px; overflow: hidden; margin-bottom: 2.5rem; background: #0e0e0e; }
.step-item { flex: 1; padding: 1rem 1.35rem; border-right: 1px solid #1e1e1e; }
.step-item:last-child { border-right: none; }
.step-num { font-family: 'DM Mono', monospace; font-size: 0.63rem; letter-spacing: 0.1em; color: #d1d5db; margin-bottom: 4px; text-transform: uppercase; }
.step-label { font-size: 0.9rem; font-weight: 500; color: #2a2a2a; }
.step-item.active .step-num { color: #e0e0e0; }
.step-item.active .step-label { color: #e0e0e0; }
.stButton > button { background: #f0f0f0 !important; color: #111 !important; border: none !important; border-radius: 8px !important; font-family: 'Geist', sans-serif !important; font-weight: 500 !important; font-size: 0.975rem !important; padding: 0.7rem 1.75rem !important; letter-spacing: -0.01em !important; transition: opacity 0.15s !important; }
.stButton > button:hover { opacity: 0.8 !important; }
.stDownloadButton > button { background: #141414 !important; color: #909090 !important; border: 1px solid #2a2a2a !important; border-radius: 8px !important; font-family: 'DM Mono', monospace !important; font-size: 0.75rem !important; letter-spacing: 0.04em !important; padding: 0.55rem 1.1rem !important; }
.stDownloadButton > button:hover { border-color: #f0f0f0 !important; color: #f0f0f0 !important; }
.stTabs [data-baseweb="tab-list"] { background: transparent !important; border-bottom: 1px solid #1e1e1e !important; gap: 0 !important; padding: 0 !important; margin-bottom: 2rem !important; box-shadow: none !important; }
.stTabs [data-baseweb="tab"] { font-family: 'DM Mono', monospace !important; font-size: 0.72rem !important; letter-spacing: 0.1em !important; text-transform: uppercase !important; color: #3a3a3a !important; padding: 0.85rem 1.35rem !important; background: transparent !important; border: none !important; border-bottom: 2px solid transparent !important; border-radius: 0 !important; }
.stTabs [aria-selected="true"] { color: #f0f0f0 !important; border-bottom-color: #f0f0f0 !important; font-weight: 600 !important; }
.stTabs [data-baseweb="tab-panel"] { padding-top: 0 !important; }
div[data-testid="stExpander"] { background: #111 !important; border: 1px solid #1e1e1e !important; border-radius: 8px !important; box-shadow: none !important; margin-bottom: 0.5rem !important; }
div[data-testid="stExpander"] summary { font-family: 'DM Mono', monospace !important; font-size: 0.72rem !important; letter-spacing: 0.07em !important; color: #444 !important; padding: 0.85rem 1.1rem !important; }
.stMultiSelect span[data-baseweb="tag"] { background: #1e1e1e !important; color: #909090 !important; border: 1px solid #2a2a2a !important; font-family: 'DM Mono', monospace !important; font-size: 0.72rem !important; border-radius: 4px !important; }
.stCheckbox label p { font-family: 'Geist', sans-serif !important; font-size: 0.975rem !important; color: #909090 !important; }
.stTextInput > div > div > input { border-radius: 8px !important; border: 1px solid #e5e7eb !important; font-family: 'DM Mono', monospace !important; font-size: 0.875rem !important; }
.stFileUploader > div { border-radius: 8px !important; border: 1px dashed #e5e7eb !important; background: #0a0a0a !important; }
.stSuccess { border-radius: 8px !important; font-family: 'Geist', sans-serif !important; font-size: 0.975rem !important; border: 1px solid !important; background: #f0fdf4 !important; border-color: #bbf7d0 !important; color: #166534 !important; }
.stError { border-radius: 8px !important; background: #fef2f2 !important; border: 1px solid #fecaca !important; color: #991b1b !important; }
[data-testid="stSidebar"] { background: #0a0a0a !important; border-right: 1px solid #1e1e1e !important; }
[data-testid="stSidebar"] * { color: #606060 !important; }
</style>
""", unsafe_allow_html=True)


# ── Helpers ────────────────────────────────────────────────────
def load_vcf_file(filename):
    path = os.path.join(BASE_DIR, "sample_data", filename)
    return open(path).read() if os.path.exists(path) else get_sample_vcf()


def run_pipeline(vcf_content, drugs, pid, groq_key, run_ix=True, gen_pdf=True):
    parsed_vcf   = parse_vcf(vcf_content)
    risk_results = run_risk_assessment(parsed_vcf, drugs)
    # Always generate explanations — static templates fire if API is unavailable
    risk_results = generate_all_explanations(groq_key, risk_results)
    all_outputs  = [
        build_output_schema(patient_id=pid, drug=r["drug"], result=r,
                            parsed_vcf=parsed_vcf, llm_exp=r.get("llm_explanation", {}))
        for r in risk_results
    ]
    ix_report = run_interaction_analysis(drugs, risk_results) if run_ix and len(drugs) > 1 else None
    pdf_bytes = None
    if gen_pdf:
        try:
            pdf_bytes = generate_pdf_report(pid, all_outputs, parsed_vcf)
        except Exception as e:
            st.warning(f"PDF generation skipped: {e}")
    return parsed_vcf, risk_results, all_outputs, ix_report, pdf_bytes


def func_css(status):
    s = (status or "").lower()
    if "no_function"  in s: return "v-nofunc"
    if "decreased"    in s: return "v-dec"
    if "increased"    in s: return "v-inc"
    return "v-norm"


SEV_PALETTE = {
    "low":      {"dot": "#f59e0b", "bg": "#fffbeb", "border": "#fde68a", "text": "#b45309"},
    "moderate": {"dot": "#f97316", "bg": "#fff7ed", "border": "#fed7aa", "text": "#c2410c"},
    "high":     {"dot": "#ef4444", "bg": "#fef2f2", "border": "#fecaca", "text": "#b91c1c"},
    "critical": {"dot": "#dc2626", "bg": "#fef2f2", "border": "#fca5a5", "text": "#991b1b"},
}
IX_PALETTE = {"low": "#f59e0b", "moderate": "#f97316", "high": "#ef4444", "critical": "#dc2626"}


def render_results(all_outputs, parsed_vcf, ix_report, pdf_bytes, pid):
    overall_sev = max(
        (o["risk_assessment"]["severity"] for o in all_outputs),
        key=lambda s: SEV_RANK.get(s, 0), default="none"
    )
    sev_label = overall_sev.title() if overall_sev != "none" else "None"
    genes_str = " · ".join(parsed_vcf.get("detected_genes", [])) or "—"

    st.markdown(f"""
    <div class="stat-grid">
      <div class="stat-cell"><div class="stat-key">Patient</div>
        <div style="font-family:'DM Mono',monospace;font-size:0.85rem;color:#e0e0e0;margin-top:4px;">{pid}</div></div>
      <div class="stat-cell"><div class="stat-key">Overall Risk</div><div class="stat-val">{sev_label}</div></div>
      <div class="stat-cell"><div class="stat-key">Drugs</div><div class="stat-val">{len(all_outputs)}</div></div>
      <div class="stat-cell"><div class="stat-key">Variants</div><div class="stat-val">{parsed_vcf['total_variants']}</div></div>
      <div class="stat-cell"><div class="stat-key">Genes</div>
        <div class="stat-val">{len(parsed_vcf['detected_genes'])}<span style="font-size:1rem;color:#d1d5db;">/6</span></div>
        <div class="stat-sub">{genes_str}</div></div>
    </div>""", unsafe_allow_html=True)

    # Downloads
    dc1, dc2, dc3 = st.columns(3)
    with dc1:
        st.download_button("Download JSON", data=json.dumps(all_outputs, indent=2),
                           file_name=f"pharmaguard_{pid}_all.json", mime="application/json",
                           use_container_width=True, key=f"dlall_{pid}")
    with dc2:
        if pdf_bytes:
            st.download_button("Download PDF", data=pdf_bytes,
                               file_name=f"pharmaguard_{pid}.pdf", mime="application/pdf",
                               use_container_width=True, key=f"dlpdf_{pid}")
    with dc3:
        if ix_report and ix_report.get("interactions_found"):
            st.download_button("Download Interactions", data=json.dumps(ix_report, indent=2),
                               file_name=f"pharmaguard_{pid}_ix.json", mime="application/json",
                               use_container_width=True, key=f"dlix_{pid}")

    st.markdown("<div style='height:1.5rem'></div>", unsafe_allow_html=True)

    # Interactions panel
    if ix_report and ix_report["interactions_found"]:
        with st.expander(f"Drug–Drug Interactions  ·  {ix_report['total_interactions']} found", expanded=True):
            for ix in ix_report["all_interactions"]:
                sev   = ix.get("severity", "low")
                color = IX_PALETTE.get(sev, "#94a3b8")
                drugs_str = " + ".join(
                    ix.get("drugs_involved", []) or
                    [ix.get("inhibitor_drug", "")] + ix.get("affected_drugs", [])
                )
                st.markdown(f"""
                <div class="ix-row">
                  <div class="ix-dot" style="background:{color};margin-top:6px;"></div>
                  <div style="flex:1;">
                    <div class="ix-title">{drugs_str}</div>
                    <div class="ix-msg">{ix.get('message', ix.get('mechanism', ''))}</div>
                    <div class="ix-rec">💡 {ix.get('recommendation', '')}</div>
                  </div>
                  <span class="ix-sev" style="background:{color}18;color:{color};">{sev.upper()}</span>
                </div>""", unsafe_allow_html=True)
    elif ix_report:
        st.markdown("""
        <div style="padding:0.75rem 1rem;background:#f0fdf4;border:1px solid #bbf7d0;
             border-radius:8px;font-size:0.8rem;color:#16a34a;margin-bottom:1rem;
             font-family:'DM Mono',monospace;font-size:0.65rem;letter-spacing:0.04em;">
          ✓ No significant drug–drug interactions detected.</div>""", unsafe_allow_html=True)

    # Per-drug cards
    for output in all_outputs:
        risk_label  = output["risk_assessment"]["risk_label"]
        drug_name   = output["drug"]
        severity    = output["risk_assessment"]["severity"]
        confidence  = output["risk_assessment"]["confidence_score"]
        gene        = output["pharmacogenomic_profile"]["primary_gene"]
        diplotype   = output["pharmacogenomic_profile"]["diplotype"]
        phenotype   = output["pharmacogenomic_profile"]["phenotype"]
        variants    = output["pharmacogenomic_profile"]["detected_variants"]
        rec         = output["clinical_recommendation"]["dosing_recommendation"]
        alts        = output["clinical_recommendation"].get("alternative_drugs", [])
        monitoring  = output["clinical_recommendation"].get("monitoring_required", "")
        exp         = output["llm_generated_explanation"]
        rc          = RISK_CONFIG.get(risk_label, RISK_CONFIG["Unknown"])
        sp          = SEV_PALETTE.get(severity, {})

        if severity in ("critical", "high") and sp:
            st.markdown(f"""
            <div class="alert" style="background:{sp['bg']};border-color:{sp['border']};">
              <span style="font-size:0.9rem;flex-shrink:0;">{'⛔' if severity=='critical' else '⚠️'}</span>
              <div>
                <div class="alert-label" style="color:{sp['text']};">{severity.title()} Alert — {drug_name}</div>
                <div class="alert-text" style="color:{sp['text']};">{rec}</div>
              </div>
            </div>""", unsafe_allow_html=True)

        st.markdown(f"""
        <div class="rcard">
          <div class="rcard-top">
            <div class="rcard-left">
              <div class="rcard-dot" style="background:{rc['dot']};"></div>
              <div>
                <div class="rcard-name">{drug_name.title()}</div>
                <div class="rcard-meta">{gene} · {diplotype} · {phenotype}</div>
              </div>
            </div>
            <span class="rcard-badge" style="color:{rc['text']};border-color:{rc['border']};background:{rc['bg']};">{rc['label']}</span>
          </div>
          <div class="rcard-body">
            <div class="mc-row">
              <div class="mc-cell"><div class="mc-key">Phenotype</div><div class="mc-val" style="color:{rc['text']};">{phenotype}</div></div>
              <div class="mc-cell"><div class="mc-key">Severity</div><div class="mc-val">{severity.title()}</div></div>
              <div class="mc-cell"><div class="mc-key">Confidence</div><div class="mc-val">{confidence:.0%}</div></div>
              <div class="mc-cell"><div class="mc-key">Variants</div><div class="mc-val">{len(variants)}</div></div>
            </div>
            <div class="conf-wrap">
              <div class="conf-header">
                <span>Prediction confidence</span>
                <span style="color:{rc['text']};font-weight:500;">{confidence:.0%}</span>
              </div>
              <div class="conf-track"><div class="conf-fill" style="width:{confidence*100:.1f}%;background:{rc['dot']};"></div></div>
            </div>
        """, unsafe_allow_html=True)

        if variants:
            rows = ""
            for v in variants:
                fc = func_css(v.get("functional_status", ""))
                fn = (v.get("functional_status") or "unknown").replace("_", " ").title()
                rows += f"""<tr>
                  <td class="v-rsid">{v.get('rsid','—')}</td>
                  <td class="v-star">{v.get('star_allele','—')}</td>
                  <td>{v.get('ref','?')} → {v.get('alt','?')}</td>
                  <td>{v.get('chrom','—')}:{v.get('pos','—')}</td>
                  <td class="{fc}">{fn}</td>
                </tr>"""
            st.markdown(f"""
            <hr class="h-rule">
            <div class="inline-label">Detected Variants ({len(variants)})</div>
            <table class="vtable">
              <thead><tr><th>rsID</th><th>Star Allele</th><th>Change</th><th>Position</th><th>Function</th></tr></thead>
              <tbody>{rows}</tbody>
            </table>""", unsafe_allow_html=True)
        else:
            st.markdown("""<div style="font-family:'DM Mono',monospace;font-size:0.68rem;color:#9ca3af;
                 padding:0.5rem 0 0.25rem;letter-spacing:0.04em;">No variants detected — wild-type (*1/*1) assumed</div>""",
                unsafe_allow_html=True)

        st.markdown(f"""
        <hr class="h-rule">
        <div class="inline-label">CPIC Recommendation</div>
        <div class="rec-box" style="background:{rc['bg']};border-color:{rc['border']};">
          <div class="rec-label" style="color:{rc['text']};">CPIC Guideline for {drug_name}</div>
          <div class="rec-text">{rec}</div>
        </div>""", unsafe_allow_html=True)

        if alts:
            chips = "".join(f'<span class="alt-chip">{a}</span>' for a in alts)
            st.markdown(f"""<div class="inline-label">Alternative Drugs</div>
            <div class="alt-chips" style="margin-bottom:1rem;">{chips}</div>""", unsafe_allow_html=True)

        if monitoring:
            st.markdown(f"""
            <div style="display:flex;gap:0.75rem;align-items:flex-start;padding:0.875rem;
                 background:#fafafa;border:1px solid #e5e7eb;border-radius:8px;margin-bottom:1rem;">
              <span style="font-size:0.9rem;">🔬</span>
              <div>
                <div class="inline-label" style="margin-bottom:0.25rem;">Monitoring</div>
                <div style="font-size:0.85rem;color:#b0b0b0;line-height:1.65;">{monitoring}</div>
              </div>
            </div>""", unsafe_allow_html=True)

        if exp.get("summary"):
            model_used = exp.get("model_used", "")
            is_static  = "static" in model_used.lower()
            model_label = model_used if model_used else "llama-3.3-70b-versatile"
            badge_class = "ai-badge-static" if is_static else "ai-badge"

            blocks = ""
            for lbl, key in [("Summary","summary"),("Biological Mechanism","biological_mechanism"),
                              ("Variant Significance","variant_significance"),("Clinical Implications","clinical_implications")]:
                if exp.get(key):
                    blocks += f"""<div class="ai-section">
                      <div class="ai-section-label">{lbl}</div>
                      <div class="ai-section-text">{exp[key]}</div>
                    </div>"""
            st.markdown(f"""
            <div class="ai-block">
              <div class="ai-header">
                <span class="{badge_class}">{model_label}</span>
                <span style="font-family:'DM Mono',monospace;font-size:0.6rem;color:#9ca3af;letter-spacing:0.06em;">
                  AI Explanation · {drug_name}
                </span>
              </div>
              {blocks}
            </div>""", unsafe_allow_html=True)

        st.markdown("</div></div>", unsafe_allow_html=True)

        with st.expander(f"Raw JSON — {drug_name}", expanded=False):
            c1, _ = st.columns([1, 3])
            with c1:
                st.download_button("Download", data=json.dumps(output, indent=2),
                                   file_name=f"pharmaguard_{pid}_{drug_name}.json",
                                   mime="application/json", key=f"dl_{pid}_{drug_name}",
                                   use_container_width=True)
            st.code(json.dumps(output, indent=2), language="json")


# ══════════════════════════════════════════════════════════════
# NAV
st.markdown("""
<div class="pg-nav">
  <div class="pg-wordmark">Pharma<em>Guard</em></div>
  <div style="display:flex;gap:0.5rem;align-items:center;">
    <span class="pg-pill">CPIC Aligned</span>
    <span class="pg-pill">RIFT 2026</span>
    <span class="pg-pill">v4.0</span>
  </div>
</div>""", unsafe_allow_html=True)

# Sidebar
with st.sidebar:
    st.markdown("""<div style="padding:1rem 0 0.5rem;">
      <div style="font-family:'DM Mono',monospace;font-size:0.6rem;letter-spacing:0.12em;
           text-transform:uppercase;color:#9ca3af;margin-bottom:0.5rem;">Groq API Key</div>
    </div>""", unsafe_allow_html=True)
    groq_api_key = st.text_input("Groq API Key", value=os.environ.get("GROQ_API_KEY", ""),
                                  type="password", label_visibility="collapsed", placeholder="gsk_...")
    st.markdown("""<div style="font-family:'DM Mono',monospace;font-size:0.6rem;color:#9ca3af;
         padding-bottom:0.5rem;">console.groq.com</div>
    <div style="font-family:'DM Mono',monospace;font-size:0.58rem;color:#3a3a3a;
         padding-bottom:1rem;line-height:1.8;">
      Primary: LLaMA 3.3 70B<br>Fallback: LLaMA 3.1 8B<br>Safety net: static templates
    </div>
    <hr style="border:none;border-top:1px solid #1e1e1e;">
    <div style="padding:0.75rem 0 0.5rem;font-family:'DM Mono',monospace;font-size:0.6rem;
         letter-spacing:0.12em;text-transform:uppercase;color:#9ca3af;">Gene Map</div>""",
        unsafe_allow_html=True)
    for g, drug in [("CYP2D6","Codeine"),("CYP2C19","Clopidogrel"),("CYP2C9","Warfarin"),
                    ("SLCO1B1","Simvastatin"),("TPMT","Azathioprine"),("DPYD","Fluorouracil")]:
        st.markdown(f"""<div style="font-family:'DM Mono',monospace;font-size:0.68rem;
             padding:3px 0;color:#374151;">{g} <span style="color:#9ca3af;">→</span> {drug}</div>""",
            unsafe_allow_html=True)

tab1, tab2 = st.tabs(["Analysis", "Test Suite"])

# ── TAB 1 ─────────────────────────────────────────────────────
with tab1:
    st.markdown("""
    <div class="steps-row">
      <div class="step-item active"><div class="step-num">01</div><div class="step-label">Upload VCF</div></div>
      <div class="step-item active"><div class="step-num">02</div><div class="step-label">Select Drugs</div></div>
      <div class="step-item"><div class="step-num">03</div><div class="step-label">Run Analysis</div></div>
      <div class="step-item"><div class="step-num">04</div><div class="step-label">Review Results</div></div>
    </div>""", unsafe_allow_html=True)

    col_l, col_r = st.columns([1.3, 1], gap="large")

    with col_l:
        st.markdown('<div class="section-label">Genomic Data</div>', unsafe_allow_html=True)
        uploaded_file = st.file_uploader("Upload VCF file (max 5 MB)", type=["vcf"],
                                          help="VCF v4.2 with GENE, STAR, FUNCTION INFO tags")
        if uploaded_file:
            sz = uploaded_file.size / (1024 * 1024)
            if sz > 5:
                st.error(f"File too large: {sz:.1f} MB — max 5 MB")
                uploaded_file = None
            else:
                peek = uploaded_file.read(500).decode("utf-8", errors="replace")
                uploaded_file.seek(0)
                if "##fileformat=VCF" not in peek and "#CHROM" not in peek:
                    st.error("Invalid VCF file — must be VCF v4.x format")
                    uploaded_file = None
                else:
                    st.success(f"✓  {uploaded_file.name}  ·  {sz:.2f} MB")

        st.markdown('<div class="section-label" style="margin-top:1rem;">Or use a scenario</div>',
                    unsafe_allow_html=True)
        scenario_opts = {
            "None": None,
            "Mixed Variants (Standard)": "sample.vcf",
            "UltraRapid Metabolizer — Codeine Toxic": "test_ultrarapid_metabolizer.vcf",
            "All Normal Wild-type": "test_all_normal_wildtype.vcf",
            "Worst Case — All Poor Metabolizers": "test_worst_case_all_pm.vcf",
        }
        chosen_label = st.selectbox("Test Scenario", list(scenario_opts.keys()), label_visibility="collapsed")
        chosen_file  = scenario_opts[chosen_label]

    with col_r:
        st.markdown('<div class="section-label">Drugs</div>', unsafe_allow_html=True)
        drug_multiselect = st.multiselect(
            "Select drugs", options=ALL_DRUGS, default=["CLOPIDOGREL"],
            format_func=lambda x: f"{x.title()}  ({GENE_DRUG_MAP.get(x, '')})",
            label_visibility="collapsed"
        )
        st.markdown('<div class="section-label" style="margin-top:0.75rem;">Or type drug names</div>',
                    unsafe_allow_html=True)
        custom_drugs = st.text_input("Custom drugs", placeholder="CODEINE, WARFARIN, ...",
                                      label_visibility="collapsed")
        st.markdown('<div class="section-label" style="margin-top:0.75rem;">Patient ID</div>',
                    unsafe_allow_html=True)
        patient_id_input = st.text_input("Patient ID", placeholder="Auto-generated if blank",
                                          label_visibility="collapsed")
        c1, c2 = st.columns(2)
        with c1: run_interactions = st.checkbox("Check interactions", value=True)
        with c2: generate_pdf    = st.checkbox("Generate PDF", value=True)

    st.markdown("<div style='height:0.5rem'></div>", unsafe_allow_html=True)
    analyze_btn = st.button("Run Analysis →", use_container_width=True)

    if analyze_btn:
        all_drugs = list(drug_multiselect)
        if custom_drugs.strip():
            all_drugs += [d.strip().upper() for d in custom_drugs.split(",") if d.strip()]
        all_drugs = list(set(all_drugs))

        if not all_drugs:
            st.error("Select at least one drug.")
            st.stop()

        vcf_content = None
        if uploaded_file:
            vcf_content = uploaded_file.read().decode("utf-8", errors="replace")
        elif chosen_file:
            vcf_content = load_vcf_file(chosen_file)
        else:
            st.error("Upload a VCF file or select a test scenario.")
            st.stop()

        pid = patient_id_input.strip() or f"PG-{str(uuid.uuid4())[:8].upper()}"

        st.markdown(f"""
        <div style="display:flex;align-items:baseline;gap:1rem;margin:2.5rem 0 1.5rem;
             padding-bottom:1rem;border-bottom:1px solid #e5e7eb;">
          <div style="font-family:'Instrument Serif',serif;font-size:1.75rem;color:#f0f0f0;
               letter-spacing:-0.02em;">Results</div>
          <div style="font-family:'DM Mono',monospace;font-size:0.68rem;color:#9ca3af;">{pid}</div>
        </div>""", unsafe_allow_html=True)

        with st.spinner("Analysing…"):
            parsed_vcf, risk_results, all_outputs, ix_report, pdf_bytes = run_pipeline(
                vcf_content, all_drugs, pid, groq_api_key, run_interactions, generate_pdf
            )

        render_results(all_outputs, parsed_vcf, ix_report, pdf_bytes, pid)

        with st.expander("VCF Parse Details", expanded=False):
            p1, p2, p3 = st.columns(3)
            p1.metric("Total Variants", parsed_vcf["total_variants"])
            p2.metric("Genes Found", len(parsed_vcf["detected_genes"]))
            p3.metric("Parse Errors", len(parsed_vcf["parse_errors"]))
            for err in parsed_vcf.get("parse_errors", []):
                st.error(err)
    else:
        st.markdown("""
        <div class="empty">
          <div class="empty-icon">⬡</div>
          <div class="empty-title">Ready for analysis</div>
          <div class="empty-hint">
            Upload a VCF · select drugs · run<br>
            CYP2D6 · CYP2C19 · CYP2C9 · SLCO1B1 · TPMT · DPYD
          </div>
        </div>""", unsafe_allow_html=True)


# ── TAB 2: TEST SUITE ──────────────────────────────────────────
TEST_SUITE = [
    {"name":"Mixed Variants","file":"sample.vcf","icon":"⬡",
     "drugs":["CLOPIDOGREL","CODEINE","AZATHIOPRINE"],
     "expected":{"CLOPIDOGREL":"Ineffective","CODEINE":"Adjust Dosage","AZATHIOPRINE":"Adjust Dosage"},
     "desc":"CYP2C19 *2/*3 · CYP2D6 *4/*1 · TPMT *3B/*1"},
    {"name":"UltraRapid Metabolizer","file":"test_ultrarapid_metabolizer.vcf","icon":"⬡",
     "drugs":["CODEINE","CLOPIDOGREL"],
     "expected":{"CODEINE":"Toxic","CLOPIDOGREL":"Safe"},
     "desc":"CYP2D6 *1xN/*1xN → URM → Codeine Toxic"},
    {"name":"All Normal Wild-type","file":"test_all_normal_wildtype.vcf","icon":"⬡",
     "drugs":ALL_DRUGS,
     "expected":{d:"Safe" for d in ALL_DRUGS},
     "desc":"Wild-type *1/*1 across all 6 genes"},
    {"name":"Worst Case — All PM","file":"test_worst_case_all_pm.vcf","icon":"⬡",
     "drugs":ALL_DRUGS,
     "expected":{"CODEINE":"Ineffective","CLOPIDOGREL":"Ineffective","WARFARIN":"Adjust Dosage",
                 "SIMVASTATIN":"Toxic","AZATHIOPRINE":"Toxic","FLUOROURACIL":"Toxic"},
     "desc":"Loss-of-function alleles across all 6 genes"},
]

RISK_DOT = {"Safe":"#22c55e","Adjust Dosage":"#f59e0b","Toxic":"#ef4444","Ineffective":"#8b5cf6","Unknown":"#9ca3af"}

with tab2:
    st.markdown("""
    <div style="margin-bottom:2rem;">
      <div style="font-family:'Instrument Serif',serif;font-size:1.75rem;color:#f0f0f0;
           letter-spacing:-0.02em;margin-bottom:0.4rem;">Test Suite</div>
      <div style="font-family:'DM Mono',monospace;font-size:0.62rem;color:#9ca3af;letter-spacing:0.06em;">
        4 scenarios · expected vs actual · pass / fail per drug
      </div>
    </div>""", unsafe_allow_html=True)

    pc = st.columns(4)
    for i, sc in enumerate(TEST_SUITE):
        with pc[i]:
            rows_html = ""
            for d, r in list(sc["expected"].items())[:4]:
                color = RISK_DOT.get(r, "#9ca3af")
                rows_html += f"""<div class="test-row">
                  <span class="test-drug">{d[:9]}</span>
                  <span class="test-result" style="color:{color};">{r}</span>
                </div>"""
            st.markdown(f"""<div class="test-card">
              <div class="test-name">{sc['name']}</div>
              <div class="test-desc">{sc['desc']}</div>
              {rows_html}
            </div>""", unsafe_allow_html=True)

    st.markdown("<div style='height:1.25rem'></div>", unsafe_allow_html=True)
    tc1, tc2 = st.columns([3, 1])
    with tc1:
        use_llm = st.checkbox("Include LLM Explanations (~30s slower, uses Groq API)", value=False)
    with tc2:
        run_all = st.button("Run All 4 Tests →", use_container_width=True)

    if run_all:
        passed, failed = 0, 0
        suite_results  = []
        for sc in TEST_SUITE:
            vcf = load_vcf_file(sc["file"])
            pid = f"TEST-{sc['name'][:6].replace(' ','').upper()}"
            key = groq_api_key if use_llm else ""
            with st.spinner(f"Running: {sc['name']}…"):
                pv, _, ao, _, _ = run_pipeline(vcf, sc["drugs"], pid, key, run_ix=False, gen_pdf=False)
            rows, sc_pass = [], True
            for out in ao:
                drug  = out["drug"]
                got   = out["risk_assessment"]["risk_label"]
                exp   = sc["expected"].get(drug, "")
                ok    = (got == exp) if exp else True
                pheno = out["pharmacogenomic_profile"]["phenotype"]
                diplo = out["pharmacogenomic_profile"]["diplotype"]
                rows.append((drug, got, exp, ok, pheno, diplo))
                if not ok: sc_pass = False
            if sc_pass: passed += 1
            else:       failed += 1
            suite_results.append({"name": sc["name"], "pass": sc_pass, "rows": rows,
                                   "outputs": ao, "file": sc["file"]})

        total    = passed + failed
        ok_color = "#16a34a" if failed == 0 else "#b45309"
        ok_bg    = "#f0fdf4" if failed == 0 else "#fffbeb"
        ok_bdr   = "#bbf7d0" if failed == 0 else "#fde68a"
        st.markdown(f"""
        <div style="background:{ok_bg};border:1px solid {ok_bdr};border-radius:10px;
             padding:1.25rem 1.5rem;margin:1.25rem 0;display:flex;align-items:center;justify-content:space-between;">
          <div style="font-family:'Instrument Serif',serif;font-size:1.4rem;color:{ok_color};">
            {'All tests passed' if failed==0 else f'{passed}/{total} tests passed'}
          </div>
          <div style="font-family:'DM Mono',monospace;font-size:0.62rem;color:{ok_color};letter-spacing:0.06em;">
            {passed} passed · {failed} failed · {int(passed/total*100)}%
          </div>
        </div>""", unsafe_allow_html=True)

        for sr in suite_results:
            sym = "✓" if sr["pass"] else "✕"
            with st.expander(f"{sym}  {sr['name']}  —  {'PASS' if sr['pass'] else 'FAIL'}", expanded=not sr["pass"]):
                rows_html = ""
                for drug, got, exp, ok, pheno, diplo in sr["rows"]:
                    rc    = RISK_CONFIG.get(got, RISK_CONFIG["Unknown"])
                    ok_c  = "#16a34a" if ok else "#dc2626"
                    ok_bg2 = "#f0fdf4" if ok else "#fef2f2"
                    rows_html += f"""<div class="rt-row">
                      <div class="rt-cell" style="font-weight:600;color:#e0e0e0;">{drug}</div>
                      <div class="rt-cell"><span style="display:inline-flex;align-items:center;gap:6px;">
                        <span style="width:6px;height:6px;border-radius:50%;background:{rc['dot']};flex-shrink:0;"></span>{got}
                      </span></div>
                      <div class="rt-cell" style="color:#9ca3af;">{exp or '—'}</div>
                      <div class="rt-cell" style="color:#6b7280;">{diplo} / {pheno}</div>
                      <div class="rt-cell" style="justify-content:center;background:{ok_bg2};color:{ok_c};font-weight:700;">{sym}</div>
                    </div>"""
                st.markdown(f"""<div class="rt-wrap">
                  <div class="rt-head">
                    <div class="rt-hcell">Drug</div><div class="rt-hcell">Result</div>
                    <div class="rt-hcell">Expected</div><div class="rt-hcell">Diplotype / Phenotype</div>
                    <div class="rt-hcell"></div>
                  </div>{rows_html}</div>""", unsafe_allow_html=True)

                d1, d2 = st.columns(2)
                with d1:
                    st.download_button("Download JSON", data=json.dumps(sr["outputs"], indent=2),
                                       file_name=f"test_{sr['file'].replace('.vcf','')}.json",
                                       mime="application/json", key=f"tsc_{sr['name'][:14]}",
                                       use_container_width=True)
                with d2:
                    st.download_button("Download VCF", data=load_vcf_file(sr["file"]),
                                       file_name=sr["file"], mime="text/plain",
                                       key=f"vcf_{sr['name'][:14]}", use_container_width=True)

        st.download_button(
            "Download Full Test Suite JSON",
            data=json.dumps([{"scenario": s["name"], "pass": s["pass"], "results": s["outputs"]}
                             for s in suite_results], indent=2),
            file_name=f"pharmaguard_tests_{datetime.now().strftime('%Y%m%d_%H%M')}.json",
            mime="application/json", use_container_width=True
        )
    else:
        st.markdown("""
        <div class="empty">
          <div class="empty-icon">▷</div>
          <div class="empty-title">One-click validation</div>
          <div class="empty-hint">
            4 scenarios · expected vs actual · pass/fail per drug<br>
            Mixed · UltraRapid · All Normal · Worst Case
          </div>
        </div>""", unsafe_allow_html=True)