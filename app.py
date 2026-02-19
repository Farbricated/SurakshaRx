"""
PharmaGuard v2.0 — Pharmacogenomic Risk Prediction System
RIFT 2026 Hackathon | Pharmacogenomics / Explainable AI Track
Redesigned UI — Clinical Precision Dashboard
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

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

SAMPLE_VCF_MAP = {
    "Mixed Variants (Standard)":              "sample.vcf",
    "UltraRapid Metabolizer - Codeine TOXIC": "test_ultrarapid_metabolizer.vcf",
    "All Normal Wild-type - All Safe":         "test_all_normal_wildtype.vcf",
    "Worst Case - All Poor Metabolizers":      "test_worst_case_all_pm.vcf",
}

ALL_DRUGS = list(DRUG_RISK_TABLE.keys())

GENE_DRUG_MAP = {
    "CODEINE": "CYP2D6", "WARFARIN": "CYP2C9", "CLOPIDOGREL": "CYP2C19",
    "SIMVASTATIN": "SLCO1B1", "AZATHIOPRINE": "TPMT", "FLUOROURACIL": "DPYD"
}

SEV_RANK  = {"none": 0, "low": 1, "moderate": 2, "high": 3, "critical": 4}

RISK_CONFIG = {
    "Safe":         {"bg": "#022c22", "border": "#10b981", "text": "#34d399", "icon": "✓", "dot": "#10b981"},
    "Adjust Dosage":{"bg": "#1c1400", "border": "#f59e0b", "text": "#fbbf24", "icon": "⚠", "dot": "#f59e0b"},
    "Toxic":        {"bg": "#1f0000", "border": "#ef4444", "text": "#f87171", "icon": "✕", "dot": "#ef4444"},
    "Ineffective":  {"bg": "#0f0a2e", "border": "#8b5cf6", "text": "#a78bfa", "icon": "∅", "dot": "#8b5cf6"},
    "Unknown":      {"bg": "#111827", "border": "#4b5563", "text": "#9ca3af", "icon": "?", "dot": "#6b7280"},
}

SEV_CONFIG = {
    "none":     {"color": "#10b981", "label": "NO RISK"},
    "low":      {"color": "#f59e0b", "label": "LOW"},
    "moderate": {"color": "#f97316", "label": "MODERATE"},
    "high":     {"color": "#ef4444", "label": "HIGH"},
    "critical": {"color": "#dc2626", "label": "CRITICAL"},
}

# ─── Page Config ──────────────────────────────────────────────
st.set_page_config(
    page_title="PharmaGuard — Genomic Risk Intelligence",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="collapsed"
)

# ─── Global CSS ───────────────────────────────────────────────
st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500;600&family=Syne:wght@400;500;600;700;800&family=DM+Sans:ital,wght@0,300;0,400;0,500;0,600;1,400&display=swap');

html, body, [class*="css"] {
    font-family: 'DM Sans', sans-serif;
    background-color: #080c14 !important;
    color: #c8d6e5;
}
.stApp { background: #080c14; }
.main .block-container { padding: 0 2rem 4rem 2rem; max-width: 1380px; }
#MainMenu, footer, header { visibility: hidden; }
[data-testid="stSidebar"] { background: #060a10 !important; border-right: 1px solid #1a2332; }

.masthead {
    display: flex; align-items: center; justify-content: space-between;
    padding: 1.6rem 0 1.4rem 0;
    border-bottom: 1px solid #111e2e;
    margin-bottom: 2rem;
}
.masthead-left { display: flex; align-items: center; gap: 1.1rem; }
.masthead-logo {
    width: 44px; height: 44px;
    background: linear-gradient(135deg, #0ea5e9 0%, #6366f1 100%);
    border-radius: 11px; display: flex; align-items: center;
    justify-content: center; font-size: 1.4rem; flex-shrink: 0;
    box-shadow: 0 0 20px #0ea5e922;
}
.masthead-title {
    font-family: 'Syne', sans-serif;
    font-size: 1.5rem; font-weight: 800; color: #f0f6ff;
    letter-spacing: -0.03em; margin: 0; line-height: 1;
}
.masthead-sub {
    font-family: 'IBM Plex Mono', monospace;
    font-size: 0.62rem; color: #2a4060;
    letter-spacing: 0.09em; text-transform: uppercase; margin-top: 4px;
}
.masthead-tags { display: flex; gap: 0.5rem; align-items: center; }
.mtag {
    font-family: 'IBM Plex Mono', monospace; font-size: 0.58rem;
    font-weight: 600; padding: 3px 10px; border-radius: 4px;
    letter-spacing: 0.08em; text-transform: uppercase;
}
.mtag-blue  { background: #0c2040; color: #0ea5e9; border: 1px solid #0ea5e933; }
.mtag-green { background: #022c22; color: #10b981; border: 1px solid #10b98133; }
.mtag-violet{ background: #1e1040; color: #8b5cf6; border: 1px solid #8b5cf633; }

.step-row {
    display: flex; margin-bottom: 2rem;
    border: 1px solid #111e2e; border-radius: 10px; overflow: hidden;
}
.step-item {
    flex: 1; padding: 0.85rem 1.2rem;
    border-right: 1px solid #111e2e; background: #080c14;
}
.step-item:last-child { border-right: none; }
.step-num {
    font-family: 'IBM Plex Mono', monospace; font-size: 0.58rem;
    color: #1a2e44; letter-spacing: 0.1em; margin-bottom: 0.2rem;
}
.step-name {
    font-family: 'Syne', sans-serif; font-size: 0.8rem;
    font-weight: 700; color: #2a4060;
}
.step-item.on { background: #0a1828; }
.step-item.on .step-num { color: #0ea5e9; }
.step-item.on .step-name { color: #7dd3fc; }

.summary-bar {
    display: flex; border-radius: 12px;
    overflow: hidden; border: 1px solid #111e2e; margin-bottom: 1.5rem;
}
.sc {
    flex: 1; background: #0a1220; padding: 1.1rem 1.4rem;
    border-right: 1px solid #111e2e;
}
.sc:last-child { border-right: none; }
.sc-lbl {
    font-family: 'IBM Plex Mono', monospace; font-size: 0.58rem;
    color: #1e3050; letter-spacing: 0.12em; text-transform: uppercase; margin-bottom: 0.35rem;
}
.sc-val {
    font-family: 'Syne', sans-serif; font-size: 1.35rem;
    font-weight: 800; color: #ddeeff; line-height: 1;
}
.sc-sub {
    font-family: 'IBM Plex Mono', monospace; font-size: 0.58rem;
    color: #1e3050; margin-top: 0.2rem;
}

.risk-card {
    border-radius: 12px; margin-bottom: 1.1rem;
    overflow: hidden; border: 1px solid;
}
.rc-head {
    display: flex; align-items: center;
    padding: 1.1rem 1.4rem; gap: 1rem;
}
.rc-icon {
    width: 38px; height: 38px; border-radius: 8px;
    display: flex; align-items: center; justify-content: center;
    font-size: 1rem; font-weight: 700; flex-shrink: 0;
    font-family: 'IBM Plex Mono', monospace; border: 1.5px solid;
}
.rc-drug {
    font-family: 'Syne', sans-serif; font-size: 1rem;
    font-weight: 800; color: #e8f2ff; letter-spacing: -0.01em;
}
.rc-gene {
    font-family: 'IBM Plex Mono', monospace; font-size: 0.62rem;
    color: #2a4060; letter-spacing: 0.06em; margin-top: 2px;
}
.rc-pill {
    margin-left: auto; font-family: 'IBM Plex Mono', monospace;
    font-size: 0.65rem; font-weight: 600; padding: 4px 12px;
    border-radius: 5px; letter-spacing: 0.06em; border: 1px solid;
}
.rc-body { padding: 0 1.4rem 1.4rem 1.4rem; }

.metric-grid {
    display: grid; grid-template-columns: repeat(4, 1fr);
    gap: 1px; background: #111e2e; border-radius: 8px;
    overflow: hidden; margin-bottom: 1rem;
}
.mg-cell { background: #080e18; padding: 0.7rem 0.9rem; }
.mg-lbl {
    font-family: 'IBM Plex Mono', monospace; font-size: 0.56rem;
    color: #1e3050; letter-spacing: 0.1em; text-transform: uppercase; margin-bottom: 0.28rem;
}
.mg-val {
    font-family: 'IBM Plex Mono', monospace; font-size: 0.8rem;
    font-weight: 600; color: #6a90b0;
}

.conf-wrap { margin-bottom: 1rem; }
.conf-labels {
    display: flex; justify-content: space-between;
    font-family: 'IBM Plex Mono', monospace; font-size: 0.58rem;
    color: #1e3050; letter-spacing: 0.08em; margin-bottom: 0.35rem;
}
.conf-track { height: 3px; background: #111e2e; border-radius: 2px; overflow: hidden; }
.conf-fill  { height: 100%; border-radius: 2px; }

.sdiv {
    display: flex; align-items: center; gap: 0.7rem;
    margin: 1rem 0 0.7rem 0;
}
.sdiv-line { flex: 1; height: 1px; background: #0d1826; }
.sdiv-lbl {
    font-family: 'IBM Plex Mono', monospace; font-size: 0.58rem;
    color: #1e3050; letter-spacing: 0.12em; text-transform: uppercase; white-space: nowrap;
}

.vt { width: 100%; border-collapse: collapse; margin-bottom: 1rem; }
.vt th {
    font-family: 'IBM Plex Mono', monospace; font-size: 0.56rem;
    color: #1e3050; letter-spacing: 0.1em; text-transform: uppercase;
    padding: 0.45rem 0.7rem; text-align: left;
    border-bottom: 1px solid #0d1826; background: #050c16;
}
.vt td {
    font-family: 'IBM Plex Mono', monospace; font-size: 0.7rem;
    color: #4a7090; padding: 0.5rem 0.7rem; border-bottom: 1px solid #090f1c;
}
.vt tr:last-child td { border-bottom: none; }
.vt tr:hover td { background: #0a1520; color: #7aaac8; }
.v-rsid { color: #0ea5e9 !important; font-weight: 600; }
.v-star { color: #a78bfa !important; }
.v-fn   { color: #f87171 !important; }
.v-fd   { color: #fbbf24 !important; }
.v-fi   { color: #60a5fa !important; }
.v-fnm  { color: #34d399 !important; }

.rec-box {
    background: #050e18; border: 1px solid #0d1e32;
    border-left: 3px solid #0ea5e9;
    border-radius: 0 8px 8px 0; padding: 0.9rem 1.1rem; margin-bottom: 0.8rem;
}
.rb-lbl {
    font-family: 'IBM Plex Mono', monospace; font-size: 0.56rem;
    color: #0ea5e9; letter-spacing: 0.1em; text-transform: uppercase; margin-bottom: 0.35rem;
}
.rb-txt {
    font-family: 'DM Sans', sans-serif; font-size: 0.8rem;
    color: #6a90b0; line-height: 1.65;
}
.rec-box.warn  { border-left-color: #f59e0b; }
.rec-box.warn .rb-lbl { color: #f59e0b; }
.rec-box.danger{ border-left-color: #ef4444; }
.rec-box.danger .rb-lbl { color: #ef4444; }
.rec-box.purple{ border-left-color: #8b5cf6; }
.rec-box.purple .rb-lbl { color: #8b5cf6; }

.ai-wrap {
    background: #050c16; border: 1px solid #0d1826;
    border-radius: 10px; overflow: hidden; margin-bottom: 0.8rem;
}
.ai-head {
    background: #080f1c; padding: 0.65rem 1.1rem;
    border-bottom: 1px solid #0d1826;
    display: flex; align-items: center; gap: 0.6rem;
}
.ai-pill {
    background: linear-gradient(135deg, #0369a1, #4338ca);
    color: #fff; font-family: 'IBM Plex Mono', monospace;
    font-size: 0.52rem; font-weight: 600; padding: 2px 8px;
    border-radius: 3px; letter-spacing: 0.08em; text-transform: uppercase;
}
.ai-head-lbl {
    font-family: 'IBM Plex Mono', monospace; font-size: 0.6rem;
    color: #2a4060; letter-spacing: 0.08em; text-transform: uppercase;
}
.ai-block { padding: 0.85rem 1.1rem; border-bottom: 1px solid #080f1c; }
.ai-block:last-child { border-bottom: none; }
.ai-blk-lbl {
    font-family: 'IBM Plex Mono', monospace; font-size: 0.56rem;
    color: #1e3a54; letter-spacing: 0.1em; text-transform: uppercase; margin-bottom: 0.3rem;
}
.ai-blk-txt {
    font-family: 'DM Sans', sans-serif; font-size: 0.8rem;
    color: #5a8090; line-height: 1.68;
}

.alert-crit {
    background: #150202; border: 1px solid #5c1010; border-left: 3px solid #ef4444;
    border-radius: 0 8px 8px 0; padding: 0.8rem 1.1rem; margin-bottom: 1rem;
    display: flex; align-items: flex-start; gap: 0.7rem;
}
.alert-warn-box {
    background: #140e00; border: 1px solid #503800; border-left: 3px solid #f59e0b;
    border-radius: 0 8px 8px 0; padding: 0.8rem 1.1rem; margin-bottom: 1rem;
    display: flex; align-items: flex-start; gap: 0.7rem;
}
.alert-icon-lbl {
    font-family: 'IBM Plex Mono', monospace; font-size: 0.58rem;
    letter-spacing: 0.1em; text-transform: uppercase; margin-bottom: 0.2rem;
}
.alert-txt {
    font-family: 'DM Sans', sans-serif; font-size: 0.78rem; line-height: 1.55;
}

.ix-card {
    border: 1px solid #1a2b40; border-radius: 9px;
    overflow: hidden; margin-bottom: 0.7rem;
}
.ix-head {
    background: #080f1c; padding: 0.7rem 1.1rem;
    display: flex; align-items: center; gap: 0.7rem;
    border-bottom: 1px solid #0d1826;
}
.ix-dot { width: 7px; height: 7px; border-radius: 50%; flex-shrink: 0; }
.ix-type-lbl {
    font-family: 'IBM Plex Mono', monospace; font-size: 0.6rem;
    color: #2a4060; letter-spacing: 0.08em; text-transform: uppercase;
}
.ix-drugs-lbl {
    font-family: 'Syne', sans-serif; font-size: 0.8rem;
    font-weight: 700; color: #6a90b0; margin-left: auto;
}
.ix-sev {
    font-family: 'IBM Plex Mono', monospace; font-size: 0.58rem;
    font-weight: 600; padding: 2px 8px; border-radius: 3px;
    letter-spacing: 0.06em; text-transform: uppercase;
}
.ix-body { padding: 0.85rem 1.1rem; }
.ix-msg {
    font-family: 'DM Sans', sans-serif; font-size: 0.78rem;
    color: #3a6080; line-height: 1.62; margin-bottom: 0.55rem;
}
.ix-rec {
    font-family: 'IBM Plex Mono', monospace; font-size: 0.65rem;
    color: #1e3a54; line-height: 1.5;
    border-top: 1px solid #0d1826; padding-top: 0.55rem;
}

.gtag {
    display: inline-block; background: #0a1826; border: 1px solid #1a2b40;
    font-family: 'IBM Plex Mono', monospace; font-size: 0.62rem;
    color: #2a5070; padding: 2px 8px; border-radius: 4px; margin: 2px;
}

.stButton > button {
    background: linear-gradient(135deg, #0369a1, #4338ca) !important;
    color: #fff !important; border: none !important; border-radius: 8px !important;
    font-family: 'Syne', sans-serif !important; font-weight: 700 !important;
    letter-spacing: 0.01em !important; padding: 0.62rem 1.5rem !important;
    font-size: 0.88rem !important; box-shadow: 0 4px 15px #0369a122 !important;
}
.stButton > button:hover { opacity: 0.92 !important; }
.stTabs [data-baseweb="tab-list"] {
    background: #080c14 !important; border-bottom: 1px solid #111e2e !important; gap: 0 !important;
}
.stTabs [data-baseweb="tab"] {
    font-family: 'IBM Plex Mono', monospace !important; font-size: 0.68rem !important;
    color: #1e3050 !important; letter-spacing: 0.1em !important;
    padding: 0.75rem 1.5rem !important; background: transparent !important;
    border: none !important; text-transform: uppercase !important;
}
.stTabs [aria-selected="true"] { color: #0ea5e9 !important; border-bottom: 2px solid #0ea5e9 !important; }
.stTabs [data-baseweb="tab-panel"] { padding-top: 1.5rem !important; }
div[data-testid="stExpander"] {
    border: 1px solid #111e2e !important; border-radius: 9px !important;
    background: #080c14 !important; margin-bottom: 0.4rem;
}
div[data-testid="stExpander"] summary {
    font-family: 'IBM Plex Mono', monospace !important;
    font-size: 0.68rem !important; color: #2a4060 !important; letter-spacing: 0.07em !important;
}
.stDownloadButton > button {
    background: #0a1220 !important; color: #2a4060 !important;
    border: 1px solid #111e2e !important; border-radius: 7px !important;
    font-family: 'IBM Plex Mono', monospace !important;
    font-size: 0.62rem !important; letter-spacing: 0.06em !important;
}
.stDownloadButton > button:hover { border-color: #0ea5e9 !important; color: #0ea5e9 !important; }
.stMultiSelect span[data-baseweb="tag"] {
    background: #0d1e30 !important; color: #7dd3fc !important;
    border: 1px solid #1a3050 !important; font-family: 'IBM Plex Mono', monospace !important;
    font-size: 0.65rem !important;
}
.stProgress > div > div > div { background: linear-gradient(90deg, #0ea5e9, #6366f1) !important; }
.stCheckbox label p {
    font-family: 'DM Sans', sans-serif !important; font-size: 0.8rem !important;
    color: #3a5870 !important;
}
</style>
""", unsafe_allow_html=True)


# ─── Helpers ──────────────────────────────────────────────────
def load_vcf_file(filename: str) -> str:
    path = os.path.join(BASE_DIR, "sample_data", filename)
    if os.path.exists(path):
        with open(path) as f:
            return f.read()
    return get_sample_vcf()


def run_pipeline(vcf_content, drugs, pid, groq_key, run_ix=True, gen_pdf=True):
    parsed_vcf   = parse_vcf(vcf_content)
    risk_results = run_risk_assessment(parsed_vcf, drugs)
    if groq_key:
        risk_results = generate_all_explanations(groq_key, risk_results)
    else:
        for r in risk_results:
            r["llm_explanation"] = {
                "summary": "Add a Groq API key in the sidebar to enable AI clinical explanations.",
                "biological_mechanism": "", "variant_significance": "",
                "clinical_implications": "", "success": False
            }
    all_outputs = [
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
            st.warning(f"PDF generation error: {e}")
    return parsed_vcf, risk_results, all_outputs, ix_report, pdf_bytes


def func_css(status: str) -> str:
    s = (status or "").lower()
    if "no_function" in s or "no function" in s: return "v-fn"
    if "decreased" in s:  return "v-fd"
    if "increased" in s:  return "v-fi"
    return "v-fnm"


def render_results(all_outputs, parsed_vcf, ix_report, pdf_bytes, pid):
    # Summary Bar
    overall_sev = max((o["risk_assessment"]["severity"] for o in all_outputs),
                      key=lambda s: SEV_RANK.get(s, 0), default="none")
    sev_cfg   = SEV_CONFIG.get(overall_sev, SEV_CONFIG["none"])
    genes_str = "  ·  ".join(parsed_vcf.get("detected_genes", [])) or "none detected"

    st.markdown(f"""
    <div class="summary-bar">
      <div class="sc">
        <div class="sc-lbl">Patient</div>
        <div class="sc-val" style="font-size:0.88rem;font-family:'IBM Plex Mono',monospace;color:#4a7090;">{pid}</div>
      </div>
      <div class="sc">
        <div class="sc-lbl">Overall Risk</div>
        <div class="sc-val" style="color:{sev_cfg['color']};">{sev_cfg['label']}</div>
      </div>
      <div class="sc">
        <div class="sc-lbl">Drugs Analyzed</div>
        <div class="sc-val">{len(all_outputs)}</div>
        <div class="sc-sub">selected</div>
      </div>
      <div class="sc">
        <div class="sc-lbl">Variants Found</div>
        <div class="sc-val">{parsed_vcf['total_variants']}</div>
        <div class="sc-sub">pharmacogenomic</div>
      </div>
      <div class="sc">
        <div class="sc-lbl">Genes Covered</div>
        <div class="sc-val">{len(parsed_vcf['detected_genes'])}<span style="font-size:0.9rem;color:#1a2e44;">/6</span></div>
        <div class="sc-sub" style="font-size:0.5rem;">{genes_str}</div>
      </div>
    </div>
    """, unsafe_allow_html=True)

    # Downloads row
    all_json = json.dumps(all_outputs, indent=2)
    dc1, dc2, dc3 = st.columns(3)
    with dc1:
        st.download_button("⬇  All Results (JSON)", data=all_json,
                           file_name=f"pharmaguard_{pid}_all.json",
                           mime="application/json", use_container_width=True, key=f"dlall_{pid}")
    with dc2:
        if pdf_bytes:
            st.download_button("⬇  Clinical PDF Report", data=pdf_bytes,
                               file_name=f"pharmaguard_{pid}_report.pdf",
                               mime="application/pdf", use_container_width=True, key=f"dlpdf_{pid}")
    with dc3:
        if ix_report and ix_report.get("interactions_found"):
            ix_json = json.dumps(ix_report, indent=2)
            st.download_button("⬇  Interaction Report", data=ix_json,
                               file_name=f"pharmaguard_{pid}_interactions.json",
                               mime="application/json", use_container_width=True, key=f"dlix_{pid}")

    st.markdown("<div style='height:0.6rem'></div>", unsafe_allow_html=True)

    # Drug-Drug Interactions
    if ix_report and ix_report["interactions_found"]:
        ix_sev = ix_report["overall_severity"]
        sev_colors = {"low":"#f59e0b","moderate":"#f97316","high":"#ef4444","critical":"#dc2626"}
        sev_bgs    = {"low":"#1c1400","moderate":"#1a0800","high":"#1a0505","critical":"#120000"}
        sev_txts   = {"low":"#fbbf24","moderate":"#fb923c","high":"#f87171","critical":"#fca5a5"}
        d_col = sev_colors.get(ix_sev,"#6b7280")

        with st.expander(f"⚡  DRUG–DRUG INTERACTIONS  ·  {ix_report['total_interactions']} ALERT(S)  ·  {ix_sev.upper()}", expanded=True):
            for ix in ix_report["all_interactions"]:
                sev  = ix.get("severity","low")
                dc   = sev_colors.get(sev,"#6b7280")
                dbg  = sev_bgs.get(sev,"#111827")
                dtxt = sev_txts.get(sev,"#9ca3af")
                drugs_str = " + ".join(
                    ix.get("drugs_involved",[]) or
                    [ix.get("inhibitor_drug","")] + ix.get("affected_drugs",[])
                )
                msg = ix.get("message", ix.get("mechanism",""))
                rec = ix.get("recommendation","")
                st.markdown(f"""
                <div class="ix-card" style="border-color:{dc}33;">
                  <div class="ix-head" style="background:{dbg};">
                    <div class="ix-dot" style="background:{dc};"></div>
                    <span class="ix-type-lbl">{ix.get('type','').replace('_',' ').title()}</span>
                    <span class="ix-drugs-lbl">{drugs_str}</span>
                    <span class="ix-sev" style="background:{dbg};color:{dtxt};border:1px solid {dc}44;">{sev.upper()}</span>
                  </div>
                  <div class="ix-body">
                    <div class="ix-msg">{msg}</div>
                    <div class="ix-rec">💡  {rec}</div>
                  </div>
                </div>""", unsafe_allow_html=True)
    elif ix_report:
        st.markdown("""
        <div style="background:#022c22;border:1px solid #065f46;border-radius:8px;
             padding:0.7rem 1.1rem;margin-bottom:1rem;font-family:'IBM Plex Mono',monospace;
             font-size:0.65rem;color:#10b981;letter-spacing:0.07em;">
          ✓  NO SIGNIFICANT DRUG–DRUG INTERACTIONS DETECTED
        </div>""", unsafe_allow_html=True)

    # Per-Drug Cards
    for output in all_outputs:
        risk_label = output["risk_assessment"]["risk_label"]
        drug_name  = output["drug"]
        severity   = output["risk_assessment"]["severity"]
        confidence = output["risk_assessment"]["confidence_score"]
        gene       = output["pharmacogenomic_profile"]["primary_gene"]
        diplotype  = output["pharmacogenomic_profile"]["diplotype"]
        phenotype  = output["pharmacogenomic_profile"]["phenotype"]
        variants   = output["pharmacogenomic_profile"]["detected_variants"]
        rec        = output["clinical_recommendation"]["dosing_recommendation"]
        alts       = output["clinical_recommendation"].get("alternative_drugs", [])
        monitoring = output["clinical_recommendation"].get("monitoring_required", "")
        exp        = output["llm_generated_explanation"]
        rc         = RISK_CONFIG.get(risk_label, RISK_CONFIG["Unknown"])
        sc         = SEV_CONFIG.get(severity, SEV_CONFIG["none"])

        # Alert banners
        if severity == "critical":
            st.markdown(f"""
            <div class="alert-crit">
              <div style="font-size:1rem;color:#ef4444;flex-shrink:0;margin-top:1px;">⛔</div>
              <div>
                <div class="alert-icon-lbl" style="color:#ef4444;">Critical Safety Alert — {drug_name}</div>
                <div class="alert-txt" style="color:#fca5a5;">{rec}</div>
              </div>
            </div>""", unsafe_allow_html=True)
        elif severity == "high":
            st.markdown(f"""
            <div class="alert-warn-box">
              <div style="font-size:1rem;color:#f59e0b;flex-shrink:0;margin-top:1px;">⚠</div>
              <div>
                <div class="alert-icon-lbl" style="color:#f59e0b;">Clinical Warning — {drug_name}</div>
                <div class="alert-txt" style="color:#fde68a;">{rec}</div>
              </div>
            </div>""", unsafe_allow_html=True)

        # Card
        st.markdown(f"""
        <div class="risk-card" style="background:{rc['bg']};border-color:{rc['border']}55;">
          <div class="rc-head">
            <div class="rc-icon" style="background:{rc['bg']};border-color:{rc['border']};color:{rc['text']};">
              {rc['icon']}
            </div>
            <div>
              <div class="rc-drug">{drug_name.title()}</div>
              <div class="rc-gene">{gene}  ·  {diplotype}  ·  {phenotype}</div>
            </div>
            <div class="rc-pill" style="background:{rc['bg']};border-color:{rc['border']};color:{rc['text']};">
              {risk_label.upper()}
            </div>
          </div>
          <div class="rc-body">
            <div class="metric-grid">
              <div class="mg-cell">
                <div class="mg-lbl">Phenotype</div>
                <div class="mg-val" style="color:{rc['text']};">{phenotype}</div>
              </div>
              <div class="mg-cell">
                <div class="mg-lbl">Severity</div>
                <div class="mg-val" style="color:{sc['color']};">{sc['label']}</div>
              </div>
              <div class="mg-cell">
                <div class="mg-lbl">Confidence</div>
                <div class="mg-val">{confidence:.0%}</div>
              </div>
              <div class="mg-cell">
                <div class="mg-lbl">Variants</div>
                <div class="mg-val">{len(variants)}</div>
              </div>
            </div>
            <div class="conf-wrap">
              <div class="conf-labels">
                <span>PREDICTION CONFIDENCE</span>
                <span style="color:{rc['dot']};">{confidence:.0%}</span>
              </div>
              <div class="conf-track">
                <div class="conf-fill" style="width:{confidence*100:.1f}%;background:{rc['dot']};"></div>
              </div>
            </div>
        """, unsafe_allow_html=True)

        # Variants table
        if variants:
            rows_html = ""
            for v in variants:
                fc = func_css(v.get("functional_status",""))
                fn_label = (v.get("functional_status") or "unknown").replace("_"," ").title()
                rows_html += f"""<tr>
                  <td class="v-rsid">{v.get('rsid','—')}</td>
                  <td class="v-star">{v.get('star_allele','—')}</td>
                  <td>{v.get('ref','?')} → {v.get('alt','?')}</td>
                  <td>{v.get('chrom','—')}:{v.get('pos','—')}</td>
                  <td class="{fc}">{fn_label}</td>
                </tr>"""
            st.markdown(f"""
            <div class="sdiv"><div class="sdiv-line"></div>
              <div class="sdiv-lbl">Detected Variants ({len(variants)})</div>
              <div class="sdiv-line"></div></div>
            <table class="vt">
              <thead><tr>
                <th>rsID</th><th>Star Allele</th><th>Nucleotide Change</th>
                <th>Genomic Position</th><th>Functional Status</th>
              </tr></thead>
              <tbody>{rows_html}</tbody>
            </table>""", unsafe_allow_html=True)
        else:
            st.markdown("""
            <div style="font-family:'IBM Plex Mono',monospace;font-size:0.65rem;
                 color:#1e3050;padding:0.4rem 0 0.8rem 0;letter-spacing:0.07em;">
              ○  No variants detected — wild-type (*1/*1) assumed — normal enzyme function
            </div>""", unsafe_allow_html=True)

        # CPIC Rec
        rec_cls = {"Toxic":"danger","Ineffective":"purple","Adjust Dosage":"warn"}.get(risk_label,"")
        st.markdown(f"""
        <div class="sdiv"><div class="sdiv-line"></div>
          <div class="sdiv-lbl">CPIC Dosing Recommendation</div>
          <div class="sdiv-line"></div></div>
        <div class="rec-box {rec_cls}">
          <div class="rb-lbl">CPIC Guideline for {drug_name}</div>
          <div class="rb-txt">{rec}</div>
        </div>""", unsafe_allow_html=True)

        if alts:
            alts_tags = "".join(f'<span class="gtag">{a}</span>' for a in alts)
            st.markdown(f"""
            <div style="margin-bottom:0.7rem;">
              <span style="font-family:'IBM Plex Mono',monospace;font-size:0.56rem;
                   color:#1e3050;letter-spacing:0.1em;text-transform:uppercase;
                   margin-right:0.5rem;">Alternative Drugs</span>
              {alts_tags}
            </div>""", unsafe_allow_html=True)

        if monitoring:
            st.markdown(f"""
            <div style="font-family:'DM Sans',sans-serif;font-size:0.76rem;color:#2a4a60;
                 padding:0.5rem 0 0.7rem 0;border-top:1px solid #0d1826;margin-top:0.3rem;">
              <span style="font-family:'IBM Plex Mono',monospace;font-size:0.56rem;
                   color:#1e3050;letter-spacing:0.1em;text-transform:uppercase;
                   margin-right:0.5rem;">Monitoring</span>
              {monitoring}
            </div>""", unsafe_allow_html=True)

        # AI Explanation
        if exp.get("summary"):
            blocks_html = ""
            for lbl, key in [("Summary","summary"),("Biological Mechanism","biological_mechanism"),
                              ("Variant Significance","variant_significance"),
                              ("Clinical Implications","clinical_implications")]:
                if exp.get(key):
                    blocks_html += f"""<div class="ai-block">
                      <div class="ai-blk-lbl">{lbl}</div>
                      <div class="ai-blk-txt">{exp[key]}</div>
                    </div>"""
            st.markdown(f"""
            <div class="sdiv"><div class="sdiv-line"></div>
              <div class="sdiv-lbl">AI Clinical Explanation</div>
              <div class="sdiv-line"></div></div>
            <div class="ai-wrap">
              <div class="ai-head">
                <span class="ai-pill">Groq LLaMA 3.3</span>
                <span class="ai-head-lbl">Pharmacogenomic Analysis · {drug_name}</span>
              </div>
              {blocks_html}
            </div>""", unsafe_allow_html=True)

        # Close card
        st.markdown("</div></div>", unsafe_allow_html=True)

        # Raw JSON expander
        json_str = json.dumps(output, indent=2)
        with st.expander(f"{{}} Raw JSON Output — {drug_name}", expanded=False):
            c1, c2 = st.columns(2)
            with c1:
                st.download_button("⬇  Download JSON", data=json_str,
                                   file_name=f"pharmaguard_{pid}_{drug_name}.json",
                                   mime="application/json", key=f"dl_{pid}_{drug_name}",
                                   use_container_width=True)
            st.code(json_str, language="json")
        st.markdown("<div style='height:0.2rem'></div>", unsafe_allow_html=True)


# ══════════════════════════════════════════════════════════════
# MASTHEAD
# ══════════════════════════════════════════════════════════════
st.markdown("""
<div class="masthead">
  <div class="masthead-left">
    <div class="masthead-logo">🧬</div>
    <div>
      <div class="masthead-title">PharmaGuard</div>
      <div class="masthead-sub">Genomic Risk Intelligence  ·  RIFT 2026  ·  CPIC Aligned  ·  Groq LLaMA 3.3</div>
    </div>
  </div>
  <div class="masthead-tags">
    <span class="mtag mtag-blue">RIFT 2026</span>
    <span class="mtag mtag-green">CPIC Aligned</span>
    <span class="mtag mtag-violet">v2.0</span>
  </div>
</div>
""", unsafe_allow_html=True)

# Sidebar
with st.sidebar:
    st.markdown("""<div style="font-family:'IBM Plex Mono',monospace;font-size:0.58rem;
        color:#1e3050;letter-spacing:0.12em;text-transform:uppercase;margin-bottom:0.4rem;">
        Groq API Key</div>""", unsafe_allow_html=True)
    groq_api_key = st.text_input("", value=os.environ.get("GROQ_API_KEY",""),
                                  type="password", label_visibility="collapsed",
                                  placeholder="gsk_...")
    st.markdown("""<div style="font-family:'IBM Plex Mono',monospace;font-size:0.56rem;
        color:#1a2a40;margin-top:0.3rem;margin-bottom:1rem;">
        Get free key at console.groq.com</div>""", unsafe_allow_html=True)
    st.markdown("---")
    for g, drug in [("CYP2D6","Codeine"),("CYP2C19","Clopidogrel"),("CYP2C9","Warfarin"),
                    ("SLCO1B1","Simvastatin"),("TPMT","Azathioprine"),("DPYD","Fluorouracil")]:
        st.markdown(f"""<div style="font-family:'IBM Plex Mono',monospace;font-size:0.62rem;
            color:#2a4060;padding:3px 0;">
            <span style="color:#0ea5e9;">{g}</span>
            <span style="color:#1a2a40;"> → </span>
            <span style="color:#3a5870;">{drug}</span></div>""", unsafe_allow_html=True)

# TABS
tab1, tab2 = st.tabs(["ANALYSIS", "TEST SUITE"])


# ── TAB 1 ─────────────────────────────────────────────────────
with tab1:
    st.markdown("""
    <div class="step-row">
      <div class="step-item on"><div class="step-num">01</div><div class="step-name">Upload VCF</div></div>
      <div class="step-item on"><div class="step-num">02</div><div class="step-name">Select Drugs</div></div>
      <div class="step-item"><div class="step-num">03</div><div class="step-name">Run Analysis</div></div>
      <div class="step-item"><div class="step-num">04</div><div class="step-name">Review Results</div></div>
    </div>
    """, unsafe_allow_html=True)

    col_l, col_r = st.columns([1.3, 1], gap="large")

    with col_l:
        st.markdown("""<span style="font-family:'IBM Plex Mono',monospace;font-size:0.58rem;
            color:#1e3050;letter-spacing:0.12em;text-transform:uppercase;display:block;
            margin-bottom:0.5rem;">STEP 01 — GENOMIC DATA</span>""", unsafe_allow_html=True)

        uploaded_file = st.file_uploader(
            "Upload VCF file (max 5 MB)",
            type=["vcf"],
            help="VCF v4.2 with GENE, STAR, FUNCTION INFO tags"
        )
        if uploaded_file is not None:
            sz = uploaded_file.size / (1024*1024)
            if sz > 5:
                st.error(f"File too large: {sz:.1f} MB — maximum is 5 MB")
                uploaded_file = None
            else:
                peek = uploaded_file.read(500).decode("utf-8", errors="replace")
                uploaded_file.seek(0)
                if "##fileformat=VCF" not in peek and "#CHROM" not in peek:
                    st.error("Invalid VCF — must contain ##fileformat=VCF header")
                    uploaded_file = None
                else:
                    st.markdown(f"""
                    <div style="background:#022c22;border:1px solid #065f4640;border-radius:7px;
                         padding:0.55rem 0.9rem;font-family:'IBM Plex Mono',monospace;
                         font-size:0.65rem;color:#10b981;margin-top:0.3rem;letter-spacing:0.06em;">
                      ✓  {uploaded_file.name}  ·  {sz:.2f} MB
                    </div>""", unsafe_allow_html=True)

        st.markdown("""<span style="font-family:'IBM Plex Mono',monospace;font-size:0.58rem;
            color:#1e3050;letter-spacing:0.12em;text-transform:uppercase;
            display:block;margin-top:1rem;margin-bottom:0.5rem;">
            OR USE TEST SCENARIO</span>""", unsafe_allow_html=True)

        scenario_opts = {
            "— none selected —": None,
            "🧬  Mixed Variants (Standard)": "sample.vcf",
            "🔴  UltraRapid Metabolizer — Codeine TOXIC": "test_ultrarapid_metabolizer.vcf",
            "🟢  All Normal Wild-type — All Safe": "test_all_normal_wildtype.vcf",
            "🚨  Worst Case — All Poor Metabolizers": "test_worst_case_all_pm.vcf",
        }
        chosen_label = st.selectbox("", list(scenario_opts.keys()), label_visibility="collapsed")
        chosen_file  = scenario_opts[chosen_label]

    with col_r:
        st.markdown("""<span style="font-family:'IBM Plex Mono',monospace;font-size:0.58rem;
            color:#1e3050;letter-spacing:0.12em;text-transform:uppercase;display:block;
            margin-bottom:0.5rem;">STEP 02 — DRUG SELECTION</span>""", unsafe_allow_html=True)

        drug_multiselect = st.multiselect(
            "Select drugs",
            options=ALL_DRUGS,
            default=["CLOPIDOGREL"],
            format_func=lambda x: f"{x.title()}  ({GENE_DRUG_MAP.get(x,'')})",
            label_visibility="collapsed"
        )

        st.markdown("""<span style="font-family:'IBM Plex Mono',monospace;font-size:0.56rem;
            color:#1a2a40;letter-spacing:0.1em;text-transform:uppercase;display:block;
            margin-top:0.7rem;margin-bottom:0.4rem;">OR TYPE MANUALLY</span>""",
                    unsafe_allow_html=True)
        custom_drugs = st.text_input("", placeholder="e.g. CODEINE, WARFARIN, AZATHIOPRINE",
                                     label_visibility="collapsed")

        st.markdown("""<span style="font-family:'IBM Plex Mono',monospace;font-size:0.56rem;
            color:#1a2a40;letter-spacing:0.1em;text-transform:uppercase;display:block;
            margin-top:0.7rem;margin-bottom:0.4rem;">PATIENT ID</span>""",
                    unsafe_allow_html=True)
        patient_id_input = st.text_input("", placeholder="Auto-generated if blank",
                                         label_visibility="collapsed")

        c1, c2 = st.columns(2)
        with c1: run_interactions = st.checkbox("Drug–Drug Interactions", value=True)
        with c2: generate_pdf    = st.checkbox("PDF Report", value=True)

    st.markdown("<div style='height:0.6rem'></div>", unsafe_allow_html=True)
    analyze_btn = st.button("▶  Run Pharmacogenomic Analysis", use_container_width=True)

    if analyze_btn:
        all_drugs = list(drug_multiselect)
        if custom_drugs.strip():
            all_drugs += [d.strip().upper() for d in custom_drugs.split(",") if d.strip()]
        all_drugs = list(set(all_drugs))

        if not all_drugs:
            st.error("Select at least one drug to analyze.")
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
        <div style="display:flex;align-items:baseline;gap:1rem;margin:1.8rem 0 1.2rem 0;
             border-bottom:1px solid #0d1826;padding-bottom:1rem;">
          <div style="font-family:'Syne',sans-serif;font-size:1.4rem;font-weight:800;
               color:#ddeeff;letter-spacing:-0.02em;">Analysis Results</div>
          <div style="font-family:'IBM Plex Mono',monospace;font-size:0.7rem;
               color:#1e3050;letter-spacing:0.08em;">{pid}</div>
        </div>
        """, unsafe_allow_html=True)

        with st.spinner("Running pharmacogenomic analysis..."):
            parsed_vcf, risk_results, all_outputs, ix_report, pdf_bytes = run_pipeline(
                vcf_content, all_drugs, pid, groq_api_key, run_interactions, generate_pdf)

        render_results(all_outputs, parsed_vcf, ix_report, pdf_bytes, pid)

        with st.expander("VCF PARSING DETAILS", expanded=False):
            p1, p2, p3 = st.columns(3)
            p1.metric("Total Variants", parsed_vcf["total_variants"])
            p2.metric("Genes Found",    len(parsed_vcf["detected_genes"]))
            p3.metric("Parse Errors",   len(parsed_vcf["parse_errors"]))
            if parsed_vcf["detected_genes"]:
                tags = "".join(f'<span class="gtag">{g}</span>' for g in parsed_vcf["detected_genes"])
                st.markdown(tags, unsafe_allow_html=True)
            for err in parsed_vcf.get("parse_errors",[]):
                st.markdown(f"""<div style="font-family:'IBM Plex Mono',monospace;
                    font-size:0.62rem;color:#ef4444;padding:2px 0;">{err}</div>""",
                            unsafe_allow_html=True)
    else:
        st.markdown("""
        <div style="text-align:center;padding:5rem 2rem;border:1px dashed #0d1826;
             border-radius:14px;margin-top:0.5rem;background:#050c16;">
          <div style="font-size:2.5rem;margin-bottom:1rem;opacity:0.15;">🧬</div>
          <div style="font-family:'Syne',sans-serif;font-size:1rem;font-weight:700;
               color:#1a2e44;margin-bottom:0.5rem;">Ready for Genomic Analysis</div>
          <div style="font-family:'IBM Plex Mono',monospace;font-size:0.62rem;
               color:#111e2e;letter-spacing:0.07em;line-height:2.2;">
            Upload a VCF file  ·  Select drugs  ·  Click Run Analysis<br>
            CYP2D6  ·  CYP2C19  ·  CYP2C9  ·  SLCO1B1  ·  TPMT  ·  DPYD
          </div>
        </div>
        """, unsafe_allow_html=True)


# ── TAB 2 — TEST SUITE ────────────────────────────────────────
TEST_SUITE = [
    {
        "name":     "Mixed Variants (Standard)",
        "file":     "sample.vcf",
        "drugs":    ["CLOPIDOGREL","CODEINE","AZATHIOPRINE"],
        "expected": {"CLOPIDOGREL":"Ineffective","CODEINE":"Adjust Dosage","AZATHIOPRINE":"Adjust Dosage"},
        "desc":     "CYP2C19 *2/*3 · CYP2D6 *4/*1 · TPMT *3B/*1",
        "icon":     "🧬",
    },
    {
        "name":     "UltraRapid Metabolizer — Codeine TOXIC",
        "file":     "test_ultrarapid_metabolizer.vcf",
        "drugs":    ["CODEINE","CLOPIDOGREL"],
        "expected": {"CODEINE":"Toxic","CLOPIDOGREL":"Safe"},
        "desc":     "CYP2D6 *1xN/*1xN duplication → URM → Toxic",
        "icon":     "🔴",
    },
    {
        "name":     "All Normal Wild-type — All Safe",
        "file":     "test_all_normal_wildtype.vcf",
        "drugs":    ALL_DRUGS,
        "expected": {d:"Safe" for d in ALL_DRUGS},
        "desc":     "Wild-type *1/*1 across all 6 genes — all Safe",
        "icon":     "🟢",
    },
    {
        "name":     "Worst Case — All Poor Metabolizers",
        "file":     "test_worst_case_all_pm.vcf",
        "drugs":    ALL_DRUGS,
        "expected": {
            "CODEINE":      "Ineffective",
            "CLOPIDOGREL":  "Ineffective",
            "WARFARIN":     "Adjust Dosage",
            "SIMVASTATIN":  "Toxic",
            "AZATHIOPRINE": "Toxic",
            "FLUOROURACIL": "Toxic",
        },
        "desc":     "Loss-of-function alleles all 6 genes — stress test",
        "icon":     "🚨",
    },
]

with tab2:
    st.markdown("""
    <div style="margin-bottom:1.5rem;">
      <div style="font-family:'Syne',sans-serif;font-size:1.2rem;font-weight:800;
           color:#ddeeff;letter-spacing:-0.02em;margin-bottom:0.3rem;">Automated Test Suite</div>
      <div style="font-family:'IBM Plex Mono',monospace;font-size:0.6rem;color:#1e3050;
           letter-spacing:0.08em;">
        4 VCF SCENARIOS  ·  VALIDATES RISK ENGINE  ·  EXPECTED vs ACTUAL  ·  JUDGE VERIFICATION
      </div>
    </div>
    """, unsafe_allow_html=True)

    # Preview cards
    pc = st.columns(4)
    for i, sc in enumerate(TEST_SUITE):
        with pc[i]:
            exp_rows = ""
            for d, r in list(sc["expected"].items())[:4]:
                rc_text = RISK_CONFIG.get(r, RISK_CONFIG["Unknown"])["text"]
                exp_rows += (
                    f"<div style='display:flex;gap:0.4rem;padding:1px 0;'>"
                    f"<span style='color:#2a4060;font-size:0.6rem;width:70px;flex-shrink:0;'>{d[:7]}</span>"
                    f"<span style='font-size:0.6rem;color:{rc_text};'>{r}</span></div>"
                )
            st.markdown(f"""
            <div style="background:#080e18;border:1px solid #111e2e;border-radius:10px;
                 padding:1rem;min-height:180px;">
              <div style="font-size:1.1rem;margin-bottom:0.5rem;">{sc['icon']}</div>
              <div style="font-family:'Syne',sans-serif;font-size:0.75rem;font-weight:700;
                   color:#5a80a0;margin-bottom:0.3rem;line-height:1.4;">{sc['name']}</div>
              <div style="font-family:'IBM Plex Mono',monospace;font-size:0.58rem;
                   color:#1a2e44;margin-bottom:0.7rem;line-height:1.5;">{sc['desc']}</div>
              <div style="font-family:'IBM Plex Mono',monospace;">{exp_rows}</div>
            </div>
            """, unsafe_allow_html=True)

    st.markdown("<div style='height:1rem'></div>", unsafe_allow_html=True)
    tc1, tc2 = st.columns([3,1])
    with tc1:
        use_llm = st.checkbox("Include LLM Explanations (uses Groq API credits, ~30s slower)", value=False)
    with tc2:
        run_all = st.button("▶  Run All 4 Tests", use_container_width=True)

    if run_all:
        st.markdown("<div style='height:0.5rem'></div>", unsafe_allow_html=True)
        passed_count, failed_count = 0, 0
        all_suite_outputs = []

        for sc in TEST_SUITE:
            vcf = load_vcf_file(sc["file"])
            pid = f"TEST_{sc['name'][:8].replace(' ','').upper()[:8]}"
            key = groq_api_key if use_llm else ""
            with st.spinner(f"Running: {sc['name']}..."):
                pv, _, ao, _, _ = run_pipeline(vcf, sc["drugs"], pid, key, run_ix=False, gen_pdf=False)
            rows, sc_pass = [], True
            for out in ao:
                drug  = out["drug"]
                got   = out["risk_assessment"]["risk_label"]
                exp   = sc["expected"].get(drug,"")
                ok    = (got==exp) if exp else True
                pheno = out["pharmacogenomic_profile"]["phenotype"]
                diplo = out["pharmacogenomic_profile"]["diplotype"]
                rows.append((drug,got,exp,ok,pheno,diplo))
                if not ok: sc_pass = False
            if sc_pass: passed_count+=1
            else:       failed_count+=1
            all_suite_outputs.append({
                "scenario":sc["name"],"pass":sc_pass,
                "rows":rows,"outputs":ao,"file":sc["file"]
            })

        # Banner
        total = passed_count + failed_count
        is_all = failed_count == 0
        b_bg  = "#022c22" if is_all else "#1c1000"
        b_bdr = "#065f46" if is_all else "#503800"
        b_col = "#10b981" if is_all else "#f59e0b"
        b_sym = "✓" if is_all else "⚠"
        st.markdown(f"""
        <div style="background:{b_bg};border:1px solid {b_bdr};border-radius:12px;
             padding:1.3rem 2rem;margin:1rem 0;text-align:center;">
          <div style="font-family:'Syne',sans-serif;font-size:1.5rem;font-weight:800;
               color:{b_col};letter-spacing:-0.02em;">
            {b_sym}  {'ALL 4 TESTS PASSED' if is_all else f'{passed_count}/{total} TESTS PASSED'}
          </div>
          <div style="font-family:'IBM Plex Mono',monospace;font-size:0.6rem;
               color:{b_col};opacity:0.55;margin-top:0.3rem;letter-spacing:0.1em;">
            {passed_count} PASSED  ·  {failed_count} FAILED  ·  {int(passed_count/total*100)}% PASS RATE
          </div>
        </div>
        """, unsafe_allow_html=True)

        for sr in all_suite_outputs:
            sc_name = sr["scenario"]
            sc_pass = sr["pass"]
            sc_icon = next((s["icon"] for s in TEST_SUITE if s["name"]==sc_name),"🧪")
            sym = "✓" if sc_pass else "✕"
            hc = "#10b981" if sc_pass else "#ef4444"

            with st.expander(f"{sym}  {sc_icon}  {sc_name}  —  {'PASS' if sc_pass else 'FAIL'}",
                             expanded=not sc_pass):
                # Header row
                st.markdown("""
                <div style="display:grid;grid-template-columns:1.1fr 1.3fr 1.3fr 1.6fr 0.35fr;
                     background:#050c16;border:1px solid #0d1826;
                     border-radius:8px 8px 0 0;border-bottom:none;">
                  <div style="font-family:'IBM Plex Mono',monospace;font-size:0.55rem;
                       color:#1a2e44;letter-spacing:0.1em;padding:0.55rem 0.85rem;">DRUG</div>
                  <div style="font-family:'IBM Plex Mono',monospace;font-size:0.55rem;
                       color:#1a2e44;letter-spacing:0.1em;padding:0.55rem 0.85rem;">RESULT</div>
                  <div style="font-family:'IBM Plex Mono',monospace;font-size:0.55rem;
                       color:#1a2e44;letter-spacing:0.1em;padding:0.55rem 0.85rem;">EXPECTED</div>
                  <div style="font-family:'IBM Plex Mono',monospace;font-size:0.55rem;
                       color:#1a2e44;letter-spacing:0.1em;padding:0.55rem 0.85rem;">DIPLOTYPE / PHENOTYPE</div>
                  <div style="padding:0.55rem 0.85rem;"></div>
                </div>
                """, unsafe_allow_html=True)

                rows_html = ""
                for drug,got,exp,ok,pheno,diplo in sr["rows"]:
                    rc   = RISK_CONFIG.get(got,RISK_CONFIG["Unknown"])
                    ok_c = "#10b981" if ok else "#ef4444"
                    ok_s = "✓" if ok else "✕"
                    rows_html += f"""
                    <div style="display:grid;grid-template-columns:1.1fr 1.3fr 1.3fr 1.6fr 0.35fr;
                         border:1px solid #0d1826;border-top:none;background:#080c14;">
                      <div style="font-family:'IBM Plex Mono',monospace;font-size:0.68rem;
                           color:#2a4060;padding:0.55rem 0.85rem;">{drug}</div>
                      <div style="padding:0.55rem 0.85rem;">
                        <span style="font-family:'IBM Plex Mono',monospace;font-size:0.65rem;
                             font-weight:600;color:{rc['text']};background:{rc['bg']};
                             border:1px solid {rc['border']}44;padding:2px 8px;border-radius:4px;">
                          {got}</span>
                      </div>
                      <div style="font-family:'IBM Plex Mono',monospace;font-size:0.68rem;
                           color:#1a2e44;padding:0.55rem 0.85rem;">{exp if exp else '—'}</div>
                      <div style="font-family:'IBM Plex Mono',monospace;font-size:0.65rem;
                           color:#1e3a54;padding:0.55rem 0.85rem;">{diplo} / {pheno}</div>
                      <div style="font-family:'IBM Plex Mono',monospace;font-size:0.78rem;
                           font-weight:700;color:{ok_c};padding:0.55rem 0.85rem;
                           text-align:center;">{ok_s}</div>
                    </div>"""
                rows_html += """<div style="border:1px solid #0d1826;border-top:none;
                    border-radius:0 0 8px 8px;background:#050c16;height:4px;"></div>"""
                st.markdown(rows_html, unsafe_allow_html=True)

                st.markdown("<div style='height:0.8rem'></div>", unsafe_allow_html=True)
                dl1, dl2 = st.columns(2)
                sc_json = json.dumps(sr["outputs"], indent=2)
                with dl1:
                    st.download_button("⬇  Download JSON", data=sc_json,
                                       file_name=f"test_{sr['file'].replace('.vcf','')}.json",
                                       mime="application/json", key=f"tsc_{sc_name[:14]}",
                                       use_container_width=True)
                with dl2:
                    vcf_raw = load_vcf_file(sr["file"])
                    st.download_button("⬇  Download VCF", data=vcf_raw,
                                       file_name=sr["file"], mime="text/plain",
                                       key=f"vcf_{sc_name[:14]}", use_container_width=True)

        st.markdown("<div style='height:0.5rem'></div>", unsafe_allow_html=True)
        full_json = json.dumps(
            [{"scenario":s["scenario"],"pass":s["pass"],"results":s["outputs"]}
             for s in all_suite_outputs], indent=2)
        st.download_button("⬇  Download Complete Test Suite JSON", data=full_json,
                           file_name=f"pharmaguard_test_suite_{datetime.now().strftime('%Y%m%d_%H%M')}.json",
                           mime="application/json", use_container_width=True)
    else:
        st.markdown("""
        <div style="text-align:center;padding:3.5rem 2rem;border:1px dashed #0d1826;
             border-radius:14px;margin-top:0.5rem;background:#050c16;">
          <div style="font-size:2rem;margin-bottom:0.8rem;opacity:0.15;
               font-family:'IBM Plex Mono',monospace;color:#0ea5e9;">▶</div>
          <div style="font-family:'Syne',sans-serif;font-size:0.95rem;font-weight:700;
               color:#1a2e44;margin-bottom:0.4rem;">One-Click Full Validation</div>
          <div style="font-family:'IBM Plex Mono',monospace;font-size:0.6rem;
               color:#0d1826;letter-spacing:0.07em;line-height:2.3;">
            4 VCF scenarios  ·  Expected vs actual risk labels  ·  Pass / fail per drug<br>
            Mixed Variants  ·  UltraRapid  ·  All Normal  ·  Worst Case All PM
          </div>
        </div>
        """, unsafe_allow_html=True)