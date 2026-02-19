"""
PharmaGuard v2.0 — Pharmacogenomic Risk Prediction System
RIFT 2026 Hackathon | Pharmacogenomics / Explainable AI Track
COMPLETE REDESIGN — Light Clinical Intelligence Theme
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

ALL_DRUGS = list(DRUG_RISK_TABLE.keys())

GENE_DRUG_MAP = {
    "CODEINE": "CYP2D6", "WARFARIN": "CYP2C9", "CLOPIDOGREL": "CYP2C19",
    "SIMVASTATIN": "SLCO1B1", "AZATHIOPRINE": "TPMT", "FLUOROURACIL": "DPYD"
}

SEV_RANK = {"none": 0, "low": 1, "moderate": 2, "high": 3, "critical": 4}

RISK_CONFIG = {
    "Safe":          {"bg": "#f0fdf4", "border": "#86efac", "accent": "#16a34a", "text": "#166534", "icon": "✓",  "badge_bg": "#dcfce7", "glow": "#22c55e20"},
    "Adjust Dosage": {"bg": "#fffbeb", "border": "#fcd34d", "accent": "#d97706", "text": "#92400e", "icon": "⚠",  "badge_bg": "#fef3c7", "glow": "#f59e0b20"},
    "Toxic":         {"bg": "#fff1f2", "border": "#fca5a5", "accent": "#dc2626", "text": "#991b1b", "icon": "✕",  "badge_bg": "#fee2e2", "glow": "#ef444420"},
    "Ineffective":   {"bg": "#f5f3ff", "border": "#c4b5fd", "accent": "#7c3aed", "text": "#5b21b6", "icon": "∅",  "badge_bg": "#ede9fe", "glow": "#8b5cf620"},
    "Unknown":       {"bg": "#f8fafc", "border": "#cbd5e1", "accent": "#64748b", "text": "#475569", "icon": "?",  "badge_bg": "#f1f5f9", "glow": "#94a3b820"},
}

SEV_CONFIG = {
    "none":     {"color": "#16a34a", "bg": "#f0fdf4", "label": "No Risk"},
    "low":      {"color": "#d97706", "bg": "#fffbeb", "label": "Low"},
    "moderate": {"color": "#ea580c", "bg": "#fff7ed", "label": "Moderate"},
    "high":     {"color": "#dc2626", "bg": "#fff1f2", "label": "High"},
    "critical": {"color": "#be123c", "bg": "#fff1f2", "label": "Critical"},
}

st.set_page_config(
    page_title="PharmaGuard — Genomic Risk Intelligence",
    page_icon="🧬", layout="wide", initial_sidebar_state="collapsed"
)

st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Bricolage+Grotesque:wght@400;500;600;700;800&family=Plus+Jakarta+Sans:ital,wght@0,300;0,400;0,500;0,600;0,700;0,800;1,400&family=JetBrains+Mono:wght@400;500;600;700&display=swap');

html, body, [class*="css"] {
    font-family: 'Plus Jakarta Sans', sans-serif !important;
    background-color: #eef2f7 !important;
    color: #1e2a3a !important;
    font-size: 15px; line-height: 1.6;
}
.stApp { background: #eef2f7 !important; }
.main .block-container { padding: 0 2.5rem 6rem !important; max-width: 1280px !important; }
#MainMenu, footer, header { visibility: hidden; }
[data-testid="stSidebar"] { background: #001a3d !important; border-right: none !important; }
[data-testid="stSidebar"] * { color: #94b4d0 !important; }

/* MASTHEAD */
.pg-masthead {
    background: #001a3d;
    margin: 0 -2.5rem 2.5rem;
    padding: 1.4rem 2.5rem;
    display: flex; align-items: center; justify-content: space-between;
    border-bottom: 3px solid #00c896;
}
.pg-brand { display: flex; align-items: center; gap: 1rem; }
.pg-logo-box {
    width: 44px; height: 44px;
    background: linear-gradient(135deg, #00c896, #0070f3);
    border-radius: 10px;
    display: flex; align-items: center; justify-content: center;
    font-size: 1.4rem; box-shadow: 0 0 30px #00c89640;
}
.pg-title {
    font-family: 'Bricolage Grotesque', sans-serif;
    font-size: 1.55rem; font-weight: 800;
    color: #ffffff; letter-spacing: -0.03em; margin: 0;
}
.pg-title span { color: #00c896; }
.pg-subtitle {
    font-family: 'JetBrains Mono', monospace;
    font-size: 0.6rem; color: #3a6080;
    letter-spacing: 0.1em; text-transform: uppercase; margin-top: 3px;
}
.pg-badges { display: flex; gap: 0.5rem; }
.pg-badge {
    font-family: 'JetBrains Mono', monospace;
    font-size: 0.6rem; font-weight: 600;
    padding: 4px 12px; border-radius: 20px;
    letter-spacing: 0.07em; text-transform: uppercase;
}
.pgb-teal  { background: #00c89620; color: #00c896; border: 1px solid #00c89640; }
.pgb-blue  { background: #0070f320; color: #60a5fa; border: 1px solid #0070f340; }
.pgb-white { background: #ffffff15; color: #94b4d0; border: 1px solid #ffffff20; }

/* STEP TRACKER */
.steps {
    display: flex; background: #fff;
    border-radius: 12px; border: 1px solid #e2e8f0;
    box-shadow: 0 1px 3px #0000000a; overflow: hidden; margin-bottom: 2rem;
}
.step { flex: 1; padding: 1rem 1.25rem; border-right: 1px solid #e2e8f0; background: #fff; }
.step:last-child { border-right: none; }
.step.on { background: #001a3d; }
.step-n {
    font-family: 'JetBrains Mono', monospace;
    font-size: 0.58rem; letter-spacing: 0.12em; color: #cbd5e1;
    margin-bottom: 0.2rem; text-transform: uppercase;
}
.step.on .step-n { color: #00c896; }
.step-t { font-size: 0.9rem; font-weight: 700; color: #cbd5e1; }
.step.on .step-t { color: #ffffff; }

/* STATS */
.stats-row {
    display: grid; grid-template-columns: repeat(5, 1fr);
    gap: 1rem; margin-bottom: 1.5rem;
}
.stat-card {
    background: #ffffff; border-radius: 12px;
    border: 1px solid #e2e8f0;
    box-shadow: 0 1px 3px #0000000a;
    padding: 1.1rem 1.25rem; position: relative; overflow: hidden;
}
.stat-card::after { content:''; position:absolute; top:0;left:0;right:0;height:3px; }
.stat-teal::after   { background: #00c896; }
.stat-blue::after   { background: #0070f3; }
.stat-amber::after  { background: #f59e0b; }
.stat-rose::after   { background: #f43f5e; }
.stat-violet::after { background: #8b5cf6; }
.stat-num {
    font-family: 'Bricolage Grotesque', sans-serif;
    font-size: 1.9rem; font-weight: 800; line-height: 1;
    margin-bottom: 0.2rem; letter-spacing: -0.03em;
}
.stat-label {
    font-family: 'JetBrains Mono', monospace;
    font-size: 0.6rem; color: #94a3b8;
    letter-spacing: 0.1em; text-transform: uppercase; margin-bottom: 0.3rem;
}
.stat-desc { font-size: 0.78rem; color: #94a3b8; line-height: 1.4; }

/* RISK CARD */
.rcard {
    border-radius: 16px; border: 1.5px solid;
    margin-bottom: 1.5rem; overflow: hidden; box-shadow: 0 4px 24px;
}
.rcard-head {
    padding: 1.25rem 1.5rem;
    display: flex; align-items: center; gap: 1rem; border-bottom: 1px solid;
}
.rcard-icon {
    width: 48px; height: 48px; border-radius: 12px;
    display: flex; align-items: center; justify-content: center;
    font-size: 1.2rem; font-weight: 700;
    font-family: 'JetBrains Mono', monospace;
    flex-shrink: 0; border: 2px solid;
}
.rcard-drug {
    font-family: 'Bricolage Grotesque', sans-serif;
    font-size: 1.25rem; font-weight: 800;
    letter-spacing: -0.02em; color: #1e2a3a;
}
.rcard-gene {
    font-family: 'JetBrains Mono', monospace;
    font-size: 0.68rem; color: #64748b; margin-top: 3px;
}
.rcard-badge {
    margin-left: auto;
    font-family: 'JetBrains Mono', monospace;
    font-size: 0.7rem; font-weight: 700;
    padding: 6px 16px; border-radius: 20px;
    letter-spacing: 0.07em; border: 1.5px solid;
}
.rcard-body { padding: 1.5rem; }

/* METRICS INSIDE CARD */
.metric-row {
    display: grid; grid-template-columns: repeat(4, 1fr);
    gap: 0.75rem; margin-bottom: 1.25rem;
}
.mc {
    background: #f8fafc; border-radius: 10px;
    border: 1px solid #e2e8f0; padding: 0.8rem 1rem;
}
.mc-label {
    font-family: 'JetBrains Mono', monospace;
    font-size: 0.6rem; color: #94a3b8;
    letter-spacing: 0.1em; text-transform: uppercase; margin-bottom: 0.3rem;
}
.mc-val {
    font-family: 'Bricolage Grotesque', sans-serif;
    font-size: 1.1rem; font-weight: 800; color: #1e2a3a; letter-spacing: -0.02em;
}

/* CONFIDENCE BAR */
.conf-outer { margin-bottom: 1.25rem; }
.conf-meta {
    display: flex; justify-content: space-between;
    font-family: 'JetBrains Mono', monospace;
    font-size: 0.63rem; color: #94a3b8;
    letter-spacing: 0.07em; margin-bottom: 0.4rem;
}
.conf-track { height: 6px; background: #e2e8f0; border-radius: 4px; overflow: hidden; }
.conf-fill  { height: 100%; border-radius: 4px; }

/* DIVIDER */
.sec-div { display: flex; align-items: center; gap: 0.8rem; margin: 1.4rem 0 1rem; }
.sec-div-line { flex: 1; height: 1px; background: #e2e8f0; }
.sec-div-label {
    font-family: 'JetBrains Mono', monospace;
    font-size: 0.62rem; color: #94a3b8;
    letter-spacing: 0.12em; text-transform: uppercase; white-space: nowrap;
}

/* VARIANTS TABLE */
.vtable { width: 100%; border-collapse: collapse; margin-bottom: 1.25rem; }
.vtable th {
    font-family: 'JetBrains Mono', monospace;
    font-size: 0.62rem; color: #94a3b8;
    letter-spacing: 0.1em; text-transform: uppercase;
    padding: 0.6rem 1rem; text-align: left;
    border-bottom: 2px solid #e2e8f0; background: #f8fafc;
}
.vtable td {
    font-family: 'JetBrains Mono', monospace;
    font-size: 0.8rem; color: #475569;
    padding: 0.65rem 1rem; border-bottom: 1px solid #f1f5f9;
}
.vtable tbody tr:hover td { background: #f8fafc; }
.vtable tbody tr:last-child td { border-bottom: none; }
.vt-rsid { color: #0070f3 !important; font-weight: 700; }
.vt-star { color: #7c3aed !important; font-weight: 700; }
.vt-nofunc  { color: #dc2626 !important; font-weight: 700; }
.vt-decfunc { color: #d97706 !important; font-weight: 700; }
.vt-incfunc { color: #0070f3 !important; font-weight: 700; }
.vt-normal  { color: #16a34a !important; font-weight: 700; }

/* RECOMMENDATION */
.rec-box { border-radius: 12px; border: 1.5px solid; padding: 1.1rem 1.25rem; margin-bottom: 1rem; }
.rec-tag {
    font-family: 'JetBrains Mono', monospace;
    font-size: 0.62rem; font-weight: 700;
    letter-spacing: 0.1em; text-transform: uppercase; margin-bottom: 0.5rem;
}
.rec-text { font-size: 0.95rem; line-height: 1.75; color: #374151; font-weight: 500; }

/* ALT DRUG TAGS */
.alt-tag {
    display: inline-block;
    background: #f0f9ff; border: 1px solid #bae6fd;
    color: #0369a1; border-radius: 20px;
    font-family: 'Plus Jakarta Sans', sans-serif;
    font-size: 0.8rem; font-weight: 600;
    padding: 3px 12px; margin: 3px;
}

/* MONITORING */
.monitoring-row {
    display: flex; align-items: flex-start; gap: 0.75rem;
    background: #f8fafc; border-radius: 10px;
    border: 1px solid #e2e8f0;
    padding: 0.9rem 1.1rem; margin-bottom: 1rem;
}
.monitoring-label {
    font-family: 'JetBrains Mono', monospace;
    font-size: 0.6rem; color: #94a3b8;
    letter-spacing: 0.1em; text-transform: uppercase; margin-bottom: 0.25rem;
}
.monitoring-text { font-size: 0.9rem; color: #4b5563; line-height: 1.6; }

/* AI EXPLANATION */
.ai-container {
    border-radius: 14px; border: 1.5px solid #e0e7ff;
    background: #fafbff; overflow: hidden; margin-bottom: 1rem;
}
.ai-header {
    background: #001a3d; padding: 0.8rem 1.25rem;
    display: flex; align-items: center; gap: 0.75rem;
}
.ai-badge {
    background: linear-gradient(135deg, #00c896, #0070f3);
    color: #fff; font-family: 'JetBrains Mono', monospace;
    font-size: 0.58rem; font-weight: 700;
    padding: 3px 10px; border-radius: 4px;
    letter-spacing: 0.08em; text-transform: uppercase;
}
.ai-header-label {
    font-family: 'JetBrains Mono', monospace;
    font-size: 0.63rem; color: #4a7090; letter-spacing: 0.08em; text-transform: uppercase;
}
.ai-section { padding: 1rem 1.25rem; border-bottom: 1px solid #e0e7ff; }
.ai-section:last-child { border-bottom: none; }
.ai-section-label {
    font-family: 'JetBrains Mono', monospace;
    font-size: 0.62rem; color: #6366f1;
    letter-spacing: 0.1em; text-transform: uppercase;
    font-weight: 700; margin-bottom: 0.45rem;
}
.ai-section-text { font-size: 0.95rem; color: #374151; line-height: 1.8; font-weight: 400; }

/* ALERT BANNERS */
.alert-critical {
    background: #fff1f2; border: 1.5px solid #fca5a5; border-left: 4px solid #dc2626;
    border-radius: 12px; padding: 1rem 1.25rem; margin-bottom: 1.25rem;
    display: flex; align-items: flex-start; gap: 0.9rem;
}
.alert-high {
    background: #fffbeb; border: 1.5px solid #fcd34d; border-left: 4px solid #f59e0b;
    border-radius: 12px; padding: 1rem 1.25rem; margin-bottom: 1.25rem;
    display: flex; align-items: flex-start; gap: 0.9rem;
}
.alert-label {
    font-family: 'JetBrains Mono', monospace;
    font-size: 0.65rem; font-weight: 700;
    letter-spacing: 0.1em; text-transform: uppercase; margin-bottom: 0.3rem;
}
.alert-body { font-size: 0.92rem; line-height: 1.65; font-weight: 500; }

/* INTERACTION CARDS */
.ix-card {
    background: #ffffff; border-radius: 12px; border: 1.5px solid;
    overflow: hidden; margin-bottom: 0.75rem;
    box-shadow: 0 1px 4px #00000008;
}
.ix-head {
    padding: 0.8rem 1.1rem; display: flex; align-items: center; gap: 0.75rem;
    border-bottom: 1px solid #f1f5f9;
}
.ix-indicator { width: 8px; height: 8px; border-radius: 50%; flex-shrink: 0; }
.ix-type {
    font-family: 'JetBrains Mono', monospace;
    font-size: 0.63rem; color: #94a3b8; letter-spacing: 0.08em; text-transform: uppercase;
}
.ix-pair { margin-left: auto; font-size: 0.9rem; font-weight: 700; color: #1e2a3a; }
.ix-sev-tag {
    font-family: 'JetBrains Mono', monospace;
    font-size: 0.62rem; font-weight: 700;
    padding: 3px 10px; border-radius: 20px;
    letter-spacing: 0.06em; text-transform: uppercase; border: 1.5px solid;
}
.ix-body { padding: 1rem 1.1rem; }
.ix-message { font-size: 0.92rem; color: #4b5563; line-height: 1.7; margin-bottom: 0.65rem; }
.ix-rec-row {
    display: flex; gap: 0.6rem; font-size: 0.85rem; color: #6b7280;
    padding-top: 0.65rem; border-top: 1px solid #f1f5f9; line-height: 1.6;
}

/* EMPTY STATE */
.empty-state {
    text-align: center; padding: 5rem 2rem;
    background: #ffffff; border-radius: 20px;
    border: 2px dashed #e2e8f0; margin-top: 0.5rem;
}
.empty-icon { font-size: 3rem; margin-bottom: 1rem; opacity: 0.3; }
.empty-title {
    font-family: 'Bricolage Grotesque', sans-serif;
    font-size: 1.2rem; font-weight: 800; color: #94a3b8; margin-bottom: 0.6rem;
    letter-spacing: -0.02em;
}
.empty-desc {
    font-family: 'JetBrains Mono', monospace;
    font-size: 0.68rem; color: #cbd5e1; letter-spacing: 0.07em; line-height: 2.2;
}

/* SUITE BANNER */
.suite-banner { border-radius: 14px; padding: 1.5rem 2rem; margin: 1rem 0; text-align: center; }
.suite-banner-title {
    font-family: 'Bricolage Grotesque', sans-serif;
    font-size: 1.6rem; font-weight: 800; letter-spacing: -0.02em; margin-bottom: 0.35rem;
}
.suite-banner-sub {
    font-family: 'JetBrains Mono', monospace;
    font-size: 0.65rem; letter-spacing: 0.1em; opacity: 0.65;
}

/* RESULT TABLE */
.result-table-wrap {
    background: #ffffff; border-radius: 12px;
    border: 1px solid #e2e8f0; overflow: hidden; margin-bottom: 1rem;
}
.rt-header {
    display: grid; grid-template-columns: 1.1fr 1.4fr 1.3fr 1.8fr 0.4fr;
    background: #f8fafc; border-bottom: 2px solid #e2e8f0;
}
.rt-hcell {
    font-family: 'JetBrains Mono', monospace;
    font-size: 0.6rem; color: #94a3b8;
    letter-spacing: 0.1em; text-transform: uppercase; padding: 0.65rem 1rem;
}
.rt-row {
    display: grid; grid-template-columns: 1.1fr 1.4fr 1.3fr 1.8fr 0.4fr;
    border-bottom: 1px solid #f1f5f9;
}
.rt-row:last-child { border-bottom: none; }
.rt-row:hover { background: #f8fafc; }
.rt-cell {
    font-family: 'JetBrains Mono', monospace;
    font-size: 0.78rem; color: #475569;
    padding: 0.7rem 1rem; display: flex; align-items: center;
}

/* GENE TAG */
.gene-tag {
    display: inline-block; background: #f1f5f9; border: 1px solid #e2e8f0;
    color: #475569; border-radius: 6px;
    font-family: 'JetBrains Mono', monospace; font-size: 0.7rem;
    padding: 2px 10px; margin: 2px;
}

/* TEST CARDS */
.test-scenario-card {
    background: #ffffff; border-radius: 12px;
    border: 1px solid #e2e8f0; box-shadow: 0 1px 3px #0000000a;
    padding: 1.1rem; min-height: 190px;
}
.test-icon { font-size: 1.3rem; margin-bottom: 0.6rem; }
.test-name { font-size: 0.85rem; font-weight: 700; color: #1e2a3a; margin-bottom: 0.4rem; line-height: 1.4; }
.test-desc {
    font-family: 'JetBrains Mono', monospace;
    font-size: 0.6rem; color: #94a3b8;
    letter-spacing: 0.05em; margin-bottom: 0.8rem; line-height: 1.6;
}
.test-result-row { display: flex; gap: 0.5rem; align-items: center; padding: 2px 0; }
.test-drug-name { font-family: 'JetBrains Mono', monospace; font-size: 0.62rem; color: #94a3b8; width: 80px; flex-shrink: 0; }
.test-result-label { font-family: 'JetBrains Mono', monospace; font-size: 0.62rem; font-weight: 700; }

/* ── STREAMLIT WIDGET OVERRIDES ── */
.stButton > button {
    background: #001a3d !important; color: #fff !important; border: none !important;
    border-radius: 10px !important;
    font-family: 'Bricolage Grotesque', sans-serif !important;
    font-weight: 700 !important; font-size: 0.95rem !important;
    padding: 0.7rem 2rem !important; letter-spacing: -0.01em !important;
    box-shadow: 0 4px 14px #001a3d25 !important;
}
.stButton > button:hover { background: #00c896 !important; box-shadow: 0 4px 20px #00c89640 !important; }
.stDownloadButton > button {
    background: #f8fafc !important; color: #475569 !important;
    border: 1.5px solid #e2e8f0 !important; border-radius: 8px !important;
    font-family: 'JetBrains Mono', monospace !important; font-size: 0.72rem !important;
}
.stDownloadButton > button:hover { border-color: #001a3d !important; color: #001a3d !important; }
.stTabs [data-baseweb="tab-list"] {
    background: #ffffff !important; border-bottom: 2px solid #e2e8f0 !important;
    border-radius: 12px 12px 0 0; gap: 0 !important;
    box-shadow: 0 1px 3px #0000000a !important; padding: 0 0.5rem !important;
}
.stTabs [data-baseweb="tab"] {
    font-family: 'JetBrains Mono', monospace !important;
    font-size: 0.72rem !important; color: #94a3b8 !important;
    letter-spacing: 0.1em !important; padding: 0.9rem 1.5rem !important;
    background: transparent !important; border: none !important;
    text-transform: uppercase !important; font-weight: 600 !important;
}
.stTabs [aria-selected="true"] {
    color: #001a3d !important; border-bottom: 2px solid #001a3d !important; font-weight: 700 !important;
}
.stTabs [data-baseweb="tab-panel"] { padding-top: 1.5rem !important; }
div[data-testid="stExpander"] {
    background: #ffffff !important; border: 1px solid #e2e8f0 !important;
    border-radius: 12px !important; box-shadow: 0 1px 3px #0000000a !important;
    margin-bottom: 0.5rem !important;
}
div[data-testid="stExpander"] summary {
    font-family: 'JetBrains Mono', monospace !important;
    font-size: 0.7rem !important; color: #475569 !important;
    letter-spacing: 0.07em !important; font-weight: 600 !important; padding: 0.9rem 1rem !important;
}
.stMultiSelect span[data-baseweb="tag"] {
    background: #e0f2fe !important; color: #0369a1 !important;
    border: 1px solid #bae6fd !important;
    font-family: 'JetBrains Mono', monospace !important; font-size: 0.7rem !important;
}
.stCheckbox label p {
    font-family: 'Plus Jakarta Sans', sans-serif !important;
    font-size: 0.92rem !important; color: #4b5563 !important; font-weight: 500 !important;
}
.stTextInput > div > div > input {
    border-radius: 8px !important; border: 1.5px solid #e2e8f0 !important;
    font-family: 'JetBrains Mono', monospace !important; font-size: 0.82rem !important;
    background: #fff !important; color: #1e2a3a !important;
}
.stTextInput > div > div > input:focus { border-color: #001a3d !important; box-shadow: 0 0 0 3px #001a3d15 !important; }
.stSelectbox > div > div > div { border-radius: 8px !important; border: 1.5px solid #e2e8f0 !important; background: #fff !important; }
.stFileUploader > div { border-radius: 10px !important; border: 2px dashed #e2e8f0 !important; background: #f8fafc !important; }
[data-testid="stMetricValue"] {
    font-family: 'Bricolage Grotesque', sans-serif !important;
    font-size: 1.8rem !important; color: #001a3d !important; font-weight: 800 !important;
}
[data-testid="stMetricLabel"] {
    font-family: 'JetBrains Mono', monospace !important;
    font-size: 0.65rem !important; color: #94a3b8 !important;
    letter-spacing: 0.08em !important; text-transform: uppercase !important;
}
.stSuccess { background: #f0fdf4 !important; border: 1px solid #86efac !important; border-radius: 10px !important; color: #166534 !important; font-family: 'Plus Jakarta Sans', sans-serif !important; }
.stError   { background: #fff1f2 !important; border: 1px solid #fca5a5 !important; border-radius: 10px !important; color: #991b1b !important; font-family: 'Plus Jakarta Sans', sans-serif !important; }
</style>
""", unsafe_allow_html=True)


# ─── Helpers ─────────────────────────────────────────────────
def load_vcf_file(filename: str) -> str:
    path = os.path.join(BASE_DIR, "sample_data", filename)
    return open(path).read() if os.path.exists(path) else get_sample_vcf()


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
    if "no_function" in s: return "vt-nofunc"
    if "decreased"  in s: return "vt-decfunc"
    if "increased"  in s: return "vt-incfunc"
    return "vt-normal"


def render_results(all_outputs, parsed_vcf, ix_report, pdf_bytes, pid):
    overall_sev = max(
        (o["risk_assessment"]["severity"] for o in all_outputs),
        key=lambda s: SEV_RANK.get(s, 0), default="none"
    )
    sc = SEV_CONFIG.get(overall_sev, SEV_CONFIG["none"])
    genes_str = " · ".join(parsed_vcf.get("detected_genes", [])) or "none detected"

    st.markdown(f"""
    <div class="stats-row">
      <div class="stat-card stat-blue">
        <div class="stat-label">Patient ID</div>
        <div class="stat-num" style="font-size:0.9rem;color:#001a3d;font-family:'JetBrains Mono',monospace;letter-spacing:0;">{pid}</div>
      </div>
      <div class="stat-card stat-teal">
        <div class="stat-label">Overall Risk</div>
        <div class="stat-num" style="color:{sc['color']};">{sc['label']}</div>
      </div>
      <div class="stat-card stat-violet">
        <div class="stat-label">Drugs Analyzed</div>
        <div class="stat-num" style="color:#7c3aed;">{len(all_outputs)}</div>
      </div>
      <div class="stat-card stat-amber">
        <div class="stat-label">Variants Found</div>
        <div class="stat-num" style="color:#d97706;">{parsed_vcf['total_variants']}</div>
      </div>
      <div class="stat-card stat-rose">
        <div class="stat-label">Genes Covered</div>
        <div class="stat-num" style="color:#f43f5e;">{len(parsed_vcf['detected_genes'])}<span style="font-size:1rem;color:#cbd5e1;">/6</span></div>
        <div class="stat-desc" style="font-size:0.65rem;">{genes_str}</div>
      </div>
    </div>
    """, unsafe_allow_html=True)

    all_json = json.dumps(all_outputs, indent=2)
    dc1, dc2, dc3 = st.columns(3)
    with dc1:
        st.download_button("⬇  All Results (JSON)", data=all_json,
                           file_name=f"pharmaguard_{pid}_all.json", mime="application/json",
                           use_container_width=True, key=f"dlall_{pid}")
    with dc2:
        if pdf_bytes:
            st.download_button("⬇  Clinical PDF Report", data=pdf_bytes,
                               file_name=f"pharmaguard_{pid}_report.pdf", mime="application/pdf",
                               use_container_width=True, key=f"dlpdf_{pid}")
    with dc3:
        if ix_report and ix_report.get("interactions_found"):
            st.download_button("⬇  Interaction Report", data=json.dumps(ix_report, indent=2),
                               file_name=f"pharmaguard_{pid}_interactions.json", mime="application/json",
                               use_container_width=True, key=f"dlix_{pid}")

    st.markdown("<div style='height:1rem'></div>", unsafe_allow_html=True)

    # Drug–Drug Interactions
    SEV_COLOR  = {"low":"#d97706","moderate":"#ea580c","high":"#dc2626","critical":"#be123c"}
    SEV_BG     = {"low":"#fffbeb","moderate":"#fff7ed","high":"#fff1f2","critical":"#fff1f2"}
    SEV_BORDER = {"low":"#fcd34d","moderate":"#fdba74","high":"#fca5a5","critical":"#fca5a5"}

    if ix_report and ix_report["interactions_found"]:
        with st.expander(
            f"⚡  Drug–Drug Interactions  ·  {ix_report['total_interactions']} Alert(s)  ·  {ix_report['overall_severity'].title()}",
            expanded=True
        ):
            for ix in ix_report["all_interactions"]:
                sev   = ix.get("severity","low")
                color = SEV_COLOR.get(sev,"#64748b")
                bg    = SEV_BG.get(sev,"#f8fafc")
                bdr   = SEV_BORDER.get(sev,"#e2e8f0")
                drugs_str = " + ".join(
                    ix.get("drugs_involved",[]) or
                    [ix.get("inhibitor_drug","")] + ix.get("affected_drugs",[])
                )
                st.markdown(f"""
                <div class="ix-card" style="border-color:{bdr};">
                  <div class="ix-head" style="background:{bg};">
                    <div class="ix-indicator" style="background:{color};"></div>
                    <span class="ix-type">{ix.get('type','').replace('_',' ').title()}</span>
                    <span class="ix-pair">{drugs_str}</span>
                    <span class="ix-sev-tag" style="color:{color};border-color:{bdr};background:#fff;">{sev.upper()}</span>
                  </div>
                  <div class="ix-body">
                    <div class="ix-message">{ix.get('message', ix.get('mechanism',''))}</div>
                    <div class="ix-rec-row"><span>💡</span><span>{ix.get('recommendation','')}</span></div>
                  </div>
                </div>""", unsafe_allow_html=True)
    elif ix_report:
        st.markdown("""
        <div style="background:#f0fdf4;border:1.5px solid #86efac;border-radius:12px;
             padding:0.9rem 1.25rem;display:flex;align-items:center;gap:0.75rem;
             margin-bottom:1rem;font-size:0.9rem;font-weight:600;color:#166534;">
          <span>✓</span>
          <span>No significant drug–drug interactions detected between selected drugs.</span>
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
        alts       = output["clinical_recommendation"].get("alternative_drugs",[])
        monitoring = output["clinical_recommendation"].get("monitoring_required","")
        exp        = output["llm_generated_explanation"]
        rc         = RISK_CONFIG.get(risk_label, RISK_CONFIG["Unknown"])
        sc         = SEV_CONFIG.get(severity, SEV_CONFIG["none"])

        if severity == "critical":
            st.markdown(f"""<div class="alert-critical">
              <span style="font-size:1.2rem;flex-shrink:0;margin-top:2px;">⛔</span>
              <div>
                <div class="alert-label" style="color:#dc2626;">Critical Safety Alert — {drug_name}</div>
                <div class="alert-body" style="color:#991b1b;">{rec}</div>
              </div>
            </div>""", unsafe_allow_html=True)
        elif severity == "high":
            st.markdown(f"""<div class="alert-high">
              <span style="font-size:1.2rem;flex-shrink:0;margin-top:2px;">⚠</span>
              <div>
                <div class="alert-label" style="color:#d97706;">Clinical Warning — {drug_name}</div>
                <div class="alert-body" style="color:#92400e;">{rec}</div>
              </div>
            </div>""", unsafe_allow_html=True)

        st.markdown(f"""
        <div class="rcard" style="border-color:{rc['border']};background:{rc['bg']};box-shadow:0 4px 24px {rc['glow']};">
          <div class="rcard-head" style="background:{rc['badge_bg']};border-bottom-color:{rc['border']};">
            <div class="rcard-icon" style="background:#fff;border-color:{rc['accent']};color:{rc['accent']};">{rc['icon']}</div>
            <div>
              <div class="rcard-drug">{drug_name.title()}</div>
              <div class="rcard-gene">{gene}  ·  {diplotype}  ·  {phenotype}</div>
            </div>
            <div class="rcard-badge" style="background:#fff;border-color:{rc['accent']};color:{rc['accent']};">{risk_label.upper()}</div>
          </div>
          <div class="rcard-body">
            <div class="metric-row">
              <div class="mc"><div class="mc-label">Phenotype</div><div class="mc-val" style="color:{rc['accent']};">{phenotype}</div></div>
              <div class="mc"><div class="mc-label">Severity</div><div class="mc-val" style="color:{sc['color']};">{sc['label']}</div></div>
              <div class="mc"><div class="mc-label">Confidence</div><div class="mc-val">{confidence:.0%}</div></div>
              <div class="mc"><div class="mc-label">Variants</div><div class="mc-val">{len(variants)}</div></div>
            </div>
            <div class="conf-outer">
              <div class="conf-meta">
                <span>Prediction Confidence</span>
                <span style="color:{rc['accent']};font-weight:700;">{confidence:.0%}</span>
              </div>
              <div class="conf-track">
                <div class="conf-fill" style="width:{confidence*100:.1f}%;background:{rc['accent']};"></div>
              </div>
            </div>
        """, unsafe_allow_html=True)

        if variants:
            rows_html = ""
            for v in variants:
                fc = func_css(v.get("functional_status",""))
                fn_label = (v.get("functional_status") or "unknown").replace("_"," ").title()
                rows_html += f"""<tr>
                  <td class="vt-rsid">{v.get('rsid','—')}</td>
                  <td class="vt-star">{v.get('star_allele','—')}</td>
                  <td>{v.get('ref','?')} → {v.get('alt','?')}</td>
                  <td>{v.get('chrom','—')}:{v.get('pos','—')}</td>
                  <td class="{fc}">{fn_label}</td>
                </tr>"""
            st.markdown(f"""
            <div class="sec-div"><div class="sec-div-line"></div>
              <div class="sec-div-label">Detected Variants ({len(variants)})</div>
              <div class="sec-div-line"></div></div>
            <table class="vtable">
              <thead><tr><th>rsID</th><th>Star Allele</th><th>Change</th><th>Position</th><th>Functional Status</th></tr></thead>
              <tbody>{rows_html}</tbody>
            </table>""", unsafe_allow_html=True)
        else:
            st.markdown("""<div style="font-family:'JetBrains Mono',monospace;font-size:0.75rem;
                 color:#94a3b8;padding:0.5rem 0 1rem;letter-spacing:0.05em;">
              ○  No variants detected — wild-type (*1/*1) assumed — normal enzyme function
            </div>""", unsafe_allow_html=True)

        style_map = {
            "Toxic":         f"background:#fff1f2;border-color:#fca5a5;",
            "Ineffective":   f"background:#f5f3ff;border-color:#c4b5fd;",
            "Adjust Dosage": f"background:#fffbeb;border-color:#fcd34d;",
            "":              f"background:#f0fdf4;border-color:#86efac;",
        }
        tag_map = {
            "Toxic": "color:#dc2626;", "Ineffective": "color:#7c3aed;",
            "Adjust Dosage": "color:#d97706;", "": "color:#16a34a;",
        }
        rk = risk_label if risk_label in style_map else ""
        st.markdown(f"""
        <div class="sec-div"><div class="sec-div-line"></div>
          <div class="sec-div-label">CPIC Dosing Recommendation</div>
          <div class="sec-div-line"></div></div>
        <div class="rec-box" style="{style_map[rk]}">
          <div class="rec-tag" style="{tag_map[rk]}">CPIC Guideline for {drug_name}</div>
          <div class="rec-text">{rec}</div>
        </div>""", unsafe_allow_html=True)

        if alts:
            tags = "".join(f'<span class="alt-tag">{a}</span>' for a in alts)
            st.markdown(f"""<div style="margin-bottom:1rem;">
              <div style="font-family:'JetBrains Mono',monospace;font-size:0.62rem;
                   color:#94a3b8;letter-spacing:0.1em;text-transform:uppercase;margin-bottom:0.5rem;">
                Alternative Drugs</div>
              {tags}
            </div>""", unsafe_allow_html=True)

        if monitoring:
            st.markdown(f"""<div class="monitoring-row">
              <span style="font-size:1rem;flex-shrink:0;">🔬</span>
              <div>
                <div class="monitoring-label">Monitoring Required</div>
                <div class="monitoring-text">{monitoring}</div>
              </div>
            </div>""", unsafe_allow_html=True)

        if exp.get("summary"):
            blocks = ""
            for lbl, key in [("Summary","summary"),("Biological Mechanism","biological_mechanism"),
                              ("Variant Significance","variant_significance"),("Clinical Implications","clinical_implications")]:
                if exp.get(key):
                    blocks += f"""<div class="ai-section">
                      <div class="ai-section-label">{lbl}</div>
                      <div class="ai-section-text">{exp[key]}</div>
                    </div>"""
            st.markdown(f"""
            <div class="sec-div"><div class="sec-div-line"></div>
              <div class="sec-div-label">AI Clinical Explanation</div>
              <div class="sec-div-line"></div></div>
            <div class="ai-container">
              <div class="ai-header">
                <span class="ai-badge">Groq LLaMA 3.3</span>
                <span class="ai-header-label">Pharmacogenomic Analysis · {drug_name}</span>
              </div>
              {blocks}
            </div>""", unsafe_allow_html=True)

        st.markdown("</div></div>", unsafe_allow_html=True)

        with st.expander(f"{{}}  Raw JSON — {drug_name}", expanded=False):
            c1, _ = st.columns([1, 3])
            with c1:
                st.download_button("⬇  Download JSON", data=json.dumps(output, indent=2),
                                   file_name=f"pharmaguard_{pid}_{drug_name}.json",
                                   mime="application/json", key=f"dl_{pid}_{drug_name}",
                                   use_container_width=True)
            st.code(json.dumps(output, indent=2), language="json")


# ══════════════════════════════════════════════════════════════
st.markdown("""
<div class="pg-masthead">
  <div class="pg-brand">
    <div class="pg-logo-box">🧬</div>
    <div>
      <div class="pg-title">Pharma<span>Guard</span></div>
      <div class="pg-subtitle">Genomic Risk Intelligence · RIFT 2026 · CPIC Aligned · Groq LLaMA 3.3</div>
    </div>
  </div>
  <div class="pg-badges">
    <span class="pg-badge pgb-teal">CPIC Aligned</span>
    <span class="pg-badge pgb-blue">RIFT 2026</span>
    <span class="pg-badge pgb-white">v2.0</span>
  </div>
</div>
""", unsafe_allow_html=True)

with st.sidebar:
    st.markdown("""<div style="padding:1.25rem 1rem 0.5rem;font-family:'JetBrains Mono',monospace;
         font-size:0.65rem;color:#4a7090;letter-spacing:0.12em;text-transform:uppercase;">Groq API Key</div>""",
        unsafe_allow_html=True)
    groq_api_key = st.text_input("", value=os.environ.get("GROQ_API_KEY",""),
                                  type="password", label_visibility="collapsed", placeholder="gsk_...")
    st.markdown("""<div style="padding:0 1rem 1rem;font-family:'JetBrains Mono',monospace;
         font-size:0.6rem;color:#2a4060;">console.groq.com → free key</div>
    <hr style="border-color:#0d2a4a;margin:0.5rem 0;">
    <div style="padding:0.75rem 1rem;font-family:'JetBrains Mono',monospace;font-size:0.62rem;
         color:#4a7090;letter-spacing:0.1em;text-transform:uppercase;margin-bottom:0.5rem;">
    Gene → Drug</div>""", unsafe_allow_html=True)
    for g, drug in [("CYP2D6","Codeine"),("CYP2C19","Clopidogrel"),("CYP2C9","Warfarin"),
                    ("SLCO1B1","Simvastatin"),("TPMT","Azathioprine"),("DPYD","Fluorouracil")]:
        st.markdown(f"""<div style="padding:4px 1rem;font-family:'JetBrains Mono',monospace;font-size:0.65rem;">
          <span style="color:#00c896;">{g}</span>
          <span style="color:#2a4060;"> → </span>
          <span style="color:#4a7090;">{drug}</span></div>""", unsafe_allow_html=True)

tab1, tab2 = st.tabs(["ANALYSIS", "TEST SUITE"])

# ─── TAB 1 ────────────────────────────────────────────────────
with tab1:
    st.markdown("""
    <div class="steps">
      <div class="step on"><div class="step-n">Step 01</div><div class="step-t">Upload VCF</div></div>
      <div class="step on"><div class="step-n">Step 02</div><div class="step-t">Select Drugs</div></div>
      <div class="step"><div class="step-n">Step 03</div><div class="step-t">Run Analysis</div></div>
      <div class="step"><div class="step-n">Step 04</div><div class="step-t">Review Results</div></div>
    </div>""", unsafe_allow_html=True)

    col_l, col_r = st.columns([1.35, 1], gap="large")

    with col_l:
        st.markdown("""<div style="font-family:'JetBrains Mono',monospace;font-size:0.65rem;
            font-weight:700;color:#001a3d;letter-spacing:0.12em;text-transform:uppercase;
            margin-bottom:0.75rem;display:flex;align-items:center;gap:0.5rem;">
            <span style="background:#001a3d;color:#00c896;width:20px;height:20px;border-radius:50%;
            display:inline-flex;align-items:center;justify-content:center;font-size:0.55rem;">1</span>
            Genomic Data Input</div>""", unsafe_allow_html=True)

        uploaded_file = st.file_uploader("Upload VCF file (max 5 MB)", type=["vcf"],
                                          help="VCF v4.2 with GENE, STAR, FUNCTION INFO tags")
        if uploaded_file is not None:
            sz = uploaded_file.size / (1024 * 1024)
            if sz > 5:
                st.error(f"File too large: {sz:.1f} MB — maximum 5 MB")
                uploaded_file = None
            else:
                peek = uploaded_file.read(500).decode("utf-8", errors="replace")
                uploaded_file.seek(0)
                if "##fileformat=VCF" not in peek and "#CHROM" not in peek:
                    st.error("Invalid VCF — must contain ##fileformat=VCF header")
                    uploaded_file = None
                else:
                    st.success(f"✓  {uploaded_file.name}  ·  {sz:.2f} MB")

        st.markdown("""<div style="font-family:'JetBrains Mono',monospace;font-size:0.65rem;
            color:#94a3b8;letter-spacing:0.1em;text-transform:uppercase;margin:1.25rem 0 0.6rem;">
            Or Use a Test Scenario</div>""", unsafe_allow_html=True)

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
        st.markdown("""<div style="font-family:'JetBrains Mono',monospace;font-size:0.65rem;
            font-weight:700;color:#001a3d;letter-spacing:0.12em;text-transform:uppercase;
            margin-bottom:0.75rem;display:flex;align-items:center;gap:0.5rem;">
            <span style="background:#001a3d;color:#00c896;width:20px;height:20px;border-radius:50%;
            display:inline-flex;align-items:center;justify-content:center;font-size:0.55rem;">2</span>
            Drug Selection</div>""", unsafe_allow_html=True)

        drug_multiselect = st.multiselect(
            "Select drugs", options=ALL_DRUGS, default=["CLOPIDOGREL"],
            format_func=lambda x: f"{x.title()}  ({GENE_DRUG_MAP.get(x,'')})",
            label_visibility="collapsed"
        )
        st.markdown("""<div style="font-family:'JetBrains Mono',monospace;font-size:0.65rem;
            color:#94a3b8;letter-spacing:0.1em;text-transform:uppercase;margin:1rem 0 0.5rem;">
            Or Type Drug Names</div>""", unsafe_allow_html=True)
        custom_drugs = st.text_input("", placeholder="e.g. CODEINE, WARFARIN, AZATHIOPRINE",
                                      label_visibility="collapsed")
        st.markdown("""<div style="font-family:'JetBrains Mono',monospace;font-size:0.65rem;
            color:#94a3b8;letter-spacing:0.1em;text-transform:uppercase;margin:1rem 0 0.5rem;">
            Patient ID</div>""", unsafe_allow_html=True)
        patient_id_input = st.text_input("", placeholder="Auto-generated if blank",
                                          label_visibility="collapsed")
        c1, c2 = st.columns(2)
        with c1: run_interactions = st.checkbox("Drug–Drug Interactions", value=True)
        with c2: generate_pdf    = st.checkbox("PDF Report", value=True)

    st.markdown("<div style='height:1rem'></div>", unsafe_allow_html=True)
    analyze_btn = st.button("▶  Run Pharmacogenomic Analysis", use_container_width=True)

    if analyze_btn:
        all_drugs = list(drug_multiselect)
        if custom_drugs.strip():
            all_drugs += [d.strip().upper() for d in custom_drugs.split(",") if d.strip()]
        all_drugs = list(set(all_drugs))

        if not all_drugs:
            st.error("Please select at least one drug to analyze.")
            st.stop()

        vcf_content = None
        if uploaded_file:
            vcf_content = uploaded_file.read().decode("utf-8", errors="replace")
        elif chosen_file:
            vcf_content = load_vcf_file(chosen_file)
        else:
            st.error("Please upload a VCF file or select a test scenario.")
            st.stop()

        pid = patient_id_input.strip() or f"PG-{str(uuid.uuid4())[:8].upper()}"

        st.markdown(f"""
        <div style="display:flex;align-items:baseline;gap:1rem;margin:2rem 0 1.5rem;
             padding-bottom:1rem;border-bottom:2px solid #e2e8f0;">
          <div style="font-family:'Bricolage Grotesque',sans-serif;font-size:1.6rem;
               font-weight:800;color:#001a3d;letter-spacing:-0.03em;">Analysis Results</div>
          <div style="font-family:'JetBrains Mono',monospace;font-size:0.75rem;
               color:#94a3b8;letter-spacing:0.08em;">{pid}</div>
        </div>""", unsafe_allow_html=True)

        with st.spinner("Running pharmacogenomic analysis..."):
            parsed_vcf, risk_results, all_outputs, ix_report, pdf_bytes = run_pipeline(
                vcf_content, all_drugs, pid, groq_api_key, run_interactions, generate_pdf
            )

        render_results(all_outputs, parsed_vcf, ix_report, pdf_bytes, pid)

        with st.expander("VCF Parsing Details", expanded=False):
            p1, p2, p3 = st.columns(3)
            p1.metric("Total Variants", parsed_vcf["total_variants"])
            p2.metric("Genes Found",    len(parsed_vcf["detected_genes"]))
            p3.metric("Parse Errors",   len(parsed_vcf["parse_errors"]))
            if parsed_vcf["detected_genes"]:
                tags = "".join(f'<span class="gene-tag">{g}</span>' for g in parsed_vcf["detected_genes"])
                st.markdown(f"<div style='margin-top:0.5rem'>{tags}</div>", unsafe_allow_html=True)
            for err in parsed_vcf.get("parse_errors", []):
                st.error(err)
    else:
        st.markdown("""
        <div class="empty-state">
          <div class="empty-icon">🧬</div>
          <div class="empty-title">Ready for Genomic Analysis</div>
          <div class="empty-desc">
            Upload a VCF file · Select drugs · Click Run Analysis<br>
            CYP2D6 · CYP2C19 · CYP2C9 · SLCO1B1 · TPMT · DPYD
          </div>
        </div>""", unsafe_allow_html=True)


# ─── TAB 2: TEST SUITE ────────────────────────────────────────
TEST_SUITE = [
    {"name":"Mixed Variants (Standard)","file":"sample.vcf","icon":"🧬",
     "drugs":["CLOPIDOGREL","CODEINE","AZATHIOPRINE"],
     "expected":{"CLOPIDOGREL":"Ineffective","CODEINE":"Adjust Dosage","AZATHIOPRINE":"Adjust Dosage"},
     "desc":"CYP2C19 *2/*3 · CYP2D6 *4/*1 · TPMT *3B/*1"},
    {"name":"UltraRapid Metabolizer — Codeine TOXIC","file":"test_ultrarapid_metabolizer.vcf","icon":"🔴",
     "drugs":["CODEINE","CLOPIDOGREL"],
     "expected":{"CODEINE":"Toxic","CLOPIDOGREL":"Safe"},
     "desc":"CYP2D6 *1xN/*1xN duplication → URM → Toxic"},
    {"name":"All Normal Wild-type — All Safe","file":"test_all_normal_wildtype.vcf","icon":"🟢",
     "drugs":ALL_DRUGS,
     "expected":{d:"Safe" for d in ALL_DRUGS},
     "desc":"Wild-type *1/*1 across all 6 genes — all Safe"},
    {"name":"Worst Case — All Poor Metabolizers","file":"test_worst_case_all_pm.vcf","icon":"🚨",
     "drugs":ALL_DRUGS,
     "expected":{"CODEINE":"Ineffective","CLOPIDOGREL":"Ineffective","WARFARIN":"Adjust Dosage",
                 "SIMVASTATIN":"Toxic","AZATHIOPRINE":"Toxic","FLUOROURACIL":"Toxic"},
     "desc":"Loss-of-function alleles across all 6 genes — stress test"},
]

RISK_COLORS_LIGHT = {
    "Safe":"#16a34a","Adjust Dosage":"#d97706",
    "Toxic":"#dc2626","Ineffective":"#7c3aed","Unknown":"#64748b",
}

with tab2:
    st.markdown("""
    <div style="margin-bottom:2rem;">
      <div style="font-family:'Bricolage Grotesque',sans-serif;font-size:1.5rem;font-weight:800;
           color:#001a3d;letter-spacing:-0.02em;margin-bottom:0.4rem;">Automated Test Suite</div>
      <div style="font-family:'JetBrains Mono',monospace;font-size:0.65rem;color:#94a3b8;letter-spacing:0.08em;">
        4 VCF scenarios · validates risk engine · expected vs actual · judge verification
      </div>
    </div>""", unsafe_allow_html=True)

    pc = st.columns(4)
    for i, sc in enumerate(TEST_SUITE):
        with pc[i]:
            exp_rows = ""
            for d, r in list(sc["expected"].items())[:4]:
                color = RISK_COLORS_LIGHT.get(r,"#64748b")
                exp_rows += f"""<div class="test-result-row">
                  <span class="test-drug-name">{d[:7]}</span>
                  <span class="test-result-label" style="color:{color};">{r}</span>
                </div>"""
            st.markdown(f"""<div class="test-scenario-card">
              <div class="test-icon">{sc['icon']}</div>
              <div class="test-name">{sc['name']}</div>
              <div class="test-desc">{sc['desc']}</div>
              {exp_rows}
            </div>""", unsafe_allow_html=True)

    st.markdown("<div style='height:1.5rem'></div>", unsafe_allow_html=True)
    tc1, tc2 = st.columns([3, 1])
    with tc1: use_llm = st.checkbox("Include LLM Explanations (uses Groq API credits, ~30s slower)", value=False)
    with tc2: run_all = st.button("▶  Run All 4 Tests", use_container_width=True)

    if run_all:
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
            if sc_pass: passed_count += 1
            else:       failed_count += 1
            all_suite_outputs.append({"scenario":sc["name"],"pass":sc_pass,"rows":rows,"outputs":ao,"file":sc["file"]})

        total  = passed_count + failed_count
        is_all = failed_count == 0
        b_col = "#16a34a" if is_all else "#d97706"
        b_bg  = "#f0fdf4" if is_all else "#fffbeb"
        b_bdr = "#86efac" if is_all else "#fcd34d"
        b_sym = "✓" if is_all else "⚠"
        st.markdown(f"""<div class="suite-banner" style="background:{b_bg};border:2px solid {b_bdr};">
          <div class="suite-banner-title" style="color:{b_col};">
            {b_sym}  {'All 4 Tests Passed' if is_all else f'{passed_count}/{total} Tests Passed'}
          </div>
          <div class="suite-banner-sub" style="color:{b_col};">
            {passed_count} passed · {failed_count} failed · {int(passed_count/total*100)}% pass rate
          </div>
        </div>""", unsafe_allow_html=True)

        for sr in all_suite_outputs:
            sc_name = sr["scenario"]
            sc_icon = next((s["icon"] for s in TEST_SUITE if s["name"]==sc_name),"🧪")
            sym = "✓" if sr["pass"] else "✕"
            with st.expander(f"{sym}  {sc_icon}  {sc_name}  —  {'PASS' if sr['pass'] else 'FAIL'}", expanded=not sr["pass"]):
                rows_html = ""
                for drug,got,exp,ok,pheno,diplo in sr["rows"]:
                    rc   = RISK_CONFIG.get(got,RISK_CONFIG["Unknown"])
                    ok_c = "#16a34a" if ok else "#dc2626"
                    ok_bg= "#f0fdf4" if ok else "#fff1f2"
                    ok_s = "✓" if ok else "✕"
                    rows_html += f"""<div class="rt-row">
                      <div class="rt-cell" style="font-weight:700;color:#1e2a3a;">{drug}</div>
                      <div class="rt-cell">
                        <span style="background:{rc['badge_bg']};color:{rc['accent']};
                             border:1.5px solid {rc['border']};border-radius:20px;
                             font-family:'JetBrains Mono',monospace;font-size:0.7rem;
                             font-weight:700;padding:3px 12px;">{got}</span>
                      </div>
                      <div class="rt-cell" style="color:#94a3b8;">{exp if exp else '—'}</div>
                      <div class="rt-cell" style="color:#64748b;">{diplo} / {pheno}</div>
                      <div class="rt-cell" style="justify-content:center;background:{ok_bg};color:{ok_c};font-weight:800;font-size:1rem;">{ok_s}</div>
                    </div>"""
                st.markdown(f"""<div class="result-table-wrap">
                  <div class="rt-header">
                    <div class="rt-hcell">Drug</div><div class="rt-hcell">Result</div>
                    <div class="rt-hcell">Expected</div><div class="rt-hcell">Diplotype / Phenotype</div>
                    <div class="rt-hcell">OK</div>
                  </div>
                  {rows_html}
                </div>""", unsafe_allow_html=True)

                dl1, dl2 = st.columns(2)
                with dl1:
                    st.download_button("⬇  Download JSON", data=json.dumps(sr["outputs"],indent=2),
                                       file_name=f"test_{sr['file'].replace('.vcf','')}.json",
                                       mime="application/json", key=f"tsc_{sc_name[:14]}", use_container_width=True)
                with dl2:
                    st.download_button("⬇  Download VCF", data=load_vcf_file(sr["file"]),
                                       file_name=sr["file"], mime="text/plain",
                                       key=f"vcf_{sc_name[:14]}", use_container_width=True)

        st.download_button(
            "⬇  Download Complete Test Suite JSON",
            data=json.dumps([{"scenario":s["scenario"],"pass":s["pass"],"results":s["outputs"]} for s in all_suite_outputs],indent=2),
            file_name=f"pharmaguard_test_suite_{datetime.now().strftime('%Y%m%d_%H%M')}.json",
            mime="application/json", use_container_width=True
        )
    else:
        st.markdown("""
        <div class="empty-state">
          <div class="empty-icon" style="font-size:2rem;">▶</div>
          <div class="empty-title">One-Click Full Validation</div>
          <div class="empty-desc">
            4 VCF scenarios · expected vs actual risk labels · pass / fail per drug<br>
            Mixed Variants · UltraRapid · All Normal · Worst Case All PM
          </div>
        </div>""", unsafe_allow_html=True)