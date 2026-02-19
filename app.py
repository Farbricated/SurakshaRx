"""
PharmaGuard v8.0 — Precision Clinical Edition
Complete UI/UX redesign: light theme, high readability, professional clinical aesthetics.
Design: "Precision Clinical" — Mayo Clinic meets Linear meets diagnostics SaaS.
Fonts: Playfair Display (display) + Plus Jakarta Sans (body) + JetBrains Mono (data)
All features from v7.0 preserved; only visual layer changed.
"""

import streamlit as st
import json, uuid, os, io
import pandas as pd
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed
from dotenv import load_dotenv

load_dotenv()

from vcf_parser import parse_vcf, get_sample_vcf
from risk_engine import run_risk_assessment, get_overall_severity, DRUG_RISK_TABLE
from llm_explainer import generate_all_explanations, generate_patient_narrative
from schema import build_output_schema
from drug_interactions import run_interaction_analysis
from pdf_report import generate_pdf_report

# ── Constants ─────────────────────────────────────────────────────────────────
BASE_DIR  = os.path.dirname(os.path.abspath(__file__))
ALL_DRUGS = list(DRUG_RISK_TABLE.keys())
GENE_DRUG_MAP = {
    "CODEINE": "CYP2D6", "WARFARIN": "CYP2C9", "CLOPIDOGREL": "CYP2C19",
    "SIMVASTATIN": "SLCO1B1", "AZATHIOPRINE": "TPMT", "FLUOROURACIL": "DPYD",
}
SEV_RANK = {"none": 0, "low": 1, "moderate": 2, "high": 3, "critical": 4}

# ── Light-theme risk config ──────────────────────────────────────────────────
RISK_CONFIG = {
    "Safe":          {"dot": "#059669", "text": "#065F46", "bg": "#ECFDF5", "border": "#A7F3D0", "emoji": "✅", "tag_bg": "#D1FAE5", "tag_text": "#065F46"},
    "Adjust Dosage": {"dot": "#D97706", "text": "#92400E", "bg": "#FFFBEB", "border": "#FDE68A", "emoji": "⚠️", "tag_bg": "#FEF3C7", "tag_text": "#92400E"},
    "Toxic":         {"dot": "#DC2626", "text": "#991B1B", "bg": "#FEF2F2", "border": "#FECACA", "emoji": "☠️", "tag_bg": "#FEE2E2", "tag_text": "#991B1B"},
    "Ineffective":   {"dot": "#7C3AED", "text": "#5B21B6", "bg": "#F5F3FF", "border": "#DDD6FE", "emoji": "❌", "tag_bg": "#EDE9FE", "tag_text": "#5B21B6"},
    "Unknown":       {"dot": "#64748B", "text": "#475569", "bg": "#F8FAFC", "border": "#E2E8F0", "emoji": "❓", "tag_bg": "#F1F5F9", "tag_text": "#475569"},
}

SEV_PALETTE = {
    "none":     {"dot": "#059669", "bg": "#ECFDF5", "border": "#A7F3D0", "text": "#065F46"},
    "low":      {"dot": "#D97706", "bg": "#FFFBEB", "border": "#FDE68A", "text": "#92400E"},
    "moderate": {"dot": "#EA580C", "bg": "#FFF7ED", "border": "#FED7AA", "text": "#9A3412"},
    "high":     {"dot": "#DC2626", "bg": "#FEF2F2", "border": "#FECACA", "text": "#991B1B"},
    "critical": {"dot": "#B91C1C", "bg": "#FFF1F1", "border": "#FCA5A5", "text": "#7F1D1D"},
}

PHENOTYPE_COLORS = {
    "PM":      {"bg": "#FEF2F2", "border": "#FECACA", "text": "#991B1B", "bar_color": "#DC2626", "label": "Poor Metabolizer",         "bar": 5},
    "IM":      {"bg": "#FFFBEB", "border": "#FDE68A", "text": "#92400E", "bar_color": "#D97706", "label": "Intermediate Metabolizer", "bar": 45},
    "NM":      {"bg": "#ECFDF5", "border": "#A7F3D0", "text": "#065F46", "bar_color": "#059669", "label": "Normal Metabolizer",       "bar": 100},
    "RM":      {"bg": "#EFF6FF", "border": "#BFDBFE", "text": "#1E40AF", "bar_color": "#3B82F6", "label": "Rapid Metabolizer",        "bar": 115},
    "URM":     {"bg": "#FFF7ED", "border": "#FED7AA", "text": "#9A3412", "bar_color": "#EA580C", "label": "Ultrarapid Metabolizer",   "bar": 130},
    "Unknown": {"bg": "#F8FAFC", "border": "#E2E8F0", "text": "#64748B", "bar_color": "#94A3B8", "label": "Unknown",                  "bar": 0},
}

POPULATION_FREQ = {
    "CYP2D6":  {"PM": 7,  "IM": 10, "NM": 77, "URM": 6},
    "CYP2C19": {"PM": 3,  "IM": 26, "NM": 52, "RM": 13, "URM": 6},
    "CYP2C9":  {"PM": 1,  "IM": 10, "NM": 89},
    "SLCO1B1": {"PM": 1,  "IM": 15, "NM": 84},
    "TPMT":    {"PM": 0.3,"IM": 10, "NM": 90},
    "DPYD":    {"PM": 0.2,"IM": 3,  "NM": 97},
}

CHROM_INFO = {
    "CYP2D6":  {"chrom": "22", "band": "q13.2",  "pos_mb": 42.5},
    "CYP2C19": {"chrom": "10", "band": "q23.33", "pos_mb": 96.7},
    "CYP2C9":  {"chrom": "10", "band": "q23.33", "pos_mb": 96.4},
    "SLCO1B1": {"chrom": "12", "band": "p12.1",  "pos_mb": 21.3},
    "TPMT":    {"chrom": "6",  "band": "p22.3",  "pos_mb": 18.1},
    "DPYD":    {"chrom": "1",  "band": "p22.1",  "pos_mb": 97.5},
}
CHROM_LEN = {"1": 248.9, "6": 170.8, "10": 133.8, "12": 133.3, "22": 50.8}

PLAIN_ENGLISH_PHENOTYPE = {
    "PM":      "Your body barely processes this medicine",
    "IM":      "Your body processes this medicine slower than average",
    "NM":      "Your body processes this medicine normally",
    "RM":      "Your body processes this medicine slightly faster than average",
    "URM":     "Your body processes this medicine dangerously fast",
    "Unknown": "Gene function unclear",
}

PLAIN_ENGLISH_RISK = {
    ("CODEINE","PM"):       "Your body can't convert codeine into a painkiller. You'd take it and feel nothing — or it could harm you.",
    ("CODEINE","URM"):      "Your body converts codeine to morphine 5× faster than normal. Even one tablet could stop your breathing.",
    ("CODEINE","IM"):       "Codeine may work less well for you. Your doctor may need to try a different painkiller.",
    ("CODEINE","NM"):       "Codeine works normally for you. Standard doses should control your pain.",
    ("WARFARIN","PM"):      "Your blood stays thin much longer than normal. Standard doses could cause dangerous bleeding.",
    ("WARFARIN","IM"):      "Warfarin lasts longer in your body than average. You'll need a lower dose.",
    ("WARFARIN","NM"):      "Warfarin works normally for you.",
    ("CLOPIDOGREL","PM"):   "This heart medication doesn't get activated in your body. It won't prevent blood clots — you need a different drug.",
    ("CLOPIDOGREL","IM"):   "This heart medication activates less than normal. You may need a stronger alternative.",
    ("CLOPIDOGREL","NM"):   "This heart medication works normally for you.",
    ("SIMVASTATIN","PM"):   "This cholesterol drug builds up in your muscles — dangerous. You need a different medication.",
    ("SIMVASTATIN","IM"):   "This cholesterol drug clears more slowly. A lower dose protects your muscles.",
    ("SIMVASTATIN","NM"):   "This cholesterol drug works normally for you.",
    ("AZATHIOPRINE","PM"):  "Your immune system drug builds up to toxic levels. Standard doses would damage your bone marrow.",
    ("AZATHIOPRINE","IM"):  "You need a lower dose of this immune drug or your bone marrow could be affected.",
    ("AZATHIOPRINE","NM"):  "This immune drug works normally for you.",
    ("FLUOROURACIL","PM"):  "Your body cannot break down this chemotherapy. Standard doses would be fatal. You need a completely different treatment.",
    ("FLUOROURACIL","IM"):  "This chemotherapy breaks down too slowly. You need half the normal dose.",
    ("FLUOROURACIL","NM"):  "This chemotherapy drug works at a normal rate in your body.",
}

PERSONAS = {
    "A": {"label": "🚨 Critical Risk",    "file": "patient_a_critical.vcf",    "drugs": ["CODEINE","FLUOROURACIL","AZATHIOPRINE"], "desc": "CYP2D6 PM · DPYD PM · TPMT PM",         "accent": "#DC2626", "bg": "#FEF2F2", "border": "#FECACA", "text": "#991B1B"},
    "B": {"label": "⚠️ Warfarin PM",      "file": "patient_b_warfarin.vcf",    "drugs": ["WARFARIN"],                              "desc": "CYP2C9 *2/*3 Poor Metabolizer",          "accent": "#D97706", "bg": "#FFFBEB", "border": "#FDE68A", "text": "#92400E"},
    "C": {"label": "💊 Drug Interaction", "file": "patient_c_interaction.vcf", "drugs": ["CLOPIDOGREL"],                           "desc": "CYP2C19 *2/*3 Poor Metabolizer",         "accent": "#7C3AED", "bg": "#F5F3FF", "border": "#DDD6FE", "text": "#5B21B6"},
    "D": {"label": "✅ All Safe",          "file": "patient_d_safe.vcf",         "drugs": ["CODEINE","WARFARIN","SIMVASTATIN"],      "desc": "Wildtype *1/*1 across all genes",        "accent": "#059669", "bg": "#ECFDF5", "border": "#A7F3D0", "text": "#065F46"},
}

TEST_SUITE = [
    {"name": "Mixed Variants",         "file": "sample.vcf",                    "drugs": ["CLOPIDOGREL","CODEINE","AZATHIOPRINE"],
     "expected": {"CLOPIDOGREL":"Ineffective","CODEINE":"Ineffective","AZATHIOPRINE":"Toxic"},
     "desc": "CYP2C19 *2/*3 · CYP2D6 *4/*4 · TPMT *3B/*3C"},
    {"name": "UltraRapid Metabolizer", "file": "test_ultrarapid_metabolizer.vcf","drugs": ["CODEINE","CLOPIDOGREL"],
     "expected": {"CODEINE":"Toxic","CLOPIDOGREL":"Safe"},
     "desc": "CYP2D6 *1xN/*1xN → URM → Codeine Toxic"},
    {"name": "All Normal Wild-type",   "file": "test_all_normal_wildtype.vcf",   "drugs": ALL_DRUGS,
     "expected": {d:"Safe" for d in ALL_DRUGS},
     "desc": "Wild-type *1/*1 across all 6 genes"},
    {"name": "Worst Case — All PM",    "file": "test_worst_case_all_pm.vcf",     "drugs": ALL_DRUGS,
     "expected": {"CODEINE":"Ineffective","CLOPIDOGREL":"Ineffective","WARFARIN":"Adjust Dosage","SIMVASTATIN":"Toxic","AZATHIOPRINE":"Toxic","FLUOROURACIL":"Toxic"},
     "desc": "Loss-of-function alleles across all 6 genes"},
]

# ── Page Config ───────────────────────────────────────────────────────────────
st.set_page_config(page_title="PharmaGuard", page_icon="🧬", layout="wide", initial_sidebar_state="collapsed")

# ══════════════════════════════════════════════════════════════════════════════
# PHASE 1–7: COMPLETE CSS SYSTEM
# ══════════════════════════════════════════════════════════════════════════════
st.markdown("""
<style>
/* ── PHASE 1: FONTS & CSS VARIABLES ───────────────────────────────────────── */
@import url('https://fonts.googleapis.com/css2?family=Playfair+Display:ital,wght@0,400;0,500;0,600;1,400;1,500&family=Plus+Jakarta+Sans:wght@300;400;500;600;700&family=JetBrains+Mono:wght@400;500;600&display=swap');

:root {
  /* Surfaces */
  --bg:          #F8F9FC;
  --surface:     #FFFFFF;
  --surface-2:   #F1F4F9;
  --surface-3:   #E8EDF5;

  /* Borders */
  --border:      #E2E8F0;
  --border-md:   #CBD5E1;
  --border-dark: #94A3B8;

  /* Text */
  --text-primary:   #0F172A;
  --text-secondary: #475569;
  --text-muted:     #94A3B8;
  --text-xmuted:    #CBD5E1;

  /* Brand */
  --brand:       #0EA5E9;
  --brand-dark:  #0284C7;
  --brand-light: #E0F2FE;
  --brand-mid:   #BAE6FD;

  /* Status */
  --safe:        #059669;
  --safe-bg:     #ECFDF5;
  --warn:        #D97706;
  --warn-bg:     #FFFBEB;
  --danger:      #DC2626;
  --danger-bg:   #FEF2F2;
  --critical:    #B91C1C;
  --violet:      #7C3AED;
  --violet-bg:   #F5F3FF;

  /* Shadows */
  --shadow-sm:   0 1px 3px rgba(15,23,42,0.06), 0 1px 2px rgba(15,23,42,0.04);
  --shadow-md:   0 4px 12px rgba(15,23,42,0.08), 0 2px 4px rgba(15,23,42,0.04);
  --shadow-lg:   0 12px 32px rgba(15,23,42,0.10), 0 4px 8px rgba(15,23,42,0.06);
  --shadow-card: 0 1px 3px rgba(15,23,42,0.05), 0 0 0 1px rgba(15,23,42,0.04);

  /* Radius */
  --r-sm:  6px;
  --r-md:  10px;
  --r-lg:  14px;
  --r-xl:  18px;
  --r-2xl: 24px;
}

/* ── PHASE 2: GLOBAL RESET ─────────────────────────────────────────────────── */
*, *::before, *::after { box-sizing: border-box; }

html, body, [class*="css"] {
  font-family: 'Plus Jakarta Sans', -apple-system, sans-serif !important;
  background: var(--bg) !important;
  color: var(--text-primary) !important;
  -webkit-font-smoothing: antialiased !important;
  text-rendering: optimizeLegibility !important;
}

.stApp { background: var(--bg) !important; }
.main .block-container {
  padding: 0 2.5rem 6rem !important;
  max-width: 1240px !important;
}
#MainMenu, footer, header { visibility: hidden; }

/* ── PHASE 2: ANIMATIONS ──────────────────────────────────────────────────── */
@keyframes fade-up   { from { opacity:0; transform:translateY(10px); } to { opacity:1; transform:translateY(0); } }
@keyframes fade-in   { from { opacity:0; } to { opacity:1; } }
@keyframes pulse-red { 0%,100%{box-shadow:0 0 0 0 rgba(220,38,38,.25)} 60%{box-shadow:0 0 0 8px rgba(220,38,38,0)} }
@keyframes shimmer   { 0%{background-position:-200% 0} 100%{background-position:200% 0} }
@keyframes bar-grow  { from{width:0} to{width:var(--w)} }
@keyframes dot-pulse { 0%,100%{opacity:1;transform:scale(1)} 50%{opacity:.6;transform:scale(1.3)} }

/* ── PHASE 2: NAVIGATION ──────────────────────────────────────────────────── */
.pg-nav {
  display: flex;
  align-items: center;
  justify-content: space-between;
  padding: 1.5rem 0 1.75rem;
  border-bottom: 1px solid var(--border);
  margin-bottom: 0;
}
.pg-logo {
  display: flex;
  align-items: baseline;
  gap: .5rem;
}
.pg-logo-text {
  font-family: 'Playfair Display', serif;
  font-size: 1.75rem;
  font-weight: 600;
  color: var(--text-primary);
  letter-spacing: -.02em;
  line-height: 1;
}
.pg-logo-text em { font-style: italic; color: var(--brand-dark); }
.pg-logo-sub {
  font-family: 'JetBrains Mono', monospace;
  font-size: .62rem;
  color: var(--text-muted);
  letter-spacing: .08em;
  text-transform: uppercase;
  font-weight: 400;
}
.pg-nav-right { display: flex; align-items: center; gap: .5rem; }
.pg-chip {
  font-family: 'JetBrains Mono', monospace;
  font-size: .6rem;
  letter-spacing: .08em;
  text-transform: uppercase;
  padding: 4px 10px;
  border-radius: 100px;
  border: 1px solid;
}
.pg-chip-default { color: var(--text-muted); border-color: var(--border); background: var(--surface); }
.pg-chip-brand   { color: var(--brand-dark);  border-color: var(--brand-mid); background: var(--brand-light); }

/* ── PHASE 2: TABS ────────────────────────────────────────────────────────── */
.stTabs [data-baseweb="tab-list"] {
  background: transparent !important;
  border-bottom: 1.5px solid var(--border) !important;
  gap: 0 !important;
  padding: 0 !important;
  margin-bottom: 2rem !important;
  box-shadow: none !important;
}
.stTabs [data-baseweb="tab"] {
  font-family: 'Plus Jakarta Sans', sans-serif !important;
  font-size: .8rem !important;
  font-weight: 500 !important;
  color: var(--text-muted) !important;
  padding: .875rem 1.5rem !important;
  background: transparent !important;
  border: none !important;
  border-bottom: 2px solid transparent !important;
  border-radius: 0 !important;
  letter-spacing: 0 !important;
  text-transform: none !important;
  transition: color .15s !important;
}
.stTabs [data-baseweb="tab"]:hover { color: var(--text-secondary) !important; }
.stTabs [aria-selected="true"] {
  color: var(--text-primary) !important;
  border-bottom-color: var(--brand) !important;
  font-weight: 600 !important;
}
.stTabs [data-baseweb="tab-panel"] { padding-top: 0 !important; }

/* ── PHASE 2: SIDEBAR ─────────────────────────────────────────────────────── */
[data-testid="stSidebar"] {
  background: var(--surface) !important;
  border-right: 1px solid var(--border) !important;
}
[data-testid="stSidebar"] * { color: var(--text-secondary) !important; }

/* ── PHASE 3: SECTION LABELS ──────────────────────────────────────────────── */
.sec-label {
  font-family: 'JetBrains Mono', monospace;
  font-size: .65rem;
  font-weight: 500;
  letter-spacing: .12em;
  text-transform: uppercase;
  color: var(--text-muted);
  margin-bottom: .75rem;
  display: flex;
  align-items: center;
  gap: .5rem;
}
.sec-label::after {
  content: '';
  flex: 1;
  height: 1px;
  background: var(--border);
}

.h-rule { border: none; border-top: 1px solid var(--border); margin: 1.25rem 0; }

/* ── PHASE 3: STEPS BAR ───────────────────────────────────────────────────── */
.steps-bar {
  display: flex;
  background: var(--surface);
  border: 1px solid var(--border);
  border-radius: var(--r-lg);
  overflow: hidden;
  margin-bottom: 2rem;
  box-shadow: var(--shadow-sm);
}
.step { flex: 1; padding: .875rem 1.25rem; border-right: 1px solid var(--border); transition: background .15s; }
.step:last-child { border-right: none; }
.step-n {
  font-family: 'JetBrains Mono', monospace;
  font-size: .58rem;
  letter-spacing: .1em;
  text-transform: uppercase;
  color: var(--brand);
  margin-bottom: 3px;
  font-weight: 500;
}
.step-l { font-size: .8rem; font-weight: 500; color: var(--text-muted); }
.step.active .step-l { color: var(--text-primary); font-weight: 600; }
.step.active { background: var(--brand-light); }
.step.active .step-n { color: var(--brand-dark); }

/* ── PHASE 3: PERSONA CARDS ───────────────────────────────────────────────── */
.persona-grid { display: grid; grid-template-columns: repeat(4,1fr); gap: .875rem; margin-bottom: 2rem; }
.persona-card {
  background: var(--surface);
  border: 1.5px solid;
  border-radius: var(--r-lg);
  padding: 1rem 1.1rem;
  cursor: pointer;
  transition: all .2s;
  box-shadow: var(--shadow-sm);
}
.persona-card:hover { transform: translateY(-2px); box-shadow: var(--shadow-md); }
.persona-label { font-size: .875rem; font-weight: 700; margin-bottom: .25rem; }
.persona-desc {
  font-family: 'JetBrains Mono', monospace;
  font-size: .6rem;
  letter-spacing: .03em;
  line-height: 1.6;
  opacity: .75;
}

/* ── PHASE 4: RISK COMMAND CENTER BANNER ─────────────────────────────────── */
.risk-banner {
  border-radius: var(--r-xl);
  padding: 1.5rem 2rem;
  margin-bottom: 1.5rem;
  border: 1.5px solid;
  animation: fade-up .35s ease;
  position: relative;
  overflow: hidden;
  box-shadow: var(--shadow-md);
}
.risk-banner::before {
  content: '';
  position: absolute;
  top: -50%;
  right: -10%;
  width: 260px;
  height: 260px;
  border-radius: 50%;
  background: currentColor;
  opacity: .04;
  pointer-events: none;
}
.rb-eyebrow {
  font-family: 'JetBrains Mono', monospace;
  font-size: .6rem;
  letter-spacing: .15em;
  text-transform: uppercase;
  opacity: .65;
  margin-bottom: .35rem;
}
.rb-headline {
  font-family: 'Playfair Display', serif;
  font-size: 2.25rem;
  font-weight: 500;
  line-height: 1.1;
  margin-bottom: .2rem;
}
.rb-stats {
  display: grid;
  grid-template-columns: repeat(4,1fr);
  gap: 1rem;
  margin-top: 1.1rem;
  padding-top: 1rem;
  border-top: 1px solid;
}
.rbs-num {
  font-family: 'Playfair Display', serif;
  font-size: 1.75rem;
  font-weight: 500;
  line-height: 1;
  margin-bottom: .1rem;
}
.rbs-key {
  font-family: 'JetBrains Mono', monospace;
  font-size: .58rem;
  letter-spacing: .1em;
  text-transform: uppercase;
  opacity: .55;
}

/* ── PHASE 4: EMERGENCY ALERT ─────────────────────────────────────────────── */
.emergency {
  background: #FFF1F1;
  border: 2px solid #FCA5A5;
  border-left: 5px solid var(--danger);
  border-radius: var(--r-lg);
  padding: 1.1rem 1.4rem;
  margin-bottom: 1rem;
  animation: pulse-red 2.5s infinite;
  box-shadow: 0 4px 16px rgba(220,38,38,.12);
}
.emergency-head { display: flex; align-items: center; gap: .65rem; margin-bottom: .4rem; }
.emergency-icon { font-size: 1.2rem; }
.emergency-drug { font-size: 1rem; font-weight: 700; color: var(--critical); letter-spacing: -.01em; }
.emergency-note { font-size: .9rem; color: #991B1B; line-height: 1.65; margin-bottom: .4rem; }
.emergency-cta {
  font-family: 'JetBrains Mono', monospace;
  font-size: .65rem;
  font-weight: 600;
  color: var(--danger);
  text-transform: uppercase;
  letter-spacing: .08em;
}

/* ── PHASE 4: GENE ACTIVITY ROW ───────────────────────────────────────────── */
.gene-row { display: grid; grid-template-columns: repeat(6,1fr); gap: .625rem; margin-bottom: 1.5rem; }
.gene-box {
  background: var(--surface);
  border: 1.5px solid;
  border-radius: var(--r-md);
  padding: .875rem .75rem;
  text-align: center;
  box-shadow: var(--shadow-sm);
  transition: box-shadow .15s, transform .15s;
}
.gene-box:hover { box-shadow: var(--shadow-md); transform: translateY(-1px); }
.gene-name {
  font-family: 'JetBrains Mono', monospace;
  font-size: .7rem;
  font-weight: 600;
  margin-bottom: .35rem;
  letter-spacing: .02em;
}
.gene-bar-track { height: 4px; border-radius: 2px; background: var(--surface-2); margin: .35rem 0; overflow: hidden; }
.gene-bar-fill { height: 100%; border-radius: 2px; transition: width .8s ease; }
.gene-pheno {
  font-family: 'JetBrains Mono', monospace;
  font-size: .6rem;
  font-weight: 500;
  letter-spacing: .03em;
}

/* ── PHASE 4: DRUG COMPARISON TABLE ──────────────────────────────────────── */
.dtable {
  background: var(--surface);
  border: 1px solid var(--border);
  border-radius: var(--r-lg);
  overflow: hidden;
  margin-bottom: 1.5rem;
  box-shadow: var(--shadow-sm);
}
.dtable-head {
  display: grid;
  grid-template-columns: 1.3fr 1.3fr 1fr 1fr 1fr 1.1fr;
  background: var(--surface-2);
  border-bottom: 1px solid var(--border);
  padding: 0 .5rem;
}
.dtable-hcell {
  font-family: 'JetBrains Mono', monospace;
  font-size: .6rem;
  letter-spacing: .1em;
  text-transform: uppercase;
  color: var(--text-muted);
  padding: .75rem .875rem;
  font-weight: 500;
}
.dtable-row {
  display: grid;
  grid-template-columns: 1.3fr 1.3fr 1fr 1fr 1fr 1.1fr;
  border-bottom: 1px solid var(--border);
  padding: 0 .5rem;
  transition: background .12s;
}
.dtable-row:last-child { border-bottom: none; }
.dtable-row:hover { background: var(--surface-2); }
.dtable-cell {
  font-family: 'Plus Jakarta Sans', sans-serif;
  font-size: .825rem;
  color: var(--text-secondary);
  padding: .75rem .875rem;
  display: flex;
  align-items: center;
}

/* ── PHASE 4: POLYGENIC RISK SCORE ───────────────────────────────────────── */
.pgx {
  background: var(--surface);
  border: 1px solid var(--border);
  border-radius: var(--r-xl);
  padding: 1.75rem 2rem;
  margin-bottom: 1.5rem;
  box-shadow: var(--shadow-md);
  position: relative;
  overflow: hidden;
}
.pgx::before {
  content: '';
  position: absolute;
  top: 0; right: 0;
  width: 200px; height: 200px;
  background: radial-gradient(circle at top right, var(--brand-light), transparent 70%);
  pointer-events: none;
}
.pgx-eye {
  font-family: 'JetBrains Mono', monospace;
  font-size: .6rem;
  letter-spacing: .15em;
  text-transform: uppercase;
  color: var(--brand-dark);
  font-weight: 500;
  margin-bottom: .35rem;
}
.pgx-score {
  font-family: 'Playfair Display', serif;
  font-size: 4rem;
  font-weight: 500;
  line-height: 1;
  letter-spacing: -.02em;
}
.pgx-label { font-size: .875rem; color: var(--text-muted); margin-bottom: 1.1rem; }
.pgx-track {
  height: 6px;
  background: var(--surface-2);
  border-radius: 3px;
  overflow: hidden;
  margin-bottom: .35rem;
}
.pgx-fill { height: 100%; border-radius: 3px; transition: width .9s cubic-bezier(.4,0,.2,1); }
.pgx-scale {
  display: flex;
  justify-content: space-between;
  font-family: 'JetBrains Mono', monospace;
  font-size: .55rem;
  color: var(--text-xmuted);
}
.pgx-pills { display: flex; flex-wrap: wrap; gap: .4rem; margin-top: 1rem; }
.pgx-pill {
  font-family: 'JetBrains Mono', monospace;
  font-size: .6rem;
  font-weight: 500;
  padding: 3px 10px;
  border-radius: 100px;
  border: 1px solid;
}

/* ── PHASE 5: DRUG x GENE HEATMAP ────────────────────────────────────────── */
.hm-wrap {
  background: var(--surface);
  border: 1px solid var(--border);
  border-radius: var(--r-lg);
  padding: 1.5rem;
  margin-bottom: 1.5rem;
  box-shadow: var(--shadow-sm);
  overflow-x: auto;
}
.hm-eye {
  font-family: 'JetBrains Mono', monospace;
  font-size: .6rem;
  letter-spacing: .12em;
  text-transform: uppercase;
  color: var(--text-muted);
  margin-bottom: 1.1rem;
  font-weight: 500;
}
.hm-grid { display: grid; gap: 4px; }
.hm-cell {
  border-radius: var(--r-sm);
  display: flex;
  flex-direction: column;
  align-items: center;
  justify-content: center;
  padding: .6rem .4rem;
  min-height: 58px;
  cursor: default;
  transition: transform .15s, box-shadow .15s;
  border: 1.5px solid;
  box-shadow: var(--shadow-sm);
}
.hm-cell:hover { transform: scale(1.07); box-shadow: var(--shadow-md); z-index: 10; position: relative; }
.hm-dname { font-family: 'JetBrains Mono', monospace; font-size: .58rem; font-weight: 600; margin-bottom: 2px; }
.hm-drisk { font-family: 'JetBrains Mono', monospace; font-size: .54rem; opacity: .8; }
.hm-header {
  font-family: 'JetBrains Mono', monospace;
  font-size: .58rem;
  letter-spacing: .06em;
  color: var(--text-muted);
  display: flex;
  align-items: center;
  justify-content: center;
  min-height: 58px;
}
.hm-legend { display: flex; gap: 1rem; margin-top: .875rem; flex-wrap: wrap; }
.hm-legend-item {
  font-family: 'JetBrains Mono', monospace;
  font-size: .58rem;
  display: flex;
  align-items: center;
  gap: 5px;
  color: var(--text-muted);
}
.hm-dot { width: 10px; height: 10px; border-radius: 3px; display: inline-block; border: 1.5px solid; }

/* ── PHASE 5: CHROMOSOME VIZ ─────────────────────────────────────────────── */
.chrom-wrap {
  background: var(--surface);
  border: 1px solid var(--border);
  border-radius: var(--r-lg);
  padding: 1.25rem 1.5rem;
  box-shadow: var(--shadow-sm);
}
.chrom-eye {
  font-family: 'JetBrains Mono', monospace;
  font-size: .6rem;
  letter-spacing: .1em;
  text-transform: uppercase;
  color: var(--text-muted);
  margin-bottom: .75rem;
  font-weight: 500;
}
.chrom-row { display: flex; align-items: center; gap: .75rem; margin-bottom: .45rem; }
.chrom-lbl {
  font-family: 'JetBrains Mono', monospace;
  font-size: .62rem;
  color: var(--text-muted);
  width: 20px;
  text-align: right;
  flex-shrink: 0;
}
.chrom-bar {
  flex: 1;
  height: 13px;
  background: var(--surface-2);
  border-radius: 7px;
  position: relative;
  overflow: visible;
  border: 1px solid var(--border);
}
.chrom-body {
  position: absolute;
  inset: 0;
  background: linear-gradient(90deg, #E2E8F0, #EFF6FF, #E2E8F0);
  border-radius: 7px;
}
.chrom-marker {
  position: absolute;
  top: -4px;
  width: 3px;
  height: 21px;
  border-radius: 2px;
  transform: translateX(-50%);
  box-shadow: 0 0 6px currentColor;
}
.chrom-gene-lbl {
  font-family: 'JetBrains Mono', monospace;
  font-size: .58rem;
  color: var(--text-secondary);
  width: 60px;
  flex-shrink: 0;
  font-weight: 500;
}
.chrom-band {
  font-family: 'JetBrains Mono', monospace;
  font-size: .55rem;
  color: var(--text-muted);
}

/* ── PHASE 5: POPULATION FREQ ────────────────────────────────────────────── */
.pop-wrap {
  background: var(--surface);
  border: 1px solid var(--border);
  border-radius: var(--r-md);
  padding: 1rem 1.25rem;
  margin-bottom: .875rem;
  box-shadow: var(--shadow-sm);
}
.pop-eye {
  font-family: 'JetBrains Mono', monospace;
  font-size: .6rem;
  letter-spacing: .1em;
  text-transform: uppercase;
  color: var(--text-muted);
  margin-bottom: .6rem;
  font-weight: 500;
}
.pop-row { display: flex; align-items: center; gap: .75rem; margin-bottom: .4rem; }
.pop-ph {
  font-family: 'JetBrains Mono', monospace;
  font-size: .65rem;
  color: var(--text-secondary);
  width: 100px;
  flex-shrink: 0;
  font-weight: 500;
}
.pop-track { flex: 1; height: 5px; background: var(--surface-2); border-radius: 3px; overflow: hidden; }
.pop-fill { height: 100%; border-radius: 3px; }
.pop-pct {
  font-family: 'JetBrains Mono', monospace;
  font-size: .62rem;
  width: 35px;
  text-align: right;
  color: var(--text-muted);
}
.pop-you {
  font-family: 'JetBrains Mono', monospace;
  font-size: .56rem;
  color: var(--brand-dark);
  font-weight: 600;
  margin-left: 3px;
}

/* ── PHASE 5: INTERACTION MATRIX ─────────────────────────────────────────── */
.ix-matrix-grid { display: grid; gap: 4px; }
.ix-cell {
  border-radius: var(--r-sm);
  display: flex;
  align-items: center;
  justify-content: center;
  min-height: 46px;
  font-family: 'JetBrains Mono', monospace;
  font-size: .58rem;
  text-align: center;
  padding: .3rem;
  transition: transform .12s;
  border: 1.5px solid;
  cursor: pointer;
  font-weight: 600;
}
.ix-cell:hover { transform: scale(1.06); z-index: 5; position: relative; box-shadow: var(--shadow-md); }
.ix-head {
  font-family: 'JetBrains Mono', monospace;
  font-size: .58rem;
  letter-spacing: .06em;
  color: var(--text-muted);
  display: flex;
  align-items: center;
  justify-content: center;
  min-height: 46px;
}

/* ── PHASE 6: DRUG RESULT CARDS ──────────────────────────────────────────── */
.rcard {
  background: var(--surface);
  border: 1px solid var(--border);
  border-radius: var(--r-xl);
  margin-bottom: 1.25rem;
  overflow: hidden;
  animation: fade-up .3s ease;
  box-shadow: var(--shadow-card);
  transition: box-shadow .2s;
}
.rcard:hover { box-shadow: var(--shadow-md); }
.rcard-top {
  padding: 1.25rem 1.5rem;
  display: flex;
  align-items: center;
  justify-content: space-between;
  border-bottom: 1px solid var(--border);
}
.rcard-left { display: flex; align-items: center; gap: .875rem; }
.rcard-dot { width: 11px; height: 11px; border-radius: 50%; flex-shrink: 0; }
.rcard-name {
  font-size: 1.05rem;
  font-weight: 700;
  letter-spacing: -.02em;
  color: var(--text-primary);
  line-height: 1;
}
.rcard-meta {
  font-family: 'JetBrains Mono', monospace;
  font-size: .65rem;
  color: var(--text-muted);
  margin-top: 3px;
  letter-spacing: .03em;
  font-weight: 400;
}
.rcard-badge {
  font-family: 'Plus Jakarta Sans', sans-serif;
  font-size: .75rem;
  font-weight: 700;
  letter-spacing: -.01em;
  padding: 5px 14px;
  border-radius: 100px;
  border: 1.5px solid;
}
.rcard-body { padding: 1.25rem 1.5rem; }

/* Metric cells */
.mc-row {
  display: grid;
  grid-template-columns: repeat(4,1fr);
  gap: 1px;
  background: var(--border);
  border-radius: var(--r-md);
  overflow: hidden;
  margin-bottom: 1.25rem;
  border: 1px solid var(--border);
}
.mc-cell { background: var(--surface-2); padding: .875rem 1rem; }
.mc-key {
  font-family: 'JetBrains Mono', monospace;
  font-size: .58rem;
  letter-spacing: .1em;
  text-transform: uppercase;
  color: var(--text-muted);
  margin-bottom: .3rem;
  font-weight: 500;
}
.mc-val {
  font-size: 1rem;
  font-weight: 700;
  color: var(--text-primary);
  letter-spacing: -.01em;
}

/* Confidence bars */
.conf-row { display: grid; grid-template-columns: 1fr 1fr; gap: 1rem; margin-bottom: 1.25rem; }
.conf-lbl {
  font-family: 'JetBrains Mono', monospace;
  font-size: .6rem;
  letter-spacing: .08em;
  text-transform: uppercase;
  color: var(--text-muted);
  margin-bottom: 5px;
  display: flex;
  justify-content: space-between;
  align-items: center;
  font-weight: 500;
}
.conf-track {
  height: 4px;
  background: var(--surface-2);
  border-radius: 2px;
  overflow: hidden;
  border: 1px solid var(--border);
}
.conf-fill { height: 100%; border-radius: 2px; transition: width .7s ease; }

/* Variant table */
.vtable { width: 100%; border-collapse: collapse; }
.vtable th {
  font-family: 'JetBrains Mono', monospace;
  font-size: .62rem;
  letter-spacing: .09em;
  text-transform: uppercase;
  color: var(--text-muted);
  padding: 0 .75rem .6rem;
  text-align: left;
  border-bottom: 1px solid var(--border);
  font-weight: 500;
}
.vtable td {
  font-family: 'JetBrains Mono', monospace;
  font-size: .78rem;
  color: var(--text-secondary);
  padding: .6rem .75rem;
  border-bottom: 1px solid var(--border);
  font-weight: 400;
}
.vtable tbody tr:last-child td { border-bottom: none; }
.vtable tbody tr:hover td { background: var(--surface-2); }
.v-rsid   { color: #2563EB !important; font-weight: 500 !important; }
.v-star   { color: #7C3AED !important; font-weight: 500 !important; }
.v-nofunc { color: var(--danger) !important; font-weight: 500 !important; }
.v-dec    { color: var(--warn) !important; font-weight: 500 !important; }
.v-norm   { color: var(--safe) !important; font-weight: 500 !important; }

/* Rec box */
.rec-box {
  border-radius: var(--r-md);
  border: 1.5px solid;
  padding: 1rem 1.25rem;
  margin-bottom: 1rem;
}
.rec-lbl {
  font-family: 'JetBrains Mono', monospace;
  font-size: .62rem;
  letter-spacing: .09em;
  text-transform: uppercase;
  font-weight: 600;
  margin-bottom: .4rem;
}
.rec-txt { font-size: .9rem; line-height: 1.75; color: var(--text-secondary); }

/* Alt chips */
.alt-chips { display: flex; flex-wrap: wrap; gap: .4rem; }
.alt-chip {
  font-family: 'JetBrains Mono', monospace;
  font-size: .7rem;
  font-weight: 500;
  color: var(--brand-dark);
  background: var(--brand-light);
  border: 1px solid var(--brand-mid);
  border-radius: 100px;
  padding: 4px 12px;
}

/* CPIC badge */
.cpic {
  font-family: 'JetBrains Mono', monospace;
  font-size: .58rem;
  letter-spacing: .05em;
  font-weight: 600;
  background: #FEF9C3;
  border: 1px solid #FDE047;
  color: #713F12;
  padding: 2px 8px;
  border-radius: 4px;
  display: inline-block;
  margin-left: .5rem;
  vertical-align: middle;
}

/* AI narrative */
.ai-narrative {
  background: linear-gradient(135deg, #F8FBFF, #F0F7FF);
  border: 1.5px solid #BFDBFE;
  border-radius: var(--r-xl);
  padding: 1.5rem;
  margin-bottom: 1.5rem;
  box-shadow: 0 2px 8px rgba(14,165,233,.06);
}
.ai-nar-head { display: flex; align-items: center; gap: .75rem; margin-bottom: .875rem; }
.ai-nar-badge {
  font-family: 'JetBrains Mono', monospace;
  font-size: .6rem;
  letter-spacing: .08em;
  text-transform: uppercase;
  background: var(--brand-light);
  border: 1px solid var(--brand-mid);
  color: var(--brand-dark);
  padding: 3px 10px;
  border-radius: 4px;
  font-weight: 600;
}
.ai-nar-title { font-size: .95rem; font-weight: 600; color: var(--brand-dark); }
.ai-nar-text { font-size: .95rem; line-height: 1.85; color: var(--text-secondary); }
.ai-block {
  border: 1px solid var(--border);
  border-radius: var(--r-lg);
  overflow: hidden;
  margin-top: 1.25rem;
  box-shadow: var(--shadow-sm);
}
.ai-block-head {
  padding: .65rem 1rem;
  background: var(--surface-2);
  border-bottom: 1px solid var(--border);
  display: flex;
  align-items: center;
  gap: .65rem;
}
.ai-badge {
  font-family: 'JetBrains Mono', monospace;
  font-size: .6rem;
  letter-spacing: .08em;
  text-transform: uppercase;
  color: var(--text-muted);
  background: var(--surface-3);
  padding: 3px 9px;
  border-radius: 4px;
  font-weight: 500;
  border: 1px solid var(--border);
}
.ai-section { padding: .9rem 1rem; border-bottom: 1px solid var(--border); }
.ai-section:last-child { border-bottom: none; }
.ai-section:hover { background: var(--surface-2); }
.ai-section-lbl {
  font-family: 'JetBrains Mono', monospace;
  font-size: .6rem;
  letter-spacing: .09em;
  text-transform: uppercase;
  color: var(--text-muted);
  margin-bottom: .4rem;
  font-weight: 500;
}
.ai-section-txt { font-size: .9rem; line-height: 1.8; color: var(--text-secondary); }

/* ── PHASE 6: CLINICAL NOTE ──────────────────────────────────────────────── */
.note-box {
  background: var(--surface-2);
  border: 1px solid var(--border);
  border-radius: var(--r-lg);
  padding: 1.25rem 1.5rem;
}
.note-box pre {
  font-family: 'JetBrains Mono', monospace;
  font-size: .8rem;
  color: var(--text-secondary);
  line-height: 1.8;
  white-space: pre-wrap;
  word-break: break-word;
  font-weight: 400;
}

/* ── PHASE 6: BEFORE / AFTER ─────────────────────────────────────────────── */
.ba-wrap {
  display: grid;
  grid-template-columns: 1fr 1fr;
  border: 1px solid var(--border);
  border-radius: var(--r-xl);
  overflow: hidden;
  margin-bottom: 1.5rem;
  box-shadow: var(--shadow-md);
}
.ba-side { padding: 1.5rem; }
.ba-label {
  font-family: 'JetBrains Mono', monospace;
  font-size: .62rem;
  letter-spacing: .1em;
  text-transform: uppercase;
  font-weight: 600;
  margin-bottom: .75rem;
}
.ba-drug { font-size: 1.05rem; font-weight: 700; margin-bottom: .3rem; }
.ba-outcome { font-size: .875rem; line-height: 1.65; }

/* ── PHASE 6: PRESCRIPTION CHECKER ──────────────────────────────────────── */
.rx-result {
  border-radius: var(--r-lg);
  padding: 1.25rem 1.5rem;
  margin-top: .875rem;
  border: 1.5px solid;
  animation: fade-up .25s ease;
  box-shadow: var(--shadow-md);
}
.rx-verdict {
  font-size: .875rem;
  font-weight: 700;
  margin-bottom: .5rem;
  letter-spacing: -.01em;
}
.rx-detail { font-size: .875rem; line-height: 1.7; color: var(--text-secondary); margin-bottom: .5rem; }
.rx-meta {
  font-family: 'JetBrains Mono', monospace;
  font-size: .6rem;
  letter-spacing: .06em;
  text-transform: uppercase;
}

/* ── PHASE 6: PATIENT MODE ───────────────────────────────────────────────── */
.patient-banner {
  border-radius: var(--r-xl);
  padding: 1.5rem;
  margin-bottom: 1.5rem;
  border: 2px solid;
  box-shadow: var(--shadow-md);
}
.patient-banner-title { font-size: 1.15rem; font-weight: 700; margin-bottom: .35rem; }
.patient-banner-sub { font-size: .9rem; line-height: 1.7; }
.patient-card {
  background: var(--surface);
  border: 1.5px solid;
  border-radius: var(--r-xl);
  padding: 1.5rem;
  margin-bottom: 1rem;
  box-shadow: var(--shadow-card);
  animation: fade-up .3s ease;
  transition: box-shadow .2s;
}
.patient-card:hover { box-shadow: var(--shadow-md); }
.patient-drug { font-size: 1.2rem; font-weight: 700; margin-bottom: .2rem; letter-spacing: -.02em; }
.patient-verdict { font-size: 1rem; font-weight: 600; line-height: 1.5; margin-bottom: .5rem; }
.patient-gene-tag {
  font-family: 'JetBrains Mono', monospace;
  font-size: .6rem;
  letter-spacing: .06em;
  color: var(--text-muted);
  margin-bottom: .6rem;
  font-weight: 500;
}
.patient-plain { font-size: .9rem; line-height: 1.8; color: var(--text-secondary); }
.patient-action {
  display: flex;
  align-items: flex-start;
  gap: .65rem;
  background: var(--surface-2);
  border: 1px solid var(--border);
  border-radius: var(--r-md);
  padding: .875rem 1rem;
  margin-top: .875rem;
}
.patient-action-text { font-size: .875rem; color: var(--text-primary); line-height: 1.65; }

/* ── PHASE 6: TEST SUITE ─────────────────────────────────────────────────── */
.test-card {
  background: var(--surface);
  border: 1px solid var(--border);
  border-radius: var(--r-lg);
  padding: 1.2rem;
  box-shadow: var(--shadow-sm);
  transition: box-shadow .15s;
}
.test-card:hover { box-shadow: var(--shadow-md); }
.test-name { font-size: .9rem; font-weight: 700; margin-bottom: .3rem; color: var(--text-primary); }
.test-desc {
  font-family: 'JetBrains Mono', monospace;
  font-size: .6rem;
  color: var(--text-muted);
  margin-bottom: .875rem;
  line-height: 1.7;
  letter-spacing: .02em;
}

/* ── PHASE 7: BUTTONS & INPUTS ───────────────────────────────────────────── */
.stButton > button {
  background: var(--text-primary) !important;
  color: #FFFFFF !important;
  border: none !important;
  border-radius: var(--r-md) !important;
  font-family: 'Plus Jakarta Sans', sans-serif !important;
  font-weight: 600 !important;
  font-size: .875rem !important;
  padding: .7rem 1.75rem !important;
  letter-spacing: -.01em !important;
  transition: all .15s !important;
  box-shadow: 0 1px 2px rgba(15,23,42,.15) !important;
}
.stButton > button:hover {
  background: #1E293B !important;
  box-shadow: 0 4px 12px rgba(15,23,42,.2) !important;
  transform: translateY(-1px) !important;
}
.stButton > button:active { transform: translateY(0) !important; }

.stDownloadButton > button {
  background: var(--surface) !important;
  color: var(--text-secondary) !important;
  border: 1px solid var(--border) !important;
  border-radius: var(--r-md) !important;
  font-family: 'JetBrains Mono', monospace !important;
  font-size: .72rem !important;
  letter-spacing: .02em !important;
  padding: .5rem 1rem !important;
  transition: all .15s !important;
  box-shadow: var(--shadow-sm) !important;
}
.stDownloadButton > button:hover {
  background: var(--text-primary) !important;
  color: #FFFFFF !important;
  border-color: var(--text-primary) !important;
  box-shadow: var(--shadow-md) !important;
}

/* File uploader */
.stFileUploader > div {
  border-radius: var(--r-lg) !important;
  border: 1.5px dashed var(--border-md) !important;
  background: var(--surface) !important;
  transition: border-color .15s, background .15s !important;
}
.stFileUploader > div:hover {
  border-color: var(--brand) !important;
  background: var(--brand-light) !important;
}

/* Inputs */
.stTextInput > div > div > input {
  border-radius: var(--r-md) !important;
  border: 1.5px solid var(--border-md) !important;
  background: var(--surface) !important;
  color: var(--text-primary) !important;
  font-family: 'Plus Jakarta Sans', sans-serif !important;
  font-size: .875rem !important;
  padding: .6rem .875rem !important;
  transition: border-color .15s, box-shadow .15s !important;
  box-shadow: var(--shadow-sm) !important;
}
.stTextInput > div > div > input:focus {
  border-color: var(--brand) !important;
  box-shadow: 0 0 0 3px rgba(14,165,233,.15) !important;
  outline: none !important;
}

/* Select / Multiselect */
.stSelectbox [data-baseweb="select"] > div,
.stMultiSelect [data-baseweb="select"] > div {
  border-radius: var(--r-md) !important;
  border: 1.5px solid var(--border-md) !important;
  background: var(--surface) !important;
  box-shadow: var(--shadow-sm) !important;
}
.stMultiSelect span[data-baseweb="tag"] {
  background: var(--brand-light) !important;
  color: var(--brand-dark) !important;
  border: 1px solid var(--brand-mid) !important;
  font-family: 'JetBrains Mono', monospace !important;
  font-size: .7rem !important;
  border-radius: 4px !important;
}

/* Checkboxes */
.stCheckbox label p {
  font-family: 'Plus Jakarta Sans', sans-serif !important;
  font-size: .875rem !important;
  color: var(--text-secondary) !important;
}

/* Toggle */
.stToggle label p {
  font-family: 'Plus Jakarta Sans', sans-serif !important;
  font-size: .875rem !important;
  color: var(--text-secondary) !important;
}

/* Expander */
div[data-testid="stExpander"] {
  background: var(--surface) !important;
  border: 1px solid var(--border) !important;
  border-radius: var(--r-md) !important;
  box-shadow: var(--shadow-sm) !important;
  margin-bottom: .5rem !important;
}
div[data-testid="stExpander"] summary {
  font-family: 'Plus Jakarta Sans', sans-serif !important;
  font-size: .8rem !important;
  font-weight: 500 !important;
  color: var(--text-secondary) !important;
  padding: .8rem 1rem !important;
  letter-spacing: 0 !important;
}
div[data-testid="stExpander"] summary:hover {
  color: var(--text-primary) !important;
}

/* Metrics */
[data-testid="stMetric"] {
  background: var(--surface) !important;
  border: 1px solid var(--border) !important;
  border-radius: var(--r-lg) !important;
  padding: 1rem 1.25rem !important;
  box-shadow: var(--shadow-sm) !important;
}
[data-testid="stMetricLabel"] { font-family: 'JetBrains Mono', monospace !important; font-size: .65rem !important; color: var(--text-muted) !important; text-transform: uppercase !important; letter-spacing: .08em !important; }
[data-testid="stMetricValue"] { font-family: 'Playfair Display', serif !important; font-size: 1.75rem !important; color: var(--text-primary) !important; }

/* Spinner */
.stSpinner > div { border-color: var(--brand) transparent transparent transparent !important; }

/* ── PHASE 7: EMPTY STATE ─────────────────────────────────────────────────── */
.empty-state {
  text-align: center;
  padding: 5rem 2rem;
  border: 1.5px dashed var(--border-md);
  border-radius: var(--r-2xl);
  background: var(--surface);
  box-shadow: var(--shadow-sm);
}
.empty-icon {
  font-size: 2.5rem;
  margin-bottom: 1rem;
  display: block;
  opacity: .35;
}
.empty-title { font-family: 'Playfair Display', serif; font-size: 1.3rem; color: var(--text-secondary); margin-bottom: .5rem; font-weight: 400; }
.empty-hint {
  font-family: 'JetBrains Mono', monospace;
  font-size: .65rem;
  color: var(--text-muted);
  letter-spacing: .06em;
  line-height: 2.2;
}

/* ── PHASE 7: SUCCESS / INFO BANNERS ─────────────────────────────────────── */
.info-strip {
  display: flex;
  align-items: flex-start;
  gap: .75rem;
  background: var(--brand-light);
  border: 1px solid var(--brand-mid);
  border-radius: var(--r-md);
  padding: .875rem 1rem;
  margin-bottom: 1rem;
}
.info-strip-text { font-size: .875rem; color: var(--brand-dark); line-height: 1.65; }

/* Code blocks */
.stCode { border-radius: var(--r-md) !important; }
pre { font-family: 'JetBrains Mono', monospace !important; }
</style>
""", unsafe_allow_html=True)


# ══════════════════════════════════════════════════════════════════════════════
# HELPERS
# ══════════════════════════════════════════════════════════════════════════════

def load_vcf_file(filename):
    path = os.path.join(BASE_DIR, "sample_data", filename)
    return open(path).read() if os.path.exists(path) else get_sample_vcf()


def run_pipeline(vcf_content, drugs, pid, groq_key, run_ix=True, gen_pdf=True, skip_llm=False):
    parsed  = parse_vcf(vcf_content)
    results = run_risk_assessment(parsed, drugs)
    results = generate_all_explanations(groq_key, results, skip_llm=skip_llm)
    outputs = [
        build_output_schema(patient_id=pid, drug=r["drug"], result=r,
                            parsed_vcf=parsed, llm_exp=r.get("llm_explanation", {}))
        for r in results
    ]
    ix  = run_interaction_analysis(drugs, results) if run_ix and len(drugs) > 1 else None
    pdf = None
    if gen_pdf:
        try:
            pdf = generate_pdf_report(pid, outputs, parsed)
        except Exception:
            pass
    return parsed, results, outputs, ix, pdf


def func_cls(status):
    s = (status or "").lower()
    if "no_function" in s: return "v-nofunc"
    if "decreased"   in s: return "v-dec"
    return "v-norm"


def sec(label):
    st.markdown(f'<div class="sec-label">{label}</div>', unsafe_allow_html=True)


# ══════════════════════════════════════════════════════════════════════════════
# COMPONENT: POLYGENIC RISK SCORE
# ══════════════════════════════════════════════════════════════════════════════

def compute_pgx_score(all_outputs):
    SEV_S  = {"none":0,"low":20,"moderate":45,"high":70,"critical":100}
    RISK_S = {"Safe":0,"Adjust Dosage":35,"Toxic":85,"Ineffective":70,"Unknown":20}
    W      = {"FLUOROURACIL":1.4,"AZATHIOPRINE":1.3,"CLOPIDOGREL":1.3,"WARFARIN":1.2,"CODEINE":1.1,"SIMVASTATIN":1.0}
    if not all_outputs: return 0,"No data",[]
    tw = ws = 0; breakdown = []
    for o in all_outputs:
        drug=o["drug"]; sev=o["risk_assessment"]["severity"]; rl=o["risk_assessment"]["risk_label"]
        gene=o["pharmacogenomic_profile"]["primary_gene"]; ph=o["pharmacogenomic_profile"]["phenotype"]
        sc=(SEV_S.get(sev,0)+RISK_S.get(rl,0))/2; wt=W.get(drug,1.0)
        ws+=sc*wt; tw+=wt; breakdown.append((gene,drug,ph,rl,sc))
    final=min(100,int(ws/tw)) if tw else 0
    label=["Low Risk","Moderate Risk","High Risk","Very High Risk","Critical Risk"][min(4,final//20)]
    return final,label,breakdown


def render_pgx_score(all_outputs):
    score, label, breakdown = compute_pgx_score(all_outputs)
    clrs = ["#059669","#D97706","#EA580C","#DC2626","#991B1B"]
    color = clrs[min(4, score//20)]
    pills = ""
    for gene,_,ph,rl,_ in breakdown:
        rc = RISK_CONFIG.get(rl, RISK_CONFIG["Unknown"])
        pills += f'<span class="pgx-pill" style="background:{rc["bg"]};border-color:{rc["border"]};color:{rc["text"]};">{gene} · {ph}</span>'
    st.markdown(f"""
    <div class="pgx">
      <div class="pgx-eye">Polygenic Risk Score</div>
      <div class="pgx-score" style="color:{color};">{score}</div>
      <div class="pgx-label">{label} — composite across {len(all_outputs)} drug{"s" if len(all_outputs)!=1 else ""}</div>
      <div class="pgx-track">
        <div class="pgx-fill" style="width:{score}%;background:linear-gradient(90deg,{color}bb,{color});"></div>
      </div>
      <div class="pgx-scale"><span>0 — No Risk</span><span>50 — High</span><span>100 — Critical</span></div>
      <div class="pgx-pills">{pills}</div>
    </div>""", unsafe_allow_html=True)


# ══════════════════════════════════════════════════════════════════════════════
# COMPONENT: RISK COMMAND CENTER BANNER
# ══════════════════════════════════════════════════════════════════════════════

def render_risk_banner(all_outputs, parsed):
    sev = max((o["risk_assessment"]["severity"] for o in all_outputs),
              key=lambda s: SEV_RANK.get(s,0), default="none")
    sp  = SEV_PALETTE.get(sev, SEV_PALETTE["none"])
    EMO   = {"none":"✅","low":"💛","moderate":"🟠","high":"🔴","critical":"🚨"}
    LABEL = {"none":"All Clear","low":"Low","moderate":"Moderate","high":"High","critical":"Critical"}
    emoji=EMO.get(sev,""); label=LABEL.get(sev,"Unknown")
    hc   = sum(1 for o in all_outputs if o["risk_assessment"]["severity"] in ("high","critical"))
    st.markdown(f"""
    <div class="risk-banner" style="background:{sp['bg']};border-color:{sp['border']};color:{sp['text']};">
      <div class="rb-eyebrow">Risk Command Center</div>
      <div class="rb-headline">{emoji} {label} Risk</div>
      <div class="rb-stats" style="border-color:{sp['border']}88;">
        <div><div class="rbs-num">{len(all_outputs)}</div><div class="rbs-key">Drugs Analysed</div></div>
        <div><div class="rbs-num" style="{'color:#DC2626' if hc else ''}">{hc}</div><div class="rbs-key">High / Critical</div></div>
        <div><div class="rbs-num">{len(parsed.get('detected_genes',[]))}</div><div class="rbs-key">Genes Detected</div></div>
        <div><div class="rbs-num">{parsed.get('total_variants',0)}</div><div class="rbs-key">Variants Found</div></div>
      </div>
    </div>""", unsafe_allow_html=True)


# ══════════════════════════════════════════════════════════════════════════════
# COMPONENT: EMERGENCY ALERTS
# ══════════════════════════════════════════════════════════════════════════════

def render_emergency_alerts(all_outputs):
    for o in all_outputs:
        if o["risk_assessment"]["severity"] == "critical":
            drug = o["drug"]
            note = o["clinical_recommendation"]["dosing_recommendation"][:220]
            st.markdown(f"""
            <div class="emergency">
              <div class="emergency-head">
                <span class="emergency-icon">🚨</span>
                <span class="emergency-drug">Critical Alert — {drug}</span>
              </div>
              <div class="emergency-note">{note}{"…" if len(o["clinical_recommendation"]["dosing_recommendation"])>220 else ""}</div>
              <div class="emergency-cta">⚡ Contact prescribing physician immediately</div>
            </div>""", unsafe_allow_html=True)


# ══════════════════════════════════════════════════════════════════════════════
# COMPONENT: GENE ACTIVITY HEATMAP
# ══════════════════════════════════════════════════════════════════════════════

def render_gene_heatmap(all_outputs):
    GENE_ORDER = ["CYP2D6","CYP2C19","CYP2C9","SLCO1B1","TPMT","DPYD"]
    gp = {o["pharmacogenomic_profile"]["primary_gene"]: o["pharmacogenomic_profile"]["phenotype"]
          for o in all_outputs}
    boxes = ""
    for g in GENE_ORDER:
        ph = gp.get(g,"Unknown")
        pc = PHENOTYPE_COLORS.get(ph, PHENOTYPE_COLORS["Unknown"])
        bar = min(100, pc["bar"])
        boxes += f"""
        <div class="gene-box" style="background:{pc['bg']};border-color:{pc['border']};">
          <div class="gene-name" style="color:{pc['text']};">{g}</div>
          <div class="gene-bar-track">
            <div class="gene-bar-fill" style="width:{bar}%;background:{pc['bar_color']};"></div>
          </div>
          <div class="gene-pheno" style="color:{pc['text']};">{ph}</div>
        </div>"""
    sec("Gene Activity Overview")
    st.markdown(f'<div class="gene-row">{boxes}</div>', unsafe_allow_html=True)


# ══════════════════════════════════════════════════════════════════════════════
# COMPONENT: DRUG COMPARISON TABLE
# ══════════════════════════════════════════════════════════════════════════════

def render_drug_table(all_outputs, pid):
    rows_html = ""; table_data = []
    for o in all_outputs:
        drug=o["drug"]; rl=o["risk_assessment"]["risk_label"]; sev=o["risk_assessment"]["severity"]
        conf=o["risk_assessment"]["confidence_score"]; gene=o["pharmacogenomic_profile"]["primary_gene"]
        ph=o["pharmacogenomic_profile"]["phenotype"]
        rc=RISK_CONFIG.get(rl,RISK_CONFIG["Unknown"]); sp=SEV_PALETTE.get(sev,SEV_PALETTE["none"])
        rows_html += f"""<div class="dtable-row">
          <div class="dtable-cell" style="font-weight:700;color:var(--text-primary);">{drug.title()}</div>
          <div class="dtable-cell">
            <span style="display:inline-flex;align-items:center;gap:7px;background:{rc['tag_bg']};
              color:{rc['tag_text']};border:1px solid {rc['border']};border-radius:100px;
              padding:3px 10px;font-size:.75rem;font-weight:600;">
              {rc['emoji']} {rl}
            </span>
          </div>
          <div class="dtable-cell"><span style="color:{sp['text']};font-weight:600;">{sev.title()}</span></div>
          <div class="dtable-cell" style="font-family:'JetBrains Mono',monospace;font-size:.75rem;color:var(--text-muted);">{gene}</div>
          <div class="dtable-cell"><span style="font-family:'JetBrains Mono',monospace;font-size:.75rem;
            color:{rc['text']};background:{rc['bg']};border:1px solid {rc['border']};
            padding:2px 8px;border-radius:4px;font-weight:600;">{ph}</span></div>
          <div class="dtable-cell">
            <div style="flex:1;height:5px;background:var(--surface-2);border-radius:3px;overflow:hidden;
              margin-right:8px;border:1px solid var(--border);">
              <div style="width:{conf*100:.0f}%;height:100%;background:{rc['dot']};border-radius:3px;"></div>
            </div>
            <span style="font-family:'JetBrains Mono',monospace;font-size:.65rem;color:var(--text-muted);
              font-weight:500;">{conf:.0%}</span>
          </div>
        </div>"""
        table_data.append({"Drug":drug,"Risk":rl,"Severity":sev,"Gene":gene,"Phenotype":ph,"Confidence":f"{conf:.0%}"})
    sec("Drug Risk Comparison")
    st.markdown(f"""
    <div class="dtable">
      <div class="dtable-head">
        <div class="dtable-hcell">Drug</div>
        <div class="dtable-hcell">Risk Label</div>
        <div class="dtable-hcell">Severity</div>
        <div class="dtable-hcell">Gene</div>
        <div class="dtable-hcell">Phenotype</div>
        <div class="dtable-hcell">Confidence</div>
      </div>
      {rows_html}
    </div>""", unsafe_allow_html=True)
    df = pd.DataFrame(table_data)
    st.download_button("⬇ Download Table as CSV", data=df.to_csv(index=False),
                       file_name=f"pharmaguard_{pid}_comparison.csv", mime="text/csv", key=f"csv_{pid}")


# ══════════════════════════════════════════════════════════════════════════════
# COMPONENT: DRUG x GENE HEATMAP
# ══════════════════════════════════════════════════════════════════════════════

def render_heatmap(all_outputs):
    DRUG_ORDER=["CODEINE","WARFARIN","CLOPIDOGREL","SIMVASTATIN","AZATHIOPRINE","FLUOROURACIL"]
    GENE_ORDER=["CYP2D6","CYP2C9","CYP2C19","SLCO1B1","TPMT","DPYD"]
    DG={"CODEINE":"CYP2D6","WARFARIN":"CYP2C9","CLOPIDOGREL":"CYP2C19","SIMVASTATIN":"SLCO1B1","AZATHIOPRINE":"TPMT","FLUOROURACIL":"DPYD"}
    rmap={o["drug"]:o for o in all_outputs}; drugs=[d for d in DRUG_ORDER if d in rmap]
    if not drugs: return
    n=len(drugs)
    hdrs='<div class="hm-header"></div>'
    for d in drugs: hdrs+=f'<div class="hm-header">{d[:5]}</div>'
    rows=""
    for gene in GENE_ORDER:
        rows+=f'<div class="hm-header" style="justify-content:flex-end;padding-right:.5rem;">{gene}</div>'
        for d in drugs:
            if DG.get(d)==gene and d in rmap:
                o=rmap[d]; rl=o["risk_assessment"]["risk_label"]; ph=o["pharmacogenomic_profile"]["phenotype"]
                rc=RISK_CONFIG.get(rl,RISK_CONFIG["Unknown"])
                sh={"Adjust Dosage":"Adjust","Ineffective":"Ineffect.","Unknown":"?"}.get(rl,rl)
                rows+=f'<div class="hm-cell" style="background:{rc["bg"]};border-color:{rc["border"]};" title="{d}×{gene}: {rl} ({ph})"><div class="hm-dname" style="color:{rc["text"]};">{sh}</div><div class="hm-drisk" style="color:{rc["text"]};">{ph}</div></div>'
            else:
                rows+=f'<div class="hm-cell" style="background:var(--surface-2);border-color:var(--border);"><div class="hm-drisk" style="color:var(--text-xmuted);">—</div></div>'
    legend="".join(f'<div class="hm-legend-item"><span class="hm-dot" style="background:{RISK_CONFIG[r]["bg"]};border-color:{RISK_CONFIG[r]["border"]};"></span>{r}</div>'
                   for r in ["Safe","Adjust Dosage","Toxic","Ineffective"])
    st.markdown(f"""
    <div class="hm-wrap">
      <div class="hm-eye">Drug × Gene Risk Matrix</div>
      <div class="hm-grid" style="grid-template-columns:78px repeat({n},1fr);">{hdrs}{rows}</div>
      <div class="hm-legend">{legend}</div>
    </div>""", unsafe_allow_html=True)


# ══════════════════════════════════════════════════════════════════════════════
# COMPONENT: CHROMOSOME VISUALIZATION
# ══════════════════════════════════════════════════════════════════════════════

def render_chromosome(all_outputs, parsed):
    detected=set(parsed.get("detected_genes",[]))
    rmap={o["pharmacogenomic_profile"]["primary_gene"]:o for o in all_outputs}
    rows=""
    for gene,info in CHROM_INFO.items():
        ch=info["chrom"]; pos=info["pos_mb"]
        pct=(pos/CHROM_LEN.get(ch,200))*100
        if gene in rmap:
            rl=rmap[gene]["risk_assessment"]["risk_label"]
            mc=RISK_CONFIG.get(rl,RISK_CONFIG["Unknown"])["dot"]
            pulse="animation:dot-pulse 2s infinite;"
        elif gene in detected:
            mc="#94A3B8"; pulse=""
        else:
            mc="#E2E8F0"; pulse=""
        rows+=f"""<div class="chrom-row">
          <div class="chrom-lbl">{ch}</div>
          <div class="chrom-bar">
            <div class="chrom-body"></div>
            <div class="chrom-marker" style="left:{pct}%;background:{mc};color:{mc};{pulse}"></div>
          </div>
          <div class="chrom-gene-lbl">{gene}</div>
          <div class="chrom-band">{info['band']}</div>
        </div>"""
    st.markdown(f"""
    <div class="chrom-wrap">
      <div class="chrom-eye">Variant Chromosome Locations</div>
      {rows}
      <div style="font-family:'JetBrains Mono',monospace;font-size:.55rem;color:var(--text-muted);margin-top:.6rem;">Coloured markers = variants detected in this VCF</div>
    </div>""", unsafe_allow_html=True)


# ══════════════════════════════════════════════════════════════════════════════
# COMPONENT: POPULATION FREQUENCY
# ══════════════════════════════════════════════════════════════════════════════

def render_pop_freq(gene, phenotype):
    freq=POPULATION_FREQ.get(gene,{})
    if not freq: return
    rows=""
    for ph, pct in sorted(freq.items(), key=lambda x:-x[1]):
        is_you=(ph==phenotype)
        pc=PHENOTYPE_COLORS.get(ph,PHENOTYPE_COLORS["Unknown"])
        you='<span class="pop-you">← You</span>' if is_you else ""
        weight="font-weight:700;" if is_you else ""
        rows+=f"""<div class="pop-row">
          <div class="pop-ph" style="{weight}{'color:'+pc['text'] if is_you else ''}">{ph}</div>
          <div class="pop-track">
            <div class="pop-fill" style="width:{min(pct,100)}%;background:{pc['bar_color'] if is_you else '#CBD5E1'};"></div>
          </div>
          <div class="pop-pct" style="{weight}{'color:'+pc['text'] if is_you else ''}">{pct}%{you}</div>
        </div>"""
    st.markdown(f"""
    <div class="pop-wrap">
      <div class="pop-eye">{gene} — Population Frequency</div>
      {rows}
    </div>""", unsafe_allow_html=True)


# ══════════════════════════════════════════════════════════════════════════════
# COMPONENT: INTERACTION MATRIX
# ══════════════════════════════════════════════════════════════════════════════

def render_ix_matrix(all_outputs, ix_report):
    if not ix_report or len(all_outputs)<2: return
    drugs=[o["drug"] for o in all_outputs]; n=len(drugs)
    sev_map={}
    for ix in ix_report.get("all_interactions",[]):
        involved=ix.get("drugs_involved",[])
        if len(involved)==2:
            sev=ix.get("severity","none")
            k1,k2=(involved[0],involved[1]),(involved[1],involved[0])
            sev_map[k1]=sev_map[k2]=sev
    MC = {
        "critical": {"bg":"#FEF2F2","text":"#991B1B","border":"#FECACA"},
        "high":     {"bg":"#FEF2F2","text":"#991B1B","border":"#FECACA"},
        "moderate": {"bg":"#FFFBEB","text":"#92400E","border":"#FDE68A"},
        "low":      {"bg":"#FEF9C3","text":"#713F12","border":"#FDE047"},
        "none":     {"bg":"#ECFDF5","text":"#065F46","border":"#A7F3D0"},
        "diag":     {"bg":"var(--surface-2)","text":"var(--text-muted)","border":"var(--border)"},
    }
    hdrs='<div class="ix-head"></div>'
    for d in drugs: hdrs+=f'<div class="ix-head">{d[:6]}</div>'
    grid=""
    for i,d1 in enumerate(drugs):
        grid+=f'<div class="ix-head" style="justify-content:flex-end;padding-right:4px;font-size:.58rem;">{d1[:6]}</div>'
        for j,d2 in enumerate(drugs):
            if i==j:
                mc=MC["diag"]
                grid+=f'<div class="ix-cell" style="background:{mc["bg"]};border-color:{mc["border"]};color:{mc["text"]};">—</div>'
            else:
                sev=sev_map.get((d1,d2),"none")
                mc=MC.get(sev,MC["none"])
                lbl=sev.upper() if sev!="none" else "OK"
                grid+=f'<div class="ix-cell" style="background:{mc["bg"]};border-color:{mc["border"]};color:{mc["text"]};">{lbl}</div>'
    sec("Drug Interaction Matrix")
    st.markdown(f"""
    <div style="background:var(--surface);border:1px solid var(--border);border-radius:var(--r-lg);
      padding:1.25rem;margin-bottom:1rem;box-shadow:var(--shadow-sm);">
      <div class="ix-matrix-grid" style="grid-template-columns:72px repeat({n},1fr);gap:4px;">
        {hdrs}{grid}
      </div>
    </div>""", unsafe_allow_html=True)
    shown=set()
    for ix in ix_report.get("all_interactions",[]):
        involved=ix.get("drugs_involved",[])
        if len(involved)==2:
            key=tuple(sorted(involved))
            if key not in shown:
                shown.add(key)
                sev=ix.get("severity","low"); sp=SEV_PALETTE.get(sev,SEV_PALETTE["low"])
                with st.expander(f"{' + '.join(involved)}  —  {sev.upper()} interaction"):
                    mech=ix.get("mechanism",ix.get("message",""))
                    rec=ix.get("recommendation","")
                    if mech: st.markdown(f'<div style="font-size:.9rem;color:var(--text-secondary);line-height:1.75;margin-bottom:.5rem;">{mech}</div>',unsafe_allow_html=True)
                    if rec:  st.markdown(f'<div style="font-family:JetBrains Mono,monospace;font-size:.73rem;color:{sp["text"]};margin-top:.5rem;font-weight:600;">→ Recommendation: {rec}</div>',unsafe_allow_html=True)


# ══════════════════════════════════════════════════════════════════════════════
# COMPONENT: AI NARRATIVE
# ══════════════════════════════════════════════════════════════════════════════

def render_ai_narrative(all_outputs, parsed, pid, groq_key, skip_llm):
    results_for_nar=[{"drug":o["drug"],"primary_gene":o["pharmacogenomic_profile"]["primary_gene"],
        "phenotype":o["pharmacogenomic_profile"]["phenotype"],"risk_label":o["risk_assessment"]["risk_label"],
        "severity":o["risk_assessment"]["severity"]} for o in all_outputs]
    with st.spinner("Generating AI clinical summary…"):
        narrative=generate_patient_narrative(pid,results_for_nar,parsed,groq_key,skip_llm)
    badge="Static Template" if (skip_llm or not groq_key) else "LLaMA 3.3 70B"
    sec("AI Clinical Summary")
    st.markdown(f"""
    <div class="ai-narrative">
      <div class="ai-nar-head">
        <span class="ai-nar-badge">{badge}</span>
        <span class="ai-nar-title">Unified Patient Summary</span>
      </div>
      <div class="ai-nar-text">{narrative}</div>
    </div>""", unsafe_allow_html=True)


# ══════════════════════════════════════════════════════════════════════════════
# COMPONENT: BEFORE / AFTER
# ══════════════════════════════════════════════════════════════════════════════

def render_before_after(all_outputs):
    dangerous=[o for o in all_outputs if o["risk_assessment"]["risk_label"] in ("Toxic","Ineffective")]
    if not dangerous: return
    o=dangerous[0]; drug=o["drug"]; rl=o["risk_assessment"]["risk_label"]
    alts=o["clinical_recommendation"].get("alternative_drugs",[])
    alt=alts[0] if alts else "Alternative medication"
    gene=o["pharmacogenomic_profile"]["primary_gene"]; ph=o["pharmacogenomic_profile"]["phenotype"]
    BEFORE={"Toxic":f"Standard {drug.lower()} dose → toxic accumulation → life-threatening outcome",
            "Ineffective":f"Standard {drug.lower()} dose → zero therapeutic effect → treatment failure"}
    sec("Before / After — PGx Impact")
    st.markdown(f"""
    <div class="ba-wrap">
      <div class="ba-side" style="background:#FFF1F1;">
        <div class="ba-label" style="color:#991B1B;">⛔ Without PharmaGuard</div>
        <div class="ba-drug" style="color:#7F1D1D;">{drug.title()} — Standard Protocol</div>
        <div class="ba-outcome" style="color:#991B1B;">{BEFORE.get(rl,"Risk undetected")}</div>
        <div style="margin-top:.75rem;font-family:'JetBrains Mono',monospace;font-size:.62rem;color:#FCA5A5;">
          {gene} {ph} phenotype undetected
        </div>
      </div>
      <div class="ba-side" style="background:#F0FDF4;">
        <div class="ba-label" style="color:#065F46;">✅ With PharmaGuard</div>
        <div class="ba-drug" style="color:#065F46;">{alt} — PGx-Guided Protocol</div>
        <div class="ba-outcome" style="color:#059669;">PGx-guided alternative → appropriate dosing → safe, effective therapy</div>
        <div style="margin-top:.75rem;font-family:'JetBrains Mono',monospace;font-size:.62rem;color:#A7F3D0;">
          {gene} {ph} phenotype identified → therapy optimised
        </div>
      </div>
    </div>""", unsafe_allow_html=True)


# ══════════════════════════════════════════════════════════════════════════════
# COMPONENT: SIDE-BY-SIDE DRUG COMPARISON
# ══════════════════════════════════════════════════════════════════════════════

def render_drug_comparison(all_outputs):
    dangerous=[o for o in all_outputs if o["risk_assessment"]["risk_label"] in ("Toxic","Ineffective")]
    if not dangerous: return
    for o in dangerous[:2]:
        drug=o["drug"]; gene=o["pharmacogenomic_profile"]["primary_gene"]; ph=o["pharmacogenomic_profile"]["phenotype"]
        alts=o["clinical_recommendation"].get("alternative_drugs",[]); rc=RISK_CONFIG.get(o["risk_assessment"]["risk_label"],RISK_CONFIG["Unknown"])
        if not alts: continue
        alt=alts[0]
        st.markdown(f"""
        <div style="margin-bottom:1.25rem;">
          <div style="display:grid;grid-template-columns:1fr 1fr;border:1.5px solid var(--border);
            border-radius:var(--r-xl);overflow:hidden;box-shadow:var(--shadow-md);">
            <div style="background:{rc['bg']};padding:1.25rem 1.5rem;border-right:1px solid var(--border);">
              <div style="font-family:'JetBrains Mono',monospace;font-size:.58rem;letter-spacing:.1em;
                text-transform:uppercase;color:{rc['text']};margin-bottom:.5rem;font-weight:600;">Current Prescription</div>
              <div style="font-size:1.05rem;font-weight:700;color:{rc['text']};margin-bottom:.3rem;">{drug.title()}</div>
              <div style="font-family:'JetBrains Mono',monospace;font-size:.65rem;color:{rc['text']};
                margin-bottom:.5rem;">{rc['emoji']} {o["risk_assessment"]["risk_label"]} · {gene} {ph}</div>
              <div style="font-size:.85rem;color:{rc['text']};opacity:.85;line-height:1.65;">
                {o["clinical_recommendation"]["dosing_recommendation"][:150]}…</div>
            </div>
            <div style="background:#F0FDF4;padding:1.25rem 1.5rem;">
              <div style="font-family:'JetBrains Mono',monospace;font-size:.58rem;letter-spacing:.1em;
                text-transform:uppercase;color:#065F46;margin-bottom:.5rem;font-weight:600;">Recommended Alternative</div>
              <div style="font-size:1.05rem;font-weight:700;color:#065F46;margin-bottom:.3rem;">{alt}</div>
              <div style="font-family:'JetBrains Mono',monospace;font-size:.65rem;color:#059669;margin-bottom:.5rem;">
                ✅ PGx-guided selection</div>
              <div style="font-size:.85rem;color:#059669;opacity:.9;line-height:1.65;">
                Selected based on {gene} {ph} phenotype per CPIC Level A guidelines.</div>
            </div>
          </div>
        </div>""", unsafe_allow_html=True)


# ══════════════════════════════════════════════════════════════════════════════
# COMPONENT: PRESCRIPTION SAFETY CHECKER
# ══════════════════════════════════════════════════════════════════════════════

def render_rx_checker(all_outputs):
    sec("Prescription Safety Checker")
    rmap={o["drug"]:o for o in all_outputs}
    drug_names=[o["drug"] for o in all_outputs]
    c1,c2=st.columns([2,1])
    with c1:
        selected=st.selectbox("Select drug",drug_names,
            format_func=lambda x:f"{x.title()}  ({GENE_DRUG_MAP.get(x,'')})",
            key="rx_drug",label_visibility="collapsed")
    with c2:
        check=st.button("Check Safety →",key="rx_check")
    if check and selected in rmap:
        o=rmap[selected]; rl=o["risk_assessment"]["risk_label"]; sev=o["risk_assessment"]["severity"]
        rec=o["clinical_recommendation"]["dosing_recommendation"]
        gene=o["pharmacogenomic_profile"]["primary_gene"]; ph=o["pharmacogenomic_profile"]["phenotype"]
        rc=RISK_CONFIG.get(rl,RISK_CONFIG["Unknown"]); sp=SEV_PALETTE.get(sev,SEV_PALETTE["none"])
        VERDICT={"Safe":"✅ Safe to Prescribe","Adjust Dosage":"⚠️ Prescribe with Dose Adjustment",
                 "Toxic":"🚨 Do Not Prescribe — High Toxicity Risk",
                 "Ineffective":"❌ Do Not Prescribe — Drug Will Be Ineffective"}
        verdict=VERDICT.get(rl,"❓ Insufficient Data")
        st.markdown(f"""
        <div class="rx-result" style="background:{rc['bg']};border-color:{rc['border']};">
          <div class="rx-verdict" style="color:{rc['text']};">{verdict}</div>
          <div class="rx-detail">{gene} {ph} phenotype detected. {rec}</div>
          <div class="rx-meta" style="color:{sp['text']};">
            Severity: {sev} · Confidence: {o["risk_assessment"]["confidence_score"]:.0%} · CPIC Level A
          </div>
        </div>""", unsafe_allow_html=True)
    elif not check:
        st.markdown("""
        <div class="info-strip">
          <span>🔍</span>
          <div class="info-strip-text">Select a drug above and click <strong>Check Safety</strong> to
          instantly validate this prescription against the patient's genotype.</div>
        </div>""", unsafe_allow_html=True)


# ══════════════════════════════════════════════════════════════════════════════
# COMPONENT: CLINICAL NOTE
# ══════════════════════════════════════════════════════════════════════════════

def render_clinical_note(all_outputs, pid):
    lines=[f"Patient {pid} — Pharmacogenomic Analysis — {datetime.utcnow().strftime('%Y-%m-%d')}",""]
    for o in all_outputs:
        gene=o["pharmacogenomic_profile"]["primary_gene"]; dip=o["pharmacogenomic_profile"]["diplotype"]
        ph=o["pharmacogenomic_profile"]["phenotype"]; drug=o["drug"]
        rl=o["risk_assessment"]["risk_label"]; rec=o["clinical_recommendation"]["dosing_recommendation"]
        alts=o["clinical_recommendation"].get("alternative_drugs",[])
        lines.append(f"Patient carries {gene} {dip} ({ph} phenotype), therefore {drug.lower()} is predicted {rl.lower()}.")
        lines.append(f"CPIC Recommendation: {rec}")
        if alts: lines.append(f"Alternatives: {', '.join(alts)}")
        lines.append("")
    lines.append("Generated by PharmaGuard v8.0 · CPIC Level A evidence · cpicpgx.org")
    note="\n".join(lines)
    sec("One-Click Clinical Note")
    st.markdown(f'<div class="note-box"><pre>{note}</pre></div>',unsafe_allow_html=True)
    st.download_button("⬇ Download Clinical Note",data=note,
        file_name=f"clinical_note_{pid}.txt",mime="text/plain",key=f"note_{pid}")


# ══════════════════════════════════════════════════════════════════════════════
# COMPONENT: PATIENT PLAIN-ENGLISH MODE
# ══════════════════════════════════════════════════════════════════════════════

def render_patient_mode(all_outputs):
    bad=any(o["risk_assessment"]["risk_label"] in ("Toxic","Ineffective") for o in all_outputs)
    if bad:
        st.markdown("""<div class="patient-banner" style="background:#FFF1F1;border-color:#FECACA;">
          <div class="patient-banner-title" style="color:#991B1B;">🚨 Important — Some medications need urgent attention</div>
          <div class="patient-banner-sub" style="color:#B91C1C;">Your genetic results show that one or more
          medications may not be safe or effective for you. Please speak with your doctor before taking these medications.</div>
        </div>""",unsafe_allow_html=True)
    else:
        st.markdown("""<div class="patient-banner" style="background:#F0FDF4;border-color:#A7F3D0;">
          <div class="patient-banner-title" style="color:#065F46;">✅ Good news — Your medications look safe</div>
          <div class="patient-banner-sub" style="color:#059669;">Based on your genetic profile, the medications
          reviewed are predicted to work normally in your body at standard doses.</div>
        </div>""",unsafe_allow_html=True)

    for o in all_outputs:
        drug=o["drug"]; rl=o["risk_assessment"]["risk_label"]
        gene=o["pharmacogenomic_profile"]["primary_gene"]; ph=o["pharmacogenomic_profile"]["phenotype"]
        alts=o["clinical_recommendation"].get("alternative_drugs",[])
        phplain=PLAIN_ENGLISH_PHENOTYPE.get(ph,ph); explain=PLAIN_ENGLISH_RISK.get((drug,ph),"")
        VERDICT={"Safe":"✅ This medicine is likely safe for you","Adjust Dosage":"⚠️ You may need a different dose",
                 "Toxic":"🚨 This medicine could be harmful to you","Ineffective":"❌ This medicine likely won't work for you"}
        verdict=VERDICT.get(rl,rl)
        bc={"Safe":"#A7F3D0","Adjust Dosage":"#FDE68A","Toxic":"#FECACA","Ineffective":"#DDD6FE"}.get(rl,"#E2E8F0")
        tc={"Safe":"#065F46","Adjust Dosage":"#92400E","Toxic":"#991B1B","Ineffective":"#5B21B6"}.get(rl,"#475569")
        action=""
        if rl in ("Toxic","Ineffective"):
            alt_text=f"They may suggest: <strong>{', '.join(alts[:3])}</strong>" if alts else "Ask about alternative medications."
            action=f'<div class="patient-action"><span style="font-size:1.1rem;">💊</span><div class="patient-action-text"><strong>Talk to your doctor before taking {drug.title()}.</strong><br>{alt_text}</div></div>'
        elif rl=="Adjust Dosage":
            action=f'<div class="patient-action"><span style="font-size:1.1rem;">📋</span><div class="patient-action-text"><strong>Tell your doctor about this result before starting {drug.title()}.</strong><br>You may need a different dose than usually prescribed.</div></div>'
        st.markdown(f"""
        <div class="patient-card" style="border-color:{bc};">
          <div class="patient-drug">{drug.title()}</div>
          <div class="patient-verdict" style="color:{tc};">{verdict}</div>
          <div class="patient-gene-tag">{gene} · {phplain}</div>
          {f'<div class="patient-plain">{explain}</div>' if explain else ''}
          {action}
        </div>""",unsafe_allow_html=True)


# ══════════════════════════════════════════════════════════════════════════════
# MASTER RESULTS RENDERER
# ══════════════════════════════════════════════════════════════════════════════

def render_results(all_outputs, parsed, ix_report, pdf_bytes, pid,
                   patient_mode=False, groq_key="", skip_llm=False):

    # 1. Risk Banner
    render_risk_banner(all_outputs, parsed)

    # 2. Emergency Alerts
    render_emergency_alerts(all_outputs)

    # 3. Gene Activity Heatmap
    render_gene_heatmap(all_outputs)

    # 4. Download row
    dc1,dc2,dc3=st.columns(3)
    with dc1:
        st.download_button("⬇ Download All JSON",data=json.dumps(all_outputs,indent=2),
            file_name=f"pharmaguard_{pid}.json",mime="application/json",
            use_container_width=True,key=f"dlall_{pid}")
    with dc2:
        if pdf_bytes:
            st.download_button("⬇ Download PDF Report",data=pdf_bytes,
                file_name=f"pharmaguard_{pid}.pdf",mime="application/pdf",
                use_container_width=True,key=f"dlpdf_{pid}")
    with dc3:
        if ix_report and ix_report.get("interactions_found"):
            st.download_button("⬇ Interactions JSON",data=json.dumps(ix_report,indent=2),
                file_name=f"pharmaguard_{pid}_ix.json",mime="application/json",
                use_container_width=True,key=f"dlix_{pid}")

    st.markdown("<div style='height:.75rem'></div>",unsafe_allow_html=True)

    if patient_mode:
        render_patient_mode(all_outputs)
        return

    # ── DOCTOR MODE ──────────────────────────────────────────

    # 5. Drug table
    render_drug_table(all_outputs, pid)

    # 6. Polygenic score
    render_pgx_score(all_outputs)

    # 7. Heatmap + Chromosome
    c1,c2=st.columns([1.4,1],gap="large")
    with c1: render_heatmap(all_outputs)
    with c2: render_chromosome(all_outputs,parsed)

    # 8. Interaction matrix
    if ix_report and len(all_outputs)>=2:
        render_ix_matrix(all_outputs,ix_report)

    # 9. AI Narrative
    render_ai_narrative(all_outputs,parsed,pid,groq_key,skip_llm)

    # 10. Before/After + Drug Comparison
    render_before_after(all_outputs)
    render_drug_comparison(all_outputs)

    # 11. Rx Checker
    render_rx_checker(all_outputs)

    # 12. Clinical Note
    render_clinical_note(all_outputs,pid)

    # 13. Individual Drug Cards
    sec("Individual Drug Analysis")

    for output in all_outputs:
        rl=output["risk_assessment"]["risk_label"]; dn=output["drug"]
        sev=output["risk_assessment"]["severity"]; conf=output["risk_assessment"]["confidence_score"]
        gene=output["pharmacogenomic_profile"]["primary_gene"]; dip=output["pharmacogenomic_profile"]["diplotype"]
        ph=output["pharmacogenomic_profile"]["phenotype"]; var=output["pharmacogenomic_profile"]["detected_variants"]
        rec=output["clinical_recommendation"]["dosing_recommendation"]
        alts=output["clinical_recommendation"].get("alternative_drugs",[])
        mon=output["clinical_recommendation"].get("monitoring_required","")
        exp=output["llm_generated_explanation"]
        rc=RISK_CONFIG.get(rl,RISK_CONFIG["Unknown"]); sp=SEV_PALETTE.get(sev,SEV_PALETTE["none"])
        cpic_lv=output.get("pharmacogenomic_profile",{}).get("cpic_evidence_level","Level A")

        st.markdown(f"""
        <div class="rcard">
          <div class="rcard-top">
            <div class="rcard-left">
              <div class="rcard-dot" style="background:{rc['dot']};box-shadow:0 0 0 3px {rc['bg']};"></div>
              <div>
                <div class="rcard-name">{dn.title()} <span class="cpic">CPIC {cpic_lv}</span></div>
                <div class="rcard-meta">{gene} · {dip} · {ph}</div>
              </div>
            </div>
            <span class="rcard-badge" style="background:{rc['tag_bg']};color:{rc['tag_text']};border-color:{rc['border']};">
              {rc['emoji']} {rl}
            </span>
          </div>
          <div class="rcard-body">
            <div class="mc-row">
              <div class="mc-cell">
                <div class="mc-key">Phenotype</div>
                <div class="mc-val" style="color:{rc['text']};">{ph}</div>
              </div>
              <div class="mc-cell">
                <div class="mc-key">Severity</div>
                <div class="mc-val" style="color:{sp['text']};">{sev.title()}</div>
              </div>
              <div class="mc-cell">
                <div class="mc-key">Confidence</div>
                <div class="mc-val">{conf:.0%}</div>
              </div>
              <div class="mc-cell">
                <div class="mc-key">Variants</div>
                <div class="mc-val">{len(var)}</div>
              </div>
            </div>""", unsafe_allow_html=True)

        # Confidence bars
        dq = min(1.0, len(var)/3.0)
        st.markdown(f"""
        <div class="conf-row">
          <div>
            <div class="conf-lbl">
              <span>Prediction Confidence</span>
              <span style="color:{rc['dot']};font-weight:700;">{conf:.0%}</span>
            </div>
            <div class="conf-track">
              <div class="conf-fill" style="width:{conf*100:.1f}%;background:{rc['dot']};"></div>
            </div>
          </div>
          <div>
            <div class="conf-lbl">
              <span>Data Quality</span>
              <span style="color:var(--text-muted);">{len(var)} variant{"s" if len(var)!=1 else ""}</span>
            </div>
            <div class="conf-track">
              <div class="conf-fill" style="width:{dq*100:.1f}%;background:#94A3B8;"></div>
            </div>
          </div>
        </div>""", unsafe_allow_html=True)

        # Variants
        if var:
            rows_html=""
            for v in var:
                fc=func_cls(v.get("functional_status",""))
                fn=(v.get("functional_status") or "unknown").replace("_"," ").title()
                rows_html+=f'<tr><td class="v-rsid">{v.get("rsid","—")}</td><td class="v-star">{v.get("star_allele","—")}</td><td class="{fc}">{fn}</td></tr>'
            st.markdown(f"""<hr class="h-rule">
            <div class="sec-label">Detected Variants ({len(var)})</div>
            <table class="vtable">
              <thead><tr><th>rsID</th><th>Star Allele</th><th>Functional Status</th></tr></thead>
              <tbody>{rows_html}</tbody>
            </table>""", unsafe_allow_html=True)

        # CPIC Rec
        st.markdown(f"""<hr class="h-rule">
        <div class="sec-label">CPIC Recommendation</div>
        <div class="rec-box" style="background:{rc['bg']};border-color:{rc['border']};">
          <div class="rec-lbl" style="color:{rc['text']};">CPIC Guideline — {dn}</div>
          <div class="rec-txt">{rec}</div>
        </div>""", unsafe_allow_html=True)

        if alts:
            chips="".join(f'<span class="alt-chip">{a}</span>' for a in alts)
            st.markdown(f'<div class="sec-label">Alternative Drugs</div><div class="alt-chips" style="margin-bottom:1rem;">{chips}</div>',unsafe_allow_html=True)

        if mon:
            st.markdown(f"""<div style="display:flex;gap:.75rem;align-items:flex-start;padding:.875rem 1rem;
                background:var(--surface-2);border:1px solid var(--border);border-radius:var(--r-md);margin-bottom:1rem;">
              <span style="font-size:1rem;margin-top:1px;">🔬</span>
              <div>
                <div class="sec-label" style="margin-bottom:.2rem;">Monitoring Protocol</div>
                <div style="font-size:.875rem;color:var(--text-secondary);line-height:1.65;">{mon}</div>
              </div>
            </div>""", unsafe_allow_html=True)

        st.markdown("</div></div>",unsafe_allow_html=True)

        render_pop_freq(gene, ph)

        if exp.get("summary"):
            is_static="static" in exp.get("model_used","").lower()
            model=exp.get("model_used","llama-3.3-70b")
            blocks=""
            for lbl,key in [("Summary","summary"),("Biological Mechanism","biological_mechanism"),
                             ("Variant Significance","variant_significance"),("Clinical Implications","clinical_implications")]:
                if exp.get(key):
                    blocks+=f'<div class="ai-section"><div class="ai-section-lbl">{lbl}</div><div class="ai-section-txt">{exp[key]}</div></div>'
            st.markdown(f"""<div class="ai-block">
              <div class="ai-block-head">
                <span class="ai-badge">{model}</span>
                <span style="font-size:.8rem;color:var(--text-muted);font-weight:500;">AI Explanation · {dn}</span>
              </div>{blocks}</div>""",unsafe_allow_html=True)

        with st.expander(f"Raw JSON — {dn}", expanded=False):
            c1,_=st.columns([1,3])
            with c1:
                st.download_button(f"⬇ {dn} JSON",data=json.dumps(output,indent=2),
                    file_name=f"pharmaguard_{pid}_{dn}.json",mime="application/json",
                    key=f"dl_{pid}_{dn}",use_container_width=True)
            st.code(json.dumps(output,indent=2),language="json")

    with st.expander("VCF Parse Details", expanded=False):
        p1,p2,p3=st.columns(3)
        p1.metric("Total Variants",parsed["total_variants"])
        p2.metric("Genes Found",len(parsed["detected_genes"]))
        p3.metric("Parse Errors",len(parsed["parse_errors"]))


# ══════════════════════════════════════════════════════════════════════════════
# NAV
# ══════════════════════════════════════════════════════════════════════════════

st.markdown("""
<div class="pg-nav">
  <div class="pg-logo">
    <div class="pg-logo-text">Pharma<em>Guard</em></div>
    <div class="pg-logo-sub">v8 · Precision Clinical</div>
  </div>
  <div class="pg-nav-right">
    <span class="pg-chip pg-chip-default">CPIC Aligned</span>
    <span class="pg-chip pg-chip-default">RIFT 2026</span>
    <span class="pg-chip pg-chip-brand">★ Precision Clinical</span>
  </div>
</div>""", unsafe_allow_html=True)

# ── Sidebar ───────────────────────────────────────────────────────────────────
with st.sidebar:
    st.markdown("""<div style="padding:1rem 0 .5rem;font-family:'JetBrains Mono',monospace;font-size:.6rem;
        letter-spacing:.12em;text-transform:uppercase;color:var(--text-muted);">Groq API Key</div>""",unsafe_allow_html=True)
    groq_api_key=st.text_input("Groq API Key",value=os.environ.get("GROQ_API_KEY",""),
        type="password",label_visibility="collapsed",placeholder="gsk_…")
    st.markdown("""<div style="font-family:'JetBrains Mono',monospace;font-size:.6rem;color:var(--text-muted);
        padding-bottom:1rem;line-height:1.9;">Model: LLaMA 3.3 70B Versatile<br>
        Fallback: static expert templates<br>Test mode: instant (no API call)</div>""",unsafe_allow_html=True)
    st.divider()
    st.markdown("""<div style="padding:.5rem 0;font-family:'JetBrains Mono',monospace;font-size:.6rem;
        letter-spacing:.12em;text-transform:uppercase;color:var(--text-muted);">Gene → Drug Map</div>""",unsafe_allow_html=True)
    for g,drug in [("CYP2D6","Codeine"),("CYP2C19","Clopidogrel"),("CYP2C9","Warfarin"),
                   ("SLCO1B1","Simvastatin"),("TPMT","Azathioprine"),("DPYD","Fluorouracil")]:
        st.markdown(f"""<div style="font-family:'JetBrains Mono',monospace;font-size:.68rem;
            padding:4px 0;color:var(--text-secondary);">{g}
            <span style="color:var(--text-muted);">→</span> {drug}</div>""",unsafe_allow_html=True)

tab1,tab2=st.tabs(["Analysis","Test Suite"])


# ══════════════════════════════════════════════════════════════════════════════
# TAB 1 — ANALYSIS
# ══════════════════════════════════════════════════════════════════════════════

with tab1:
    # Persona buttons
    sec("Quick Demo — Select Patient Persona")
    p_cols=st.columns(4)
    for i,(key,p) in enumerate(PERSONAS.items()):
        with p_cols[i]:
            st.markdown(f"""
            <div class="persona-card" style="background:{p['bg']};border-color:{p['border']};">
              <div class="persona-label" style="color:{p['text']};">{p['label']}</div>
              <div class="persona-desc" style="color:{p['accent']};">{p['desc']}</div>
            </div>""",unsafe_allow_html=True)
            if st.button("Load →",key=f"persona_{key}",use_container_width=True):
                st.session_state["persona_file"]=p["file"]
                st.session_state["persona_drugs"]=p["drugs"]

    st.markdown("<div style='height:1.25rem'></div>",unsafe_allow_html=True)

    # Steps
    st.markdown("""<div class="steps-bar">
      <div class="step active"><div class="step-n">01</div><div class="step-l">Upload VCF</div></div>
      <div class="step active"><div class="step-n">02</div><div class="step-l">Select Drugs</div></div>
      <div class="step"><div class="step-n">03</div><div class="step-l">Run Analysis</div></div>
      <div class="step"><div class="step-n">04</div><div class="step-l">Review Results</div></div>
    </div>""",unsafe_allow_html=True)

    col_l,col_r=st.columns([1.3,1],gap="large")

    with col_l:
        sec("Genomic Data")
        uploaded=st.file_uploader("Upload VCF",type=["vcf"])
        if uploaded:
            sz=uploaded.size/(1024*1024)
            if sz>5:
                st.error(f"File too large ({sz:.1f} MB). Max 5 MB."); uploaded=None
            else:
                peek=uploaded.read(400).decode("utf-8",errors="replace"); uploaded.seek(0)
                if "##fileformat=VCF" not in peek and "#CHROM" not in peek:
                    st.error("Invalid VCF file."); uploaded=None
                else:
                    st.success(f"✓  {uploaded.name}  ·  {sz:.2f} MB")

        sec("Or select a test scenario")
        scenario_opts={
            "None":None,
            "Mixed Variants (Standard)":"sample.vcf",
            "UltraRapid Metabolizer":"test_ultrarapid_metabolizer.vcf",
            "All Normal Wild-type":"test_all_normal_wildtype.vcf",
            "Worst Case — All Poor Metabolizers":"test_worst_case_all_pm.vcf",
            "Patient A — Critical Risk":"patient_a_critical.vcf",
            "Patient B — Warfarin PM":"patient_b_warfarin.vcf",
            "Patient C — Clopidogrel PM":"patient_c_interaction.vcf",
            "Patient D — All Safe":"patient_d_safe.vcf",
        }
        chosen_label=st.selectbox("Scenario",list(scenario_opts.keys()),label_visibility="collapsed")
        chosen_file=scenario_opts[chosen_label]

    with col_r:
        sec("Medications")
        default_drugs=st.session_state.get("persona_drugs",["CLOPIDOGREL"])
        drugs_selected=st.multiselect("Select drugs",options=ALL_DRUGS,default=default_drugs,
            format_func=lambda x:f"{x.title()}  ({GENE_DRUG_MAP.get(x,'')})",
            label_visibility="collapsed")

        sec("Custom drugs (comma-separated)")
        custom_drugs=st.text_input("Custom",placeholder="CODEINE, WARFARIN…",label_visibility="collapsed")

        sec("Patient ID")
        pid_input=st.text_input("Patient ID",placeholder="Auto-generated if blank",label_visibility="collapsed")

        c1,c2=st.columns(2)
        with c1: do_ix=st.checkbox("Check drug interactions",value=True)
        with c2: do_pdf=st.checkbox("Generate PDF report",value=True)

        sec("View Mode")
        patient_mode=st.toggle("Patient Plain-English Mode",
            help="Converts clinical jargon into plain language any patient can understand")

    st.markdown("<div style='height:.5rem'></div>",unsafe_allow_html=True)
    run_btn=st.button("Run Analysis →",use_container_width=True)

    if run_btn:
        all_drugs=list(drugs_selected)
        if custom_drugs.strip():
            all_drugs+=[d.strip().upper() for d in custom_drugs.split(",") if d.strip()]
        all_drugs=list(set(all_drugs))
        if not all_drugs:
            st.error("Please select at least one drug."); st.stop()

        vcf_content=None
        if "persona_file" in st.session_state and not uploaded and not chosen_file:
            vcf_content=load_vcf_file(st.session_state["persona_file"])
        elif uploaded:
            vcf_content=uploaded.read().decode("utf-8",errors="replace")
        elif chosen_file:
            vcf_content=load_vcf_file(chosen_file)
        else:
            st.error("Please upload a VCF or select a test scenario."); st.stop()

        pid=pid_input.strip() or f"PG-{str(uuid.uuid4())[:8].upper()}"

        st.markdown(f"""
        <div style="display:flex;align-items:baseline;gap:1rem;margin:2.5rem 0 1.5rem;
            padding-bottom:1rem;border-bottom:1.5px solid var(--border);">
          <div style="font-family:'Playfair Display',serif;font-size:1.9rem;font-weight:500;
            color:var(--text-primary);letter-spacing:-.02em;">Analysis Results</div>
          <div style="font-family:'JetBrains Mono',monospace;font-size:.7rem;
            color:var(--text-muted);font-weight:400;">{pid}</div>
        </div>""",unsafe_allow_html=True)

        with st.spinner("Analysing genomic data…"):
            parsed,risk_results,all_outputs,ix_report,pdf_bytes=run_pipeline(
                vcf_content,all_drugs,pid,groq_api_key,do_ix,do_pdf)

        render_results(all_outputs,parsed,ix_report,pdf_bytes,pid,
            patient_mode=patient_mode,groq_key=groq_api_key,skip_llm=(not groq_api_key))

    else:
        st.markdown("""
        <div class="empty-state">
          <span class="empty-icon">🧬</span>
          <div class="empty-title">Ready for analysis</div>
          <div class="empty-hint">
            Click a Patient Persona above for an instant demo<br>
            Or upload a VCF file · select medications · click Run Analysis<br><br>
            Genes: CYP2D6 · CYP2C19 · CYP2C9 · SLCO1B1 · TPMT · DPYD<br><br>
            v8.0 Precision Clinical · Light Theme · WCAG AA Accessible
          </div>
        </div>""",unsafe_allow_html=True)


# ══════════════════════════════════════════════════════════════════════════════
# TAB 2 — TEST SUITE
# ══════════════════════════════════════════════════════════════════════════════

RISK_DOT={"Safe":"#059669","Adjust Dosage":"#D97706","Toxic":"#DC2626","Ineffective":"#7C3AED","Unknown":"#94A3B8"}

with tab2:
    st.markdown("""
    <div style="margin-bottom:2rem;">
      <div style="font-family:'Playfair Display',serif;font-size:1.9rem;font-weight:500;
        color:var(--text-primary);letter-spacing:-.02em;margin-bottom:.4rem;">Test Suite</div>
      <div style="font-family:'JetBrains Mono',monospace;font-size:.65rem;
        color:var(--text-muted);letter-spacing:.04em;">
        4 clinical scenarios · parallel execution · pass/fail validation per drug
      </div>
    </div>""",unsafe_allow_html=True)

    tc=st.columns(4)
    for i,sc in enumerate(TEST_SUITE):
        with tc[i]:
            rh="".join(
                f'<div style="display:flex;align-items:center;gap:.5rem;margin-bottom:3px;">'
                f'<span style="width:7px;height:7px;border-radius:50%;background:{RISK_DOT.get(r,"#94A3B8")};flex-shrink:0;"></span>'
                f'<span style="font-family:JetBrains Mono,monospace;font-size:.65rem;color:var(--text-secondary);">{d[:10]}</span>'
                f'<span style="font-family:JetBrains Mono,monospace;font-size:.65rem;font-weight:600;color:{RISK_DOT.get(r,"#94A3B8);")}margin-left:auto;">{r}</span>'
                f'</div>'
                for d,r in list(sc["expected"].items())[:4]
            )
            st.markdown(f'<div class="test-card"><div class="test-name">{sc["name"]}</div>'
                        f'<div class="test-desc">{sc["desc"]}</div>{rh}</div>',unsafe_allow_html=True)

    st.markdown("<div style='height:1.25rem'></div>",unsafe_allow_html=True)
    tb1,tb2=st.columns([3,1])
    with tb1: use_llm=st.checkbox("Include LLM Explanations (requires Groq API key)",value=False)
    with tb2: run_all=st.button("Run All 4 Tests →",use_container_width=True)

    if run_all:
        def run_one_test(sc):
            vcf=load_vcf_file(sc["file"])
            pid=f"TEST-{sc['name'][:6].replace(' ','').upper()}"
            pv,_,ao,_,_=run_pipeline(vcf,sc["drugs"],pid,
                groq_api_key if use_llm else "",run_ix=False,gen_pdf=False,skip_llm=not use_llm)
            rows,ok=[],True
            for out in ao:
                drug=out["drug"]; got=out["risk_assessment"]["risk_label"]
                exp=sc["expected"].get(drug,""); passed=(got==exp) if exp else True
                ph=out["pharmacogenomic_profile"]["phenotype"]; dp=out["pharmacogenomic_profile"]["diplotype"]
                rows.append((drug,got,exp,passed,ph,dp))
                if not passed: ok=False
            return {"name":sc["name"],"pass":ok,"rows":rows,"outputs":ao,"file":sc["file"]}

        prog=st.empty(); prog.info("Running all 4 scenarios in parallel…")
        results=[None]*4
        with ThreadPoolExecutor(max_workers=4) as ex:
            futs={ex.submit(run_one_test,sc):i for i,sc in enumerate(TEST_SUITE)}
            done=0
            for f in as_completed(futs):
                results[futs[f]]=f.result(); done+=1
                prog.info(f"Completed {done}/4 scenarios…")
        prog.empty()

        passed=sum(1 for r in results if r["pass"]); failed=4-passed
        oc="#059669" if failed==0 else "#D97706"
        ob="#ECFDF5" if failed==0 else "#FFFBEB"
        od="#A7F3D0" if failed==0 else "#FDE68A"
        st.markdown(f"""
        <div style="background:{ob};border:1.5px solid {od};border-radius:var(--r-lg);
            padding:1.25rem 1.5rem;margin:1.25rem 0;display:flex;align-items:center;
            justify-content:space-between;box-shadow:var(--shadow-sm);">
          <div style="font-family:'Playfair Display',serif;font-size:1.35rem;font-weight:500;color:{oc};">
            {'All tests passed ✓' if failed==0 else f'{passed}/4 tests passed'}
          </div>
          <div style="font-family:'JetBrains Mono',monospace;font-size:.65rem;font-weight:600;color:{oc};">
            {passed} passed · {failed} failed · {int(passed/4*100)}%
          </div>
        </div>""",unsafe_allow_html=True)

        for sr in results:
            sym="✓ Pass" if sr["pass"] else "✗ Fail"
            with st.expander(f"{sym}  —  {sr['name']}",expanded=not sr["pass"]):
                for drug,got,exp,ok,ph,dp in sr["rows"]:
                    rc=RISK_CONFIG.get(got,RISK_CONFIG["Unknown"])
                    ok_bg="#ECFDF5" if ok else "#FEF2F2"
                    ok_cl="#059669" if ok else "#DC2626"
                    st.markdown(f"""
                    <div style="display:grid;grid-template-columns:1fr 1.2fr 1.2fr 1.5fr 60px;
                        border:1px solid var(--border);border-radius:var(--r-md);
                        padding:.5rem .875rem;margin-bottom:.4rem;background:var(--surface);
                        box-shadow:var(--shadow-sm);align-items:center;">
                      <div style="font-weight:700;color:var(--text-primary);font-size:.85rem;">{drug}</div>
                      <div style="display:inline-flex;align-items:center;gap:6px;">
                        <span style="width:7px;height:7px;border-radius:50%;background:{rc['dot']};flex-shrink:0;"></span>
                        <span style="font-size:.82rem;color:{rc['text']};font-weight:600;">{got}</span>
                      </div>
                      <div style="font-family:'JetBrains Mono',monospace;font-size:.75rem;color:var(--text-muted);">{exp or '—'}</div>
                      <div style="font-family:'JetBrains Mono',monospace;font-size:.72rem;color:var(--text-muted);">{dp} / {ph}</div>
                      <div style="font-family:'JetBrains Mono',monospace;font-size:.72rem;font-weight:700;
                        color:{ok_cl};background:{ok_bg};padding:3px 8px;border-radius:4px;text-align:center;">
                        {'PASS' if ok else 'FAIL'}
                      </div>
                    </div>""",unsafe_allow_html=True)
                d1,d2=st.columns(2)
                with d1:
                    st.download_button("⬇ JSON",data=json.dumps(sr["outputs"],indent=2),
                        file_name=f"test_{sr['file'].replace('.vcf','')}.json",mime="application/json",
                        key=f"tsc_{sr['name'][:14]}",use_container_width=True)
                with d2:
                    st.download_button("⬇ VCF",data=load_vcf_file(sr["file"]),
                        file_name=sr["file"],mime="text/plain",
                        key=f"vcf_{sr['name'][:14]}",use_container_width=True)

        st.download_button("⬇ Download Full Test Suite JSON",
            data=json.dumps([{"scenario":s["name"],"pass":s["pass"],"results":s["outputs"]} for s in results],indent=2),
            file_name=f"pharmaguard_tests_{datetime.now().strftime('%Y%m%d_%H%M')}.json",
            mime="application/json",use_container_width=True)
    else:
        st.markdown("""
        <div class="empty-state">
          <span class="empty-icon">▷</span>
          <div class="empty-title">One-click validation</div>
          <div class="empty-hint">
            4 clinical scenarios · ThreadPoolExecutor parallel execution<br>
            Mixed · UltraRapid · All Normal · Worst Case All PM<br>
            Pass/fail per drug · expected vs actual comparison
          </div>
        </div>""",unsafe_allow_html=True)