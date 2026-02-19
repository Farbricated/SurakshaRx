"""
PharmaGuard v9.1 — Precision Clinical Pharmacogenomics
Fixes applied vs v9.0:
  BUG 1: AI badge showed raw "static-template-v5 (rate-limited)" → now "Static Template"
  BUG 2: Clinical Note pre text was invisible → explicit color on pre and .note-box pre
  BUG 3: File uploader filename/size text invisible → color overrides for inner elements
"""

import streamlit as st
import json, uuid, os, re, io
import pandas as pd
from datetime import datetime
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

RISK_CFG = {
    "Safe":         {"color":"#16A34A","bg":"#F0FDF4","border":"#BBF7D0","text":"#14532D","tag_bg":"#DCFCE7","tag_text":"#15803D","shape":"●","severity_dot":"#16A34A"},
    "Adjust Dosage":{"color":"#D97706","bg":"#FFFBEB","border":"#FDE68A","text":"#78350F","tag_bg":"#FEF3C7","tag_text":"#92400E","shape":"▲","severity_dot":"#D97706"},
    "Toxic":        {"color":"#B91C1C","bg":"#FEF2F2","border":"#FECACA","text":"#7F1D1D","tag_bg":"#FEE2E2","tag_text":"#991B1B","shape":"⬛","severity_dot":"#DC2626"},
    "Ineffective":  {"color":"#6D28D9","bg":"#F5F3FF","border":"#DDD6FE","text":"#4C1D95","tag_bg":"#EDE9FE","tag_text":"#5B21B6","shape":"◆","severity_dot":"#7C3AED"},
    "Unknown":      {"color":"#475569","bg":"#F8FAFC","border":"#E2E8F0","text":"#334155","tag_bg":"#F1F5F9","tag_text":"#475569","shape":"?","severity_dot":"#64748B"},
}

SEV_CFG = {
    "none":     {"color":"#16A34A","bg":"#F0FDF4","border":"#BBF7D0","text":"#14532D","label":"None"},
    "low":      {"color":"#D97706","bg":"#FFFBEB","border":"#FDE68A","text":"#78350F","label":"Low"},
    "moderate": {"color":"#EA580C","bg":"#FFF7ED","border":"#FED7AA","text":"#7C2D12","label":"Moderate"},
    "high":     {"color":"#DC2626","bg":"#FEF2F2","border":"#FECACA","text":"#7F1D1D","label":"High"},
    "critical": {"color":"#B91C1C","bg":"#FFF1F1","border":"#FCA5A5","text":"#450A0A","label":"Critical"},
}

PHENO_CFG = {
    "PM":      {"bg":"#FEF2F2","border":"#FECACA","text":"#7F1D1D","bar":"#DC2626","label":"Poor Metabolizer","pct":5},
    "IM":      {"bg":"#FFFBEB","border":"#FDE68A","text":"#78350F","bar":"#D97706","label":"Intermediate Metabolizer","pct":45},
    "NM":      {"bg":"#F0FDF4","border":"#BBF7D0","text":"#14532D","bar":"#16A34A","label":"Normal Metabolizer","pct":100},
    "RM":      {"bg":"#EFF6FF","border":"#BFDBFE","text":"#1E3A8A","bar":"#2563EB","label":"Rapid Metabolizer","pct":115},
    "URM":     {"bg":"#FFF7ED","border":"#FED7AA","text":"#7C2D12","bar":"#EA580C","label":"Ultrarapid Metabolizer","pct":130},
    "Unknown": {"bg":"#F8FAFC","border":"#E2E8F0","text":"#475569","bar":"#94A3B8","label":"Unknown","pct":0},
}

POP_FREQ = {
    "CYP2D6":  {"PM":7,"IM":10,"NM":77,"URM":6},
    "CYP2C19": {"PM":3,"IM":26,"NM":52,"RM":13,"URM":6},
    "CYP2C9":  {"PM":1,"IM":10,"NM":89},
    "SLCO1B1": {"PM":1,"IM":15,"NM":84},
    "TPMT":    {"PM":0.3,"IM":10,"NM":90},
    "DPYD":    {"PM":0.2,"IM":3,"NM":97},
}

CHROM_INFO = {
    "CYP2D6":  {"chrom":"22","band":"q13.2","pos_mb":42.5},
    "CYP2C19": {"chrom":"10","band":"q23.33","pos_mb":96.7},
    "CYP2C9":  {"chrom":"10","band":"q23.33","pos_mb":96.4},
    "SLCO1B1": {"chrom":"12","band":"p12.1","pos_mb":21.3},
    "TPMT":    {"chrom":"6","band":"p22.3","pos_mb":18.1},
    "DPYD":    {"chrom":"1","band":"p22.1","pos_mb":97.5},
}
CHROM_LEN = {"1":248.9,"6":170.8,"10":133.8,"12":133.3,"22":50.8}

PLAIN_PHENO = {
    "PM":"Your body barely processes this medicine",
    "IM":"Your body processes this medicine slower than average",
    "NM":"Your body processes this medicine normally",
    "RM":"Your body processes this medicine slightly faster than average",
    "URM":"Your body processes this medicine dangerously fast",
    "Unknown":"Gene function unclear",
}

PLAIN_RISK = {
    ("CODEINE","PM"):      "Your body can't convert codeine into a painkiller — it won't help your pain.",
    ("CODEINE","URM"):     "Your body converts codeine to morphine extremely fast. Even one tablet could be life-threatening.",
    ("CODEINE","IM"):      "Codeine may be less effective. Your doctor may need to try a different painkiller.",
    ("CODEINE","NM"):      "Codeine works normally for you. Standard doses should manage pain safely.",
    ("WARFARIN","PM"):     "Warfarin stays in your body much longer than normal. Standard doses could cause dangerous bleeding.",
    ("WARFARIN","IM"):     "Warfarin clears more slowly. You'll likely need a lower dose.",
    ("WARFARIN","NM"):     "Warfarin works normally for you. Standard INR monitoring applies.",
    ("CLOPIDOGREL","PM"):  "This heart medication won't activate properly, leaving you unprotected against blood clots.",
    ("CLOPIDOGREL","IM"):  "This heart medication activates less than normal. A stronger alternative may be needed.",
    ("CLOPIDOGREL","NM"):  "This heart medication works normally for you.",
    ("SIMVASTATIN","PM"):  "This cholesterol drug can't be cleared properly and may build up in your muscles, causing serious damage.",
    ("SIMVASTATIN","IM"):  "This cholesterol drug clears more slowly. A lower dose will protect your muscles.",
    ("SIMVASTATIN","NM"):  "This cholesterol drug works normally for you.",
    ("AZATHIOPRINE","PM"): "This immune drug builds up to dangerous levels. Standard doses would seriously harm your bone marrow.",
    ("AZATHIOPRINE","IM"): "You need a lower dose of this immune drug to stay safe.",
    ("AZATHIOPRINE","NM"): "This immune drug works normally for you.",
    ("FLUOROURACIL","PM"): "Your body cannot break down this chemotherapy. Standard doses would be life-threatening.",
    ("FLUOROURACIL","IM"): "This chemotherapy breaks down too slowly. You need a significantly reduced dose.",
    ("FLUOROURACIL","NM"): "This chemotherapy works at a normal rate in your body.",
}

PERSONAS = {
    "A":{"label":"Critical Risk","file":"patient_a_critical.vcf","drugs":["CODEINE","FLUOROURACIL","AZATHIOPRINE"],"desc":"CYP2D6 PM · DPYD PM · TPMT PM","sev":"critical"},
    "B":{"label":"Warfarin PM","file":"patient_b_warfarin.vcf","drugs":["WARFARIN"],"desc":"CYP2C9 *2/*3 Poor Metabolizer","sev":"high"},
    "C":{"label":"Drug Interaction","file":"patient_c_interaction.vcf","drugs":["CLOPIDOGREL"],"desc":"CYP2C19 *2/*3 Poor Metabolizer","sev":"high"},
    "D":{"label":"All Safe","file":"patient_d_safe.vcf","drugs":["CODEINE","WARFARIN","SIMVASTATIN"],"desc":"Wildtype *1/*1 all genes","sev":"none"},
}

TEST_SUITE = [
    {"name":"Mixed Variants","file":"sample.vcf","drugs":["CLOPIDOGREL","CODEINE","AZATHIOPRINE"],
     "expected":{"CLOPIDOGREL":"Ineffective","CODEINE":"Ineffective","AZATHIOPRINE":"Toxic"},
     "desc":"CYP2C19 *2/*3 · CYP2D6 *4/*4 · TPMT *3B/*3C"},
    {"name":"UltraRapid Metabolizer","file":"test_ultrarapid_metabolizer.vcf","drugs":["CODEINE","CLOPIDOGREL"],
     "expected":{"CODEINE":"Toxic","CLOPIDOGREL":"Safe"},"desc":"CYP2D6 *1xN/*1xN → URM → Codeine Toxic"},
    {"name":"All Normal Wild-type","file":"test_all_normal_wildtype.vcf","drugs":ALL_DRUGS,
     "expected":{d:"Safe" for d in ALL_DRUGS},"desc":"Wild-type *1/*1 across all 6 genes"},
    {"name":"Worst Case — All PM","file":"test_worst_case_all_pm.vcf","drugs":ALL_DRUGS,
     "expected":{"CODEINE":"Ineffective","CLOPIDOGREL":"Ineffective","WARFARIN":"Adjust Dosage","SIMVASTATIN":"Toxic","AZATHIOPRINE":"Toxic","FLUOROURACIL":"Toxic"},
     "desc":"Loss-of-function alleles across all 6 genes"},
]

# ── Page Config ───────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="PharmaGuard — Pharmacogenomic Risk",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="collapsed"
)

# ══════════════════════════════════════════════════════════════════════════════
# CSS — v9.1 (all 3 visibility bugs fixed inline)
# ══════════════════════════════════════════════════════════════════════════════
st.markdown("""<style>
@import url('https://fonts.googleapis.com/css2?family=DM+Sans:ital,opsz,wght@0,9..40,300;0,9..40,400;0,9..40,500;0,9..40,600;0,9..40,700;1,9..40,400&family=JetBrains+Mono:wght@400;500;600&display=swap');

:root {
  --sp-1:4px;--sp-2:8px;--sp-3:12px;--sp-4:16px;--sp-5:20px;--sp-6:24px;--sp-8:32px;--sp-10:40px;--sp-12:48px;--sp-16:64px;
  --bg:#FAFBFC;--surface:#FFFFFF;--surface-sub:#F1F5F9;--surface-sub2:#E8EEF5;
  --border-light:#E8EDF5;--border:#D1D9E6;--border-dark:#94A3B8;
  --text-primary:#0F172A;--text-secondary:#334155;--text-muted:#64748B;--text-xmuted:#94A3B8;
  --brand:#1D4ED8;--brand-dark:#1E3A8A;--brand-hover:#1E40AF;--brand-light:#EFF6FF;--brand-border:#BFDBFE;
  --safe:#16A34A;--safe-bg:#F0FDF4;--safe-border:#BBF7D0;
  --warn:#D97706;--warn-bg:#FFFBEB;--danger:#B91C1C;--danger-bg:#FEF2F2;--danger-light:#DC2626;
  --shadow-xs:0 1px 2px rgba(15,23,42,.04);--shadow-sm:0 1px 3px rgba(15,23,42,.06),0 1px 2px rgba(15,23,42,.04);
  --shadow-md:0 4px 12px rgba(15,23,42,.08),0 2px 4px rgba(15,23,42,.04);
  --shadow-lg:0 12px 32px rgba(15,23,42,.10),0 4px 8px rgba(15,23,42,.04);
  --r-sm:6px;--r-md:8px;--r-lg:12px;--r-xl:16px;--r-2xl:20px;--r-full:9999px;
  --font-body:'DM Sans',-apple-system,BlinkMacSystemFont,sans-serif;
  --font-mono:'JetBrains Mono','Fira Code',monospace;
}

*,*::before,*::after{box-sizing:border-box;margin:0;padding:0;}
html,body,[class*="css"]{font-family:var(--font-body)!important;font-size:16px!important;background:var(--bg)!important;color:var(--text-primary)!important;-webkit-font-smoothing:antialiased!important;}
.stApp{background:var(--bg)!important;}
.main .block-container{padding:0 var(--sp-10) var(--sp-16)!important;max-width:1280px!important;}
#MainMenu,footer,header{visibility:hidden;}

@keyframes fade-up{from{opacity:0;transform:translateY(12px)}to{opacity:1;transform:translateY(0)}}
@keyframes fade-in{from{opacity:0}to{opacity:1}}
@keyframes slide-in{from{opacity:0;transform:translateX(-8px)}to{opacity:1;transform:translateX(0)}}
@keyframes pulse-once{0%{transform:scale(1);box-shadow:0 0 0 0 rgba(185,28,28,.3)}40%{transform:scale(1.02);box-shadow:0 0 0 8px rgba(185,28,28,0)}100%{transform:scale(1);box-shadow:0 0 0 0 rgba(185,28,28,0)}}
@keyframes bar-fill{from{width:0!important}}
@keyframes score-count{from{opacity:0;transform:scale(.85)}to{opacity:1;transform:scale(1)}}

.reveal-card{animation:fade-up .32s cubic-bezier(.4,0,.2,1) both;}
.reveal-card:nth-child(1){animation-delay:.04s}.reveal-card:nth-child(2){animation-delay:.10s}
.reveal-card:nth-child(3){animation-delay:.16s}.reveal-card:nth-child(4){animation-delay:.22s}
.reveal-card:nth-child(5){animation-delay:.28s}.reveal-card:nth-child(6){animation-delay:.34s}

/* Navigation */
.pg-nav{display:flex;align-items:center;justify-content:space-between;padding:var(--sp-6) 0;border-bottom:1px solid var(--border-light);margin-bottom:var(--sp-8);animation:fade-in .4s ease;}
.pg-brand-name{font-size:1.5rem;font-weight:700;color:var(--text-primary);letter-spacing:-.03em;line-height:1;}
.pg-brand-name span{color:var(--brand);}
.pg-brand-sub{font-family:var(--font-mono);font-size:.8rem;color:var(--text-xmuted);letter-spacing:.1em;text-transform:uppercase;}
.pg-nav-badges{display:flex;align-items:center;gap:var(--sp-2);}
.pg-badge{font-family:var(--font-mono);font-size:.8rem;font-weight:500;letter-spacing:.06em;text-transform:uppercase;padding:4px 10px;border-radius:var(--r-full);border:1px solid;white-space:nowrap;}
.pg-badge-default{color:var(--text-muted);border-color:var(--border);background:var(--surface);}
.pg-badge-brand{color:var(--brand);border-color:var(--brand-border);background:var(--brand-light);font-weight:600;}

/* Trust strip */
.trust-strip{display:flex;align-items:center;gap:var(--sp-6);padding:var(--sp-3) var(--sp-4);background:var(--brand-light);border:1px solid var(--brand-border);border-radius:var(--r-md);margin-bottom:var(--sp-8);}
.trust-item{display:flex;align-items:center;gap:var(--sp-2);font-size:.85rem;color:var(--brand-dark);font-weight:500;white-space:nowrap;}
.trust-sep{width:1px;height:16px;background:var(--brand-border);flex-shrink:0;}

/* Tabs */
.stTabs [data-baseweb="tab-list"]{background:transparent!important;border-bottom:1px solid var(--border-light)!important;gap:0!important;padding:0!important;margin-bottom:var(--sp-8)!important;box-shadow:none!important;}
.stTabs [data-baseweb="tab"]{font-family:var(--font-body)!important;font-size:1rem!important;font-weight:500!important;color:var(--text-muted)!important;padding:var(--sp-3) var(--sp-5)!important;background:transparent!important;border:none!important;border-bottom:2.5px solid transparent!important;border-radius:0!important;transition:color .15s!important;}
.stTabs [aria-selected="true"]{color:var(--brand)!important;border-bottom-color:var(--brand)!important;font-weight:600!important;}
.stTabs [data-baseweb="tab-panel"]{padding-top:0!important;}

/* Sidebar */
[data-testid="stSidebar"]{background:var(--surface)!important;border-right:1px solid var(--border-light)!important;}

/* Section label */
.sec-label{display:flex;align-items:center;gap:var(--sp-3);font-size:.8rem;font-weight:600;letter-spacing:.1em;text-transform:uppercase;color:var(--text-muted);margin-bottom:var(--sp-4);}
.sec-label::after{content:'';flex:1;height:1px;background:var(--border-light);}

/* Steps */
.steps{display:flex;background:var(--surface);border:1px solid var(--border-light);border-radius:var(--r-xl);overflow:hidden;margin-bottom:var(--sp-8);box-shadow:var(--shadow-xs);}
.step{flex:1;padding:var(--sp-4) var(--sp-5);border-right:1px solid var(--border-light);}
.step:last-child{border-right:none;}
.step-num{font-family:var(--font-mono);font-size:.7rem;font-weight:600;letter-spacing:.12em;text-transform:uppercase;color:var(--text-xmuted);margin-bottom:3px;}
.step-lbl{font-size:.875rem;font-weight:500;color:var(--text-muted);}
.step.done .step-num{color:var(--brand);}.step.done .step-lbl{color:var(--text-primary);font-weight:600;}.step.done{background:#F8FAFF;}

/* Persona cards */
.persona-grid{display:grid;grid-template-columns:repeat(4,1fr);gap:var(--sp-3);margin-bottom:var(--sp-8);}
.persona-card{background:var(--surface);border:1.5px solid var(--border-light);border-radius:var(--r-lg);padding:var(--sp-4);transition:all .2s cubic-bezier(.4,0,.2,1);box-shadow:var(--shadow-xs);cursor:pointer;}
.persona-card:hover{transform:translateY(-2px);box-shadow:var(--shadow-md);border-color:var(--brand-border);}
.pc-sev{display:inline-flex;align-items:center;gap:5px;font-size:.75rem;font-weight:600;padding:3px 10px;border-radius:var(--r-full);border:1px solid;margin-bottom:var(--sp-2);}
.pc-name{font-size:.875rem;font-weight:700;color:var(--text-primary);margin-bottom:3px;}
.pc-desc{font-family:var(--font-mono);font-size:.7rem;color:var(--text-muted);line-height:1.7;}

/* Risk command center */
.risk-center{border-radius:var(--r-2xl);padding:var(--sp-8);margin-bottom:var(--sp-6);border:1.5px solid;position:relative;overflow:hidden;box-shadow:var(--shadow-md);}
.rc-eyebrow{font-family:var(--font-mono);font-size:.7rem;font-weight:600;letter-spacing:.14em;text-transform:uppercase;opacity:.6;margin-bottom:var(--sp-1);}
.rc-headline{font-size:2.5rem;font-weight:700;letter-spacing:-.03em;line-height:1.1;margin-bottom:var(--sp-1);}
.rc-sub{font-size:.9rem;opacity:.75;margin-bottom:var(--sp-5);}
.rc-stats{display:grid;grid-template-columns:repeat(4,1fr);gap:var(--sp-5);padding-top:var(--sp-5);border-top:1px solid;border-color:inherit;opacity:.8;}
.rc-stat-num{font-size:2rem;font-weight:700;letter-spacing:-.03em;line-height:1;margin-bottom:3px;}
.rc-stat-lbl{font-family:var(--font-mono);font-size:.65rem;font-weight:500;letter-spacing:.1em;text-transform:uppercase;opacity:.6;}

/* Critical alert */
.crit-alert{display:flex;gap:var(--sp-4);background:#FFF1F2;border:1px solid #FECDD3;border-left:4px solid var(--danger);border-radius:var(--r-lg);padding:var(--sp-4) var(--sp-5);margin-bottom:var(--sp-4);animation:pulse-once .8s ease .3s both;box-shadow:var(--shadow-sm);}
.crit-title{font-size:.95rem;font-weight:700;color:var(--danger);margin-bottom:3px;}
.crit-note{font-size:.875rem;color:#7F1D1D;line-height:1.65;margin-bottom:var(--sp-2);}
.crit-action{font-family:var(--font-mono);font-size:.7rem;font-weight:600;color:var(--danger-light);letter-spacing:.08em;text-transform:uppercase;}

/* Gene row */
.gene-row{display:grid;grid-template-columns:repeat(6,1fr);gap:var(--sp-3);margin-bottom:var(--sp-6);}
.gene-box{background:var(--surface);border:1.5px solid var(--border-light);border-radius:var(--r-lg);padding:var(--sp-4) var(--sp-3);text-align:center;box-shadow:var(--shadow-xs);transition:box-shadow .15s,transform .15s,border-color .15s;}
.gene-box:hover{box-shadow:var(--shadow-md);transform:translateY(-2px);}
.gene-nm{font-family:var(--font-mono);font-size:.75rem;font-weight:600;margin-bottom:var(--sp-2);color:var(--text-secondary);}
.gene-track{height:3px;border-radius:2px;background:var(--surface-sub);margin:var(--sp-2) 0;overflow:hidden;}
.gene-fill{height:100%;border-radius:2px;}
.gene-ph{font-family:var(--font-mono);font-size:.8rem;font-weight:600;letter-spacing:.03em;}

/* Drug table */
.dtab{background:var(--surface);border:1px solid var(--border-light);border-radius:var(--r-xl);overflow:hidden;margin-bottom:var(--sp-6);box-shadow:var(--shadow-sm);}
.dtab-head{display:grid;grid-template-columns:1.4fr 1.2fr .9fr 1fr .9fr 1.1fr;background:var(--surface-sub);border-bottom:1px solid var(--border-light);}
.dtab-hcell{font-family:var(--font-mono);font-size:.65rem;font-weight:600;letter-spacing:.1em;text-transform:uppercase;color:var(--text-muted);padding:var(--sp-3) var(--sp-4);}
.dtab-row{display:grid;grid-template-columns:1.4fr 1.2fr .9fr 1fr .9fr 1.1fr;border-bottom:1px solid var(--border-light);transition:background .12s;}
.dtab-row:last-child{border-bottom:none;}.dtab-row:hover{background:var(--surface-sub);}
.dtab-cell{font-size:.9rem;color:var(--text-secondary);padding:var(--sp-3) var(--sp-4);display:flex;align-items:center;}

/* Risk badge */
.risk-badge{display:inline-flex;align-items:center;gap:6px;font-size:.8rem;font-weight:600;padding:4px 12px;border-radius:var(--r-full);border:1.5px solid;letter-spacing:-.01em;}

/* PGx card */
.pgx-card{background:var(--surface);border:1px solid var(--border-light);border-radius:var(--r-2xl);padding:var(--sp-8);margin-bottom:var(--sp-6);box-shadow:var(--shadow-md);position:relative;overflow:hidden;}
.pgx-card::before{content:'';position:absolute;top:0;right:0;width:280px;height:280px;background:radial-gradient(circle at top right,var(--brand-light) 0%,transparent 65%);pointer-events:none;}
.pgx-eyebrow{font-family:var(--font-mono);font-size:.7rem;font-weight:600;letter-spacing:.14em;text-transform:uppercase;color:var(--brand);margin-bottom:var(--sp-2);}
.pgx-score{font-size:4.5rem;font-weight:700;letter-spacing:-.04em;line-height:1;margin-bottom:4px;animation:score-count .5s cubic-bezier(.4,0,.2,1) .1s both;}
.pgx-label{font-size:.9rem;color:var(--text-muted);margin-bottom:var(--sp-5);}
.pgx-marker{position:relative;height:6px;background:var(--surface-sub);border-radius:3px;overflow:visible;margin-bottom:var(--sp-5);}
.pgx-fill{position:absolute;top:0;left:0;height:100%;border-radius:3px;transition:width .9s cubic-bezier(.4,0,.2,1);}
.pgx-indicator{position:absolute;top:-4px;width:14px;height:14px;border-radius:50%;background:var(--surface);border:3px solid;transform:translateX(-50%);box-shadow:var(--shadow-sm);transition:left .9s cubic-bezier(.4,0,.2,1);}
.pgx-thresh-labels{display:flex;justify-content:space-between;font-family:var(--font-mono);font-size:.65rem;color:var(--text-xmuted);margin-bottom:var(--sp-3);}
.pgx-pills{display:flex;flex-wrap:wrap;gap:var(--sp-2);}
.pgx-pill{font-family:var(--font-mono);font-size:.7rem;font-weight:600;padding:3px 10px;border-radius:var(--r-full);border:1px solid;letter-spacing:.03em;}

/* Heatmap */
.hm-wrap{background:var(--surface);border:1px solid var(--border-light);border-radius:var(--r-xl);padding:var(--sp-6);margin-bottom:var(--sp-6);box-shadow:var(--shadow-sm);overflow-x:auto;}
.hm-eyebrow{font-family:var(--font-mono);font-size:.7rem;font-weight:600;letter-spacing:.12em;text-transform:uppercase;color:var(--text-muted);margin-bottom:var(--sp-5);}
.hm-grid{display:grid;gap:3px;}
.hm-cell{border-radius:var(--r-sm);display:flex;flex-direction:column;align-items:center;justify-content:center;padding:var(--sp-3) var(--sp-2);min-height:56px;border:1.5px solid;transition:transform .12s,box-shadow .12s;cursor:default;}
.hm-cell:hover{transform:scale(1.06);box-shadow:var(--shadow-md);z-index:5;position:relative;}
.hm-cell-name{font-family:var(--font-mono);font-size:.7rem;font-weight:600;margin-bottom:2px;}
.hm-cell-risk{font-family:var(--font-mono);font-size:.65rem;opacity:.8;}
.hm-header{font-family:var(--font-mono);font-size:.65rem;letter-spacing:.05em;color:var(--text-muted);display:flex;align-items:center;justify-content:center;min-height:56px;}
.hm-legend{display:flex;gap:var(--sp-5);margin-top:var(--sp-4);flex-wrap:wrap;}
.hm-legend-item{font-family:var(--font-mono);font-size:.7rem;display:flex;align-items:center;gap:5px;color:var(--text-muted);}
.hm-dot{width:10px;height:10px;border-radius:3px;display:inline-block;border:1.5px solid;}

/* Chromosome */
.chrom-wrap{background:var(--surface);border:1px solid var(--border-light);border-radius:var(--r-xl);padding:var(--sp-5) var(--sp-6);box-shadow:var(--shadow-sm);}
.chrom-eyebrow{font-family:var(--font-mono);font-size:.7rem;font-weight:600;letter-spacing:.12em;text-transform:uppercase;color:var(--text-muted);margin-bottom:var(--sp-4);}
.chrom-row{display:flex;align-items:center;gap:var(--sp-3);margin-bottom:var(--sp-2);}
.chrom-chr{font-family:var(--font-mono);font-size:.75rem;color:var(--text-muted);width:18px;text-align:right;flex-shrink:0;}
.chrom-bar{flex:1;height:11px;background:var(--surface-sub);border-radius:6px;position:relative;overflow:visible;border:1px solid var(--border-light);}
.chrom-body{position:absolute;inset:0;background:linear-gradient(90deg,#DDE3EE,#EEF2F8,#DDE3EE);border-radius:6px;}
.chrom-marker{position:absolute;top:-5px;width:3px;height:21px;border-radius:2px;transform:translateX(-50%);}
.chrom-gene{font-family:var(--font-mono);font-size:.75rem;color:var(--text-secondary);width:56px;flex-shrink:0;font-weight:500;}
.chrom-band{font-family:var(--font-mono);font-size:.65rem;color:var(--text-xmuted);}

/* Population freq */
.pop-wrap{background:var(--surface);border:1px solid var(--border-light);border-radius:var(--r-lg);padding:var(--sp-4) var(--sp-5);margin-bottom:var(--sp-4);box-shadow:var(--shadow-xs);}
.pop-eyebrow{font-family:var(--font-mono);font-size:.65rem;font-weight:600;letter-spacing:.1em;text-transform:uppercase;color:var(--text-muted);margin-bottom:var(--sp-3);}
.pop-row{display:flex;align-items:center;gap:var(--sp-3);margin-bottom:var(--sp-2);}
.pop-ph{font-family:var(--font-mono);font-size:.75rem;color:var(--text-secondary);width:96px;flex-shrink:0;font-weight:500;}
.pop-track{flex:1;height:4px;background:var(--surface-sub);border-radius:2px;overflow:hidden;}
.pop-fill{height:100%;border-radius:2px;animation:bar-fill .7s cubic-bezier(.4,0,.2,1) both;}
.pop-pct{font-family:var(--font-mono);font-size:.7rem;width:32px;text-align:right;color:var(--text-muted);}
.pop-you{font-family:var(--font-mono);font-size:.65rem;color:var(--brand);font-weight:700;margin-left:3px;}

/* Interaction matrix */
.ix-grid{display:grid;gap:3px;}
.ix-cell{border-radius:var(--r-sm);display:flex;align-items:center;justify-content:center;min-height:44px;font-family:var(--font-mono);font-size:.65rem;text-align:center;padding:var(--sp-1);font-weight:700;border:1.5px solid;transition:transform .12s,box-shadow .12s;}
.ix-cell:hover{transform:scale(1.06);z-index:5;position:relative;box-shadow:var(--shadow-md);}
.ix-head{font-family:var(--font-mono);font-size:.65rem;letter-spacing:.05em;color:var(--text-muted);display:flex;align-items:center;justify-content:center;min-height:44px;}

/* Drug card */
.dcard{background:var(--surface);border:1px solid var(--border-light);border-radius:var(--r-2xl);margin-bottom:var(--sp-6);overflow:hidden;box-shadow:var(--shadow-sm);transition:box-shadow .2s;}
.dcard:hover{box-shadow:var(--shadow-md);}
.dcard-header{display:flex;align-items:center;justify-content:space-between;padding:var(--sp-5) var(--sp-6);border-bottom:1px solid var(--border-light);}
.dcard-left{display:flex;align-items:center;gap:var(--sp-4);}
.dcard-indicator{width:10px;height:10px;border-radius:50%;flex-shrink:0;}
.dcard-drug{font-size:1.125rem;font-weight:700;letter-spacing:-.02em;color:var(--text-primary);}
.dcard-meta{font-family:var(--font-mono);font-size:.75rem;color:var(--text-muted);margin-top:3px;}
.dcard-body{padding:var(--sp-6);}

/* Metrics row */
.metrics-row{display:grid;grid-template-columns:repeat(4,1fr);gap:1px;background:var(--border-light);border-radius:var(--r-lg);overflow:hidden;border:1px solid var(--border-light);margin-bottom:var(--sp-5);}
.metric-cell{background:var(--surface-sub);padding:var(--sp-4);}
.metric-key{font-family:var(--font-mono);font-size:.65rem;font-weight:600;letter-spacing:.1em;text-transform:uppercase;color:var(--text-muted);margin-bottom:4px;}
.metric-val{font-size:1.125rem;font-weight:700;color:var(--text-primary);letter-spacing:-.02em;}

.conf-grid{display:grid;grid-template-columns:1fr 1fr;gap:var(--sp-5);margin-bottom:var(--sp-5);}
.conf-label{font-family:var(--font-mono);font-size:.65rem;font-weight:600;letter-spacing:.08em;text-transform:uppercase;color:var(--text-muted);display:flex;justify-content:space-between;margin-bottom:5px;}
.conf-track{height:4px;background:var(--surface-sub);border-radius:2px;overflow:hidden;}
.conf-fill{height:100%;border-radius:2px;animation:bar-fill .7s cubic-bezier(.4,0,.2,1) both;}

/* Variant table */
.vtable{width:100%;border-collapse:collapse;}
.vtable th{font-family:var(--font-mono);font-size:.65rem;font-weight:600;letter-spacing:.1em;text-transform:uppercase;color:var(--text-muted);padding:0 var(--sp-3) var(--sp-3);text-align:left;border-bottom:1px solid var(--border-light);}
.vtable td{font-family:var(--font-mono);font-size:.85rem;color:var(--text-secondary);padding:var(--sp-2) var(--sp-3);border-bottom:1px solid var(--border-light);}
.vtable tbody tr:last-child td{border-bottom:none;}.vtable tbody tr:hover td{background:var(--surface-sub);}
.v-rsid{color:#2563EB!important;font-weight:500!important;}.v-star{color:#7C3AED!important;font-weight:500!important;}
.v-nofunc{color:var(--danger)!important;font-weight:500!important;}.v-dec{color:var(--warn)!important;font-weight:500!important;}.v-norm{color:var(--safe)!important;font-weight:500!important;}

/* Recommendation box */
.rec-box{border-radius:var(--r-lg);border:1.5px solid;padding:var(--sp-4) var(--sp-5);margin-bottom:var(--sp-4);}
.rec-label{font-family:var(--font-mono);font-size:.65rem;font-weight:600;letter-spacing:.1em;text-transform:uppercase;margin-bottom:var(--sp-2);}
.rec-text{font-size:.95rem;line-height:1.75;color:var(--text-secondary);}
.alt-chips{display:flex;flex-wrap:wrap;gap:var(--sp-2);}
.alt-chip{font-family:var(--font-mono);font-size:.75rem;font-weight:500;color:var(--brand);background:var(--brand-light);border:1px solid var(--brand-border);border-radius:var(--r-full);padding:4px 12px;}
.cpic-badge{font-family:var(--font-mono);font-size:.65rem;font-weight:700;background:#FEFCE8;border:1px solid #FDE047;color:#713F12;padding:2px 8px;border-radius:4px;display:inline-block;margin-left:var(--sp-2);vertical-align:middle;letter-spacing:.05em;}

/* AI block */
.ai-block{background:linear-gradient(135deg,#F8FBFF 0%,#EFF6FF 100%);border:1.5px solid var(--brand-border);border-radius:var(--r-xl);overflow:hidden;margin-bottom:var(--sp-5);box-shadow:0 2px 8px rgba(29,78,216,.06);}
.ai-header{display:flex;align-items:center;gap:var(--sp-3);padding:var(--sp-3) var(--sp-5);background:rgba(255,255,255,.7);border-bottom:1px solid var(--brand-border);}
/* BUG 1 FIX: badge always shows clean label — raw model strings now cleaned in Python */
.ai-badge-pill{font-family:var(--font-mono);font-size:.7rem;font-weight:600;letter-spacing:.08em;text-transform:uppercase;background:var(--brand-light);border:1px solid var(--brand-border);color:var(--brand);padding:3px 9px;border-radius:var(--r-sm);}
.ai-title{font-size:.9rem;font-weight:600;color:var(--brand-dark);}
.ai-section{padding:var(--sp-4) var(--sp-5);border-bottom:1px solid rgba(191,219,254,.5);}
.ai-section:last-child{border-bottom:none;}.ai-section:hover{background:rgba(255,255,255,.5);}
.ai-sec-label{font-family:var(--font-mono);font-size:.65rem;font-weight:600;letter-spacing:.1em;text-transform:uppercase;color:var(--brand);margin-bottom:var(--sp-2);}
.ai-sec-text{font-size:.9rem;line-height:1.8;color:var(--text-secondary);}

/* Narrative */
.narrative-box{background:var(--brand-light);border:1.5px solid var(--brand-border);border-radius:var(--r-xl);padding:var(--sp-6);margin-bottom:var(--sp-6);box-shadow:var(--shadow-sm);}
.narrative-header{display:flex;align-items:center;gap:var(--sp-3);margin-bottom:var(--sp-4);}
.narrative-text{font-size:.95rem;line-height:1.85;color:var(--text-secondary);}

/* Before/After */
.ba-grid{display:grid;grid-template-columns:1fr 1fr;border:1.5px solid var(--border-light);border-radius:var(--r-xl);overflow:hidden;margin-bottom:var(--sp-6);box-shadow:var(--shadow-md);}
.ba-side{padding:var(--sp-6);}
.ba-side-lbl{font-family:var(--font-mono);font-size:.7rem;font-weight:700;letter-spacing:.1em;text-transform:uppercase;margin-bottom:var(--sp-3);}
.ba-drug{font-size:.9rem;font-weight:700;margin-bottom:3px;}
.ba-text{font-size:.875rem;line-height:1.65;}
.ba-gene{font-family:var(--font-mono);font-size:.7rem;margin-top:var(--sp-3);opacity:.65;}

/* Rx checker */
.rx-result{border-radius:var(--r-lg);padding:var(--sp-5) var(--sp-6);margin-top:var(--sp-4);border:1.5px solid;animation:fade-up .25s ease;box-shadow:var(--shadow-md);}
.rx-verdict{font-size:.95rem;font-weight:700;margin-bottom:var(--sp-2);}
.rx-detail{font-size:.875rem;line-height:1.7;color:var(--text-secondary);margin-bottom:var(--sp-2);}
.rx-meta{font-family:var(--font-mono);font-size:.7rem;letter-spacing:.06em;text-transform:uppercase;}

/* BUG 2 FIX: Clinical note — pre text was invisible.
   Root cause: global pre only set font-family with no color. Streamlit's dark
   container backgrounds caused inherited text to be near-invisible.
   Fix: explicit color on .note-box pre AND the global pre rule. */
.note-box{background:var(--surface-sub);border:1px solid var(--border-light);border-radius:var(--r-xl);padding:var(--sp-6);}
.note-box pre,.note-box pre *{
  font-family:var(--font-mono)!important;
  font-size:.85rem!important;
  color:var(--text-secondary)!important;
  background:transparent!important;
  line-height:1.85!important;
  white-space:pre-wrap!important;
  word-break:break-word!important;
  font-weight:400!important;
}

/* Patient mode */
.patient-banner{border-radius:var(--r-xl);padding:var(--sp-6);margin-bottom:var(--sp-6);border:1.5px solid;box-shadow:var(--shadow-md);}
.patient-banner-title{font-size:1.125rem;font-weight:700;margin-bottom:var(--sp-2);}
.patient-banner-sub{font-size:.9rem;line-height:1.7;}
.pcard{background:var(--surface);border:1.5px solid;border-radius:var(--r-xl);padding:var(--sp-6);margin-bottom:var(--sp-4);box-shadow:var(--shadow-sm);animation:fade-up .3s ease both;transition:box-shadow .2s;}
.pcard:hover{box-shadow:var(--shadow-md);}
.pcard-drug{font-size:1.1rem;font-weight:700;letter-spacing:-.02em;margin-bottom:3px;}
.pcard-verdict{font-size:.9rem;font-weight:600;line-height:1.5;margin-bottom:var(--sp-2);}
.pcard-gene{font-family:var(--font-mono);font-size:.7rem;letter-spacing:.06em;color:var(--text-muted);margin-bottom:var(--sp-3);font-weight:600;text-transform:uppercase;}
.pcard-plain{font-size:.9rem;line-height:1.8;color:var(--text-secondary);}
.pcard-action{display:flex;align-items:flex-start;gap:var(--sp-3);background:var(--surface-sub);border:1px solid var(--border-light);border-radius:var(--r-lg);padding:var(--sp-4);margin-top:var(--sp-4);}
.pcard-action-text{font-size:.875rem;color:var(--text-primary);line-height:1.65;}

/* Disclaimer */
.disclaimer-box{display:flex;gap:var(--sp-4);background:#FFFBEB;border:1px solid #FDE68A;border-left:3px solid var(--warn);border-radius:var(--r-lg);padding:var(--sp-4) var(--sp-5);margin-bottom:var(--sp-6);}
.disclaimer-text{font-size:.85rem;color:#78350F;line-height:1.7;}

/* Test suite card */
.tc-card{background:var(--surface);border:1px solid var(--border-light);border-radius:var(--r-xl);padding:var(--sp-5);box-shadow:var(--shadow-xs);}
.tc-name{font-size:.95rem;font-weight:700;color:var(--text-primary);margin-bottom:var(--sp-1);}
.tc-desc{font-family:var(--font-mono);font-size:.7rem;color:var(--text-muted);margin-bottom:var(--sp-4);line-height:1.7;}

/* Empty state */
.empty-state{text-align:center;padding:5rem 2rem;border:1.5px dashed var(--border);border-radius:var(--r-2xl);background:var(--surface);box-shadow:var(--shadow-xs);}
.empty-icon{font-size:2.5rem;display:block;margin-bottom:var(--sp-4);opacity:.3;}
.empty-title{font-size:1.125rem;font-weight:600;color:var(--text-secondary);margin-bottom:var(--sp-2);}
.empty-hint{font-family:var(--font-mono);font-size:.7rem;color:var(--text-muted);letter-spacing:.04em;line-height:2.4;}

.info-strip{display:flex;align-items:flex-start;gap:var(--sp-3);background:var(--brand-light);border:1px solid var(--brand-border);border-radius:var(--r-md);padding:var(--sp-4);margin-bottom:var(--sp-4);}
.info-strip-text{font-size:.875rem;color:var(--brand-dark);line-height:1.65;}

/* Buttons */
.stButton>button{background:var(--brand)!important;color:#FFF!important;border:none!important;border-radius:var(--r-md)!important;font-family:var(--font-body)!important;font-weight:600!important;font-size:.95rem!important;padding:.6875rem 1.75rem!important;height:48px!important;transition:all .15s!important;box-shadow:0 1px 2px rgba(29,78,216,.2)!important;min-height:44px!important;}
.stButton>button:hover{background:var(--brand-hover)!important;box-shadow:0 4px 12px rgba(29,78,216,.3)!important;transform:translateY(-1px)!important;}
.stDownloadButton>button{background:var(--surface)!important;color:var(--text-secondary)!important;border:1.5px solid var(--border)!important;border-radius:var(--r-md)!important;font-family:var(--font-mono)!important;font-size:.8rem!important;padding:.5rem 1rem!important;transition:all .15s!important;box-shadow:var(--shadow-xs)!important;min-height:44px!important;}
.stDownloadButton>button:hover{background:var(--brand)!important;color:#FFF!important;border-color:var(--brand)!important;}

/* BUG 3 FIX: File uploader filename + size text was invisible.
   Root cause: .stFileUploader>div set background:var(--surface) with no text color.
   Streamlit's injected dark label styles bled through making text invisible.
   Fix: explicit color overrides on every text element inside the uploader. */
.stFileUploader>div{border-radius:var(--r-xl)!important;border:2px dashed var(--border)!important;background:var(--surface)!important;transition:border-color .15s,background .15s!important;padding:var(--sp-8)!important;}
.stFileUploader>div:hover{border-color:var(--brand)!important;background:var(--brand-light)!important;}
.stFileUploader label,.stFileUploader p,.stFileUploader small,.stFileUploader span,
[data-testid="stFileUploader"] label,[data-testid="stFileUploader"] p,
[data-testid="stFileUploader"] small,[data-testid="stFileUploader"] span{
  color:var(--text-secondary)!important;font-family:var(--font-body)!important;
}
[data-testid="stFileUploaderFile"] span,
[data-testid="stFileUploaderFile"] p,
[data-testid="stFileUploaderFile"] small{color:var(--text-primary)!important;font-weight:500!important;}

/* Inputs */
.stTextInput>div>div>input{border-radius:var(--r-md)!important;border:1.5px solid var(--border)!important;background:var(--surface)!important;color:var(--text-primary)!important;font-family:var(--font-body)!important;font-size:.95rem!important;padding:.6875rem .875rem!important;height:48px!important;transition:border-color .15s,box-shadow .15s!important;box-shadow:var(--shadow-xs)!important;}
.stTextInput>div>div>input:focus{border-color:var(--brand)!important;box-shadow:0 0 0 3px rgba(147,197,253,.4)!important;outline:none!important;}
.stSelectbox [data-baseweb="select"]>div,.stMultiSelect [data-baseweb="select"]>div{border-radius:var(--r-md)!important;border:1.5px solid var(--border)!important;background:var(--surface)!important;min-height:48px!important;}
.stMultiSelect span[data-baseweb="tag"]{background:var(--brand-light)!important;color:var(--brand-dark)!important;border:1px solid var(--brand-border)!important;font-family:var(--font-mono)!important;font-size:.75rem!important;border-radius:5px!important;}
div[data-testid="stExpander"]{background:var(--surface)!important;border:1px solid var(--border-light)!important;border-radius:var(--r-lg)!important;box-shadow:var(--shadow-xs)!important;margin-bottom:var(--sp-2)!important;}
div[data-testid="stExpander"] summary{font-family:var(--font-body)!important;font-size:.9rem!important;font-weight:500!important;color:var(--text-secondary)!important;padding:var(--sp-4) var(--sp-5)!important;}
[data-testid="stMetric"]{background:var(--surface)!important;border:1px solid var(--border-light)!important;border-radius:var(--r-xl)!important;padding:var(--sp-5)!important;box-shadow:var(--shadow-xs)!important;}
[data-testid="stMetricLabel"]{font-family:var(--font-mono)!important;font-size:.7rem!important;color:var(--text-muted)!important;text-transform:uppercase!important;letter-spacing:.1em!important;font-weight:600!important;}
[data-testid="stMetricValue"]{font-size:1.875rem!important;color:var(--text-primary)!important;font-weight:700!important;letter-spacing:-.02em!important;}

/* BUG 2 FIX also: global pre must set explicit color (was font-family only) */
.stCode{border-radius:var(--r-lg)!important;}
pre{font-family:var(--font-mono)!important;color:var(--text-secondary)!important;white-space:pre-wrap!important;word-break:break-word!important;}
.stSpinner>div{border-color:var(--brand) transparent transparent transparent!important;}
</style>""", unsafe_allow_html=True)


# ══════════════════════════════════════════════════════════════════════════════
# UTILITY FUNCTIONS
# ══════════════════════════════════════════════════════════════════════════════

def load_vcf(filename):
    p = os.path.join(BASE_DIR, "sample_data", filename)
    return open(p).read() if os.path.exists(p) else get_sample_vcf()

def run_pipeline(vcf, drugs, pid, key, run_ix=True, gen_pdf=True, skip_llm=False):
    parsed  = parse_vcf(vcf)
    results = run_risk_assessment(parsed, drugs)
    results = generate_all_explanations(key, results, skip_llm=skip_llm)
    outputs = [build_output_schema(patient_id=pid, drug=r["drug"], result=r,
                parsed_vcf=parsed, llm_exp=r.get("llm_explanation", {})) for r in results]
    ix  = run_interaction_analysis(drugs, results) if run_ix and len(drugs) > 1 else None
    pdf = None
    if gen_pdf:
        try:
            pdf = generate_pdf_report(pid, outputs, parsed)
        except:
            pass
    return parsed, results, outputs, ix, pdf

def func_cls(status):
    s = (status or "").lower()
    if "no_function" in s or "no function" in s:
        return "v-nofunc"
    if any(x in s for x in ["decreased","splice","missense","frame","stop","pathogenic"]) and "synonymous" not in s:
        return "v-dec"
    return "v-norm"

def sec(label):
    st.markdown(f'<div class="sec-label">{label}</div>', unsafe_allow_html=True)

def risk_badge_html(rl):
    rc = RISK_CFG.get(rl, RISK_CFG["Unknown"])
    return (f'<span class="risk-badge" style="background:{rc["tag_bg"]};color:{rc["tag_text"]};'
            f'border-color:{rc["border"]};">'
            f'<span style="font-size:.8rem;">{rc["shape"]}</span>{rl}</span>')

def clean_model_label(raw_model: str):
    """
    BUG 1 FIX: Return (display_label, is_static) with a clean human-readable label.
    Strips parenthetical suffixes like '(rate-limited)', '(no API key)', '(error: ...)'.
    """
    is_static = "static" in raw_model.lower()
    if is_static:
        return "Static Template", True
    clean = re.sub(r"\s*\(.*?\)$", "", raw_model).strip()
    return (clean or raw_model), False


# ══════════════════════════════════════════════════════════════════════════════
# VISUAL COMPONENTS
# ══════════════════════════════════════════════════════════════════════════════

def compute_pgx(outputs):
    SEV_S  = {"none": 0, "low": 20, "moderate": 45, "high": 70, "critical": 100}
    RISK_S = {"Safe": 0, "Adjust Dosage": 35, "Toxic": 85, "Ineffective": 70, "Unknown": 20}
    W      = {"FLUOROURACIL": 1.4, "AZATHIOPRINE": 1.3, "CLOPIDOGREL": 1.3,
               "WARFARIN": 1.2, "CODEINE": 1.1, "SIMVASTATIN": 1.0}
    if not outputs:
        return 0, "No data", []
    tw = ws = 0
    bd = []
    for o in outputs:
        drug = o["drug"]
        sev  = o["risk_assessment"]["severity"]
        rl   = o["risk_assessment"]["risk_label"]
        gene = o["pharmacogenomic_profile"]["primary_gene"]
        ph   = o["pharmacogenomic_profile"]["phenotype"]
        sc   = (SEV_S.get(sev, 0) + RISK_S.get(rl, 0)) / 2
        wt   = W.get(drug, 1.0)
        ws  += sc * wt
        tw  += wt
        bd.append((gene, drug, ph, rl, sc))
    final = min(100, int(ws / tw)) if tw else 0
    labels = ["Low Risk", "Moderate Risk", "High Risk", "Very High Risk", "Critical Risk"]
    label  = labels[min(4, final // 20)]
    return final, label, bd


def render_pgx(outputs):
    score, label, bd = compute_pgx(outputs)
    SCORE_COLORS = ["#16A34A", "#D97706", "#EA580C", "#DC2626", "#B91C1C"]
    color = SCORE_COLORS[min(4, score // 20)]
    pills = ""
    for gene, _, ph, rl, _ in bd:
        rc = RISK_CFG.get(rl, RISK_CFG["Unknown"])
        pills += (f'<span class="pgx-pill" style="background:{rc["tag_bg"]};border-color:{rc["border"]};'
                  f'color:{rc["tag_text"]};">{gene} · {ph}</span>')
    st.markdown(f"""
    <div class="pgx-card">
      <div class="pgx-eyebrow">Polygenic Risk Score</div>
      <div class="pgx-score" style="color:{color};">{score}</div>
      <div class="pgx-label">{label} — composite across {len(outputs)} drug{"s" if len(outputs)!=1 else ""}</div>
      <div class="pgx-marker">
        <div class="pgx-fill" style="width:{score}%;background:linear-gradient(90deg,{color}99,{color});"></div>
        <div class="pgx-indicator" style="left:{score}%;border-color:{color};"></div>
      </div>
      <div class="pgx-thresh-labels">
        <span>0 — No Risk</span><span>25</span><span>50 — High</span><span>75</span><span>100 — Critical</span>
      </div>
      <div class="pgx-pills">{pills}</div>
    </div>""", unsafe_allow_html=True)


def render_risk_center(outputs, parsed):
    sev = max((o["risk_assessment"]["severity"] for o in outputs),
              key=lambda s: SEV_RANK.get(s, 0), default="none")
    sp  = SEV_CFG.get(sev, SEV_CFG["none"])
    EMO = {"none": "✓", "low": "△", "moderate": "⚠", "high": "⛔", "critical": "🚨"}
    hc  = sum(1 for o in outputs if o["risk_assessment"]["severity"] in ("high", "critical"))
    st.markdown(f"""
    <div class="risk-center" style="background:{sp['bg']};border-color:{sp['border']};color:{sp['text']};">
      <div class="rc-eyebrow">Risk Command Center</div>
      <div class="rc-headline">{EMO.get(sev,"")} {sp['label']} Risk Profile</div>
      <div class="rc-sub">Patient pharmacogenomic assessment across {len(outputs)} medication{"s" if len(outputs)!=1 else ""}</div>
      <div class="rc-stats" style="border-top-color:{sp['border']}88;">
        <div><div class="rc-stat-num">{len(outputs)}</div><div class="rc-stat-lbl">Drugs Assessed</div></div>
        <div><div class="rc-stat-num" style="{'color:#B91C1C' if hc else ''}">{hc}</div><div class="rc-stat-lbl">High / Critical</div></div>
        <div><div class="rc-stat-num">{len(parsed.get('detected_genes',[]))}</div><div class="rc-stat-lbl">Genes Detected</div></div>
        <div><div class="rc-stat-num">{parsed.get('total_variants',0)}</div><div class="rc-stat-lbl">Variants Found</div></div>
      </div>
    </div>""", unsafe_allow_html=True)


def render_critical_alerts(outputs):
    for o in outputs:
        if o["risk_assessment"]["severity"] == "critical":
            drug = o["drug"]
            note = o["clinical_recommendation"]["dosing_recommendation"][:240]
            st.markdown(f"""
            <div class="crit-alert">
              <div style="font-size:1.25rem;flex-shrink:0;padding-top:1px;">🚨</div>
              <div>
                <div class="crit-title">Critical Safety Alert — {drug}</div>
                <div class="crit-note">{note}{"…" if len(o["clinical_recommendation"]["dosing_recommendation"])>240 else ""}</div>
                <div class="crit-action">⚡ Contact prescribing physician immediately</div>
              </div>
            </div>""", unsafe_allow_html=True)


def render_disclaimer():
    st.markdown("""
    <div class="disclaimer-box">
      <span style="font-size:1rem;flex-shrink:0;">📋</span>
      <div class="disclaimer-text">
        <strong>For informational purposes only.</strong> These results require review by a qualified
        clinical pharmacologist or geneticist before any medication changes. All recommendations are
        based on CPIC Level A evidence — verify at <strong>cpicpgx.org</strong>.
      </div>
    </div>""", unsafe_allow_html=True)


def render_gene_row(outputs):
    GENE_ORDER = ["CYP2D6", "CYP2C19", "CYP2C9", "SLCO1B1", "TPMT", "DPYD"]
    gp = {o["pharmacogenomic_profile"]["primary_gene"]: o["pharmacogenomic_profile"]["phenotype"]
          for o in outputs}
    boxes = ""
    for g in GENE_ORDER:
        ph = gp.get(g, "Unknown")
        pc = PHENO_CFG.get(ph, PHENO_CFG["Unknown"])
        bar = min(100, pc["pct"])
        active = "active" if g in gp else ""
        boxes += f"""
        <div class="gene-box {active}" style="{'border-color:'+pc['border']+';' if active else ''}">
          <div class="gene-nm" style="{'color:'+pc['text']+';' if active else ''}">{g}</div>
          <div class="gene-track">
            <div class="gene-fill" style="width:{bar}%;background:{pc['bar']};"></div>
          </div>
          <div class="gene-ph" style="color:{pc['text'] if active else 'var(--text-xmuted)'};">{ph}</div>
        </div>"""
    sec("Gene Activity Overview")
    st.markdown(f'<div class="gene-row">{boxes}</div>', unsafe_allow_html=True)


def render_drug_table(outputs, pid):
    rows = ""
    data = []
    for o in outputs:
        drug = o["drug"]
        rl   = o["risk_assessment"]["risk_label"]
        sev  = o["risk_assessment"]["severity"]
        conf = o["risk_assessment"]["confidence_score"]
        gene = o["pharmacogenomic_profile"]["primary_gene"]
        ph   = o["pharmacogenomic_profile"]["phenotype"]
        rc   = RISK_CFG.get(rl, RISK_CFG["Unknown"])
        sp   = SEV_CFG.get(sev, SEV_CFG["none"])
        badge = risk_badge_html(rl)
        rows += f"""<div class="dtab-row">
          <div class="dtab-cell" style="font-weight:700;color:var(--text-primary);">{drug.title()}</div>
          <div class="dtab-cell">{badge}</div>
          <div class="dtab-cell"><span style="color:{sp['text']};font-weight:600;">{sp['label']}</span></div>
          <div class="dtab-cell" style="font-family:var(--font-mono);font-size:.8rem;color:var(--text-muted);">{gene}</div>
          <div class="dtab-cell"><span style="font-family:var(--font-mono);font-size:.8rem;color:{rc['tag_text']};background:{rc['tag_bg']};border:1px solid {rc['border']};padding:2px 8px;border-radius:4px;font-weight:600;">{ph}</span></div>
          <div class="dtab-cell">
            <div style="flex:1;height:4px;background:var(--surface-sub);border-radius:2px;overflow:hidden;margin-right:8px;">
              <div style="width:{conf*100:.0f}%;height:100%;background:{rc['severity_dot']};border-radius:2px;"></div>
            </div>
            <span style="font-family:var(--font-mono);font-size:.75rem;color:var(--text-muted);font-weight:600;">{conf:.0%}</span>
          </div>
        </div>"""
        data.append({"Drug": drug, "Risk": rl, "Severity": sev, "Gene": gene,
                      "Phenotype": ph, "Confidence": f"{conf:.0%}"})
    sec("Drug Risk Summary")
    st.markdown(f"""
    <div class="dtab">
      <div class="dtab-head">
        <div class="dtab-hcell">Drug</div><div class="dtab-hcell">Risk Label</div>
        <div class="dtab-hcell">Severity</div><div class="dtab-hcell">Gene</div>
        <div class="dtab-hcell">Phenotype</div><div class="dtab-hcell">Confidence</div>
      </div>{rows}
    </div>""", unsafe_allow_html=True)
    df = pd.DataFrame(data)
    st.download_button("⬇ Download CSV", data=df.to_csv(index=False),
        file_name=f"pharmaguard_{pid}.csv", mime="text/csv", key=f"csv_{pid}")


def render_heatmap(outputs):
    DRUG_ORD = ["CODEINE","WARFARIN","CLOPIDOGREL","SIMVASTATIN","AZATHIOPRINE","FLUOROURACIL"]
    GENE_ORD = ["CYP2D6","CYP2C9","CYP2C19","SLCO1B1","TPMT","DPYD"]
    DG = {"CODEINE":"CYP2D6","WARFARIN":"CYP2C9","CLOPIDOGREL":"CYP2C19",
          "SIMVASTATIN":"SLCO1B1","AZATHIOPRINE":"TPMT","FLUOROURACIL":"DPYD"}
    rmap  = {o["drug"]: o for o in outputs}
    drugs = [d for d in DRUG_ORD if d in rmap]
    if not drugs:
        return
    n = len(drugs)
    hdrs = '<div class="hm-header"></div>'
    for d in drugs:
        hdrs += f'<div class="hm-header">{d[:5]}</div>'
    rows = ""
    for gene in GENE_ORD:
        rows += f'<div class="hm-header" style="justify-content:flex-end;padding-right:6px;">{gene}</div>'
        for d in drugs:
            if DG.get(d) == gene and d in rmap:
                o  = rmap[d]
                rl = o["risk_assessment"]["risk_label"]
                ph = o["pharmacogenomic_profile"]["phenotype"]
                rc = RISK_CFG.get(rl, RISK_CFG["Unknown"])
                sh = {"Adjust Dosage":"Adjust","Ineffective":"Ineffect.","Unknown":"?"}.get(rl, rl)
                rows += (f'<div class="hm-cell" style="background:{rc["bg"]};border-color:{rc["border"]};" '
                         f'title="{d}×{gene}: {rl} ({ph})">'
                         f'<div class="hm-cell-name" style="color:{rc["text"]};">{sh}</div>'
                         f'<div class="hm-cell-risk" style="color:{rc["text"]};">{ph}</div></div>')
            else:
                rows += '<div class="hm-cell" style="background:var(--surface-sub);border-color:var(--border-light);"><div class="hm-cell-risk" style="color:var(--text-xmuted);">—</div></div>'
    legend = "".join(
        f'<div class="hm-legend-item"><span class="hm-dot" style="background:{RISK_CFG[r]["bg"]};border-color:{RISK_CFG[r]["border"]};"></span><span>{RISK_CFG[r]["shape"]} {r}</span></div>'
        for r in ["Safe", "Adjust Dosage", "Toxic", "Ineffective"])
    st.markdown(f"""
    <div class="hm-wrap">
      <div class="hm-eyebrow">Drug × Gene Risk Matrix</div>
      <div class="hm-grid" style="grid-template-columns:80px repeat({n},1fr);">{hdrs}{rows}</div>
      <div class="hm-legend">{legend}</div>
    </div>""", unsafe_allow_html=True)


def render_chromosome(outputs, parsed):
    det  = set(parsed.get("detected_genes", []))
    rmap = {o["pharmacogenomic_profile"]["primary_gene"]: o for o in outputs}
    rows = ""
    for gene, info in CHROM_INFO.items():
        ch  = info["chrom"]
        pos = info["pos_mb"]
        pct = (pos / CHROM_LEN.get(ch, 200)) * 100
        if gene in rmap:
            rl = rmap[gene]["risk_assessment"]["risk_label"]
            mc = RISK_CFG.get(rl, RISK_CFG["Unknown"])["severity_dot"]
        elif gene in det:
            mc = "#94A3B8"
        else:
            mc = "#DDE3EE"
        rows += f"""<div class="chrom-row">
          <div class="chrom-chr">{ch}</div>
          <div class="chrom-bar">
            <div class="chrom-body"></div>
            <div class="chrom-marker" style="left:{pct}%;background:{mc};box-shadow:0 0 5px {mc}88;"></div>
          </div>
          <div class="chrom-gene">{gene}</div>
          <div class="chrom-band">{info['band']}</div>
        </div>"""
    st.markdown(f"""
    <div class="chrom-wrap">
      <div class="chrom-eyebrow">Variant Chromosome Locations</div>
      {rows}
      <div style="font-family:var(--font-mono);font-size:.65rem;color:var(--text-xmuted);margin-top:var(--sp-3);">
        Coloured markers = variants detected · Grey = undetected
      </div>
    </div>""", unsafe_allow_html=True)


def render_pop_freq(gene, ph):
    freq = POP_FREQ.get(gene, {})
    if not freq:
        return
    rows = ""
    for p, pct in sorted(freq.items(), key=lambda x: -x[1]):
        you = (p == ph)
        pc  = PHENO_CFG.get(p, PHENO_CFG["Unknown"])
        you_tag = f'<span class="pop-you">← You</span>' if you else ""
        w = "font-weight:700;" if you else ""
        rows += f"""<div class="pop-row">
          <div class="pop-ph" style="{w}{'color:'+pc['text']+';' if you else ''}">{pc['label']}</div>
          <div class="pop-track"><div class="pop-fill" style="width:{min(pct,100)}%;background:{pc['bar'] if you else '#CBD5E1'};"></div></div>
          <div class="pop-pct" style="{w}{'color:'+pc['text']+';' if you else ''}">{pct}%{you_tag}</div>
        </div>"""
    st.markdown(f"""
    <div class="pop-wrap">
      <div class="pop-eyebrow">{gene} — Population Distribution</div>{rows}
    </div>""", unsafe_allow_html=True)


def render_ix_matrix(outputs, ix):
    if not ix or len(outputs) < 2:
        return
    drugs = [o["drug"] for o in outputs]
    n     = len(drugs)
    sm    = {}
    for x in ix.get("all_interactions", []):
        inv = x.get("drugs_involved", [])
        if len(inv) == 2:
            sv = x.get("severity", "none")
            sm[(inv[0], inv[1])] = sm[(inv[1], inv[0])] = sv
    MC = {
        "critical": {"bg":"#FEF2F2","text":"#7F1D1D","border":"#FECACA"},
        "high":     {"bg":"#FEF2F2","text":"#7F1D1D","border":"#FECACA"},
        "moderate": {"bg":"#FFFBEB","text":"#78350F","border":"#FDE68A"},
        "low":      {"bg":"#FEFCE8","text":"#713F12","border":"#FDE047"},
        "none":     {"bg":"#F0FDF4","text":"#14532D","border":"#BBF7D0"},
        "diag":     {"bg":"var(--surface-sub)","text":"var(--text-muted)","border":"var(--border-light)"},
    }
    hdrs = '<div class="ix-head"></div>'
    for d in drugs:
        hdrs += f'<div class="ix-head">{d[:6]}</div>'
    grid = ""
    for i, d1 in enumerate(drugs):
        grid += f'<div class="ix-head" style="justify-content:flex-end;padding-right:4px;">{d1[:6]}</div>'
        for j, d2 in enumerate(drugs):
            if i == j:
                mc = MC["diag"]
                grid += f'<div class="ix-cell" style="background:{mc["bg"]};border-color:{mc["border"]};color:{mc["text"]};">—</div>'
            else:
                sv = sm.get((d1, d2), "none")
                mc = MC.get(sv, MC["none"])
                lbl = sv.upper() if sv != "none" else "OK"
                grid += f'<div class="ix-cell" style="background:{mc["bg"]};border-color:{mc["border"]};color:{mc["text"]};">{lbl}</div>'
    sec("Drug Interaction Matrix")
    st.markdown(f"""
    <div style="background:var(--surface);border:1px solid var(--border-light);border-radius:var(--r-xl);padding:var(--sp-5);margin-bottom:var(--sp-4);box-shadow:var(--shadow-sm);">
      <div class="ix-grid" style="grid-template-columns:76px repeat({n},1fr);gap:3px;">{hdrs}{grid}</div>
    </div>""", unsafe_allow_html=True)
    shown = set()
    for x in ix.get("all_interactions", []):
        inv = x.get("drugs_involved", [])
        key = tuple(sorted(inv))
        if len(inv) == 2 and key not in shown:
            shown.add(key)
            sv = x.get("severity", "low")
            sp = SEV_CFG.get(sv, SEV_CFG["low"])
            with st.expander(f"{' + '.join(inv)}  —  {sv.upper()} interaction"):
                mech = x.get("mechanism", x.get("message", ""))
                rec  = x.get("recommendation", "")
                if mech:
                    st.markdown(f'<div style="font-size:.9rem;color:var(--text-secondary);line-height:1.75;margin-bottom:var(--sp-2);">{mech}</div>', unsafe_allow_html=True)
                if rec:
                    st.markdown(f'<div style="font-family:var(--font-mono);font-size:.8rem;color:{sp["text"]};margin-top:var(--sp-2);font-weight:600;">→ {rec}</div>', unsafe_allow_html=True)


def render_narrative(outputs, parsed, pid, key, skip_llm):
    results_for = [{"drug": o["drug"],
                    "primary_gene": o["pharmacogenomic_profile"]["primary_gene"],
                    "phenotype": o["pharmacogenomic_profile"]["phenotype"],
                    "risk_label": o["risk_assessment"]["risk_label"],
                    "severity": o["risk_assessment"]["severity"]}
                   for o in outputs]
    with st.spinner("Generating AI clinical summary…"):
        nar = generate_patient_narrative(pid, results_for, parsed, key, skip_llm)
    model_label = "Static Template" if (skip_llm or not key) else "LLaMA 3.3 70B"
    sec("AI Clinical Summary")
    st.markdown(f"""
    <div class="narrative-box">
      <div class="narrative-header">
        <span class="ai-badge-pill">{model_label}</span>
        <span style="font-size:.9rem;font-weight:600;color:var(--brand-dark);">Unified Patient Summary</span>
      </div>
      <div class="narrative-text">{nar}</div>
    </div>""", unsafe_allow_html=True)


def render_before_after(outputs):
    bad = [o for o in outputs if o["risk_assessment"]["risk_label"] in ("Toxic", "Ineffective")]
    if not bad:
        return
    o    = bad[0]
    drug = o["drug"]
    rl   = o["risk_assessment"]["risk_label"]
    alts = o["clinical_recommendation"].get("alternative_drugs", [])
    alt  = alts[0] if alts else "Alternative medication"
    gene = o["pharmacogenomic_profile"]["primary_gene"]
    ph   = o["pharmacogenomic_profile"]["phenotype"]
    BEFORE = {
        "Toxic":      f"Standard {drug.lower()} dose → toxic accumulation → serious harm",
        "Ineffective": f"Standard {drug.lower()} dose → zero therapeutic effect → treatment failure",
    }
    sec("Clinical Impact — Before & After PGx")
    st.markdown(f"""
    <div class="ba-grid">
      <div class="ba-side" style="background:#FFF1F2;border-right:1px solid var(--border-light);">
        <div class="ba-side-lbl" style="color:#B91C1C;">⛔ Without PharmaGuard</div>
        <div class="ba-drug" style="color:#7F1D1D;">{drug.title()} — Standard Protocol</div>
        <div class="ba-text" style="color:#7F1D1D;">{BEFORE.get(rl,"Risk undetected")}</div>
        <div class="ba-gene" style="color:#FECACA;">{gene} {ph} phenotype undetected</div>
      </div>
      <div class="ba-side" style="background:#F0FDF4;">
        <div class="ba-side-lbl" style="color:#14532D;">✓ With PharmaGuard</div>
        <div class="ba-drug" style="color:#15803D;">{alt} — PGx-Guided Protocol</div>
        <div class="ba-text" style="color:#16A34A;">Appropriate alternative selected → safe, effective therapy</div>
        <div class="ba-gene" style="color:#BBF7D0;">{gene} {ph} phenotype identified → therapy optimised</div>
      </div>
    </div>""", unsafe_allow_html=True)


def render_rx_checker(outputs):
    sec("Prescription Safety Checker")
    rmap  = {o["drug"]: o for o in outputs}
    drugs = [o["drug"] for o in outputs]
    c1, c2 = st.columns([2, 1])
    with c1:
        sel = st.selectbox("Select drug", drugs,
              format_func=lambda x: f"{x.title()}  ({GENE_DRUG_MAP.get(x,'')})",
              key="rx_drug", label_visibility="collapsed")
    with c2:
        check = st.button("Check Safety →", key="rx_check")
    if check and sel in rmap:
        o    = rmap[sel]
        rl   = o["risk_assessment"]["risk_label"]
        sev  = o["risk_assessment"]["severity"]
        rec  = o["clinical_recommendation"]["dosing_recommendation"]
        gene = o["pharmacogenomic_profile"]["primary_gene"]
        ph   = o["pharmacogenomic_profile"]["phenotype"]
        rc   = RISK_CFG.get(rl, RISK_CFG["Unknown"])
        sp   = SEV_CFG.get(sev, SEV_CFG["none"])
        VERDICT = {
            "Safe":         "✓ Safe to Prescribe",
            "Adjust Dosage":"△ Prescribe with Dose Adjustment",
            "Toxic":        "⛔ Do Not Prescribe — Toxicity Risk",
            "Ineffective":  "◆ Do Not Prescribe — Drug Ineffective",
        }
        st.markdown(f"""
        <div class="rx-result" style="background:{rc['bg']};border-color:{rc['border']};">
          <div class="rx-verdict" style="color:{rc['text']};">{VERDICT.get(rl, rl)}</div>
          <div class="rx-detail">{gene} {ph} phenotype detected. {rec}</div>
          <div class="rx-meta" style="color:{sp['text']};">Severity: {sp['label']} · Confidence: {o["risk_assessment"]["confidence_score"]:.0%} · CPIC Level A</div>
        </div>""", unsafe_allow_html=True)
    elif not check:
        st.markdown("""<div class="info-strip"><span>🔍</span>
          <div class="info-strip-text">Select a drug and click <strong>Check Safety</strong> to validate against this patient's genotype.</div>
        </div>""", unsafe_allow_html=True)


def render_clinical_note(outputs, pid):
    lines = [f"PharmaGuard Clinical Note — Patient {pid} — {datetime.utcnow().strftime('%Y-%m-%d')}",
             "=" * 60, ""]
    for o in outputs:
        gene = o["pharmacogenomic_profile"]["primary_gene"]
        dip  = o["pharmacogenomic_profile"]["diplotype"]
        ph   = o["pharmacogenomic_profile"]["phenotype"]
        drug = o["drug"]
        rl   = o["risk_assessment"]["risk_label"]
        rec  = o["clinical_recommendation"]["dosing_recommendation"]
        alts = o["clinical_recommendation"].get("alternative_drugs", [])
        lines.append(f"DRUG: {drug}")
        lines.append(f"Gene: {gene} | Diplotype: {dip} | Phenotype: {ph} | Risk: {rl}")
        lines.append(f"CPIC: {rec}")
        if alts:
            lines.append(f"Alternatives: {', '.join(alts)}")
        lines.append("")
    lines += ["", "—" * 60,
              "Generated by PharmaGuard v9.1 · CPIC Level A evidence · cpicpgx.org",
              "NOT FOR CLINICAL USE WITHOUT VALIDATION BY A QUALIFIED CLINICIAN"]
    note = "\n".join(lines)
    sec("One-Click Clinical Note")
    st.markdown(f'<div class="note-box"><pre>{note}</pre></div>', unsafe_allow_html=True)
    st.download_button("⬇ Download Clinical Note", data=note,
        file_name=f"clinical_note_{pid}.txt", mime="text/plain", key=f"note_{pid}")


def render_patient_mode(outputs):
    bad = any(o["risk_assessment"]["risk_label"] in ("Toxic","Ineffective") for o in outputs)
    if bad:
        st.markdown("""<div class="patient-banner" style="background:#FFF1F2;border-color:#FECACA;">
          <div class="patient-banner-title" style="color:#B91C1C;">🚨 Important — Some medications need urgent attention</div>
          <div class="patient-banner-sub" style="color:#7F1D1D;">Your genetic results show that one or more medications may not be safe or effective for you. Please speak with your doctor before taking these medications.</div>
        </div>""", unsafe_allow_html=True)
    else:
        st.markdown("""<div class="patient-banner" style="background:#F0FDF4;border-color:#BBF7D0;">
          <div class="patient-banner-title" style="color:#14532D;">✓ Good news — Your medications look safe</div>
          <div class="patient-banner-sub" style="color:#16A34A;">Based on your genetic profile, the medications reviewed are predicted to work normally at standard doses.</div>
        </div>""", unsafe_allow_html=True)
    for o in outputs:
        drug    = o["drug"]
        rl      = o["risk_assessment"]["risk_label"]
        gene    = o["pharmacogenomic_profile"]["primary_gene"]
        ph      = o["pharmacogenomic_profile"]["phenotype"]
        alts    = o["clinical_recommendation"].get("alternative_drugs", [])
        phplain = PLAIN_PHENO.get(ph, ph)
        explain = PLAIN_RISK.get((drug, ph), "")
        VERDICT = {
            "Safe":         "✓ This medicine is likely safe for you",
            "Adjust Dosage":"△ You may need a different dose",
            "Toxic":        "⛔ This medicine could be harmful to you",
            "Ineffective":  "◆ This medicine likely won't work for you",
        }
        rc = RISK_CFG.get(rl, RISK_CFG["Unknown"])
        action = ""
        if rl in ("Toxic", "Ineffective"):
            alt_text = f"They may suggest: <strong>{', '.join(alts[:3])}</strong>" if alts else "Ask about alternative medications."
            action = f'<div class="pcard-action"><span style="font-size:1rem;">💊</span><div class="pcard-action-text"><strong>Talk to your doctor before taking {drug.title()}.</strong><br>{alt_text}</div></div>'
        elif rl == "Adjust Dosage":
            action = f'<div class="pcard-action"><span style="font-size:1rem;">📋</span><div class="pcard-action-text"><strong>Tell your doctor about this result before starting {drug.title()}.</strong><br>You may need a different dose than usually prescribed.</div></div>'
        st.markdown(f"""
        <div class="pcard" style="border-color:{rc['border']};">
          <div class="pcard-drug">{drug.title()}</div>
          <div class="pcard-verdict" style="color:{rc['text']};">{VERDICT.get(rl, rl)}</div>
          <div class="pcard-gene">{gene} · {phplain}</div>
          {f'<div class="pcard-plain">{explain}</div>' if explain else ''}
          {action}
        </div>""", unsafe_allow_html=True)


# ══════════════════════════════════════════════════════════════════════════════
# MASTER RESULTS RENDERER
# ══════════════════════════════════════════════════════════════════════════════

def render_results(outputs, parsed, ix, pdf_bytes, pid, patient_mode=False, key="", skip_llm=False):
    render_disclaimer()
    render_risk_center(outputs, parsed)
    render_critical_alerts(outputs)

    dc1, dc2, dc3 = st.columns(3)
    with dc1:
        st.download_button("⬇ Download All JSON", data=json.dumps(outputs, indent=2),
            file_name=f"pharmaguard_{pid}.json", mime="application/json",
            use_container_width=True, key=f"dlall_{pid}")
    with dc2:
        if pdf_bytes:
            st.download_button("⬇ Download PDF Report", data=pdf_bytes,
                file_name=f"pharmaguard_{pid}.pdf", mime="application/pdf",
                use_container_width=True, key=f"dlpdf_{pid}")
    with dc3:
        if ix and ix.get("interactions_found"):
            st.download_button("⬇ Interactions JSON", data=json.dumps(ix, indent=2),
                file_name=f"pharmaguard_{pid}_ix.json", mime="application/json",
                use_container_width=True, key=f"dlix_{pid}")

    st.markdown("<div style='height:var(--sp-3)'></div>", unsafe_allow_html=True)

    if patient_mode:
        render_patient_mode(outputs)
        return

    render_gene_row(outputs)
    render_drug_table(outputs, pid)
    render_pgx(outputs)

    c1, c2 = st.columns([1.4, 1], gap="large")
    with c1: render_heatmap(outputs)
    with c2: render_chromosome(outputs, parsed)

    if ix and len(outputs) >= 2:
        render_ix_matrix(outputs, ix)

    render_narrative(outputs, parsed, pid, key, skip_llm)
    render_before_after(outputs)
    render_rx_checker(outputs)
    render_clinical_note(outputs, pid)

    # ── Individual Drug Cards ─────────────────────────────────────────────────
    sec("Individual Drug Analysis")
    for output in outputs:
        rl   = output["risk_assessment"]["risk_label"]
        drug = output["drug"]
        sev  = output["risk_assessment"]["severity"]
        conf = output["risk_assessment"]["confidence_score"]
        gene = output["pharmacogenomic_profile"]["primary_gene"]
        dip  = output["pharmacogenomic_profile"]["diplotype"]
        ph   = output["pharmacogenomic_profile"]["phenotype"]
        var  = output["pharmacogenomic_profile"]["detected_variants"]
        rec  = output["clinical_recommendation"]["dosing_recommendation"]
        alts = output["clinical_recommendation"].get("alternative_drugs", [])
        mon  = output["clinical_recommendation"].get("monitoring_required", "")
        exp  = output["llm_generated_explanation"]
        rc   = RISK_CFG.get(rl, RISK_CFG["Unknown"])
        sp   = SEV_CFG.get(sev, SEV_CFG["none"])
        cpic_lv = output.get("pharmacogenomic_profile", {}).get("cpic_evidence_level", "Level A")

        st.markdown(f"""
        <div class="dcard reveal-card">
          <div class="dcard-header">
            <div class="dcard-left">
              <div class="dcard-indicator" style="background:{rc['severity_dot']};box-shadow:0 0 0 3px {rc['bg']};"></div>
              <div>
                <div class="dcard-drug">{drug.title()}
                  <span class="cpic-badge">CPIC {cpic_lv}</span>
                </div>
                <div class="dcard-meta">{gene} · {dip} · {ph}</div>
              </div>
            </div>
            {risk_badge_html(rl)}
          </div>
          <div class="dcard-body">
            <div class="metrics-row">
              <div class="metric-cell"><div class="metric-key">Phenotype</div><div class="metric-val" style="color:{rc['text']};font-size:.95rem;">{ph}</div></div>
              <div class="metric-cell"><div class="metric-key">Severity</div><div class="metric-val" style="color:{sp['text']};font-size:.95rem;">{sp['label']}</div></div>
              <div class="metric-cell"><div class="metric-key">Confidence</div><div class="metric-val">{conf:.0%}</div></div>
              <div class="metric-cell"><div class="metric-key">Variants</div><div class="metric-val">{len(var)}</div></div>
            </div>""", unsafe_allow_html=True)

        dq = min(1.0, len(var) / 3.0)
        st.markdown(f"""
        <div class="conf-grid">
          <div>
            <div class="conf-label"><span>Prediction Confidence</span><span style="color:{rc['severity_dot']};font-weight:700;">{conf:.0%}</span></div>
            <div class="conf-track"><div class="conf-fill" style="width:{conf*100:.1f}%;background:{rc['severity_dot']};"></div></div>
          </div>
          <div>
            <div class="conf-label"><span>Data Quality</span><span style="color:var(--text-muted);">{len(var)} variant{"s" if len(var)!=1 else ""}</span></div>
            <div class="conf-track"><div class="conf-fill" style="width:{dq*100:.1f}%;background:#94A3B8;"></div></div>
          </div>
        </div>""", unsafe_allow_html=True)

        if var:
            rows_html = ""
            for v in var:
                fc = func_cls(v.get("functional_status", ""))
                fn = (v.get("functional_status") or "unknown").replace("_", " ").title()
                rows_html += (f'<tr><td class="v-rsid">{v.get("rsid","—")}</td>'
                              f'<td class="v-star">{v.get("star_allele","—")}</td>'
                              f'<td class="{fc}">{fn}</td></tr>')
            st.markdown(f"""
            <div style="margin-bottom:var(--sp-4);">
              <div class="conf-label" style="margin-bottom:var(--sp-3);">Detected Variants ({len(var)})</div>
              <table class="vtable">
                <thead><tr><th>rsID</th><th>Star Allele</th><th>Functional Status</th></tr></thead>
                <tbody>{rows_html}</tbody>
              </table>
            </div>""", unsafe_allow_html=True)

        # CPIC Recommendation
        rec_color = rc["bg"]; rec_border = rc["border"]; rec_text_c = rc["text"]
        st.markdown(f"""
        <div class="rec-box" style="background:{rec_color};border-color:{rec_border};">
          <div class="rec-label" style="color:{rec_text_c};">CPIC Recommendation — {drug}</div>
          <div class="rec-text">{rec}</div>
        </div>""", unsafe_allow_html=True)

        if mon:
            st.markdown(f"""
            <div class="rec-box" style="background:var(--surface-sub);border-color:var(--border-light);">
              <div class="rec-label" style="color:var(--text-muted);">🔬 Monitoring Protocol</div>
              <div class="rec-text">{mon}</div>
            </div>""", unsafe_allow_html=True)

        if alts:
            chips = "".join(f'<span class="alt-chip">{a}</span>' for a in alts)
            st.markdown(f"""
            <div style="margin-bottom:var(--sp-4);">
              <div class="conf-label" style="margin-bottom:var(--sp-2);">Alternative Medications</div>
              <div class="alt-chips">{chips}</div>
            </div>""", unsafe_allow_html=True)

        # Population frequency
        render_pop_freq(gene, ph)

        # ── AI Explanation ────────────────────────────────────────────────────
        if exp.get("summary"):
            raw_model = exp.get("model_used", "llama-3.3-70b")
            # BUG 1 FIX: clean_model_label() strips "(rate-limited)" etc.
            model, is_static = clean_model_label(raw_model)
            blocks = ""
            for lbl, k in [("Summary","summary"), ("Biological Mechanism","biological_mechanism"),
                           ("Variant Significance","variant_significance"), ("Clinical Implications","clinical_implications")]:
                if exp.get(k):
                    blocks += (f'<div class="ai-section">'
                               f'<div class="ai-sec-label">{lbl}</div>'
                               f'<div class="ai-sec-text">{exp[k]}</div>'
                               f'</div>')
            st.markdown(f"""
            <div class="ai-block">
              <div class="ai-header">
                <span class="ai-badge-pill">{model}</span>
                <span class="ai-title">AI Explanation · {drug}</span>
              </div>{blocks}
            </div>""", unsafe_allow_html=True)

        with st.expander(f"Raw JSON — {drug}"):
            st.json(output)

        st.markdown('</div></div>', unsafe_allow_html=True)


# ══════════════════════════════════════════════════════════════════════════════
# NAVIGATION + LAYOUT
# ══════════════════════════════════════════════════════════════════════════════

def render_nav(key_present):
    model_badge = "LLaMA 3.3 70B" if key_present else "Static Mode"
    model_cls   = "pg-badge-brand" if key_present else "pg-badge-default"
    st.markdown(f"""
    <div class="pg-nav">
      <div class="pg-brand">
        <div class="pg-brand-name">Pharma<span>Guard</span></div>
        <div class="pg-brand-sub">v9.1 · Precision Clinical · RIFT 2026</div>
      </div>
      <div class="pg-nav-badges">
        <span class="pg-badge pg-badge-brand">CPIC Level A</span>
        <span class="pg-badge {model_cls}">{model_badge}</span>
        <span class="pg-badge pg-badge-default">RIFT 2026</span>
      </div>
    </div>
    <div class="trust-strip">
      <div class="trust-item">🔒 Genetic data analyzed locally — never stored</div>
      <div class="trust-sep"></div>
      <div class="trust-item">✓ CPIC Level A Evidence Guidelines</div>
      <div class="trust-sep"></div>
      <div class="trust-item">🧬 6 Pharmacogenes · 6 High-Risk Drugs</div>
      <div class="trust-sep"></div>
      <div class="trust-item">⚕ For review by qualified clinicians only</div>
    </div>""", unsafe_allow_html=True)


def render_steps(has_vcf, has_drugs, has_results):
    steps = [
        ("01", "Upload VCF", has_vcf),
        ("02", "Select Drugs", has_drugs),
        ("03", "Run Analysis", has_results),
        ("04", "Review Results", has_results),
    ]
    html = '<div class="steps">'
    for num, lbl, done in steps:
        cls = "step done" if done else "step"
        html += f'<div class="{cls}"><div class="step-num">{num}</div><div class="step-lbl">{lbl}</div></div>'
    html += "</div>"
    st.markdown(html, unsafe_allow_html=True)


def render_persona_demo(key):
    SEV_COLORS = {
        "critical": {"sev_bg":"#FEF2F2","sev_border":"#FECACA","sev_text":"#7F1D1D","sev_label":"Critical"},
        "high":     {"sev_bg":"#FFF7ED","sev_border":"#FED7AA","sev_text":"#7C2D12","sev_label":"High Risk"},
        "moderate": {"sev_bg":"#FFFBEB","sev_border":"#FDE68A","sev_text":"#78350F","sev_label":"Moderate"},
        "none":     {"sev_bg":"#F0FDF4","sev_border":"#BBF7D0","sev_text":"#14532D","sev_label":"All Safe"},
    }
    st.markdown('<div style="font-size:.8rem;font-weight:600;letter-spacing:.1em;text-transform:uppercase;color:var(--text-muted);margin-bottom:var(--sp-3);">Quick Demo — Select Patient Persona</div>', unsafe_allow_html=True)
    cols = st.columns(4)
    for i, (pid, p) in enumerate(PERSONAS.items()):
        sc = SEV_COLORS.get(p["sev"], SEV_COLORS["none"])
        with cols[i]:
            st.markdown(f"""
            <div class="persona-card">
              <div class="pc-sev" style="background:{sc['sev_bg']};border-color:{sc['sev_border']};color:{sc['sev_text']};">
                {sc['sev_label']}
              </div>
              <div class="pc-name">{p['label']}</div>
              <div class="pc-desc">{p['desc']}</div>
            </div>""", unsafe_allow_html=True)
            if st.button(f"Load {p['label']}", key=f"persona_{pid}", use_container_width=True):
                vcf = load_vcf(p["file"])
                pid_gen = f"PG-{uuid.uuid4().hex[:8].upper()}"
                with st.spinner(f"Running {p['label']} analysis…"):
                    parsed, results, outputs, ix, pdf = run_pipeline(
                        vcf, p["drugs"], pid_gen, key, skip_llm=not bool(key))
                st.session_state["results"]      = outputs
                st.session_state["parsed"]       = parsed
                st.session_state["ix"]           = ix
                st.session_state["pdf"]          = pdf
                st.session_state["patient_id"]   = pid_gen
                st.session_state["results_key"]  = key
                st.session_state["results_skip"] = not bool(key)
                st.rerun()


def render_test_suite(key):
    st.markdown("### Test Suite")
    for i, tc in enumerate(TEST_SUITE):
        with st.container():
            st.markdown(f"""
            <div class="tc-card">
              <div class="tc-name">{tc['name']}</div>
              <div class="tc-desc">{tc['desc']}</div>
            </div>""", unsafe_allow_html=True)
            if st.button(f"▶ Run: {tc['name']}", key=f"tc_{i}", use_container_width=True):
                vcf = load_vcf(tc["file"])
                pid = f"TC-{uuid.uuid4().hex[:6].upper()}"
                with st.spinner(f"Running {tc['name']}…"):
                    parsed, results, outputs, ix, pdf = run_pipeline(
                        vcf, tc["drugs"], pid, key, skip_llm=True)
                # Show pass/fail
                passed = all(o["risk_assessment"]["risk_label"] == tc["expected"].get(o["drug"])
                             for o in outputs if o["drug"] in tc["expected"])
                if passed:
                    st.success(f"✓ PASS — {tc['name']}")
                else:
                    st.error(f"✗ FAIL — {tc['name']}")
                    for o in outputs:
                        drug = o["drug"]
                        got  = o["risk_assessment"]["risk_label"]
                        want = tc["expected"].get(drug, "—")
                        icon = "✓" if got == want else "✗"
                        st.markdown(f"`{icon} {drug}: got={got} want={want}`")
                st.session_state["results"]      = outputs
                st.session_state["parsed"]       = parsed
                st.session_state["ix"]           = ix
                st.session_state["pdf"]          = pdf
                st.session_state["patient_id"]   = pid
                st.session_state["results_key"]  = key
                st.session_state["results_skip"] = True
                st.rerun()


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════

def main():
    # ── Sidebar ───────────────────────────────────────────────────────────────
    with st.sidebar:
        st.markdown("### ⚙ Settings")
        groq_key  = st.text_input("Groq API Key", type="password",
                        placeholder="gsk_…", help="Required for LLM explanations")
        model_sel = st.selectbox("Model", ["LLaMA 3.3 70B Versatile"],
                        help="Model used for AI explanations")
        use_static = st.checkbox("Test mode: instant (no API call)", value=not bool(groq_key))
        st.markdown("---")
        st.markdown("**Gene → Drug Map**")
        for drug, gene in GENE_DRUG_MAP.items():
            st.markdown(f"`{gene}` → {drug.title()}")

    key       = groq_key.strip() if groq_key else ""
    skip_llm  = use_static or not key

    render_nav(bool(key))

    # ── Tabs ──────────────────────────────────────────────────────────────────
    tab_analysis, tab_suite = st.tabs(["Analysis", "Test Suite"])

    # ── Analysis Tab ──────────────────────────────────────────────────────────
    with tab_analysis:
        has_results = "results" in st.session_state and st.session_state["results"]
        render_steps(has_vcf=True, has_drugs=True, has_results=has_results)
        render_persona_demo(key)

        col_input, col_results = st.columns([1, 2], gap="large")

        with col_input:
            sec("Genomic Data")

            # Persona shortcut info
            persona_sel = st.selectbox("Or select a test scenario",
                options=["None"] + [p["label"] for p in PERSONAS.values()],
                key="persona_sel")

            vcf_file = st.file_uploader(
                "Upload VCF file", type=["vcf"],
                help="Drag and drop file here\nLimit 200MB per file • VCF",
                label_visibility="collapsed")

            # Determine VCF source
            vcf_text = None
            if vcf_file:
                vcf_text = vcf_file.read().decode("utf-8")
            elif persona_sel != "None":
                for p in PERSONAS.values():
                    if p["label"] == persona_sel:
                        vcf_text = load_vcf(p["file"])
                        break

            if vcf_text:
                st.success(f"✓ {getattr(vcf_file,'name','Scenario VCF')} · {len(vcf_text)/1024:.2f} MB — Ready for analysis")

            sec("Medications to Analyse")
            default_drugs = ALL_DRUGS
            if persona_sel != "None":
                for p in PERSONAS.values():
                    if p["label"] == persona_sel:
                        default_drugs = p["drugs"]
                        break
            selected_drugs = st.multiselect("Select drugs", ALL_DRUGS,
                default=default_drugs, label_visibility="collapsed")
            custom_raw = st.text_input("Custom drugs (comma-separated)", placeholder="CODEINE, WARFARIN…")
            if custom_raw:
                extras = [d.strip().upper() for d in custom_raw.split(",") if d.strip()]
                selected_drugs = list(dict.fromkeys(selected_drugs + extras))

            sec("Patient ID")
            patient_id_input = st.text_input("Patient ID", placeholder="Auto-generated if blank",
                                              label_visibility="collapsed")
            pid = patient_id_input.strip() or f"PG-{uuid.uuid4().hex[:8].upper()}"

            sec("View Mode")
            patient_mode = st.checkbox("Patient-friendly view (plain language)", value=False)

            run_btn = st.button("Run Analysis →", use_container_width=True,
                                disabled=not vcf_text or not selected_drugs)

            if run_btn and vcf_text and selected_drugs:
                with st.spinner("Analysing pharmacogenomic profile…"):
                    parsed, results, outputs, ix, pdf = run_pipeline(
                        vcf_text, selected_drugs, pid, key,
                        run_ix=len(selected_drugs) > 1,
                        gen_pdf=True,
                        skip_llm=skip_llm)
                st.session_state["results"]      = outputs
                st.session_state["parsed"]       = parsed
                st.session_state["ix"]           = ix
                st.session_state["pdf"]          = pdf
                st.session_state["patient_id"]   = pid
                st.session_state["results_key"]  = key
                st.session_state["results_skip"] = skip_llm
                st.rerun()

        with col_results:
            sec("Analysis Results")
            if "results" in st.session_state and st.session_state["results"]:
                res_pid  = st.session_state["patient_id"]
                res_outs = st.session_state["results"]
                res_par  = st.session_state["parsed"]
                res_ix   = st.session_state.get("ix")
                res_pdf  = st.session_state.get("pdf")
                res_key  = st.session_state.get("results_key", key)
                res_skip = st.session_state.get("results_skip", skip_llm)
                st.markdown(f"""
                <div style="display:flex;align-items:center;gap:var(--sp-3);margin-bottom:var(--sp-4);">
                  <span style="font-family:var(--font-mono);font-size:1rem;font-weight:700;
                    color:var(--brand);background:var(--brand-light);border:1px solid var(--brand-border);
                    padding:4px 12px;border-radius:var(--r-full);">{res_pid}</span>
                  <span style="font-size:.75rem;color:var(--text-muted);cursor:pointer;" title="Copy">📋</span>
                </div>""", unsafe_allow_html=True)
                render_results(res_outs, res_par, res_ix, res_pdf, res_pid,
                               patient_mode=patient_mode, key=res_key, skip_llm=res_skip)
            else:
                st.markdown("""
                <div class="empty-state">
                  <span class="empty-icon">🧬</span>
                  <div class="empty-title">No analysis results yet</div>
                  <div class="empty-hint">
                    Upload a VCF file<br>
                    Select medications<br>
                    Click Run Analysis
                  </div>
                </div>""", unsafe_allow_html=True)

    # ── Test Suite Tab ────────────────────────────────────────────────────────
    with tab_suite:
        render_test_suite(key)


if __name__ == "__main__":
    main()