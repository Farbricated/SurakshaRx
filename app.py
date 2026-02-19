"""
PharmaGuard v1.0 — Pharmacogenomic Risk Prediction System
RIFT 2026 Hackathon | Pharmacogenomics / Explainable AI Track

Streamlit Application Entry Point
"""

import streamlit as st
import json
import uuid
import os
from datetime import datetime

from vcf_parser import parse_vcf, get_sample_vcf
from risk_engine import run_risk_assessment, get_overall_severity, DRUG_RISK_TABLE
from llm_explainer import generate_all_explanations
from schema import build_output_schema

# ─── Page Config ──────────────────────────────────────────────
st.set_page_config(
    page_title="PharmaGuard — Pharmacogenomic Risk System",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ─── Styling ──────────────────────────────────────────────────
st.markdown("""
<style>
    @import url('https://fonts.googleapis.com/css2?family=Space+Grotesk:wght@300;400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap');

    html, body, [class*="css"] {
        font-family: 'Space Grotesk', sans-serif;
    }

    .main-header {
        background: linear-gradient(135deg, #0f172a 0%, #1e293b 50%, #0f4c75 100%);
        padding: 2rem 2.5rem;
        border-radius: 16px;
        margin-bottom: 2rem;
        border: 1px solid #1e40af44;
    }
    .main-header h1 {
        color: #e0f2fe;
        font-size: 2.4rem;
        font-weight: 700;
        margin: 0;
        letter-spacing: -0.5px;
    }
    .main-header p {
        color: #94a3b8;
        margin: 0.4rem 0 0 0;
        font-size: 1rem;
    }

    .risk-badge-safe {
        background: #065f46; color: #6ee7b7;
        padding: 6px 18px; border-radius: 999px;
        font-weight: 700; font-size: 1rem;
        border: 1.5px solid #10b981;
        display: inline-block;
    }
    .risk-badge-adjust {
        background: #78350f; color: #fde68a;
        padding: 6px 18px; border-radius: 999px;
        font-weight: 700; font-size: 1rem;
        border: 1.5px solid #f59e0b;
        display: inline-block;
    }
    .risk-badge-toxic {
        background: #7f1d1d; color: #fca5a5;
        padding: 6px 18px; border-radius: 999px;
        font-weight: 700; font-size: 1rem;
        border: 1.5px solid #ef4444;
        display: inline-block;
    }
    .risk-badge-ineffective {
        background: #312e81; color: #c4b5fd;
        padding: 6px 18px; border-radius: 999px;
        font-weight: 700; font-size: 1rem;
        border: 1.5px solid #8b5cf6;
        display: inline-block;
    }
    .risk-badge-unknown {
        background: #1f2937; color: #9ca3af;
        padding: 6px 18px; border-radius: 999px;
        font-weight: 700; font-size: 1rem;
        border: 1.5px solid #6b7280;
        display: inline-block;
    }

    .metric-card {
        background: #1e293b;
        border: 1px solid #334155;
        border-radius: 12px;
        padding: 1.2rem 1.5rem;
        margin: 0.5rem 0;
    }
    .metric-card h4 {
        color: #64748b;
        font-size: 0.75rem;
        text-transform: uppercase;
        letter-spacing: 1px;
        margin: 0 0 0.4rem 0;
    }
    .metric-card p {
        color: #e2e8f0;
        font-size: 1.1rem;
        font-weight: 600;
        margin: 0;
        font-family: 'JetBrains Mono', monospace;
    }

    .explanation-box {
        background: #0f172a;
        border: 1px solid #1e40af44;
        border-left: 4px solid #3b82f6;
        border-radius: 8px;
        padding: 1rem 1.5rem;
        margin: 0.8rem 0;
    }
    .explanation-box h5 {
        color: #60a5fa;
        font-size: 0.8rem;
        text-transform: uppercase;
        letter-spacing: 1px;
        margin: 0 0 0.5rem 0;
    }
    .explanation-box p {
        color: #cbd5e1;
        line-height: 1.7;
        margin: 0;
    }

    .warning-box {
        background: #451a03;
        border: 1px solid #92400e;
        border-left: 4px solid #f97316;
        border-radius: 8px;
        padding: 1rem 1.5rem;
        margin: 0.8rem 0;
    }

    .critical-box {
        background: #450a0a;
        border: 1px solid #991b1b;
        border-left: 4px solid #ef4444;
        border-radius: 8px;
        padding: 1rem 1.5rem;
        margin: 0.8rem 0;
    }

    .stButton>button {
        background: linear-gradient(135deg, #1d4ed8, #2563eb);
        color: white;
        border: none;
        border-radius: 8px;
        padding: 0.6rem 1.5rem;
        font-weight: 600;
        font-family: 'Space Grotesk', sans-serif;
    }
    .stButton>button:hover {
        background: linear-gradient(135deg, #1e40af, #1d4ed8);
    }

    div[data-testid="stSidebar"] {
        background: #0f172a;
    }

    .json-output {
        font-family: 'JetBrains Mono', monospace;
        font-size: 0.8rem;
        background: #0f172a;
        border: 1px solid #334155;
        border-radius: 8px;
        padding: 1rem;
        max-height: 400px;
        overflow-y: auto;
    }
</style>
""", unsafe_allow_html=True)

# ─── Sidebar ──────────────────────────────────────────────────
with st.sidebar:
    st.markdown("### ⚙️ Configuration")
    groq_api_key = st.text_input(
        "Groq API Key",
        value=os.environ.get("GROQ_API_KEY", ""),
        type="password",
        help="Enter your Groq API key for LLM-generated explanations"
    )
    
    st.markdown("---")
    st.markdown("### 🧬 Supported Genes")
    for g in ["CYP2D6", "CYP2C19", "CYP2C9", "SLCO1B1", "TPMT", "DPYD"]:
        st.markdown(f"• `{g}`")
    
    st.markdown("---")
    st.markdown("### 💊 Supported Drugs")
    for d in DRUG_RISK_TABLE.keys():
        st.markdown(f"• {d.title()}")
    
    st.markdown("---")
    st.markdown("### ℹ️ About")
    st.markdown("""
    **PharmaGuard v1.0**  
    RIFT 2026 Hackathon  
    Pharmacogenomics / XAI Track  
    
    Based on [CPIC Guidelines](https://cpicpgx.org)
    """)

# ─── Main Header ─────────────────────────────────────────────
st.markdown("""
<div class="main-header">
    <h1>🧬 PharmaGuard</h1>
    <p>Pharmacogenomic Risk Prediction System · RIFT 2026 · Powered by Groq LLaMA 3.3</p>
</div>
""", unsafe_allow_html=True)

# ─── Input Section ───────────────────────────────────────────
col1, col2 = st.columns([1.2, 1])

with col1:
    st.markdown("#### 📁 Upload VCF File")
    uploaded_file = st.file_uploader(
        "Drag and drop or browse",
        type=["vcf"],
        help="Upload a standard VCF v4.2 file (max 5MB)"
    )
    
    use_sample = st.checkbox("Use sample VCF for demo", value=False)
    
    if use_sample:
        st.info("🧪 Sample VCF loaded — includes variants for CYP2C19, CYP2D6, SLCO1B1, TPMT, DPYD")

with col2:
    st.markdown("#### 💊 Drug Selection")
    drug_multiselect = st.multiselect(
        "Select drugs to analyze",
        options=list(DRUG_RISK_TABLE.keys()),
        default=["CLOPIDOGREL"],
        format_func=lambda x: x.title()
    )
    
    custom_drugs = st.text_input(
        "Or type additional drugs (comma-separated)",
        placeholder="e.g., CODEINE, WARFARIN",
        help="Supported: CODEINE, WARFARIN, CLOPIDOGREL, SIMVASTATIN, AZATHIOPRINE, FLUOROURACIL"
    )
    
    patient_id = st.text_input(
        "Patient ID (optional)",
        placeholder="Auto-generated if empty"
    )

# ─── Analyze Button ──────────────────────────────────────────
st.markdown("---")
analyze_btn = st.button("🔬 Run Pharmacogenomic Analysis", use_container_width=True)

# ─── Analysis Logic ──────────────────────────────────────────
if analyze_btn:
    # Validate inputs
    all_drugs = list(drug_multiselect)
    if custom_drugs.strip():
        all_drugs += [d.strip().upper() for d in custom_drugs.split(",") if d.strip()]
    all_drugs = list(set(all_drugs))
    
    if not all_drugs:
        st.error("⚠️ Please select at least one drug.")
        st.stop()
    
    if not uploaded_file and not use_sample:
        st.error("⚠️ Please upload a VCF file or enable the sample VCF.")
        st.stop()
    
    if not groq_api_key:
        st.warning("⚠️ No Groq API key provided. LLM explanations will be skipped.")
    
    # Read VCF
    if use_sample:
        vcf_content = get_sample_vcf()
        vcf_filename = "sample_pharmaguard.vcf"
    else:
        vcf_content = uploaded_file.read().decode("utf-8", errors="replace")
        vcf_filename = uploaded_file.name
    
    pid = patient_id.strip() if patient_id.strip() else f"PATIENT_{str(uuid.uuid4())[:8].upper()}"
    
    with st.spinner("🧬 Parsing VCF file..."):
        parsed_vcf = parse_vcf(vcf_content)
    
    with st.spinner("⚙️ Running risk assessment..."):
        risk_results = run_risk_assessment(parsed_vcf, all_drugs)
    
    if groq_api_key:
        with st.spinner("🤖 Generating LLM clinical explanations via Groq..."):
            risk_results = generate_all_explanations(groq_api_key, risk_results)
    else:
        for r in risk_results:
            r["llm_explanation"] = {
                "summary": "LLM explanation not available (no API key provided).",
                "biological_mechanism": "",
                "variant_significance": "",
                "clinical_implications": "",
                "success": False,
            }
    
    # Build JSON outputs
    all_outputs = []
    for result in risk_results:
        output = build_output_schema(
            patient_id=pid,
            drug=result["drug"],
            result=result,
            parsed_vcf=parsed_vcf,
            llm_exp=result.get("llm_explanation", {}),
        )
        all_outputs.append(output)
    
    # ─── Results Display ────────────────────────────────────
    st.markdown("---")
    st.markdown(f"## 📊 Analysis Results — `{pid}`")
    
    # VCF Parse Summary
    with st.expander("📋 VCF Parsing Summary", expanded=False):
        pcol1, pcol2, pcol3 = st.columns(3)
        pcol1.metric("Total Variants Detected", parsed_vcf["total_variants"])
        pcol2.metric("Target Genes Found", len(parsed_vcf["detected_genes"]))
        pcol3.metric("Parse Errors", len(parsed_vcf["parse_errors"]))
        
        if parsed_vcf["detected_genes"]:
            st.markdown(f"**Detected Genes:** `{'`, `'.join(parsed_vcf['detected_genes'])}`")
        if parsed_vcf["parse_errors"]:
            for err in parsed_vcf["parse_errors"][:5]:
                st.warning(err)
    
    # Drug Results
    RISK_BADGE_CLASS = {
        "Safe": "risk-badge-safe",
        "Adjust Dosage": "risk-badge-adjust",
        "Toxic": "risk-badge-toxic",
        "Ineffective": "risk-badge-ineffective",
        "Unknown": "risk-badge-unknown",
    }
    
    RISK_EMOJI = {
        "Safe": "🟢",
        "Adjust Dosage": "🟡",
        "Toxic": "🔴",
        "Ineffective": "🟣",
        "Unknown": "⚪",
    }
    
    for output in all_outputs:
        risk_label = output["risk_assessment"]["risk_label"]
        emoji = RISK_EMOJI.get(risk_label, "⚪")
        badge_class = RISK_BADGE_CLASS.get(risk_label, "risk-badge-unknown")
        drug_name = output["drug"]
        severity = output["risk_assessment"]["severity"].upper()
        confidence = output["risk_assessment"]["confidence_score"]
        
        with st.expander(f"{emoji} {drug_name.title()} — {risk_label} (Severity: {severity})", expanded=True):
            
            # Risk alert for critical cases
            if severity in ("CRITICAL", "HIGH"):
                st.markdown(f"""
                <div class="critical-box">
                    <strong>⚠️ CLINICAL ALERT:</strong> {output['clinical_recommendation']['dosing_recommendation']}
                </div>
                """, unsafe_allow_html=True)
            elif severity == "MODERATE":
                st.markdown(f"""
                <div class="warning-box">
                    <strong>⚠️ CLINICAL NOTE:</strong> {output['clinical_recommendation']['dosing_recommendation']}
                </div>
                """, unsafe_allow_html=True)
            
            # Metric cards
            m1, m2, m3, m4 = st.columns(4)
            with m1:
                st.markdown(f"""<div class="metric-card">
                    <h4>Risk Label</h4>
                    <p><span class="{badge_class}">{risk_label}</span></p>
                </div>""", unsafe_allow_html=True)
            with m2:
                st.markdown(f"""<div class="metric-card">
                    <h4>Primary Gene</h4>
                    <p>{output['pharmacogenomic_profile']['primary_gene']}</p>
                </div>""", unsafe_allow_html=True)
            with m3:
                st.markdown(f"""<div class="metric-card">
                    <h4>Diplotype</h4>
                    <p>{output['pharmacogenomic_profile']['diplotype']}</p>
                </div>""", unsafe_allow_html=True)
            with m4:
                st.markdown(f"""<div class="metric-card">
                    <h4>Phenotype</h4>
                    <p>{output['pharmacogenomic_profile']['phenotype']}</p>
                </div>""", unsafe_allow_html=True)
            
            # Confidence
            st.markdown(f"**Confidence Score:** `{confidence:.0%}`")
            st.progress(confidence)
            
            # Detected Variants
            variants = output["pharmacogenomic_profile"]["detected_variants"]
            if variants:
                st.markdown("**🔬 Detected Variants:**")
                vcols = st.columns([1, 1, 1, 1, 2])
                vcols[0].markdown("**rsID**")
                vcols[1].markdown("**Gene**")
                vcols[2].markdown("**Star Allele**")
                vcols[3].markdown("**REF>ALT**")
                vcols[4].markdown("**Functional Status**")
                for v in variants:
                    vcols[0].markdown(f"`{v['rsid']}`")
                    vcols[1].markdown(f"`{v.get('gene','N/A')}`")
                    vcols[2].markdown(f"`{v.get('star_allele','N/A')}`")
                    vcols[3].markdown(f"`{v.get('ref','?')}>{v.get('alt','?')}`")
                    vcols[4].markdown(v.get("functional_status", "Unknown"))
            else:
                st.info("No pharmacogenomic variants detected for this gene. Wild-type (*1/*1) assumed.")
            
            # Clinical Recommendation
            st.markdown("---")
            st.markdown("**📋 Clinical Recommendation (CPIC)**")
            st.markdown(f"> {output['clinical_recommendation']['dosing_recommendation']}")
            
            if output["clinical_recommendation"]["alternative_drugs"]:
                st.markdown(f"**Alternatives:** {', '.join(output['clinical_recommendation']['alternative_drugs'])}")
            
            if output["clinical_recommendation"]["monitoring_required"]:
                st.markdown(f"**Monitoring:** {output['clinical_recommendation']['monitoring_required']}")
            
            # LLM Explanation
            exp = output["llm_generated_explanation"]
            if exp.get("summary"):
                st.markdown("---")
                st.markdown("**🤖 AI Clinical Explanation (Groq LLaMA 3.3)**")
                
                if exp.get("summary"):
                    st.markdown(f"""<div class="explanation-box">
                        <h5>📌 Summary</h5>
                        <p>{exp['summary']}</p>
                    </div>""", unsafe_allow_html=True)
                
                if exp.get("biological_mechanism"):
                    st.markdown(f"""<div class="explanation-box">
                        <h5>🔬 Biological Mechanism</h5>
                        <p>{exp['biological_mechanism']}</p>
                    </div>""", unsafe_allow_html=True)
                
                if exp.get("variant_significance"):
                    st.markdown(f"""<div class="explanation-box">
                        <h5>🧬 Variant Significance</h5>
                        <p>{exp['variant_significance']}</p>
                    </div>""", unsafe_allow_html=True)
                
                if exp.get("clinical_implications"):
                    st.markdown(f"""<div class="explanation-box">
                        <h5>🏥 Clinical Implications</h5>
                        <p>{exp['clinical_implications']}</p>
                    </div>""", unsafe_allow_html=True)
            
            # JSON Output
            st.markdown("---")
            st.markdown("**📄 JSON Output**")
            json_str = json.dumps(output, indent=2)
            st.download_button(
                label=f"⬇️ Download JSON ({drug_name})",
                data=json_str,
                file_name=f"pharmaguard_{pid}_{drug_name}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json",
                mime="application/json",
                key=f"dl_{drug_name}"
            )
            with st.expander("View Raw JSON", expanded=False):
                st.code(json_str, language="json")
    
    # Download all results
    if len(all_outputs) > 1:
        st.markdown("---")
        all_json = json.dumps(all_outputs, indent=2)
        st.download_button(
            label="⬇️ Download All Results (JSON)",
            data=all_json,
            file_name=f"pharmaguard_{pid}_all_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json",
            mime="application/json",
        )

else:
    # Landing state
    st.markdown("""
    <div style="text-align:center; padding: 3rem; color: #475569;">
        <div style="font-size: 4rem; margin-bottom: 1rem;">🧬</div>
        <h3 style="color: #94a3b8;">Ready to Analyze</h3>
        <p>Upload a VCF file, select drugs, and click <strong>Run Pharmacogenomic Analysis</strong></p>
        <p style="font-size:0.9rem;">Analyzes: CYP2D6 · CYP2C19 · CYP2C9 · SLCO1B1 · TPMT · DPYD</p>
    </div>
    """, unsafe_allow_html=True)