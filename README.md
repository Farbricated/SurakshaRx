# 🧬 PharmaGuard — Pharmacogenomic Risk Prediction System

> RIFT 2026 Hackathon | Pharmacogenomics / Explainable AI Track

[![Live Demo](https://img.shields.io/badge/Live-Demo-blue)](YOUR_STREAMLIT_URL_HERE)
[![LinkedIn](https://img.shields.io/badge/LinkedIn-Demo_Video-blue)](YOUR_LINKEDIN_VIDEO_URL_HERE)

## 🔗 Links
- **Live Demo:** `YOUR_STREAMLIT_URL_HERE`
- **LinkedIn Demo Video:** `YOUR_LINKEDIN_VIDEO_URL_HERE`

---

## 📋 Project Overview

PharmaGuard analyzes patient genetic data (VCF files) and drug names to predict personalized pharmacogenomic risks, providing clinically actionable recommendations with LLM-generated explanations powered by Groq's LLaMA 3.3.

Adverse drug reactions kill over 100,000 Americans annually. Many are preventable through pharmacogenomic testing. PharmaGuard makes this accessible, explainable, and fast.

---

## 🏗️ Architecture

```
User uploads VCF file + selects drug(s)
        ↓
VCF Parser (vcf_parser.py)
  • Parses VCF v4.2 format
  • Extracts variants for 6 target genes
  • Returns diplotype data
        ↓
Risk Engine (risk_engine.py)
  • Diplotype → Phenotype lookup (CPIC)
  • Phenotype × Drug → Risk Label
  • Returns: risk label, severity, confidence
        ↓
Groq LLM Explainer (llm_explainer.py)
  • LLaMA 3.3 70B via Groq API
  • Generates clinical explanation sections
  • Biological mechanism, variant significance
        ↓
Schema Builder (schema.py)
  • Assembles validated JSON output
  • Pydantic schema compliance
        ↓
Streamlit UI (app.py)
  • Color-coded risk display
  • Expandable sections
  • JSON download
```

---

## 🧬 Supported Genes & Drugs

| Gene | Drug(s) |
|------|---------|
| CYP2D6 | Codeine |
| CYP2C19 | Clopidogrel |
| CYP2C9 | Warfarin |
| SLCO1B1 | Simvastatin |
| TPMT | Azathioprine |
| DPYD | Fluorouracil |

---

## 🛠️ Tech Stack

| Layer | Technology |
|-------|-----------|
| Frontend + Backend | Streamlit (Python) |
| LLM API | Groq — LLaMA 3.3 70B Versatile |
| Data Validation | Pydantic v2 |
| VCF Parsing | Pure Python |
| Deployment | Streamlit Community Cloud |
| Guidelines | CPIC (cpicpgx.org) |

---

## 🚀 Installation & Local Setup

```bash
# 1. Clone the repository
git clone https://github.com/YOUR_USERNAME/pharma-guard.git
cd pharma-guard

# 2. Create virtual environment
python -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\activate

# 3. Install dependencies
pip install -r requirements.txt

# 4. Set up environment variables
cp .env.example .env
# Edit .env and add your GROQ_API_KEY

# 5. Run the app
streamlit run app.py
```

---

## 🔑 Environment Variables

```env
GROQ_API_KEY=your_groq_api_key_here
```

Get your free Groq API key at: https://console.groq.com

---

## 📤 API Output Schema

```json
{
  "patient_id": "PATIENT_XXX",
  "drug": "DRUG_NAME",
  "timestamp": "ISO8601_timestamp",
  "risk_assessment": {
    "risk_label": "Safe|Adjust Dosage|Toxic|Ineffective|Unknown",
    "confidence_score": 0.95,
    "severity": "none|low|moderate|high|critical"
  },
  "pharmacogenomic_profile": {
    "primary_gene": "CYP2D6",
    "diplotype": "*4/*4",
    "phenotype": "PM",
    "detected_variants": [{ "rsid": "rs1065852", ... }]
  },
  "clinical_recommendation": {
    "cpic_guideline": "...",
    "dosing_recommendation": "...",
    "alternative_drugs": [...],
    "monitoring_required": "...",
    "contraindicated": false
  },
  "llm_generated_explanation": {
    "summary": "...",
    "biological_mechanism": "...",
    "variant_significance": "...",
    "clinical_implications": "..."
  },
  "quality_metrics": {
    "vcf_parsing_success": true,
    "variants_detected": 3,
    "genes_analyzed": ["CYP2D6"],
    "parse_errors": [],
    "explanation_generated": true
  }
}
```

---

## 📁 Project Structure

```
pharma-guard/
├── app.py              # Main Streamlit application
├── vcf_parser.py       # VCF file parsing & variant extraction
├── risk_engine.py      # CPIC-based risk prediction engine
├── llm_explainer.py    # Groq API LLM explanation generator
├── schema.py           # Pydantic output schema builder
├── requirements.txt    # Python dependencies
├── .env.example        # Environment variable template
├── sample_data/
│   └── sample.vcf      # Test VCF file with multiple gene variants
└── README.md
```

---

## 👥 Team Members

- Member 1 — Name | Role
- Member 2 — Name | Role

---

## 🏆 Hackathon

**RIFT 2026** | Pharmacogenomics / Explainable AI Track  
Hashtags: #RIFT2026 #PharmaGuard #Pharmacogenomics #AIinHealthcare

---

## ⚠️ Disclaimer

PharmaGuard is a hackathon prototype for educational purposes only. It is NOT a medical device and should NOT be used for clinical decision-making without proper validation and regulatory approval.