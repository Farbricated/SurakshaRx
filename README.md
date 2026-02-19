# 🧬 SurakshaRx — Pharmacogenomic Risk Prediction System

> RIFT 2026 Hackathon | Pharmacogenomics / Explainable AI Track | **v5.0 ★**

[![Live Demo](https://img.shields.io/badge/Live-Demo-brightgreen)](YOUR_STREAMLIT_URL_HERE)
[![LinkedIn](https://img.shields.io/badge/LinkedIn-Demo_Video-blue)](YOUR_LINKEDIN_VIDEO_URL_HERE)
[![CPIC](https://img.shields.io/badge/CPIC-Aligned-orange)](https://cpicpgx.org)
[![Python](https://img.shields.io/badge/Python-3.10%2B-blue)](https://python.org)
[![Groq](https://img.shields.io/badge/LLM-LLaMA_3.3_70B-purple)](https://groq.com)

## 🔗 Links
- **Live Demo:** `YOUR_STREAMLIT_URL_HERE`
- **LinkedIn Demo Video:** `YOUR_LINKEDIN_VIDEO_URL_HERE`

---

## 📋 What is SurakshaRx?

Adverse drug reactions kill **over 100,000 Americans annually** — and many are preventable. Most people don't know their genes affect how medicines work in their body. SurakshaRx changes that.

Upload a genetic data file (VCF), select medications, and SurakshaRx instantly tells you:
- Which drugs are **safe** for you
- Which need a **different dose**
- Which are **toxic** or **won't work** based on your DNA
- **Why** — with full AI-generated clinical explanations

All powered by CPIC guidelines (the gold standard in clinical pharmacogenomics) and LLaMA 3.3 70B via Groq.

---

## ✨ v5.0 New Features (Hackathon Edition)

| Feature | What it does | Wow Factor |
|---------|-------------|-----------|
| **🎯 Polygenic Risk Score** | Combines all gene results into a single 0–100 patient risk number | Judges see one clear number, instantly understand severity |
| **🔥 Drug × Gene Heatmap** | Color-coded matrix: all drugs × all genes, hover for details | Most visually striking feature — judges will screenshot it |
| **🧬 Chromosome Visualization** | Shows WHERE on each chromosome your variants live, with glowing animated markers | Nobody else will have this |
| **👥 Population Frequency** | Shows how rare your genotype is vs global population (e.g. "Only 7% of people have this") | Personal and emotional — makes results tangible |
| **💬 Patient Plain-English Mode** | Toggle converts all clinical jargon to language any patient understands | Judges love patient-centric design |
| **⚡ Parallel Test Execution** | ThreadPoolExecutor runs all 4 test scenarios simultaneously | Test suite completes in ~3s instead of ~12s |
| **🔒 Thread-Safe LLM Cache** | RLock-protected cache eliminates duplicate API calls | No more rate-limiting in parallel tests |

---

## 🏗️ Architecture

```
User uploads VCF file + selects drug(s)
        ↓
VCF Parser (vcf_parser.py)
  • Parses VCF v4.2 format
  • Extracts variants for 6 target genes
  • Returns diplotype data per gene
        ↓
Risk Engine (risk_engine.py)
  • Diplotype → Phenotype lookup (CPIC tables)
  • Phenotype × Drug → Risk Label + Severity
  • Returns: risk_label, severity, confidence, CPIC recommendation
        ↓
┌─────────────────────────────────────────────────────┐
│  LLM Explainer (llm_explainer.py) — Thread-safe     │
│  • Thread-safe RLock cache: same drug/phenotype      │
│    served instantly, zero extra API calls            │
│  • LLaMA 3.3 70B via Groq, 3 retries + backoff      │
│  • Static expert templates as fallback               │
│  • All templates cite rsIDs (rubric requirement)     │
│  • skip_llm=True path for instant test mode          │
└─────────────────────────────────────────────────────┘
        ↓
Schema Builder (schema.py)
  • Assembles fully validated JSON output
  • All fields populated — no nulls
        ↓
Drug Interactions (drug_interactions.py)
  • Shared-gene risk detection
  • Known dangerous combination alerts
  • CYP inhibitor phenocopying detection
        ↓
Streamlit UI (app.py v5.0) — 5 NEW VISUAL FEATURES
  • Polygenic Risk Score (0-100 composite gauge)
  • Drug × Gene Heatmap (color matrix)
  • Chromosome Visualization (animated variant markers)
  • Population Frequency bars (your genotype vs world)
  • Patient Plain-English Mode toggle
  • PDF Report download
  • Parallel test suite (ThreadPoolExecutor)
```

---

## 🧬 Supported Genes & Drugs

| Gene | Drug | Clinical Significance |
|------|------|-----------------------|
| **CYP2D6** | Codeine | URM = respiratory death risk; PM = zero pain relief |
| **CYP2C19** | Clopidogrel | PM = heart attack risk — drug won't work at all |
| **CYP2C9** | Warfarin | PM = fatal bleeding risk at standard doses |
| **SLCO1B1** | Simvastatin | Poor Function = rhabdomyolysis (muscle destruction) |
| **TPMT** | Azathioprine | PM = bone marrow destruction at standard doses |
| **DPYD** | Fluorouracil | PM = fatal chemotherapy toxicity |

---

## 📊 Test Suite — All 4 Scenarios

| Scenario | Drugs | Key Variants | Expected Results |
|----------|-------|-------------|-----------------|
| **Mixed Variants** | Clopidogrel, Codeine, Azathioprine | CYP2C19 *2/*3, CYP2D6 *4/*4, TPMT *3B/*3C | Ineffective, Ineffective, Toxic ✅ |
| **UltraRapid Metabolizer** | Codeine, Clopidogrel | CYP2D6 *1xN/*1xN | Toxic, Safe ✅ |
| **All Normal Wild-type** | All 6 drugs | *1/*1 across all genes | All Safe ✅ |
| **Worst Case — All PM** | All 6 drugs | LOF alleles all 6 genes | Ineffective×2, Adjust×1, Toxic×3 ✅ |

All tests run **in parallel** via ThreadPoolExecutor — completes in ~3 seconds.

---

## 🛠️ Tech Stack

| Layer | Technology |
|-------|-----------|
| Frontend + Backend | Streamlit (Python) |
| LLM | Groq — LLaMA 3.3 70B Versatile |
| LLM Fallback | Static expert templates (rsID-cited) |
| Parallel Execution | `concurrent.futures.ThreadPoolExecutor` |
| Thread Safety | `threading.RLock` for LLM cache |
| PDF Reports | fpdf2 |
| Data Validation | Pydantic v2 |
| VCF Parsing | Pure Python (no dependencies) |
| Guidelines | CPIC (cpicpgx.org) — Level A evidence |
| Deployment | Streamlit Community Cloud |

---

## 🚀 Installation & Local Setup

```bash
# 1. Clone the repository
git clone https://github.com/YOUR_USERNAME/pharma-guard.git
cd pharma-guard

# 2. Create virtual environment (recommended)
python3 -m venv venv
source venv/bin/activate        # macOS/Linux
# venv\Scripts\activate         # Windows

# 3. Install dependencies
pip install -r requirements.txt

# 4. Set up environment variables
cp .env.example .env
# Edit .env and add: GROQ_API_KEY=your_key_here

# 5. Run the app
streamlit run app.py
```

The app will open at `http://localhost:8501`

---

## 🔑 Environment Variables

```env
GROQ_API_KEY=your_groq_api_key_here
```

Get your **free** Groq API key at: https://console.groq.com

> **No API key?** The app works fine without one — it falls back to expert-written static templates that cite all relevant rsIDs and star alleles.

---

## 📤 Output JSON Schema

Every analysis produces a fully-validated JSON file:

```json
{
  "patient_id": "PG-A1B2C3D4",
  "drug": "CLOPIDOGREL",
  "timestamp": "2026-02-19T14:30:00Z",
  "risk_assessment": {
    "risk_label": "Ineffective",
    "confidence_score": 0.94,
    "severity": "high"
  },
  "pharmacogenomic_profile": {
    "primary_gene": "CYP2C19",
    "diplotype": "*2/*3",
    "phenotype": "PM",
    "detected_variants": [
      {
        "rsid": "rs4244285",
        "chrom": "chr22",
        "pos": "42522613",
        "ref": "G",
        "alt": "A",
        "gene": "CYP2C19",
        "star_allele": "*2",
        "functional_status": "no_function"
      }
    ]
  },
  "clinical_recommendation": {
    "cpic_guideline": "CPIC Guideline for CLOPIDOGREL",
    "dosing_recommendation": "Use alternative antiplatelet: prasugrel or ticagrelor.",
    "alternative_drugs": ["Prasugrel", "Ticagrelor"],
    "monitoring_required": "Platelet function testing. Monitor for ischaemic events.",
    "contraindicated": false
  },
  "llm_generated_explanation": {
    "summary": "...",
    "biological_mechanism": "...",
    "variant_significance": "...",
    "clinical_implications": "...",
    "model_used": "llama-3.3-70b-versatile"
  },
  "quality_metrics": {
    "vcf_parsing_success": true,
    "variants_detected": 2,
    "genes_analyzed": ["CYP2C19", "CYP2D6", "CYP2C9", "SLCO1B1", "TPMT", "DPYD"],
    "parse_errors": [],
    "explanation_generated": true
  }
}
```

---

## 📁 Project Structure

```
pharma-guard/
├── app.py                          # Main Streamlit app (v5.0 — all new features)
├── vcf_parser.py                   # VCF v4.2 parsing & variant extraction
├── risk_engine.py                  # CPIC-based risk prediction engine
├── llm_explainer.py                # Groq LLM explainer (thread-safe, skip_llm mode)
├── schema.py                       # Output JSON schema builder
├── drug_interactions.py            # Drug-drug interaction checker
├── pdf_report.py                   # Clinical PDF report generator
├── requirements.txt                # Python dependencies
├── .env.example                    # Environment variable template
├── .gitignore                      # Git ignore (excludes .env, patient data)
├── sample_data/
│   ├── sample.vcf                  # Mixed variants test case
│   ├── test_ultrarapid_metabolizer.vcf
│   ├── test_all_normal_wildtype.vcf
│   └── test_worst_case_all_pm.vcf
└── README.md                       # This file
```

---

## 📖 Usage Examples

### Example 1 — Run the sample analysis
1. Open the app at your live URL
2. Select **"Mixed Variants (Standard)"** from the scenario dropdown
3. Select drugs: **Clopidogrel + Codeine + Azathioprine**
4. Click **Run Analysis →**
5. See: Polygenic Risk Score, Drug×Gene Heatmap, Chromosome Visualization

### Example 2 — Patient Mode
1. Run any analysis
2. Toggle **"🧬 → 💬 Patient Plain-English Mode"**
3. All clinical jargon converts to language any patient understands
4. Perfect for patient-facing applications

### Example 3 — Run all tests in parallel
1. Go to the **Test Suite** tab
2. Click **Run All 4 Tests →** (no API key needed — uncheck LLM)
3. All 4 scenarios run simultaneously, results in ~3 seconds
4. Download full test suite JSON for submission

### Example 4 — API integration
```python
import json

result = json.load(open("SurakshaRx_PATIENT_001_WARFARIN.json"))
risk      = result["risk_assessment"]["risk_label"]        # "Adjust Dosage"
gene      = result["pharmacogenomic_profile"]["primary_gene"]  # "CYP2C9"
diplotype = result["pharmacogenomic_profile"]["diplotype"]     # "*2/*3"
rec       = result["clinical_recommendation"]["dosing_recommendation"]
alts      = result["clinical_recommendation"]["alternative_drugs"]
```

---

## 🆕 Changelog

### v5.0 (RIFT 2026 Hackathon Final)
- ✨ **Polygenic Risk Score** — composite 0-100 patient risk score across all genes/drugs
- ✨ **Drug × Gene Heatmap** — interactive color matrix with hover details
- ✨ **Chromosome Visualization** — animated variant location markers on chromosome diagrams
- ✨ **Population Frequency** — how rare is this patient's genotype vs global population
- ✨ **Patient Plain-English Mode** — toggle converts all clinical jargon to plain language
- ✨ **Parallel test execution** — ThreadPoolExecutor, 4x faster test suite
- ✨ **Thread-safe LLM cache** — RLock prevents duplicate API calls in parallel mode

### v4.2
- Thread-safe RLock cache for LLM explainer
- Reduced backoff: 5s/15s/30s → 1s/3s/8s
- All static templates now cite rsIDs

### v4.1
- Fixed TEST_SUITE expected values (CODEINE *4/*4 → Ineffective, AZATHIOPRINE *3B/*3C → Toxic)
- Version bump

### v4.0
- LLM retry + fallback to static templates
- Drug-drug interaction checker
- PDF report generation
- Drug-specific monitoring recommendations

---

## 👥 Team

| Name | Role |
|------|------|
| Member 1 | Lead Developer |
| Member 2 | Clinical Domain Expert |

---

## 🏆 Hackathon Submission

**Event:** RIFT 2026  
**Track:** Pharmacogenomics / Explainable AI  
**Hashtags:** `#RIFT2026` `#SurakshaRx` `#Pharmacogenomics` `#AIinHealthcare` `#ExplainableAI`

---

## ⚠️ Disclaimer

SurakshaRx is a hackathon prototype for educational purposes only. It is **NOT** a medical device and should **NOT** be used for clinical decision-making without validation by a qualified clinical pharmacologist or geneticist. All recommendations should be verified against current CPIC guidelines at [cpicpgx.org](https://cpicpgx.org).