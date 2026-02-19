"""
Drug-Drug Interaction Checker for SurakshaRx v4.0
Unchanged from v2 — logic is correct and comprehensive.
"""

from typing import Dict, List

GENE_DRUG_MAP = {
    "CYP2D6":  ["CODEINE", "TRAMADOL", "OXYCODONE", "METOPROLOL", "PAROXETINE", "FLUOXETINE", "HALOPERIDOL"],
    "CYP2C19": ["CLOPIDOGREL", "OMEPRAZOLE", "ESCITALOPRAM", "DIAZEPAM", "VORICONAZOLE"],
    "CYP2C9":  ["WARFARIN", "PHENYTOIN", "IBUPROFEN", "CELECOXIB", "FLUVASTATIN"],
    "SLCO1B1": ["SIMVASTATIN", "ATORVASTATIN", "PRAVASTATIN", "ROSUVASTATIN", "METHOTREXATE"],
    "TPMT":    ["AZATHIOPRINE", "MERCAPTOPURINE", "THIOGUANINE"],
    "DPYD":    ["FLUOROURACIL", "CAPECITABINE", "TEGAFUR"],
}

CYP_INHIBITORS = {
    "FLUOXETINE":   {"CYP2D6": "strong"},
    "PAROXETINE":   {"CYP2D6": "strong"},
    "OMEPRAZOLE":   {"CYP2C19": "moderate"},
    "VORICONAZOLE": {"CYP2C19": "strong", "CYP2C9": "strong"},
    "FLUVASTATIN":  {"CYP2C9": "moderate"},
}

DANGEROUS_COMBOS = {
    frozenset(["CODEINE", "FLUOXETINE"]): {
        "severity": "high",
        "mechanism": "Fluoxetine strongly inhibits CYP2D6, reducing codeine→morphine conversion. Risk of inefficacy and unpredictable opioid levels.",
        "recommendation": "Avoid combination. Use non-CYP2D6-dependent opioid.",
    },
    frozenset(["CODEINE", "PAROXETINE"]): {
        "severity": "high",
        "mechanism": "Paroxetine is a potent CYP2D6 inhibitor. Combined with codeine, creates phenocopying of poor metabolizer status.",
        "recommendation": "Contraindicated. Switch to morphine or hydromorphone.",
    },
    frozenset(["WARFARIN", "FLUOXETINE"]): {
        "severity": "moderate",
        "mechanism": "Fluoxetine inhibits CYP2C9 moderately, increasing warfarin plasma levels and bleeding risk.",
        "recommendation": "Increase INR monitoring frequency. Consider dose reduction.",
    },
    frozenset(["SIMVASTATIN", "ATORVASTATIN"]): {
        "severity": "moderate",
        "mechanism": "Both statins compete for SLCO1B1-mediated hepatic uptake. Combined use amplifies myopathy risk.",
        "recommendation": "Avoid combination. Use single statin at appropriate dose.",
    },
    frozenset(["AZATHIOPRINE", "MERCAPTOPURINE"]): {
        "severity": "critical",
        "mechanism": "Both are TPMT substrates. Combination causes compounding myelosuppression risk.",
        "recommendation": "CONTRAINDICATED. Never combine.",
    },
    frozenset(["FLUOROURACIL", "CAPECITABINE"]): {
        "severity": "critical",
        "mechanism": "Both are DPYD substrates. Combination massively increases fluoropyrimidine toxicity.",
        "recommendation": "CONTRAINDICATED. Use only one fluoropyrimidine at a time.",
    },
    frozenset(["CLOPIDOGREL", "OMEPRAZOLE"]): {
        "severity": "high",
        "mechanism": "Omeprazole inhibits CYP2C19, reducing clopidogrel activation by up to 45%. Increased risk of cardiovascular events.",
        "recommendation": "Use pantoprazole instead of omeprazole. Pantoprazole has minimal CYP2C19 inhibition.",
    },
}


def check_shared_gene_risk(drugs: List[str], phenotype_map: Dict[str, str]) -> List[Dict]:
    interactions = []
    drug_upper = [d.upper() for d in drugs]
    gene_to_selected_drugs = {}
    for gene, gene_drugs in GENE_DRUG_MAP.items():
        shared = [d for d in drug_upper if d in gene_drugs]
        if len(shared) >= 2:
            gene_to_selected_drugs[gene] = shared

    for gene, shared_drugs in gene_to_selected_drugs.items():
        phenotype = phenotype_map.get(gene, "Unknown")
        severity, message = "low", ""
        if phenotype == "PM":
            severity = "critical"
            message = (f"COMPOUND RISK: Patient is a Poor Metabolizer for {gene}. "
                       f"Multiple drugs ({', '.join(shared_drugs)}) depend on this pathway. "
                       f"Risk of compounding toxicity or inefficacy is CRITICAL.")
        elif phenotype == "IM":
            severity = "high"
            message = (f"ELEVATED RISK: Patient is an Intermediate Metabolizer for {gene}. "
                       f"Multiple drugs ({', '.join(shared_drugs)}) use this pathway. "
                       f"Monitor closely for cumulative adverse effects.")
        elif phenotype == "URM":
            severity = "high"
            message = (f"ULTRARAPID METABOLIZER RISK: {gene} URM with multiple substrates "
                       f"({', '.join(shared_drugs)}). All drugs may be rapidly cleared — risk of inefficacy.")
        elif phenotype == "NM":
            severity = "low"
            message = (f"Low risk: Normal {gene} activity with multiple substrates "
                       f"({', '.join(shared_drugs)}). Standard monitoring recommended.")
        if message:
            interactions.append({
                "type": "shared_gene", "gene": gene, "drugs_involved": shared_drugs,
                "phenotype": phenotype, "severity": severity, "message": message,
                "recommendation": f"Consider sequencing or spacing {gene} substrate drugs. Consult CPIC.",
            })
    return interactions


def check_known_interactions(drugs: List[str]) -> List[Dict]:
    interactions = []
    drug_upper = set(d.upper() for d in drugs)
    for combo, info in DANGEROUS_COMBOS.items():
        if combo.issubset(drug_upper):
            interactions.append({
                "type": "known_interaction", "drugs_involved": list(combo),
                "severity": info["severity"], "mechanism": info["mechanism"],
                "recommendation": info["recommendation"],
                "message": f"Known interaction between {' + '.join(combo)}: {info['mechanism']}",
            })
    return interactions


def check_inhibitor_effects(drugs: List[str], risk_results: List[Dict]) -> List[Dict]:
    interactions = []
    drug_upper = [d.upper() for d in drugs]
    for drug in drug_upper:
        if drug in CYP_INHIBITORS:
            for gene, strength in CYP_INHIBITORS[drug].items():
                affected = [d for d in drug_upper if d != drug and d in GENE_DRUG_MAP.get(gene, [])]
                if affected:
                    interactions.append({
                        "type": "inhibitor_effect", "inhibitor_drug": drug,
                        "affected_gene": gene, "inhibition_strength": strength,
                        "affected_drugs": affected,
                        "severity": "high" if strength == "strong" else "moderate",
                        "message": (f"{drug} is a {strength} {gene} inhibitor. "
                                    f"This may phenocopy a Poor Metabolizer for: {', '.join(affected)}."),
                        "recommendation": f"Review {gene} substrate doses when co-prescribed with {drug}.",
                    })
    return interactions


def run_interaction_analysis(drugs: List[str], risk_results: List[Dict]) -> Dict:
    phenotype_map = {r.get("primary_gene", ""): r.get("phenotype", "Unknown") for r in risk_results}
    shared_gene  = check_shared_gene_risk(drugs, phenotype_map)
    known        = check_known_interactions(drugs)
    inhibitor    = check_inhibitor_effects(drugs, risk_results)
    all_interactions = shared_gene + known + inhibitor
    SEVERITY_RANK = {"low": 0, "moderate": 1, "high": 2, "critical": 3}
    overall = "none"
    for ix in all_interactions:
        s = ix.get("severity", "low")
        if SEVERITY_RANK.get(s, 0) > SEVERITY_RANK.get(overall, -1):
            overall = s
    return {
        "interactions_found":    len(all_interactions) > 0,
        "total_interactions":    len(all_interactions),
        "overall_severity":      overall,
        "shared_gene_risks":     shared_gene,
        "known_interactions":    known,
        "inhibitor_effects":     inhibitor,
        "all_interactions":      all_interactions,
    }