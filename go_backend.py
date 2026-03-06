from flask import Flask, request, jsonify
from flask_cors import CORS
from Bio import Entrez
import requests
import re
import os
import json

app = Flask(__name__)
CORS(app)

Entrez.email = "akinarli99@gmail.com"

# ─── NCBI Protein Search ──────────────────────────────────────────────────────

def search_ncbi_go(organism, go_term):
    try:
        go_clean = go_term.strip()

        queries = [
            f'"{organism}"[Organism] AND "{go_clean}"[All Fields]',
            f'"{organism}"[Organism] AND {go_clean}[All Fields]',
        ]

        ids = []
        for query in queries:
            handle = Entrez.esearch(db="protein", term=query, retmax="20")
            record = Entrez.read(handle)
            handle.close()
            ids = record.get("IdList", [])
            print(f"[NCBI] Query: {query} → {len(ids)} results")
            if ids:
                break

        if not ids:
            return []

        # Fetch GenBank records to get GO evidence codes
        handle = Entrez.efetch(db="protein", id=",".join(ids[:20]), rettype="gb", retmode="text")
        gb_text = handle.read()
        handle.close()

        # Also get summaries for basic info
        handle = Entrez.esummary(db="protein", id=",".join(ids[:20]), retmode="json")
        data = handle.read()
        handle.close()

        summary = json.loads(data)

        # Parse evidence codes from GenBank flat file
        evidence_map = {}  # accession -> list of evidence codes
        current_acc = None
        for line in gb_text.split("\n"):
            if line.startswith("ACCESSION"):
                current_acc = line.split()[1].strip()
                if current_acc not in evidence_map:
                    evidence_map[current_acc] = []
            if current_acc and "/GO_" in line:
                # Extract evidence code e.g. [Evidence IEA]
                ev_match = re.findall(r'\[Evidence\s+([A-Z]+)\]', line)
                for ev in ev_match:
                    if ev not in evidence_map[current_acc]:
                        evidence_map[current_acc].append(ev)
                # Also extract GO term label
        
        # Parse GO annotations per accession from GenBank
        go_annotations = {}  # accession -> list of GO strings
        current_acc = None
        for line in gb_text.split("\n"):
            if line.startswith("ACCESSION"):
                current_acc = line.split()[1].strip()
                if current_acc not in go_annotations:
                    go_annotations[current_acc] = []
            if current_acc and "/GO_" in line:
                go_match = re.search(r'/(GO_\w+)="(GO:\d+[^"]*)"', line)
                if go_match:
                    go_type = go_match.group(1)  # GO_function, GO_process, GO_component
                    go_val = go_match.group(2)
                    ev_match = re.findall(r'\[Evidence\s+([A-Z]+)\]', line)
                    ev_str = ", ".join(ev_match) if ev_match else ""
                    go_annotations[current_acc].append({
                        "type": go_type,
                        "term": go_val,
                        "evidence": ev_str
                    })

        results = []
        for uid in ids[:20]:
            entry = summary.get("result", {}).get(uid, {})
            if not entry:
                continue
            acc = entry.get("accessionversion", entry.get("caption", uid))
            title = entry.get("title", "")
            organism_name = entry.get("organism", "")

            # Get base accession (without version) for matching
            base_acc = acc.split(".")[0]
            ev_codes = evidence_map.get(acc, evidence_map.get(base_acc, []))
            go_anns = go_annotations.get(acc, go_annotations.get(base_acc, []))

            results.append({
                "accession": acc,
                "name": title,
                "organism": organism_name,
                "evidence_codes": ev_codes,
                "go_annotations": go_anns,
                "link": f"https://www.ncbi.nlm.nih.gov/protein/{acc}"
            })

        return results

    except Exception as e:
        print(f"[NCBI ERROR] {e}")
        import traceback
        traceback.print_exc()
        return []


# ─── UniProt Search ───────────────────────────────────────────────────────────

def search_uniprot_go(organism, go_term):
    try:
        go_clean = go_term.strip()
        is_go_id = bool(re.match(r"GO:\d+", go_clean, re.IGNORECASE))

        if is_go_id:
            go_num = re.sub(r"(?i)go:", "", go_clean)
            query = f'(organism_name:"{organism}") AND (go:{go_num})'
        else:
            query = f'(organism_name:"{organism}") AND (go:"{go_clean}")'

        print(f"[UniProt] Query: {query}")

        url = "https://rest.uniprot.org/uniprotkb/search"
        params = {
            "query": query,
            "format": "json",
            "size": "25",
            "fields": "accession,protein_name,gene_names,organism_name,go_p,go_f,go_c,protein_existence,reviewed"
        }

        r = requests.get(url, params=params, timeout=30)
        print(f"[UniProt] Status: {r.status_code}")

        if not r.ok:
            print(f"[UniProt ERROR] {r.text[:300]}")
            return []

        data = r.json()
        results = []

        # Protein existence mapping
        pe_map = {
            1: "Experimental evidence",
            2: "Transcript level",
            3: "Inferred from homology",
            4: "Predicted",
            5: "Uncertain"
        }

        for entry in data.get("results", []):
            acc = entry.get("primaryAccession", "")

            # Protein name
            pn = entry.get("proteinDescription", {})
            name = ""
            rec = pn.get("recommendedName", {})
            if rec:
                name = rec.get("fullName", {}).get("value", "")
            if not name:
                subs = pn.get("submissionNames", [])
                if subs:
                    name = subs[0].get("fullName", {}).get("value", "")

            # Protein existence
            pe_val = entry.get("proteinExistence", "")
            # UniProt API returns string like "1: Evidence at protein level" or integer
            if isinstance(pe_val, int):
                protein_existence = pe_map.get(pe_val, str(pe_val))
            elif isinstance(pe_val, str):
                # Try to extract number
                pe_num_match = re.match(r"(\d+)", pe_val)
                if pe_num_match:
                    pe_num = int(pe_num_match.group(1))
                    protein_existence = pe_map.get(pe_num, pe_val)
                else:
                    protein_existence = pe_val
            else:
                protein_existence = "Unknown"

            # Reviewed status
            is_reviewed = entry.get("entryType", "") == "UniProtKB reviewed (Swiss-Prot)"

            # GO annotations with evidence
            go_terms = []
            for go_field in ["goBiologicalProcess", "goMolecularFunction", "goCellularComponent"]:
                for go_entry in entry.get(go_field, []):
                    go_id = go_entry.get("goId", "")
                    go_name = go_entry.get("name", "")
                    # Evidence from properties
                    evidences = go_entry.get("evidences", [])
                    ev_codes = list(set([ev.get("evidenceCode", "") for ev in evidences if ev.get("evidenceCode")]))
                    if go_id:
                        go_terms.append({
                            "id": go_id,
                            "name": go_name,
                            "evidence": ", ".join(ev_codes) if ev_codes else ""
                        })

            # Fallback: cross references
            if not go_terms:
                for ref in entry.get("uniProtKBCrossReferences", []):
                    if ref.get("database") == "GO":
                        go_id = ref.get("id", "")
                        props = {p["key"]: p["value"] for p in ref.get("properties", [])}
                        go_name = props.get("GoTerm", "")
                        go_ev = props.get("GoEvidenceType", "")
                        go_terms.append({
                            "id": go_id,
                            "name": go_name,
                            "evidence": go_ev
                        })

            results.append({
                "accession": acc,
                "name": name or "Unknown protein",
                "organism": entry.get("organism", {}).get("scientificName", organism),
                "protein_existence": protein_existence,
                "reviewed": is_reviewed,
                "go_terms": go_terms,
                "link": f"https://www.uniprot.org/uniprotkb/{acc}"
            })

        print(f"[UniProt] Found {len(results)} results")
        return results

    except Exception as e:
        print(f"[UNIPROT ERROR] {e}")
        import traceback
        traceback.print_exc()
        return []


# ─── Endpoint ─────────────────────────────────────────────────────────────────

@app.route("/search", methods=["POST"])
def search():
    data = request.get_json()
    organism = (data or {}).get("organism", "").strip()
    go_term = (data or {}).get("go_term", "").strip()

    if not organism or not go_term:
        return jsonify({"error": "organism ve go_term gerekli"}), 400

    print(f"[INFO] Arama: {organism} | {go_term}")

    ncbi_results = search_ncbi_go(organism, go_term)
    uniprot_results = search_uniprot_go(organism, go_term)

    return jsonify({
        "organism": organism,
        "go_term": go_term,
        "ncbi": ncbi_results,
        "uniprot": uniprot_results,
        "ncbi_total": len(ncbi_results),
        "uniprot_total": len(uniprot_results)
    })


if __name__ == "__main__":
    port = int(os.environ.get("PORT", 5051))
    print(f"GO Search backend http://localhost:{port} adresinde calisiyor...")
    app.run(host="0.0.0.0", port=port, debug=False)
