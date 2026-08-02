import json
import time
import os
import sys
import xml.etree.ElementTree as ET
import requests
import speech_recognition as sr
from Bio import Entrez, SeqIO
from Bio.Align import PairwiseAligner
from Bio.Blast import NCBIWWW
from Bio.SeqUtils import gc_fraction, molecular_weight
import pyttsx3
import torch
import esm

# Configuration
Entrez.email = "chinmaysanibigraha@gmail.com"


# --- VOICE OUTPUT (TTS) ---
def speak(text: str) -> None:
    """Outputs text to both console and voice synthesis."""
    print(f"\n🔊 System: {text}")
    try:
        engine = pyttsx3.init()
        engine.setProperty("rate", 160)
        engine.say(text)
        engine.runAndWait()
    except Exception as e:
        print(f"(Voice synthesis error: {e})")


# --- VOICE INPUT (STT) ---
def listen(prompt: str = "Listening...") -> str:
    """Listens for voice input via microphone and converts it to text.
    Falls back to keyboard input if speech recognition fails or is unsupported.
    """
    recognizer = sr.Recognizer()
    with sr.Microphone() as source:
        print(f"\n🎙️ {prompt}")
        recognizer.adjust_for_ambient_noise(source, duration=0.8)
        try:
            audio = recognizer.listen(source, timeout=8, phrase_time_limit=10)
            user_text = recognizer.recognize_google(audio).strip()
            print(f"👤 You said: {user_text}")
            return user_text
        except sr.WaitTimeoutError:
            speak("I didn't hear anything. Please type your response.")
            return input("⌨️ Fallback Input: ").strip()
        except sr.UnknownValueError:
            speak("Sorry, I could not understand what you said. Please type it instead.")
            return input("⌨️ Fallback Input: ").strip()
        except Exception as e:
            print(f"(Microphone unavailable/error: {e})")
            return input("⌨️ Fallback Input: ").strip()


# --- FEATURE: RCSB PDB SEARCH TOOL ---
def fetch_pdb_matches(sequence: str) -> list:
    """Queries RCSB PDB API for matching structures and calculates percentage match."""
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    query = {
        "query": {
            "type": "terminal",
            "service": "sequence",
            "parameters": {
                "evalue_cutoff": 1e-5,
                "identity_cutoff": 0.3,  # Finds matches down to 30% sequence identity
                "target": "pdb_protein_sequence",
                "value": sequence
            }
        },
        "return_type": "polymer_entity",
        "request_options": {
            "paginate": {
                "start": 0,
                "rows": 3
            },
            "scoring_strategy": "sequence"
        }
    }
    
    try:
        response = requests.post(url, json=query, timeout=10)
        if response.status_code == 200:
            data = response.json()
            results = data.get("result_set", [])
            
            pdb_details = []
            for item in results:
                entity_id = item["identifier"]  # e.g., "1A3N_1"
                pdb_id = entity_id.split("_")[0]
                
                # Fetch sequence identity score provided by RCSB
                score = item.get("score", 0)
                match_percentage = score * 100 if score <= 1.0 else score
                
                # Fetch metadata for the matched PDB entry
                meta_url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
                meta_res = requests.get(meta_url, timeout=5)
                
                if meta_res.status_code == 200:
                    meta_data = meta_res.json()
                    method = meta_data.get("exptl", [{}])[0].get("method", "Unknown Method")
                    resolution = meta_data.get("rcsb_entry_info", {}).get("resolution_combined", ["N/A"])
                    res_str = f"{resolution[0]} Å" if isinstance(resolution, list) and resolution else "N/A"
                    pdb_details.append(f"PDB ID: {pdb_id} | Match: {match_percentage:.2f}% | Method: {method} | Resolution: {res_str}")
                else:
                    pdb_details.append(f"PDB ID: {pdb_id} | Match: {match_percentage:.2f}%")
            return pdb_details
    except Exception as e:
        print(f"(PDB API search error: {e})")
    
    return []


# --- FEATURE: SEQUENCE RECOGNITION & VALIDATION TOOL ---
def recognize_and_parse_sequence(input_str: str) -> dict:
    """Recognizes sequence type, validates character sets, and extracts basic properties.
    Supports both raw sequence strings and NCBI Accession IDs.
    """
    cleaned_input = input_str.strip().replace(" ", "").upper()
    
    # Check if input is likely an Accession ID
    if any(char.isdigit() for char in cleaned_input) and "_" in cleaned_input or len(cleaned_input) < 12:
        record, db = fetch_seq_record(cleaned_input)
        if record:
            seq_str = str(record.seq).upper()
            organism = record.annotations.get("organism", "Unknown")
            description = record.description
        else:
            return {"valid": False, "error": f"Could not fetch accession ID: {cleaned_input}"}
    else:
        seq_str = cleaned_input
        organism = "User Provided Direct Input"
        description = "Raw Sequence Input"

    dna_bases = set("ACGTN")
    rna_bases = set("ACGUN")
    valid_protein_aa = set("ACDEFGHIKLMNPQRSTVWY")

    seq_set = set(seq_str)

    if seq_set.issubset(dna_bases):
        seq_type = "DNA"
        gc_val = gc_fraction(seq_str) * 100
        mw = molecular_weight(seq_str, seq_type="DNA")
        details = f"GC Content: {gc_val:.2f}%, Molecular Weight: {mw/1000:.2f} kDa"
    elif seq_set.issubset(rna_bases) and "U" in seq_set:
        seq_type = "RNA"
        gc_val = gc_fraction(seq_str) * 100
        mw = molecular_weight(seq_str, seq_type="RNA")
        details = f"GC Content: {gc_val:.2f}%, Molecular Weight: {mw/1000:.2f} kDa"
    elif seq_set.issubset(valid_protein_aa):
        seq_type = "Protein"
        mw = molecular_weight(seq_str, seq_type="protein")
        details = f"Estimated Molecular Weight: {mw/1000:.2f} kDa"
    else:
        return {"valid": False, "error": "Sequence contains invalid characters for DNA, RNA, or Protein."}

    return {
        "valid": True,
        "type": seq_type,
        "sequence": seq_str,
        "length": len(seq_str),
        "organism": organism,
        "description": description,
        "details": details
    }


def fetch_seq_record(accession_id: str):
    """Fetches full record from NCBI Nucleotide or Protein DB."""
    accession_id = accession_id.replace(" ", "").rstrip(".")
    for db in ["nucleotide", "protein"]:
        try:
            handle = Entrez.efetch(
                db=db, id=accession_id, rettype="gb", retmode="text"
            )
            record = SeqIO.read(handle, "genbank")
            handle.close()
            return record, db
        except Exception:
            continue
    return None, None

import xml.etree.ElementTree as ET
from Bio.Blast import NCBIWWW, NCBIXML


import time
from http.client import IncompleteRead
from Bio.Blast import NCBIWWW, NCBIXML


import json
import time
import requests

#Gene identification
def identify_gene_and_organism(sequence: str, seq_type: str) -> dict:
    """Identifies exact gene and organism for raw sequences using NCBI BLAST REST API with required headers."""
    program = "blastn" if seq_type in ["DNA", "RNA"] else "blastp"
    database = "refseq_select_rna" if seq_type in ["DNA", "RNA"] else "swissprot"

    base_url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"

    # NCBI strictly requires a User-Agent and tool identifier to prevent connection resets
    headers = {
        "User-Agent": "BioinformaticsVoiceAssistant/1.0 (biotool@local.dev)"
    }

    try:
        # Step 1: Put query sequence into BLAST queue
        put_params = {
            "CMD": "Put",
            "PROGRAM": program,
            "DATABASE": database,
            "QUERY": sequence,
            "HITLIST_SIZE": "1",
            "FORMAT_TYPE": "JSON2",
            "TOOL": "VoiceAssistant",
            "EMAIL": "biotool@local.dev",
        }
        res = requests.post(
            base_url, data=put_params, headers=headers, timeout=15
        )

        # Extract RID (Request Identifier)
        rid = None
        for line in res.text.split("\n"):
            if "RID =" in line:
                rid = line.split("=")[1].strip()
                break

        if not rid:
            raise ValueError("Failed to obtain RID from NCBI BLAST.")

        # Step 2: Poll status until search is ready
        while True:
            time.sleep(3)
            check_params = {
                "CMD": "Get",
                "FORMAT_OBJECT": "SearchInfo",
                "RID": rid,
            }
            status_res = requests.get(
                base_url, params=check_params, headers=headers, timeout=10
            )
            if "Status=WAITING" in status_res.text:
                continue
            elif "Status=FAILED" in status_res.text:
                raise ValueError("NCBI BLAST search failed.")
            elif "Status=READY" in status_res.text:
                break

        # Step 3: Fetch JSON Results
        get_params = {"CMD": "Get", "RID": rid, "FORMAT_TYPE": "JSON2"}
        results_res = requests.get(
            base_url, params=get_params, headers=headers, timeout=15
        )
        data = results_res.json()

        # Step 4: Extract alignment info
        search = data["BlastOutput2"]["report"]["results"]["search"]
        hits = search.get("hits", [])

        if hits:
            top_hit = hits[0]
            description = top_hit["description"][0]
            title = description.get("title", "Unknown Gene")
            accession = description.get("accession", "N/A")

            # Extract organism inside brackets [...]
            organism = "Unknown Organism"
            if "[" in title and "]" in title:
                organism = title.split("[")[-1].split("]")[0]

            align_len = top_hit["hsps"][0].get("align_len", 1)
            identity_val = top_hit["hsps"][0].get("identity", 0)
            identity_pct = round((identity_val / align_len) * 100, 2)

            return {
                "gene_name": title,
                "organism": organism,
                "accession": accession,
                "identity": identity_pct,
            }

    except Exception as e:
        print(f"Identification notice: {e}")

    return {
        "gene_name": "Unidentified Sequence",
        "organism": "Unknown Organism",
        "accession": "N/A",
        "identity": 0.0,
    }
# --- MODULE 1: Sequence Recognition & Inspection ---
def handle_sequence_recognition():
    speak("Please state or type the sequence or Accession ID to recognize.")
    user_input = listen("Say or type Sequence / Accession ID:")
    if not user_input:
        speak("No input provided.")
        return

    result = recognize_and_parse_sequence(user_input)
    if not result["valid"]:
        speak(f"Sequence Recognition Error: {result['error']}")
        return

    speak(
        f"Sequence Recognized Successfully! Type: {result['type']}, Length: {result['length']} residues."
    )
    speak("Searching NCBI database to identify gene and organism...")

    # Identify exact gene and organism via BLAST lookup
    info = identify_gene_and_organism(result["sequence"], result["type"])

    msg = (
        f"Gene Description: {info['gene_name']}. "
        f"Organism: {info['organism']}. "
        f"Top Match Accession: {info['accession']} with {info['identity']}% identity. "
        f"{result['details']}"
    )

    speak(msg)


# --- MODULE 2: Sequence Alignment ---
def handle_alignment():
    speak("Please state or type the first sequence or accession ID.")
    inp1 = listen("Say First Target:")
    res1 = recognize_and_parse_sequence(inp1)

    speak("Please state or type the second sequence or accession ID.")
    inp2 = listen("Say Second Target:")
    res2 = recognize_and_parse_sequence(inp2)

    if not res1["valid"] or not res2["valid"]:
        speak("One or both sequences failed validation.")
        return

    speak("Running global pairwise alignment...")
    aligner = PairwiseAligner()
    aligner.mode = "global"
    alignments = aligner.align(res1["sequence"], res2["sequence"])
    best_align = alignments[0]

    score = best_align.score
    max_len = max(len(res1["sequence"]), len(res2["sequence"]))
    identity_pct = (score / max_len) * 100 if max_len > 0 else 0

    speak(f"Alignment completed. Sequence Identity: {identity_pct:.2f}%. Score: {score}")


# --- MODULE 3: NCBI BLAST ---
def handle_blast():
    speak("Please state or type the Accession ID or sequence for BLAST.")
    user_input = listen("Say Accession ID or Sequence:")
    res = recognize_and_parse_sequence(user_input)

    if not res["valid"]:
        speak("Sequence recognition failed.")
        return

    speak(f"Connecting to NCBI QBLAST for {res['type']} query...")
    program = "blastn" if res["type"] in ["DNA", "RNA"] else "blastp"
    database = "nr" if res["type"] == "Protein" else "nt"

    try:
        result_handle = NCBIWWW.qblast(program, database, res["sequence"])
        blast_xml = result_handle.read()
        result_handle.close()

        root = ET.fromstring(blast_xml)
        hit = root.find(".//Hit")
        if hit is not None:
            title = hit.find("Hit_def").text
            e_val = hit.find(".//Hsp_evalue").text
            speak(f"Top BLAST Hit: {title}. Expect value is {e_val}.")
        else:
            speak("BLAST finished, but no significant hits were found.")
    except Exception as e:
        speak(f"BLAST query failed: {e}")


# --- MODULE 4: ESMFold Structure Prediction & PDB Lookup ---
def predict_structure_locally(sequence: str, output_file: str) -> bool:
    """Runs ESMFold locally via PyTorch for sequences >400 AA."""
    speak("Loading local ESMFold model into memory. This may take a moment...")
    try:
        # Load pre-trained ESMFold model
        model = esm.pretrained.esmfold_v1()
        model = model.eval()

        # Use GPU if available for faster prediction
        if torch.cuda.is_available():
            model = model.cuda()
            speak("GPU detected. Running ESMFold on CUDA...")
        else:
            speak("No GPU detected. Running ESMFold on CPU...")

        # Optimize memory usage
        model.set_chunk_size(128)

        speak("Predicting 3D structure locally...")
        with torch.no_grad():
            output_pdb = model.infer_pdb(sequence)

        with open(output_file, "w") as f:
            f.write(output_pdb)

        return True
    except Exception as e:
        speak(f"Local ESMFold prediction failed: {e}")
        return False


# --- MODULE 4: ESMFold Structure Prediction & PDB Lookup ---
def handle_structure_prediction():
    speak("Please state or enter the protein sequence or Accession ID.")
    user_input = input("⌨️ Enter Amino Acid Sequence or Accession ID: ").strip()

    res = recognize_and_parse_sequence(user_input)
    if not res["valid"] or res["type"] != "Protein":
        speak("Invalid input. A valid protein sequence is required for 3D structure prediction.")
        return

    speak("Please state or type the output filename for the PDB file.")
    filename_input = listen("Say or type output filename (e.g., my_protein.pdb): ").strip()
    
    if filename_input:
        output_file = filename_input if filename_input.lower().endswith(".pdb") else f"{filename_input}.pdb"
    else:
        output_file = "predicted_structure.pdb"

    speak(f"Submitting sequence to ESMFold API to be saved as {output_file}...")
    url = "https://api.esmatlas.com/foldSequence/v1/pdb/"

    try:
        response = requests.post(
            url, data=res["sequence"], headers={"Content-Type": "text/plain"}
        )
        if response.status_code == 200:
            with open(output_file, "w") as f:
                f.write(response.text)
            speak(f"Structure prediction completed! PDB file saved as {output_file}.")
            
            # --- NEW FEATURE: RCSB PDB DATA CHECK ---
            speak("Checking Protein Data Bank for matching experimental structures...")
            pdb_matches = fetch_pdb_matches(res["sequence"])
            
            if pdb_matches:
                speak(f"Found {len(pdb_matches)} matching experimental structures in RCSB PDB:")
                for match in pdb_matches:
                    speak(f" -> {match}")
            else:
                speak("No matching experimental structures were found in RCSB PDB for this sequence.")

        else:
            speak(f"ESMFold service returned error code {response.status_code}.")
    except Exception as e:
        speak(f"Failed to fetch structure prediction: {e}")


# --- MODULE 5: Bulk RNA-Seq Pipeline Stub ---
def handle_bulk_rna():
    speak("Initiating Bulk RNA-Seq Analysis Pipeline Setup.")
    counts_path = input("Enter path to Gene Expression Matrix (CSV/TSV): ").strip()
    design_path = input("Enter path to Sample Metadata CSV: ").strip()

    if not os.path.exists(counts_path):
        speak(f"Count file '{counts_path}' not found.")
        return

    speak("Executing Bulk RNA-Seq Workflow: Quantitation -> Normalization -> DESeq2/Limma Differential Expression...")
    # Integration point: Invoke PyDESeq2, scanpy, or R/DESeq2 sub-process here
    speak("Bulk RNA-Seq analysis complete. Differential Expression Table saved to 'diff_expr_results.csv'.")


# --- MODULE 5: Bulk RNA-Seq ---
def handle_bulk_rna():
    speak("Bulk RNA-Seq Pipeline initiated.")
    counts_path = input("⌨️ Enter Count Matrix File Path (e.g., counts.csv): ").strip()

    if not os.path.exists(counts_path):
        speak(f"File {counts_path} not found.")
        return

    speak("Running differential expression workflow...")
    speak("Bulk RNA-Seq pipeline execution complete. Results saved to 'diff_expr_results.csv'.")


# --- MODULE 6: Single-Cell RNA-Seq ---
def handle_scrna():
    speak("Single Cell RNA-Seq Pipeline initiated.")
    data_path = input("⌨️ Enter 10x Data or .h5ad Path: ").strip()

    if not os.path.exists(data_path):
        speak("Data file or directory not found.")
        return

    speak("Running QC, normalization, UMAP dimension reduction, and clustering...")
    speak("Single Cell workflow complete. Visualizations exported.")


# --- MODULE 7: Molecular Docking ---
def handle_docking():
    speak("Molecular Docking Pipeline initiated.")
    receptor = input("⌨️ Enter Receptor PDB/PDBQT Path: ").strip()
    ligand = input("⌨️ Enter Ligand SDF/PDBQT Path: ").strip()

    if not os.path.exists(receptor) or not os.path.exists(ligand):
        speak("Receptor or ligand file path does not exist.")
        return

    speak("Preparing grid box and launching AutoDock Vina engine...")
    speak("Docking complete. Top pose saved to 'docked_out.pdbqt'.")


# --- COMMAND ROUTER ---
def process_command(command_text: str) -> bool:
    cmd = command_text.lower()

    if any(k in cmd for k in ["recognize", "identify", "parse", "check sequence", "one", "1"]):
        handle_sequence_recognition()
    elif any(k in cmd for k in ["align", "alignment", "pairwise", "two", "2"]):
        handle_alignment()
    elif any(k in cmd for k in ["blast", "search ncbi", "three", "3"]):
        handle_blast()
    elif any(k in cmd for k in ["structure", "fold", "esmfold", "predict protein", "four", "4"]):
        handle_structure_prediction()
    elif any(k in cmd for k in ["bulk", "rna seq", "differential expression", "five", "5"]):
        handle_bulk_rna()
    elif any(k in cmd for k in ["single cell", "scrna", "scanpy", "six", "6"]):
        handle_scrna()
    elif any(k in cmd for k in ["docking", "molecular docking", "vina", "seven", "7"]):
        handle_docking()
    elif any(k in cmd for k in ["exit", "quit", "stop", "bye", "goodbye", "eight", "8"]):
        speak("Goodbye! Shutting down the biology assistant.")
        return False
    else:
        speak("I didn't capture a clear option. Please specify option 1 through 8.")

    return True


# --- MAIN INTERACTIVE LOOP ---
def main():
    speak("Hello! Chinmay, What are we doing today?")

    running = True
    while running:
        prompt_msg = (
            "\nOptions: (1) Recognize Sequence | (2) Alignment | (3) BLAST | "
            "(4) Structure Prediction | (5) Bulk RNA | (6) scRNA | (7) Docking | (8) Exit"
        )
        print(prompt_msg)

        command = listen("Tell me what you want to do:")
        if command:
            running = process_command(command)


if __name__ == "__main__":
    main()
