import os
import time
import requests
import pandas as pd
from collections import Counter
import streamlit as st
import py3Dmol
from stmol import showmol
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction


st.set_page_config(page_title="Bioinformatics Sequence Pipeline", layout="wide")


AA_NAMES = {
    'A': 'Alanine', 'C': 'Cysteine', 'D': 'Aspartic Acid', 'E': 'Glutamic Acid',
    'F': 'Phenylalanine', 'G': 'Glycine', 'H': 'Histidine', 'I': 'Isoleucine',
    'K': 'Lysine', 'L': 'Leucine', 'M': 'Methionine', 'N': 'Asparagine',
    'P': 'Proline', 'Q': 'Glutamine', 'R': 'Arginine', 'S': 'Serine',
    'T': 'Threonine', 'V': 'Valine', 'W': 'Tryptophan', 'Y': 'Tyrosine'
}



def fetch_sequence(input_query: str, uploaded_file) -> str:
    """Parses an uploaded FASTA file, fetches an Accession ID, or processes raw text."""
    if uploaded_file is not None:
        string_data = uploaded_file.getvalue().decode("utf-8")
        from io import StringIO
        record = SeqIO.read(StringIO(string_data), "fasta")
        return str(record.seq).upper()
    
    input_query = input_query.strip()
    if not input_query:
        return ""
        
    
    if any(char.isdigit() for char in input_query) and len(input_query) < 20:
        for db in ["nucleotide", "protein"]:
            try:
                handle = Entrez.efetch(db=db, id=input_query, rettype="fasta", retmode="text")
                record = SeqIO.read(handle, "fasta")
                handle.close()
                return str(record.seq).upper()
            except Exception:
                continue
        st.error("Could not fetch sequence for the given Accession ID.")
        return ""
    
    
    return input_query.replace(" ", "").replace("\n", "").upper()


def identify_sequence(seq: str):
    """Determines sequence type, calculates GC content, and identifies gene via NCBI BLAST."""
    seq_set = set(seq)
    dna_bases = set("ACGTN")
    rna_bases = set("ACGUN")
    protein_bases = set("ACDEFGHIKLMNPQRSTVWY")

    if seq_set.issubset(dna_bases):
        seq_type = "DNA"
    elif seq_set.issubset(rna_bases) and "U" in seq_set:
        seq_type = "RNA"
    elif seq_set.issubset(protein_bases):
        seq_type = "Protein"
    else:
        st.error("Invalid sequence characters detected.")
        return None, None, None

    gc_content = gc_fraction(seq) * 100 if seq_type in ["DNA", "RNA"] else None
    
    program = "blastn" if seq_type in ["DNA", "RNA"] else "blastp"
    database = "nt" if seq_type in ["DNA", "RNA"] else "nr"
    base_url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
    headers = {"User-Agent": "StructuredBioPipeline/1.0"}
    
    put_params = {
        "CMD": "Put", "PROGRAM": program, "DATABASE": database,
        "QUERY": seq, "HITLIST_SIZE": "1", "FORMAT_TYPE": "JSON2"
    }
    
    gene_name = "Unknown Gene"
    try:
        res = requests.post(base_url, data=put_params, headers=headers)
        rid = next((line.split("=")[1].strip() for line in res.text.split("\n") if "RID =" in line), None)
        
        if rid:
            with st.spinner("Waiting for NCBI BLAST alignment results..."):
                while True:
                    time.sleep(4)
                    status_res = requests.get(base_url, params={"CMD": "Get", "FORMAT_OBJECT": "SearchInfo", "RID": rid}, headers=headers)
                    if "Status=READY" in status_res.text:
                        break
                    elif "Status=FAILED" in status_res.text:
                        gene_name = "Unknown Gene (BLAST Failed)"
                        break
            
            results_res = requests.get(base_url, params={"CMD": "Get", "RID": rid, "FORMAT_TYPE": "JSON2"}, headers=headers)
            try:
                data = results_res.json()
                hits = data.get("BlastOutput2", {}).get("report", {}).get("results", {}).get("search", {}).get("hits", [])
                gene_name = hits[0]["description"][0]["title"] if hits else "Unknown Gene (No hits found)"
            except Exception:
                gene_name = "Unknown Gene (NCBI API Error)"
        else:
            gene_name = "Unknown Gene (Failed to obtain RID)"
            
    except Exception as e:
        gene_name = f"Unknown Gene (Connection Error: {e})"

    return seq_type, gc_content, gene_name


def central_dogma_pipeline(seq: str, seq_type: str):
    """Handles transcription and translation."""
    bio_seq = Seq(seq)
    transcript = None
    
    if seq_type == "DNA":
        transcript = str(bio_seq.transcribe())
        protein = str(bio_seq.transcribe().translate(to_stop=True))
    elif seq_type == "RNA":
        transcript = seq
        protein = str(bio_seq.translate(to_stop=True))
    else:
        protein = seq

    return transcript, protein


def fetch_pdb_similar(protein_seq: str):
    """Fetches top 5 similar structures from RCSB PDB."""
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    query = {
        "query": {
            "type": "terminal",
            "service": "sequence",
            "parameters": {
                "evalue_cutoff": 1,
                "identity_cutoff": 0.3,
                "target": "pdb_protein_sequence",
                "value": protein_seq
            }
        },
        "return_type": "polymer_entity",
        "request_options": {
            "paginate": {"start": 0, "rows": 5},
            "scoring_strategy": "sequence"
        }
    }
    
    pdb_matches = []
    
    try:
        response = requests.post(url, json=query, timeout=15)
        
        # RCSB returns HTTP 204 when there are no matching results
        if response.status_code == 204 or not response.text.strip():
            return []
            
        response.raise_for_status()
        
        # Safely parse JSON payload
        data = response.json()
        results = data.get("result_set", [])
        
        for item in results:
            pdb_id = item["identifier"].split("_")[0]
            score = item.get("score", 0.0)
            match_pct = score * 100 if score <= 1.0 else score
            
            pdb_matches.append({
                "PDB ID": pdb_id,
                "Sequence Identity (%)": f"{match_pct:.2f}"
            })
            
    except requests.exceptions.Timeout:
        st.error("RCSB PDB API request timed out. Please try again.")
    except requests.exceptions.ConnectionError:
        st.error("Unable to connect to RCSB PDB. Please check your network connection.")
    except requests.exceptions.HTTPError as err:
        st.error(f"RCSB PDB API returned an error: {err}")
    except requests.exceptions.JSONDecodeError:
        st.error("Received invalid JSON response from RCSB PDB API.")
    except Exception as e:
        st.error(f"An unexpected error occurred: {e}")
        
    return pdb_matches

def analyze_amino_acids(protein_seq: str):
    """Calculates top 10 most frequent amino acids."""
    total_aa = len(protein_seq)
    counts = Counter(protein_seq)
    valid_counts = {aa: count for aa, count in counts.items() if aa in AA_NAMES}
    top_10 = Counter(valid_counts).most_common(10)
    
    data = []
    for aa, count in top_10:
        pct = (count / total_aa) * 100
        data.append({
            "Amino Acid": AA_NAMES[aa],
            "Code": aa,
            "Count": count,
            "Percentage (%)": round(pct, 2)
        })
    return pd.DataFrame(data)


def predict_structure_esm(protein_seq: str):
    """Predicts 3D structure using ESMFold API."""
    url = "https://api.esmatlas.com/foldSequence/v1/pdb/"
    res = requests.post(url, data=protein_seq, headers={"Content-Type": "text/plain"})
    if res.status_code == 200:
        return res.text
    return None



st.title("Welcome to The ProtCraft Wizard🧙🏻‍♂️")


st.sidebar.header("Settings")
user_email = st.sidebar.text_input("NCBI Entrez Email", value="your.email@example.com")
Entrez.email = user_email


st.header("1. Input Sequence")
input_option = st.radio("Choose Input Method:", ("Raw Sequence / Accession ID", "Upload FASTA File"))

raw_input = ""
uploaded_file = None

if input_option == "Raw Sequence / Accession ID":
    raw_input = st.text_area("Enter Sequence or Accession ID:", value="", placeholder="e.g. NM_000518 or ATGCG...")
else:
    uploaded_file = st.file_uploader("Upload FASTA file", type=["fasta", "fas", "fa"])

if st.button("Run Pipeline", type="primary"):
    with st.spinner("Processing input sequence..."):
        sequence = fetch_sequence(raw_input, uploaded_file)
        
    if not sequence:
        st.error("No valid sequence input detected. Please provide a sequence, accession ID, or FASTA file.")
    else:
        st.success("Sequence successfully loaded!")
        
        
        st.header("2. Identification & BLAST Analysis")
        with st.spinner("Analyzing sequence type and querying BLAST..."):
            seq_type, gc_content, gene_name = identify_sequence(sequence)
        
        if seq_type:
            col1, col2, col3 = st.columns(3)
            col1.metric("Sequence Type", seq_type)
            col2.metric("GC Content", f"{gc_content:.2f}%" if gc_content is not None else "N/A")
            col3.metric("Sequence Length", f"{len(sequence)} bp/aa")
            st.info(f"**Predicted Gene Name (BLAST):** {gene_name}")

            
            st.header("3. Transcription & Translation")
            transcript, protein_seq = central_dogma_pipeline(sequence, seq_type)
            
            if transcript:
                with st.expander("View mRNA Transcript"):
                    st.text_area("RNA Sequence", transcript, height=100)
            
            with st.expander("View Translated Protein Sequence", expanded=True):
                st.text_area("Protein Sequence", protein_seq, height=100)

            
            st.header("4. Protein Analysis")
            col_a, col_b = st.columns(2)
            
            with col_a:
                st.subheader("Top 10 Amino Acids Frequency")
                df_aa = analyze_amino_acids(protein_seq)
                if not df_aa.empty:
                    st.dataframe(df_aa, use_container_width=True)
                else:
                    st.write("No amino acid data available.")

            with col_b:
                st.subheader("Top PDB Sequence Matches")
                with st.spinner("Searching RCSB PDB..."):
                    matches = fetch_pdb_similar(protein_seq)
                if matches:
                    st.dataframe(pd.DataFrame(matches), use_container_width=True)
                else:
                    st.write("No significant RCSB PDB matches found.")

            
            st.header("5. 3D Structure Prediction")
            if len(protein_seq) > 400:
                st.warning("Protein length exceeds 400 amino acids. ESMFold API predictions are restricted to sequences ≤ 400 residues.")
            else:
                with st.spinner("Predicting 3D structure using ESMFold..."):
                    pdb_data = predict_structure_esm(protein_seq)

                if pdb_data:
                    st.success("3D Structure predicted successfully!")
                    
                    
                    view = py3Dmol.view(width=800, height=500)
                    view.addModel(pdb_data, "pdb")
                    view.setStyle({'cartoon': {'color': 'spectrum'}})
                    view.zoomTo()
                    st.subheader("Interactive 3D Structure Viewer")
                    showmol(view, height=500, width=800)
                    
                    
                    st.download_button(
                        label="Download PDB File",
                        data=pdb_data,
                        file_name="predicted_structure.pdb",
                        mime="chemical/x-pdb"
                    )
                else:
                    st.error("Failed to predict 3D structure using ESMFold API.")
