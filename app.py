import os
import time
import requests
import pandas as pd
import xml.etree.ElementTree as ET
from collections import Counter
import streamlit as st
import py3Dmol
from stmol import showmol
from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW
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
    """Determines sequence type, calculates GC content, and identifies top 5 gene matches via NCBI BLAST."""
    seq = seq.strip().upper()
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
        return None, None, []

    gc_content = gc_fraction(seq) * 100 if seq_type in ["DNA", "RNA"] else None

    
    program = "blastn" if seq_type in ["DNA", "RNA"] else "blastp"
    database = "nt" if seq_type in ["DNA", "RNA"] else "nr"
    
    gene_matches = []
    
    try:
        
        result_handle = NCBIWWW.qblast(
            program=program, 
            database=database, 
            sequence=seq,
            hitlist_size=5  
        )
        blast_xml = result_handle.read()
        result_handle.close()

        root = ET.fromstring(blast_xml)
        query_len = float(root.find(".//BlastOutput_query-len").text) if root.find(".//BlastOutput_query-len") is not None else len(seq)
        
        
        for hit in root.findall(".//Hit")[:5]:
            
            accession_elem = hit.find("Hit_accession")
            accession_id = accession_elem.text if accession_elem is not None else "N/A"
            
            
            title_elem = hit.find("Hit_def")
            title = title_elem.text if title_elem is not None else "Unknown Gene"
            
            
            hsp = hit.find(".//Hsp")
            hit_gc = None
            pct_match = 0.0
            
            if hsp is not None:
                
                identity_elem = hsp.find("Hsp_identity")
                align_len_elem = hsp.find("Hsp_align-len")
                if identity_elem is not None and align_len_elem is not None:
                    identity = float(identity_elem.text)
                    align_len = float(align_len_elem.text)
                    pct_match = (identity / align_len) * 100 if align_len > 0 else 0.0
                
                
                hseq_elem = hsp.find("Hsp_hseq")
                if hseq_elem is not None and hseq_elem.text and seq_type in ["DNA", "RNA"]:
                    target_seq = hseq_elem.text.upper().replace("-", "")
                    if target_seq:
                        hit_gc = gc_fraction(target_seq) * 100
            
            gene_matches.append({
                "Gene Name": title,
                "Accession ID": accession_id,
                "GC Content (%)": f"{hit_gc:.2f}" if hit_gc is not None else "N/A",
                "Match Percentage (%)": f"{pct_match:.2f}"
            })
            
    except Exception as e:
        st.warning(f"NCBI BLAST query encounter: {e}")

    return seq_type, gc_content, gene_matches


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
            "type": "terminal", "service": "sequence",
            "parameters": {"evalue_cutoff": 1, "identity_cutoff": 0.3, "target": "pdb_protein_sequence", "value": protein_seq}
        },
        "return_type": "polymer_entity",
        "request_options": {"paginate": {"start": 0, "rows": 5}, "scoring_strategy": "sequence"}
    }
    
    response = requests.post(url, json=query)
    pdb_matches = []
    if response.status_code == 200:
        results = response.json().get("result_set", [])
        for item in results:
            pdb_id = item["identifier"].split("_")[0]
            match_pct = item.get("score", 0) * 100
            pdb_matches.append({"PDB ID": pdb_id, "Sequence Identity (%)": f"{match_pct:.2f}"})
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

def color_protein_sequence_block(seq: str) -> str:
    """Renders sequence with solid colored background blocks matching standard MSA/Clustal tools."""
    
    bg_colors = {
        'A': '#80a0f0', 'I': '#80a0f0', 'L': '#80a0f0', 'M': '#80a0f0', 'F': '#80a0f0', 'W': '#80a0f0', 'V': '#80a0f0', # Blue: Hydrophobic
        'R': '#f01505', 'K': '#f01505',                                                                               # Red: Basic / Positive
        'N': '#00ff00', 'Q': '#00ff00',                                                                               # Green: Polar
        'D': '#c000c0', 'E': '#c000c0',                                                                               # Pink/Magenta: Acidic / Negative
        'C': '#f08080',                                                                                               # Pink-Red: Cysteine
        'G': '#f09040',                                                                                               # Orange: Glycine
        'P': '#ffff00',                                                                                               # Yellow: Proline
        'H': '#15a4a4', 'Y': '#15a4a4',                                                                               # Cyan: Aromatic
        'S': '#15a400', 'T': '#15a400'                                                                                # Dark Green: Hydroxylated
    }

    styled_html = "<div style='font-family: monospace; font-size: 15px; word-break: break-all; line-height: 2.0; background-color: #222; padding: 14px; border-radius: 6px; letter-spacing: 1px;'>"
    
    for aa in seq:
        bg = bg_colors.get(aa, "#ffffff")
       
        text_color = "#ffffff" if aa in ['R', 'K', 'S', 'T', 'D', 'E'] else "#000000"
        
        styled_html += (
            f"<span style='background-color: {bg}; color: {text_color}; "
            f"font-weight: bold; padding: 2px 5px; margin: 1px 0px; "
            f"display: inline-block; text-align: center; border-radius: 2px;'>{aa}</span>"
        )
        
    styled_html += "</div>"
    return styled_html


st.title("Welcome to ProtCraft Wizard🧙‍♂️")


st.sidebar.header("Settings")
user_email = st.sidebar.text_input("NCBI Entrez Email", value="your.email@example.com")
Entrez.email = user_email


st.header("1. Input Sequence🧬")
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
            seq_type, gc_content, gene_matches = identify_sequence(sequence)
        
        if seq_type:
            col1, col2, col3 = st.columns(3)
            col1.metric("Sequence Type", seq_type)
            col2.metric("GC Content", f"{gc_content:.2f}%" if gc_content is not None else "N/A")
            col3.metric("Sequence Length", f"{len(sequence)} bp/aa")
            
            st.subheader("Top 5 Gene Matches (NCBI BLAST)")
            if gene_matches:
                df_matches = pd.DataFrame(gene_matches)
                st.dataframe(df_matches, use_container_width=True)
            else:
                st.info("No significant BLAST hits found.")

           
            st.header("3. Transcription & Translation")
            transcript, protein_seq = central_dogma_pipeline(sequence, seq_type)
            
            if transcript:
                with st.expander("View mRNA Transcript"):
                    st.text_area("RNA Sequence", transcript, height=100)
            
            with st.expander("View Translated Protein Sequence", expanded=True):
                st.markdown(color_protein_sequence(protein_seq), unsafe_allow_html=True)
                st.caption("🟦 Hydrophobic | 🟥 Basic | 🟩 Polar | 🟪 Acidic | 🟧 Glycine | 🟨 Proline")

            
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
