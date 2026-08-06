from collections import Counter
import io
import urllib.parse
import xml.etree.ElementTree as ET

from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
import pandas as pd
import requests
import streamlit as st

st.set_page_config(page_title="ProtCraft Wizard", layout="wide", initial_sidebar_state="expanded")

AA_NAMES = {
    "A": "Alanine", "C": "Cysteine", "D": "Aspartic Acid", "E": "Glutamic Acid",
    "F": "Phenylalanine", "G": "Glycine", "H": "Histidine", "I": "Isoleucine",
    "K": "Lysine", "L": "Leucine", "M": "Methionine", "N": "Asparagine",
    "P": "Proline", "Q": "Glutamine", "R": "Arginine", "S": "Serine",
    "T": "Threonine", "V": "Valine", "W": "Tryptophan", "Y": "Tyrosine"
}

def inject_theme():
    st.markdown("""
    <style>
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;600;800&family=JetBrains+Mono:wght@400;700&display=swap');
    html, body, [class*="st-"] { font-family: 'Inter', sans-serif !important; }
    code, pre { font-family: 'JetBrains Mono', monospace !important; }
    
    @keyframes bgMove {
        0% { background-position: 0% 50%; }
        50% { background-position: 100% 50%; }
        100% { background-position: 0% 50%; }
    }

    .stApp {
        background: linear-gradient(135deg, #020617, #090d1a, #030712);
        background-size: 300% 300%;
        animation: bgMove 15s ease infinite;
        color: #f1f5f9;
    }

    .stApp::before {
        content: "";
        position: fixed;
        top: -30%; left: -30%; width: 160%; height: 160%;
        background: radial-gradient(circle at 50% 50%, rgba(56,189,248,0.18), rgba(129,140,248,0.22), transparent 70%);
        filter: blur(90px);
        animation: bgMove 12s ease infinite alternate;
        z-index: 0;
        pointer-events: none;
    }

    .block-container { position: relative; z-index: 1; padding-top: 2rem; }

    div[data-testid="stExpander"], div[data-testid="stMetric"] {
        background: rgba(255, 255, 255, 0.02) !important;
        backdrop-filter: blur(16px);
        border: 1px solid rgba(255, 255, 255, 0.07) !important;
        border-radius: 14px !important;
        box-shadow: 0 8px 30px rgba(0,0,0,0.4);
    }

    .stButton > button {
        background: linear-gradient(135deg, #0ea5e9 0%, #2563eb 100%) !important;
        color: #fff !important; border-radius: 12px !important;
        font-weight: 600 !important; border: none !important;
        box-shadow: 0 4px 20px rgba(14, 165, 233, 0.4);
        transition: all 0.3s ease;
    }
    .stButton > button:hover { transform: translateY(-2px); box-shadow: 0 6px 25px rgba(14, 165, 233, 0.6); }

    section[data-testid="stSidebar"] { background: rgba(3, 7, 18, 0.9) !important; backdrop-filter: blur(20px); }
    header[data-testid="stHeader"] { background: transparent !important; }

    div[data-testid="stRadio"] div[role="radiogroup"] { display: flex; gap: 16px; }
    div[data-testid="stRadio"] div[role="radiogroup"] label {
        flex: 1; background: rgba(255,255,255,0.02); border: 1px solid rgba(255,255,255,0.08);
        border-radius: 12px; padding: 12px 18px; cursor: pointer; transition: all 0.3s;
    }
    div[data-testid="stRadio"] div[role="radiogroup"] label:has(input:checked) {
        background: rgba(56, 189, 248, 0.15); border-color: #38bdf8;
        box-shadow: 0 0 15px rgba(56,189,248,0.25);
    }
    </style>
    """, unsafe_allow_html=True)

inject_theme()

st.sidebar.header("Settings")
Entrez.email = st.sidebar.text_input("NCBI Email", value="user@example.com")

st.markdown("<h1 style='font-size: 2.5rem; font-weight: 800; margin-bottom: 0;'>ProtCraft <span style='color: #38bdf8;'>Wizard</span> 🧙‍♂️</h1><p style='color: #94a3b8;'>Next-Gen Bioinformatics Pipeline</p>", unsafe_allow_html=True)

input_option = st.radio("Input Method", ["Raw Sequence / Accession ID", "Upload FASTA File"], horizontal=True, label_visibility="collapsed")
raw_input = st.text_area("Input", placeholder="Enter sequence or ID...", height=100) if input_option == "Raw Sequence / Accession ID" else None
uploaded_file = st.file_uploader("Upload FASTA", type=["fasta", "fa"]) if input_option == "Upload FASTA File" else None

if st.button("Run Pipeline", type="primary"):
    seq = ""
    if uploaded_file:
        seq = str(SeqIO.read(io.StringIO(uploaded_file.getvalue().decode("utf-8")), "fasta").seq).upper()
    elif raw_input:
        q = raw_input.strip()
        if any(c.isdigit() for c in q) and len(q) < 20:
            for db in ["nucleotide", "protein"]:
                try:
                    h = Entrez.efetch(db=db, id=q, rettype="fasta", retmode="text")
                    seq = str(SeqIO.read(h, "fasta").seq).upper()
                    h.close()
                    break
                except: continue
        if not seq: seq = q.replace(" ", "").replace("\n", "").upper()

    if not seq:
        st.error("Please provide a valid sequence.")
    else:
        st.success("Sequence loaded successfully!")
        seq_set = set(seq)
        seq_type = "DNA" if seq_set.issubset(set("ACGTN")) else ("RNA" if seq_set.issubset(set("ACGUN")) and "U" in seq_set else "Protein" if seq_set.issubset(set("ACDEFGHIKLMNPQRSTVWY")) else None)
        
        if not seq_type:
            st.error("Invalid sequence characters.")
        else:
            gc = gc_fraction(seq) * 100 if seq_type in ["DNA", "RNA"] else None
            c1, c2, c3 = st.columns(3)
            c1.metric("Type", seq_type)
            c2.metric("GC Content", f"{gc:.2f}%" if gc is not None else "N/A")
            c3.metric("Length", f"{len(seq)} bp/aa")

            prog, db = ("blastn", "nt") if seq_type in ["DNA", "RNA"] else ("blastp", "nr")
            matches = []
            try:
                handle = NCBIWWW.qblast(prog, db, seq, hitlist_size=3)
                root = ET.fromstring(handle.read())
                handle.close()
                for hit in root.findall(".//Hit")[:3]:
                    matches.append({"Gene": hit.find("Hit_def").text, "ID": hit.find("Hit_accession").text})
            except: pass

            if matches:
                st.subheader("Top BLAST Matches")
                st.dataframe(pd.DataFrame(matches), use_container_width=True)

            bio_seq = Seq(seq)
            transcript = str(bio_seq.transcribe()) if seq_type == "DNA" else (seq if seq_type == "RNA" else None)
            protein_seq = str(bio_seq.translate()) if seq_type == "DNA" else (Seq(seq).translate() if seq_type == "RNA" else seq)

            if transcript:
                with st.expander("mRNA Transcript"):
                    st.text_area("RNA", transcript, height=80)

            with st.expander("Translated Protein Sequence", expanded=True):
                st.code(protein_seq, language=None)

            st.subheader("SWISS-MODEL Structure Prediction")
            encoded = urllib.parse.quote(protein_seq)
            st.markdown(f'<a href="https://swissmodel.expasy.org/interactive?target={encoded}" target="_blank"><button style="background:linear-gradient(135deg,#0ea5e9,#2563eb);color:white;padding:10px 24px;border:none;border-radius:8px;font-weight:600;cursor:pointer;">🚀 Open in SWISS-MODEL</button></a>', unsafe_allow_html=True)
