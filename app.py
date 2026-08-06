from collections import Counter
import io
import os
import time
import urllib.parse
import xml.etree.ElementTree as ET
from textwrap import dedent

from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
import pandas as pd
import requests
import streamlit as st
import streamlit.components.v1 as components


st.set_page_config(
    page_title="ProtCraft Wizard | Bioinformatics Suite", 
    layout="wide",
    initial_sidebar_state="expanded"
)


AA_NAMES = {
    "A": "Alanine",
    "C": "Cysteine",
    "D": "Aspartic Acid",
    "E": "Glutamic Acid",
    "F": "Phenylalanine",
    "G": "Glycine",
    "H": "Histidine",
    "I": "Isoleucine",
    "K": "Lysine",
    "L": "Leucine",
    "M": "Methionine",
    "N": "Asparagine",
    "P": "Proline",
    "Q": "Glutamine",
    "R": "Arginine",
    "S": "Serine",
    "T": "Threonine",
    "V": "Valine",
    "W": "Tryptophan",
    "Y": "Tyrosine",
}


def inject_custom_ui_theme():
    """Injects high-end animated neon UI, floating particles, and interactive cards."""
    css = """
    <style>
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700;800&family=JetBrains+Mono:wght@400;500;700&display=swap');

    html, body, [class*="st-"] {
        font-family: 'Inter', sans-serif !important;
    }

    code, pre, .stDataFrame {
        font-family: 'JetBrains Mono', monospace !important;
    }

    @keyframes breathingGradient {
        0% { background-position: 0% 50%; }
        50% { background-position: 100% 50%; }
        100% { background-position: 0% 50%; }
    }

    .stApp {
        background: linear-gradient(-45deg, #020617, #090d1a, #040711, #020308);
        background-size: 400% 400%;
        animation: breathingGradient 20s ease infinite;
        color: #f1f5f9;
    }

    /* Floating Neon Glow Animation Background */
    @keyframes floatOrb1 {
        0% { transform: translate(0px, 0px) scale(1); opacity: 0.25; }
        50% { transform: translate(80px, 100px) scale(1.25); opacity: 0.45; }
        100% { transform: translate(-40px, 50px) scale(0.9); opacity: 0.25; }
    }

    @keyframes floatOrb2 {
        0% { transform: translate(0px, 0px) scale(1); opacity: 0.2; }
        50% { transform: translate(-100px, -70px) scale(1.3); opacity: 0.4; }
        100% { transform: translate(60px, -40px) scale(1); opacity: 0.2; }
    }

    .stApp::before {
        content: "";
        position: fixed;
        top: 10%; left: 15%; width: 450px; height: 450px;
        background: radial-gradient(circle, rgba(14, 165, 233, 0.22), transparent 70%);
        filter: blur(90px);
        pointer-events: none;
        z-index: 0;
        animation: floatOrb1 14s ease-in-out infinite alternate;
    }

    .stApp::after {
        content: "";
        position: fixed;
        bottom: 10%; right: 15%; width: 500px; height: 500px;
        background: radial-gradient(circle, rgba(168, 85, 247, 0.2), transparent 70%);
        filter: blur(100px);
        pointer-events: none;
        z-index: 0;
        animation: floatOrb2 18s ease-in-out infinite alternate;
    }

    .block-container {
        position: relative;
        z-index: 1;
        padding-top: 2.5rem;
        padding-bottom: 3.5rem;
    }

    /* Glassmorphism Cards with smooth hover glow */
    div[data-testid="stExpander"], div[data-testid="stMetric"] {
        background: rgba(255, 255, 255, 0.02) !important;
        backdrop-filter: blur(18px);
        -webkit-backdrop-filter: blur(18px);
        border: 1px solid rgba(255, 255, 255, 0.07) !important;
        border-radius: 16px !important;
        padding: 16px !important;
        box-shadow: 0 10px 30px 0 rgba(0, 0, 0, 0.4);
        position: relative;
        z-index: 2;
        transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
    }

    div[data-testid="stExpander"]:hover {
        border-color: rgba(56, 189, 248, 0.35) !important;
        box-shadow: 0 12px 35px 0 rgba(56, 189, 248, 0.12);
        transform: translateY(-2px);
    }

    h1, h2, h3 {
        color: #f8fafc !important;
        font-weight: 700 !important;
        letter-spacing: -0.03em;
    }

    /* Animated Neon Buttons */
    .stButton > button {
        background: linear-gradient(135deg, #0ea5e9 0%, #2563eb 100%) !important;
        color: #ffffff !important;
        border: 1px solid rgba(255, 255, 255, 0.2) !important;
        border-radius: 12px !important;
        padding: 0.75rem 2rem !important;
        font-weight: 600 !important;
        letter-spacing: 0.02em;
        transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1) !important;
        box-shadow: 0 4px 20px rgba(14, 165, 233, 0.4);
        position: relative;
        z-index: 2;
    }

    .stButton > button:hover {
        transform: translateY(-3px) scale(1.02);
        box-shadow: 0 8px 30px rgba(14, 165, 233, 0.7);
        border-color: rgba(255, 255, 255, 0.4) !important;
    }

    /* Sidebar Glassmorphism */
    section[data-testid="stSidebar"] {
        background-color: rgba(4, 6, 12, 0.92) !important;
        backdrop-filter: blur(22px);
        border-right: 1px solid rgba(255, 255, 255, 0.06);
        z-index: 10;
    }
    
    header[data-testid="stHeader"] {
        background: transparent !important;
    }

    /* Sleek Radio Navigation Cards */
    div[data-testid="stRadio"] div[role="radiogroup"] {
        display: flex !important;
        flex-direction: row !important;
        gap: 16px !important;
        width: 100% !important;
        align-items: stretch !important;
    }

    div[data-testid="stRadio"] div[role="radiogroup"] label {
        flex: 1 1 0px !important;
        min-height: 78px !important;
        background: rgba(255, 255, 255, 0.02) !important;
        backdrop-filter: blur(12px) !important;
        border: 1px solid rgba(255, 255, 255, 0.08) !important;
        border-radius: 14px !important;
        padding: 14px 20px !important;
        transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1) !important;
        cursor: pointer !important;
        display: flex !important;
        align-items: center !important;
        justify-content: flex-start !important;
        box-shadow: 0 4px 16px rgba(0,0,0,0.2);
    }

    div[data-testid="stRadio"] div[role="radiogroup"] label:hover {
        background: rgba(56, 189, 248, 0.08) !important;
        border-color: rgba(56, 189, 248, 0.4) !important;
        transform: translateY(-2px);
    }

    div[data-testid="stRadio"] div[role="radiogroup"] label:has(input:checked) {
        background: linear-gradient(135deg, rgba(14, 165, 233, 0.18) 0%, rgba(37, 99, 235, 0.18) 100%) !important;
        border-color: #38bdf8 !important;
        box-shadow: 0 0 22px rgba(56, 189, 248, 0.3) !important;
    }

    textarea, input {
        background-color: rgba(255, 255, 255, 0.02) !important;
        border: 1px solid rgba(255, 255, 255, 0.08) !important;
        color: #f8fafc !important;
        border-radius: 12px !important;
    }
    textarea:focus, input:focus {
        border-color: #38bdf8 !important;
        box-shadow: 0 0 12px rgba(56, 189, 248, 0.25) !important;
    }
    </style>
    """
    st.markdown(css, unsafe_allow_html=True)


def fetch_sequence(input_query: str, uploaded_file) -> str:
    if uploaded_file is not None:
        string_data = uploaded_file.getvalue().decode("utf-8")
        record = SeqIO.read(io.StringIO(string_data), "fasta")
        return str(record.seq).upper()
    input_query = input_query.strip()
    if not input_query:
        return ""
    if any(char.isdigit() for char in input_query) and len(input_query) < 20:
        for db in ["nucleotide", "protein"]:
            try:
                handle = Entrez.efetch(
                    db=db, id=input_query, rettype="fasta", retmode="text"
                )
                record = SeqIO.read(handle, "fasta")
                handle.close()
                return str(record.seq).upper()
            except Exception:
                continue
    return input_query.replace(" ", "").replace("\n", "").upper()


def identify_sequence(seq: str):
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
            program=program, database=database, sequence=seq, hitlist_size=5
        )
        blast_xml = result_handle.read()
        result_handle.close()
        root = ET.fromstring(blast_xml)
        for hit in root.findall(".//Hit")[:5]:
            accession_elem = hit.find("Hit_accession")
            accession_id = (
                accession_elem.text if accession_elem is not None else "N/A"
            )
            title_elem = hit.find("Hit_def")
            title = (
                title_elem.text if title_elem is not None else "Unknown Gene"
            )
            hsp = hit.find(".//Hsp")
            hit_gc, pct_match = None, 0.0
            if hsp is not None:
                identity_elem = hsp.find("Hsp_identity")
                align_len_elem = hsp.find("Hsp_align-len")
                if (
                    identity_elem is not None
                    and align_len_elem is not None
                    and float(align_len_elem.text) > 0
                ):
                    pct_match = (
                        float(identity_elem.text)
                        / float(align_len_elem.text)
                    ) * 100
                hseq_elem = hsp.find("Hsp_hseq")
                if (
                    hseq_elem is not None
                    and hseq_elem.text
                    and seq_type in ["DNA", "RNA"]
                ):
                    target_seq = hseq_elem.text.upper().replace("-", "")
                    if target_seq:
                        hit_gc = gc_fraction(target_seq) * 100
            gene_matches.append(
                {
                    "Gene Name": title,
                    "Accession ID": accession_id,
                    "GC Content (%)": f"{hit_gc:.2f}"
                    if hit_gc is not None
                    else "N/A",
                    "Match Percentage (%)": f"{pct_match:.2f}",
                }
            )
    except Exception as e:
        st.warning(f"NCBI BLAST query encounter: {e}")
    return seq_type, gc_content, gene_matches


def central_dogma_pipeline(seq: str, seq_type: str):
    bio_seq = Seq(seq)
    if seq_type == "DNA":
        return str(bio_seq.transcribe()), str(bio_seq.translate())
    elif seq_type == "RNA":
        return seq, str(bio_seq.translate())
    return None, seq


def find_orfs(dna_seq: str, min_protein_length: int = 15):
    orfs = []
    seq_obj = Seq(dna_seq)
    for strand, target_seq in [
        ("+", seq_obj),
        ("-", seq_obj.reverse_complement()),
    ]:
        target_str = str(target_seq)
        total_len = len(target_str)
        for frame in range(3):
            translated = target_seq[frame:].translate(to_stop=False)
            translated_str = str(translated)
            i = 0
            while i < len(translated_str):
                start_idx = translated_str.find("M", i)
                if start_idx == -1:
                    break
                stop_idx = translated_str.find("*", start_idx)
                if stop_idx == -1:
                    protein_len = len(translated_str) - start_idx
                    if protein_len >= min_protein_length:
                        start_nt = (start_idx * 3) + frame + 1
                        orfs.append(
                            {
                                "Strand": strand,
                                "Frame": f"Frame +{frame+1}"
                                if strand == "+"
                                else f"Frame -{frame+1}",
                                "Start (nt)": start_nt,
                                "End (nt)": total_len,
                                "Length (aa)": protein_len,
                                "Protein Sequence": translated_str[start_idx:],
                            }
                        )
                    break
                else:
                    protein_len = stop_idx - start_idx
                    if protein_len >= min_protein_length:
                        start_nt = (start_idx * 3) + frame + 1
                        end_nt = (stop_idx * 3) + frame + 3
                        orfs.append(
                            {
                                "Strand": strand,
                                "Frame": f"Frame +{frame+1}"
                                if strand == "+"
                                else f"Frame -{frame+1}",
                                "Start (nt)": start_nt,
                                "End (nt)": end_nt,
                                "Length (aa)": protein_len,
                                "Protein Sequence": translated_str[
                                    start_idx:stop_idx
                                ],
                            }
                        )
                    i = stop_idx + 1
    return orfs


def render_orf_diagram(orfs, total_seq_len):
    """Renders a clean graphical map of all 6 reading frames and detected ORFs."""
    frames = [
        "Frame +1",
        "Frame +2",
        "Frame +3",
        "Frame -1",
        "Frame -2",
        "Frame -3",
    ]
    frame_colors = {
        "Frame +1": "#38bdf8",
        "Frame +2": "#818cf8",
        "Frame +3": "#c084fc",
        "Frame -1": "#f43f5e",
        "Frame -2": "#fb923c",
        "Frame -3": "#facc15",
    }

    html = f"""
    <div style="background: rgba(10, 14, 23, 0.85); border: 1px solid rgba(255, 255, 255, 0.08); border-radius: 16px; padding: 22px; font-family: 'JetBrains Mono', monospace; box-shadow: 0 10px 35px rgba(0,0,0,0.5);">
        <div style="font-size: 13px; color: #94a3b8; margin-bottom: 16px; display: flex; justify-content: space-between; align-items: center;">
            <span style="color: #f1f5f9; font-weight: 600;"><b>6-Frame ORF Graphic Map</b></span>
            <span style="background: rgba(56, 189, 248, 0.1); color: #38bdf8; padding: 4px 10px; border-radius: 20px; font-size: 11px; border: 1px solid rgba(56, 189, 248, 0.2);">Total Length: {total_seq_len} nt</span>
        </div>
        <div style="position: relative; width: 100%;">
    """

    for f in frames:
        html += f"""
        <div style="display: flex; align-items: center; margin-bottom: 10px;">
            <div style="width: 85px; font-size: 11px; color: #cbd5e1; font-weight: 600;">{f}</div>
            <div style="position: relative; flex-grow: 1; height: 18px; background: rgba(255,255,255,0.02); border-radius: 6px; border: 1px solid rgba(255,255,255,0.05);">
        """

        frame_orfs = [o for o in orfs if o["Frame"] == f]
        for orf in frame_orfs:
            left_pct = max(0, min(100, (orf["Start (nt)"] / total_seq_len) * 100))
            width_pct = max(
                1,
                min(
                    100 - left_pct,
                    ((orf["End (nt)"] - orf["Start (nt)"]) / total_seq_len)
                    * 100,
                ),
            )
            color = frame_colors.get(f, "#38bdf8")

            html += f"""
            <div title="ORF: {orf['Start (nt)']} - {orf['End (nt)']} nt ({orf['Length (aa)']} aa)" 
                 style="position: absolute; left: {left_pct}%; width: {width_pct}%; height: 100%; background: {color}; opacity: 0.9; border-radius: 4px; cursor: pointer; box-shadow: 0 0 12px {color};">
            </div>
            """

        html += """
            </div>
        </div>
        """

    html += """
        </div>
        <div style="display: flex; gap: 18px; font-size: 11px; color: #94a3b8; margin-top: 18px; justify-content: flex-end;">
            <span style="display: flex; align-items: center; gap: 6px;"><span style="width:10px; height:10px; background:#38bdf8; display:inline-block; border-radius:3px; box-shadow: 0 0 8px #38bdf8;"></span> + Frames</span>
            <span style="display: flex; align-items: center; gap: 6px;"><span style="width:10px; height:10px; background:#f43f5e; display:inline-block; border-radius:3px; box-shadow: 0 0 8px #f43f5e;"></span> - Frames</span>
        </div>
    </div>
    """
    components.html(html, height=270, scrolling=False)


def fetch_pdb_similar(protein_seq: str):
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    query = {
        "query": {
            "type": "terminal",
            "service": "sequence",
            "parameters": {
                "evalue_cutoff": 1,
                "identity_cutoff": 0.3,
                "target": "pdb_protein_sequence",
                "value": protein_seq,
            },
        },
        "return_type": "polymer_entity",
        "request_options": {
            "paginate": {"start": 0, "rows": 5},
            "scoring_strategy": "sequence",
        },
    }
    response = requests.post(url, json=query)
    pdb_matches = []
    if response.status_code == 200:
        for item in response.json().get("result_set", []):
            full_id = item["identifier"]
            pdb_id = full_id.split("_")[0][:4]
            match_pct = item.get("score", 0) * 100
            if not any(d["PDB ID"] == pdb_id for d in pdb_matches):
                pdb_matches.append(
                    {
                        "PDB ID": pdb_id,
                        "Sequence Identity (%)": f"{match_pct:.2f}",
                    }
                )
    return pdb_matches


def analyze_amino_acids(protein_seq: str):
    if not protein_seq:
        return pd.DataFrame()
    total_aa = len(protein_seq)
    counts = Counter(protein_seq)
    valid_counts = {aa: count for aa, count in counts.items() if aa in AA_NAMES}
    data = [
        {
            "Amino Acid": AA_NAMES[aa],
            "Code": aa,
            "Count": count,
            "Percentage (%)": round((count / total_aa) * 100, 2),
        }
        for aa, count in Counter(valid_counts).most_common(10)
    ]
    return pd.DataFrame(data)


def color_protein_sequence_block(seq: str) -> str:
    bg_colors = {
        "A": "#38bdf8", "I": "#38bdf8", "L": "#38bdf8", "M": "#38bdf8", "F": "#38bdf8", "W": "#38bdf8", "V": "#38bdf8",
        "R": "#f43f5e", "K": "#f43f5e",
        "N": "#22c55e", "Q": "#22c55e",
        "D": "#c084fc", "E": "#c084fc",
        "C": "#fb923c",
        "G": "#facc15",
        "P": "#e879f9",
        "H": "#2dd4bf", "Y": "#2dd4bf",
        "S": "#4ade80", "T": "#4ade80",
    }
    styled_html = "<div style='font-family: \"JetBrains Mono\", monospace; font-size: 14px; word-break: break-all; line-height: 2.2; background: rgba(3, 7, 18, 0.75); border: 1px solid rgba(255, 255, 255, 0.08); padding: 18px; border-radius: 14px; letter-spacing: 1.5px; box-shadow: inset 0 2px 10px rgba(0,0,0,0.6);'>"
    for aa in seq:
        bg = bg_colors.get(aa, "#334155")
        text_color = "#030712" if aa in ["A", "I", "L", "M", "F", "W", "V", "N", "Q", "C", "G", "P", "H", "Y", "S", "T"] else "#ffffff"
        styled_html += f"<span style='background-color: {bg}; color: {text_color}; font-weight: 700; padding: 3px 6px; margin: 2px 1px; display: inline-block; text-align: center; border-radius: 4px; box-shadow: 0 2px 6px rgba(0,0,0,0.25);'>{aa}</span>"
    styled_html += "</div>"
    return styled_html


inject_custom_ui_theme()

st.markdown(
    """<style>
.header-container {
    display: flex;
    flex-direction: row;
    justify-content: space-between;
    align-items: flex-start;
    width: 100%;
    margin-bottom: 2rem;
    position: relative;
}
@keyframes titleTextGlow {
    0% { filter: drop-shadow(0px 0px 10px rgba(56, 189, 248, 0.4)); }
    50% { filter: drop-shadow(0px 0px 30px rgba(56, 189, 248, 0.95)); }
    100% { filter: drop-shadow(0px 0px 10px rgba(56, 189, 248, 0.4)); }
}
.title-glow-text {
    background: linear-gradient(135deg, #38bdf8 0%, #60a5fa 50%, #2563eb 100%);
    -webkit-background-clip: text;
    -webkit-text-fill-color: transparent;
    animation: titleTextGlow 4s ease-in-out infinite;
    display: inline-block;
}

.central-dogma-anim-vertical {
    position: absolute;
    top: -15px;
    right: 15px;
    width: 100px;  
    height: 320px; 
    overflow: visible;
    z-index: 100;
}

@keyframes spinDna {
    0% { transform: scaleX(1); }
    50% { transform: scaleX(-0.85); filter: drop-shadow(0 0 14px #38bdf8); }
    100% { transform: scaleX(1); }
}
.dna-layer {
    transform-origin: 60px 80px;
    animation: spinDna 5s linear infinite;
}

@keyframes swayRna {
    0%, 100% { transform: translateX(0px); }
    50% { transform: translateX(4px); }
}
.rna-strand {
    animation: swayRna 3s ease-in-out infinite;
}

@keyframes pulseProtein {
    0%, 100% { transform: scale(1) rotate(0deg); filter: drop-shadow(0 0 6px rgba(192, 132, 252, 0.6)); }
    50% { transform: scale(1.12) rotate(4deg); filter: drop-shadow(0 0 20px rgba(192, 132, 252, 1)); }
}
.protein-cluster {
    transform-origin: 60px 345px;
    animation: pulseProtein 3.5s ease-in-out infinite;
}

@keyframes processArrow {
    0%, 100% { opacity: 0.3; }
    50% { opacity: 1; filter: drop-shadow(0 0 10px #38bdf8); }
}
.process-arrow {
    animation: processArrow 2s infinite;
}
</style>
<div class="header-container">
<div style="flex: 1; min-width: 320px;">
    <h1 style='font-size: 3rem; font-weight: 800; margin: 0; color: #f8fafc; white-space: nowrap;'>
        ProtCraft <span class='title-glow-text'>Wizard</span> 🧙‍♂️
    </h1>
    <p style="color: #94a3b8; font-size: 1.05rem; margin-top: 8px; font-weight: 400;">Next-Gen Bioinformatics Sequence Pipeline & Structure Prediction Suite</p>
</div>
<div style="flex-shrink: 0; text-align: right;">
<svg class="central-dogma-anim-vertical" viewBox="0 0 120 400" fill="none" xmlns="http://www.w3.org/2000/svg">
<defs>
<linearGradient id="unifiedGrad" x1="0" y1="0" x2="0" y2="400" gradientUnits="userSpaceOnUse">
<stop offset="0%" stop-color="#38bdf8"/>
<stop offset="50%" stop-color="#818cf8"/>
<stop offset="100%" stop-color="#c084fc"/>
</linearGradient>
</defs>
<g class="dna-layer">
<path d="M 40 30 C 40 55, 80 65, 80 90 C 80 115, 40 125, 40 150" stroke="url(#unifiedGrad)" stroke-width="4.5" stroke-linecap="round"/>
<path d="M 80 30 C 80 55, 40 65, 40 90 C 40 115, 80 125, 80 150" stroke="url(#unifiedGrad)" stroke-width="4.5" stroke-linecap="round" opacity="0.85"/>
<line x1="42" y1="30" x2="78" y2="30" stroke="url(#unifiedGrad)" stroke-width="2.5"/>
<line x1="50" y1="45" x2="70" y2="45" stroke="url(#unifiedGrad)" stroke-width="2.5"/>
<line x1="65" y1="75" x2="55" y2="75" stroke="url(#unifiedGrad)" stroke-width="2.5"/>
<line x1="78" y1="90" x2="42" y2="90" stroke="url(#unifiedGrad)" stroke-width="2.5"/>
<line x1="65" y1="105" x2="55" y2="105" stroke="url(#unifiedGrad)" stroke-width="2.5"/>
<line x1="50" y1="135" x2="70" y2="135" stroke="url(#unifiedGrad)" stroke-width="2.5"/>
<line x1="42" y1="150" x2="78" y2="150" stroke="url(#unifiedGrad)" stroke-width="2.5"/>
</g>
<g class="process-arrow">
<line x1="60" y1="158" x2="60" y2="178" stroke="url(#unifiedGrad)" stroke-width="2" stroke-dasharray="4 2"/>
<polygon points="55,174 65,174 60,182" fill="url(#unifiedGrad)"/>
</g>
<g class="rna-strand">
<path d="M 50 195 C 75 220, 45 260, 70 285" stroke="url(#unifiedGrad)" stroke-width="4.5" fill="none" stroke-linecap="round"/>
</g>
<g class="process-arrow">
<line x1="60" y1="298" x2="60" y2="318" stroke="url(#unifiedGrad)" stroke-width="2" stroke-dasharray="4 2"/>
<polygon points="55,314 65,314 60,322" fill="url(#unifiedGrad)"/>
</g>
<g class="protein-cluster">
<path d="M 40 345 L 60 335 L 80 345 L 70 365 L 50 365 Z" stroke="url(#unifiedGrad)" stroke-width="3" fill="none" stroke-linejoin="round"/>
<path d="M 60 335 L 50 365" stroke="url(#unifiedGrad)" stroke-width="3" fill="none" opacity="0.6"/>
<circle cx="40" cy="345" r="9" fill="url(#unifiedGrad)"/>
<circle cx="60" cy="335" r="9" fill="url(#unifiedGrad)"/>
<circle cx="80" cy="345" r="9" fill="url(#unifiedGrad)"/>
<circle cx="70" cy="365" r="9" fill="url(#unifiedGrad)"/>
<circle cx="50" cy="365" r="9" fill="url(#unifiedGrad)"/>
</g>
</svg>
</div>
</div>""",
    unsafe_allow_html=True,
)

st.sidebar.header("Pipeline Settings")
user_email = st.sidebar.text_input(
    "NCBI Entrez Email", value="your.email@example.com"
)
Entrez.email = user_email

st.header("Input Sequence")
st.markdown(
    "<p style='font-size: 0.95rem; font-weight: 500; margin-bottom: 0.6rem; color: #cbd5e1;'>Choose Input Method:</p>",
    unsafe_allow_html=True,
)

input_option = st.radio(
    "Choose Input Method:",
    options=["Raw Sequence / Accession ID", "Upload FASTA File"],
    horizontal=True,
    label_visibility="collapsed",
)

raw_input = ""
uploaded_file = None

if input_option == "Raw Sequence / Accession ID":
    raw_input = st.text_area(
        "Enter Sequence or Accession ID:",
        value="",
        placeholder="e.g., NM_000518 or ATGCG...",
        height=140,
    )
else:
    uploaded_file = st.file_uploader(
        "Upload FASTA file", type=["fasta", "fas", "fa"]
    )

st.markdown("<div style='margin-top: 1rem;'></div>", unsafe_allow_html=True)

if st.button("Run Pipeline", type="primary"):
    with st.spinner("Processing input sequence..."):
        sequence = fetch_sequence(raw_input, uploaded_file)

    if not sequence:
        st.error(
            "No valid sequence input detected. Please provide a sequence, accession ID, or FASTA file."
        )
    else:
        st.success("Sequence successfully loaded!")

        st.header("Identification & BLAST Analysis")
        with st.spinner("Analyzing sequence type and querying BLAST..."):
            seq_type, gc_content, gene_matches = identify_sequence(sequence)

        if seq_type:
            col1, col2, col3 = st.columns(3)
            col1.metric("Sequence Type", seq_type)
            col2.metric(
                "GC Content",
                f"{gc_content:.2f}%" if gc_content is not None else "N/A",
            )
            col3.metric("Sequence Length", f"{len(sequence)} bp/aa")

            st.subheader("Top 5 Gene Matches (NCBI BLAST)")
            if gene_matches:
                df_matches = pd.DataFrame(gene_matches)
                st.dataframe(df_matches, use_container_width=True)
            else:
                st.info("No significant BLAST hits found.")

            st.header("Transcription & Translation")
            transcript, protein_seq = central_dogma_pipeline(
                sequence, seq_type
            )

            if transcript:
                with st.expander("View mRNA Transcript"):
                    st.text_area("RNA Sequence", transcript, height=100)

            with st.expander(
                "View Translated Protein Sequence", expanded=True
            ):
                st.markdown(
                    color_protein_sequence_block(protein_seq),
                    unsafe_allow_html=True,
                )
                st.caption(
                    "🟦 Hydrophobic | 🟥 Basic | 🟩 Polar | 🟪 Acidic | 🟧 Glycine | 🟨 Proline"
                )

            st.header("Open Reading Frame (ORF) Diagram Map & Sequences")
            if seq_type in ["DNA", "RNA"]:
                dna_for_orf = sequence.replace("U", "T")

                with st.spinner("Scanning 6 reading frames for ORFs..."):
                    orf_list = find_orfs(dna_for_orf, min_protein_length=15)

                if orf_list:
                    render_orf_diagram(orf_list, len(dna_for_orf))

                    st.subheader("Detected ORFs Sequence Data")
                    for idx, orf in enumerate(orf_list):
                        with st.expander(
                            f"ORF #{idx+1} | Strand: {orf['Strand']} | Frame: {orf['Frame']} | Coordinates: {orf['Start (nt)']} - {orf['End (nt)']} nt | Length: {orf['Length (aa)']} aa"
                        ):
                            st.markdown(
                                f"**Protein Sequence:**", unsafe_allow_html=True
                            )
                            st.markdown(
                                color_protein_sequence_block(
                                    orf["Protein Sequence"]
                                ),
                                unsafe_allow_html=True,
                            )
                            st.text_area(
                                f"Raw Protein Sequence #{idx+1}",
                                orf["Protein Sequence"],
                                height=80,
                                key=f"orf_seq_{idx}",
                            )
                else:
                    st.info(
                        "No Open Reading Frames found meeting the minimum length criteria."
                    )
            else:
                st.info(
                    "ORF scanning is available for DNA and RNA input sequence types."
                )

            st.header("Protein Analysis")
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
                    st.dataframe(
                        pd.DataFrame(matches), use_container_width=True
                    )
                else:
                    st.write("No significant RCSB PDB matches found.")

            # --- SWISS-MODEL URL-BASED PREDICTION INTEGRATION ---
            st.header("Protein Structure Prediction (SWISS-MODEL)")
            
            if protein_seq:
                encoded_seq = urllib.parse.quote(protein_seq)
                swiss_url = f"https://swissmodel.expasy.org/interactive?target={encoded_seq}"
                
                st.markdown(
                    f"""
                    <div style="background: rgba(10, 14, 23, 0.85); border: 1px solid rgba(56, 189, 248, 0.25); border-radius: 16px; padding: 24px; box-shadow: 0 10px 35px rgba(0,0,0,0.5); text-align: center;">
                        <h3 style="color: #38bdf8; margin-top: 0; margin-bottom: 10px;">Ready for Homology Modeling</h3>
                        <p style="color: #cbd5e1; font-size: 0.95rem; margin-bottom: 20px; line-height: 1.5;">
                            Click below to seamlessly transfer your translated protein sequence ({len(protein_seq)} aa) directly into SWISS-MODEL for automated template search and 3D structure prediction.
                        </p>
                        <a href="{swiss_url}" target="_blank" style="display: inline-block; background: linear-gradient(135deg, #0ea5e9 0%, #2563eb 100%); color: #ffffff; padding: 0.8rem 2.2rem; border-radius: 12px; font-weight: 700; text-decoration: none; box-shadow: 0 4px 20px rgba(14, 165, 233, 0.45); transition: all 0.3s ease;">
                            🚀 Open in SWISS-MODEL
                        </a>
                    </div>
                    """,
                    unsafe_allow_html=True
                )
            else:
                st.warning("No protein sequence available for modeling.")
