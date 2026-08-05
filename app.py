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
import streamlit.components.v1 as components
from textwrap import dedent

st.set_page_config(page_title="Bioinformatics Sequence Pipeline", layout="wide")


AA_NAMES = {
    'A': 'Alanine', 'C': 'Cysteine', 'D': 'Aspartic Acid', 'E': 'Glutamic Acid',
    'F': 'Phenylalanine', 'G': 'Glycine', 'H': 'Histidine', 'I': 'Isoleucine',
    'K': 'Lysine', 'L': 'Leucine', 'M': 'Methionine', 'N': 'Asparagine',
    'P': 'Proline', 'Q': 'Glutamine', 'R': 'Arginine', 'S': 'Serine',
    'T': 'Threonine', 'V': 'Valine', 'W': 'Tryptophan', 'Y': 'Tyrosine'
}

def inject_custom_ui_theme():
    """Injects dynamic breathing background, glassmorphism UI, and equal-sized side-by-side radio cards."""
    css = """
    <style>
    /* 1. Base Deep Black Background */
    @keyframes breathingGradient {
        0% { background-position: 0% 50%; }
        50% { background-position: 100% 50%; }
        100% { background-position: 0% 50%; }
    }

    .stApp {
        background: linear-gradient(-45deg, #000000, #0a0a0c, #050811, #000000);
        background-size: 400% 400%;
        animation: breathingGradient 18s ease infinite;
        color: #f8fafc;
    }

    /* 2. Dynamic Blurry Light Zone */
    @keyframes lightSweep {
        0% { opacity: 0.2; transform: scale(1); }
        50% { opacity: 0.6; transform: scale(1.15); }
        100% { opacity: 0.2; transform: scale(1); }
    }

    .stApp::after {
        content: "";
        position: fixed;
        top: -20%; left: -20%; width: 140%; height: 140%;
        background: linear-gradient(
            125deg, 
            transparent 15%, 
            rgba(56, 189, 248, 0.12) 35%, 
            rgba(37, 99, 235, 0.22) 50%, 
            rgba(56, 189, 248, 0.12) 65%, 
            transparent 85%
        );
        background-size: 200% 200%;
        filter: blur(100px);
        pointer-events: none;
        z-index: 0;
        animation: lightSweep 14s ease-in-out infinite alternate;
    }

    /* Elevate content above background light layer */
    .block-container {
        position: relative;
        z-index: 1;
    }

    /* 3. Glassmorphism Containers */
    div[data-testid="stExpander"], div[data-testid="stMetric"] {
        background: rgba(255, 255, 255, 0.02) !important;
        backdrop-filter: blur(12px);
        -webkit-backdrop-filter: blur(12px);
        border: 1px solid rgba(255, 255, 255, 0.08) !important;
        border-radius: 12px !important;
        padding: 10px !important;
        position: relative;
        z-index: 2;
    }

    /* Headers & Text */
    h1, h2, h3 {
        color: #f1f5f9 !important;
        font-weight: 700 !important;
        letter-spacing: -0.02em;
    }

    /* Primary Interactive Buttons */
    .stButton > button {
        background: linear-gradient(135deg, #0284c7 0%, #2563eb 100%) !important;
        color: #ffffff !important;
        border: none !important;
        border-radius: 8px !important;
        padding: 0.6rem 1.5rem !important;
        font-weight: 600 !important;
        transition: all 0.3s ease !important;
        box-shadow: 0 4px 15px rgba(37, 99, 235, 0.3);
        position: relative;
        z-index: 2;
    }

    .stButton > button:hover {
        transform: translateY(-2px);
        box-shadow: 0 6px 20px rgba(37, 99, 235, 0.5);
    }

    /* Sidebar Styling */
    section[data-testid="stSidebar"] {
        background-color: rgba(10, 10, 12, 0.85) !important;
        backdrop-filter: blur(16px);
        border-right: 1px solid rgba(255, 255, 255, 0.08);
        z-index: 10;
    }
    
    header[data-testid="stHeader"] {
        background: transparent !important;
    }

    /* 4. Side-by-Side Selectable Glass Cards (Equal Dimensions) */
    div[data-testid="stRadio"] div[role="radiogroup"] {
        display: flex !important;
        flex-direction: row !important;
        gap: 16px !important;
        width: 100% !important;
        align-items: stretch !important;
    }

    div[data-testid="stRadio"] div[role="radiogroup"] label {
        flex: 1 1 0px !important;
        min-height: 72px !important;
        background: rgba(255, 255, 255, 0.03) !important;
        backdrop-filter: blur(10px) !important;
        border: 1px solid rgba(255, 255, 255, 0.1) !important;
        border-radius: 10px !important;
        padding: 12px 20px !important;
        transition: all 0.3s ease !important;
        cursor: pointer !important;
        display: flex !important;
        align-items: center !important;
        justify-content: flex-start !important;
    }

    /* Card Hover State */
    div[data-testid="stRadio"] div[role="radiogroup"] label:hover {
        background: rgba(56, 189, 248, 0.08) !important;
        border-color: rgba(56, 189, 248, 0.4) !important;
        transform: translateY(-1px);
    }

    /* Active Selected Card Highlight */
    div[data-testid="stRadio"] div[role="radiogroup"] label:has(input:checked) {
        background: rgba(56, 189, 248, 0.12) !important;
        border-color: #38bdf8 !important;
        box-shadow: 0 0 15px rgba(56, 189, 248, 0.2) !important;
    }
    </style>
    """
    st.markdown(css, unsafe_allow_html=True)

def render_protein_3d_viewer(
    pdb_input: str, is_pdb_id: bool = True, height: int = 480
):
    """Renders a touch-optimized, mobile-friendly 3D protein viewer using 3Dmol.js.

    Parameters:
    - pdb_input: PDB ID string (e.g. '1A2C') OR raw PDB file contents.
    - is_pdb_id: True if passing a 4-letter PDB ID, False if passing raw PDB string (e.g., from ESMFold).
    - height: Height of the canvas in pixels.
    """
    # Sanitize string input for JavaScript injection if raw PDB data is passed
    if is_pdb_id:
        fetch_js = f"v.addModelAsPdbId('{pdb_input.strip()}');"
    else:
        escaped_pdb = (
            pdb_input.replace("\\", "\\\\")
            .replace("`", "\\`")
            .replace("${", "\\${")
        )
        fetch_js = f"v.addModel(`{escaped_pdb}`, 'pdb');"

    html_code = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <meta name="viewport" content="width=device-width, initial-scale=1.0, maximum-scale=1.0, user-scalable=no" />
        <script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.0.4/3Dmol-min.js"></script>
        <style>
            * {{
                box-sizing: border-box;
                margin: 0;
                padding: 0;
            }}
            body, html {{
                width: 100%;
                height: 100%;
                overflow: hidden;
                background-color: transparent;
                font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
            }}
            
            /* Glassmorphism Outer Wrapper */
            .viewer-wrapper {{
                position: relative;
                width: 100%;
                height: {height}px;
                background: rgba(10, 10, 12, 0.6);
                backdrop-filter: blur(12px);
                border: 1px solid rgba(255, 255, 255, 0.1);
                border-radius: 14px;
                overflow: hidden;
            }}

            /* Canvas Container - touch-action: none prevents parent page scrolling during 3D rotation */
            #viewport {{
                width: 100%;
                height: 100%;
                touch-action: none;
            }}

            /* Mobile Floating Quick Action Toolbar */
            .controls-bar {{
                position: absolute;
                bottom: 12px;
                left: 50%;
                transform: translateX(-50%);
                display: flex;
                gap: 6px;
                background: rgba(0, 0, 0, 0.75);
                backdrop-filter: blur(8px);
                padding: 6px 12px;
                border-radius: 30px;
                border: 1px solid rgba(255, 255, 255, 0.15);
                z-index: 10;
                width: max-content;
                max-width: 92%;
                overflow-x: auto;
            }}

            .control-btn {{
                background: rgba(255, 255, 255, 0.08);
                color: #e2e8f0;
                border: 1px solid rgba(255, 255, 255, 0.12);
                border-radius: 20px;
                padding: 6px 12px;
                font-size: 11px;
                font-weight: 600;
                cursor: pointer;
                white-space: nowrap;
                transition: all 0.2s ease;
                -webkit-tap-highlight-color: transparent;
            }}

            .control-btn:active, .control-btn.active {{
                background: rgba(56, 189, 248, 0.25);
                border-color: #38bdf8;
                color: #38bdf8;
            }}

            /* Mobile Gesture Hint Badge */
            .gesture-hint {{
                position: absolute;
                top: 10px;
                right: 12px;
                background: rgba(0, 0, 0, 0.5);
                color: rgba(255, 255, 255, 0.5);
                font-size: 10px;
                padding: 4px 8px;
                border-radius: 6px;
                pointer-events: none;
                backdrop-filter: blur(4px);
            }}
        </style>
    </head>
    <body>
        <div class="viewer-wrapper">
            <div class="gesture-hint">👆 1-Finger Rotate | 🤏 Pinch Zoom</div>
            <div id="viewport"></div>
            
            <div class="controls-bar">
                <button class="control-btn active" onclick="setStyle('cartoon')">Cartoon</button>
                <button class="control-btn" onclick="setStyle('stick')">Sticks</button>
                <button class="control-btn" onclick="setStyle('sphere')">Sphere</button>
                <button class="control-btn" onclick="toggleSurface()">Surface</button>
                <button class="control-btn" onclick="resetView()">Reset</button>
            </div>
        </div>

        <script>
            let viewer = null;
            let currentModel = null;
            let surfaceObj = null;
            let currentStyle = 'cartoon';

            document.addEventListener("DOMContentLoaded", function() {{
                let element = document.getElementById('viewport');
                let config = {{ backgroundColor: '0x000000', backgroundAlpha: 0.0 }};
                
                viewer = $3Dmol.createViewer(element, config);
                let v = viewer;

                {fetch_js}
                
                currentModel = v.getModel();
                setStyle('cartoon');
                v.zoomTo();
                v.render();
            }});

            function setStyle(type) {{
                if (!viewer) return;
                currentStyle = type;
                
                // Update active state on buttons
                document.querySelectorAll('.control-btn').forEach(btn => {{
                    if(['Cartoon', 'Sticks', 'Sphere'].includes(btn.innerText)) {{
                        btn.classList.remove('active');
                    }}
                }});

                viewer.setStyle({{}}, {{}}); // Clear existing styles

                if (type === 'cartoon') {{
                    viewer.setStyle({{}}, {{cartoon: {{color: 'spectrum'}}}});
                }} else if (type === 'stick') {{
                    viewer.setStyle({{}}, {{stick: {{colorscheme: 'amino'}}}});
                }} else if (type === 'sphere') {{
                    viewer.setStyle({{}}, {{sphere: {{scale: 0.28, colorscheme: 'spectrum'}}}});
                }}

                viewer.render();
            }}

            function toggleSurface() {{
                if (!viewer) return;
                let btn = event.target;
                
                if (surfaceObj) {{
                    viewer.removeSurface(surfaceObj);
                    surfaceObj = null;
                    btn.classList.remove('active');
                }} else {{
                    btn.classList.add('active');
                    surfaceObj = viewer.addSurface($3Dmol.SurfaceType.MS, {{
                        opacity: 0.6,
                        color: 'white'
                    }});
                }}
                viewer.render();
            }}

            function resetView() {{
                if (!viewer) return;
                viewer.zoomTo();
                viewer.render();
            }}
        </script>
    </body>
    </html>
    """
    components.html(html_code, height=height + 10, scrolling=False)

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

def predict_structure_colabfold(
    sequence: str, api_url: str = "YOUR_COLABFOLD_API_ENDPOINT"
):
    """Sends sequence to ColabFold GPU backend and returns PDB string."""
    try:
        payload = {"sequence": sequence.strip()}
        response = requests.post(
            f"{api_url}/predict", json=payload, timeout=300
        )

        if response.status_code == 200:
            res_json = response.json()
            return res_json.get("pdb_data")
        else:
            st.error(f"ColabFold server error (Status {response.status_code})")
            return None
    except Exception as e:
        st.error(f"Failed to connect to ColabFold service: {str(e)}")
        return None

def color_protein_sequence_block(seq: str) -> str:
    """Renders sequence with solid colored background blocks matching standard MSA/Clustal tools."""
    bg_colors = {
        'A': '#80a0f0', 'I': '#80a0f0', 'L': '#80a0f0', 'M': '#80a0f0', 'F': '#80a0f0', 'W': '#80a0f0', 'V': '#80a0f0',
        'R': '#f01505', 'K': '#f01505',
        'N': '#00ff00', 'Q': '#00ff00',
        'D': '#c000c0', 'E': '#c000c0',
        'C': '#f08080',
        'G': '#f09040',
        'P': '#ffff00',
        'H': '#15a4a4', 'Y': '#15a4a4',
        'S': '#15a400', 'T': '#15a400'
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

st.set_page_config(
    page_title="Bioinformatics Sequence Pipeline", layout="wide"
)
inject_custom_ui_theme()

st.markdown("""<style>
.header-container {
display: flex;
flex-direction: row;
justify-content: space-between;
align-items: flex-start;
width: 100%;
margin-bottom: 1.5rem;
position: relative;
}
@keyframes titleTextGlow {
0% { filter: drop-shadow(0px 0px 8px rgba(56, 189, 248, 0.35)); }
50% { filter: drop-shadow(0px 0px 24px rgba(56, 189, 248, 0.85)); }
100% { filter: drop-shadow(0px 0px 8px rgba(56, 189, 248, 0.35)); }
}
.title-glow-text {
background: linear-gradient(90deg, #38bdf8 0%, #60a5fa 50%, #3b82f6 100%);
-webkit-background-clip: text;
-webkit-text-fill-color: transparent;
animation: titleTextGlow 4s ease-in-out infinite;
display: inline-block;
}

/* NEW SLIGHTLY SMALLER SVG ANIMATION STYLING */
.central-dogma-anim-vertical {
position: absolute;
top: -10px;
right: 15px;
width: 100px;  
height: 320px; 
overflow: visible;
z-index: 100;
}

@keyframes spinDna {
0% { transform: scaleX(1); }
50% { transform: scaleX(-0.85); filter: drop-shadow(0 0 10px #38bdf8); }
100% { transform: scaleX(1); }
}
.dna-layer {
transform-origin: 60px 80px;
animation: spinDna 5s linear infinite;
}

@keyframes swayRna {
0%, 100% { transform: translateX(0px); }
50% { transform: translateX(3px); }
}
.rna-strand {
animation: swayRna 3s ease-in-out infinite;
}

@keyframes pulseProtein {
0%, 100% { transform: scale(1) rotate(0deg); filter: drop-shadow(0 0 5px rgba(192, 132, 252, 0.5)); }
50% { transform: scale(1.08) rotate(3deg); filter: drop-shadow(0 0 15px rgba(192, 132, 252, 1)); }
}
.protein-cluster {
transform-origin: 60px 345px;
animation: pulseProtein 3.5s ease-in-out infinite;
}

@keyframes processArrow {
0%, 100% { opacity: 0.3; }
50% { opacity: 1; filter: drop-shadow(0 0 6px white); }
}
.process-arrow {
animation: processArrow 2s infinite;
}
</style>
<div class="header-container">
<div style="flex: 1; min-width: 320px;">
<h1 style='font-size: 2.8rem; font-weight: 800; margin: 0; color: #f8fafc; white-space: nowrap;'>
Welcome to <span class='title-glow-text'>ProtCraft Wizard</span> 🧙‍♂️
</h1>
</div>
<div style="flex-shrink: 0; text-align: right;">
<svg class="central-dogma-anim-vertical" viewBox="0 0 120 400" fill="none" xmlns="http://www.w3.org/2000/svg">
<defs>
<!-- Single Unified Gradient for the entire vertical structure -->
<linearGradient id="unifiedGrad" x1="0" y1="0" x2="0" y2="400" gradientUnits="userSpaceOnUse">
<stop offset="0%" stop-color="#38bdf8"/>
<stop offset="50%" stop-color="#818cf8"/>
<stop offset="100%" stop-color="#c084fc"/>
</linearGradient>
</defs>

<!-- 1. DISTINCT DNA DOUBLE HELIX -->
<g class="dna-layer">
<!-- Intersecting Double Helix Paths -->
<path d="M 40 30 C 40 55, 80 65, 80 90 C 80 115, 40 125, 40 150" stroke="url(#unifiedGrad)" stroke-width="4.5" stroke-linecap="round"/>
<path d="M 80 30 C 80 55, 40 65, 40 90 C 40 115, 80 125, 80 150" stroke="url(#unifiedGrad)" stroke-width="4.5" stroke-linecap="round" opacity="0.85"/>

<!-- Connecting Base Pairs (Rungs) -->
<line x1="42" y1="30" x2="78" y2="30" stroke="url(#unifiedGrad)" stroke-width="2.5"/>
<line x1="50" y1="45" x2="70" y2="45" stroke="url(#unifiedGrad)" stroke-width="2.5"/>
<line x1="65" y1="75" x2="55" y2="75" stroke="url(#unifiedGrad)" stroke-width="2.5"/>
<line x1="78" y1="90" x2="42" y2="90" stroke="url(#unifiedGrad)" stroke-width="2.5"/>
<line x1="65" y1="105" x2="55" y2="105" stroke="url(#unifiedGrad)" stroke-width="2.5"/>
<line x1="50" y1="135" x2="70" y2="135" stroke="url(#unifiedGrad)" stroke-width="2.5"/>
<line x1="42" y1="150" x2="78" y2="150" stroke="url(#unifiedGrad)" stroke-width="2.5"/>
</g>

<!-- TRANSCRIPTION ARROW -->
<g class="process-arrow">
<line x1="60" y1="158" x2="60" y2="178" stroke="url(#unifiedGrad)" stroke-width="2" stroke-dasharray="4 2"/>
<polygon points="55,174 65,174 60,182" fill="url(#unifiedGrad)"/>
</g>

<!-- 2. DISTINCT RNA SINGLE STRAND WITH SOLID BACKBONE -->
<g class="rna-strand">
<!-- Solid Wavy Single Backbone -->
<path d="M 50 195 C 75 220, 45 260, 70 285" stroke="url(#unifiedGrad)" stroke-width="4.5" fill="none" stroke-linecap="round"/>
</g>

<!-- TRANSLATION ARROW -->
<g class="process-arrow">
<line x1="60" y1="298" x2="60" y2="318" stroke="url(#unifiedGrad)" stroke-width="2" stroke-dasharray="4 2"/>
<polygon points="55,314 65,314 60,322" fill="url(#unifiedGrad)"/>
</g>

<!-- 3. DISTINCT FOLDED PROTEIN POLYPEPTIDE CHAIN -->
<g class="protein-cluster">
<!-- Tangled Polypeptide Backbone Connectors -->
<path d="M 40 345 L 60 335 L 80 345 L 70 365 L 50 365 Z" stroke="url(#unifiedGrad)" stroke-width="3" fill="none" stroke-linejoin="round"/>
<path d="M 60 335 L 50 365" stroke="url(#unifiedGrad)" stroke-width="3" fill="none" opacity="0.6"/>

<!-- Amino Acid Spheres -->
<circle cx="40" cy="345" r="9" fill="url(#unifiedGrad)"/>
<circle cx="60" cy="335" r="9" fill="url(#unifiedGrad)"/>
<circle cx="80" cy="345" r="9" fill="url(#unifiedGrad)"/>
<circle cx="70" cy="365" r="9" fill="url(#unifiedGrad)"/>
<circle cx="50" cy="365" r="9" fill="url(#unifiedGrad)"/>
</g>
</svg>
</div>
</div>""", unsafe_allow_html=True)
st.sidebar.header("Settings")
user_email = st.sidebar.text_input(
    "NCBI Entrez Email", value="your.email@example.com"
)
Entrez.email = user_email


# --- SECTION 1: INPUT SEQUENCE ---
st.header("Input Sequence")

# Normal plain text label
st.markdown(
    "<p style='font-size: 1rem; font-weight: 500; margin-bottom: 0.5rem; color: #e2e8f0;'>Choose Input Method:</p>",
    unsafe_allow_html=True,
)

# Radio options with hidden internal label
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

            
            st.header("Protein Structure Visualization")

            engine = st.radio(
                "Select Prediction Engine:",
                options=[
                    "ESMFold (Ultra-Fast | <400 aa)",
                    "ColabFold / AlphaFold2 (High-Accuracy | Up to 1200 aa)",
                ],
                horizontal=True,
                key="structure_prediction_engine",
            )

            pdb_data = None

            if "ESMFold" in engine:
                if len(protein_seq) > 400:
                    st.warning(
                        "Protein length exceeds 400 amino acids. ESMFold API predictions are restricted to shorter sequences."
                    )
                else:
                    with st.spinner("Predicting 3D structure using ESMFold..."):
                        pdb_data = predict_structure_esm(protein_seq)
            else:
                with st.spinner(
                    "Generating MSAs and predicting structure using ColabFold (1–3 mins)..."
                ):
                    # Replace URL below with your active ngrok/Modal endpoint
                    COLABFOLD_API = (
                        "https://banshee-remedy-oblong.ngrok-free.dev/"
                    )
                    pdb_data = predict_structure_colabfold(
                        protein_seq, api_url=COLABFOLD_API
                    )

            if pdb_data:
                st.success("3D Structure predicted successfully!")

                st.subheader("Interactive 3D Structure Viewer")
                render_protein_3d_viewer(
                    pdb_input=pdb_data, is_pdb_id=False, height=500
                )

                st.download_button(
                    label="Download PDB File",
                    data=pdb_data,
                    file_name="predicted_structure.pdb",
                    mime="chemical/x-pdb",
                )
            elif "ESMFold" in engine and len(protein_seq) <= 400:
                st.error("Failed to predict 3D structure using ESMFold API.")
            elif "ColabFold" in engine:
                st.error(
                    "Failed to predict 3D structure using ColabFold API. Check backend connection."
                )
