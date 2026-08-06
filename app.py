import io
import time
import xml.etree.ElementTree as ET

from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
import pandas as pd
import requests
import streamlit as st
import streamlit.components.v1 as components

st.set_page_config(page_title="ProtCraft Wizard", layout="wide", initial_sidebar_state="collapsed")

def inject_theme():
    st.markdown("""
    <style>
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;600;800&family=JetBrains+Mono:wght@400;700&display=swap');
    
    h1, h2, h3, p, label, span { font-family: 'Inter', sans-serif; }
    code, pre { font-family: 'JetBrains Mono', monospace !important; }
    
    header[data-testid="stHeader"], [data-testid="collapsedControl"] { display: none !important; }
    
    @keyframes titleTextGlow {
        0% { filter: drop-shadow(0px 0px 8px rgba(56, 189, 248, 0.4)); }
        50% { filter: drop-shadow(0px 0px 25px rgba(56, 189, 248, 0.95)); }
        100% { filter: drop-shadow(0px 0px 8px rgba(56, 189, 248, 0.4)); }
    }
    .title-glow { color: #38bdf8; animation: titleTextGlow 3.5s ease-in-out infinite; display: inline-block; }

    @keyframes orbit1 {
        0% { transform: rotate(0deg) translateX(15vw) scale(1); }
        50% { transform: rotate(180deg) translateX(15vw) scale(1.2); }
        100% { transform: rotate(360deg) translateX(15vw) scale(1); }
    }
    @keyframes orbit2 {
        0% { transform: rotate(180deg) translateX(15vw) scale(1); }
        50% { transform: rotate(360deg) translateX(15vw) scale(1.2); }
        100% { transform: rotate(540deg) translateX(15vw) scale(1); }
    }

    .stApp { background: #020617; color: #f1f5f9; overflow-x: hidden; }
    
    .stApp::before {
        content: ""; position: fixed; top: 20%; left: 25%; width: 50vw; height: 50vw;
        background: radial-gradient(circle, rgba(56,189,248,0.25), transparent 60%);
        filter: blur(80px); animation: orbit1 12s linear infinite; z-index: 0; pointer-events: none;
    }
    .stApp::after {
        content: ""; position: fixed; top: 20%; left: 25%; width: 50vw; height: 50vw;
        background: radial-gradient(circle, rgba(168,85,247,0.25), transparent 60%);
        filter: blur(80px); animation: orbit2 12s linear infinite; z-index: 0; pointer-events: none;
    }

    .block-container { position: relative; z-index: 1; padding-top: 3rem; }

    div[data-testid="stExpander"], div[data-testid="stMetric"] {
        background: rgba(255, 255, 255, 0.03) !important; backdrop-filter: blur(20px);
        border: 1px solid rgba(255, 255, 255, 0.08) !important; border-radius: 14px !important;
    }

    div[data-testid="stRadio"] div[role="radiogroup"] { display: flex; gap: 16px; }
    div[data-testid="stRadio"] div[role="radiogroup"] label {
        flex: 1; background: rgba(255,255,255,0.03); border: 1px solid rgba(255,255,255,0.1);
        border-radius: 12px; padding: 0 20px; cursor: pointer; transition: all 0.3s;
        min-height: 75px; display: flex; align-items: center; 
    }
    div[data-testid="stRadio"] div[role="radiogroup"] label:has(input:checked) {
        background: rgba(56, 189, 248, 0.15); border-color: #38bdf8; box-shadow: 0 0 20px rgba(56,189,248,0.2);
    }

    button[kind="primary"] {
        background: linear-gradient(135deg, #0ea5e9 0%, #2563eb 100%) !important;
        color: #fff !important; border-radius: 12px !important; padding: 0.6rem 2rem !important;
        font-weight: 600 !important; border: none !important; box-shadow: 0 4px 20px rgba(14, 165, 233, 0.4);
    }
    </style>
    """, unsafe_allow_html=True)

def render_protein_3d_viewer(pdb_data: str, height: int = 500):
    escaped_pdb = pdb_data.replace("\\", "\\\\").replace("`", "\\`").replace("${", "\\${").replace("\n", "\\n")
    html_code = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.0.4/3Dmol-min.js"></script>
        <style>
            body, html {{ width: 100%; height: 100%; margin: 0; overflow: hidden; background-color: transparent; }}
            .viewer-wrapper {{ position: relative; width: 100%; height: {height}px; background: rgba(10, 10, 12, 0.6); backdrop-filter: blur(12px); border: 1px solid rgba(255, 255, 255, 0.1); border-radius: 14px; overflow: hidden; }}
            #viewport {{ width: 100%; height: 100%; }}
            .controls {{ position: absolute; bottom: 12px; left: 50%; transform: translateX(-50%); display: flex; gap: 8px; background: rgba(0, 0, 0, 0.8); padding: 8px 16px; border-radius: 30px; border: 1px solid rgba(255,255,255,0.15); z-index: 10; }}
            .btn {{ background: rgba(255,255,255,0.1); color: #fff; border: none; border-radius: 20px; padding: 6px 14px; cursor: pointer; font-family: sans-serif; font-size: 12px; font-weight: 600; transition: 0.3s; }}
            .btn:hover {{ background: #38bdf8; }}
        </style>
    </head>
    <body>
        <div class="viewer-wrapper">
            <div id="viewport"></div>
            <div class="controls">
                <button class="btn" onclick="v.setStyle({{}},{{cartoon:{{color:'spectrum'}}}}); v.render();">Cartoon</button>
                <button class="btn" onclick="v.setStyle({{}},{{stick:{{}}}}); v.render();">Sticks</button>
                <button class="btn" onclick="v.setStyle({{}},{{sphere:{{scale: 0.3, color:'spectrum'}}}}); v.render();">Sphere</button>
            </div>
        </div>
        <script>
            let v = null;
            document.addEventListener("DOMContentLoaded", function() {{
                v = $3Dmol.createViewer(document.getElementById('viewport'), {{backgroundColor: '0x000000', backgroundAlpha: 0.0}});
                v.addModel(`{escaped_pdb}`, "pdb");
                v.setStyle({{}}, {{cartoon: {{color: 'spectrum'}}}});
                v.zoomTo(); v.render();
            }});
        </script>
    </body>
    </html>
    """
    components.html(html_code, height=height + 15, scrolling=False)


inject_theme()

st.sidebar.header("API Configurations")
Entrez.email = st.sidebar.text_input("NCBI Email", value="protcraft@example.com")
swiss_token = st.sidebar.text_input("SWISS-MODEL API Token", type="password", help="Required for background prediction.")

st.markdown("<h1 style='font-size: 3rem; font-weight: 800; margin-bottom: 0;'>ProtCraft <span class='title-glow'>Wizard</span> 🧙‍♂️</h1><p style='color: #94a3b8; margin-bottom: 2rem;'>Next-Gen Bioinformatics Pipeline</p>", unsafe_allow_html=True)

input_option = st.radio("Input Method", ["Raw Sequence / Accession ID", "Upload FASTA File"], horizontal=True, label_visibility="collapsed")
raw_input = st.text_area("Input", placeholder="Enter sequence (e.g., ATGC...) or ID (e.g., NM_000518)...", height=120) if "Raw" in input_option else None
uploaded_file = st.file_uploader("Upload FASTA", type=["fasta", "fa"]) if "Upload" in input_option else None

if st.button("Run Pipeline", type="primary"):
    seq = ""
    if uploaded_file:
        seq = str(SeqIO.read(io.StringIO(uploaded_file.getvalue().decode("utf-8")), "fasta").seq).upper()
    elif raw_input:
        q = raw_input.strip()
        if any(c.isdigit() for c in q) and len(q) < 20:
            for db in ["nucleotide", "protein"]:
                try:
                    seq = str(SeqIO.read(Entrez.efetch(db=db, id=q, rettype="fasta", retmode="text"), "fasta").seq).upper()
                    break
                except: continue
        if not seq: seq = q.replace(" ", "").replace("\n", "").upper()

    if not seq:
        st.error("Please provide a valid sequence.")
    else:
        seq_set = set(seq)
        seq_type = "DNA" if seq_set.issubset(set("ACGTN")) else ("RNA" if seq_set.issubset(set("ACGUN")) and "U" in seq_set else "Protein" if seq_set.issubset(set("ACDEFGHIKLMNPQRSTVWY")) else None)
        
        if not seq_type:
            st.error("Invalid sequence characters detected.")
        else:
            gc = gc_fraction(seq) * 100 if seq_type in ["DNA", "RNA"] else None
            c1, c2, c3 = st.columns(3)
            c1.metric("Type", seq_type)
            c2.metric("GC Content", f"{gc:.2f}%" if gc else "N/A")
            c3.metric("Length", f"{len(seq)}")

            matches = []
            try:
                prog, db = ("blastn", "nt") if seq_type in ["DNA", "RNA"] else ("blastp", "nr")
                root = ET.fromstring(NCBIWWW.qblast(prog, db, seq, hitlist_size=3).read())
                matches = [{"Gene": h.find("Hit_def").text, "ID": h.find("Hit_accession").text} for h in root.findall(".//Hit")[:3]]
            except: pass

            if matches:
                st.subheader("Top BLAST Matches")
                st.dataframe(pd.DataFrame(matches), use_container_width=True)

            bio_seq = Seq(seq)
            transcript = str(bio_seq.transcribe()) if seq_type == "DNA" else (seq if seq_type == "RNA" else None)
            protein = str(bio_seq.translate()) if seq_type == "DNA" else (str(bio_seq.translate()) if seq_type == "RNA" else seq)

            if transcript:
                with st.expander("mRNA Transcript"):
                    st.text_area("RNA", transcript, height=100)

            with st.expander("Translated Protein Sequence", expanded=True):
                st.code(protein, language=None)

            st.subheader("SWISS-MODEL Automated Prediction")
            if not swiss_token:
                st.warning("⚠️ Please provide a SWISS-MODEL API Token in the sidebar to run background predictions.")
            else:
                with st.spinner("Submitting to SWISS-MODEL and building model (this may take a few minutes)..."):
                    try:
                        headers = {"Authorization": f"Token {swiss_token}"}
                        data = {"target_sequences": protein, "project_title": "ProtCraft_Automated"}
                        res = requests.post("https://swissmodel.expasy.org/project/", headers=headers, json=data)
                        
                        if res.status_code in [200, 201, 202]:
                            proj_url = res.json().get("url")
                            status = "RUNNING"
                            
                            status_placeholder = st.empty()
                            while status in ["RUNNING", "PENDING", "QUEUED"]:
                                status_placeholder.info(f"SWISS-MODEL Status: {status}... Polling server.")
                                time.sleep(10)
                                req = requests.get(proj_url, headers=headers).json()
                                status = req.get("status", "UNKNOWN")
                            
                            status_placeholder.empty()
                            
                            if status == "COMPLETED" and req.get("models"):
                                st.success("Model successfully generated!")
                                pdb_url = req["models"][0]["coordinates_url"]
                                pdb_text = requests.get(pdb_url, headers=headers).text
                                render_protein_3d_viewer(pdb_text)
                            else:
                                st.error("SWISS-MODEL finished but could not generate a valid model for this sequence.")
                        else:
                            st.error(f"Failed to submit job. Check API Token. (Status: {res.status_code})")
                    except Exception as e:
                        st.error(f"API Error: {e}")
