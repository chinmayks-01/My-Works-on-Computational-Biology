from collections import Counter
import io
import os
import time
import urllib.parse
import xml.etree.ElementTree as ET
from textwrap import dedent
import gzip

from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
import pandas as pd
import requests
import streamlit as st
import streamlit.components.v1 as components


st.set_page_config(page_title="Bioinformatics Sequence Pipeline", layout="wide")


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
    """Injects dynamic breathing background, glassmorphism UI, and equal-sized side-by-side radio cards."""
    css = """
    <style>
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

    .block-container {
        position: relative;
        z-index: 1;
    }

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

    h1, h2, h3 {
        color: #f1f5f9 !important;
        font-weight: 700 !important;
        letter-spacing: -0.02em;
    }

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

    section[data-testid="stSidebar"] {
        background-color: rgba(10, 10, 12, 0.85) !important;
        backdrop-filter: blur(16px);
        border-right: 1px solid rgba(255, 255, 255, 0.08);
        z-index: 10;
    }
    
    header[data-testid="stHeader"] {
        background: transparent !important;
    }

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

    div[data-testid="stRadio"] div[role="radiogroup"] label:hover {
        background: rgba(56, 189, 248, 0.08) !important;
        border-color: rgba(56, 189, 248, 0.4) !important;
        transform: translateY(-1px);
    }

    div[data-testid="stRadio"] div[role="radiogroup"] label:has(input:checked) {
        background: rgba(56, 189, 248, 0.12) !important;
        border-color: #38bdf8 !important;
        box-shadow: 0 0 15px rgba(56, 189, 248, 0.2) !important;
    }
    </style>
    """
    st.markdown(css, unsafe_allow_html=True)


def render_header_with_horizontal_dna():
    """Renders page header alongside an interactive 3D particle horizontal rotating DNA strand that distorts on hover."""
    header_html = """
    <style>
    .header-wrapper {
        display: flex;
        flex-direction: row;
        align-items: center;
        justify-content: space-between;
        width: 100%;
        margin-bottom: 1rem;
        gap: 20px;
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
    .dna-container {
        flex: 1;
        max-width: 580px;
        height: 120px;
        display: flex;
        justify-content: flex-end;
        align-items: center;
    }
    #horizontalDnaCanvas {
        width: 100%;
        height: 100%;
        cursor: pointer;
    }
    </style>
    <div class="header-wrapper">
        <div style="flex-shrink: 0;">
            <h1 style='font-size: 2.6rem; font-weight: 800; margin: 0; color: #f8fafc; white-space: nowrap;'>
                Welcome to <span class='title-glow-text'>ProtCraft Wizard</span> 🧙‍♂️
            </h1>
        </div>
        <div class="dna-container">
            <canvas id="horizontalDnaCanvas"></canvas>
        </div>
    </div>

    <script>
    (function() {
        const canvas = document.getElementById('horizontalDnaCanvas');
        if (!canvas) return;
        const ctx = canvas.getContext('2d');

        function setCanvasDimensions() {
            canvas.width = canvas.parentElement.clientWidth || 550;
            canvas.height = 120;
        }
        setCanvasDimensions();
        window.addEventListener('resize', setCanvasDimensions);

        let rotationAngle = 0;
        const mouse = { x: -1000, y: -1000, active: false };

        canvas.addEventListener('mousemove', (e) => {
            const rect = canvas.getBoundingClientRect();
            mouse.x = e.clientX - rect.left;
            mouse.y = e.clientY - rect.top;
            mouse.active = true;
        });

        canvas.addEventListener('mouseleave', () => {
            mouse.x = -1000;
            mouse.y = -1000;
            mouse.active = false;
        });

        const numSteps = 38;
        const strandRadius = 28;
        const rungDotsCount = 5;

        // Initialize particle set
        class DNAParticle {
            constructor(type, indexRatio, rungFraction, color) {
                this.type = type; // 'strand1', 'strand2', or 'rung'
                this.indexRatio = indexRatio; // 0.0 to 1.0 along horizontal axis
                this.rungFraction = rungFraction; // 0.0 to 1.0 between strands
                this.color = color;
                
                this.x = 0;
                this.y = 0;
                this.vx = 0;
                this.vy = 0;
                this.z = 0;
                this.scale = 1;
            }

            calculateTarget(angle, width, height) {
                const length = Math.min(width - 20, 520);
                const startX = (width - length) / 2;
                const cy = height / 2;

                const currentX = startX + this.indexRatio * length;
                const nodeAngle = angle + this.indexRatio * Math.PI * 7;

                if (this.type === 'strand1') {
                    const ty = cy + Math.sin(nodeAngle) * strandRadius;
                    const tz = Math.cos(nodeAngle) * strandRadius;
                    return { tx: currentX, ty: ty, tz: tz };
                } else if (this.type === 'strand2') {
                    const ty = cy + Math.sin(nodeAngle + Math.PI) * strandRadius;
                    const tz = Math.cos(nodeAngle + Math.PI) * strandRadius;
                    return { tx: currentX, ty: ty, tz: tz };
                } else {
                    // Base-pair rung particle interpolation
                    const y1 = cy + Math.sin(nodeAngle) * strandRadius;
                    const z1 = Math.cos(nodeAngle) * strandRadius;

                    const y2 = cy + Math.sin(nodeAngle + Math.PI) * strandRadius;
                    const z2 = Math.cos(nodeAngle + Math.PI) * strandRadius;

                    const ty = y1 + (y2 - y1) * this.rungFraction;
                    const tz = z1 + (z2 - z1) * this.rungFraction;
                    return { tx: currentX, ty: ty, tz: tz };
                }
            }

            update(angle, width, height) {
                const target = this.calculateTarget(angle, width, height);
                this.z = target.tz;
                this.scale = 1 + this.z / 140;

                // Distortion physics when mouse approaches
                const dx = this.x - mouse.x;
                const dy = this.y - mouse.y;
                const dist = Math.hypot(dx, dy);
                const distortRadius = 75;

                if (mouse.active && dist < distortRadius && dist > 0) {
                    const force = (1 - dist / distortRadius) * 14;
                    const anglePush = Math.atan2(dy, dx);
                    this.vx += Math.cos(anglePush) * force;
                    this.vy += Math.sin(anglePush) * force;
                }

                // Elastic spring force returning particle to target helix coordinate
                const spring = 0.08;
                const friction = 0.82;

                this.vx += (target.tx - this.x) * spring;
                this.vy += (target.ty - this.y) * spring;

                this.vx *= friction;
                this.vy *= friction;

                this.x += this.vx;
                this.y += this.vy;
            }

            draw(ctx) {
                const alpha = (this.z + strandRadius) / (2 * strandRadius) * 0.65 + 0.35;
                const baseRadius = this.type === 'rung' ? 2.2 : 3.6;
                const radius = Math.max(0.6, baseRadius * this.scale);

                ctx.save();
                ctx.shadowBlur = (this.type === 'rung' ? 6 : 10) * this.scale;
                ctx.shadowColor = this.color;
                ctx.fillStyle = this.color;
                ctx.globalAlpha = Math.max(0.15, alpha);
                ctx.beginPath();
                ctx.arc(this.x, this.y, radius, 0, Math.PI * 2);
                ctx.fill();
                ctx.restore();
            }
        }

        const particles = [];

        // Build DNA Particle Mesh
        for (let i = 0; i <= numSteps; i++) {
            const ratio = i / numSteps;

            // Strand 1 Particle (Cyan)
            particles.push(new DNAParticle('strand1', ratio, 0, '#38bdf8'));

            // Strand 2 Particle (Purple)
            particles.push(new DNAParticle('strand2', ratio, 0, '#c084fc'));

            // Nucleotide Rung Dots (Magenta / Pink)
            for (let j = 1; j <= rungDotsCount; j++) {
                const rungFraction = j / (rungDotsCount + 1);
                particles.push(new DNAParticle('rung', ratio, rungFraction, '#f472b6'));
            }
        }

        // Initialize particle positions
        particles.forEach(p => {
            const initTarget = p.calculateTarget(0, canvas.width, canvas.height);
            p.x = initTarget.tx;
            p.y = initTarget.ty;
        });

        function animate() {
            ctx.clearRect(0, 0, canvas.width, canvas.height);
            rotationAngle += 0.022;

            // Update particle physics
            particles.forEach(p => p.update(rotationAngle, canvas.width, canvas.height));

            // Z-sorting for realistic 3D depth perception
            particles.sort((a, b) => a.z - b.z);

            // Draw particles
            particles.forEach(p => p.draw(ctx));

            requestAnimationFrame(animate);
        }

        animate();
    })();
    </script>
    """
    components.html(header_html, height=130, scrolling=False)


def render_header_with_horizontal_dna():
    """Renders page header alongside a wide interactive 3D particle horizontal rotating DNA strand that distorts freely on hover."""
    header_html = """
    <style>
    .header-wrapper {
        display: flex;
        flex-direction: row;
        align-items: center;
        justify-content: space-between;
        width: 100%;
        margin-bottom: 0.5rem;
        gap: 20px;
        overflow: visible;
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
    .dna-container {
        flex: 1;
        max-width: 620px;
        height: 170px;
        display: flex;
        justify-content: flex-end;
        align-items: center;
        overflow: visible;
    }
    #horizontalDnaCanvas {
        width: 100%;
        height: 100%;
        cursor: pointer;
        overflow: visible;
    }
    </style>
    <div class="header-wrapper">
        <div style="flex-shrink: 0;">
            <h1 style='font-size: 2.6rem; font-weight: 800; margin: 0; color: #f8fafc; white-space: nowrap;'>
                Welcome to <span class='title-glow-text'>ProtCraft Wizard</span> 🧙‍♂️
            </h1>
        </div>
        <div class="dna-container">
            <canvas id="horizontalDnaCanvas"></canvas>
        </div>
    </div>

    <script>
    (function() {
        const canvas = document.getElementById('horizontalDnaCanvas');
        if (!canvas) return;
        const ctx = canvas.getContext('2d');

        function setCanvasDimensions() {
            canvas.width = canvas.parentElement.clientWidth || 580;
            canvas.height = 170;
        }
        setCanvasDimensions();
        window.addEventListener('resize', setCanvasDimensions);

        let rotationAngle = 0;
        const mouse = { x: -1000, y: -1000, active: false };

        canvas.addEventListener('mousemove', (e) => {
            const rect = canvas.getBoundingClientRect();
            mouse.x = e.clientX - rect.left;
            mouse.y = e.clientY - rect.top;
            mouse.active = true;
        });

        canvas.addEventListener('mouseleave', () => {
            mouse.x = -1000;
            mouse.y = -1000;
            mouse.active = false;
        });

        const numSteps = 38;
        const strandRadius = 46; // Increased width for a bolder helix
        const rungDotsCount = 6;

        class DNAParticle {
            constructor(type, indexRatio, rungFraction, color) {
                this.type = type; // 'strand1', 'strand2', or 'rung'
                this.indexRatio = indexRatio;
                this.rungFraction = rungFraction;
                this.color = color;
                
                this.x = 0;
                this.y = 0;
                this.vx = 0;
                this.vy = 0;
                this.z = 0;
                this.scale = 1;
            }

            calculateTarget(angle, width, height) {
                const length = Math.min(width - 20, 540);
                const startX = (width - length) / 2;
                const cy = height / 2;

                const currentX = startX + this.indexRatio * length;
                const nodeAngle = angle + this.indexRatio * Math.PI * 6.5;

                if (this.type === 'strand1') {
                    const ty = cy + Math.sin(nodeAngle) * strandRadius;
                    const tz = Math.cos(nodeAngle) * strandRadius;
                    return { tx: currentX, ty: ty, tz: tz };
                } else if (this.type === 'strand2') {
                    const ty = cy + Math.sin(nodeAngle + Math.PI) * strandRadius;
                    const tz = Math.cos(nodeAngle + Math.PI) * strandRadius;
                    return { tx: currentX, ty: ty, tz: tz };
                } else {
                    const y1 = cy + Math.sin(nodeAngle) * strandRadius;
                    const z1 = Math.cos(nodeAngle) * strandRadius;

                    const y2 = cy + Math.sin(nodeAngle + Math.PI) * strandRadius;
                    const z2 = Math.cos(nodeAngle + Math.PI) * strandRadius;

                    const ty = y1 + (y2 - y1) * this.rungFraction;
                    const tz = z1 + (z2 - z1) * this.rungFraction;
                    return { tx: currentX, ty: ty, tz: tz };
                }
            }

            update(angle, width, height) {
                const target = this.calculateTarget(angle, width, height);
                this.z = target.tz;
                this.scale = 1 + this.z / 160;

                // Expanded repulsion zone so particles fly freely on hover
                const dx = this.x - mouse.x;
                const dy = this.y - mouse.y;
                const dist = Math.hypot(dx, dy);
                const distortRadius = 110;

                if (mouse.active && dist < distortRadius && dist > 0) {
                    const force = (1 - dist / distortRadius) * 22;
                    const anglePush = Math.atan2(dy, dx);
                    this.vx += Math.cos(anglePush) * force;
                    this.vy += Math.sin(anglePush) * force;
                }

                // Elastic spring return force
                const spring = 0.075;
                const friction = 0.84;

                this.vx += (target.tx - this.x) * spring;
                this.vy += (target.ty - this.y) * spring;

                this.vx *= friction;
                this.vy *= friction;

                this.x += this.vx;
                this.y += this.vy;
            }

            draw(ctx) {
                const alpha = (this.z + strandRadius) / (2 * strandRadius) * 0.65 + 0.35;
                const baseRadius = this.type === 'rung' ? 2.4 : 3.8;
                const radius = Math.max(0.6, baseRadius * this.scale);

                ctx.save();
                ctx.shadowBlur = (this.type === 'rung' ? 6 : 12) * this.scale;
                ctx.shadowColor = this.color;
                ctx.fillStyle = this.color;
                ctx.globalAlpha = Math.max(0.15, alpha);
                ctx.beginPath();
                ctx.arc(this.x, this.y, radius, 0, Math.PI * 2);
                ctx.fill();
                ctx.restore();
            }
        }

        const particles = [];

        // Construct Mesh
        for (let i = 0; i <= numSteps; i++) {
            const ratio = i / numSteps;

            particles.push(new DNAParticle('strand1', ratio, 0, '#38bdf8'));
            particles.push(new DNAParticle('strand2', ratio, 0, '#c084fc'));

            for (let j = 1; j <= rungDotsCount; j++) {
                const rungFraction = j / (rungDotsCount + 1);
                particles.push(new DNAParticle('rung', ratio, rungFraction, '#f472b6'));
            }
        }

        particles.forEach(p => {
            const initTarget = p.calculateTarget(0, canvas.width, canvas.height);
            p.x = initTarget.tx;
            p.y = initTarget.ty;
        });

        function animate() {
            ctx.clearRect(0, 0, canvas.width, canvas.height);
            rotationAngle += 0.022;

            particles.forEach(p => p.update(rotationAngle, canvas.width, canvas.height));
            particles.sort((a, b) => a.z - b.z);
            particles.forEach(p => p.draw(ctx));

            requestAnimationFrame(animate);
        }

        animate();
    })();
    </script>
    """
    components.html(header_html, height=180, scrolling=False)


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
    <div style="background: rgba(15, 23, 42, 0.7); border: 1px solid rgba(255, 255, 255, 0.1); border-radius: 12px; padding: 20px; font-family: monospace;">
        <div style="font-size: 13px; color: #94a3b8; margin-bottom: 12px; display: flex; justify-content: space-between;">
            <span><b>6-Frame ORF Graphic Map</b></span>
            <span>Total Length: {total_seq_len} nt</span>
        </div>
        <div style="position: relative; width: 100%;">
    """

    for f in frames:
        html += f"""
        <div style="display: flex; align-items: center; margin-bottom: 8px;">
            <div style="width: 80px; font-size: 11px; color: #cbd5e1; font-weight: bold;">{f}</div>
            <div style="position: relative; flex-grow: 1; height: 16px; background: rgba(255,255,255,0.04); border-radius: 4px; border: 1px solid rgba(255,255,255,0.06);">
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
                 style="position: absolute; left: {left_pct}%; width: {width_pct}%; height: 100%; background: {color}; opacity: 0.85; border-radius: 3px; cursor: pointer; box-shadow: 0 0 8px {color};">
            </div>
            """

        html += """
            </div>
        </div>
        """

    html += """
        </div>
        <div style="display: flex; gap: 15px; font-size: 11px; color: #94a3b8; margin-top: 15px; justify-content: flex-end;">
            <span style="display: flex; align-items: center; gap: 5px;"><span style="width:10px; height:10px; background:#38bdf8; display:inline-block; border-radius:2px;"></span> + Frames</span>
            <span style="display: flex; align-items: center; gap: 5px;"><span style="width:10px; height:10px; background:#f43f5e; display:inline-block; border-radius:2px;"></span> - Frames</span>
        </div>
    </div>
    """
    components.html(html, height=260, scrolling=False)


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
        "A": "#80a0f0",
        "I": "#80a0f0",
        "L": "#80a0f0",
        "M": "#80a0f0",
        "F": "#80a0f0",
        "W": "#80a0f0",
        "V": "#80a0f0",
        "R": "#f01505",
        "K": "#f01505",
        "N": "#00ff00",
        "Q": "#00ff00",
        "D": "#c000c0",
        "E": "#c000c0",
        "C": "#f08080",
        "G": "#f09040",
        "P": "#ffff00",
        "H": "#15a4a4",
        "Y": "#15a4a4",
        "S": "#15a400",
        "T": "#15a400",
    }
    styled_html = "<div style='font-family: monospace; font-size: 15px; word-break: break-all; line-height: 2.0; background-color: #222; padding: 14px; border-radius: 6px; letter-spacing: 1px;'>"
    for aa in seq:
        bg = bg_colors.get(aa, "#ffffff")
        text_color = (
            "#ffffff" if aa in ["R", "K", "S", "T", "D", "E"] else "#000000"
        )
        styled_html += f"<span style='background-color: {bg}; color: {text_color}; font-weight: bold; padding: 2px 5px; margin: 1px 0px; display: inline-block; text-align: center; border-radius: 2px;'>{aa}</span>"
    styled_html += "</div>"
    return styled_html


# Apply Custom Theme CSS
inject_custom_ui_theme()

# Render Header with Horizontal Interactive Particle DNA Strand
render_header_with_horizontal_dna()

st.sidebar.header("Settings")
user_email = st.sidebar.text_input(
    "NCBI Entrez Email", value="your.email@example.com"
)
Entrez.email = user_email
swiss_token = st.sidebar.text_input("SWISS-MODEL API Token", type="password", help="Required for background protein structure prediction.")

st.header("Input Sequence")
st.markdown(
    "<p style='font-size: 1rem; font-weight: 500; margin-bottom: 0.5rem; color: #e2e8f0;'>Choose Input Method:</p>",
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

            # --- SWISS-MODEL AUTOMATED BACKGROUND PREDICTION ---
            st.header("Protein Structure Prediction (SWISS-MODEL API)")
            
            if protein_seq:
                if not swiss_token:
                    st.warning("⚠️ Please provide your SWISS-MODEL API Token in the sidebar to enable background prediction.")
                else:
                    with st.spinner("Submitting sequence to SWISS-MODEL and generating homology model (this may take a few minutes)..."):
                        try:
                            headers = {
                                "Authorization": f"Token {swiss_token}",
                                "Content-Type": "application/json",
                                "Accept": "application/json"
                            }
                            
                            clean_protein_seq = protein_seq.replace("*", "").strip()
                            
                            data = {
                                "target_sequences": [clean_protein_seq], 
                                "project_title": "ProtCraft Automated Job"
                            }
                            res = requests.post("https://swissmodel.expasy.org/automodel/", headers=headers, json=data)
                            
                            if res.status_code in [200, 201, 202]:
                                resp_json = res.json()
                                project_id = resp_json.get("project_id")
                                status = resp_json.get("status", "QUEUED")
                                
                                if not project_id:
                                    st.error(f"API did not return a project_id. Full response: {resp_json}")
                                else:
                                    status_placeholder = st.empty()
                                    poll_url = f"https://swissmodel.expasy.org/project/{project_id}/models/summary/"
                                    
                                    poll_req = None
                                    while status in ["RUNNING", "PENDING", "QUEUED"]:
                                        status_placeholder.info(f"SWISS-MODEL API Status: {status} (ID: {project_id})... Polling server (Please wait).")
                                        time.sleep(10)
                                        
                                        poll_req = requests.get(poll_url, headers=headers)
                                        if poll_req.status_code == 200:
                                            status = poll_req.json().get("status", "UNKNOWN")
                                        else:
                                            st.error(f"Polling error {poll_req.status_code}: {poll_req.text}")
                                            status = "API_ERROR"
                                            break
                                    
                                    status_placeholder.empty()
                                    
                                    if status == "COMPLETED" and poll_req:
                                        models = poll_req.json().get("models", [])
                                        if models:
                                            st.success("Homology model successfully generated!")
                                            
                                            pdb_url = models[0].get("coordinates_url") 
                                            if not pdb_url:
                                                pdb_url = f"https://swissmodel.expasy.org/project/{project_id}/models/01.pdb"
                                                
                                            # Download the file
                                            pdb_res = requests.get(pdb_url, headers=headers)
                                            if pdb_res.status_code == 200:
                                                # Try to decompress the response content if it is gzipped
                                                try:
                                                    pdb_bytes = gzip.decompress(pdb_res.content)
                                                    pdb_text = pdb_bytes.decode('utf-8')
                                                except Exception:
                                                    # If decompression fails, assume it's already plain text
                                                    pdb_text = pdb_res.text
                                                
                                                if "ATOM" in pdb_text or "HEADER" in pdb_text:
                                                    render_protein_3d_viewer(pdb_text, height=500)
                                                else:
                                                    st.error("SWISS-MODEL returned an invalid file (not a PDB).")
                                                    st.code(pdb_text[:500])
                                            else:
                                                st.error(f"Failed to download the generated PDB coordinate file. (Code: {pdb_res.status_code})")
                                        else:
                                            st.error("SWISS-MODEL job finished, but no suitable template was found to generate a valid model.")
                                    elif status in ["FAILED", "UNKNOWN", "API_ERROR"]:
                                        st.error(f"SWISS-MODEL job stopped with status: {status}")
                            else:
                                st.error(f"Failed to submit job. Please check your API Token and sequence. (Status Code: {res.status_code})\n\n{res.text}")
                        except Exception as e:
                            st.error(f"Failed to connect to SWISS-MODEL API: {e}")
            else:
                st.warning("No protein sequence available for modeling.")
