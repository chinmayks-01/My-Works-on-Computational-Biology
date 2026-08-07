def inject_custom_ui_theme():
    """Injects dynamic breathing background, glassmorphism UI, equal-sized side-by-side radio cards, and cursor interactive animation."""
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

    # Executing JavaScript via Streamlit components to ensure script execution
    cursor_js = """
    <script>
    (function() {
        try {
            const parentDoc = window.parent.document;
            let canvas = parentDoc.getElementById('cursor-particle-canvas');
            
            if (!canvas) {
                canvas = parentDoc.createElement('canvas');
                canvas.id = 'cursor-particle-canvas';
                canvas.style.position = 'fixed';
                canvas.style.top = '0';
                canvas.style.left = '0';
                canvas.style.width = '100vw';
                canvas.style.height = '100vh';
                canvas.style.pointerEvents = 'none';
                canvas.style.zIndex = '999999';
                parentDoc.body.appendChild(canvas);
            }

            const ctx = canvas.getContext('2d');
            let width = canvas.width = parentDoc.documentElement.clientWidth || window.innerWidth;
            let height = canvas.height = parentDoc.documentElement.clientHeight || window.innerHeight;

            function updateBounds() {
                width = canvas.width = parentDoc.documentElement.clientWidth || window.innerWidth;
                height = canvas.height = parentDoc.documentElement.clientHeight || window.innerHeight;
            }

            parentDoc.defaultView.addEventListener('resize', updateBounds);

            const particles = [];
            const colors = ['#38bdf8', '#818cf8', '#c084fc', '#3b82f6', '#0284c7'];

            parentDoc.addEventListener('mousemove', function(e) {
                for (let i = 0; i < 2; i++) {
                    particles.push({
                        x: e.clientX,
                        y: e.clientY,
                        vx: (Math.random() - 0.5) * 1.5,
                        vy: (Math.random() - 0.5) * 1.5,
                        radius: Math.random() * 2.5 + 1,
                        color: colors[Math.floor(Math.random() * colors.length)],
                        alpha: 0.85,
                        decay: Math.random() * 0.02 + 0.015
                    });
                }
            });

            function renderParticles() {
                ctx.clearRect(0, 0, width, height);

                for (let i = 0; i < particles.length; i++) {
                    let p = particles[i];
                    p.x += p.vx;
                    p.y += p.vy;
                    p.alpha -= p.decay;

                    if (p.alpha <= 0) {
                        particles.splice(i, 1);
                        i--;
                        continue;
                    }

                    ctx.save();
                    ctx.globalAlpha = p.alpha;
                    ctx.beginPath();
                    ctx.arc(p.x, p.y, p.radius, 0, Math.PI * 2);
                    ctx.fillStyle = p.color;
                    ctx.shadowBlur = 8;
                    ctx.shadowColor = p.color;
                    ctx.fill();
                    ctx.restore();
                }

                requestAnimationFrame(renderParticles);
            }

            renderParticles();
        } catch(err) {
            console.error("Cursor particle error:", err);
        }
    })();
    </script>
    """
    components.html(cursor_js, height=0, width=0)
