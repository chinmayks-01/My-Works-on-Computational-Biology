def render_header_with_horizontal_dna():
    """Renders page header alongside an unconstrained 3D particle horizontally-stretched rotating DNA strand in gradient blue with soft blurry ends."""
    # Main page title container with DOM mounting target for canvas
    st.markdown(
        """
        <style>
        .header-wrapper {
            display: flex;
            flex-direction: row;
            align-items: center;
            justify-content: space-between;
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
        #dna-canvas-target {
            flex: 1;
            max-width: 650px;
            height: 120px;
            position: relative;
        }
        </style>
        <div class="header-wrapper">
            <div style="flex-shrink: 0;">
                <h1 style='font-size: 2.6rem; font-weight: 800; margin: 0; color: #f8fafc; white-space: nowrap;'>
                    Welcome to <span class='title-glow-text'>ProtCraft Wizard</span> 🧙‍♂️
                </h1>
            </div>
            <div id="dna-canvas-target"></div>
        </div>
        """,
        unsafe_allow_html=True,
    )

    # Injected script with gradient blue color scheme and progressive end-blurring
    dna_js = """
    <script>
    (function() {
        try {
            const parentDoc = window.parent.document;
            const targetContainer = parentDoc.getElementById('dna-canvas-target');
            if (!targetContainer) return;

            let canvas = parentDoc.getElementById('unboundedDnaCanvas');
            if (!canvas) {
                canvas = parentDoc.createElement('canvas');
                canvas.id = 'unboundedDnaCanvas';
                canvas.style.position = 'absolute';
                canvas.style.top = '-80px';
                canvas.style.left = '-80px';
                canvas.style.pointerEvents = 'auto';
                canvas.style.cursor = 'pointer';
                canvas.style.zIndex = '5';
                targetContainer.appendChild(canvas);
            }

            const ctx = canvas.getContext('2d');

            function setCanvasDimensions() {
                const rect = targetContainer.getBoundingClientRect();
                const containerWidth = rect.width > 50 ? rect.width : 580;
                canvas.width = containerWidth + 160;
                canvas.height = 280;
                canvas.style.width = canvas.width + 'px';
                canvas.style.height = canvas.height + 'px';
            }
            setCanvasDimensions();
            parentDoc.defaultView.addEventListener('resize', setCanvasDimensions);

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

            const numSteps = 42;
            const strandRadius = 24;
            const rungDotsCount = 5;

            // Helper to interpolate gradient blue colors across the strand
            function getGradientBlueColor(ratio, type, rungFraction) {
                // Electric Sky Blue -> Royal Blue -> Cobalt Gradient
                if (type === 'strand1') {
                    return ratio < 0.5 ? '#38bdf8' : '#60a5fa'; // Bright Cyan-Blue to Sky Blue
                } else if (type === 'strand2') {
                    return ratio < 0.5 ? '#2563eb' : '#1d4ed8'; // Royal Blue to Sapphire
                } else {
                    // Nucleotide rungs: soft intermediate blue gradient
                    return rungFraction < 0.5 ? '#0284c7' : '#3b82f6';
                }
            }

            class DNAParticle {
                constructor(type, indexRatio, rungFraction) {
                    this.type = type;
                    this.indexRatio = indexRatio; // 0.0 (left) to 1.0 (right)
                    this.rungFraction = rungFraction;
                    this.color = getGradientBlueColor(indexRatio, type, rungFraction);
                    
                    this.x = 0;
                    this.y = 0;
                    this.vx = 0;
                    this.vy = 0;
                    this.z = 0;
                    this.scale = 1;
                }

                calculateTarget(angle, width, height) {
                    const length = width - 180;
                    const startX = 90;
                    const cy = height / 2;

                    const currentX = startX + this.indexRatio * length;
                    const nodeAngle = angle + this.indexRatio * Math.PI * 3.6;

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
                    this.scale = 1 + this.z / 150;

                    // Mouse scatter physics
                    const dx = this.x - mouse.x;
                    const dy = this.y - mouse.y;
                    const dist = Math.hypot(dx, dy);
                    const distortRadius = 100;

                    if (mouse.active && dist < distortRadius && dist > 0) {
                        const force = (1 - dist / distortRadius) * 20;
                        const anglePush = Math.atan2(dy, dx);
                        this.vx += Math.cos(anglePush) * force;
                        this.vy += Math.sin(anglePush) * force;
                    }

                    // Elastic spring return force
                    const spring = 0.08;
                    const friction = 0.83;

                    this.vx += (target.tx - this.x) * spring;
                    this.vy += (target.ty - this.y) * spring;

                    this.vx *= friction;
                    this.vy *= friction;

                    this.x += this.vx;
                    this.y += this.vy;
                }

                draw(ctx) {
                    // Edge Factor: 1 at center, smoothly tapering down towards 0 at both ends
                    const edgeFactor = Math.sin(this.indexRatio * Math.PI);

                    // Alpha gradually softens towards both ends (from 100% to ~38%)
                    const depthAlpha = (this.z + strandRadius) / (2 * strandRadius) * 0.65 + 0.35;
                    const edgeAlphaMultiplier = 0.38 + 0.62 * edgeFactor;
                    const alpha = depthAlpha * edgeAlphaMultiplier;

                    // Blur increases gradually towards both ends for a soft out-of-focus glow
                    const endBlurAddition = (1 - edgeFactor) * 8.5;
                    const baseBlur = this.type === 'rung' ? 5 : 10;

                    const baseRadius = this.type === 'rung' ? 2.2 : 3.6;
                    const radius = Math.max(0.6, baseRadius * this.scale);

                    ctx.save();
                    ctx.shadowBlur = (baseBlur + endBlurAddition) * this.scale;
                    ctx.shadowColor = this.color;
                    ctx.fillStyle = this.color;
                    ctx.globalAlpha = Math.max(0.12, alpha);
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

                particles.push(new DNAParticle('strand1', ratio, 0));
                particles.push(new DNAParticle('strand2', ratio, 0));

                for (let j = 1; j <= rungDotsCount; j++) {
                    const rungFraction = j / (rungDotsCount + 1);
                    particles.push(new DNAParticle('rung', ratio, rungFraction));
                }
            }

            particles.forEach(p => {
                const initTarget = p.calculateTarget(0, canvas.width, canvas.height);
                p.x = initTarget.tx;
                p.y = initTarget.ty;
            });

            function animate() {
                ctx.clearRect(0, 0, canvas.width, canvas.height);
                rotationAngle += 0.02;

                particles.forEach(p => p.update(rotationAngle, canvas.width, canvas.height));
                particles.sort((a, b) => a.z - b.z);
                particles.forEach(p => p.draw(ctx));

                requestAnimationFrame(animate);
            }

            animate();
        } catch(err) {
            console.error("DNA Canvas Error:", err);
        }
    })();
    </script>
    """
    components.html(dna_js, height=0, width=0)
