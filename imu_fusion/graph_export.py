''' Export the REAL GTSAM factor graph as a typed, time-stamped structure and an
    interactive + animated, self-contained HTML viewer.

    The static `run_study.plot_factorgraph` is a hand-drawn schematic.  This module
    instead walks the SAME wiring as `celestial_factor_graph.build_graph` and emits
    every variable and factor with its type, connected variables, sigma, and the
    streaming `step` at which it appears / marginalises out of the fixed-lag window
    (`realtime.StreamingEstimator` semantics).  `graph_structure` cross-checks
    itself against the real `build_graph` graph (`size()` / `keys()`), so the export
    cannot silently drift from what GTSAM actually solves.

    `write_graph_viewer` renders a single self-contained HTML (inlined JSON + vanilla
    JS/SVG force layout, no external hosts) suitable for publishing as an Artifact.

    (c) 2026.  MIT License (see LICENSE file).
'''

import json

# Factor-type -> (colour, human label).  Palette consistent with plot_factorgraph.
FACTOR_STYLE = {
    "prior":   ("#d1495b", "attitude + height prior"),
    "alt_Sun": ("#c8992e", "altitude ☉"),
    "alt_Moon": ("#5b6b86", "altitude ☾"),
    "az_Sun":  ("#e0ad4b", "azimuth ☉"),
    "az_Moon": ("#8aa0bf", "azimuth ☾"),
    "q_Sun":   ("#b07d1a", "parallactic q ☉"),
    "q_Moon":  ("#6b7f9e", "parallactic q ☾"),
    "diff":    ("#3ddc97", "Δq (Sun−Moon) — horizon-free"),
    "imu":     ("#edae49", "IMU preintegration"),
    "vprior":  ("#9bb0bd", "velocity prior"),
    "bprior":  ("#c9a978", "bias prior"),
}
VAR_STYLE = {"pose": ("#e8f1f3", "pose X (position+attitude)"),
             "vel": ("#e7efe8", "velocity V"),
             "bias": ("#f2e2c2", "IMU bias B")}


def _marg_step(times, i, lag_s):
    ''' The keyframe index at which keyframe i falls out of the fixed-lag window
        (first later shot more than lag_s newer), or None if it never does. '''
    for j in range(i + 1, len(times)):
        if times[j] - times[i] > lag_s:
            return j
    return None


def graph_structure(scenario, use_imu=True, use_azimuth=True,
                    use_parallactic=True, use_differential=True,
                    lag_s=30.0, validate=True):
    ''' Typed, time-stamped node/factor structure of the real graph.  Mirrors
        `celestial_factor_graph.build_graph`; when `validate`, asserts the factor
        count and variable ids match the actual GTSAM graph. '''
    kfs = scenario.keyframes
    times = [kf.time_s for kf in kfs]
    nodes, factors = [], []

    def marg(i):
        return _marg_step(times, i, lag_s)

    # --- variables ---
    for kf in kfs:
        nodes.append(dict(id=f"x{kf.index}", label=f"X{kf.index}", kind="pose",
                          step=kf.index, marg=marg(kf.index)))
    if use_imu:
        for kf in kfs:
            nodes.append(dict(id=f"v{kf.index}", label=f"V{kf.index}", kind="vel",
                              step=kf.index, marg=marg(kf.index)))
        nodes.append(dict(id="b0", label="B", kind="bias", step=0, marg=None))

    def add(fid, ftype, vs, step, mstep, sigma):
        factors.append(dict(id=fid, type=ftype, vars=vs, step=step,
                            marg=mstep, sigma=sigma))

    # --- factors, in build_graph order ---
    for kf in kfs:
        i, x = kf.index, f"x{kf.index}"
        m = marg(i)
        add(f"prior{i}", "prior", [x], i, m, "attitude 1e-3 rad / E,N loose")
        for o in kf.observations:
            b = o.body
            add(f"alt_{b}_{i}", f"alt_{b}", [x], i, m,
                f"{o.alt_sigma_arcmin:.2f}′")
            if use_azimuth:
                add(f"az_{b}_{i}", f"az_{b}", [x], i, m,
                    f"{o.az_sigma_arcmin:.2f}′")
            if use_parallactic and o.par_valid:
                add(f"q_{b}_{i}", f"q_{b}", [x], i, m,
                    f"{o.par_sigma_deg:.3f}°")
        if (use_parallactic and use_differential
                and getattr(kf, "diff_valid", False)):
            add(f"diff{i}", "diff", [x], i, m, f"{kf.diff_q_sigma:.3f}° horizon-free")
    if use_imu:
        for k in range(len(kfs) - 1):
            add(f"imu{k}", "imu",
                [f"x{k}", f"v{k}", f"x{k+1}", f"v{k+1}", "b0"], k + 1,
                marg(k), "preintegrated accel+gyro")
        add("vprior", "vprior", ["v0"], 0, marg(0), "5 m/s")
        add("bprior", "bprior", ["b0"], 0, None, "0.1")

    struct = dict(
        nodes=nodes, factors=factors,
        n_shots=len(kfs), lag_s=lag_s,
        epoch=kfs[0].time_iso if kfs else "",
        factor_style=FACTOR_STYLE, var_style=VAR_STYLE,
        counts=dict(variables=len(nodes), factors=len(factors)))

    if validate:
        import gtsam
        from .celestial_factor_graph import build_graph
        g, init = build_graph(scenario, use_imu, use_azimuth, use_parallactic,
                              use_differential=use_differential)
        real_vars = {gtsam.DefaultKeyFormatter(k) for k in init.keys()}
        got_vars = {n["id"] for n in nodes}
        assert len(factors) == g.size(), \
            f"factor count {len(factors)} != real graph {g.size()}"
        assert got_vars == real_vars, \
            f"variable ids {got_vars ^ real_vars} differ from real graph"
        for f in factors:
            for v in f["vars"]:
                assert v in got_vars, f"factor {f['id']} references missing {v}"
    return struct


def write_graph_json(structure, path):
    with open(path, "w") as fh:
        json.dump(structure, fh)


def write_graph_viewer(structure, path):
    ''' Render a self-contained interactive + animated HTML viewer (no external
        hosts) of the given graph structure. '''
    data = json.dumps(structure)
    html = _HTML_TEMPLATE.replace("/*__DATA__*/", data)
    with open(path, "w") as fh:
        fh.write(html)


_HTML_TEMPLATE = r"""<title>Live factor graph — iPhone Sun+Moon fix</title>
<style>
:root{--bg:#0f1419;--panel:#1a222c;--fg:#e6edf3;--mut:#9db0c0;--line:#2b3946}
@media(prefers-color-scheme:light){:root{--bg:#f5f7fa;--panel:#fff;--fg:#12212e;--mut:#4a5b6b;--line:#dfe6ee}}
:root[data-theme=dark]{--bg:#0f1419;--panel:#1a222c;--fg:#e6edf3;--mut:#9db0c0;--line:#2b3946}
:root[data-theme=light]{--bg:#f5f7fa;--panel:#fff;--fg:#12212e;--mut:#4a5b6b;--line:#dfe6ee}
*{box-sizing:border-box}html,body{margin:0;height:100%;background:var(--bg);color:var(--fg);
font:14px/1.5 -apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif}
#wrap{display:flex;flex-direction:column;height:100vh}
header{padding:.6rem .9rem;border-bottom:1px solid var(--line)}
h1{font-size:1rem;margin:0}.sub{color:var(--mut);font-size:.8rem;margin-top:.15rem}
#main{flex:1;position:relative;overflow:hidden}
svg{width:100%;height:100%;display:block;cursor:grab}svg.drag{cursor:grabbing}
.edge{stroke:var(--line);stroke-width:1}
.varlbl{font-size:9px;fill:var(--fg);pointer-events:none;text-anchor:middle;dominant-baseline:central}
.fade{opacity:.12}
#controls{padding:.5rem .9rem;border-top:1px solid var(--line);display:flex;
gap:.8rem;align-items:center;flex-wrap:wrap;background:var(--panel)}
button{background:#00b4d8;border:0;color:#012;border-radius:6px;padding:.3rem .7rem;
font-weight:600;cursor:pointer}button.sec{background:var(--line);color:var(--fg)}
input[type=range]{flex:1;min-width:120px}
#legend{position:absolute;top:.6rem;right:.6rem;background:var(--panel);border:1px solid var(--line);
border-radius:8px;padding:.5rem .6rem;max-width:230px;font-size:.78rem}
.chip{display:flex;align-items:center;gap:.4rem;padding:.12rem 0;cursor:pointer;user-select:none}
.chip.off{opacity:.35}.sw{width:11px;height:11px;border-radius:2px;flex:0 0 auto}
.sw.circle{border-radius:50%}
#tip{position:absolute;pointer-events:none;background:var(--panel);border:1px solid var(--line);
border-radius:6px;padding:.35rem .5rem;font-size:.76rem;display:none;max-width:220px;z-index:5}
.mono{font-variant-numeric:tabular-nums}
</style>
<div id="wrap">
<header><h1>Live factor graph — <span id="stitle"></span></h1>
<div class="sub" id="ssub"></div></header>
<div id="main"><svg id="svg"></svg><div id="legend"></div><div id="tip"></div></div>
<div id="controls">
<button id="play">▶ Play</button>
<button class="sec" id="mode">Streaming: on</button>
<span class="mono" id="steplbl" style="min-width:8em"></span>
<input type="range" id="slider" min="0" value="0" step="1">
<button class="sec" id="reset">Reset view</button>
</div></div>
<script>
const D=/*__DATA__*/;
const NS='http://www.w3.org/2000/svg';
const svg=document.getElementById('svg'), tip=document.getElementById('tip');
document.getElementById('stitle').textContent=`${D.counts.variables} variables · ${D.counts.factors} factors`;
document.getElementById('ssub').textContent=`Real GTSAM graph · canonical epoch ${D.epoch} · ${D.n_shots} Sun+Moon shots · fixed-lag ${D.lag_s}s. Drag nodes, scroll to zoom, toggle factor types, scrub the timeline.`;
// build layout nodes (variables + factors) and edges
const vmap={}; D.nodes.forEach(n=>{n.kind_='var';vmap[n.id]=n});
const fmap={}; D.factors.forEach(f=>{f.kind_='fac';fmap[f.id]=f});
const L=[...D.nodes,...D.factors];
const edges=[]; D.factors.forEach(f=>f.vars.forEach(v=>{if(vmap[v])edges.push({a:f,b:vmap[v]})}));
// initial positions: poses along a line, others near their pose
let W=svg.clientWidth||900,H=svg.clientHeight||600;
L.forEach((n,i)=>{n.x=W/2+Math.cos(i*2.4)*160*Math.random()+ (n.step||0)*40-D.n_shots*20;
 n.y=H/2+Math.sin(i*2.4)*120;n.vx=0;n.vy=0;n.pin=false;});
D.nodes.forEach(n=>{if(n.kind==='pose'){n.x=120+n.step*(Math.max(W-240,300)/Math.max(D.n_shots-1,1));n.y=H*0.5;}});
// force sim
function tick(){
 const vis=L.filter(n=>n._vis);
 for(let i=0;i<vis.length;i++){const a=vis[i];
  for(let j=i+1;j<vis.length;j++){const b=vis[j];
   let dx=a.x-b.x,dy=a.y-b.y,d2=dx*dx+dy*dy+.01,d=Math.sqrt(d2);
   let f=1400/d2; if(d<1){dx=Math.random();dy=Math.random();d=1;}
   const fx=f*dx/d,fy=f*dy/d;
   if(!a.pin){a.vx+=fx;a.vy+=fy;} if(!b.pin){b.vx-=fx;b.vy-=fy;}}}
 edges.forEach(e=>{if(!e.a._vis||!e.b._vis)return;
  let dx=e.b.x-e.a.x,dy=e.b.y-e.a.y,d=Math.sqrt(dx*dx+dy*dy)+.01,k=(d-46)*0.02;
  const fx=k*dx/d,fy=k*dy/d;
  if(!e.a.pin){e.a.vx+=fx;e.a.vy+=fy;} if(!e.b.pin){e.b.vx-=fx;e.b.vy-=fy;}});
 vis.forEach(n=>{if(n.pin)return; if(n.kind==='pose'){n.vy+=(H*0.5-n.y)*0.02;}
  n.vx*=0.82;n.vy*=0.82;n.x+=Math.max(-30,Math.min(30,n.vx));n.y+=Math.max(-30,Math.min(30,n.vy));});
}
// rendering
const gE=el('g'),gF=el('g'),gN=el('g'),gT=el('g');
svg.append(gE,gF,gN,gT);
function el(t){return document.createElementNS(NS,t);}
const shapes=new Map();
L.forEach(n=>{
 let s;
 if(n.kind_==='var'){s=el('circle');s.setAttribute('r',n.kind==='pose'?15:n.kind==='bias'?12:11);
  s.setAttribute('fill',D.var_style[n.kind][0]);s.setAttribute('stroke','#23303a');s.setAttribute('stroke-width',1.4);gN.append(s);
  const tx=el('text');tx.setAttribute('class','varlbl');tx.textContent=n.label;gT.append(tx);n._t=tx;
 }else{s=el('rect');const w=n.type==='imu'?13:10;s.setAttribute('width',w);s.setAttribute('height',w);
  s.setAttribute('rx',2);s.setAttribute('fill',D.factor_style[n.type][0]);s.setAttribute('stroke','#23303a');gF.append(s);}
 s.style.cursor='pointer';shapes.set(n.id,s);
 s.addEventListener('mousemove',ev=>showTip(ev,n));s.addEventListener('mouseleave',()=>tip.style.display='none');
 s.addEventListener('mousedown',ev=>startDrag(ev,n));
});
const elines=edges.map(e=>{const l=el('line');l.setAttribute('class','edge');gE.append(l);return l;});
function draw(){
 L.forEach(n=>{const s=shapes.get(n.id);const on=n._vis, faded=n._fade;
  s.style.display=on?'':'none'; s.setAttribute('class',faded?'fade':'');
  if(n.kind_==='var'){s.setAttribute('cx',n.x);s.setAttribute('cy',n.y);
   n._t.style.display=on?'':'none';n._t.setAttribute('class','varlbl'+(faded?' fade':''));
   n._t.setAttribute('x',n.x);n._t.setAttribute('y',n.y);}
  else{const w=n.type==='imu'?13:10;s.setAttribute('x',n.x-w/2);s.setAttribute('y',n.y-w/2);}});
 edges.forEach((e,i)=>{const l=elines[i];const on=e.a._vis&&e.b._vis;l.style.display=on?'':'none';
  l.setAttribute('class','edge'+((e.a._fade||e.b._fade)?' fade':''));
  l.setAttribute('x1',e.a.x);l.setAttribute('y1',e.a.y);l.setAttribute('x2',e.b.x);l.setAttribute('y2',e.b.y);});
}
// visibility from timeline + legend + mode
const hidden=new Set(); let streaming=true; let step=D.n_shots;
function applyVis(){
 L.forEach(n=>{
  const appeared=(n.step||0)<=step;
  const typeHidden=n.kind_==='fac'&&hidden.has(n.type);
  const marged=streaming&&n.marg!=null&&n.marg<=step;
  n._vis=appeared&&!typeHidden&&!(marged&&false); // keep marged but faded
  n._fade=marged;
  if(marged)n._vis=appeared&&!typeHidden;
 });
}
// legend
const lg=document.getElementById('legend');
const seen=new Set();
Object.keys(D.var_style).forEach(k=>{const [c,lab]=D.var_style[k];
 lg.append(chip(c,lab,null,true));});
D.factors.forEach(f=>{if(seen.has(f.type))return;seen.add(f.type);
 const[c,lab]=D.factor_style[f.type];lg.append(chip(c,lab,f.type,false));});
function chip(color,label,type,circle){const d=document.createElement('div');d.className='chip';
 const sw=document.createElement('span');sw.className='sw'+(circle?' circle':'');sw.style.background=color;
 const tx=document.createElement('span');tx.textContent=label;d.append(sw,tx);
 if(type!=null)d.addEventListener('click',()=>{if(hidden.has(type)){hidden.delete(type);d.classList.remove('off');}
  else{hidden.add(type);d.classList.add('off');}applyVis();draw();});
 return d;}
// tooltip
function showTip(ev,n){tip.style.display='block';
 const r=svg.getBoundingClientRect();tip.style.left=(ev.clientX-r.left+12)+'px';tip.style.top=(ev.clientY-r.top+12)+'px';
 if(n.kind_==='var')tip.innerHTML=`<b>${n.label}</b><br>${D.var_style[n.kind][1]}<br>appears @ shot ${n.step}`+(n.marg!=null?`<br>marginalises @ shot ${n.marg}`:'');
 else tip.innerHTML=`<b>${D.factor_style[n.type][1]}</b><br>σ = ${n.sigma}<br>on ${n.vars.join(', ')}<br>added @ shot ${n.step}`;
}
// drag + pan + zoom via viewBox
let vb={x:0,y:0,w:W,h:H};function setVB(){svg.setAttribute('viewBox',`${vb.x} ${vb.y} ${vb.w} ${vb.h}`);}
function toWorld(ev){const r=svg.getBoundingClientRect();return{x:vb.x+(ev.clientX-r.left)/r.width*vb.w,y:vb.y+(ev.clientY-r.top)/r.height*vb.h};}
let drag=null,pan=null;
function startDrag(ev,n){ev.stopPropagation();drag=n;n.pin=true;svg.classList.add('drag');}
svg.addEventListener('mousedown',ev=>{if(drag)return;pan={x:ev.clientX,y:ev.clientY,vx:vb.x,vy:vb.y};svg.classList.add('drag');});
window.addEventListener('mousemove',ev=>{if(drag){const p=toWorld(ev);drag.x=p.x;drag.y=p.y;drag.vx=drag.vy=0;}
 else if(pan){const r=svg.getBoundingClientRect();vb.x=pan.vx-(ev.clientX-pan.x)/r.width*vb.w;vb.y=pan.vy-(ev.clientY-pan.y)/r.height*vb.h;setVB();}});
window.addEventListener('mouseup',()=>{if(drag){drag.pin=false;}drag=null;pan=null;svg.classList.remove('drag');});
svg.addEventListener('wheel',ev=>{ev.preventDefault();const p=toWorld(ev);const f=ev.deltaY>0?1.12:0.9;
 vb.x=p.x-(p.x-vb.x)*f;vb.y=p.y-(p.y-vb.y)*f;vb.w*=f;vb.h*=f;setVB();},{passive:false});
document.getElementById('reset').onclick=()=>{vb={x:0,y:0,w:W,h:H};setVB();};
// timeline
const slider=document.getElementById('slider');slider.max=D.n_shots;slider.value=D.n_shots;
const steplbl=document.getElementById('steplbl');
function setStep(s){step=+s;slider.value=s;steplbl.textContent=`shot ${s} / ${D.n_shots}`;applyVis();draw();}
slider.oninput=e=>setStep(e.target.value);
let playing=false,timer=null;
document.getElementById('play').onclick=function(){playing=!playing;this.textContent=playing?'⏸ Pause':'▶ Play';
 if(playing){if(step>=D.n_shots)setStep(0);timer=setInterval(()=>{if(step>=D.n_shots){playing=false;this.textContent='▶ Play';clearInterval(timer);return;}setStep(step+1);},900);}
 else clearInterval(timer);};
document.getElementById('mode').onclick=function(){streaming=!streaming;this.textContent='Streaming: '+(streaming?'on':'off');applyVis();draw();};
// main loop
setVB();setStep(D.n_shots);
function loop(){for(let k=0;k<2;k++)tick();draw();requestAnimationFrame(loop);}
loop();
</script>
"""
