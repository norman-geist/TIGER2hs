#!/usr/bin/env python3
"""Generate a self-contained, offline SWARM pool viewer from a NAMD log.

Usage: python3 swarm_pool_view.py [output | file | "output/job%i.out"] [--open]
Python 3.8+; standard library only. No installation or network access required.
Snapshots are delimited by SWARM_POOL headers (or SWARM_POOL_META lines).
SWARM_POOL columns 1 and 4 (zero-indexed, including the tag) supply cycle
and CV. Job files are read in numeric order; playback follows log order.
If a snapshot reports multiple CVs, optional SWARM_POOL_META can identify
which was active. Otherwise the sole reported CV is used directly.
Draw probabilities use reported snapshot weights. Recorded best scores
come from NEW_SEED / OLD_SEED and may reset. Score direction is inferred
when possible and can be adjusted separately for each CV in the viewer.
"""

import argparse
import gzip
import glob
import html
import json
import math
import re
import sys
import webbrowser
from pathlib import Path


def parse_log(lines):
    snapshots, warnings, cv_names = [], [], set()
    pending_meta, context, offers, records, current = None, {}, {}, {}, None
    source = None

    def warn(line, message):
        warnings.append("Line {}: {}".format(line, message))

    def start(line):
        nonlocal current, pending_meta
        current = {
            "line": line, "cycle": None, "run": context.get("run"),
            "cv_id": context.get("cv_id"), "active": None, "nodes": {},
            "records": dict(records), "offers": dict(offers),
            "source": source,
        }
        if pending_meta:
            current["cycle"] = pending_meta.get("cycle")
            current["active"] = pending_meta.get("active")
            if "run" in pending_meta:
                current["run"] = pending_meta["run"]
        pending_meta = None

    def finish():
        nonlocal current
        if current and current["nodes"]:
            reported = {cv for n in current["nodes"].values() for cv in n["cvs"]}
            if current["active"] is None and len(reported) == 1:
                current["active"] = next(iter(reported))
            current["nodes"] = sorted(current["nodes"].values(), key=lambda n: n["id"])
            snapshots.append(current)
        current = None

    for line_no, line in enumerate(lines, 1):
        parts = line.split()
        if not parts:
            continue
        tag = parts[0]
        if tag == "SWARM_VIEW_SOURCE":
            finish()
            source = json.loads(line.split(None, 1)[1])
            context, offers, pending_meta = {}, {}, None
        elif tag == "SWARM_CYCLE":
            finish()
            offers = {}
            context = {}
            try:
                context = {"cv_id": parts[1], "cycle": int(parts[2])}
                if "RUN" in parts:
                    context["run"] = int(parts[parts.index("RUN") + 1])
            except (ValueError, IndexError):
                warn(line_no, "malformed SWARM_CYCLE; run metadata unavailable")
        elif tag == "SWARM_POOL_META":
            finish()
            #The logged Tcl list normally consists of plain CV identifiers.
            tokens = re.findall(r'\{([^{}]*)\}|"([^"]*)"|(\S+)', line.strip())
            tokens = [next((v for v in token if v != ""), "") for token in tokens]
            fields = dict(zip(tokens[1::2], tokens[2::2]))
            pending_meta = {"active": fields.get("CV")}
            try:
                if "Cycle" in fields:
                    pending_meta["cycle"] = int(fields["Cycle"])
                if "RUN" in fields or "Run" in fields:
                    pending_meta["run"] = int(fields.get("RUN", fields.get("Run")))
            except ValueError:
                warn(line_no, "invalid integer in SWARM_POOL_META")
        elif tag in ("NEW_SEED", "OLD_SEED"):
            try:
                value = float(parts[4])
                if not math.isfinite(value):
                    raise ValueError()
                records[parts[1]] = value
            except (ValueError, IndexError):
                warn(line_no, "invalid seed record")
        elif tag == "SEED_OFFER":
            try:
                value = float(parts[4])
                if not math.isfinite(value):
                    raise ValueError()
                offers[int(parts[3])] = value
            except (ValueError, IndexError):
                warn(line_no, "invalid replica offer")
        elif tag == "SWARM_POOL":
            if len(parts) > 1 and parts[1] == "Cycle":
                finish()
                start(line_no)
                if "Weight" not in parts:
                    warn(line_no, "Weight column absent; weights will be unavailable")
                continue
            if current is None:
                start(line_no)
            try:
                if len(parts) not in (8, 9):
                    raise ValueError("expected 7 or 8 columns after SWARM_POOL")
                cycle, node_id, parent = map(int, parts[1:4])
                cv = parts[4]
                score, attempted, successful = map(float, parts[5:8])
                weight = float(parts[8]) if len(parts) == 9 else None
                numbers = [score, attempted, successful] + ([] if weight is None else [weight])
                if not all(math.isfinite(v) for v in numbers):
                    raise ValueError("non-finite numeric value")
                if attempted < 0 or successful < 0 or (weight is not None and weight < 0):
                    raise ValueError("negative counter or weight")
                if current["cycle"] is not None and current["cycle"] != cycle:
                    #Also support reports that omit repeated headers.
                    finish()
                    start(line_no)
                current["cycle"] = cycle
                cv_names.add(cv)
                node = current["nodes"].setdefault(node_id, {"id": node_id, "parent": parent, "cvs": {}})
                if node["parent"] != parent:
                    warn(line_no, "conflicting parents for node {}".format(node_id))
                if cv in node["cvs"]:
                    warn(line_no, "duplicate row for node {}, CV {}".format(node_id, cv))
                node["cvs"][cv] = {"score": score, "attempts": attempted, "successes": successful, "weight": weight}
            except (ValueError, IndexError) as exc:
                warn(line_no, "invalid SWARM_POOL row ({})".format(exc))
    finish()

    cvs = sorted(cv_names)
    if not snapshots:
        raise ValueError("No valid SWARM_POOL snapshots found in the input.")

    #Associate numeric SWARM CV IDs with the CV reported in each snapshot.
    cv_map = {}
    for s in snapshots:
        if s["active"] and s["cv_id"] is not None:
            old = cv_map.get(s["cv_id"])
            if old and old != s["active"]:
                warn(s["line"], "inconsistent CV ID mapping; seed record association omitted")
                cv_map[s["cv_id"]] = False
            elif old is not False:
                cv_map[s["cv_id"]] = s["active"]
    if len(cvs) == 1:
        for s in snapshots:
            s["active"] = s["active"] or cvs[0]
            for cv_id in s["records"]:
                cv_map.setdefault(cv_id, cvs[0])
    elif any(s["active"] is None for s in snapshots):
        warnings.append("Some snapshots report several CVs in column 4. Those rows identify each score's CV, but not which CV was active. SWARM_POOL_META is optional for identifying the active CV in this format.")

    for index, s in enumerate(snapshots):
        s["record"] = {cv_map[k]: v for k, v in s.pop("records").items() if cv_map.get(k)}
        s.pop("cv_id", None)
        nodes = {n["id"]: n for n in s["nodes"]}
        snapshot_cvs = {cv for n in s["nodes"] for cv in n["cvs"]}
        prefix = "Snapshot {} (cycle {})".format(index + 1, s["cycle"])
        if 0 not in nodes:
            warnings.append(prefix + ": virtual root 0 is missing")
        elif nodes[0]["parent"] != -1:
            warnings.append(prefix + ": root parent is not -1")
        for node in s["nodes"]:
            if node["id"] and node["parent"] not in nodes:
                warnings.append(prefix + ": node {} has missing parent {}".format(node["id"], node["parent"]))
            missing = snapshot_cvs - set(node["cvs"])
            if missing:
                warnings.append(prefix + ": node {} is missing CVs {}".format(node["id"], ", ".join(sorted(missing))))
            seen, cursor = set(), node["id"]
            while cursor in nodes and cursor != 0:
                if cursor in seen:
                    warnings.append(prefix + ": parent cycle involving node {}".format(cursor))
                    break
                seen.add(cursor)
                cursor = nodes[cursor]["parent"]
            if node["id"] == 0:
                for cv, row in node["cvs"].items():
                    if row["weight"]:
                        warnings.append(prefix + ": root has nonzero weight for " + cv)
        if s["active"] and s["active"] not in cvs:
            warnings.append(prefix + ": active CV has no reported rows")
    directions = {}
    for cv in cvs:
        values = [s["record"][cv] for s in snapshots if cv in s["record"]]
        changes = [b - a for a, b in zip(values, values[1:]) if b != a]
        directions[cv] = "max" if changes and all(v > 0 for v in changes) else "min"
    return {"snapshots": snapshots, "cvs": cvs, "warnings": list(dict.fromkeys(warnings)), "directions": directions}


TEMPLATE = r'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>SWARM pool — __TITLE__</title>
<style>
:root {color-scheme:light dark;--bg:#f8fafc;--fg:#17212e;--panel:#fff;--line:#cbd5e1;--muted:#536173;--series:#176b9b;--other:#bb6127;--tip:#8c3455;--highlight:#e6f2fb}
@media(prefers-color-scheme:dark){:root{--bg:#111820;--fg:#e4eaf2;--panel:#1c2632;--line:#495667;--muted:#aebaca;--series:#6bb9e5;--other:#efaa70;--tip:#ef91b6;--highlight:#263f52}}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--fg);font:15px/1.5 system-ui,sans-serif}main{max-width:1250px;margin:auto;padding:24px}h1{font-size:24px;margin:0 0 8px}h2{font-size:18px;margin:24px 0 10px}.muted{color:var(--muted)}button,select{font:inherit;color:var(--fg);background:var(--panel);border:1px solid var(--line);border-radius:5px;padding:6px 10px}button{cursor:pointer}button:disabled{opacity:.45;cursor:default}button:focus-visible,select:focus-visible,input:focus-visible{outline:2px solid var(--series);outline-offset:3px}.controls{display:flex;gap:12px;align-items:center;flex-wrap:wrap;margin:16px 0}.controls label{display:flex;align-items:center;gap:7px}input[type=range]{flex:1;min-width:130px;accent-color:var(--series)}#status{font-variant-numeric:tabular-nums}#graph-wrap{overflow:auto;max-height:640px;border:1px solid var(--line);margin-top:12px;background:var(--panel)}#graph{display:block;width:100%;min-width:520px}svg text{fill:var(--fg);font:12px system-ui,sans-serif}#selected{padding:10px 0;min-height:42px;font-variant-numeric:tabular-nums}.legend{display:flex;flex-wrap:wrap;gap:16px;font-size:13px;color:var(--muted)}.legend .one{color:var(--series)}.legend .two{color:var(--other)}.table-wrap{overflow:auto}table{border-collapse:collapse;width:100%;font-variant-numeric:tabular-nums}th,td{text-align:right;padding:7px 10px;border-bottom:1px solid var(--line);white-space:nowrap}th:first-child,td:first-child{text-align:left}tbody tr{cursor:pointer}tbody tr:hover,tbody tr.selected{background:var(--highlight)}#trend{width:100%;display:block}details{margin-top:18px}#warnings{overflow-wrap:anywhere;padding-left:22px}.node{cursor:pointer}.node:focus{outline:none}.node:focus circle{stroke:var(--fg);stroke-width:3}#snapshot-label{min-width:120px}.empty{fill:var(--muted)}@media(max-width:600px){main{padding:14px}h1{font-size:21px}.controls{gap:8px}#status{font-size:14px}}
</style></head><body><main>
<h1>SWARM pool development</h1><div class="muted">__TITLE__</div>
<div class="controls"><button id="prev" type="button">Previous</button><button id="play" type="button">Play</button><button id="next" type="button">Next</button><label>Delay (s) <input id="delay" type="number" required min="0.05" max="60" step="any" value="0.9" style="width:90px;font:inherit;color:var(--fg);background:var(--panel);border:1px solid var(--line);border-radius:5px;padding:6px"></label><label for="cycle" id="snapshot-label">Snapshot <span id="position"></span></label><input type="range" id="cycle" min="0" value="0" aria-label="Snapshot"></div>
<div class="controls"><label>CV <select id="cv"></select></label><label>Best score <select id="direction"><option value="min">Lower is better</option><option value="max">Higher is better</option></select></label><label>Node color <select id="color"><option value="id">Node ID</option><option value="score">CV score</option></select></label><label>Node size <select id="size"><option value="weight">Weight</option><option value="probability">Draw probability</option><option value="uniform">Uniform</option></select></label><label>Edge width <select id="edges"><option value="uniform">Uniform</option><option value="weight">Smaller endpoint weight</option></select></label></div>
<div id="status" aria-live="polite"></div>
<div id="graph-wrap"><svg id="graph" role="img" aria-label="Pool parent-child graph"></svg></div>
<div class="legend"><span id="color-legend">Color identifies node ID; neighboring IDs use contrasting hues</span><span>Arrow = parent → child</span><span id="edge-legend">Edges have uniform width</span><span>Root 0 is virtual</span></div>
<div id="selected" aria-live="polite">Select a node to inspect its counters.</div>
<div class="table-wrap"><table><thead><tr><th>Node</th><th>Parent</th><th>Score</th><th>Attempts</th><th>Successes</th><th>Weight</th><th>Draw probability</th></tr></thead><tbody id="rows"></tbody></table></div>
<h2>Best score over snapshots</h2><div class="legend"><span class="one">━ Recorded swarming best</span><span class="two">━ Best retained in pool</span></div>
<svg id="trend" role="img" aria-label="Best scores over snapshots"></svg>
<div class="muted" style="font-size:13px">Probabilities use the weights at the reported snapshot.</div>
<details id="validation"><summary id="validation-title"></summary><ul id="warnings"></ul></details>
</main><script id="pool-data" type="application/json">__DATA__</script><script>
'use strict';
const report=JSON.parse(document.getElementById('pool-data').textContent),data=report.snapshots;
const byId=id=>document.getElementById(id),slider=byId('cycle'),cvSelect=byId('cv'),direction=byId('direction'),sizeSelect=byId('size'),delayInput=byId('delay'),edgeSelect=byId('edges'),colorSelect=byId('color');
const colorTheme=window.matchMedia?window.matchMedia('(prefers-color-scheme: dark)'):null;
function idColor(id){return `hsl(${((id*137.508)%360).toFixed(3)},72%,${colorTheme&&colorTheme.matches?64:46}%)`}
delayInput.value=report.delay??0.9;
let selected=null,timer=null,lastCv=null,trendCacheKey=null;
const cvDirections={...report.directions},scoreRanges=Object.create(null),histories=Object.create(null);
for(const cv of report.cvs)histories[cv]=[];
//Compute score ranges and pool extrema once; redraws inspect only current nodes.
for(let i=0;i<data.length;i++){
const s=data[i],extrema=Object.create(null);
for(const n of s.nodes){if(!n.id)continue;for(const [cv,row] of Object.entries(n.cvs)){
const range=scoreRanges[cv]||(scoreRanges[cv]={min:Infinity,max:-Infinity});
range.min=Math.min(range.min,row.score);range.max=Math.max(range.max,row.score);
const pair=extrema[cv]||(extrema[cv]={min:Infinity,max:-Infinity});
pair.min=Math.min(pair.min,row.score);pair.max=Math.max(pair.max,row.score);
}}
for(const cv of report.cvs)histories[cv].push({x:i+1,min:extrema[cv]?.min??null,max:extrema[cv]?.max??null,record:s.record[cv]??null,s});
}
slider.max=data.length-1;slider.value=data.length-1;
function option(value,label){const o=document.createElement('option');o.value=value;o.textContent=label;cvSelect.append(o)}
if(report.cvs.length>1 && data.some(s=>s.active))option('__active__','Follow active CV');
for(const cv of report.cvs)option(cv,cv);
if(report.cvs.length>1 && !data.some(s=>s.active))cvSelect.value=report.cvs[0];
byId('validation-title').textContent=report.warnings.length?`${report.warnings.length} validation warning(s)`:'Validation: no structural warnings';
for(const w of report.warnings){const li=document.createElement('li');li.textContent=w;byId('warnings').append(li)}
const fmt=x=>x===null||x===undefined?'—':Number(x.toPrecision(6)).toString();
const NS='http://www.w3.org/2000/svg';
function el(tag,attrs={},text){const n=document.createElementNS(NS,tag);for(const [k,v] of Object.entries(attrs))n.setAttribute(k,v);if(text!==undefined)n.textContent=text;return n}
function activeCv(s){return cvSelect.value==='__active__'?s.active:cvSelect.value}
function stop(){if(timer)clearInterval(timer);timer=null;byId('play').textContent='Play'}
function schedule(){if(timer)clearInterval(timer);timer=setInterval(()=>{if(+slider.value>=data.length-1){stop();return}slider.value=+slider.value+1;draw();if(+slider.value===data.length-1)stop()},Number(delayInput.value)*1000)}
function selectNode(id){selected=id;draw()}
function drawGraph(s,cv,total){
const svg=byId('graph');svg.replaceChildren();
const nodes=s.nodes.map(n=>({...n,children:[],row:n.cvs[cv]})),map=new Map(nodes.map(n=>[n.id,n]));
const visited=new Set(),roots=[];
for(const n of nodes){if(n.parent!==n.id&&map.has(n.parent))map.get(n.parent).children.push(n);else roots.push(n)}
//Guard traversal so malformed parent cycles remain inspectable.
let leaf=0;function place(n,depth){if(visited.has(n.id))return;visited.add(n.id);n.depth=depth;const children=n.children.filter(c=>!visited.has(c.id));children.forEach(c=>place(c,depth+1));const positioned=children.filter(c=>c.y!==undefined);n.y=positioned.length?positioned.reduce((v,c)=>v+c.y,0)/positioned.length:55+leaf++*72;}
roots.sort((a,b)=>a.id-b.id).forEach(n=>place(n,0));
for(const n of nodes)if(!visited.has(n.id))place(n,0);
const depth=Math.max(1,...nodes.map(n=>n.depth)),width=Math.max(520,byId('graph-wrap').clientWidth,depth*160+140),height=Math.max(230,leaf*72+40);
svg.setAttribute('viewBox',`0 0 ${width} ${height}`);svg.style.width=width+'px';svg.style.height=height+'px';
svg.append(el('title',{},`Snapshot ${+slider.value+1}, CV ${cv||'unknown'}, ${nodes.length-Number(map.has(0))} real nodes`));
for(const n of nodes)n.x=65+n.depth*(width-140)/depth;
const radius=n=>{if(!n.id)return 10;if(sizeSelect.value==='uniform'||!n.row)return 14;const w=n.row.weight;if(w===null)return 14;const v=sizeSelect.value==='probability'?(total>0?w/total:0)*nodes.length:w;return Math.min(27,Math.sqrt(130+170*Math.max(0,v)))};
const defs=el('defs'),marker=el('marker',{id:'arrow',viewBox:'0 -3 6 6',refX:0,markerWidth:4,markerHeight:4,markerUnits:'strokeWidth',orient:'auto'});marker.append(el('path',{d:'M0,-3L6,0L0,3',fill:'var(--muted)'}));defs.append(marker);svg.append(defs);
for(const n of nodes){const p=map.get(n.parent);if(!p||p===n)continue;const dx=n.x-p.x,dy=n.y-p.y,len=Math.hypot(dx,dy);if(!len)continue;
const weighted=edgeSelect.value==='weight'&&p.id!==0&&p.row&&n.row&&p.row.weight!==null&&n.row.weight!==null;
const strength=weighted?Math.min(p.row.weight,n.row.weight):null;
//Fixed compressed scale, unchanged across snapshots and CV rotations.
const thickness=weighted?1+4*Math.sqrt(strength/(1+strength)):1.5;
//Anchor the arrow at its base; the shaft ends before the triangular head.
const endGap=radius(n)+5+4*thickness;
const line=el('line',{'data-parent':p.id,'data-child':n.id,x1:p.x+dx*(radius(p)+3)/len,y1:p.y+dy*(radius(p)+3)/len,x2:n.x-dx*endGap/len,y2:n.y-dy*endGap/len,stroke:'var(--muted)','stroke-width':thickness,'stroke-linecap':'butt','marker-end':'url(#arrow)'});
line.append(el('title',{},weighted?`Node ${p.id} → ${n.id}; smaller endpoint weight ${fmt(strength)}`:`Node ${p.id} → ${n.id}; uniform width`));svg.append(line);}
const {min:lo,max:hi}=scoreRanges[cv]||{min:Infinity,max:-Infinity};
for(const n of nodes){const row=n.row;let fraction=row&&hi>lo?(row.score-lo)/(hi-lo):.5;if(direction.value==='min')fraction=1-fraction;
const g=el('g',{transform:`translate(${n.x},${n.y})`,class:'node',tabindex:0,role:'button','aria-label':`Node ${n.id}, score ${row?fmt(row.score):'unavailable'}`});
const fill=!n.id?'var(--line)':colorSelect.value==='id'?idColor(n.id):row?'var(--series)':'var(--line)';
g.append(el('circle',{r:radius(n),fill,'fill-opacity':colorSelect.value==='score'&&n.id&&row?0.25+0.7*fraction:1,stroke:n.id===selected?'var(--fg)':'none','stroke-width':2}));
g.append(el('text',{'text-anchor':'middle',y:-radius(n)-7},`Node ${n.id}`));g.append(el('text',{'text-anchor':'middle',y:radius(n)+17},n.id?(row?fmt(row.score):'missing CV'):'virtual'));
g.append(el('title',{},`Node ${n.id}; parent ${n.parent}; attempts ${row?fmt(row.attempts):'—'}; successes ${row?fmt(row.successes):'—'}; weight ${row?fmt(row.weight):'—'}`));
g.addEventListener('click',()=>selectNode(n.id));g.addEventListener('keydown',e=>{if(e.key==='Enter'||e.key===' '){e.preventDefault();selectNode(n.id)}});svg.append(g);}
}
function drawTrend(cv){
const svg=byId('trend'),width=Math.max(300,svg.clientWidth),height=265,m={left:66,right:24,top:22,bottom:52};
const x=v=>m.left+(v-1)/Math.max(1,data.length-1)*(width-m.left-m.right);
const key=JSON.stringify([cv,direction.value,width]),current=svg.querySelector('[data-current-snapshot]');
if(trendCacheKey===key&&current){current.setAttribute('x1',x(+slider.value+1));current.setAttribute('x2',x(+slider.value+1));return}
trendCacheKey=key;svg.replaceChildren();svg.setAttribute('viewBox',`0 0 ${width} ${height}`);
const points=(histories[cv]||[]).map(p=>({...p,pool:direction.value==='min'?p.min:p.max}));
let lo=Infinity,hi=-Infinity;for(const p of points)for(const v of [p.pool,p.record])if(v!==null){lo=Math.min(lo,v);hi=Math.max(hi,v)}
if(!Number.isFinite(lo)){svg.append(el('text',{x:width/2,y:80,'text-anchor':'middle'},'No scores available for this CV'));return}
const pad=(hi-lo)*.1||Math.max(1,Math.abs(lo)*.05);lo-=pad;hi+=pad;
const y=v=>height-m.bottom-(v-lo)/(hi-lo)*(height-m.top-m.bottom);
svg.append(el('rect',{x:m.left,y:m.top,width:width-m.left-m.right,height:height-m.top-m.bottom,fill:'none',stroke:'var(--line)'}));
for(let i=0;i<5;i++){const v=lo+(hi-lo)*i/4,yy=y(v);svg.append(el('line',{x1:m.left,x2:width-m.right,y1:yy,y2:yy,stroke:'var(--line)','stroke-opacity':.35}));svg.append(el('text',{x:m.left-9,y:yy+4,'text-anchor':'end'},fmt(v)));}
const count=width<450?3:6;const ticks=[...new Set(Array.from({length:count},(_,i)=>Math.round(1+(data.length-1)*i/(count-1))))];for(const tick of ticks)svg.append(el('text',{x:x(tick),y:height-m.bottom+21,'text-anchor':'middle'},tick));
svg.append(el('text',{x:(m.left+width-m.right)/2,y:height-8,'text-anchor':'middle'},'Snapshot in log order'));
svg.append(el('text',{transform:`translate(17,${(m.top+height-m.bottom)/2}) rotate(-90)`,'text-anchor':'middle'},`CV score`));
for(const [key,color] of [['record','var(--series)'],['pool','var(--other)']]){let d='',prev=null;for(const p of points){if(p[key]===null){prev=null;continue}d+=prev?`H${x(p.x)}V${y(p[key])}`:`M${x(p.x)},${y(p[key])}`;prev=p;}svg.append(el('path',{d,fill:'none',stroke:color,'stroke-width':2}));
for(const p of points)if(p[key]!==null){const g=el('g',{role:'button',tabindex:0,'aria-label':`Snapshot ${p.x}, ${key} score ${p[key]}`});g.append(el('circle',{cx:x(p.x),cy:y(p[key]),r:10,fill:'transparent'}));g.append(el('circle',{cx:x(p.x),cy:y(p[key]),r:3,fill:color}));g.append(el('title',{},`Snapshot ${p.x}; cycle ${p.s.cycle}; ${key==='record'?'recorded best':'pool best'} ${p[key]}`));const choose=()=>{slider.value=p.x-1;draw()};g.addEventListener('click',choose);g.addEventListener('keydown',e=>{if(e.key==='Enter'){choose()}});svg.append(g);}}
svg.append(el('line',{'data-current-snapshot':'',x1:x(+slider.value+1),x2:x(+slider.value+1),y1:m.top,y2:height-m.bottom,stroke:'var(--fg)','stroke-opacity':.4}));
}
function draw(){const i=+slider.value,s=data[i],cv=activeCv(s);if(cv!==lastCv){direction.value=cvDirections[cv]||'min';lastCv=cv}
byId('color-legend').textContent=colorSelect.value==='id'?'Color identifies node ID; neighboring IDs use contrasting hues':'Darker fill = better score on a fixed scale per CV';
byId('edge-legend').textContent=edgeSelect.value==='weight'?'Thicker edges = both endpoints more strongly weighted (compressed scale); root edges uniform':'Edges have uniform width';
byId('position').textContent=`${i+1} / ${data.length}`;byId('prev').disabled=i===0;byId('next').disabled=i===data.length-1;
byId('status').textContent=`${s.source?s.source+' · ':''}Cycle ${s.cycle} · TIGER2 run ${s.run??'unknown'} · CV ${cv||'unknown'} · Active CV ${s.active||'unknown'} · ${s.nodes.filter(n=>n.id).length} real nodes`;
let total=0,complete=true;for(const n of s.nodes){if(!n.id)continue;const row=n.cvs[cv];if(!row||row.weight===null){complete=false;continue}total+=row.weight}
if(!complete)total=null;
const tbody=byId('rows');tbody.replaceChildren();for(const n of s.nodes){const r=n.cvs[cv],tr=document.createElement('tr');tr.className=n.id===selected?'selected':'';const cells=[n.id,n.parent,n.id?(r?fmt(r.score):'—'):'virtual',r?fmt(r.attempts):'—',r?fmt(r.successes):'—',r?fmt(r.weight):'—',r&&total>0?(100*r.weight/total).toFixed(2)+'%':'—'];for(const value of cells){const td=document.createElement('td');td.textContent=value;tr.append(td)}tr.addEventListener('click',()=>selectNode(n.id));tbody.append(tr)}
const node=s.nodes.find(n=>n.id===selected),r=node?.cvs[cv];byId('selected').textContent=node?`Node ${node.id} · Parent ${node.parent} · Attempts ${r?fmt(r.attempts):'—'} · Successes ${r?fmt(r.successes):'—'} · Weight ${r?fmt(r.weight):'—'}`:selected!==null?`Node ${selected} is absent from this snapshot.`:'Select a node to inspect its counters.';
drawGraph(s,cv,total);drawTrend(cv);
}
byId('prev').addEventListener('click',()=>{stop();slider.value=Math.max(0,+slider.value-1);draw()});byId('next').addEventListener('click',()=>{stop();slider.value=Math.min(data.length-1,+slider.value+1);draw()});slider.addEventListener('input',()=>{stop();draw()});
byId('play').addEventListener('click',()=>{if(timer){stop();return}if(!delayInput.checkValidity()){delayInput.reportValidity();return}if(+slider.value===data.length-1){slider.value=0;draw()}byId('play').textContent='Pause';schedule()});
delayInput.addEventListener('input',()=>{if(delayInput.checkValidity()&&timer)schedule()});
delayInput.addEventListener('change',()=>{if(!delayInput.checkValidity()){stop();delayInput.reportValidity()}});
cvSelect.addEventListener('change',draw);direction.addEventListener('change',()=>{const cv=activeCv(data[+slider.value]);if(cv)cvDirections[cv]=direction.value;draw()});sizeSelect.addEventListener('change',draw);
edgeSelect.addEventListener('change',draw);colorSelect.addEventListener('change',draw);
if(colorTheme&&colorTheme.addEventListener)colorTheme.addEventListener('change',draw);
window.addEventListener('resize',draw);draw();
</script></body></html>'''


def build_report(report, title):
    payload = json.dumps(report, separators=(",", ":"), ensure_ascii=True, allow_nan=False)
    payload = payload.replace("<", "\\u003c").replace(">", "\\u003e").replace("&", "\\u0026")
    return TEMPLATE.replace("__TITLE__", html.escape(title)).replace("__DATA__", payload)


def resolve_inputs(inputs):
    if inputs == ["-"]:
        return []
    if "-" in inputs:
        raise ValueError("Standard input cannot be combined with file inputs.")

    def job_order(path):
        match = re.fullmatch(r"job(\d+)\.out(?:\.gz)?", path.name)
        return (str(path.parent), 0, int(match.group(1)), path.name) if match else (str(path.parent), 1, 0, path.name)

    paths = []
    for name in inputs:
        path = Path(name)
        if path.is_dir():
            matches = [p for p in path.iterdir() if p.is_file() and re.fullmatch(r"job\d+\.out(?:\.gz)?", p.name)]
        elif "%i" in name or glob.has_magic(name):
            matches = [Path(p) for p in glob.glob(name.replace("%i", "*")) if Path(p).is_file()]
        elif path.is_file():
            matches = [path]
        else:
            raise ValueError("Input not found: {}".format(name))
        if not matches:
            raise ValueError("No matching job output files: {}".format(name))
        paths.extend(matches)
    unique = {p.resolve(): p for p in paths}
    return sorted(unique.values(), key=job_order)


def iter_logs(paths):
    for path in paths:
        yield "SWARM_VIEW_SOURCE " + json.dumps(path.name)
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rt", encoding="utf-8-sig", errors="replace") as stream:
            yield from stream


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("logs", nargs="*", help="directory (default: output), file(s), quoted glob, job%%i.out pattern, or - for stdin")
    parser.add_argument("--output", "-o", type=Path, help="HTML destination (default: directory/swarm.pool.html or <file stem>.pool.html)")
    parser.add_argument("--open", action="store_true", help="open the generated report in the default browser")
    parser.add_argument("--delay", type=float, default=0.9, help="initial playback delay in seconds, adjustable in the report (0.05–60; default: 0.9)")
    args = parser.parse_args(argv)
    try:
        if not math.isfinite(args.delay) or not 0.05 <= args.delay <= 60:
            raise ValueError("Playback delay must be between 0.05 and 60 seconds.")
        inputs = args.logs or ["output"]
        paths = resolve_inputs(inputs)
        output = args.output
        if output is None:
            if len(inputs) == 1 and Path(inputs[0]).is_dir():
                output = Path(inputs[0]) / "swarm.pool.html"
            elif len(paths) == 1:
                source = paths[0]
                name = source.stem if source.suffix == ".gz" else source.name
                output = source.parent / (Path(name).stem + ".pool.html")
            else:
                output = Path("swarm.pool.html")
        if any(output.resolve() == path.resolve() for path in paths):
            raise ValueError("Output must differ from every input log.")
        if inputs == ["-"]:
            report = parse_log(sys.stdin)
        else:
            report = parse_log(iter_logs(paths))
        title = paths[0].name if len(paths) == 1 else "{} job files".format(len(paths)) if paths else "Standard input"
        report["delay"] = args.delay
        output.write_text(build_report(report, title), encoding="utf-8")
    except (OSError, ValueError, EOFError) as exc:
        parser.exit(1, "Error: {}\n".format(exc))
    print("Wrote {} ({} snapshots; CVs: {})".format(output, len(report["snapshots"]), ", ".join(report["cvs"])))
    for warning in report["warnings"][:20]:
        print("Warning: " + warning, file=sys.stderr)
    if len(report["warnings"]) > 20:
        print("Warning: {} additional warnings are listed in the report.".format(len(report["warnings"]) - 20), file=sys.stderr)
    if args.open:
        try:
            if not webbrowser.open(output.resolve().as_uri()):
                print("Warning: browser could not be opened; open the HTML file manually.", file=sys.stderr)
        except webbrowser.Error as exc:
            print("Warning: {}. Open the HTML file manually.".format(exc), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
