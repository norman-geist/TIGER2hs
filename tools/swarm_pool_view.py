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
Click nodes or table rows to toggle persistent tracking outlines; colors and legend also export in GIFs.
Diagnostic plots show tracked-node probabilities, frontier probability and pool topology.
Graph views fit the page width or screen without internal scrollbars.
Pool CV spread shows the retained score range and median, excluding virtual root 0.
Replica CV spread uses SEED_OFFER scores before reseeding, showing individual
replicas, range and median alongside the best retained pool score.
All time-series plots support drag-to-select rectangle zoom, directional
pan buttons and reset; zoom is retained per CV during playback. Offers
are matched by CV and sampling event; legacy cycle indices are also accepted.
Render GIF exports the full graph sequence using the current display settings,
playback delay, and smooth-motion preference; rendering can be cancelled.
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
        nonlocal current, pending_meta, offers
        current = {
            "line": line, "cycle": None, "run": context.get("run"),
            "cv_id": context.get("cv_id"), "active": None, "nodes": {},
            "records": dict(records), "offers": offers,
            "source": source,
        }
        if pending_meta:
            current["cycle"] = pending_meta.get("cycle")
            current["active"] = pending_meta.get("active")
            if "run" in pending_meta:
                current["run"] = pending_meta["run"]
        pending_meta = None
        #Offers describe this sampling event; never carry them into another snapshot.
        offers = {}
        if current["cv_id"] is None and len(current["offers"]) == 1:
            current["cv_id"] = next(iter(current["offers"]))

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
                cv_id, run, replica_id = parts[1], int(parts[2]), int(parts[3])
                if replica_id < 0:
                    raise ValueError()
                batch = offers.get(cv_id)
                if batch is None or batch["run"] != run:
                    batch = offers[cv_id] = {"run": run, "scores": {}}
                batch["scores"][replica_id] = value
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
            for cv_id in set(s["records"]) | set(s["offers"]):
                cv_map.setdefault(cv_id, cvs[0])
    elif any(s["active"] is None for s in snapshots):
        warnings.append("Some snapshots report several CVs in column 4. Those rows identify each score's CV, but not which CV was active. SWARM_POOL_META is optional for identifying the active CV in this format.")

    for index, s in enumerate(snapshots):
        s["record"] = {cv_map[k]: v for k, v in s.pop("records").items() if cv_map.get(k)}
        s["replicas"] = {}
        for cv_id, batch in s.pop("offers").items():
            cv = cv_map.get(cv_id, cv_id if cv_id in cvs else None)
            if cv and (s["run"] is None or batch["run"] in (s["run"], s["cycle"])):
                s["replicas"][cv] = batch
        if s["run"] is None and s["active"] in s["replicas"]:
            s["run"] = s["replicas"][s["active"]]["run"]
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
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--fg);font:15px/1.5 system-ui,sans-serif}main{max-width:1250px;margin:auto;padding:24px}h1{font-size:24px;margin:0 0 8px}h2{font-size:18px;margin:24px 0 10px}.muted{color:var(--muted)}button,select{font:inherit;color:var(--fg);background:var(--panel);border:1px solid var(--line);border-radius:5px;padding:6px 10px}button{cursor:pointer}button:disabled{opacity:.45;cursor:default}button:focus-visible,select:focus-visible,input:focus-visible{outline:2px solid var(--series);outline-offset:3px}.controls{display:flex;gap:12px;align-items:center;flex-wrap:wrap;margin:16px 0}.controls label{display:flex;align-items:center;gap:7px}input[type=range]{flex:1;min-width:130px;accent-color:var(--series)}input[type=number]{appearance:textfield;-moz-appearance:textfield}input[type=number]::-webkit-inner-spin-button,input[type=number]::-webkit-outer-spin-button{-webkit-appearance:none;margin:0}#status{font-variant-numeric:tabular-nums}#graph-wrap{border:1px solid var(--line);margin-top:12px;background:var(--panel)}#graph{display:block;width:100%}svg text{fill:var(--fg);font:12px system-ui,sans-serif}#selected{padding:10px 0;min-height:42px;font-variant-numeric:tabular-nums}.legend{display:flex;flex-wrap:wrap;gap:16px;font-size:13px;color:var(--muted)}.legend .one{color:var(--series)}.legend .two{color:var(--other)}.table-wrap{overflow:auto}table{border-collapse:collapse;width:100%;font-variant-numeric:tabular-nums}th,td{text-align:right;padding:7px 10px;border-bottom:1px solid var(--line);white-space:nowrap}th:first-child,td:first-child{text-align:left}tbody tr{cursor:pointer}tbody tr:hover,tbody tr.selected{background:var(--highlight)}#trend,#spread,#replica-spread,#tracked-probability,#frontier-probability,#topology{width:100%;display:block}details{margin-top:18px}#warnings{overflow-wrap:anywhere;padding-left:22px}.node{cursor:pointer}.node:focus{outline:none}.node:focus-visible text{font-weight:700}#tracked-panel{position:fixed;right:16px;bottom:16px;z-index:5;padding:12px;background:var(--panel);border:1px solid var(--line);border-radius:8px;box-shadow:0 3px 16px #0002;max-width:260px;max-height:45vh;overflow:auto}#tracked-panel[hidden]{display:none}#tracked-panel strong{display:block}.tracked-item{display:flex;gap:8px;align-items:center}.tracked-item input{width:30px;height:26px;padding:0;border:0;background:none}.tracked-list{display:flex;gap:8px;flex-wrap:wrap;align-items:center;margin:8px 0}.tracked-list button{border-color:#dc2626}.tracked-list button.absent{opacity:.6}#snapshot-label{min-width:120px}.empty{fill:var(--muted)}@media(max-width:600px){main{padding:14px}h1{font-size:21px}.controls{gap:8px}#status{font-size:14px}}

.plot-controls{display:flex;align-items:center;gap:8px;flex-wrap:wrap;margin:8px 0;font-size:12px}.plot-controls button{min-width:34px}.zoom-plot{touch-action:none;cursor:crosshair;user-select:none}.zoom-plot:focus-visible{outline:2px solid var(--series);outline-offset:3px}
.viewer-controls{display:grid;gap:12px;margin:22px 0 16px}.toolbar-panel{padding:16px 18px;background:var(--panel);border:1px solid var(--line);border-radius:10px}.toolbar-heading{font-size:12px;font-weight:650;letter-spacing:.06em;text-transform:uppercase;color:var(--muted);margin-bottom:12px}.playback-row{display:flex;align-items:center;gap:20px;flex-wrap:wrap}.button-group{display:flex;gap:6px;align-items:center}.viewer-controls button{min-height:38px}.viewer-controls #play{min-width:76px;background:var(--series);border-color:var(--series);color:var(--panel);font-weight:600}.delay-control{display:flex;gap:8px;align-items:center}.delay-control label{font-size:13px;color:var(--muted);white-space:nowrap}.delay-control input{width:76px;height:38px;text-align:center;font:inherit;color:var(--fg);background:var(--panel);border:1px solid var(--line);border-radius:5px}.delay-control button{width:34px;padding:6px}.motion-control{display:flex;align-items:center;gap:7px;font-size:14px;white-space:nowrap}.motion-control input{accent-color:var(--series)}#render-gif{margin-left:auto}.timeline-row{display:flex;align-items:center;gap:18px;margin-top:16px;padding-top:14px;border-top:1px solid var(--line)}.timeline-row #snapshot-label{display:flex;gap:8px;align-items:center;font-size:13px;min-width:145px;color:var(--muted)}#position{font-variant-numeric:tabular-nums;color:var(--fg);font-weight:600}.timeline-row #cycle{width:100%;min-width:0;margin:0}.display-grid{display:grid;grid-template-columns:repeat(6,minmax(0,1fr));gap:14px}.display-field{display:flex;flex-direction:column;gap:6px;min-width:0}.display-field>span{font-size:13px;color:var(--muted)}.display-field select{width:100%;min-width:0;height:38px}.tracking-row{display:flex;align-items:center;gap:18px;justify-content:space-between}.tracking-row p{margin:0;font-size:13px;color:var(--muted)}.tracking-row button{white-space:nowrap}#status{font-size:13px;color:var(--muted);padding:0 2px}#export-status:empty{display:none}#export-status{margin-top:6px}
@media(max-width:1050px){.display-grid{grid-template-columns:repeat(3,minmax(0,1fr))}.playback-row{gap:12px}}
@media(max-width:600px){.toolbar-panel{padding:14px}.display-grid{grid-template-columns:repeat(2,minmax(0,1fr));gap:12px}.playback-row{gap:12px 16px}#render-gif{margin-left:0}.timeline-row{flex-direction:column;align-items:stretch;gap:10px}.timeline-row #snapshot-label{justify-content:space-between}.tracking-row{flex-wrap:wrap;gap:10px}.viewer-controls{margin-top:18px}}

</style></head><body><main>
<h1>SWARM pool development</h1><div class="muted">__TITLE__</div>
<section class="viewer-controls" aria-label="Viewer controls">
<div class="toolbar-panel">
<div class="toolbar-heading">Graph display</div>
<div class="display-grid">
<label class="display-field"><span>CV</span><select id="cv"></select></label>
<label class="display-field"><span>Best score</span><select id="direction"><option value="min">Lower is better</option><option value="max">Higher is better</option></select></label>
<label class="display-field"><span>Graph view</span><select id="graph-view"><option value="width">Fit width</option><option value="screen">Fit screen</option></select></label>
<label class="display-field"><span>Node color</span><select id="color"><option value="id">Node ID</option><option value="score">CV score</option></select></label>
<label class="display-field"><span>Node size</span><select id="size"><option value="weight">Weight</option><option value="probability">Draw probability</option><option value="uniform">Uniform</option></select></label>
<label class="display-field"><span>Edge width</span><select id="edges"><option value="uniform">Uniform</option><option value="weight">Endpoint weight</option></select></label>
</div>
</div>
<div class="toolbar-panel">
<div class="toolbar-heading">Playback</div>
<div class="playback-row">
<div class="button-group" role="group" aria-label="Playback navigation"><button id="prev" type="button">Previous</button><button id="play" type="button">Play</button><button id="next" type="button">Next</button></div>
<div class="delay-control"><label for="delay">Delay (s)</label><div class="button-group"><button id="delay-down" type="button" aria-label="Decrease playback delay">−</button><input id="delay" type="number" required min="0.05" max="60" step="any" value="0.9"><button id="delay-up" type="button" aria-label="Increase playback delay">+</button></div></div>
<label class="motion-control"><input id="motion" type="checkbox" checked> Smooth motion</label>
<button id="render-gif" type="button">Render GIF</button>
</div>
<div class="timeline-row"><label for="cycle" id="snapshot-label">Snapshot <span id="position"></span></label><input type="range" id="cycle" min="0" value="0" aria-label="Snapshot"></div>
</div>

</section>
<div id="status" aria-live="polite"></div><div id="export-status" role="status" aria-live="polite" class="muted"></div>
<div id="graph-wrap"><svg id="graph" role="img" aria-label="Pool parent-child graph"></svg></div>
<div class="legend"><span id="color-legend">Color identifies node ID; neighboring IDs use contrasting hues</span><span>Arrow = parent → child</span><span id="edge-legend">Edges have uniform width</span><span>Root 0 is virtual</span></div>


<div class="toolbar-panel"><div class="toolbar-heading">Node tracking</div><div class="tracking-row"><p>Click a node or table row to track it. Choose outline colors in the corner legend; tracking is included in GIF exports.</p><button id="clear-tracked" type="button">Clear tracked nodes</button></div></div>
<div id="selected" aria-live="polite">Select a node to inspect its counters.</div>
<div class="table-wrap"><table><thead><tr><th>Node</th><th>Parent</th><th>Score</th><th>Attempts</th><th>Successes</th><th>Weight</th><th>Draw probability</th></tr></thead><tbody id="rows"></tbody></table></div>
<h2>Pool CV spread</h2>
<div id="spread-summary" aria-live="polite"></div>
<div class="legend"><span class="two">Shaded band: minimum–maximum</span><span class="one">━ Median</span><span>Virtual root excluded; gaps indicate unreported CV scores</span></div>
<svg id="spread" role="img" aria-label="Retained pool CV range and median over snapshots"></svg>
<h2>Replica CV spread</h2>
<div id="replica-spread-summary" aria-live="polite"></div>
<div class="legend"><span class="two">Shaded band: replica minimum–maximum</span><span class="one">━ Replica median</span><span style="color:#8b5cf6">┄ Best retained in pool</span></div>
<svg id="replica-spread" role="img" aria-label="Replica CV scores, range and median compared with the pool frontier over snapshots"></svg>
<h2>Best score over snapshots</h2><div class="legend"><span class="one">━ Recorded swarming best</span><span class="two">━ Best retained in pool</span></div>
<svg id="trend" role="img" aria-label="Best scores over snapshots"></svg>
<div class="muted" style="font-size:13px">Probabilities use the weights at the reported snapshot.</div>
<h2>Tracked-node draw probabilities</h2>
<div class="muted">Colors match the tracking legend. Gaps mean a node is absent or its CV probability is unavailable.</div>
<svg id="tracked-probability" role="img" aria-label="Tracked-node draw probabilities over snapshots"></svg>
<h2>Frontier draw probability</h2>
<div class="controls"><label>Frontier width <input id="frontier-width" type="number" required value="5" min="0" max="100" step="any" style="width:70px"> % of observed CV range</label></div>
<div class="muted">Combined probability of nodes within this distance of the best retained score. The observed range is fixed per CV across the report.</div>
<div id="frontier-summary"></div><svg id="frontier-probability" role="img" aria-label="Combined frontier draw probability over snapshots"></svg>
<h2>Pool topology over time</h2>
<div class="legend"><span class="one">━ Real nodes</span><span class="two">━ Tips</span><span style="color:#8b5cf6">━ Branch points</span></div>
<div id="topology-summary"></div><svg id="topology" role="img" aria-label="Pool size, tip count and branch-point count over snapshots"></svg>
<details id="validation"><summary id="validation-title"></summary><ul id="warnings"></ul></details>
</main><aside id="tracked-panel" hidden aria-label="Tracking legend"><strong>Tracked nodes</strong><div id="tracked-list" class="tracked-list"></div></aside><script id="pool-data" type="application/json">__DATA__</script><script>
'use strict';
const report=JSON.parse(document.getElementById('pool-data').textContent),data=report.snapshots;
const byId=id=>document.getElementById(id),slider=byId('cycle'),cvSelect=byId('cv'),direction=byId('direction'),sizeSelect=byId('size'),delayInput=byId('delay'),edgeSelect=byId('edges'),colorSelect=byId('color');
const colorTheme=window.matchMedia?window.matchMedia('(prefers-color-scheme: dark)'):null;
function idColor(id){return `hsl(${((id*137.508)%360).toFixed(3)},72%,${colorTheme&&colorTheme.matches?64:46}%)`}
delayInput.value=report.delay??0.9;
const trackedNodes=new Set(),trackedColors=new Map(),trackPalette=['#dc2626','#2563eb','#16a34a','#9333ea','#ea580c','#0891b2','#db2777'];
const trackColor=id=>trackedColors.get(id)||'#dc2626';
const zoomPlotIds=['spread','replica-spread','trend','tracked-probability','frontier-probability','topology'];
const spreadViews=new Map(),spreadGeometry=new Map();
let spreadDrag=null,spreadClickBlockedUntil=0;
const topologyHistory=[],replicaHistories=Object.create(null),spreadCaches=Object.create(null);
let selected=null,timer=null,lastCv=null,trendCacheKey=null;
const cvDirections={...report.directions},scoreRanges=Object.create(null),histories=Object.create(null);
for(const cv of report.cvs){histories[cv]=[];replicaHistories[cv]=[]}
//Compute score ranges and pool extrema once; redraws inspect only current nodes.
for(let i=0;i<data.length;i++){
const s=data[i],extrema=Object.create(null),real=s.nodes.filter(n=>n.id),children=new Map(s.nodes.map(n=>[n.id,0]));
for(const n of real)if(children.has(n.parent)&&n.parent!==n.id)children.set(n.parent,children.get(n.parent)+1);
topologyHistory.push({nodes:real.length,tips:real.filter(n=>children.get(n.id)===0).length,branches:real.filter(n=>children.get(n.id)>1).length});
for(const n of s.nodes){if(!n.id)continue;for(const [cv,row] of Object.entries(n.cvs)){
const range=scoreRanges[cv]||(scoreRanges[cv]={min:Infinity,max:-Infinity});
range.min=Math.min(range.min,row.score);range.max=Math.max(range.max,row.score);
const pair=extrema[cv]||(extrema[cv]={min:Infinity,max:-Infinity,values:[]});
pair.values.push(row.score);pair.min=Math.min(pair.min,row.score);pair.max=Math.max(pair.max,row.score);
}}
for(const cv of report.cvs){
const values=extrema[cv]?.values||[];values.sort((a,b)=>a-b);
const count=values.length,middle=Math.floor(count/2),median=count?(count%2?values[middle]:values[middle-1]/2+values[middle]/2):null;
const complete=real.length>0&&real.every(n=>n.cvs[cv]&&n.cvs[cv].weight!==null),total=complete?real.reduce((sum,n)=>sum+n.cvs[cv].weight,0):0;
const probabilities=complete&&total>0?Object.fromEntries(real.map(n=>[n.id,n.cvs[cv].weight/total])):null;
const batch=s.replicas?.[cv],replicas=Object.entries(batch?.scores||{}).map(([id,score])=>({id:Number(id),score})).sort((a,b)=>a.score-b.score||a.id-b.id);
const rc=replicas.length,rm=Math.floor(rc/2);
replicaHistories[cv].push({x:i+1,min:rc?replicas[0].score:null,max:rc?replicas[rc-1].score:null,median:rc?(rc%2?replicas[rm].score:replicas[rm-1].score/2+replicas[rm].score/2):null,count:rc,replicas,run:s.run??batch?.run,s,poolMin:extrema[cv]?.min??null,poolMax:extrema[cv]?.max??null});
histories[cv].push({probabilities,x:i+1,min:extrema[cv]?.min??null,max:extrema[cv]?.max??null,median,count,record:s.record[cv]??null,s});
}
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
function selectNode(id){if(exportRunning)return;selected=id;if(trackedNodes.has(id))trackedNodes.delete(id);else{trackedNodes.add(id);if(!trackedColors.has(id))trackedColors.set(id,trackPalette[(trackedColors.size)%trackPalette.length])}draw()}
byId('clear-tracked').addEventListener('click',()=>{if(exportRunning)return;trackedNodes.clear();selected=null;draw()});
byId('graph-view').addEventListener('change',draw);
const graphNodes=new Map(),graphEdges=new Map();
let graphFrame=null,graphSnapshot=null,graphSize=null,graphTransition=null;
let exportRunning=false,exportAbort=false,exportStepping=false;
const reducedMotion=window.matchMedia?window.matchMedia('(prefers-reduced-motion: reduce)'):null;
const motionInput=byId('motion');
function drawGraph(s,cv,total){
const svg=byId('graph'),snapshot=+slider.value;
const animate=graphSnapshot!==null&&motionInput.checked&&!(reducedMotion&&reducedMotion.matches)&&(snapshot!==graphSnapshot||graphFrame!==null);
if(graphFrame!==null){cancelAnimationFrame(graphFrame);graphFrame=null}
if(!svg.querySelector('[data-node-layer]')){
svg.replaceChildren();svg.append(el('title'));
const defs=el('defs'),marker=el('marker',{id:'arrow',viewBox:'0 -3 6 6',refX:0,markerWidth:4,markerHeight:4,markerUnits:'strokeWidth',orient:'auto'});
marker.append(el('path',{d:'M0,-3L6,0L0,3',fill:'var(--muted)'}));defs.append(marker);svg.append(defs);
svg.append(el('g',{'data-edge-layer':''}),el('g',{'data-node-layer':''}));
}
const edgeLayer=svg.querySelector('[data-edge-layer]'),nodeLayer=svg.querySelector('[data-node-layer]');
const nodes=s.nodes.map(n=>({...n,children:[],row:n.cvs[cv]})),map=new Map(nodes.map(n=>[n.id,n]));
const visited=new Set(),roots=[];
for(const n of nodes){if(n.parent!==n.id&&map.has(n.parent))map.get(n.parent).children.push(n);else roots.push(n)}
let leaf=0;
function place(n,depth){
if(visited.has(n.id))return;visited.add(n.id);n.depth=depth;
const children=n.children.filter(c=>!visited.has(c.id));children.forEach(c=>place(c,depth+1));
const positioned=children.filter(c=>c.y!==undefined);
n.y=positioned.length?positioned.reduce((v,c)=>v+c.y,0)/positioned.length:55+leaf++*88;
}
roots.sort((a,b)=>a.id-b.id).forEach(n=>place(n,0));
for(const n of nodes)if(!visited.has(n.id))place(n,0);
const depth=Math.max(1,...nodes.map(n=>n.depth)),width=Math.max(520,byId('graph-wrap').clientWidth,depth*160+140),height=Math.max(230,leaf*88+40);
svg.querySelector('title').textContent=`Snapshot ${snapshot+1}, CV ${cv||'unknown'}, ${nodes.length-Number(map.has(0))} real nodes`;
for(const n of nodes)n.x=65+n.depth*(width-140)/depth;
const radius=n=>{
if(!n.id)return 10;if(sizeSelect.value==='uniform'||!n.row)return 14;
const w=n.row.weight;if(w===null)return 14;
const v=sizeSelect.value==='probability'?(total>0?w/total:0)*nodes.length:w;
return Math.min(27,Math.sqrt(130+170*Math.max(0,v)));
};
const before=new Map([...graphNodes].map(([id,n])=>[id,{x:n.x,y:n.y,r:n.r,opacity:n.opacity}]));
for(const entry of graphNodes.values())entry.active=false;
const {min:lo,max:hi}=scoreRanges[cv]||{min:Infinity,max:-Infinity};
for(const n of nodes){
let entry=graphNodes.get(n.id);
if(!entry){
//New nodes start at the nearest ancestor present in the preceding frame.
let cursor=n.parent,ancestor=null;const seen=new Set();
while(cursor!==undefined&&!seen.has(cursor)){
seen.add(cursor);if(before.has(cursor)){ancestor=before.get(cursor);break}cursor=map.get(cursor)?.parent;
}
const start=animate&&ancestor?ancestor:{x:n.x,y:n.y};
const g=el('g',{class:'node','data-node-id':n.id,role:'button'});
g.append(el('circle'),el('text',{'text-anchor':'middle'}),el('text',{'text-anchor':'middle'}),el('title'));
g.addEventListener('click',()=>selectNode(n.id));
g.addEventListener('keydown',e=>{if(e.key==='Enter'||e.key===' '){e.preventDefault();selectNode(n.id)}});
nodeLayer.append(g);
entry={g,x:start.x,y:start.y,r:radius(n),opacity:animate?0:1};graphNodes.set(n.id,entry);
}
entry.active=true;entry.node=n;entry.target={x:n.x,y:n.y,r:radius(n),opacity:1};
entry.g.setAttribute('tabindex','0');entry.g.removeAttribute('aria-hidden');
entry.g.setAttribute('aria-label',`Node ${n.id}, score ${n.row?fmt(n.row.score):'unavailable'}`);
entry.g.style.pointerEvents='';
const row=n.row;let fraction=row&&hi>lo?(row.score-lo)/(hi-lo):.5;if(direction.value==='min')fraction=1-fraction;
const fill=!n.id?'var(--line)':colorSelect.value==='id'?idColor(n.id):row?'var(--series)':'var(--line)';
const circle=entry.g.querySelector('circle');circle.setAttribute('fill',fill);
circle.setAttribute('fill-opacity',colorSelect.value==='score'&&n.id&&row?0.25+0.7*fraction:1);
circle.setAttribute('stroke',trackedNodes.has(n.id)?trackColor(n.id):'none');circle.setAttribute('stroke-width',3);circle.setAttribute('vector-effect','non-scaling-stroke');entry.g.setAttribute('aria-pressed',String(trackedNodes.has(n.id)));
const texts=entry.g.querySelectorAll('text');texts[0].textContent=`Node ${n.id}`;texts[1].textContent=n.id?(row?fmt(row.score):'missing CV'):'virtual';
entry.g.querySelector('title').textContent=`Node ${n.id}; parent ${n.parent}; attempts ${row?fmt(row.attempts):'—'}; successes ${row?fmt(row.successes):'—'}; weight ${row?fmt(row.weight):'—'}`;
}
for(const entry of graphNodes.values()){
entry.from={x:entry.x,y:entry.y,r:entry.r,opacity:entry.opacity};
if(!entry.active){entry.target={...entry.from,opacity:0};entry.g.setAttribute('tabindex','-1');entry.g.setAttribute('aria-hidden','true');entry.g.style.pointerEvents='none'}
}
for(const edge of graphEdges.values())edge.active=false;
for(const n of nodes){
const p=map.get(n.parent);if(!p||p===n)continue;
const weighted=edgeSelect.value==='weight'&&p.id!==0&&p.row&&n.row&&p.row.weight!==null&&n.row.weight!==null;
const strength=weighted?Math.min(p.row.weight,n.row.weight):null;
const thickness=weighted?1+4*Math.sqrt(strength/(1+strength)):1.5;
let edge=graphEdges.get(n.id);
if(!edge){
const line=el('line',{stroke:'var(--muted)','stroke-linecap':'butt','marker-end':'url(#arrow)'});line.append(el('title'));edgeLayer.append(line);
edge={line,parent:n.parent,child:n.id,width:thickness,opacity:animate?0:1,attach:null};graphEdges.set(n.id,edge);
}
edge.fromWidth=edge.width;edge.targetWidth=thickness;edge.fromOpacity=edge.opacity;edge.targetOpacity=1;
//A reconnection moves the parent end from its currently displayed attachment.
edge.reconnect=edge.parent!==n.parent||(animate&&edge.reconnect);
edge.fromAttach=edge.attach?{...edge.attach}:null;
edge.parent=n.parent;edge.active=true;
edge.line.setAttribute('data-parent',p.id);edge.line.setAttribute('data-child',n.id);
edge.line.querySelector('title').textContent=weighted?`Node ${p.id} → ${n.id}; smaller endpoint weight ${fmt(strength)}`:`Node ${p.id} → ${n.id}; uniform width`;
}
for(const edge of graphEdges.values())if(!edge.active){
edge.fromWidth=edge.width;edge.targetWidth=edge.width;edge.fromOpacity=edge.opacity;edge.targetOpacity=0;edge.reconnect=false;
}
const startSize=graphSize||{width,height},lerp=(a,b,t)=>a+(b-a)*t;
function render(t){
graphSize={width:lerp(startSize.width,width,t),height:lerp(startSize.height,height,t)};
svg.setAttribute('viewBox',`0 0 ${graphSize.width} ${graphSize.height}`);svg.style.width='100%';const fittedHeight=byId('graph-wrap').clientWidth*graphSize.height/graphSize.width;svg.style.height=(byId('graph-view').value==='screen'?Math.min(fittedHeight,Math.max(230,window.innerHeight*.7)):fittedHeight)+'px';svg.setAttribute('preserveAspectRatio','xMidYMin meet');
for(const entry of graphNodes.values()){
for(const key of ['x','y','r','opacity'])entry[key]=lerp(entry.from[key],entry.target[key],t);
entry.g.setAttribute('transform',`translate(${entry.x},${entry.y})`);entry.g.setAttribute('opacity',entry.opacity);
entry.g.querySelector('circle').setAttribute('r',entry.r);
const texts=entry.g.querySelectorAll('text');texts[0].setAttribute('y',-entry.r-7);texts[1].setAttribute('y',entry.r+17);
}
for(const edge of graphEdges.values()){
const p=graphNodes.get(edge.parent),n=graphNodes.get(edge.child);if(!p||!n){edge.line.setAttribute('visibility','hidden');continue}
edge.width=lerp(edge.fromWidth,edge.targetWidth,t);edge.opacity=lerp(edge.fromOpacity,edge.targetOpacity,t);
const attachment=edge.reconnect&&edge.fromAttach?{x:lerp(edge.fromAttach.x,p.x,t),y:lerp(edge.fromAttach.y,p.y,t),r:lerp(edge.fromAttach.r,p.r,t)}:{x:p.x,y:p.y,r:p.r};
edge.attach=attachment;
const dx=n.x-attachment.x,dy=n.y-attachment.y,len=Math.hypot(dx,dy),startGap=attachment.r+3,endGap=n.r+5+4*edge.width;
edge.line.setAttribute('stroke-width',edge.width);
edge.line.setAttribute('opacity',Math.min(edge.opacity,p.opacity,n.opacity));
//Hide an edge while emerging nodes overlap, rather than drawing it backwards.
if(len<=startGap+endGap){edge.line.setAttribute('visibility','hidden');continue}
edge.line.removeAttribute('visibility');
edge.line.setAttribute('x1',attachment.x+dx*startGap/len);edge.line.setAttribute('y1',attachment.y+dy*startGap/len);
edge.line.setAttribute('x2',n.x-dx*endGap/len);edge.line.setAttribute('y2',n.y-dy*endGap/len);
}
}
function finish(){
render(1);graphFrame=null;
for(const [id,entry] of graphNodes)if(!entry.active){entry.g.remove();graphNodes.delete(id)}
for(const [id,edge] of graphEdges){if(!edge.active){edge.line.remove();graphEdges.delete(id)}else edge.reconnect=false}
}
graphSnapshot=snapshot;
const seconds=Number(delayInput.value),duration=Math.min(500,700*(Number.isFinite(seconds)&&seconds>0?seconds:0.9));
graphTransition={render,finish,duration};
if(!animate){finish();return}
if(exportRunning){render(0);return}
const started=performance.now();render(0);
function frame(now){
const progress=Math.min(1,Math.max(0,(now-started)/duration)),eased=progress*progress*(3-2*progress);
render(eased);
if(progress<1)graphFrame=requestAnimationFrame(frame);else finish();
}
graphFrame=requestAnimationFrame(frame);
}

const diagnosticCache=new Map();
let frontierPercent=5;
for(const id of zoomPlotIds){
const svg=byId(id),controls=document.createElement('div');controls.className='plot-controls';svg.classList.add('zoom-plot');svg.setAttribute('tabindex','0');svg.setAttribute('aria-keyshortcuts','ArrowLeft ArrowRight ArrowUp ArrowDown + -');
for(const [action,label,title] of [['in','+','Zoom in'],['out','−','Zoom out'],['left','←','Move window left'],['right','→','Move window right'],['up','↑','Move window up'],['down','↓','Move window down'],['reset','Reset zoom','Reset zoom']]){
const button=document.createElement('button');button.type='button';button.textContent=label;button.title=title;button.setAttribute('aria-label',title);button.disabled=true;
if(action==='reset')button.dataset.resetPlot=id;else if(action==='in'||action==='out'){button.dataset.zoomPlot=id;button.dataset.zoomAction=action}else{button.dataset.panPlot=id;button.dataset.panDirection=action}controls.append(button);
}
const hint=document.createElement('span');hint.className='muted';hint.textContent='Drag to zoom · Arrow keys: move window · + / −: zoom';controls.append(hint);svg.before(controls);
const keyboardPlot=e=>{const action={ArrowLeft:'left',ArrowRight:'right',ArrowUp:'up',ArrowDown:'down'}[e.key],zoom=e.key==='+'||e.code==='NumpadAdd'?'in':e.key==='-'||e.key==='−'||e.code==='NumpadSubtract'?'out':null,g=spreadGeometry.get(id);
if((!action&&!zoom)||e.altKey||e.ctrlKey||e.metaKey||exportRunning||!g)return;
if(action&&!spreadViews.has(g.viewKey))return;
e.preventDefault();svg.focus({preventScroll:true});if(zoom)zoomPlot(id,zoom);else panPlot(id,action);
};svg.addEventListener('keydown',keyboardPlot);controls.addEventListener('keydown',keyboardPlot);
}
function syncPlotControls(id){
const g=spreadGeometry.get(id),view=g&&spreadViews.has(g.viewKey),epsilon=1e-9;
document.querySelector(`[data-reset-plot="${id}"]`).disabled=!view;
for(const button of document.querySelectorAll(`[data-zoom-plot="${id}"]`))button.disabled=button.dataset.zoomAction==='out'?!view:!g||(g.x[1]-g.x[0]<=(g.fullX[1]-g.fullX[0])/1000+epsilon&&g.y[1]-g.y[0]<=(g.fullY[1]-g.fullY[0])/1000+epsilon);
for(const button of document.querySelectorAll(`[data-pan-plot="${id}"]`)){
const action=button.dataset.panDirection;
button.disabled=!view||(action==='left'?g.x[0]<=g.fullX[0]+epsilon:action==='right'?g.x[1]>=g.fullX[1]-epsilon:action==='up'?g.y[1]>=g.fullY[1]-epsilon:g.y[0]<=g.fullY[0]+epsilon);
}
}
function plotLayer(svg,id,width,height,m){
const defs=el('defs'),clip=el('clipPath',{id:id+'-clip'});clip.append(el('rect',{x:m.left,y:m.top,width:width-m.left-m.right,height:height-m.top-m.bottom}));defs.append(clip);svg.append(defs);
const layer=el('g',{'clip-path':`url(#${id}-clip)`});svg.append(layer);return layer;
}
function linePlot(id,series,key,label,percent=false,empty='No data available'){
const svg=byId(id),width=Math.max(300,svg.clientWidth),height=265,m={left:66,right:24,top:22,bottom:52};
const cv=id==='topology'?'__topology__':activeCv(data[+slider.value]),viewKey=JSON.stringify([id,cv]),view=spreadViews.get(viewKey),xDomain=view?.x||[1,Math.max(2,data.length)];
const x=i=>m.left+(i+1-xDomain[0])/(xDomain[1]-xDomain[0])*(width-m.left-m.right),cacheKey=JSON.stringify([key,width,view]);
const current=svg.querySelector('[data-current-snapshot]');
if(diagnosticCache.get(id)===cacheKey&&current){current.setAttribute('x1',x(+slider.value));current.setAttribute('x2',x(+slider.value));return}
diagnosticCache.set(id,cacheKey);spreadGeometry.delete(id);syncPlotControls(id);svg.replaceChildren();svg.setAttribute('viewBox',`0 0 ${width} ${height}`);
const values=series.flatMap(s=>s.values.filter(v=>v!==null));
if(!values.length){svg.append(el('text',{x:width/2,y:80,'text-anchor':'middle'},empty));return}
const fullY=[0,percent?100:Math.max(1,...values)*1.1],yDomain=view?.y||fullY,lo=yDomain[0],hi=yDomain[1],y=v=>height-m.bottom-(v-lo)/(hi-lo)*(height-m.top-m.bottom);
spreadGeometry.set(id,{cv,viewKey,x:[...xDomain],y:[lo,hi],fullX:[1,Math.max(2,data.length)],fullY,width,height,m});syncPlotControls(id);
svg.append(el('rect',{x:m.left,y:m.top,width:width-m.left-m.right,height:height-m.top-m.bottom,fill:'none',stroke:'var(--line)'}));
for(let i=0;i<5;i++){const v=lo+(hi-lo)*i/4,yy=y(v);svg.append(el('line',{x1:m.left,x2:width-m.right,y1:yy,y2:yy,stroke:'var(--line)','stroke-opacity':.35}));svg.append(el('text',{x:m.left-9,y:yy+4,'text-anchor':'end'},percent?fmt(v)+'%':fmt(v)))}
const count=width<450?3:6,ticks=[...new Set(Array.from({length:count},(_,i)=>Math.round(xDomain[0]+(xDomain[1]-xDomain[0])*i/(count-1))-1))].filter(i=>i>=0&&i<data.length&&i+1>=xDomain[0]&&i+1<=xDomain[1]);
for(const i of ticks)svg.append(el('text',{x:x(i),y:height-m.bottom+21,'text-anchor':'middle'},data[i].cycle??i+1));
svg.append(el('text',{x:(m.left+width-m.right)/2,y:height-8,'text-anchor':'middle'},'Cycle (snapshots in log order)'));
svg.append(el('text',{transform:`translate(17,${(m.top+height-m.bottom)/2}) rotate(-90)`,'text-anchor':'middle'},label));
const layer=plotLayer(svg,id,width,height,m);
for(const s of series){let d='',previous=false;
s.values.forEach((v,i)=>{if(v===null){previous=false;return}d+=(previous?'L':'M')+x(i)+','+y(v);previous=true});
layer.append(el('path',{d,fill:'none',stroke:s.color,'stroke-width':2}));
s.values.forEach((v,i)=>{if(v===null)return;const g=el('g',{role:'button',tabindex:0,'aria-label':`${s.name}, cycle ${data[i].cycle}, ${fmt(v)}${percent?' percent':''}`});
g.append(el('circle',{cx:x(i),cy:y(v),r:9,fill:'transparent'}),el('circle',{cx:x(i),cy:y(v),r:3,fill:s.color}),el('title',{},`${s.name}; snapshot ${i+1}; cycle ${data[i].cycle}; ${fmt(v)}${percent?'%':''}`));
const choose=()=>{if(exportRunning||performance.now()<spreadClickBlockedUntil)return;stop();slider.value=i;draw()};g.addEventListener('click',choose);g.addEventListener('keydown',e=>{if(e.key==='Enter'||e.key===' '){e.preventDefault();choose()}});layer.append(g);
});}
layer.append(el('line',{'data-current-snapshot':'',x1:x(+slider.value),x2:x(+slider.value),y1:m.top,y2:height-m.bottom,stroke:'var(--fg)','stroke-opacity':.4}));
}
function drawDiagnostics(cv){
const points=histories[cv]||[],ids=[...trackedNodes].filter(id=>id!==0).sort((a,b)=>a-b),colors=ids.map(trackColor);
linePlot('tracked-probability',ids.map(id=>({name:`Node ${id}`,color:trackColor(id),values:points.map(p=>p.probabilities?.[id]===undefined?null:100*p.probabilities[id])})),[cv,ids,colors],'Draw probability',true,ids.length?'No reported draw probabilities for tracked nodes':'Track nodes to show their draw probabilities');
const fraction=frontierPercent/100,range=scoreRanges[cv],gap=range?fraction*(range.max-range.min):0;
const frontier=points.map(p=>{if(!p.probabilities)return null;let probability=0;for(const n of p.s.nodes){if(!n.id)continue;const row=n.cvs[cv];if(!row)continue;const near=direction.value==='min'?row.score<=p.min+gap:row.score>=p.max-gap;if(near)probability+=p.probabilities[n.id]??0}return 100*probability});
linePlot('frontier-probability',[{name:'Frontier',color:'var(--series)',values:frontier}],[cv,direction.value,fraction],'Draw probability',true);
const value=frontier[+slider.value];byId('frontier-summary').textContent=`${cv||'Selected CV'} · Frontier width ${fmt(gap)} CV units · Current probability ${value===null||value===undefined?'unavailable':fmt(value)+'%'}`;
linePlot('topology',[
{name:'Pool nodes',color:'var(--series)',values:topologyHistory.map(p=>p.nodes)},
{name:'Tips',color:'var(--other)',values:topologyHistory.map(p=>p.tips)},
{name:'Branch points',color:'#8b5cf6',values:topologyHistory.map(p=>p.branches)}
],'topology','Node count');
const t=topologyHistory[+slider.value];byId('topology-summary').textContent=`${t.nodes} real nodes · ${t.tips} tips · ${t.branches} branch points (virtual root excluded)`;
}
byId('frontier-width').addEventListener('change',()=>{const input=byId('frontier-width');if(input.value.trim()&&input.checkValidity()){frontierPercent=Number(input.value);draw()}else input.reportValidity()});

function drawSpread(cv,replicas=false){
const plotId=replicas?'replica-spread':'spread',svg=byId(plotId),points=(replicas?replicaHistories:histories)[cv]||[],p=points[+slider.value];
byId(plotId+'-summary').textContent=p&&p.count?`${cv} · Minimum ${fmt(p.min)} · Maximum ${fmt(p.max)} · Range ${fmt(p.max-p.min)} · Median ${fmt(p.median)} · ${p.count} scored ${replicas?'replicas · Run '+(p.run??'unknown'):'nodes'}`:`${cv||'Selected CV'} · ${replicas?'No replica offers':'No retained scores'} reported in this snapshot`;
const width=Math.max(300,svg.clientWidth),height=265,m={left:66,right:24,top:22,bottom:52};
const viewKey=JSON.stringify([plotId,cv]),view=spreadViews.get(viewKey);
const xDomain=view?.x||[1,Math.max(2,data.length)],x=v=>m.left+(v-xDomain[0])/(xDomain[1]-xDomain[0])*(width-m.left-m.right);
byId(plotId).classList.add('zoom-plot');
document.querySelector(`[data-reset-plot="${plotId}"]`).disabled=!view;
const key=JSON.stringify([cv,width,replicas?direction.value:null,view]),current=svg.querySelector('[data-current-snapshot]');
if(spreadCaches[plotId]===key&&current){current.setAttribute('x1',x(+slider.value+1));current.setAttribute('x2',x(+slider.value+1));return}
spreadCaches[plotId]=key;spreadGeometry.delete(plotId);syncPlotControls(plotId);svg.replaceChildren();svg.setAttribute('viewBox',`0 0 ${width} ${height}`);
let lo=Infinity,hi=-Infinity;for(const p of points)if(p.count){lo=Math.min(lo,p.min);hi=Math.max(hi,p.max)}
if(!Number.isFinite(lo)){svg.append(el('text',{x:width/2,y:80,'text-anchor':'middle'},replicas?'No replica offers available for this CV':'No retained scores available for this CV'));return}
if(replicas)for(const p of points){const v=direction.value==='min'?p.poolMin:p.poolMax;if(v!==null){lo=Math.min(lo,v);hi=Math.max(hi,v)}}
const pad=(hi-lo)*.1||Math.max(1,Math.abs(lo)*.05);lo-=pad;hi+=pad;
const fullY=[lo,hi];if(view?.y){[lo,hi]=view.y}
spreadGeometry.set(plotId,{cv,viewKey,x:[...xDomain],y:[lo,hi],fullX:[1,Math.max(2,data.length)],fullY,width,height,m});syncPlotControls(plotId);
const y=v=>height-m.bottom-(v-lo)/(hi-lo)*(height-m.top-m.bottom);
svg.append(el('rect',{x:m.left,y:m.top,width:width-m.left-m.right,height:height-m.top-m.bottom,fill:'none',stroke:'var(--line)'}));
for(let i=0;i<5;i++){const v=lo+(hi-lo)*i/4,yy=y(v);svg.append(el('line',{x1:m.left,x2:width-m.right,y1:yy,y2:yy,stroke:'var(--line)','stroke-opacity':.35}));svg.append(el('text',{x:m.left-9,y:yy+4,'text-anchor':'end'},fmt(v)))}
const count=width<450?3:6,ticks=[...new Set(Array.from({length:count},(_,i)=>Math.round(xDomain[0]+(xDomain[1]-xDomain[0])*i/(count-1))))].filter(v=>v>=1&&v<=data.length&&v>=xDomain[0]&&v<=xDomain[1]);
for(const tick of ticks)svg.append(el('text',{x:x(tick),y:height-m.bottom+21,'text-anchor':'middle'},data[tick-1].cycle??tick));
svg.append(el('text',{x:(m.left+width-m.right)/2,y:height-8,'text-anchor':'middle'},'Cycle (snapshots in log order)'));
svg.append(el('text',{transform:`translate(17,${(m.top+height-m.bottom)/2}) rotate(-90)`,'text-anchor':'middle'},'CV score'));
const defs=el('defs'),clip=el('clipPath',{id:plotId+'-clip'});clip.append(el('rect',{x:m.left,y:m.top,width:width-m.left-m.right,height:height-m.top-m.bottom}));defs.append(clip);svg.append(defs);
const layer=el('g',{'clip-path':`url(#${plotId}-clip)`});svg.append(layer);
//Separate segments prevent interpolation across snapshots with missing CV data.
let segment=[];
function flush(){if(!segment.length)return;
const upper=segment.map(p=>`${x(p.x)},${y(p.max)}`),lower=segment.slice().reverse().map(p=>`${x(p.x)},${y(p.min)}`);
layer.append(el('path',{d:'M'+upper.join('L')+'L'+lower.join('L')+'Z',fill:'var(--other)','fill-opacity':.18,stroke:'none'}));
for(const field of ['min','max','median'])layer.append(el('path',{d:'M'+segment.map(p=>`${x(p.x)},${y(p[field])}`).join('L'),fill:'none',stroke:field==='median'?'var(--series)':'var(--other)','stroke-width':field==='median'?2:1}));segment=[];
}
for(const p of points){if(p.count)segment.push(p);else flush()}flush();
if(replicas){let d='',previous=false;for(const p of points){const v=direction.value==='min'?p.poolMin:p.poolMax;if(v===null){previous=false;continue}d+=(previous?'L':'M')+x(p.x)+','+y(v);previous=true}layer.append(el('path',{d,fill:'none',stroke:'#8b5cf6','stroke-width':2,'stroke-dasharray':'5 4'}))}
for(const p of points)if(p.count){
const g=el('g',{role:'button',tabindex:0,'aria-label':`Cycle ${p.s.cycle}, minimum ${p.min}, maximum ${p.max}, median ${p.median}`});
g.append(el('line',{x1:x(p.x),x2:x(p.x),y1:y(p.min),y2:y(p.max),stroke:'var(--other)','stroke-width':2}));
g.append(el('circle',{cx:x(p.x),cy:y(p.median),r:10,fill:'transparent'}));g.append(el('circle',{cx:x(p.x),cy:y(p.median),r:3,fill:'var(--series)'}));
if(replicas)for(const r of p.replicas){const dot=el('circle',{cx:x(p.x),cy:y(r.score),r:3,fill:'var(--other)','fill-opacity':.65,'data-replica-id':r.id});dot.append(el('title',{},`Replica ${r.id}; run ${p.run}; ${cv}: ${fmt(r.score)}`));g.append(dot)}
g.append(el('title',{},`Snapshot ${p.x}; cycle ${p.s.cycle}; ${cv}: minimum ${fmt(p.min)}, maximum ${fmt(p.max)}, median ${fmt(p.median)}; ${p.count} ${replicas?'replicas':'nodes'}`));
const choose=()=>{if(exportRunning||performance.now()<spreadClickBlockedUntil)return;stop();slider.value=p.x-1;draw()};g.addEventListener('click',choose);g.addEventListener('keydown',e=>{if(e.key==='Enter'||e.key===' '){e.preventDefault();choose()}});layer.append(g);
}
layer.append(el('line',{'data-current-snapshot':'',x1:x(+slider.value+1),x2:x(+slider.value+1),y1:m.top,y2:height-m.bottom,stroke:'var(--fg)','stroke-opacity':.4}));
}

//Zoom in data coordinates so labels, line widths and replica dots remain legible.
function limitSpreadDomain(domain,full){
const range=full[1]-full[0],span=Math.min(range,Math.max(range/1000,domain[1]-domain[0]));
const low=Math.max(full[0],Math.min(full[1]-span,domain[0]));return [low,low+span];
}
function redrawPlot(id,cv){
if(id==='spread'||id==='replica-spread')drawSpread(cv,id==='replica-spread');else if(id==='trend')drawTrend(cv);else drawDiagnostics(activeCv(data[+slider.value]));
}
function updateSpreadView(id,x,y){
const g=spreadGeometry.get(id);if(!g)return;
x=limitSpreadDomain(x,g.fullX);y=limitSpreadDomain(y,g.fullY);
const full=(d,f)=>Math.abs(d[0]-f[0])<1e-9&&Math.abs(d[1]-f[1])<1e-9;
if(full(x,g.fullX)&&full(y,g.fullY))spreadViews.delete(g.viewKey);else spreadViews.set(g.viewKey,{x,y});
redrawPlot(id,g.cv);
}
for(const id of zoomPlotIds){
const svg=byId(id);
const coordinates=(e,g=spreadGeometry.get(id))=>{const r=svg.getBoundingClientRect();if(!g)return null;
return {g,fx:(e.clientX-r.left)/r.width*g.width,fy:(e.clientY-r.top)/r.height*g.height}};
const clamp=(p,g)=>({x:Math.max(g.m.left,Math.min(g.width-g.m.right,p.fx)),y:Math.max(g.m.top,Math.min(g.height-g.m.bottom,p.fy))});
svg.addEventListener('pointerdown',e=>{if(e.button!==0||exportRunning)return;const p=coordinates(e);if(!p)return;const {g,fx,fy}=p;
if(fx<g.m.left||fx>g.width-g.m.right||fy<g.m.top||fy>g.height-g.m.bottom)return;
svg.focus({preventScroll:true});stop();spreadDrag={id,pointer:e.pointerId,startX:e.clientX,startY:e.clientY,start:clamp(p,g),g:{...g,x:[...g.x],y:[...g.y]},moved:false,box:null};
});
svg.addEventListener('pointermove',e=>{const d=spreadDrag;if(!d||d.id!==id||d.pointer!==e.pointerId)return;
if(!d.moved&&Math.hypot(e.clientX-d.startX,e.clientY-d.startY)<4)return;
if(!d.moved){d.moved=true;svg.setPointerCapture(e.pointerId);d.box=el('rect',{'data-zoom-selection':'',fill:'var(--series)','fill-opacity':.15,stroke:'var(--series)','stroke-width':1.5,'stroke-dasharray':'4 3','pointer-events':'none'});svg.append(d.box)}
e.preventDefault();spreadClickBlockedUntil=performance.now()+350;
const p=clamp(coordinates(e,d.g),d.g);
for(const [k,v] of Object.entries({x:Math.min(d.start.x,p.x),y:Math.min(d.start.y,p.y),width:Math.abs(p.x-d.start.x),height:Math.abs(p.y-d.start.y)}))d.box.setAttribute(k,v);
});
const finish=e=>{const d=spreadDrag;if(!d||d.id!==id||d.pointer!==e.pointerId)return;
spreadDrag=null;d.box?.remove();if(d.moved)spreadClickBlockedUntil=performance.now()+350;
if(svg.hasPointerCapture(e.pointerId))svg.releasePointerCapture(e.pointerId);
if(!d.moved||e.type==='pointercancel')return;
const g=d.g,p=clamp(coordinates(e,g),g),r=svg.getBoundingClientRect();
//Tiny selections are treated as canceled drags, never as point clicks.
if(Math.abs(p.x-d.start.x)*r.width/g.width<6||Math.abs(p.y-d.start.y)*r.height/g.height<6)return;
const xv=v=>g.x[0]+(v-g.m.left)/(g.width-g.m.left-g.m.right)*(g.x[1]-g.x[0]);
const yv=v=>g.y[1]-(v-g.m.top)/(g.height-g.m.top-g.m.bottom)*(g.y[1]-g.y[0]);
updateSpreadView(id,[xv(Math.min(d.start.x,p.x)),xv(Math.max(d.start.x,p.x))],[yv(Math.max(d.start.y,p.y)),yv(Math.min(d.start.y,p.y))]);
};
svg.addEventListener('pointerup',finish);svg.addEventListener('pointercancel',finish);
svg.addEventListener('dblclick',()=>{if(exportRunning)return;const g=spreadGeometry.get(id);if(g){spreadViews.delete(g.viewKey);redrawPlot(id,g.cv)}});
}
function panPlot(id,action){
if(exportRunning)return;const g=spreadGeometry.get(id);if(!g||!spreadViews.has(g.viewKey))return;stop();
const dx=(g.x[1]-g.x[0])*.2*(action==='left'?-1:action==='right'?1:0),dy=(g.y[1]-g.y[0])*.2*(action==='down'?-1:action==='up'?1:0);
updateSpreadView(id,g.x.map(v=>v+dx),g.y.map(v=>v+dy));
}
function zoomPlot(id,action){
if(exportRunning)return;const g=spreadGeometry.get(id);if(!g)return;stop();
const factor=action==='in'?.8:1.25,scale=d=>{const middle=(d[0]+d[1])/2,half=(d[1]-d[0])*factor/2;return [middle-half,middle+half]};
updateSpreadView(id,scale(g.x),scale(g.y));
}
document.querySelectorAll('[data-zoom-plot]').forEach(button=>button.addEventListener('click',()=>zoomPlot(button.dataset.zoomPlot,button.dataset.zoomAction)));
document.querySelectorAll('[data-pan-plot]').forEach(button=>button.addEventListener('click',()=>panPlot(button.dataset.panPlot,button.dataset.panDirection)));
document.querySelectorAll('[data-reset-plot]').forEach(button=>button.addEventListener('click',()=>{if(exportRunning)return;const id=button.dataset.resetPlot,g=spreadGeometry.get(id);if(g){spreadViews.delete(g.viewKey);redrawPlot(id,g.cv)}}));

function drawTrend(cv){
const svg=byId('trend'),width=Math.max(300,svg.clientWidth),height=265,m={left:66,right:24,top:22,bottom:52};
const viewKey=JSON.stringify(['trend',cv]),view=spreadViews.get(viewKey),xDomain=view?.x||[1,Math.max(2,data.length)];
const x=v=>m.left+(v-xDomain[0])/(xDomain[1]-xDomain[0])*(width-m.left-m.right);
const key=JSON.stringify([cv,direction.value,width,view]),current=svg.querySelector('[data-current-snapshot]');
if(trendCacheKey===key&&current){current.setAttribute('x1',x(+slider.value+1));current.setAttribute('x2',x(+slider.value+1));return}
trendCacheKey=key;spreadGeometry.delete('trend');syncPlotControls('trend');svg.replaceChildren();svg.setAttribute('viewBox',`0 0 ${width} ${height}`);
const points=(histories[cv]||[]).map(p=>({...p,pool:direction.value==='min'?p.min:p.max}));
let lo=Infinity,hi=-Infinity;for(const p of points)for(const v of [p.pool,p.record])if(v!==null){lo=Math.min(lo,v);hi=Math.max(hi,v)}
if(!Number.isFinite(lo)){svg.append(el('text',{x:width/2,y:80,'text-anchor':'middle'},'No scores available for this CV'));return}
const pad=(hi-lo)*.1||Math.max(1,Math.abs(lo)*.05);lo-=pad;hi+=pad;
const fullY=[lo,hi];if(view?.y){[lo,hi]=view.y}
spreadGeometry.set('trend',{cv,viewKey,x:[...xDomain],y:[lo,hi],fullX:[1,Math.max(2,data.length)],fullY,width,height,m});syncPlotControls('trend');
const y=v=>height-m.bottom-(v-lo)/(hi-lo)*(height-m.top-m.bottom);
svg.append(el('rect',{x:m.left,y:m.top,width:width-m.left-m.right,height:height-m.top-m.bottom,fill:'none',stroke:'var(--line)'}));
for(let i=0;i<5;i++){const v=lo+(hi-lo)*i/4,yy=y(v);svg.append(el('line',{x1:m.left,x2:width-m.right,y1:yy,y2:yy,stroke:'var(--line)','stroke-opacity':.35}));svg.append(el('text',{x:m.left-9,y:yy+4,'text-anchor':'end'},fmt(v)));}
const count=width<450?3:6;const ticks=[...new Set(Array.from({length:count},(_,i)=>Math.round(xDomain[0]+(xDomain[1]-xDomain[0])*i/(count-1))))].filter(v=>v>=1&&v<=data.length&&v>=xDomain[0]&&v<=xDomain[1]);for(const tick of ticks)svg.append(el('text',{x:x(tick),y:height-m.bottom+21,'text-anchor':'middle'},tick));
svg.append(el('text',{x:(m.left+width-m.right)/2,y:height-8,'text-anchor':'middle'},'Snapshot in log order'));
svg.append(el('text',{transform:`translate(17,${(m.top+height-m.bottom)/2}) rotate(-90)`,'text-anchor':'middle'},`CV score`));
const layer=plotLayer(svg,'trend',width,height,m);
for(const [key,color] of [['record','var(--series)'],['pool','var(--other)']]){let d='',prev=null;for(const p of points){if(p[key]===null){prev=null;continue}d+=prev?`H${x(p.x)}V${y(p[key])}`:`M${x(p.x)},${y(p[key])}`;prev=p;}layer.append(el('path',{d,fill:'none',stroke:color,'stroke-width':2}));
for(const p of points)if(p[key]!==null){const g=el('g',{role:'button',tabindex:0,'aria-label':`Snapshot ${p.x}, ${key} score ${p[key]}`});g.append(el('circle',{cx:x(p.x),cy:y(p[key]),r:10,fill:'transparent'}));g.append(el('circle',{cx:x(p.x),cy:y(p[key]),r:3,fill:color}));g.append(el('title',{},`Snapshot ${p.x}; cycle ${p.s.cycle}; ${key==='record'?'recorded best':'pool best'} ${p[key]}`));const choose=()=>{if(exportRunning||performance.now()<spreadClickBlockedUntil)return;stop();slider.value=p.x-1;draw()};g.addEventListener('click',choose);g.addEventListener('keydown',e=>{if(e.key==='Enter'){choose()}});layer.append(g);}}
layer.append(el('line',{'data-current-snapshot':'',x1:x(+slider.value+1),x2:x(+slider.value+1),y1:m.top,y2:height-m.bottom,stroke:'var(--fg)','stroke-opacity':.4}));
}
function draw(){if(exportRunning&&!exportStepping)return;const i=+slider.value,s=data[i],cv=activeCv(s);if(cv!==lastCv){direction.value=cvDirections[cv]||'min';lastCv=cv}
byId('color-legend').textContent=colorSelect.value==='id'?'Color identifies node ID; neighboring IDs use contrasting hues':'Darker fill = better score on a fixed scale per CV';
byId('edge-legend').textContent=edgeSelect.value==='weight'?'Thicker edges = both endpoints more strongly weighted (compressed scale); root edges uniform':'Edges have uniform width';
byId('position').textContent=`${i+1} / ${data.length}`;byId('prev').disabled=i===0;byId('next').disabled=i===data.length-1;
byId('status').textContent=`${s.source?s.source+' · ':''}Cycle ${s.cycle} · TIGER2 run ${s.run??'unknown'} · CV ${cv||'unknown'} · Active CV ${s.active||'unknown'} · ${s.nodes.filter(n=>n.id).length} real nodes`;
let total=0,complete=true;for(const n of s.nodes){if(!n.id)continue;const row=n.cvs[cv];if(!row||row.weight===null){complete=false;continue}total+=row.weight}
if(!complete)total=null;
const tbody=byId('rows');tbody.replaceChildren();for(const n of s.nodes){const r=n.cvs[cv],tr=document.createElement('tr');tr.className=trackedNodes.has(n.id)?'selected':'';tr.setAttribute('data-node-id',n.id);tr.setAttribute('tabindex','0');tr.setAttribute('aria-label',`Node ${n.id}, ${trackedNodes.has(n.id)?'tracked':'not tracked'}`);tr.addEventListener('keydown',e=>{if(e.key==='Enter'||e.key===' '){e.preventDefault();selectNode(n.id)}});const cells=[n.id,n.parent,n.id?(r?fmt(r.score):'—'):'virtual',r?fmt(r.attempts):'—',r?fmt(r.successes):'—',r?fmt(r.weight):'—',r&&total>0?(100*r.weight/total).toFixed(2)+'%':'—'];for(const value of cells){const td=document.createElement('td');td.textContent=value;tr.append(td)}tr.addEventListener('click',()=>selectNode(n.id));tbody.append(tr)}
const node=s.nodes.find(n=>n.id===selected),r=node?.cvs[cv];byId('selected').textContent=node?`Node ${node.id} · Parent ${node.parent} · Attempts ${r?fmt(r.attempts):'—'} · Successes ${r?fmt(r.successes):'—'} · Weight ${r?fmt(r.weight):'—'}`:selected!==null?`Node ${selected} is absent from this snapshot.`:'Select a node to inspect its counters.';
const trackedList=byId('tracked-list');
byId('tracked-panel').hidden=!trackedNodes.size;
const legendKey=JSON.stringify([...trackedNodes].sort((a,b)=>a-b).map(id=>[id,trackColor(id),s.nodes.some(n=>n.id===id)]));
if(trackedList.dataset.key!==legendKey){trackedList.dataset.key=legendKey;trackedList.replaceChildren();
for(const id of [...trackedNodes].sort((a,b)=>a-b)){
const present=s.nodes.some(n=>n.id===id),item=document.createElement('div'),input=document.createElement('input'),button=document.createElement('button');item.className='tracked-item';input.type='color';input.value=trackColor(id);input.setAttribute('aria-label',`Tracking color for node ${id}`);input.dataset.nodeId=id;
input.addEventListener('change',()=>{if(exportRunning)return;trackedColors.set(id,input.value);draw()});
button.type='button';button.textContent=`Node ${id}${present?'':' (absent)'} ×`;button.className=present?'':'absent';button.style.borderColor=trackColor(id);button.setAttribute('aria-label',`Stop tracking node ${id}${present?'':', absent from this snapshot'}`);button.addEventListener('click',()=>selectNode(id));item.append(input,button);trackedList.append(item);
}}
for(const control of trackedList.querySelectorAll('input,button'))control.disabled=exportRunning;
byId('clear-tracked').disabled=exportRunning||!trackedNodes.size;
drawGraph(s,cv,total);drawSpread(cv);drawSpread(cv,true);drawTrend(cv);drawDiagnostics(cv);
}
byId('prev').addEventListener('click',()=>{stop();slider.value=Math.max(0,+slider.value-1);draw()});byId('next').addEventListener('click',()=>{stop();slider.value=Math.min(data.length-1,+slider.value+1);draw()});slider.addEventListener('input',()=>{stop();draw()});
byId('play').addEventListener('click',()=>{if(timer){stop();return}if(!delayInput.checkValidity()){delayInput.reportValidity();return}if(+slider.value===data.length-1){slider.value=0;draw()}byId('play').textContent='Pause';schedule()});
function adjustDelay(change){
const value=Number(delayInput.value),base=delayInput.value.trim()&&Number.isFinite(value)?value:report.delay??0.9;
delayInput.value=String(Math.round(Math.max(0.05,Math.min(60,base+change))*1000000)/1000000);
delayInput.dispatchEvent(new Event('input'));
}
byId('delay-down').addEventListener('click',()=>adjustDelay(-0.05));
byId('delay-up').addEventListener('click',()=>adjustDelay(0.05));
delayInput.addEventListener('input',()=>{if(delayInput.checkValidity()&&timer)schedule()});
delayInput.addEventListener('change',()=>{if(!delayInput.checkValidity()){stop();delayInput.reportValidity()}});
cvSelect.addEventListener('change',draw);direction.addEventListener('change',()=>{const cv=activeCv(data[+slider.value]);if(cv)cvDirections[cv]=direction.value;draw()});sizeSelect.addEventListener('change',draw);
motionInput.addEventListener('change',draw);
if(reducedMotion&&reducedMotion.addEventListener)reducedMotion.addEventListener('change',draw);
edgeSelect.addEventListener('change',draw);colorSelect.addEventListener('change',draw);
if(colorTheme&&colorTheme.addEventListener)colorTheme.addEventListener('change',draw);
//GIF89a encoder: fixed 256-color RGB palette and GIF LZW compression.
//All frames cover the full canvas; no external libraries or services are used.
function gifLzw(pixels){
const bytes=[];let buffer=0,bits=0,size=9,next=258,dict=new Map();
function emit(code){
buffer|=code<<bits;bits+=size;
while(bits>=8){bytes.push(buffer&255);buffer>>>=8;bits-=8}
//The decoder builds its dictionary one code later than the encoder.
if(next>=(1<<size)&&size<12)size++;
}
emit(256);let prefix=pixels[0];
for(let i=1;i<pixels.length;i++){
const value=pixels[i],key=prefix*256+value,found=dict.get(key);
if(found!==undefined){prefix=found;continue}
emit(prefix);
if(next<4096)dict.set(key,next++);
else{emit(256);dict=new Map();next=258;size=9}
prefix=value;
}
emit(prefix);emit(257);if(bits)bytes.push(buffer&255);
return Uint8Array.from(bytes);
}
function gifEncoder(width,height){
const chunks=[],word=n=>[n&255,(n>>8)&255],ascii=s=>Array.from(s,c=>c.charCodeAt(0));
const palette=[];
for(let i=0;i<256;i++)palette.push(Math.round((i>>5)*255/7),Math.round(((i>>2)&7)*255/7),Math.round((i&3)*255/3));
chunks.push(Uint8Array.from([...ascii('GIF89a'),...word(width),...word(height),247,0,0,...palette,33,255,11,...ascii('NETSCAPE2.0'),3,1,0,0,0]));
return {
add(rgba,delay){
const pixels=new Uint8Array(width*height);
for(let i=0,j=0;i<pixels.length;i++,j+=4)pixels[i]=(rgba[j]&224)|((rgba[j+1]>>3)&28)|(rgba[j+2]>>6);
const compressed=gifLzw(pixels),header=[33,249,4,4,...word(delay),0,0,44,0,0,0,0,...word(width),...word(height),0,8];
chunks.push(Uint8Array.from(header));
for(let i=0;i<compressed.length;i+=255){const block=compressed.subarray(i,i+255);chunks.push(Uint8Array.of(block.length),block)}
chunks.push(Uint8Array.of(0));
},
finish(){return new Blob([...chunks,Uint8Array.of(59)],{type:'image/gif'})}
};
}
function gifBounds(){
let width=Math.max(520,byId('graph-wrap').clientWidth),height=230;
for(const s of data){
const nodes=new Map(s.nodes.map(n=>[n.id,n])),children=new Map(s.nodes.map(n=>[n.id,0]));let depth=0;
for(const n of s.nodes){if(nodes.has(n.parent)&&n.parent!==n.id)children.set(n.parent,children.get(n.parent)+1);
let cursor=n,steps=0,seen=new Set([n.id]);
while(nodes.has(cursor.parent)&&!seen.has(cursor.parent)){cursor=nodes.get(cursor.parent);seen.add(cursor.id);steps++}
depth=Math.max(depth,steps);
}
const leaves=Math.max(1,[...children.values()].filter(n=>n===0).length);
width=Math.max(width,Math.max(1,depth)*160+140);height=Math.max(height,leaves*88+40);
}
const scale=Math.min(1,1200/Math.max(width,height));
return {width,height,pixelsWide:Math.max(1,Math.round(width*scale)),pixelsHigh:Math.max(1,Math.round(height*scale))};
}
async function gifPixels(canvas,context,bounds){
const original=byId('graph'),clone=original.cloneNode(true);
const originals=[original,...original.querySelectorAll('*')],copies=[clone,...clone.querySelectorAll('*')];
const properties=['fill','stroke','stroke-width','fill-opacity','stroke-opacity','opacity','font-family','font-size','font-weight','text-anchor','visibility','vector-effect'];
for(let i=0;i<originals.length;i++){
const computed=getComputedStyle(originals[i]);
for(const property of properties)copies[i].style.setProperty(property,computed.getPropertyValue(property));
}
clone.setAttribute('xmlns',NS);clone.setAttribute('viewBox',`0 0 ${bounds.width} ${bounds.height}`);
clone.setAttribute('width',canvas.width);clone.setAttribute('height',bounds.pixelsHigh);
clone.style.width=canvas.width+'px';clone.style.height=bounds.pixelsHigh+'px';
const blob=new Blob([new XMLSerializer().serializeToString(clone)],{type:'image/svg+xml;charset=utf-8'}),url=URL.createObjectURL(blob);
try{
const image=await new Promise((resolve,reject)=>{const img=new Image();img.onload=()=>resolve(img);img.onerror=()=>reject(new Error('Could not rasterize the graph.'));img.src=url});
context.fillStyle=getComputedStyle(document.body).backgroundColor;context.fillRect(0,0,canvas.width,canvas.height);context.drawImage(image,0,bounds.legendHeight,canvas.width,bounds.pixelsHigh);
if(trackedNodes.size){
const ids=[...trackedNodes].sort((a,b)=>a-b),font=12,pad=8,line=19,legendWidth=Math.min(190,canvas.width-16),legendHeight=pad*2+line*(ids.length+1),left=canvas.width-legendWidth-8,top=8;
context.save();context.globalAlpha=.94;context.fillStyle=getComputedStyle(byId('graph-wrap')).backgroundColor;context.fillRect(left,top,legendWidth,legendHeight);context.globalAlpha=1;context.strokeStyle='#94a3b8';context.lineWidth=1;context.strokeRect(left,top,legendWidth,legendHeight);context.fillStyle=getComputedStyle(document.body).color;context.font=`bold ${font}px sans-serif`;context.fillText('Tracked nodes',left+pad,top+pad+font);if(left>180)context.fillText(`Cycle ${data[+slider.value].cycle} · CV ${activeCv(data[+slider.value])||'unknown'}`,8,top+pad+font);context.font=`${font}px sans-serif`;
const snapshot=data[+slider.value];ids.forEach((id,i)=>{const y=top+pad+line*(i+1)+font;context.fillStyle=trackColor(id);context.fillRect(left+pad,y-9,10,10);context.fillStyle=getComputedStyle(document.body).color;context.fillText(`Node ${id}${snapshot.nodes.some(n=>n.id===id)?'':' (absent)'}`,left+pad+17,y)});context.restore();
}
return context.getImageData(0,0,canvas.width,canvas.height).data;
}finally{URL.revokeObjectURL(url)}
}
async function renderGif(){
if(exportRunning){exportAbort=true;byId('export-status').textContent='Cancelling GIF rendering…';return}
if(!delayInput.checkValidity()){delayInput.reportValidity();return}
stop();if(graphFrame!==null){cancelAnimationFrame(graphFrame);graphFrame=null}if(graphTransition)graphTransition.finish();
const savedSnapshot=slider.value,savedMotion=motionInput.checked,controls=[...document.querySelectorAll('main button, main input, main select, #tracked-panel button, #tracked-panel input')].filter(n=>n.id!=='render-gif'),disabled=controls.map(n=>n.disabled);
const button=byId('render-gif'),status=byId('export-status');
exportRunning=true;exportAbort=false;controls.forEach(n=>n.disabled=true);button.textContent='Cancel render';
const delay=Math.max(1,Math.round(Number(delayInput.value)*100));let frames=0;
function step(index){slider.value=index;exportStepping=true;try{draw()}finally{exportStepping=false}}
try{
const bounds=gifBounds(),canvas=document.createElement('canvas');bounds.legendHeight=trackedNodes.size?32+19*(trackedNodes.size+1):0;canvas.width=bounds.pixelsWide;canvas.height=bounds.pixelsHigh+bounds.legendHeight;
const context=canvas.getContext('2d',{willReadFrequently:true}),encoder=gifEncoder(canvas.width,canvas.height);
async function capture(hold){if(exportAbort)return;const pixels=await gifPixels(canvas,context,bounds);if(exportAbort)return;encoder.add(pixels,hold);frames++}
if(!context)throw new Error('Canvas rendering is unavailable in this browser.');
if(document.fonts&&document.fonts.ready)await document.fonts.ready;
for(let i=0;i<data.length&&!exportAbort;i++){
status.textContent=`Rendering GIF: snapshot ${i+1} / ${data.length} · ${frames} frames`;
//The first snapshot appears immediately, with no transition from the viewer.
motionInput.checked=i===0?false:savedMotion;step(i);
const transition=graphTransition;
const smooth=i>0&&savedMotion&&!(reducedMotion&&reducedMotion.matches);
if(smooth){
const motionDelay=Math.max(1,Math.min(delay,Math.round(transition.duration/10))),count=Math.min(motionDelay,Math.max(1,Math.ceil(transition.duration/80)));
for(let j=0;j<count&&!exportAbort;j++){
const t=(j+1)/count;transition.render(t*t*(3-2*t));
const hold=Math.floor((j+1)*motionDelay/count)-Math.floor(j*motionDelay/count)+(j===count-1?delay-motionDelay:0);
await capture(hold);
}
transition.finish();
}else{transition.finish();await capture(delay)}
await new Promise(resolve=>setTimeout(resolve,0));
}
if(exportAbort){status.textContent='GIF rendering cancelled.';return}
const gif=encoder.finish(),url=URL.createObjectURL(gif),link=document.createElement('a');link.href=url;link.download='swarm-pool.gif';document.body.append(link);link.click();link.remove();setTimeout(()=>URL.revokeObjectURL(url),60000);
status.textContent=`GIF ready · ${frames} frames · ${canvas.width} × ${canvas.height} pixels`;
}catch(error){status.textContent='GIF rendering failed: '+error.message}
finally{
if(graphFrame!==null){cancelAnimationFrame(graphFrame);graphFrame=null}
exportRunning=false;exportStepping=false;motionInput.checked=savedMotion;slider.value=savedSnapshot;
graphSnapshot=null;controls.forEach((n,i)=>n.disabled=disabled[i]);button.textContent='Render GIF';draw();
}
}
byId('render-gif').addEventListener('click',renderGif);

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
