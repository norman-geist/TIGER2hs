#!/usr/bin/env python3
"""Coordinate-PCA/MSM prototype for SHOOTING-TIGER2h baseline segments.

Install:
    python -m pip install numpy scipy scikit-learn matplotlib MDAnalysis deeptime networkx
Run:
    python shooting_msm.py JOB.psf --dt-ps 20 --lag-ps 100 \
        --lags-ps 20 40 100 200 500 --output output --out msm_results
Input: output/REPLICA/JOB.MILESTONE.SEGMENT.dcd (numeric IDs).
Only closed, immutable, regularly spaced baseline segments belong here.
Temporary *.tmp.dcd files are excluded. Never concatenate segments for counting.
--dt-ps is the actual saved-frame spacing, NOT necessarily a TIGER2h cycle.

Coordinates must be whole and have consistent atom order. Optional --unwrap
makes bonded fragments whole; it does not assemble separate chains into a complex.
Alignment removes global translation/rotation but not domain movements. A stable
reference-domain --align-selection can be useful for multidomain proteins.

Default estimator: nonreversible row-normalized counts without pseudocounts.
--reversible uses deeptime's reversible maximum-likelihood estimator. Neither
estimator automatically corrects non-equilibrium starting-state or stopping bias.
The model is restricted to the most populated strongly connected component.
Excluded states cannot be assigned relative equilibrium populations.

FES: -RT log of a histogram weighted by MSM stationary populations, using the
observed within-state coordinate distribution. Conditional on the active component
and local equilibration; NOT a validated equilibrium landscape by default.
The raw PCA histogram is labelled sampling density and is not a free energy.

Validation: implied timescales, available-pair counts, segment lengths, and a
segment-held-out CK diagnostic. Exchange-coupled segments need not be independent;
this diagnostic supplies neither independent-run validation nor confidence bounds.
Frames are held in RAM (float32). --stride/--selection reduce memory; --fit-max
limits PCA/clustering training frames only, not data used for transition counts.
"""
import argparse
import csv
import json
import re
import sys
import time
import logging
import hashlib
from types import SimpleNamespace
from sklearn.metrics import pairwise_distances_argmin
from pathlib import Path

import numpy as np
from scipy.sparse.csgraph import connected_components
from sklearn.decomposition import PCA
from sklearn.cluster import MiniBatchKMeans
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


_STARTED = time.perf_counter()


def progress(message):
    logging.getLogger('shooting_msm').info('[%7.1fs] %s', time.perf_counter()-_STARTED, message)


def write_csv(path, header, rows):
    progress(f'Writing {path.name}')
    with open(path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)


def lag_frames(ps, dt):
    if not np.isfinite(dt) or dt <= 0:
        raise ValueError('Effective frame spacing must be finite and positive')
    if not np.isfinite(ps) or ps <= 0:
        raise ValueError(f'Lag {ps} ps must be finite and positive')
    ratio = ps / dt
    if not np.isfinite(ratio):
        raise ValueError('Lag/frame-spacing ratio is too large')
    n = int(round(ratio))
    if n < 1 or not np.isclose(n * dt, ps, rtol=1e-7, atol=1e-8):
        lower = max(1, int(np.floor(ratio)))
        upper = max(lower + 1, int(np.ceil(ratio)))
        raise ValueError(
            f'Lag {ps:g} ps must be a positive multiple of the effective frame '
            f'spacing {dt:g} ps (--dt-ps * --stride). '
            f'Nearby valid lags: {lower*dt:g} or {upper*dt:g} ps')
    return n


def counts(trajectories, lag, nstates):
    """Every lagged pair stays inside one input trajectory."""
    C = np.zeros((nstates, nstates), dtype=np.int64)
    for d in trajectories:
        if len(d) > lag:
            np.add.at(C, (d[:-lag], d[lag:]), 1)
    return C


def components(C):
    _, labels = connected_components(C > 0, directed=True, connection='strong')
    return [np.flatnonzero(labels == k) for k in np.unique(labels)]


def estimate(C, reversible=False):
    if len(C) < 2 or np.any(C.sum(axis=1) == 0) or len(components(C)) != 1:
        raise ValueError('At least two strongly connected states with outgoing counts required')
    if reversible:
        from deeptime.markov.tools.estimation import transition_matrix
        P, pi = transition_matrix(C.astype(float), reversible=True, return_statdist=True)
        P, pi = np.asarray(P), np.asarray(pi)
    else:
        P = C / C.sum(axis=1, keepdims=True)
        # Constrained stationary linear system (also handles periodic chains).
        A = np.vstack((P.T - np.eye(len(P)), np.ones(len(P))))
        b = np.r_[np.zeros(len(P)), 1.0]
        pi = np.linalg.lstsq(A, b, rcond=None)[0]
    if np.min(pi) < -1e-8 or not np.allclose(P.sum(axis=1), 1):
        raise ValueError('Invalid estimated stationary distribution or transition matrix')
    pi = np.maximum(pi, 0)
    pi /= pi.sum()
    if not np.allclose(pi @ P, pi, atol=1e-6):
        raise ValueError('Stationarity check failed')
    return P, pi


def relaxation(P, lag_ps, n=5):
    eig = np.linalg.eigvals(P)
    eig = np.delete(eig, np.argmin(abs(eig - 1)))
    eig = eig[np.argsort(-abs(eig))]
    out = []
    for v in eig[:n]:
        a = abs(v)
        t = -lag_ps / np.log(a) if 0 < a < 1 - 1e-12 else np.nan
        out.append((float(t), float(v.real), float(v.imag)))
    return out


def selected_coordinates(u, atoms, fit, args):
    """Yield (source frame, PCA xyz, fit xyz), retaining only needed atom batches.

    libdcd still reads full coordinate records into an internal single-frame
    buffer. This avoids full-system Python Timestep copies for every frame;
    it is not sparse disk I/O. Unwrapping requires the full bonded system.
    """
    if args.unwrap:
        for ts in u.trajectory[args.discard_frames::args.stride]:
            u.atoms.unwrap(compound='fragments', reference=None, inplace=True)
            yield ts.frame, atoms.positions, fit.positions
        return

    indices = np.union1d(atoms.indices, fit.indices)
    needed = u.atoms[indices]
    feature_map = np.searchsorted(indices, atoms.indices)
    fit_map = np.searchsorted(indices, fit.indices)
    # Bound the selected-coordinate batch to 128 frames / about 32 MiB.
    batch_frames = max(1, min(128, (32 * 1024**2) // (12 * len(indices))))
    total = len(u.trajectory)
    for start in range(args.discard_frames, total, batch_frames * args.stride):
        stop = min(total, start + batch_frames * args.stride)
        block = u.trajectory.timeseries(atomgroup=needed, start=start, stop=stop,
                                       step=args.stride, order='fac')
        expected = range(start, stop, args.stride)
        if len(block) != len(expected):
            raise ValueError('Incomplete coordinate batch read from DCD')
        for frame, xyz in zip(expected, block):
            yield frame, xyz[feature_map], xyz[fit_map]


def load_segments(args, out):
    import MDAnalysis as mda
    from MDAnalysis.analysis.align import rotation_matrix
    progress(f'Discovering segments under {args.output}')
    job = args.jobname or args.psf.stem
    rx = re.compile(re.escape(job) + r'\.(\d+)\.(\d+)\.dcd$')
    found = []
    for path in args.output.glob('*/*.dcd'):
        m = rx.fullmatch(path.name)
        if m and path.parent.name.isdigit():
            found.append((int(path.parent.name), int(m[1]), int(m[2]), path))
    found.sort(key=lambda x: x[:3])
    if not found:
        raise ValueError(f'No files match {args.output}/REPLICA/{job}.MILESTONE.SEGMENT.dcd')
    total_bytes = sum(item[3].stat().st_size for item in found)
    progress(f'Found {len(found)} segments; total DCD size {total_bytes/1024**2:.1f} MiB')
    arrays, metadata = [], []
    loading_start = time.perf_counter()
    loaded_frames = 0
    saved = getattr(args, "_model", None)
    ref = saved["reference"].copy() if saved is not None else None
    reference_center = saved["reference_center"].copy() if saved is not None else None
    u = None
    atoms = fit = None
    for number, (rep, milestone, segment, path) in enumerate(found, 1):
        file_start = time.perf_counter()
        progress(f'[{number}/{len(found)}] Opening DCD: {path}')
        before = path.stat()
        if u is None:
            progress(f'Loading PSF once: {args.psf}')
            u = mda.Universe(str(args.psf), str(path))
            fit = u.select_atoms(args.align_selection or args.selection)
            atoms = u.select_atoms(args.selection)
            if saved is not None:
                if not (np.array_equal(atoms.indices, saved['feature_indices']) and
                        np.array_equal(fit.indices, saved['alignment_indices'])):
                    raise ValueError('Saved model atom selections do not match topology')
                np.savez(out/'alignment_reference.npz', centered_coordinates=ref,
                         center=reference_center, alignment_indices=fit.indices, feature_indices=atoms.indices)
            mode = 'full-system frame processing for unwrapping' if args.unwrap else 'selected-atom batch extraction'
            progress(f'Reader mode: {mode}; topology and atom selections reused across segments')
        else:
            u.load_new(str(path))
        try:
            if len(fit) < 3 or len(atoms) < 1:
                raise ValueError('Need >=3 alignment atoms and >=1 PCA atom')
            progress(f'  {len(u.atoms)} total atoms; {len(atoms)} PCA atoms; '
                     f'{len(fit)} alignment atoms; {len(u.trajectory)} frames')
            frames = []
            last_update = time.perf_counter()
            frame_ids = []
            for frame, feature_xyz, fit_xyz in selected_coordinates(u, atoms, fit, args):
                xyz = fit_xyz.astype(np.float64)
                center = xyz.mean(axis=0)
                centered = xyz - center
                if ref is None:
                    if np.linalg.matrix_rank(centered) < 2:
                        raise ValueError('Alignment atoms are collinear')
                    ref = centered.copy()
                    reference_center = center.copy()
                    np.savez(out / 'alignment_reference.npz', centered_coordinates=ref,
                             center=reference_center, alignment_indices=fit.indices,
                             feature_indices=atoms.indices)
                R, _ = rotation_matrix(centered, ref)
                coords = (feature_xyz - center) @ R.T
                if not np.all(np.isfinite(coords)):
                    raise ValueError(f'Nonfinite aligned coordinates: {path}, frame {frame}')
                frames.append(coords.ravel().astype(np.float32))
                frame_ids.append(frame)
                if time.perf_counter()-last_update >= 5:
                    progress(f'  Reading/alignment: {len(frames)} retained frames; '
                             f'source frame {frame+1}/{len(u.trajectory)}')
                    last_update = time.perf_counter()
            if not frames:
                progress(f'Skipping {path}: no frames after discard')
                continue
            after = path.stat()
            if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
                raise ValueError(f'File changed during reading: {path}; use completed segments only')
            arrays.append(np.asarray(frames))
            loaded_frames += len(frames)
            elapsed = time.perf_counter()-loading_start
            eta = elapsed/number*(len(found)-number)
            progress(f'  Finished in {time.perf_counter()-file_start:.2f}s; '
                     f'{loaded_frames} frames loaded; approximate loading ETA {eta:.0f}s')
            metadata.append(dict(path=str(path.resolve()), replica=rep, milestone=milestone,
                                 segment=segment, frames=len(frames), raw_frames=len(u.trajectory),
                                 source_frames=frame_ids,
                                 duration_ps=(len(frames)-1)*args.dt_ps*args.stride))
        finally:
            u.trajectory.close()
    if not arrays:
        raise ValueError('No usable frames')
    progress(f'Loading complete in {time.perf_counter()-loading_start:.1f}s; '
             f'aligned coordinate arrays {sum(a.nbytes for a in arrays)/1024**2:.1f} MiB')
    return arrays, metadata


def plot_landscapes(out, Z, labels, active, pi, temperature, bins):
    H, xe, ye = np.histogram2d(Z[:, 0], Z[:, 1], bins=bins)
    fig, ax = plt.subplots(figsize=(7, 5))
    im = ax.pcolormesh(xe, ye, np.ma.masked_where(H.T == 0, np.log10(1+H.T)), shading='auto')
    fig.colorbar(im, ax=ax, label='log10(1 + frame count)')
    ax.set(xlabel='PC1 (Å)', ylabel='PC2 (Å)', title='PCA sampling density — not equilibrium probability')
    fig.tight_layout(); fig.savefig(out/'pca_sampling.png', dpi=180); plt.close(fig)
    if pi is None:
        return
    occupancy = np.bincount(labels, minlength=int(labels.max())+1)
    state_weights = np.zeros(len(occupancy))
    state_weights[active] = pi / occupancy[active]
    H, _, _ = np.histogram2d(Z[:, 0], Z[:, 1], bins=(xe, ye), weights=state_weights[labels])
    F = np.full_like(H, np.nan)
    mask = H > 0
    F[mask] = -0.008314462618 * temperature * np.log(H[mask]/H[mask].max())
    fig, ax = plt.subplots(figsize=(7, 5))
    im = ax.pcolormesh(xe, ye, np.ma.masked_invalid(F.T), shading='auto', cmap='viridis')
    fig.colorbar(im, ax=ax, label='Relative free energy (kJ/mol)')
    ax.set(xlabel='PC1 (Å)', ylabel='PC2 (Å)', title='MSM-weighted landscape — active component only')
    fig.tight_layout(); fig.savefig(out/'pca_free_energy.png', dpi=180); plt.close(fig)
    np.savez(out/'landscape.npz', x_edges=xe, y_edges=ye, probability=H, free_energy_kj_mol=F)


def plot_network(out, centers, C, active, pi, threshold, max_edges, routes=None):
    import networkx as nx
    G = nx.DiGraph()
    G.add_nodes_from(range(len(C)))
    edges = [(i,j,int(C[i,j])) for i,j in zip(*np.nonzero(C))
             if i != j and C[i,j] >= threshold]
    edges.sort(key=lambda e:e[2], reverse=True)
    edges = edges[:max_edges]
    G.add_weighted_edges_from(edges)
    pos = {i:centers[i,:2] for i in range(len(C))}
    population = np.zeros(len(C))
    if pi is not None:
        population[active] = pi
    colors = ['#237a9b' if i in set(active) else '#bbbbbb' for i in G.nodes]
    sizes = 90 + 1400*population/max(population.max(), 1e-12)
    fig, ax = plt.subplots(figsize=(9,7))
    nx.draw_networkx_nodes(G, pos, node_size=sizes, node_color=colors, ax=ax, alpha=.85)
    nx.draw_networkx_edges(G, pos, ax=ax, arrows=True, arrowsize=10,
        edgelist=[(i,j) for i,j,_ in edges],
        width=[.4 + np.log1p(e[2])*.35 for e in edges], alpha=.4,
        connectionstyle='arc3,rad=0.08', node_size=sizes)
    if routes is not None:
        from matplotlib.lines import Line2D
        import textwrap
        handles=[]
        for n,route in enumerate(routes['paths']):
            color=plt.get_cmap('tab10')(n % 10)
            pairs=list(zip(route['states'][:-1],route['states'][1:]))
            # Always show route edges, independently of background display filters.
            G.add_edges_from(pairs)
            nx.draw_networkx_edges(G,pos,edgelist=pairs,ax=ax,arrows=True,
                arrowsize=16,width=2.5,edge_color=[color],alpha=.9,node_size=sizes,
                connectionstyle=f'arc3,rad={.12+.09*n}')
            description=' -> '.join(map(str,route['states']))
            label=f"Path {n+1}: {100*route['fraction']:.1f}% of total flux\n"+textwrap.fill(description,width=45)
            handles.append(Line2D([0],[0],color=color,lw=2.5,label=label))
        for ids,shape,color,label in [(routes['start_states'],'s','#35a853','Start A'),
                                      (routes['end_states'],'D','#e6534c','End B')]:
            ids=[i for i in ids if i in G]
            if ids:
                nx.draw_networkx_nodes(G,pos,nodelist=ids,node_shape=shape,
                    node_color=color,node_size=sizes[ids],edgecolors='black',linewidths=1,ax=ax)
                handles.append(Line2D([0],[0],marker=shape,color='none',markerfacecolor=color,
                    markeredgecolor='black',markersize=8,label=label))
        if handles:
            ax.legend(handles=handles,loc='upper right',fontsize=8,framealpha=.95,
                title=f"Routes shown: {100*routes['covered_fraction']:.1f}% of total flux",
                title_fontsize=9)
    if len(C) <= 80:
        nx.draw_networkx_labels(G, pos, font_size=7, ax=ax)
    elif routes is not None:
        nodes=set(routes['start_states']+routes['end_states'])
        for route in routes['paths']:nodes.update(route['states'])
        nx.draw_networkx_labels(G,pos,labels={i:str(i) for i in nodes if i in G},font_size=7,ax=ax)
    ax.set_axis_on(); ax.tick_params(left=True, bottom=True, labelleft=True, labelbottom=True)
    ax.set(xlabel='PC1 (Å)', ylabel='PC2 (Å)',
           title='Observed transition network\nArrows: counts; node size: active-component MSM population')
    filename='transition_network.png'
    if routes is not None:
        filename='transition_paths.png'
        ax.set_title('Dominant A -> B reactive pathways (TPT)\nSame PCA nodes; coloured arrows show flux-decomposition routes')
        if routes['status'] != 'ok':
            ax.text(.02,.02,routes['message'],transform=ax.transAxes,fontsize=9,
                    bbox=dict(facecolor='white',alpha=.9),wrap=True)
    fig.tight_layout(); fig.savefig(out/filename, dpi=180); plt.close(fig)
    return len(edges)


def transition_paths(P, pi, active, C, start_states, end_states, fraction, max_paths, lag_ps):
    """TPT on the primary MSM; retain original cluster IDs in every output."""
    result=dict(status='unavailable',message='',start_states=sorted(set(start_states)),
                end_states=sorted(set(end_states)),paths=[],covered_fraction=0.0,
                requested_fraction=fraction,max_paths=max_paths,
                lag_ps=lag_ps,total_flux_per_lag=None,total_flux_per_ps=None)
    if P is None:
        result['message']='No connected MSM available for pathway analysis.'
        return result,None
    mapping={int(state):i for i,state in enumerate(active)}
    missing=sorted((set(start_states)|set(end_states))-set(mapping))
    if missing:
        result['message']=f'Endpoint states outside active component: {missing}. No routes estimated.'
        return result,None
    from deeptime.markov import reactive_flux
    flux=reactive_flux(P,[mapping[i] for i in result['start_states']],
                      [mapping[i] for i in result['end_states']],stationary_distribution=pi)
    total=float(flux.total_flux)
    if not np.isfinite(total) or total <= 0:
        result['message']='No positive finite A -> B reactive flux.'
        return result,None
    paths,capacities=flux.pathways(fraction=fraction,maxiter=max_paths)
    # Enforce the user-facing cap independently of backend iteration semantics.
    paths,capacities=paths[:max_paths],capacities[:max_paths]
    for path,capacity in zip(paths,capacities):
        states=[int(active[i]) for i in path]
        result['paths'].append(dict(states=states,flux_per_lag=float(capacity),
            flux_per_ps=float(capacity)/lag_ps,fraction=float(capacity)/total,
            unobserved_directed_edges=[[i,j] for i,j in zip(states[:-1],states[1:]) if C[i,j]==0]))
    result.update(status='ok',message='Dominant net-reactive-flux decomposition',
                  covered_fraction=float(np.sum(capacities)/total),
                  total_flux_per_lag=total,total_flux_per_ps=total/lag_ps)
    return result,flux


def export_paths(out, args, P, pi, active, C, centers):
    progress(f'TPT routes: states {args.start_states} -> {args.end_states}')
    result,flux=transition_paths(P,pi,active,C,args.start_states,args.end_states,
                                args.path_fraction,args.max_paths,args.lag_ps)
    result['interpretation']=(
        'Path fractions are shares of total A-to-B net reactive flux in this decomposition, '
        'not probabilities of exact time-resolved trajectories. Decompositions can be nonunique '
        'when routes have equal capacity. Applies to the active component only. '
        'No coarse-graining of microstates or automatic endpoint identification is performed.')
    result['state_id_note']='IDs are PCA-cluster states, NOT injection milestone IDs; they may change when fitting settings/data change.'
    (out/'paths_summary.json').write_text(json.dumps(result,indent=2)+'\n')
    cumulative=0.0;rows=[]
    for i,route in enumerate(result['paths'],1):
        cumulative+=route['fraction']
        rows.append([i,' -> '.join(map(str,route['states'])),route['flux_per_lag'],
                     route['flux_per_ps'],route['fraction'],cumulative,
                     '; '.join(f'{a}->{b}' for a,b in route['unobserved_directed_edges'])])
    write_csv(out/'paths.csv',['path','states','flux_per_lag','flux_per_ps',
              'fraction_of_total_flux','cumulative_fraction','edges_with_zero_observed_count'],rows)
    if flux is not None:
        write_csv(out/'committors.csv',['state','forward_committor','backward_committor'],
                  zip(active,flux.forward_committor,flux.backward_committor))
        write_csv(out/'reactive_edges.csv',['source','destination','net_flux_per_lag',
                  'net_flux_per_ps','observed_count'],
                  ((active[i],active[j],flux.net_flux[i,j],flux.net_flux[i,j]/args.lag_ps,
                    C[active[i],active[j]]) for i,j in zip(*np.nonzero(flux.net_flux))))
    plot_network(out,centers,C,active,pi,args.min_edge_count,args.max_edges,routes=result)
    progress(f"Path analysis: {result['status']}; {len(result['paths'])} routes; "
             f"{result['covered_fraction']:.1%} of total reactive flux displayed")
    if result['status']!='ok':progress(result['message'])
    return result


def load_fitted_model(args):
    """Portable numeric archive: never deserialize executable pickle objects."""
    with np.load(args.load_model, allow_pickle=False) as archive:
        model = {name: archive[name] for name in archive.files}
    if int(model['schema_version']) != 1:
        raise ValueError('Unsupported fitted-model schema version')
    if str(model['psf_sha256']) != hashlib.sha256(args.psf.read_bytes()).hexdigest():
        raise ValueError('PSF differs from saved model; identical topology file is required')
    pcs = model['pca_components']; centers = model['cluster_centers']
    if pcs.ndim != 2 or centers.ndim != 2 or centers.shape[1] != len(pcs):
        raise ValueError('Invalid saved PCA/cluster dimensions')
    for name in ['pca_components','pca_mean','cluster_centers','reference','reference_center']:
        if not np.isfinite(model[name]).all():
            raise ValueError('Nonfinite values in fitted model: ' + name)
    args.selection = str(model['selection'])
    args.align_selection = str(model['align_selection'])
    args.unwrap = bool(model['unwrap'])
    args.pcs, args.clusters = len(pcs), len(centers)
    args._model = model


def export_structures(args, out, representatives, routes):
    """Actual sampled representatives, aligned into the fitted PCA reference."""
    import MDAnalysis as mda
    from MDAnalysis.analysis.align import rotation_matrix
    progress('Exporting representative node and pathway PDBs')
    with np.load(out/'alignment_reference.npz', allow_pickle=False) as archive:
        ref = archive['centered_coordinates']
        fit_indices = archive['alignment_indices']
    node_dir = out/'node_pdbs'; node_dir.mkdir()
    u = mda.Universe(str(args.psf))
    selected = u.select_atoms(args.pdb_selection)
    if not len(selected):
        raise ValueError('--pdb-selection matches no atoms')
    current = None; node_files = {}
    for state, dcd, frame in representatives:
        if current != dcd:
            u.load_new(dcd); current = dcd
        u.trajectory[int(frame)]
        if args.unwrap:
            u.atoms.unwrap(compound='fragments', reference=None, inplace=True)
        xyz = u.atoms.positions.astype(np.float64)
        center = xyz[fit_indices].mean(axis=0)
        R, _ = rotation_matrix(xyz[fit_indices]-center, ref)
        u.atoms.positions = (xyz-center) @ R.T
        filename = node_dir/f'state_{state:04d}.pdb'
        with mda.Writer(str(filename), multiframe=False, bonds=None) as writer:
            writer.write(selected)
        node_files[state] = filename
    if current is not None:
        u.trajectory.close()
    manifest = []
    if routes is not None and routes['paths']:
        path_dir = out/'path_pdbs'; path_dir.mkdir()
        for number, route in enumerate(routes['paths'], 1):
            filename = path_dir/f'path_{number:03d}.pdb'
            # Reuse the exact node structures; each MODEL is a state, not a time frame.
            with filename.open('w') as handle:
                handle.write('REMARK 900 ORDERED STATE REPRESENTATIVES; NOT A CONTINUOUS TRAJECTORY\n')
                handle.write(f"REMARK 900 FRACTION OF TOTAL REACTIVE FLUX {route['fraction']:.8g}\n")
                for model, state in enumerate(route['states'], 1):
                    handle.write(f'REMARK 900 MODEL {model} STATE {state}\nMODEL     {model:4d}\n')
                    for line in node_files[state].read_text().splitlines():
                        if line.startswith(('ATOM  ', 'HETATM', 'TER   ')):
                            handle.write(line+'\n')
                    handle.write('ENDMDL\n')
                    manifest.append([number, model, state, str(filename.relative_to(out)),
                                     str(node_files[state].relative_to(out))])
                handle.write('END\n')
        write_csv(out/'path_structures.csv', ['path','model','state','path_pdb','node_pdb'], manifest)
    (node_dir/'README.txt').write_text(
        'Each PDB is the sampled frame nearest its cluster center in retained PCA space.\n'
        'Source DCD and zero-based frame index: ../representatives.csv.\n'
        'Coordinates use the centered alignment reference; selection: '+args.pdb_selection+'\n'
        'Unsampled states have no representative. Path MODEL order is state order, not time.\n')


def assess_quality(out, args, C, active, pi, dtrajs, lengths, occupancy, spectral, ck_rows, routes):
    """Descriptive screening rules, not calibrated confidence or hypothesis tests."""
    lag = lag_frames(args.lag_ps, args.dt_ps * args.stride)
    thresholds = dict(min_active_frame_fraction=0.90, min_pairs_per_state=100,
        min_exit_pairs_per_state=20, min_segments_with_exits_per_state=3,
        min_valid_lags=3, min_lag_span_ratio=2.0, max_relative_timescale_range=0.20,
        min_ck_pairs_per_row=50, min_ck_population_coverage=0.80,
        min_ck_longer_horizons=2, max_ck_row_total_variation=0.15,
        min_path_edge_supporting_segments=3)
    findings=[]; concerns=[]; missing=[]; actions=[]
    def flag(text, action):
        concerns.append(text)
        if action not in actions: actions.append(action)
    exists = pi is not None and len(active)>1
    coverage=float(occupancy[active].sum()/max(1,occupancy.sum()))
    contributes=int(np.sum(lengths>lag))
    metrics=dict(active_frame_fraction=coverage, contributing_segments=contributes,
                 total_segments=len(lengths), selected_lag_ps=args.lag_ps)
    findings.append(f'{contributes}/{len(lengths)} segments are long enough to contribute transition pairs at {args.lag_ps:g} ps.')
    findings.append(f'The active connected model contains {len(active)}/{len(C)} states and {coverage:.1%} of sampled frames. This is sampling coverage, not equilibrium coverage.')
    if not exists:
        verdict='NO ESTIMABLE CONNECTED KINETIC MODEL'
        summary=('Do not trust kinetic times, equilibrium free energies or transition-route weights from this dataset yet. '
                 'The structural projection and representative PDBs can still describe what was sampled.')
        actions.append('Collect continuous baseline segments long enough to observe transitions connecting multiple states; inspect lag_coverage.csv.')
    else:
        restricted=C[np.ix_(active,active)]
        pairs=restricted.sum(axis=1); exits=pairs-np.diag(restricted)
        support=np.zeros(len(active),int)
        edge_support=np.zeros_like(C)
        for d in dtrajs:
            if len(d)<=lag: continue
            source,target=d[:-lag],d[lag:]
            observed_edges=np.unique(np.column_stack((source,target)),axis=0)
            edge_support[observed_edges[:,0],observed_edges[:,1]]+=1
            for j,state in enumerate(active):
                support[j]+=bool(np.any((source==state)&(target!=state)&np.isin(target,active)))
        outgoing=int(C[active].sum()-restricted.sum())
        weak=active[(pairs<thresholds['min_pairs_per_state'])|
                     (exits<thresholds['min_exit_pairs_per_state'])|
                     (support<thresholds['min_segments_with_exits_per_state'])].tolist()
        metrics.update(excluded_outgoing_pairs=outgoing,weakly_supported_states=weak,
            state_support=[dict(state=int(state),pairs=int(pairs[j]),exit_pairs=int(exits[j]),
                                segments_with_exits=int(support[j])) for j,state in enumerate(active)])
        findings.append(f'The least-supported active state has {int(pairs.min())} outgoing pairs, {int(exits.min())} exit pairs and exits in {int(support.min())} segments (these minima can concern different states).')
        if coverage<thresholds['min_active_frame_fraction']:
            flag('A substantial part of the sampled data lies outside the active component; the MSM cannot describe the whole sampled ensemble.',
                 'Improve connectivity between sampled regions; interpret populations and free energies only within the active component.')
        if outgoing:
            flag(f'{outgoing} observed outgoing pairs leave the active component and are omitted from the fitted model.',
                 'Sample returns from excluded regions before interpreting the restricted model as an equilibrium ensemble.')
        if weak:
            flag(f'Active states with limited transition support: {weak}. Sliding pairs may reflect repeated observations of the same event.',
                 'Collect more baseline transitions from the weakly supported states, preferably across independent runs.')
        # Magnitude-ranked spectra provide a screen only, not tracked physical modes.
        by_lag={}
        for t,mode,scale,real,imag in spectral:
            if t>=args.lag_ps and mode<=min(3,len(active)-1):
                by_lag.setdefault(t,[]).append((scale,real,imag))
        valid={t:v for t,v in by_lag.items() if len(v)==min(3,len(active)-1) and
               all(np.isfinite(x[0]) and x[0]>0 and x[1]>0 and abs(x[2])<1e-8 for x in v)}
        metrics['usable_lags_for_timescale_screen']=sorted(valid)
        if (len(valid)<thresholds['min_valid_lags'] or
                max(valid,default=0)/max(min(valid,default=1),1e-30)<thresholds['min_lag_span_ratio']):
            missing.append('Lag stability is not established: need at least three usable lags at/above the selected lag spanning a factor of two, with positive real slow modes.')
            actions.append('Add compatible --lags-ps values at and above the chosen lag; ensure segments are long enough for them.')
        else:
            values=np.array([[x[0] for x in valid[t]] for t in sorted(valid)])
            spreads=np.ptp(values,axis=0)/np.median(values,axis=0)
            metrics['relative_slow_timescale_ranges']=spreads.tolist()
            findings.append(f'The largest relative range of the first {values.shape[1]} slow timescales across usable lags is {spreads.max():.1%}. Modes are ranked, not structurally matched.')
            if spreads.max()>thresholds['max_relative_timescale_range']:
                flag('Slow timescales vary appreciably with lag; a stable kinetic regime has not been demonstrated.',
                     'Inspect implied_timescales.png and try longer segments, alternative lags or a better state representation.')
        ck_qualified=[]; ck_details=[]
        for multiple,horizon,_,_,_ in ck_rows:
            if multiple<=1: continue
            with np.load(out/f'ck_{multiple}.npz',allow_pickle=False) as ck:
                enough=ck['row_counts']>=thresholds['min_ck_pairs_per_row']
                pop=float(pi[enough].sum())
                tv=0.5*np.nansum(abs(ck['predicted']-ck['observed']),axis=1)
                worst=float(tv[enough].max()) if enough.any() else None
                ck_details.append(dict(horizon_ps=float(horizon),supported_population=pop,
                                       max_supported_row_total_variation=worst))
                if pop>=thresholds['min_ck_population_coverage'] and worst is not None: ck_qualified.append(worst)
        metrics['ck_longer_horizons']=ck_details
        if len(ck_qualified)<thresholds['min_ck_longer_horizons']:
            missing.append('CK validation lacks two longer prediction horizons with at least 50 observed pairs per evaluated row covering 80% of active-model population.')
            actions.append('Collect longer held-out segments so CK can test predictions beyond the fitting lag.')
        if ck_qualified:
            findings.append(f'Maximum supported-row CK total-variation distance is {max(ck_qualified):.3f} (0 means agreement; 1 means disjoint distributions).')
            if max(ck_qualified)>thresholds['max_ck_row_total_variation']:
                flag('Held-out predictions disagree appreciably with observed transitions under the screening rule; sampling noise and model error are not separated.',
                     'Inspect ck_validation.png; increase independent sampling and test lag/state choices before reporting kinetic times.')
        if routes is not None:
            if routes['status']!='ok':
                findings.append('Requested pathway analysis is unavailable: '+routes['message'])
            else:
                weak_edges=[]
                for route in routes['paths']:
                    for i,j in zip(route['states'][:-1],route['states'][1:]):
                        if edge_support[i,j]<thresholds['min_path_edge_supporting_segments']: weak_edges.append((int(i),int(j)))
                weak_edges=sorted(set(weak_edges))
                metrics['path_edges_with_fewer_than_three_supporting_segments']=weak_edges
                findings.append(f'Displayed routes cover {routes["covered_fraction"]:.1%} of fitted reactive flux; this is not a confidence level.')
                if weak_edges:
                    findings.append(f'Route edges seen in fewer than three segments: {weak_edges}. Route rankings are particularly tentative; zero-count edges can arise from reversible estimation.')
                    actions.append('Sample the weakly supported pathway edges before interpreting route rankings quantitatively.')
        if concerns:
            verdict='LIMITED RELIABILITY — EXPLORATORY USE ONLY'
            summary='The model can organise sampled structures and suggest candidate routes, but the diagnostics flag limitations. Do not yet treat kinetic times, free-energy differences or pathway rankings as quantitatively reliable.'
        elif missing:
            verdict='VALIDATION INCOMPLETE — TRUST NOT YET ESTABLISHED'
            summary='No warning threshold was crossed in the available checks, but key validation evidence is missing. Use the model for exploration; quantitative kinetics remain provisional.'
        else:
            verdict='ENCOURAGING DIAGNOSTICS — STILL PROVISIONAL'
            summary='Connectivity, transition support, lag sensitivity and CK agreement meet the screening rules. This supports further interpretation within the sampled component, but does not establish convergence or validated physical kinetics.'
            actions.append('Check uncertainty, independent-run reproducibility and sensitivity to clustering/PCA before quantitative claims.')
    limitations=[
        'The thresholds below are transparent heuristic screening rules, not statistical confidence levels or universal acceptance criteria.',
        'Sliding transition pairs and exchange-coupled segments are correlated; segment counts are not independent-run counts.',
        'These diagnostics cannot establish unbiased within-state sampling or rule out injection and exchange-conditioned stopping bias.',
        'No confidence intervals, convergence-versus-runtime test or independent-run validation are currently available. High PCA variance capture alone does not validate kinetics.']
    if args.load_model is not None:
        limitations.append('The reused PCA/clustering model may have seen held-out frames; CK is not fully independent representation validation.')
    result=dict(verdict=verdict,summary=summary,findings=findings,concerns=concerns,
                missing_evidence=missing,recommended_actions=actions,limitations=limitations,
                screening_thresholds=thresholds,metrics=metrics)
    paragraphs=['QUALITY ASSESSMENT',verdict,summary]
    for title,items in [('Evidence',findings),('Concerns',concerns),('Missing evidence',missing),
                        ('Next steps',actions),('Limits of this assessment',limitations)]:
        if items: paragraphs.extend(['',title]+['- '+x for x in items])
    paragraphs.extend(['','Screening thresholds (heuristics):'])
    paragraphs.extend('- '+name.replace('_',' ')+': '+str(value) for name,value in thresholds.items())
    text='\n'.join(paragraphs)+'\n'
    (out/'quality_report.txt').write_text(text)
    (out/'quality_report.json').write_text(json.dumps(result,indent=2)+'\n')
    return result,text


def analyze(args, arrays, metadata, out):
    rng = np.random.default_rng(args.seed)
    lengths = np.array([len(x) for x in arrays])
    offsets = np.r_[0, np.cumsum(lengths)]
    progress('Assembling aligned coordinate matrix')
    X = np.concatenate(arrays)
    progress(f'Coordinate matrix: {X.shape[0]} frames x {X.shape[1]} features')
    dt = args.dt_ps * args.stride
    target_lag = lag_frames(args.lag_ps, dt)
    lags = sorted(set([target_lag]+[lag_frames(v, dt) for v in args.lags_ps]))
    # Whole-segment holdout. Not an independent-run split.
    ids = rng.permutation(len(arrays))
    ntest = max(1, int(round(len(ids)*args.test_fraction))) if len(ids)>=3 and args.test_fraction>0 else 0
    test_ids = set(ids[:ntest].tolist())
    training = np.concatenate([np.arange(offsets[i], offsets[i+1])
                              for i in range(len(arrays)) if i not in test_ids])
    if len(training)>args.fit_max:
        training = np.sort(rng.choice(training,args.fit_max,replace=False))
    saved = getattr(args, '_model', None)
    if saved is None:
        ncomp = min(args.pcs, len(training)-1, X.shape[1])
        if ncomp < 2:
            raise ValueError('Need enough training frames/features for >=2 PCs')
        pca = PCA(n_components=ncomp, svd_solver='randomized', random_state=args.seed)
        stage_start = time.perf_counter()
        progress(f'Fitting PCA: {len(training)} training frames, {ncomp} components')
        pca.fit(X[training])
        progress(f'PCA fit complete in {time.perf_counter()-stage_start:.1f}s; projecting all frames')
        Z = pca.transform(X)
        progress(f'PCA projection complete; retained variance {pca.explained_variance_ratio_.sum():.1%}')
        if not np.isfinite(pca.explained_variance_ratio_).all() or pca.explained_variance_.sum()<=0:
            raise ValueError('No coordinate variance after alignment')
        k = min(args.clusters,len(training))
        if k < 2:
            raise ValueError('Need >=2 clusters')
        km = MiniBatchKMeans(n_clusters=k, random_state=args.seed, n_init=10,
                            batch_size=min(4096,max(256,len(training))), reassignment_ratio=0)
        stage_start = time.perf_counter()
        progress(f'Clustering: {k} centers in {ncomp} PCs, {len(training)} training frames')
        km.fit(Z[training])
        progress(f'Clustering fit complete in {time.perf_counter()-stage_start:.1f}s; assigning all frames')
        all_labels = km.predict(Z)
    else:
        progress('Reusing saved alignment, PCA and cluster centers; no refitting')
        pca = SimpleNamespace(mean_=saved['pca_mean'], components_=saved['pca_components'],
                              explained_variance_ratio_=saved['explained_variance_ratio'])
        km = SimpleNamespace(cluster_centers_=saved['cluster_centers'])
        ncomp, k = len(pca.components_), len(km.cluster_centers_)
        Z = (X - pca.mean_) @ pca.components_.T
        all_labels = pairwise_distances_argmin(Z, km.cluster_centers_)
    with np.load(out/'alignment_reference.npz', allow_pickle=False) as reference:
        np.savez_compressed(out/'fitted_model.npz', schema_version=np.array(1),
            psf_sha256=np.array(hashlib.sha256(args.psf.read_bytes()).hexdigest()),
            selection=np.array(args.selection), align_selection=np.array(args.align_selection or args.selection),
            unwrap=np.array(args.unwrap), reference=reference['centered_coordinates'],
            reference_center=reference['center'], alignment_indices=reference['alignment_indices'],
            feature_indices=reference['feature_indices'], pca_mean=pca.mean_,
            pca_components=pca.components_, explained_variance_ratio=pca.explained_variance_ratio_,
            cluster_centers=km.cluster_centers_)
    progress('State assignments complete; saving PCA and trajectory tables')
    dtrajs = [all_labels[offsets[i]:offsets[i+1]] for i in range(len(arrays))]
    for i,m in enumerate(metadata):
        m['split'] = 'test' if i in test_ids else 'train'
    write_csv(out/'segments.csv', ['path','replica','milestone','segment','frames','duration_ps','split'],
              [[m[h] for h in ['path','replica','milestone','segment','frames','duration_ps','split']] for m in metadata])
    write_csv(out/'frames.csv', ['segment_index','source_frame','time_since_segment_start_ps','state','PC1_A','PC2_A'],
              ((i,f,f*args.dt_ps,int(dtrajs[i][j]),*Z[offsets[i]+j,:2])
               for i,m in enumerate(metadata) for j,f in enumerate(m['source_frames'])))
    np.savez_compressed(out/'pca_and_states.npz', scores=Z, state_labels=all_labels, offsets=offsets,
        pca_mean=pca.mean_, pca_components=pca.components_,
        explained_variance_ratio=pca.explained_variance_ratio_, cluster_centers=km.cluster_centers_)
    progress('Plotting PCA variance and segment-length diagnostics')
    fig,axes=plt.subplots(1,2,figsize=(10,4))
    axes[0].plot(np.arange(1,ncomp+1), np.cumsum(pca.explained_variance_ratio_), 'o-')
    axes[0].set(xlabel='Number of PCs',ylabel='Cumulative explained variance',ylim=(0,1.02),
                title='Original fitted PCA variance' if saved is not None else 'Fitted PCA variance')
    axes[1].hist((lengths-1)*dt,bins=min(30,len(lengths)))
    axes[1].axvline(args.lag_ps,color='red',ls='--',label='Selected lag')
    axes[1].set(xlabel='Retained segment duration (ps)',ylabel='Segments');axes[1].legend()
    fig.tight_layout();fig.savefig(out/'sampling_diagnostics.png',dpi=180);plt.close(fig)
    progress(f'Counting within-segment transitions at lag {args.lag_ps:g} ps')
    C=counts(dtrajs,target_lag,k)
    occupancy=np.bincount(all_labels,minlength=k)
    comps=components(C)
    valid=[s for s in comps if len(s)>=2 and C[np.ix_(s,s)].sum()>0]
    active=max(valid,key=lambda s:occupancy[s].sum()) if valid else np.array([],dtype=int)
    component_ids=np.empty(k,dtype=int)
    for c,s in enumerate(comps):component_ids[s]=c
    progress(f'Observed pairs: {C.sum()}; strongly connected components: {len(comps)}; '
             f'selected active states: {len(active)}/{k}')
    progress('Estimating selected-lag MSM' if len(active) else 'No connected MSM available; retaining sampling diagnostics')
    P=pi=None
    if len(active):P,pi=estimate(C[np.ix_(active,active)],args.reversible)
    np.savez_compressed(out/'msm.npz', counts=C, active_states=active,
        transition_matrix=P if P is not None else np.empty((0,0)),
        stationary_distribution=pi if pi is not None else np.empty(0), lag_ps=args.lag_ps)
    stat=np.full(k,np.nan)
    if pi is not None:stat[active]=pi
    write_csv(out/'states.csv',['state','frames','observed_component','active','stationary_probability','PC1_A','PC2_A'],
        ((i,occupancy[i],component_ids[i],i in active,stat[i],*km.cluster_centers_[i,:2]) for i in range(k)))
    write_csv(out/'transitions.csv',['source','destination','observed_pairs','supporting_segments'],
        ((i,j,C[i,j],sum(bool(np.any((d[:-target_lag]==i)&(d[target_lag:]==j)))
                        for d in dtrajs if len(d)>target_lag)) for i,j in zip(*np.nonzero(C))))
    progress('Plotting PCA density and MSM-weighted free energy (if identifiable)')
    plot_landscapes(out,Z,all_labels,active,pi,args.temperature,args.bins)
    progress('Plotting transition network')
    shown=plot_network(out,km.cluster_centers_,C,active,pi,args.min_edge_count,args.max_edges)
    path_result=None
    if args.start_states is not None:
        path_result=export_paths(out,args,P,pi,active,C,km.cluster_centers_)
    # Use the same active states at every lag; do not quietly switch components.
    lag_rows=[]; spectral=[]
    for l in lags:
        progress(f'Lag scan: {l*dt:g} ps')
        Cl=counts(dtrajs,l,k)
        row=[l*dt,int(np.sum(lengths>l)),int(Cl.sum()),int(Cl.sum()-np.trace(Cl))]
        modes=[]
        if len(active):
            try:
                Pl,_=estimate(Cl[np.ix_(active,active)],args.reversible)
                modes=relaxation(Pl,l*dt)
            except ValueError as exc:
                progress(f'Lag {l*dt:g} ps: model unavailable ({exc})')
        for m in range(5):
            t,re_,im_=modes[m] if m<len(modes) else (np.nan,np.nan,np.nan)
            spectral.append([l*dt,m+1,t,re_,im_])
        lag_rows.append(row)
    write_csv(out/'lag_coverage.csv',['lag_ps','contributing_segments','pairs','off_diagonal_pairs'],lag_rows)
    write_csv(out/'implied_timescales.csv',['lag_ps','mode_by_magnitude','timescale_ps','eigenvalue_real','eigenvalue_imag'],spectral)
    fig,ax=plt.subplots(figsize=(7,5));sp=np.array(spectral)
    for m in range(1,6):
        v=sp[sp[:,1]==m];ax.plot(v[:,0],v[:,2],'o-',label=f'Mode {m}')
    ax.plot([lags[0]*dt,lags[-1]*dt],[lags[0]*dt,lags[-1]*dt],'k:',label='timescale = lag')
    ax.set(xlabel='Lag (ps)',ylabel='Relaxation / decay-envelope timescale (ps)',
           title='Fixed active states; gaps indicate insufficient connectivity')
    ax.set_xscale('log');ax.set_yscale('log');ax.legend(fontsize=8)
    fig.tight_layout();fig.savefig(out/'implied_timescales.png',dpi=180);plt.close(fig)
    # Save nearest sampled conformations, rather than unphysical PCA centroids.
    representative_rows=[]
    for state in range(k):
        members=np.flatnonzero(all_labels==state)
        if not len(members):continue
        f=int(members[np.argmin(np.sum((Z[members]-km.cluster_centers_[state])**2,axis=1))])
        si=int(np.searchsorted(offsets,f,side='right')-1)
        representative_rows.append([state,metadata[si]['path'],metadata[si]['source_frames'][f-offsets[si]]])
    write_csv(out/'representatives.csv',['state','dcd','source_frame'],representative_rows)
    export_structures(args, out, representative_rows, path_result)
    # Held-out CK: train on train segments, compare directly observed test probabilities.
    ck_status='Not available: insufficient training connectivity or no held-out segments'
    progress('Starting held-out Chapman-Kolmogorov diagnostic')
    ck_rows=[]
    if len(active) and test_ids:
        train=[d for i,d in enumerate(dtrajs) if i not in test_ids]
        test=[d for i,d in enumerate(dtrajs) if i in test_ids]
        try:
            Pt,_=estimate(counts(train,target_lag,k)[np.ix_(active,active)],args.reversible)
            for mult in [1,2,3,5]:
                progress(f'  CK prediction horizon: {args.lag_ps*mult:g} ps')
                obs=counts(test,target_lag*mult,k)[np.ix_(active,active)]
                totals=obs.sum(axis=1);ok=totals>0
                if not ok.any():continue
                predicted=np.linalg.matrix_power(Pt,mult)
                empirical=np.full_like(predicted,np.nan)
                empirical[ok]=obs[ok]/totals[ok,None]
                rmse=float(np.sqrt(np.mean((predicted[ok]-empirical[ok])**2)))
                ck_rows.append([mult,args.lag_ps*mult,int(obs.sum()),int(ok.sum()),rmse])
                np.savez(out/f'ck_{mult}.npz',states=active,predicted=predicted,observed=empirical,row_counts=totals)
            if ck_rows:
                ck_status='Segment-held-out diagnostic available; heuristic assessment only, no confidence intervals'
                fig,axes=plt.subplots(1,2,figsize=(11,4));r=np.array(ck_rows)
                axes[0].plot(r[:,1],r[:,4],'o-');axes[0].set(xlabel='Prediction horizon (ps)',ylabel='Probability RMSE',
                    title='Held-out CK diagnostic (no uncertainty bounds)')
                for si in np.argsort(-pi)[:3]:
                    pred=[];obs=[]
                    for multiple in r[:,0].astype(int):
                        with np.load(out/f'ck_{multiple}.npz') as v:
                            pred.append(v['predicted'][si,si]);obs.append(v['observed'][si,si])
                    line,=axes[1].plot(r[:,1],pred,'-',label=f'State {active[si]} predicted')
                    axes[1].plot(r[:,1],obs,'o--',color=line.get_color(),label=f'State {active[si]} observed')
                axes[1].set(xlabel='Prediction horizon (ps)',ylabel='Return probability',ylim=(-.03,1.03))
                axes[1].legend(fontsize=7)
                fig.tight_layout();fig.savefig(out/'ck_validation.png',dpi=180);plt.close(fig)
        except ValueError as exc:
            ck_status=f'Not available: {exc}'
            progress(f'CK diagnostic unavailable: {exc}')
    write_csv(out/'ck_validation.csv',['multiple','horizon_ps','pairs','observed_rows','rmse'],ck_rows)
    outgoing=int(C[active].sum()-C[np.ix_(active,active)].sum()) if len(active) else 0
    report=dict(segments=len(arrays),frames=len(X),retained_baseline_ps=float(np.sum((lengths-1)*dt)),
        frame_interval_ps=dt,selected_lag_ps=args.lag_ps,pca_components=ncomp,
        pca_variance_captured=float(pca.explained_variance_ratio_.sum()),clusters=k,
        observed_strong_components=len(comps),active_states=active.tolist(),
        active_frame_fraction=float(occupancy[active].sum()/len(X)),
        outgoing_pairs_excluded_from_active_model=outgoing,
        observed_pairs=int(C.sum()),off_diagonal_pairs=int(C.sum()-np.trace(C)),
        estimator='reversible MLE' if args.reversible else 'nonreversible row normalization',
        ck_status=ck_status,network_edges_displayed=shown,path_analysis=path_result,
        options={k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items() if not k.startswith('_')})
    if args.load_model is not None:
        ck_status += '; reused feature model: representation may have seen held-out frames'
        report['ck_status'] = ck_status
    progress('Assessing result quality and remaining validation gaps')
    quality,quality_text=assess_quality(out,args,C,active,pi,dtrajs,lengths,occupancy,spectral,ck_rows,path_result)
    report['quality_assessment']=quality
    progress('Writing model report')
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    status='MSM ESTIMATED — NOT AUTOMATICALLY VALIDATED' if P is not None else 'INSUFFICIENT OBSERVED CONNECTIVITY FOR AN MSM'
    lines=[status,'',f'Segments: {len(arrays)}; frames: {len(X)}; lag: {args.lag_ps:g} ps',
           f'Active states: {len(active)}/{k}; active frame coverage: {report["active_frame_fraction"]:.1%}',
           f'Observed lagged pairs: {C.sum()}; off-diagonal pairs: {report["off_diagonal_pairs"]}',
           f'Excluded outgoing pairs from active component: {outgoing}',f'CK: {ck_status}','',
           'INTERPRETATION',
           '- Raw PCA density reflects the injection/sampling policy, not equilibrium populations.',
           '- FES is conditional on the active component and assumes adequate within-state sampling.',
           '- Excluded components have unknown relative weights; omitted exits can bias the restricted model.',
           '- Sliding pairs and exchange-coupled segments are correlated; counts are not independent samples.',
           '- Eigenmodes are sorted by magnitude, not tracked by structural identity across lags.',
           '- Nonreal/negative eigenvalues describe oscillatory decay; inspect eigenvalues in the CSV.',
           '- CK is a noisy segment-held-out diagnostic, not independent-run validation; no CIs are provided.',
           '- Non-equilibrium starts and exchange-conditioned termination can bias estimates.',
           '- Inspect lag stability, CK predictions, connectivity and independent runs before kinetic claims.',
           '- No automatic macrostate assignment or endpoint identification; TPT requires explicit endpoint states.',
           '', 'FILES',
           'pca_sampling.png; pca_free_energy.png (if identifiable); transition_network.png;',
           'sampling_diagnostics.png; implied_timescales.png; ck_validation.png (if testable);',
           'states.csv; transitions.csv; segments.csv; frames.csv; lag_coverage.csv;',
           'fitted_model.npz; node_pdbs/; path_pdbs/ (when routes requested);',
           'representatives.csv; msm.npz; pca_and_states.npz; landscape.npz (if identifiable);',
           'quality_report.txt; quality_report.json; report.json; analysis.log.']
    if path_result is not None:
        lines += ['', 'PATHWAY ANALYSIS',path_result['message'],
                  f"Displayed {len(path_result['paths'])} paths covering {path_result['covered_fraction']:.1%} of total reactive flux.",
                  'See transition_paths.png, paths.csv and paths_summary.json; committors.csv and reactive_edges.csv when available.',
                  'Flux shares are decomposition weights, not probabilities of exact microscopic trajectories.']
    lines += ['',quality_text]
    (out/'report.txt').write_text('\n'.join(lines)+'\n')
    print('\n'.join(lines[:8]).rstrip())
    print('\n'+quality['verdict']+'\n'+quality['summary'])
    # The first two findings repeat the run statistics printed immediately above.
    for item in quality['findings'][2:]+quality['concerns']+quality['missing_evidence']:
        print('- '+item)
    for action in quality['recommended_actions']:
        print('Next: '+action)
    print('Full assessment and heuristic thresholds: '+str(out/'quality_report.txt'))


def build_parser():
    p = argparse.ArgumentParser(
        description=(
            'Build a coordinate-PCA Markov state model from separate baseline DCD segments.\n'
            'Counts never cross segment boundaries. Input names must be:\n'
            '  INPUT/REPLICA/JOB.MILESTONE.SEGMENT.dcd\n'
            'REPLICA, MILESTONE and SEGMENT are numeric; temporary DCDs are excluded.'),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Examples:
  python shooting_msm.py ng_beta36_orig.psf --jobname ng_beta36 \\
      --output output --out msm_results --dt-ps 16 --lag-ps 96 \\
      --lags-ps 16 32 64 96 160 320 --pcs 5 --clusters 20

  # Fit globally using one stable domain, but analyse all protein C-alpha atoms:
  python shooting_msm.py job.psf --dt-ps 20 --lag-ps 100 \\
      --align-selection "protein and name CA and resid 1:60"

Pathways:
  Add --start-states 0 --end-states 7 --max-paths 5 --path-fraction 0.9
  to overlay dominant A-to-B routes on a second copy of the network graph.
  You may specify sets, e.g. --start-states 0 1 --end-states 7 8.
  Choose IDs from states.csv/representatives.csv, not milestone filenames.
  Use --load-model previous_results/fitted_model.npz to preserve state IDs as
  data accumulate. Without it, refitting can change their meaning.
  Legends give fractions of TOTAL reactive flux, not percentages renormalized
  over displayed routes. The route limit may prevent reaching the target fraction.
  Every intermediate microstate is retained; no automatic macrostate merging.
  Original count-network filters do not remove coloured route edges.

Reusable model and PDBs:
  Every run saves fitted_model.npz (numeric data; no executable pickle).
  Add --load-model msm_results/fitted_model.npz --out msm_updated to reuse it.
  Supply all completed DCDs for cumulative analysis; old counts are not appended.
  The identical PSF is required. Saved selections, unwrapping, alignment, PCA and
  centers override fitting settings. Timing and MSM estimation remain configurable.
  Variance percentages describe the original PCA fit, not the new dataset.
  node_pdbs/state_XXXX.pdb contains the nearest currently sampled representative.
  Default export selection is protein; --pdb-selection "all" includes solvent.
  path_pdbs/path_XXX.pdb contains one MODEL per ordered node representative.
  These are structural route illustrations, NOT continuous transition trajectories.
  path_structures.csv maps each MODEL to its state; representatives.csv gives
  source DCDs and zero-based source frames. Unsampled nodes have no PDB.

Timing:
  --dt-ps describes the original SAVED FRAME spacing, not the integrator timestep.
  Effective spacing = --dt-ps * --stride; every lag must be a positive multiple.
  With --dt-ps 16 --stride 2, valid lags include 32, 64, 96, ... ps.
  Initial frames are discarded BEFORE striding each segment. Sampling restarts
  at original frame --discard-frames, then advances by --stride frames.
  Variable saved-frame spacing is not supported.

Quality assessment:
  quality_report.txt and quality_report.json explain diagnostic support, concerns
  and next steps; the main reports include them too. Console output gives a
  concise verdict. Rules are heuristic,
  not confidence levels; no automatic guarantee of convergence or kinetic validity.
  Lag screening requires >=3 usable lags at/above the chosen lag spanning >=2x.

Interpretation:
  The primary lag controls the saved MSM, stationary weights, FES and network.
  Extra lags test implied timescales on the SAME active state set.
  The original network shows observed counts; the pathway overlay shows routes
  from reactive flux decomposition. Neither arrow style represents a rate constant.
  Grey nodes are outside the selected active component; display filters affect
  arrows only. No pseudocounts are added to create artificial connections.
  Raw PCA density is sampling density, not equilibrium probability. The FES uses
  MSM stationary weights and is conditional on the active connected component
  and adequate within-state sampling. Disconnected components have unknown
  relative equilibrium weights.
  New PCA/clustering fits use training-segment frames only; the final MSM uses
  all segments. Reused models may already have seen held-out frames. CK uses a
  separate training MSM and held-out segment counts.
  Held-out segments may still be correlated through replica exchange; no
  confidence intervals or statistical convergence guarantee are supplied.

Practical notes:
  Load full-system PSF/DCDs with matching atom order. Coordinates must be whole;
  --unwrap reconstructs bonded fragments, not disconnected chains as a complex.
  Atom selections are fixed using the first input frame and reused throughout.
  PSF topology is parsed once. Without --unwrap, selected PCA/alignment atoms
  are extracted in small batches; libdcd still reads full frame records internally.
  With --unwrap, the full bonded system is processed frame by frame.
  Only supply completed immutable DCDs. Each rerun requires a new/empty --out.
  All selected coordinates are stored in RAM. --fit-max limits fitting work,
  not total memory; --stride and --selection reduce the coordinate data size.
  Diagnostics go to stderr and OUT/analysis.log. Results include PNG plots,
  CSV tables, NPZ model arrays, report.txt and report.json.
""")
    p.add_argument('psf', type=Path, help='Required PSF topology matching every DCD atom order.')

    g = p.add_argument_group('Input and output')
    g.add_argument('--jobname', metavar='NAME',
        help='DCD filename stem JOB. Default: PSF basename without .psf; override if the names differ.')
    g.add_argument('--output', type=Path, default=Path('output'), metavar='DIR',
        help='Input root containing numeric replica subdirectories. Default: %(default)s.')
    g.add_argument('--out', type=Path, default=Path('msm_results'), metavar='DIR',
        help='Results directory; must be new or empty. Default: %(default)s.')

    g = p.add_argument_group('Reusable state definition and structures')
    g.add_argument('--load-model', type=Path, default=None, metavar='NPZ',
        help='Reuse fitted_model.npz; fixes alignment reference, selections, unwrapping, PCA and cluster IDs. Overrides corresponding fitting options. Requires identical PSF. Default: fit anew. Every run writes fitted_model.npz; transition counts/MSM are always re-estimated from current inputs, not appended.')
    g.add_argument('--pdb-selection', default='protein', metavar='SELECTION',
        help='Atoms exported in representative node/path PDBs, independent of PCA selection. Use "all" for the full system. Default: %(default)s. Each path MODEL is a state representative, not a continuous trajectory.')

    g = p.add_argument_group('Timing and frame sampling')
    g.add_argument('--dt-ps', type=float, required=True, metavar='PS',
        help='Required constant original DCD frame interval in picoseconds; no default.')
    g.add_argument('--lag-ps', type=float, required=True, metavar='PS',
        help='Required primary MSM lag in picoseconds; no default. Must match the effective frame grid.')
    g.add_argument('--lags-ps', type=float, nargs='+', default=[], metavar='PS',
        help='Additional lags in picoseconds for implied-timescale diagnostics. The primary lag is always included. Default: none.')
    g.add_argument('--stride', type=int, default=1, metavar='N',
        help='Keep every Nth original frame within each segment (N >= 1). Effective interval becomes dt-ps * N. Default: %(default)s.')
    g.add_argument('--discard-frames', type=int, default=0, metavar='N',
        help='Discard the first N original frames of EACH segment before striding (N >= 0). Useful for a start-relaxation sensitivity test. Default: %(default)s.')

    g = p.add_argument_group('Coordinates, PCA and clustering')
    g.add_argument('--selection', default='protein and name CA', metavar='EXPR',
        help='MDAnalysis atom selection whose aligned Cartesian coordinates enter PCA; quote expressions containing spaces. Default: "%(default)s".')
    g.add_argument('--align-selection', default=None, metavar='EXPR',
        help='MDAnalysis selection for rotational/translational fitting (at least 3 noncollinear atoms). Default: same as --selection. Reference is the first retained input frame, or the saved reference with --load-model.')
    g.add_argument('--unwrap', action='store_true',
        help='Reconstruct whole bonded fragments per frame using PSF bonds and valid DCD cell dimensions. Does not assemble disconnected chains. Default: off.')
    g.add_argument('--pcs', type=int, default=10, metavar='N',
        help='Number of unwhitened coordinate PCs used for clustering (N >= 2). Capped by feature count and training frames minus one. Plots show PC1/PC2 only. Default: %(default)s.')
    g.add_argument('--clusters', type=int, default=50, metavar='N',
        help='Requested MiniBatchKMeans microstates in retained PCA space (N >= 2); capped by available fitting frames. These are not automatically metastable macrostates. Default: %(default)s.')
    g.add_argument('--fit-max', type=int, default=50000, metavar='N',
        help='Maximum randomly sampled training frames for PCA and clustering fits (N >= 3). All retained frames are projected, assigned and counted afterward. Does not cap total memory. Default: %(default)s.')

    g = p.add_argument_group('MSM estimation and validation')
    g.add_argument('--reversible', action='store_true',
        help='Use deeptime reversible maximum-likelihood estimation with detailed balance. Default: off (nonreversible row-normalized counts). Neither option automatically removes reseeding bias.')
    g.add_argument('--test-fraction', type=float, default=.2, metavar='F',
        help='Fraction of whole segments held out from PCA/clustering fitting and CK training (0 <= F < 0.5); final MSM uses all segments. Zero disables holdout. With fewer than 3 segments holdout is disabled. Default: %(default)s.')
    g.add_argument('--seed', type=int, default=42, metavar='N',
        help='Nonnegative random seed for segment splitting, fitting-frame subsampling, PCA and clustering. Default: %(default)s.')

    g = p.add_argument_group('Transition pathways (optional)')
    g.add_argument('--start-states',type=int,nargs='+',default=None,metavar='ID',
        help='Source state ID(s) A from states.csv or network node labels, NOT injection milestone IDs. Requires --end-states. Default: none (pathway analysis disabled).')
    g.add_argument('--end-states',type=int,nargs='+',default=None,metavar='ID',
        help='Target state ID(s) B, disjoint from A and in the same active MSM component. Requires --start-states. Default: none.')
    g.add_argument('--max-paths',type=int,default=5,metavar='N',
        help='Maximum dominant TPT flux-decomposition routes to compute and overlay (N >= 1). All route nodes are retained. Default: %(default)s.')
    g.add_argument('--path-fraction',type=float,default=.9,metavar='F',
        help='Target fraction of total A-to-B reactive flux to explain (0 < F <= 1). Stop when reached or --max-paths is exhausted; actual coverage is reported. Default: %(default)s.')

    g = p.add_argument_group('Landscapes and network display')
    g.add_argument('--temperature', type=float, default=300, metavar='K',
        help='Baseline temperature in kelvin for F = -RT ln(p), reported in kJ/mol. Does not change transition counts or kinetic timescales. Default: %(default)s.')
    g.add_argument('--bins', type=int, default=40, metavar='N',
        help='Histogram bins per axis for PC1/PC2 density and free-energy plots (N >= 1); N squared cells. Does not change MSM clustering. Default: %(default)s.')
    g.add_argument('--min-edge-count', type=int, default=2, metavar='N',
        help='Minimum observed lagged-pair count for a directed off-diagonal network arrow (N >= 1). Display only; all counts remain in the MSM. Default: %(default)s.')
    g.add_argument('--max-edges', type=int, default=200, metavar='N',
        help='Maximum displayed network arrows, keeping those with the largest counts (N >= 1). Display only; does not alter the MSM. Default: %(default)s.')
    return p


def main():
    p = build_parser()
    args=p.parse_args()
    if not args.psf.is_file():p.error('PSF does not exist')
    if args.load_model is not None:
        try:
            load_fitted_model(args)
        except (OSError, ValueError, KeyError, TypeError) as exc:
            p.error(f'Cannot load fitted model: {exc}')
    if not np.isfinite(args.temperature) or args.temperature<=0:p.error('Temperature must be finite and positive')
    if min(args.stride,args.pcs,args.clusters,args.fit_max,args.bins,args.min_edge_count,args.max_edges)<1:p.error('Positive integer options required')
    if args.pcs<2 or args.clusters<2 or args.discard_frames<0:p.error('Need >=2 PCs/clusters and nonnegative discard')
    if not 0<=args.test_fraction<.5:p.error('--test-fraction must be >=0 and <0.5')
    if args.fit_max < 3:p.error('--fit-max must be at least 3')
    if args.seed < 0:p.error('--seed must be nonnegative')
    if not args.output.is_dir():p.error(f'Input directory does not exist: {args.output}')

    if (args.start_states is None) != (args.end_states is None):
        p.error('--start-states and --end-states must be supplied together')
    if args.max_paths<1:p.error('--max-paths must be positive')
    if not np.isfinite(args.path_fraction) or not 0<args.path_fraction<=1:
        p.error('--path-fraction must be >0 and <=1')
    if args.start_states is not None:
        if set(args.start_states)&set(args.end_states):p.error('Start/end states must be disjoint')
        if any(i<0 or i>=args.clusters for i in args.start_states+args.end_states):
            p.error('Endpoint IDs must be >=0 and < --clusters; use cluster IDs, not milestone IDs')
        try:
            from deeptime.markov import reactive_flux
        except ImportError:
            p.error('Pathway analysis requires deeptime; install it before running')

    # Validate every requested lag before DCD I/O or creating output files.
    effective_dt = args.dt_ps * args.stride
    lag_errors = []
    for option, value in [('--lag-ps', args.lag_ps)] + [('--lags-ps', v) for v in args.lags_ps]:
        try:
            lag_frames(value, effective_dt)
        except ValueError as exc:
            lag_errors.append(f'{option}: {exc}')
    if lag_errors:
        p.error('Invalid lag settings:\n  ' + '\n  '.join(lag_errors))
    if args.out.exists() and not args.out.is_dir():p.error('--out must be a directory')
    if args.out.exists() and any(args.out.iterdir()):p.error('Use a new/empty --out directory to avoid stale results')
    args.out.mkdir(parents=True,exist_ok=True)
    logger = logging.getLogger('shooting_msm')
    logger.setLevel(logging.INFO)
    logger.propagate = False
    formatter = logging.Formatter('%(asctime)s %(message)s', datefmt='%H:%M:%S')
    for handler in (logging.StreamHandler(sys.stderr), logging.FileHandler(args.out/'analysis.log')):
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    if args.load_model is not None:
        progress(f'Loaded state definition {args.load_model}: stored selections, unwrap, PCA and clusters override fitting options')
    progress(f'Starting analysis; PSF: {args.psf}; results: {args.out}')
    progress(f'Preflight passed: effective frame spacing {effective_dt:g} ps; all requested lags compatible')
    arrays,metadata=load_segments(args,args.out)
    analyze(args,arrays,metadata,args.out)
    progress('Analysis complete')
    print(f'Results: {args.out.resolve()}')

if __name__=='__main__':
    try:main()
    except (ValueError,OSError,ImportError) as e:
        sys.exit(f'ERROR: {e}')
