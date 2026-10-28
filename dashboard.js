/* ------------------------------------------------------------------ utils */
const $ = s => document.querySelector(s);
let CVD = false;
const el = (t, a = {}, ...kids) => { const n = document.createElement(t);
  for (const [k, v] of Object.entries(a)) { if (k === 'class') n.className = v;
    else if (k === 'html') n.innerHTML = v; else if (k.startsWith('on')) n[k] = v;
    else n.setAttribute(k, v); }
  kids.flat().forEach(c => { if (c === null || c === undefined) return;
    n.append(c.nodeType ? c : document.createTextNode(c)); });
  return n; };
const fmt = (x, d = 2) => (x === null || x === undefined || Number.isNaN(x)) ? '–'
  : (typeof x === 'number' ? (Number.isInteger(x) ? x : x.toFixed(d)) : x);
/* Display name for a species. Sample ids carry the suffix of the file they
   came from (Takifugu_rubripes_pep), which is provenance, not a name, and it
   repeats on every axis label. COMMON[] is filled from NCBI's taxdump when the
   metadata has a common_name column. */
const SUFFIX = /_(pep|protein|prot|aa|cds|longest|primary)$/i;
const nice = s => String(s).replace(SUFFIX, '').replace(/_/g, ' ');
let COMMON = {};
/* When a sample sheet is present it holds the real species name, which beats
   anything a regex can recover from a file name such as
   X_alpha_020202_okay. The method is appended only where several proteomes of
   the same species would otherwise carry identical labels. */
const SAMP = D.samples || {};
const SPCOUNT = {};
Object.values(SAMP).forEach(m => { if (m && m.species)
  SPCOUNT[m.species] = (SPCOUNT[m.species] || 0) + 1; });
function sciName(s) {
  const m = SAMP[s];
  return (m && m.species) ? nice(m.species) : nice(s);
}
function qualifier(s) {
  const m = SAMP[s];
  return (m && m.species && SPCOUNT[m.species] > 1 && m.method) ? m.method : '';
}
function commonName(s) { const c = COMMON[s]; return c ? String(c) : ''; }
/* "Scientific name (Common name)", or just the scientific name when unknown.
   `max` trims the common name rather than the scientific one, since the latter
   is the identifier. */
function fullName(s, max) {
  const sci = sciName(s), c = commonName(s), q = qualifier(s);
  const tail = q ? ` · ${q}` : '';
  if (!c) return sci + tail;
  let out = `${sci} (${c})${tail}`;
  if (max && out.length > max) {
    const room = max - sci.length - tail.length - 3;
    out = room > 4 ? `${sci} (${c.slice(0, room - 1)}…)${tail}` : sci + tail;
  }
  return out;
}
/* Right-aligned species label on a canvas: the binomial in italic, the
   vernacular upright in brackets after it, which is the convention and is also
   what stops a reader taking "zebrafish" for part of the name. Trimming eats
   into the common name only; the binomial is the identifier. `budget` is in
   pixels, not characters, because the two fonts differ in width. */
function drawSpeciesLabel(g, s, x, y, budget, size = 11) {
  const sci = sciName(s), c = commonName(s);
  const q = qualifier(s);
  const IT = `italic ${size}px sans-serif`, RO = `${size}px sans-serif`;
  g.textAlign = 'right';
  if (!c && !q) { g.font = IT; g.fillText(sci, x, y); return; }
  g.font = IT; const wSci = g.measureText(sci).width;
  g.font = RO; let txt = (c ? ` (${c})` : '') + (q ? ` \u00b7 ${q}` : '');
  if (c && wSci + g.measureText(txt).width > budget) {
    let cc = c;
    const tail = q ? ` \u00b7 ${q}` : '';
    while (cc.length > 1 &&
           wSci + g.measureText(` (${cc}…)${tail}`).width > budget) cc = cc.slice(0, -1);
    txt = ` (${cc}…)${tail}`;
  }
  const wC = g.measureText(txt).width;
  g.fillText(txt, x, y);
  g.font = IT; g.fillText(sci, x - wC, y);
}
/* The same label as two SVG tspans. */
function svgSpeciesLabel(s, attrs, size = 11) {
  const t = svgEl('text', {...attrs, 'font-size': size});
  const sci = sciName(s), c = commonName(s);
  const q = qualifier(s);
  const a = svgEl('tspan', {'font-style': 'italic'}); a.textContent = sci; t.append(a);
  if (c || q) { const b = svgEl('tspan', {});
    b.textContent = (c ? ` (${c})` : '') + (q ? ` \u00b7 ${q}` : ''); t.append(b); }
  return t;
}
/* Text is measured on a canvas to size gutters, but SVG labels render in the
   page font. Measuring in canvas's default sans-serif under-sized gutters by
   several pixels and clipped the longest label, so measure in the real font. */
const UIFONT = (() => { try { return getComputedStyle(document.body).fontFamily || 'sans-serif'; }
                        catch (e) { return 'sans-serif'; } })();
const measurer = (px, weight = '') => { const g = el('canvas').getContext('2d');
  g.font = `${weight} ${px}px ${UIFONT}`.trim(); return g; };
/* Secondary detail (summary tables, long lists) sits behind a one-line toggle,
   closed, so the figure is what a reader meets first. */
function collapsible(title, ...kids) {
  return el('details', {class: 'more'}, el('summary', {}, title), ...kids);
}
/* Draw one label in k when a categorical axis has more rows than it has room
   for. Returns 1 below the limit, so small runs are unaffected. */
const labelStep = (n, max = 60) => Math.max(1, Math.ceil(n / max));
/* A <select> reports the first option's value in a browser, but code that
   assumes it can crash on an empty value. Always fall back to a known-good
   member of the list. */
const pick = (sel, list) => (sel && list.includes(sel.value)) ? sel.value : list[0];
const dl = (name, text) => { const b = new Blob([text], {type: 'text/tab-separated-values'});
  const u = URL.createObjectURL(b); const a = el('a', {href: u, download: name}); a.click();
  URL.revokeObjectURL(u); };

/* Taxonomic filtering works by rebuilding a VIEW of the payload rather than
   teaching every tab about a filter: each tab already reads D.species[i] and
   D.counts[fam][i] at the current index, so subsetting those arrays and the
   SP/SPI lookups propagates everywhere, including figures and exports. D0 is
   the untouched payload the view is derived from. */
const D0 = {species: D.species, counts: D.counts, per10k: D.per10k,
            rescued: D.rescued || {}};
const TAX_COLS = D.tax_cols || [];
let FILTER = {};                  // {rank: value}, empty means everything
let SP, SPI, TRAITS, GROUP_BY = D.trait_name || null;

function speciesPasses(s) {
  for (const k in FILTER) {
    const want = FILTER[k]; if (!want) continue;
    const have = (s.tax && s.tax[k]) || (k === (D.trait_name || '') ? s.trait : null);
    if (String(have ?? '—') !== want) return false;
  }
  return true;
}
function groupValue(s) {
  if (!GROUP_BY) return s.trait ?? null;
  if (s.tax && GROUP_BY in s.tax) return s.tax[GROUP_BY];
  return GROUP_BY === D.trait_name ? (s.trait ?? null) : null;
}
function rebuildView() {
  const keep = [];
  D0.species.forEach((s, i) => { if (speciesPasses(s)) keep.push(i); });
  D.species = keep.map(i => ({...D0.species[i], trait: groupValue(D0.species[i])}));
  const sub = m => { if (!m) return m; const o = {};
    for (const f in m) o[f] = keep.map(i => m[f][i]); return o; };
  D.counts = sub(D0.counts); D.per10k = sub(D0.per10k); D.rescued = sub(D0.rescued);
  D.trait_name = GROUP_BY || D.trait_name;
  SP = D.species.map(s => s.name);
  SPI = Object.fromEntries(SP.map((s, i) => [s, i]));
  TRAITS = [...new Set(D.species.map(s => s.trait).filter(x => x))].sort();
  return keep.length;
}
COMMON = Object.fromEntries(D0.species
  .map(s => [s.name, (s.tax && (s.tax.common_name || s.tax.common)) || ''])
  .filter(([, c]) => c));
rebuildView();
const FAM = D.families.map(f => f.family);
const FMETA = Object.fromEntries(D.families.map(f => [f.family, f]));
const TCOL = ['#bc4749','#e9a13b','#2a9d8f','#1b4965','#7b2cbf','#6a994e','#9c6644'];
let traitColour = t => t == null ? '#c9ced4' : TCOL[TRAITS.indexOf(t) % TCOL.length];
const CATS = [...new Set(D.families.map(f => f.category))];
const CCOL = ['#1b4965','#2a9d8f','#8ecae6','#bc4749','#e9a13b','#7b2cbf','#6a994e','#9c6644'];
const catColour = c => CCOL[CATS.indexOf(c) % CCOL.length];

const tip = $('#tip');
function showTip(e, html) { tip.innerHTML = html; tip.style.opacity = 1;
  const pad = 14; let x = e.pageX + pad, y = e.pageY + pad;
  if (x + tip.offsetWidth > window.innerWidth + window.scrollX) x = e.pageX - tip.offsetWidth - pad;
  tip.style.left = x + 'px'; tip.style.top = y + 'px'; }
function hideTip() { tip.style.opacity = 0; }

/* matrix access: 'raw' | 'per10k' */
function matrix(kind) { return (kind === 'per10k' && D.per10k) ? D.per10k : D.counts; }
function col(kind, fam) { const m = matrix(kind); return m[fam] || SP.map(() => 0); }
function stats(v) { const n = v.length, mean = v.reduce((a, b) => a + b, 0) / n;
  const sd = Math.sqrt(v.reduce((a, b) => a + (b - mean) ** 2, 0) / (n - 1 || 1));
  const s = [...v].sort((a, b) => a - b);
  const q = p => { const i = (s.length - 1) * p, lo = Math.floor(i), hi = Math.ceil(i);
    return s[lo] + (s[hi] - s[lo]) * (i - lo); };
  return {mean, sd, min: s[0], max: s[s.length - 1], q1: q(.25), med: q(.5), q3: q(.75)}; }
function zcol(kind, fam) { const v = col(kind, fam), st = stats(v);
  return st.sd === 0 ? v.map(() => 0) : v.map(x => (x - st.mean) / st.sd); }

/* families that carry no variance cannot be z-scored; say so rather than hide them */
const INVARIANT = {};
FAM.forEach(f => { const raw = D.counts[f]; if (!raw) return;
  const sum = raw.reduce((a, b) => a + b, 0);
  const allSame = raw.every(x => x === raw[0]);
  if (sum === 0) INVARIANT[f] = 'never detected in any species';
  else if (allSame) INVARIANT[f] = 'constant at ' + raw[0] + ' in every species'; });
let VARFAM = FAM.filter(f => !INVARIANT[f]);

/* diverging + sequential colour ramps, no dependencies */
function rdbu(t) { t = Math.max(0, Math.min(1, t));
  const a = [[5,48,97],[33,102,172],[67,147,195],[146,197,222],[209,229,240],
             [247,247,247],[253,219,199],[244,165,130],[214,96,77],[178,24,43],[103,0,31]];
  const i = Math.min(a.length - 2, Math.floor(t * (a.length - 1))), f = t * (a.length - 1) - i;
  const c = a[i].map((v, k) => Math.round(v + (a[i + 1][k] - v) * f));
  return `rgb(${c[0]},${c[1]},${c[2]})`; }
function viridis(t) { t = Math.max(0, Math.min(1, t));
  const a = [[68,1,84],[72,40,120],[62,74,137],[49,104,142],[38,130,142],
             [31,158,137],[53,183,121],[110,206,88],[181,222,43],[253,231,37]];
  const i = Math.min(a.length - 2, Math.floor(t * (a.length - 1))), f = t * (a.length - 1) - i;
  const c = a[i].map((v, k) => Math.round(v + (a[i + 1][k] - v) * f));
  return `rgb(${c[0]},${c[1]},${c[2]})`; }

/* sortable table with search and TSV export */
function makeTable(host, cols, rows, opts = {}) {
  /* Drop columns that are empty in every row. A column of dashes is not
     information, it is a reminder that a step was not run, and the Methods tab
     is where that belongs. `keepAll` opts out for tables where the reader is
     checking that something is absent. */
  if (!opts.keepAll && rows.length) {
    const empty = v => v === null || v === undefined || v === '' ||
      v === '\u2013' || v === '\u2014' || (typeof v === 'number' && Number.isNaN(v));
    const keep = cols.map((_, i) => !rows.every(r => empty(r[i])));
    if (keep.some(k => !k) && keep.some(k => k)) {
      const at = [];                       // old column index -> new, or -1
      let j = 0; keep.forEach(k => at.push(k ? j++ : -1));
      const remap = a => (a || []).map(i => at[i]).filter(i => i >= 0);
      cols = cols.filter((_, i) => keep[i]);
      rows = rows.map(r => r.filter((_, i) => keep[i]));
      const sb = at[opts.sortBy ?? 0];
      opts = {...opts, sortBy: sb >= 0 ? sb : 0, desc: sb >= 0 ? opts.desc : true,
              mono: remap(opts.mono), ital: remap(opts.ital), sci: remap(opts.sci)};
    }
  }
  let sortI = opts.sortBy ?? 0, dir = opts.desc ? -1 : 1, filter = '';
  const search = el('input', {type: 'search', placeholder: 'filter…',
    oninput: e => { filter = e.target.value.toLowerCase(); draw(); }});
  const dlBtn = el('button', {class: 'act', onclick: () => dl((opts.name || 'table') + '.tsv',
    [cols.join('\t'), ...view().map(r => r.join('\t'))].join('\n'))}, 'download TSV');
  const count = el('span', {class: 'hint'});
  const wrap = el('div', {class: 'scroll'});
  host.append(el('div', {class: 'row'}, el('div', {class: 'ctl'},
    el('label', {}, 'search'), search), dlBtn, count), wrap);
  const view = () => { let r = rows;
    if (filter) r = r.filter(x => x.some(v => String(v).toLowerCase().includes(filter)));
    return [...r].sort((a, b) => { const x = a[sortI], y = b[sortI];
      if (typeof x === 'number' && typeof y === 'number') return (x - y) * dir;
      return String(x).localeCompare(String(y)) * dir; }); };
  function draw() { const r = view(); count.textContent = r.length + ' rows';
    const cap = opts.cap || 800;
    const t = el('table', {},
      el('thead', {}, el('tr', {}, cols.map((c, i) => el('th', {
        onclick: () => { if (sortI === i) dir = -dir; else { sortI = i; dir = 1; } draw(); }},
        c + (sortI === i ? (dir > 0 ? ' ▲' : ' ▼') : ''))))),
      el('tbody', {}, r.slice(0, cap).map(row => el('tr', {}, row.map((v, i) =>
        el('td', {class: (opts.mono && opts.mono.includes(i) ? 'mono ' : '') +
          (opts.ital && opts.ital.includes(i) ? 'sci' : '')},
          (opts.sci && opts.sci.includes(i) && typeof v === 'number' && v !== 0)
            ? v.toPrecision(3) : fmt(v)))))));
    wrap.innerHTML = ''; wrap.append(t);
    if (r.length > cap) wrap.append(el('div', {class: 'note'},
      `showing the first ${cap} of ${r.length} rows. Narrow with the filter box, or download the full TSV.`)); }
  draw();
}

/* ------------------------------------------------------------- OVERVIEW */
/* Horizontal stacked bars: one row per group, one segment per key. Used for
   the composition view, which is the question the Overview exists to answer:
   what is each group's defensome made of. */
function stackedBars(rows, keys, colourOf, W, opts = {}) {
  const mg = measurer(12);
  const labs = rows.map(r => `${opts.raw ? r.name : nice(r.name)} (${r.n.toLocaleString()})`);
  const ML = Math.round(Math.min(W * 0.34, Math.max(90, Math.max(...labs.map(t =>
    mg.measureText(t).width), 10) + 14)));
  const MT = 8, rh = opts.rh || 28;
  // the right margin is sized to the end labels, so "100% complete" is not cut
  const ends = rows.map(r => opts.share ? (opts.endLabel ? opts.endLabel(r) : '100%')
                                        : fmt(keys.reduce((a, k) => a + (r.v[k] || 0), 0), 1));
  mg.font = `11px ${UIFONT}`;
  const MR = Math.ceil(Math.max(40, ...ends.map(t => mg.measureText(String(t)).width)) + 14);
  const H = MT + rows.length * rh + 26;
  const sv = svgEl('svg', {width: W, height: H, viewBox: `0 0 ${W} ${H}`});
  const span = W - ML - MR;
  const mx = Math.max(...rows.map(r => keys.reduce((a, k) => a + (r.v[k] || 0), 0)), 1e-9);
  rows.forEach((r, i) => {
    const y = MT + i * rh, tot = keys.reduce((a, k) => a + (r.v[k] || 0), 0);
    const scale = opts.share ? (tot ? span / tot : 0) : span / mx;
    let x = ML;
    keys.forEach(k => {
      const v = r.v[k] || 0; if (!v) return;
      const w = v * scale;
      const rect = svgEl('rect', {x, y: y + 4, width: Math.max(0.5, w), height: rh - 12,
        fill: colourOf(k)});
      rect.addEventListener('mousemove', e => showTip(e,
        `<b>${opts.raw ? r.name : nice(r.name)}</b> · ${k}<br>${fmt(v, 2)} ${opts.unit || ''}` +
        `<br>${(tot ? v / tot * 100 : 0).toFixed(1)}% of the defensome` +
        `<br>${r.n.toLocaleString()} ${opts.nUnit || 'species'}`));
      rect.addEventListener('mouseleave', hideTip);
      sv.append(rect); x += w;
    });
    const t = svgEl('text', {x: ML - 8, y: y + rh / 2 + 1, 'text-anchor': 'end',
      'font-size': 12, 'dominant-baseline': 'middle'});
    t.textContent = labs[i]; sv.append(t);
    const t2 = svgEl('text', {x: x + 6, y: y + rh / 2 + 1, 'font-size': 11,
      fill: '#6b7280', 'dominant-baseline': 'middle'});
    t2.textContent = opts.share ? (opts.endLabel ? opts.endLabel(r) : '100%') : fmt(tot, 1); sv.append(t2);
  });
  return sv;
}

function tabOverview(host) {
  /* Home: what this is, where to go, and one picture. Nothing about how the
     calls were made; that is on the About page. */
  const grp = D.trait_name || (RANK_COLS.includes('class') ? 'class' : TAX_COLS[0]);
  const nPhyla = RANK_COLS.includes('phylum')
    ? new Set(D0.species.map(s => (s.tax || {}).phylum).filter(x => x)).size : 0;
  const hero = el('div', {class: 'hero'});
  hero.append(el('h1', {}, 'Chemical defensome'),
    el('div', {class: 'stats'},
      el('a', {href: '#/species'}, el('b', {}, D0.species.length.toLocaleString()), ' species'),
      el('a', {href: '#/families'}, el('b', {}, FAM.length), ' gene families'),
      ...(nPhyla ? [el('span', {}, el('b', {}, nPhyla), ' phyla')] : [])),
    el('div', {class: 'row', style: 'margin:14px 0 0;gap:10px'},
      el('a', {class: 'btn', href: '#/species'}, 'Browse species'),
      el('a', {class: 'btn', href: '#/families'}, 'Browse gene families'),
      ...(TRAIT_COLS.length ? [el('a', {class: 'btn ghost', href: '#/compare'}, `Compare by ${TRAIT_COLS[0]}`)] : [])));
  host.append(hero);

  /* browse by group: one card per value of the grouping, largest first */
  if (grp) {
    const n = {};
    D.species.forEach(sp => { const v = sp.trait ?? (sp.tax || {})[grp]; if (v) n[v] = (n[v] || 0) + 1; });
    const vals = Object.entries(n).sort((a, b) => b[1] - a[1]);
    if (vals.length > 1) {
      const p = el('div', {class: 'panel'});
      p.append(el('h2', {}, `By ${grp}`));
      const cards = el('div', {class: 'cards'});
      vals.slice(0, 16).forEach(([v, k]) => cards.append(el('a', {href: '#/species', class: 'card',
        onclick: e => { e.preventDefault(); showTaxon(grp, v); }},
        el('i', {class: 'sw', style: 'background:' + traitColour(v)}),
        el('span', {class: 'nm'}, nice(v)), el('span', {class: 'n'}, k.toLocaleString()))));
      if (vals.length > 16) cards.append(el('a', {href: '#/species', class: 'card more'},
        `+ ${vals.length - 16} more`));
      p.append(cards); host.append(p);
    }
  }

  /* ---- composition by group ---- */
  const S = {kind: D.per10k ? 'per10k' : 'raw', share: true, level: 'category'};
  const p = el('div', {class: 'panel'});
  p.append(el('h2', {}, 'What each group\u2019s defensome is made of'));
  const seg = (opts, cur, set) => el('div', {class: 'seg'}, ...opts.map(([v, l]) =>
    el('button', {class: v === cur() ? 'on' : '', onclick: e => { set(v);
      [...e.target.parentNode.children].forEach(b => b.classList.toggle('on', b === e.target)); draw(); }}, l)));
  p.append(el('div', {class: 'row', style: 'align-items:center'},
    seg([['share', 'share'], ['abs', 'copies per species']], () => S.share ? 'share' : 'abs', v => { S.share = v === 'share'; }),
    seg([['category', 'by category'], ['family', 'by family']], () => S.level, v => { S.level = v; })));
  const body = el('div'); const leg = el('div', {class: 'legend'});
  p.append(body, leg); host.append(p);

  function draw() {
    body.innerHTML = ''; leg.innerHTML = '';
    const keys = S.level === 'category' ? CATS : FAM;
    const colourOf = S.level === 'category' ? catColour
      : (f => catColour(FMETA[f] ? FMETA[f].category : CATS[0]));
    const famsOf = k => S.level === 'category'
      ? FAM.filter(f => FMETA[f].category === k) : [k];
    const byGroup = {};
    D.species.forEach((sp, i) => {
      const g = sp.trait ?? '\u2014';
      (byGroup[g] = byGroup[g] || []).push(i); });
    const rows = Object.entries(byGroup)
      .map(([name, idx]) => { const v = {};
        keys.forEach(k => { const fs = famsOf(k);
          v[k] = idx.reduce((a, i) =>
            a + fs.reduce((b, f) => b + col(S.kind, f)[i], 0), 0) / idx.length; });
        return {name, n: idx.length, v, tot: keys.reduce((a, k) => a + v[k], 0)}; })
      .sort((a, b) => b.n - a.n).slice(0, 30);
    if (!rows.length) { body.append(el('div', {class: 'empty'}, 'No species in view.')); return; }
    const W = Math.max(520, Math.min(1180, (p.clientWidth || 1100) - 36));
    const unit = S.kind === 'per10k' ? `per 10k ${D.size_unit || 'proteins'}` : 'copies';
    body.append(figure(stackedBars(rows, keys, colourOf, W, {share: S.share, unit}),
      'defensome_composition'));
    if (S.level === 'category')
      keys.forEach(k => leg.append(el('span', {},
        el('i', {class: 'sw', style: 'background:' + colourOf(k)}), k)));
    tabOverview._rows = rows;
  }
  draw();
}

/* Quality control, for the About page. */
function qcPanels(host) {
  const nInv = Object.keys(INVARIANT).length;
  const flagged = D0.species.filter(x => x.flag && x.flag !== 'ok');
  if (flagged.length) {
    const f = el('div', {class: 'panel'});
    f.append(el('h2', {}, `Flagged species (${flagged.length})`), el('p', {class: 'hint'},
      'Each species against the median of the control families. INFLATED: one control at 3x the median or two at 2x; SPARSE: two controls missing.'));
    dataGrid(f, [
      {label: 'species', get: r => sciName(r.name), cell: r => spLink(r.name)},
      {label: 'flag', get: r => r.flag}, {label: 'worst family', get: r => r.worst_family || ''},
      {label: 'max ratio', num: true, get: r => r.max_ratio}], flagged,
      {sortBy: 3, desc: true, pageSize: 20, name: 'qc_flagged', unit: 'species', units: 'species'});
    host.append(f);
  }
  if (D.qc_families) {
    const q = el('div', {class: 'panel'});
    q.append(el('h2', {}, 'Are the control families behaving?'),
      el('p', {class: 'hint'}, 'Observed median across species against the literature expectation for insect genomes. A high ratio means the map rule is catching more than the family name implies, which is a map problem rather than a genome problem.'));
    const qb = el('div');
    makeTable(qb, ['family', 'expected', 'observed median', 'min', 'max', 'ratio', 'verdict'],
      D.qc_families.map(r => [r.family, r.expected, r.observed_median, r.observed_min,
        r.observed_max, r.median_vs_expected, r.verdict]), {sortBy: 5, desc: true, name: 'qc_families'});
    q.append(collapsible(`control families (${D.qc_families.length})`, qb));
    host.append(q);
  }

  if (nInv) {
    const z = el('div', {class: 'panel'});
    z.append(el('h2', {}, nInv === 1 ? 'One family cannot be z-scored'
                                     : `${nInv} families cannot be z-scored`),
      el('p', {class: 'hint'}, 'A family with no variance has no z-score, so it never appears in a heatmap. Both cases below are findings, not gaps, so they are listed rather than dropped.'));
    makeTable(z, ['family', 'Pfam', 'why', 'what to do'],
      Object.entries(INVARIANT).map(([f, why]) => [f, FMETA[f] ? FMETA[f].pfam : '', why,
        why.startsWith('never')
          ? (D.qc_zero && (D.qc_zero.find(x => x.family === f) || {}).verdict
              ? 'verdict: ' + D.qc_zero.find(x => x.family === f).verdict.replace(/_/g, ' ').toLowerCase()
              : 'run qc: it writes qc_zero_families.tsv with a verdict for every empty family')
          : 'single-copy control: the pipeline is calibrated']),
      {sortBy: 0, mono: [1], name: 'invariant_families', keepAll: true});
    host.append(z);
  }
}
const kpiEmpty = v => v === null || v === undefined || v === '' ||
  v === '\u2013' || v === '\u2014' || (typeof v === 'number' && Number.isNaN(v));
const kpi = (v, l) => kpiEmpty(v) ? null : el('div', {class: 'kpi'}, el('div', {class: 'v'}, String(v)),
  el('div', {class: 'l'}, l));

/* -------------------------------------------------------------- HEATMAP */
function tabHeatmap(host) {
  const S = {kind: D.per10k ? 'per10k' : 'raw', scale: 'z', tiers: new Set(['CORE']),
    cats: new Set(CATS), sort: 'trait', sortFam: null};
  const p = el('div', {class: 'panel'});
  p.append(el('h2', {}, 'Species by family'), el('p', {class: 'hint'},
    'Click a column header to sort species by that family. Click a cell for the family detail. Raw counts track proteome size; per-10k asks what share of the gene budget goes to defence.'));
  const ctrls = el('div', {class: 'row'});
  const mk = (label, node) => el('div', {class: 'ctl'}, el('label', {}, label), node);
  const kindSel = el('select', {onchange: e => { S.kind = e.target.value; draw(); }},
    ...(D.per10k ? [el('option', {value: 'per10k'}, `per 10k ${D.size_unit || 'proteins'}`)] : []),
    el('option', {value: 'raw'}, 'raw copies'));
  const scaleSel = el('select', {onchange: e => { S.scale = e.target.value; draw(); }},
    el('option', {value: 'z'}, 'z-score per family'),
    el('option', {value: 'abs'}, 'absolute value'),
    el('option', {value: 'log'}, 'log10(x+1)'));
  // Sorting is offered on every metadata column, not only the grouping one, so
  // you can group by class and still read the rows in habitat order.
  const sortSel = el('select', {onchange: e => { S.sort = e.target.value; S.sortFam = null; draw(); }},
    el('option', {value: 'trait'}, (D.trait_name || 'trait') + ' (group)'),
    el('option', {value: 'name'}, 'species name'),
    el('option', {value: 'total'}, 'total defensome'),
    el('option', {value: 'size'}, 'proteome size'),
    ...TAX_COLS.filter(c => c !== D.trait_name)
      .map(c => el('option', {value: 'tax:' + c}, c)));
  ctrls.append(mk('values', kindSel), mk('scale', scaleSel), mk('sort species by', sortSel));
  const tierChips = el('div', {class: 'chips'});
  ['CORE', 'BROAD'].forEach(t => tierChips.append(el('span', {
    class: 'chip' + (S.tiers.has(t) ? ' on' : ''),
    onclick: e => { S.tiers.has(t) ? S.tiers.delete(t) : S.tiers.add(t);
      e.target.classList.toggle('on'); draw(); }}, t)));
  const catChips = el('div', {class: 'chips'});
  CATS.forEach(c => catChips.append(el('span', {class: 'chip on',
    onclick: e => { S.cats.has(c) ? S.cats.delete(c) : S.cats.add(c);
      e.target.classList.toggle('on'); draw(); }}, c)));
  ctrls.append(mk('tier', tierChips), mk('category', catChips));
  p.append(ctrls);
  /* Two canvases: the column headers stay put while the rows scroll, which is
     what makes a 1,575-row matrix usable. They share one horizontal scroll
     position so the headers never drift off their columns. */
  const hdr = el('canvas', {style: 'display:block'});
  const hbox = el('div', {style: 'overflow:hidden'}, hdr);
  const cv = el('canvas');
  const body = el('div', {style: 'overflow:auto;max-height:620px'}, cv);
  // One frame around headers and rows, sized to the matrix rather than to the
  // panel, so the border hugs the data instead of enclosing empty space.
  const frame = el('div', {class: 'matrix'}, hbox, body);
  body.onscroll = () => { hbox.scrollLeft = body.scrollLeft; };
  // The export has to be the whole matrix, headers included, not the part of
  // the body that happens to be scrolled into view.
  const saveWhole = () => { const out = el('canvas');
    out.width = hdr.width; out.height = hdr.height + cv.height;
    const o = out.getContext('2d');
    o.fillStyle = '#fff'; o.fillRect(0, 0, out.width, out.height);
    o.drawImage(hdr, 0, 0); o.drawImage(cv, 0, hdr.height);
    savePNG(out, 'defensome_heatmap'); };
  p.append(frame, el('div', {class: 'row', style: 'margin:6px 0 0'},
    el('button', {class: 'act', onclick: saveWhole}, 'PNG')));
  const cbar = el('div', {style: 'margin-top:10px'}); p.append(cbar);
  const leg = el('div', {class: 'legend'}); p.append(leg);
  host.append(p);

  function draw() {
    const fams = FAM.filter(f => S.tiers.has(FMETA[f].tier) && S.cats.has(FMETA[f].category)
      && !(S.scale === 'z' && INVARIANT[f]));
    let order = SP.map((s, i) => i);
    const coreV = fams.length ? fams : VARFAM;
    const tot = i => coreV.reduce((a, f) => a + col(S.kind, f)[i], 0);
    const taxSort = S.sort.startsWith('tax:') ? S.sort.slice(4) : null;
    // The value whose contiguous runs get a band label: whatever the rows are
    // ordered by, if that is categorical.
    const bandOf = taxSort ? (i => (D.species[i].tax || {})[taxSort])
                           : (i => D.species[i].trait);
    if (S.sortFam) { const v = col(S.kind, S.sortFam); order.sort((a, b) => v[b] - v[a]); }
    else if (S.sort === 'name') order.sort((a, b) => sciName(SP[a]).localeCompare(sciName(SP[b])));
    else if (S.sort === 'total') order.sort((a, b) => tot(b) - tot(a));
    else if (S.sort === 'size') order.sort((a, b) => D.species[b].n_proteins - D.species[a].n_proteins);
    else order.sort((a, b) => String(bandOf(a) ?? '￿').localeCompare(String(bandOf(b) ?? '￿'))
      || tot(b) - tot(a));
    // With hundreds of rows a label per row is unreadable and the ink hides the
    // data. Show every k-th name, and when rows are grouped put the group name
    // once against its band instead.
    const grouped = !S.sortFam && (S.sort === 'trait' || taxSort);
    tabHeatmap._lastOrder = order.map(i => [SP[i], bandOf(i)]);
    const LBL_MAX = 60;
    const step = Math.max(1, Math.ceil(order.length / LBL_MAX));

    // Band labels replace per-row names once rows are grouped and dense: at
    // that density a name per row is unreadable and the ink hides the data.
    const bandLabels = grouped && order.length > LBL_MAX;
    // Size the label gutter to the labels actually drawn, rather than
    // reserving a fixed 210px that is half empty for short names and too
    // narrow for long ones.
    const mg = el('canvas').getContext('2d');
    const IT = 'italic 11px sans-serif', RO = '11px sans-serif';
    let labW = 0;
    if (bandLabels) {
      mg.font = '600 11.5px sans-serif';
      new Set(order.map(i => String(bandOf(i) ?? '\u2014'))).forEach(v => {
        labW = Math.max(labW, mg.measureText(`${nice(v)} (0000)`).width); });
    } else {
      order.forEach((si, r) => {
        if (!(order.length <= LBL_MAX || r % step === 0)) return;
        mg.font = IT; let w = mg.measureText(sciName(SP[si])).width;
        const c = commonName(SP[si]), q = qualifier(SP[si]);
        if (c || q) { mg.font = RO;
          w += mg.measureText((c ? ` (${c})` : '') + (q ? ` \u00b7 ${q}` : '')).width; }
        labW = Math.max(labW, w); });
    }
    const CH = 17, R = 30;
    const L = Math.round(Math.min(330, Math.max(95, labW + 30)));
    // The header block is rotated 75 degrees, so the room it needs is set by
    // the longest family name, not by a constant that is too tall for eight
    // families and too short for forty-four.
    mg.font = '11px sans-serif';
    const hdrW = Math.max(...fams.map(f => mg.measureText(f).width), 10);
    const T = Math.round(Math.min(230, Math.max(46, hdrW * 0.97 + 18)));
    // Spread the cells across the width that is actually there. Eight families
    // in a 1,400px panel should not be a 200px thumbnail.
    const avail = (p.clientWidth || 1134) - 50;   // panel width less padding and scrollbar
    const CW = Math.max(26, Math.min(60, Math.floor((avail - L - R) / Math.max(fams.length, 1))));
    const WIDE = L + fams.length * CW + R;
    hdr.width = WIDE; hdr.height = T;
    cv.width = WIDE; cv.height = order.length * CH + 2;
    const hg = hdr.getContext('2d'); hg.clearRect(0, 0, WIDE, T);
    const g = cv.getContext('2d'); g.clearRect(0, 0, cv.width, cv.height);
    const vals = {}; fams.forEach(f => {
      vals[f] = S.scale === 'z' ? zcol(S.kind, f)
        : S.scale === 'log' ? col(S.kind, f).map(x => Math.log10(x + 1))
        : col(S.kind, f); });
    let lo, hi;
    if (S.scale === 'z') { lo = -2.5; hi = 2.5; }
    else { const all = fams.flatMap(f => vals[f]); lo = 0; hi = Math.max(...all) || 1; }
    hg.font = '11px sans-serif';
    fams.forEach((f, j) => { hg.save(); hg.translate(L + j * CW + CW / 2, T - 8);
      hg.rotate(-Math.PI / 2.4); hg.textAlign = 'left';
      hg.fillStyle = INVARIANT[f] ? '#9aa1a8' : '#1c1f23'; hg.fillText(f, 0, 0); hg.restore(); });
    order.forEach((si, r) => {
      g.fillStyle = traitColour(D.species[si].trait);
      g.fillRect(L - 12, r * CH + 2, 7, CH - 4);
      if (!bandLabels && (order.length <= LBL_MAX || r % step === 0)) {
        g.fillStyle = '#1c1f23'; g.textBaseline = 'alphabetic';
        drawSpeciesLabel(g, SP[si], L - 18, r * CH + CH - 5, L - 24);
      }
      fams.forEach((f, j) => { const v = vals[f][si];
        g.fillStyle = S.scale === 'z' ? rdbu((v - lo) / (hi - lo)) : viridis((v - lo) / (hi - lo));
        g.fillRect(L + j * CW, r * CH, CW - 1, CH - 1); });
    });
    tabHeatmap._geom = () => ({L, T, CW, CH, rows: order.length, fams});
    tabHeatmap._lastBands = [];
    if (bandLabels) {       // one label per contiguous run of the grouping variable
      const tr = i => String(bandOf(order[i]) ?? '—');
      let run0 = 0;
      for (let r = 1; r <= order.length; r++) {
        if (r === order.length || tr(r) !== tr(run0)) {
          const y0 = run0 * CH, y1 = r * CH;
          if (y1 - y0 > 14) {
            g.save(); g.fillStyle = '#1c1f23'; g.font = '600 11.5px sans-serif';
            g.textAlign = 'right'; g.textBaseline = 'middle';
            g.fillText(`${nice(tr(run0))} (${r - run0})`, L - 18, (y0 + y1) / 2);
            g.restore();
            tabHeatmap._lastBands.push([tr(run0), r - run0]);
          }
          run0 = r;
        }
      }
    }
    cv.onmousemove = e => { const b = cv.getBoundingClientRect();
      const j = Math.floor((e.clientX - b.left - L) / CW), r = Math.floor((e.clientY - b.top) / CH);
      if (j < 0 || r < 0 || j >= fams.length || r >= order.length) return hideTip();
      const f = fams[j], si = order[r];
      showTip(e, `<b>${fullName(SP[si])}</b><br>${f} (${FMETA[f].category}, ${FMETA[f].tier})<br>
        raw ${D.counts[f][si]}${D.per10k ? ' · per10k ' + D.per10k[f][si].toFixed(2) : ''}<br>
        z ${zcol(S.kind, f)[si].toFixed(2)}${INVARIANT[f] ? '<br><i>' + INVARIANT[f] + '</i>' : ''}
        ${D.species[si].trait ? '<br>' + D.trait_name + ': ' + D.species[si].trait : ''}`); };
    cv.onmouseleave = hideTip;
    cv.onclick = e => { const b = cv.getBoundingClientRect();
      const j = Math.floor((e.clientX - b.left - L) / CW), r = Math.floor((e.clientY - b.top) / CH);
      if (j >= 0 && j < fams.length && r >= 0 && r < order.length) {
        location.hash = famHref(fams[j]); } };
    hdr.style.cursor = 'pointer';
    hdr.onclick = e => { const b = hdr.getBoundingClientRect();
      const j = Math.floor((e.clientX - b.left - L) / CW);
      if (j >= 0 && j < fams.length) { S.sortFam = fams[j]; draw(); } };
    hdr.onmousemove = e => { const b = hdr.getBoundingClientRect();
      const j = Math.floor((e.clientX - b.left - L) / CW);
      if (j < 0 || j >= fams.length) return hideTip();
      showTip(e, `<b>${fams[j]}</b><br>${FMETA[fams[j]].category}, ${FMETA[fams[j]].tier}` +
        '<br>click to sort the species by this family'); };
    hdr.onmouseleave = hideTip;
    cbar.innerHTML = '';
    cbar.append(colourbar(S.scale === 'z' ? rdbu : viridis, lo, hi,
      S.scale === 'z' ? `z-score, ${S.kind === 'per10k' ? `copies per 10k ${D.size_unit || 'proteins'}` : 'raw copies'}`
                      : (S.scale === 'log' ? 'log10(x+1)' : (S.kind === 'per10k' ? `copies per 10k ${D.size_unit || 'proteins'}` : 'raw copies')),
      280, 12));
    leg.innerHTML = '';
    TRAITS.forEach(t => leg.append(el('span', {}, el('i', {class: 'sw',
      style: 'background:' + traitColour(t)}), t)));
    const hidden = FAM.filter(f => S.scale === 'z' && INVARIANT[f]
      && S.tiers.has(FMETA[f].tier) && S.cats.has(FMETA[f].category));
    if (hidden.length) leg.append(el('span', {style: 'color:#6b7280'},
      '· hidden (no variance, cannot be z-scored): ' + hidden.join(', ')));
  }
  draw();
  // test hooks: the drawn row order and the live control state
  tabHeatmap._state = S; tabHeatmap._canvas = () => cv;
  tabHeatmap._header = () => hdr;
  tabHeatmap._sortSel = () => sortSel;
  tabHeatmap._setSort = v => { sortSel.value = v;
    sortSel.dispatchEvent(new Event('change')); };
}

/* --------------------------------------------------------------- FAMILY */
/* A family page, like a species page: what the family is, how it is spread
   across the species in view, its members with links to each gene, and the
   trees built for it. The filter bar applies; the switcher changes family. */
let famSel, FAM_CUR = null;
function tabFamily(host) {
  const fam = FAM.includes(FAM_CUR) ? FAM_CUR : (FAM.includes('CYP') ? 'CYP' : FAM[0]);
  FAM_CUR = fam;
  const m = FMETA[fam];
  const head = el('div', {class: 'ehead'});
  famSel = el('select', {}, ...FAM.map(f => el('option', {value: f}, f)));
  famSel.value = fam;
  const kindSel = el('select', {}, ...(D.per10k ? [el('option', {value: 'per10k'}, `per 10k ${D.size_unit || 'proteins'}`)] : []),
    el('option', {value: 'raw'}, 'raw copies'));
  head.append(el('div', {style: 'display:flex;justify-content:space-between;gap:12px;flex-wrap:wrap'},
    el('div', {}, crumb(link('#/families', 'Families'), el('span', {}, m.category), el('span', {}, fam)),
      el('h1', {}, fam), el('div', {class: 'sub'}, (FAMDESC[fam] ? FAMDESC[fam] + ' · ' : '') + m.category)),
    el('div', {class: 'row', style: 'margin:0'}, el('div', {class: 'ctl'}, el('label', {}, 'family'), famSel),
      el('div', {class: 'ctl'}, el('label', {}, 'values'), kindSel))));
  head.append(el('div', {class: 'meta'}, ...String(m.pfam).split(',').map(x =>
    extLink(interpro(x.trim()), 'InterPro ' + x.trim()))));
  const jump = el('div', {class: 'jump'}, ...[['fam-dist', 'Distribution'], ['fam-species', 'Species'],
    ['fam-domains', 'Domains'], ['fam-members', 'Members'], ['fam-trees', 'Trees'], ['fam-downloads', 'Downloads']]
    .map(([id, t]) => el('a', {href: '#' + id, onclick: e => { e.preventDefault();
      const n = document.getElementById ? document.getElementById(id) : null; if (n && n.scrollIntoView) n.scrollIntoView({behavior: 'smooth'}); }}, t)));
  head.append(jump);
  host.append(head);
  const body = el('div'); host.append(body);
  const draw = () => { body.innerHTML = ''; renderFamily(body, fam, pick(kindSel, ['per10k', 'raw'])); };
  famSel.onchange = () => { FAM_CUR = famSel.value; go(famHref(famSel.value)); };
  kindSel.onchange = draw; draw();
  familyExtras(host, fam);
}

function renderFamily(host, fam, kind) {
  const m = FMETA[fam];
  if (!m) { host.append(el('div', {class: 'note'}, 'Unknown family: ' + fam)); return; }
  const v = col(kind, fam), st = stats(v);
  const nPres = D.counts[fam].filter(x => x > 0).length;
  host.append(el('div', {class: 'grid g4'},
    kpi(`${nPres.toLocaleString()} of ${SP.length.toLocaleString()}`, 'species with a copy'),
    kpi(fmt(st.med), 'median ' + (kind === 'per10k' ? `per 10k ${D.size_unit || 'proteins'}` : 'copies')),
    kpi(fmt(st.min) + ' – ' + fmt(st.max), 'range'),
    kpi(D.counts[fam].reduce((a, b) => a + b, 0).toLocaleString(), 'genes in view')));

  /* distribution by group: full width, so it can grow to as many groups as
     the filter leaves in view; the numbers behind it are one click away */
  if (TRAITS.length) {
    const c = el('div', {class: 'panel', id: 'fam-dist'});
    c.append(el('h2', {}, 'By ' + D.trait_name), el('p', {class: 'hint'},
      'Species are not independent, so treat any separation here as a screen and not a result.'));
    const groups = TRAITS.map(t => ({t, v: SP.map((s, i) => [D.species[i].trait, v[i]])
      .filter(x => x[0] === t).map(x => x[1])}));
    const W = Math.max(520, Math.min(1100, (host.clientWidth || 1000) - 40));
    c.append(figure(boxPlot(groups, groups.length > BOX_VERTICAL_MAX ? W : 560, 300), fam + '_by_trait'));
    const rows = groups.filter(g => g.v.length).map(g => { const q = stats(g.v);
      return [g.t, g.v.length, +fmt(q.mean), +fmt(q.med), q.min, q.max]; });
    const box = el('div');
    makeTable(box, [D.trait_name, 'n', 'mean', 'median', 'min', 'max'], rows,
      {sortBy: 3, desc: true, name: fam + '_by_' + D.trait_name});
    c.append(collapsible(`numbers by ${D.trait_name} (${rows.length})`, box));
    host.append(c);
  }
  /* domain completeness: one bar and one line; the detail on request */
  if (D.completeness && D.completeness[fam]) {
    const c = el('div', {class: 'panel', id: 'fam-domains'});
    const cc = D.completeness[fam], tot = Object.values(cc).reduce((a, b) => a + b, 0);
    const cols = {COMPLETE: '#2a9d8f', PARTIAL: '#e9a13b', FRAGMENT: '#bc4749', MISSING_DOMAIN: '#adb5bd'};
    const pc = k => tot ? (100 * (cc[k] || 0) / tot).toFixed(1) : '0';
    c.append(el('h2', {}, `Domain completeness: ${pc('COMPLETE')}% complete`),
      el('p', {class: 'hint'}, 'PARTIAL: a required domain is truncated below the family threshold. FRAGMENT: far shorter than the family norm.'));
    const bar = el('div', {style: 'display:flex;height:14px;border-radius:3px;overflow:hidden;margin:6px 0'});
    Object.entries(cc).forEach(([k, n]) => { if (!n) return;
      bar.append(el('div', {style: `width:${n / tot * 100}%;background:${cols[k] || '#ccc'}`,
        title: `${k}: ${n}`})); });
    c.append(bar, el('div', {class: 'legend', style: 'margin-top:4px'},
      Object.entries(cc).filter(([, n]) => n).map(([k, n]) =>
        el('span', {}, el('i', {class: 'sw', style: 'background:' + (cols[k] || '#ccc')}),
          `${k.toLowerCase().replace('_', ' ')} ${n.toLocaleString()}`))));
    if (D.comp_sf) {
      const rows = D.comp_sf.filter(r => r[1] === fam && SPI[r[0]] !== undefined)
        .map(r => [sciName(r[0]), r[2], r[3], r[4], r[5],
                   +(r[2] / Math.max(1, r[2] + r[3] + r[4] + r[5]) * 100).toFixed(1)]);
      if (rows.length) { const box = el('div');
        makeTable(box, ['species', 'complete', 'partial', 'fragment', 'missing', '% complete'],
          rows, {sortBy: 5, ital: [0], name: fam + '_completeness'});
        c.append(collapsible(`per species (${rows.length.toLocaleString()}, least complete first)`, box)); }
    }
    if (D.architectures) { const ar = D.architectures.filter(r => r[0] === fam).slice(0, 12);
      if (ar.length) c.append(collapsible(`domain architectures (top ${ar.length})`,
        el('table', {}, el('thead', {}, el('tr', {},
          ['architecture', 'proteins'].map(h => el('th', {}, h)))),
          el('tbody', {}, ar.map(r => el('tr', {}, el('td', {class: 'mono'}, r[1]),
            el('td', {}, r[2]))))))); }
    host.append(c);
  }

  /* per-species bar */
  const c3 = el('div', {class: 'panel', id: 'fam-species'});
  c3.append(el('h2', {}, 'Per species'), el('p', {class: 'hint'}, 'Ranked; click a name for the species page.'));
  const order = SP.map((s, i) => i).sort((a, b) => v[b] - v[a]);
  // A bar per species is a 27,000-pixel scroll at catalogue scale. Show the top
  // of the ranking by default and let the reader open the rest; the heatmap and
  // the exported TSV are the places to read all of it.
  const BAR_TOP = 20;
  const barHold = el('div');
  const drawBars = all => { barHold.innerHTML = '';
    const rows = all ? order : order.slice(0, BAR_TOP);
    barHold.append(figure(barChart(rows.map(i => ({label: fullName(SP[i], 32), value: v[i], href: spHref(SP[i]),
      colour: traitColour(D.species[i].trait),
      extra: `raw ${D.counts[fam][i]} · ${D.species[i].n_proteins.toLocaleString()} ${D.size_unit || 'proteins'}`})),
      880, Math.max(220, rows.length * 17)), fam + '_per_species'));
    if (TRAITS.length > 1 && TRAITS.length <= 12) barHold.append(el('div', {class: 'legend'},
      TRAITS.map(t => el('span', {}, el('i', {class: 'sw', style: 'background:' + traitColour(t)}), nice(t)))));
    if (order.length > BAR_TOP) barHold.append(el('button', {class: 'act',
      style: 'margin-top:8px', onclick: () => drawBars(!all)},
      all ? `show top ${BAR_TOP} only` : `show all ${order.length} species`));
  };
  drawBars(order.length <= BAR_TOP);
  c3.append(barHold);
  host.append(c3);

}

/* --------------------------------------------------------- SVG PRIMITIVES */
const svgEl = (t, a = {}) => { const n = document.createElementNS('http://www.w3.org/2000/svg', t);
  for (const [k, v] of Object.entries(a)) n.setAttribute(k, v); return n; };

/* Box plot that stays readable at any number of groups.
   Up to 8 groups: the vertical layout. Beyond that, one horizontal row per
   group sorted by median, so names read left to right with no rotation and the
   order itself carries the comparison. Groups with fewer than MIN_N species
   are left out, with a line saying how many: a "box" of one or two species is
   a point pretending to be a distribution. Individual species are drawn only
   while there are few enough to see; past that they are a grey smear that
   hides the boxes. */
const BOX_MIN_N = 3, BOX_MAX_DOTS = 300, BOX_VERTICAL_MAX = 8;
function boxPlot(groups, W, H) {
  const keep = groups.filter(g => g.v.length >= BOX_MIN_N);
  const dropped = groups.length - keep.length;
  const wrap = node => {
    if (!dropped) return node;
    const n = groups.filter(g => g.v.length < BOX_MIN_N).reduce((a, g) => a + g.v.length, 0);
    return el('div', {}, node, el('div', {class: 'hint', style: 'margin:4px 0 0'},
      `${dropped} group${dropped === 1 ? '' : 's'} with fewer than ${BOX_MIN_N} species not shown (${n} species).`));
  };
  if (!keep.length) return el('div', {class: 'note'}, `No group has ${BOX_MIN_N} or more species.`);
  const dots = keep.reduce((a, g) => a + g.v.length, 0) <= BOX_MAX_DOTS;
  if (keep.length > BOX_VERTICAL_MAX) return wrap(hBoxPlot(keep, Math.max(W, 560), dots));
  return wrap(vBoxPlot(keep, W, H, dots));
}

function hBoxPlot(groups, W, dots) {
  const gs = groups.map(g => ({...g, s: stats(g.v)}))
    .sort((a, b) => b.s.med - a.s.med || b.v.length - a.v.length);
  // past the palette, colour would repeat across unrelated groups and suggest
  // a link that is not there; the label already says which group a row is
  const one = gs.length > TCOL.length;
  const all = gs.flatMap(g => g.v);
  const lo = Math.min(...all), hi = Math.max(...all), pad = (hi - lo) * .05 || 1;
  const x0 = lo - pad, x1 = hi + pad;
  const mg = measurer(11.5);
  const labs = gs.map(g => `${g.t} (${g.v.length})`);
  const ML = Math.ceil(Math.min(W * .4, Math.max(...labs.map(t => mg.measureText(t).width)) + 16));
  const MR = 16, MT = 6, MB = 30, rh = 22;
  const H = MT + gs.length * rh + MB;
  const X = v => ML + (W - ML - MR) * (v - x0) / (x1 - x0);
  const sv = svgEl('svg', {width: W, height: H, viewBox: `0 0 ${W} ${H}`});
  for (let k = 0; k <= 4; k++) { const v = x0 + (x1 - x0) * k / 4;
    sv.append(svgEl('line', {x1: X(v), x2: X(v), y1: MT, y2: H - MB, stroke: '#eceff2'}));
    const t = svgEl('text', {x: X(v), y: H - MB + 16, 'text-anchor': 'middle',
      'font-size': 10, fill: '#6b7280'}); t.textContent = fmt(v, 1); sv.append(t); }
  gs.forEach((g, i) => { const y = MT + i * rh + rh / 2, col = one ? '#1b4965' : traitColour(g.t), s = g.s;
    sv.append(svgEl('line', {x1: X(s.min), x2: X(s.max), y1: y, y2: y, stroke: '#9aa1a8'}));
    const r = svgEl('rect', {x: X(s.q1), y: y - 7, width: Math.max(1, X(s.q3) - X(s.q1)),
      height: 14, fill: col, 'fill-opacity': .35, stroke: col});
    r.addEventListener('mousemove', e => showTip(e, `<b>${g.t}</b> · ${g.v.length} species<br>` +
      `median ${fmt(s.med, 2)} · IQR ${fmt(s.q1, 2)} to ${fmt(s.q3, 2)}<br>range ${fmt(s.min, 2)} to ${fmt(s.max, 2)}`));
    r.addEventListener('mouseleave', hideTip); sv.append(r);
    sv.append(svgEl('line', {x1: X(s.med), x2: X(s.med), y1: y - 7, y2: y + 7,
      stroke: col, 'stroke-width': 2.5}));
    if (dots) g.v.forEach((v, q) => sv.append(svgEl('circle', {cx: X(v),
      cy: y + ((q * 2654435761 % 100) / 100 - .5) * 10, r: 1.8, fill: col, 'fill-opacity': .6})));
    const t = svgEl('text', {x: ML - 8, y: y + 4, 'text-anchor': 'end', 'font-size': 11.5,
      class: 'axlab'}); t.textContent = labs[i]; sv.append(t); });
  return sv;
}

function vBoxPlot(groups, W, H, dots = true) {
  const all = groups.flatMap(g => g.v); if (!all.length) return el('div', {}, 'no data');
  const lo = Math.min(...all), hi = Math.max(...all), pad = (hi - lo) * .1 || 1;
  const y0 = lo - pad, y1 = hi + pad;
  const MT = 12, MR = 12;
  let ML = 56;
  let bw = (W - ML - MR) / groups.length;
  /* Category names are as long as they are; the axis has to adapt to them.
     Measure the longest, and if it will not fit upright in one band, lay the
     labels at an angle and give the axis the depth that needs. Past the point
     where even angled labels collide, show every k-th. */
  const mg = measurer(11);
  const labels = groups.map(g => `${g.t} (n=${g.v.length})`);
  const maxLab = Math.max(...labels.map(t => mg.measureText(t).width), 0);
  const rotate = maxLab > bw * 0.92;
  if (rotate) {
    // A label laid at -38 degrees runs left of its band, so the leftmost one
    // needs the axis gutter widened or it is cut off at the edge.
    ML = Math.max(ML, Math.ceil(maxLab * 0.79 - bw / 2 + 8));
    ML = Math.min(ML, Math.round(W * 0.35));
    bw = (W - ML - MR) / groups.length;
  }
  const gstep = rotate ? Math.max(1, Math.ceil(13 / Math.max(bw, 1))) : 1;
  const MB = rotate ? Math.round(Math.min(150, Math.max(54, maxLab * 0.63 + 30))) : 44;
  const sv = svgEl('svg', {width: W, height: H, viewBox: `0 0 ${W} ${H}`});
  const Y = v => MT + (H - MT - MB) * (1 - (v - y0) / (y1 - y0));
  for (let k = 0; k <= 4; k++) { const v = y0 + (y1 - y0) * k / 4;
    sv.append(svgEl('line', {x1: ML, x2: W - MR, y1: Y(v), y2: Y(v), stroke: '#eceff2'}));
    const t = svgEl('text', {x: ML - 8, y: Y(v) + 4, 'text-anchor': 'end',
      'font-size': 10, fill: '#6b7280'}); t.textContent = v.toFixed(1); sv.append(t); }
  groups.forEach((g, j) => { const cx = ML + bw * (j + .5), col = traitColour(g.t);
    if (g.v.length) { const s = stats(g.v);
      sv.append(svgEl('line', {x1: cx, x2: cx, y1: Y(s.min), y2: Y(s.max), stroke: '#556'}));
      sv.append(svgEl('rect', {x: cx - bw * .22, y: Y(s.q3), width: bw * .44,
        height: Math.max(1, Y(s.q1) - Y(s.q3)), fill: col, 'fill-opacity': .3,
        stroke: col}));
      sv.append(svgEl('line', {x1: cx - bw * .22, x2: cx + bw * .22, y1: Y(s.med),
        y2: Y(s.med), stroke: col, 'stroke-width': 2}));
      if (dots) g.v.forEach((v, q) => sv.append(svgEl('circle', {
        cx: cx + ((q * 2654435761 % 100) / 100 - .5) * bw * .30, cy: Y(v), r: 2.6,
        fill: col, 'fill-opacity': .8, stroke: '#fff', 'stroke-width': .5}))); }
    if (j % gstep) return;
    const ly = H - MB + (rotate ? 14 : 16);
    const t = svgEl('text', {x: cx, y: ly, 'font-size': 11, class: 'axlab',
      'text-anchor': rotate ? 'end' : 'middle',
      ...(rotate ? {transform: `rotate(-38 ${cx} ${ly})`} : {})});
    t.textContent = rotate ? labels[j] : g.t; sv.append(t);
    if (!rotate) { const t2 = svgEl('text', {x: cx, y: H - MB + 30, 'text-anchor': 'middle',
      'font-size': 10, fill: '#6b7280'}); t2.textContent = 'n=' + g.v.length; sv.append(t2); } });
  return sv;
}

function barChart(items, W, H, signed = false) {
  // the label gutter fits the longest label, so names are never cut at the edge
  const mg = measurer(11);
  const ML = Math.round(Math.min(330, Math.max(80, ...items.map(d => mg.measureText(String(d.label)).width)) + 16));
  const MR = 60, MT = 8;
  // a few rows should not be stretched into 50-pixel slabs to fill a height
  // computed for many; cap the row height and let the figure be shorter
  const rows = items.length, rh = Math.min(22, Math.max(9, (H - MT) / rows));
  const hh = MT + rows * rh;
  const mx = Math.max(...items.map(d => Math.abs(d.value))) || 1;
  const PADL = signed ? 46 : 0;
  const zero = signed ? ML + PADL + (W - ML - MR - PADL) / 2 : ML;
  const span = signed ? (W - ML - MR - PADL) / 2 : (W - ML - MR);
  const sv = svgEl('svg', {width: W, height: hh, viewBox: `0 0 ${W} ${hh}`});
  if (signed) sv.append(svgEl('line', {x1: zero, x2: zero, y1: MT, y2: hh, stroke: '#ccc'}));
  items.forEach((d, i) => { const y = MT + i * rh, w = Math.abs(d.value) / mx * span;
    const x = signed ? (d.value >= 0 ? zero : zero - w) : ML;
    const r = svgEl('rect', {x, y: y + 1, width: Math.max(1, w), height: rh - 2,
      fill: d.colour || '#1b4965'});
    r.addEventListener('mousemove', e => showTip(e,
      `<b>${d.label}</b><br>${fmt(d.value, 3)}${d.extra ? '<br>' + d.extra : ''}`));
    r.addEventListener('mouseleave', hideTip); sv.append(r);
    if (rh >= 9) { const t = svgEl('text', {x: ML - 8, y: y + rh / 2 + 3.5,
      'text-anchor': 'end', 'font-size': Math.min(11, rh - 1), fill: '#1c1f23'});
      t.textContent = d.label;
      if (d.href) { t.style.cursor = 'pointer'; t.setAttribute('fill', '#1b4965');
        t.addEventListener('click', () => { location.hash = d.href; }); }
      sv.append(t); }
    const t2 = svgEl('text', {x: (signed && d.value < 0 ? x - 5 : x + w + 5),
      y: y + rh / 2 + 3.5, 'text-anchor': signed && d.value < 0 ? 'end' : 'start',
      'font-size': Math.min(10, rh - 1), fill: '#6b7280'});
    t2.textContent = d.vlabel != null ? d.vlabel : fmt(d.value, 2); sv.append(t2); });
  return sv;
}

/* ------------------------------------------------------------ ORDINATION */
function jacobiEig(A, n) {   /* symmetric eigendecomposition, enough for 40x40 */
  const V = Array.from({length: n}, (_, i) => Array.from({length: n}, (_, j) => i === j ? 1 : 0));
  const a = A.map(r => r.slice());
  for (let sweep = 0; sweep < 60; sweep++) {
    let off = 0; for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) off += a[i][j] ** 2;
    if (off < 1e-12) break;
    for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) {
      if (Math.abs(a[i][j]) < 1e-14) continue;
      const th = (a[j][j] - a[i][i]) / (2 * a[i][j]);
      const t = Math.sign(th || 1) / (Math.abs(th) + Math.sqrt(th * th + 1));
      const c = 1 / Math.sqrt(t * t + 1), s = t * c;
      for (let k = 0; k < n; k++) { const aik = a[i][k], ajk = a[j][k];
        a[i][k] = c * aik - s * ajk; a[j][k] = s * aik + c * ajk; }
      for (let k = 0; k < n; k++) { const aki = a[k][i], akj = a[k][j];
        a[k][i] = c * aki - s * akj; a[k][j] = s * aki + c * akj; }
      for (let k = 0; k < n; k++) { const vki = V[k][i], vkj = V[k][j];
        V[k][i] = c * vki - s * vkj; V[k][j] = s * vki + c * vkj; } }
  }
  const idx = Array.from({length: n}, (_, i) => i).sort((x, y) => a[y][y] - a[x][x]);
  return {values: idx.map(i => a[i][i]), vectors: idx.map(i => V.map(r => r[i]))};
}

function tabOrdination(host) {
  const S = {kind: D.per10k ? 'per10k' : 'raw', fams: new Set(
    VARFAM.filter(f => FMETA[f].tier === 'CORE'))};
  const p = el('div', {class: 'panel'});
  p.append(el('h2', {}, 'Ordination'), el('p', {class: 'hint'},
    'Each family is z-scored first so a 150-copy family does not outweigh a 5-copy one, then the covariance matrix is decomposed. The loadings panel says which families drive the separation. If PC1 explains only 30 to 40 percent, do not read this as a single defence gradient.'));
  const kindSel = el('select', {onchange: e => { S.kind = e.target.value; draw(); }},
    ...(D.per10k ? [el('option', {value: 'per10k'}, `per 10k ${D.size_unit || 'proteins'}`)] : []),
    el('option', {value: 'raw'}, 'raw copies'));
  const chips = el('div', {class: 'chips'});
  VARFAM.forEach(f => chips.append(el('span', {class: 'chip' + (S.fams.has(f) ? ' on' : ''),
    onclick: e => { S.fams.has(f) ? S.fams.delete(f) : S.fams.add(f);
      e.target.classList.toggle('on'); draw(); }}, f)));
  p.append(el('div', {class: 'row'}, el('div', {class: 'ctl'}, el('label', {}, 'values'), kindSel),
    el('button', {class: 'act', onclick: () => { S.fams = new Set(VARFAM.filter(f =>
      FMETA[f].tier === 'CORE')); refreshChips(); draw(); }}, 'CORE only'),
    el('button', {class: 'act', onclick: () => { S.fams = new Set(VARFAM); refreshChips(); draw(); }}, 'all'),
    el('button', {class: 'act', onclick: () => { S.fams = new Set(); refreshChips(); draw(); }}, 'none')));
  p.append(el('div', {class: 'ctl'}, el('label', {}, 'families included'), chips));
  const out = el('div'); p.append(out); host.append(p);
  const refreshChips = () => [...chips.children].forEach(c =>
    c.classList.toggle('on', S.fams.has(c.textContent)));
  function draw() { out.innerHTML = '';
    const fams = VARFAM.filter(f => S.fams.has(f));
    if (fams.length < 3) { out.append(el('div', {class: 'note'}, 'Select at least three families.')); return; }
    const Z = SP.map((_, i) => fams.map(f => zcol(S.kind, f)[i]));
    const n = fams.length, N = SP.length;
    const C = Array.from({length: n}, () => new Array(n).fill(0));
    for (let i = 0; i < n; i++) for (let j = i; j < n; j++) {
      let s = 0; for (let k = 0; k < N; k++) s += Z[k][i] * Z[k][j];
      C[i][j] = C[j][i] = s / (N - 1); }
    const {values, vectors} = jacobiEig(C, n);
    const tot = values.reduce((a, b) => a + Math.max(0, b), 0);
    const scores = Z.map(r => [0, 1].map(k => r.reduce((a, v, q) => a + v * vectors[k][q], 0)));
    out.append(el('div', {class: 'grid g2'},
      wrapPanel('Species scores', figure(scatter(scores.map((s, i) => ({x: s[0], y: s[1],
        label: nice(SP[i]), colour: traitColour(D.species[i].trait),
        extra: (D.species[i].trait ?? '') + ' · ' + D.species[i].n_proteins.toLocaleString() + ' proteins'})),
        620, 520, `PC1 (${(values[0] / tot * 100).toFixed(1)}%)`,
        `PC2 (${(values[1] / tot * 100).toFixed(1)}%)`, true), 'pca_scores')),
      wrapPanel('Family loadings', figure(scatter(fams.map((f, q) => ({x: vectors[0][q],
        y: vectors[1][q], label: f, colour: catColour(FMETA[f].category),
        extra: FMETA[f].category})), 620, 520, 'PC1 loading', 'PC2 loading', fams.length),
        'pca_loadings'))));
    const varRows = values.slice(0, 8).map((v, k) => [`PC${k + 1}`,
      +(v / tot * 100).toFixed(2), +(values.slice(0, k + 1).reduce((a, b) => a + b, 0) / tot * 100).toFixed(2)]);
    const vp = el('div', {class: 'panel'});
    vp.append(el('h2', {}, 'Variance explained'));
    makeTable(vp, ['component', '% variance', 'cumulative %'], varRows, {name: 'pca_variance'});
    out.append(vp);
  }
  draw();
}
const wrapPanel = (title, node) => el('div', {class: 'panel'}, el('h2', {}, title), node);

function scatter(pts, W, H, xl, yl, labels) {
  const ML = 56, MB = 44, MT = 12, MR = 16;
  const xs = pts.map(p => p.x), ys = pts.map(p => p.y);
  const pad = v => { const lo = Math.min(...v), hi = Math.max(...v), d = (hi - lo) * .12 || 1;
    return [lo - d, hi + d]; };
  const [x0, x1] = pad(xs), [y0, y1] = pad(ys);
  const X = v => ML + (W - ML - MR) * (v - x0) / (x1 - x0);
  const Y = v => MT + (H - MT - MB) * (1 - (v - y0) / (y1 - y0));
  const sv = svgEl('svg', {width: W, height: H, viewBox: `0 0 ${W} ${H}`});
  sv.append(svgEl('line', {x1: X(0), x2: X(0), y1: MT, y2: H - MB, stroke: '#e3e6ea'}));
  sv.append(svgEl('line', {x1: ML, x2: W - MR, y1: Y(0), y2: Y(0), stroke: '#e3e6ea'}));
  pts.forEach(p => { const c = svgEl('circle', {cx: X(p.x), cy: Y(p.y), r: 5,
    fill: p.colour, 'fill-opacity': .85, stroke: '#fff', 'stroke-width': 1});
    c.addEventListener('mousemove', e => showTip(e, `<b>${p.label}</b>` +
      (p.extra ? '<br>' + p.extra : '') + `<br>${p.x.toFixed(3)}, ${p.y.toFixed(3)}`));
    c.addEventListener('mouseleave', hideTip); sv.append(c); });
  /* Labelling every point turns a thousand-species ordination into a grey
     smear. Label the points furthest from the centroid, which are the ones
     worth naming, and only where no label has been placed already; the rest
     are still on hover. */
  if (labels) {
    /* Labelling every point turns a hundred-species ordination into a grey
       smear. Label the points furthest from the centroid, which are the ones
       worth naming, skip any whose box would land on a label already placed or
       run off the panel, and leave the rest to the tooltip. Boxes are measured,
       not guessed, because a guessed width lets long names overlap. */
    const mg = el('canvas').getContext('2d');
    mg.font = '8.5px sans-serif';
    const FS = 8.5, PAD = 2;
    const mx = xs.reduce((a, b) => a + b, 0) / pts.length;
    const my = ys.reduce((a, b) => a + b, 0) / pts.length;
    const cap = typeof labels === 'number' ? labels : (pts.length <= 40 ? pts.length : 30);
    const rank = pts.map((p, i) => [i, (p.x - mx) ** 2 + (p.y - my) ** 2])
      .sort((a, b) => b[1] - a[1]).slice(0, cap).map(r => r[0]);
    const placed = [];
    const clear = (x0, y0, x1, y1) => !placed.some(b =>
      x0 < b[2] + PAD && x1 + PAD > b[0] && y0 < b[3] + PAD && y1 + PAD > b[1]);
    let drawn = 0;
    rank.forEach(i => {
      const p = pts[i], px = X(p.x), py = Y(p.y);
      const w = mg.measureText(p.label).width, h = FS;
      // right of the point by default, left when that would leave the panel
      let x0 = px + 7;
      if (x0 + w > W - MR) x0 = px - 7 - w;
      if (x0 < ML) return;
      const y0 = py - h / 2, y1 = py + h / 2;
      if (y0 < MT || y1 > H - MB) return;
      if (!clear(x0, y0, x0 + w, y1)) return;
      placed.push([x0, y0, x0 + w, y1]); drawn++;
      const t = svgEl('text', {x: x0.toFixed(1), y: (py + 3).toFixed(1),
        'font-size': FS, fill: '#444'});
      t.textContent = p.label; sv.append(t);
    });
    if (drawn < pts.length) {
      const n = svgEl('text', {x: W - MR, y: MT + 10, 'text-anchor': 'end',
        'font-size': 9.5, fill: '#9aa1a8'});
      n.textContent = `${drawn} of ${pts.length} labelled, hover for the rest`;
      sv.append(n);
    }
  }
  const tx = svgEl('text', {x: (ML + W - MR) / 2, y: H - 8, 'text-anchor': 'middle',
    'font-size': 11, fill: '#6b7280'}); tx.textContent = xl; sv.append(tx);
  const ty = svgEl('text', {x: 14, y: (MT + H - MB) / 2, 'font-size': 11, fill: '#6b7280',
    transform: `rotate(-90 14 ${(MT + H - MB) / 2})`, 'text-anchor': 'middle'});
  ty.textContent = yl; sv.append(ty);
  return sv;
}

/* ---------------------------------------------------------- CORRELATION */
function tabCorrelation(host) {
  const S = {kind: D.per10k ? 'per10k' : 'raw', tier: 'CORE'};
  const p = el('div', {class: 'panel'});
  p.append(el('h2', {}, 'Do families expand together?'), el('p', {class: 'hint'},
    'Pearson correlation across species between every pair of families. Blocks of red are families that expand and contract as a unit. This is the same information the ordination rotates, without the rotation.'));
  const kindSel = el('select', {onchange: e => { S.kind = e.target.value; draw(); }},
    ...(D.per10k ? [el('option', {value: 'per10k'}, `per 10k ${D.size_unit || 'proteins'}`)] : []),
    el('option', {value: 'raw'}, 'raw copies'));
  const tierSel = el('select', {onchange: e => { S.tier = e.target.value; draw(); }},
    el('option', {value: 'CORE'}, 'CORE only'), el('option', {value: 'ALL'}, 'all families'));
  p.append(el('div', {class: 'row'}, el('div', {class: 'ctl'}, el('label', {}, 'values'), kindSel),
    el('div', {class: 'ctl'}, el('label', {}, 'tier'), tierSel)));
  const cv = el('canvas'); p.append(figure(cv, 'family_correlation'));
  const cbar = el('div', {style: 'margin-top:10px'}); p.append(cbar);
  const pairs = el('div'); p.append(pairs); host.append(p);
  function draw() {
    const fams = VARFAM.filter(f => S.tier === 'ALL' || FMETA[f].tier === 'CORE');
    const Zs = fams.map(f => zcol(S.kind, f)), n = fams.length, N = SP.length;
    const R = Array.from({length: n}, () => new Array(n).fill(0));
    for (let i = 0; i < n; i++) for (let j = i; j < n; j++) {
      let s = 0; for (let k = 0; k < N; k++) s += Zs[i][k] * Zs[j][k];
      R[i][j] = R[j][i] = s / (N - 1); }
    const ord = Array.from({length: n}, (_, i) => i)
      .sort((a, b) => R[b].reduce((x, y) => x + y, 0) - R[a].reduce((x, y) => x + y, 0));
    const CS = 20;
    // Gutters measured from the names, so a long family name is never clipped
    // and a short set is not given 130px of empty margin.
    const mgc = el('canvas').getContext('2d'); mgc.font = '10.5px sans-serif';
    const wMax = Math.max(...fams.map(f => mgc.measureText(f).width), 10);
    const L = Math.round(Math.min(220, wMax + 14)), T = L;
    cv.width = L + n * CS + 20; cv.height = T + n * CS + 20;
    const g = cv.getContext('2d'); g.clearRect(0, 0, cv.width, cv.height);
    g.font = '10.5px sans-serif';
    ord.forEach((fi, r) => { g.textAlign = 'right'; g.fillStyle = '#1c1f23';
      g.fillText(fams[fi], L - 6, T + r * CS + CS - 6);
      g.save(); g.translate(L + r * CS + CS - 6, T - 6); g.rotate(-Math.PI / 2);
      g.textAlign = 'left'; g.fillText(fams[fi], 0, 0); g.restore(); });
    ord.forEach((fi, r) => ord.forEach((fj, c) => {
      g.fillStyle = rdbu((R[fi][fj] + 1) / 2);
      g.fillRect(L + c * CS, T + r * CS, CS - 1, CS - 1); }));
    cv.onmousemove = e => { const b = cv.getBoundingClientRect();
      const c = Math.floor((e.clientX - b.left - L) / CS), r = Math.floor((e.clientY - b.top - T) / CS);
      if (c < 0 || r < 0 || c >= n || r >= n) return hideTip();
      showTip(e, `<b>${fams[ord[r]]}</b> vs <b>${fams[ord[c]]}</b><br>r = ${R[ord[r]][ord[c]].toFixed(3)}`); };
    cv.onmouseleave = hideTip;
    cbar.innerHTML = '';
    cbar.append(colourbar(rdbu, -1, 1, 'Pearson r', 280, 12));
    const rows = [];
    for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++)
      rows.push([fams[i], fams[j], +R[i][j].toFixed(3)]);
    pairs.innerHTML = ''; pairs.append(el('h2', {style: 'margin-top:16px'}, 'Pairs'));
    makeTable(pairs, ['family A', 'family B', 'r'], rows, {sortBy: 2, desc: true, name: 'family_correlations'});
  }
  draw();
}

/* -------------------------------------------------------------- DOMAINS */
function tabDomains(host) {
  if (!D.completeness) { host.append(el('div', {class: 'panel'},
    el('h2', {}, 'Domains'), el('p', {class: 'hint'},
    'No domain output found. Run: defensome.py domains --out <results>'))); return; }
  const p = el('div', {class: 'panel'});
  p.append(el('h2', {}, 'Architecture completeness'), el('p', {class: 'hint'},
    'A protein is COMPLETE when every domain its family requires is present above the family coverage threshold, using merged non-overlapping HMM segments. Families whose model matches in two pieces would look truncated if scored on the best single segment.'));
  const fams = Object.keys(D.completeness).sort((a, b) => {
    const ta = D.completeness[a], tb = D.completeness[b];
    const pa = ta.COMPLETE / Object.values(ta).reduce((x, y) => x + y, 0);
    const pb = tb.COMPLETE / Object.values(tb).reduce((x, y) => x + y, 0);
    return pa - pb; });
  const cols = {COMPLETE: '#2a9d8f', PARTIAL: '#e9a13b', FRAGMENT: '#bc4749', MISSING_DOMAIN: '#adb5bd'};
  const box = el('div');
  fams.forEach(f => { const cc = D.completeness[f], tot = Object.values(cc).reduce((a, b) => a + b, 0);
    const row = el('div', {style: 'display:flex;align-items:center;gap:10px;margin:3px 0'});
    row.append(el('div', {style: 'width:170px;font-size:12px'}, `${f} (${tot})`));
    const bar = el('div', {style: 'flex:1;display:flex;height:16px;border-radius:3px;overflow:hidden'});
    Object.entries(cc).forEach(([k, v]) => { if (!v) return;
      const d = el('div', {style: `width:${v / tot * 100}%;background:${cols[k]}`});
      d.addEventListener('mousemove', e => showTip(e, `<b>${f}</b><br>${k}: ${v} (${(v / tot * 100).toFixed(1)}%)`));
      d.addEventListener('mouseleave', hideTip); bar.append(d); });
    row.append(bar, el('div', {style: 'width:56px;text-align:right;font-size:12px;color:#6b7280'},
      (cc.COMPLETE / tot * 100).toFixed(1) + '%'));
    box.append(row); });
  p.append(box, el('div', {class: 'legend'}, Object.entries(cols).map(([k, v]) =>
    el('span', {}, el('i', {class: 'sw', style: 'background:' + v}), k))));
  host.append(p);

  if (D.architectures) { const p2 = el('div', {class: 'panel'});
    p2.append(el('h2', {}, 'Domain architectures'), el('p', {class: 'hint'},
      'Ordered Pfam domain strings observed per family. Unexpected architectures are usually fused or split gene models.'));
    makeTable(p2, ['family', 'architecture', 'proteins'], D.architectures,
      {sortBy: 2, desc: true, mono: [1], name: 'architectures', cap: 1200});
    host.append(p2); }

  if (D.flagged_proteins) { const p3 = el('div', {class: 'panel'});
    p3.append(el('h2', {}, 'Proteins that are not COMPLETE'), el('p', {class: 'hint'},
      'The individual gene models worth eyeballing. domain_coverages gives the merged HMM coverage per required domain.'));
    makeTable(p3, ['species', 'protein', 'family', 'status', 'length', 'length ratio',
      'min domain cov', 'domain coverages'], D.flagged_proteins,
      {sortBy: 6, mono: [1, 7], name: 'flagged_proteins', cap: 1000});
    host.append(p3); }
}

/* ---------------------------------------------------------------- CLANS */
function tabClans(host) {
  const keys = Object.keys(D.clans || {});
  if (!keys.length) { host.append(el('div', {class: 'panel'}, el('h2', {}, 'Clans'),
    el('p', {class: 'hint'}, 'No clan output found. Run: defensome.py cyp --out <results> --refs <clan-labelled FASTA>, then genetree.'))); return; }
  const p = el('div', {class: 'panel'});
  const sel = el('select', {}, ...keys.map(k => el('option', {value: k}, k)));
  const modeSel = el('select', {}, el('option', {value: 'frac'}, 'fraction of family'),
    el('option', {value: 'per10k'}, `per 10k ${D.size_unit || 'proteins'}`), el('option', {value: 'raw'}, 'raw counts'));
  p.append(el('h2', {}, 'Clan composition'), el('p', {class: 'hint'},
    'Fractions answer whether the balance between clans shifts; per-10k answers whether a clan expanded. They can disagree, and the disagreement is usually the interesting part.'),
    el('div', {class: 'row'}, el('div', {class: 'ctl'}, el('label', {}, 'family'), sel),
      el('div', {class: 'ctl'}, el('label', {}, 'scale'), modeSel)));
  const body = el('div'); p.append(body); host.append(p);
  const draw = () => { body.innerHTML = ''; renderClans(body, pick(sel, keys), pick(modeSel, ['frac','per10k','raw'])); };
  sel.onchange = draw; modeSel.onchange = draw; draw();
}

function renderClans(host, fam, mode) {
  const C = D.clans[fam];
  if (!C) { host.append(el('div', {class: 'note'}, 'No clan data for ' + fam)); return; }
  const clans = C.clans;
  const ccol = c => CCOL[clans.indexOf(c) % CCOL.length];
  const rows = SP.filter(s => C.data[s]).map(s => { const i = SPI[s], raw = C.data[s];
    const tot = raw.reduce((a, b) => a + b, 0);
    const v = mode === 'frac' ? raw.map(x => tot ? x / tot : 0)
      : mode === 'per10k' ? raw.map(x => D.species[i].n_proteins ? x / D.species[i].n_proteins * 1e4 : 0)
      : raw;
    return {sp: s, i, raw, v, tot}; });
  rows.sort((a, b) => (b.v[0] || 0) - (a.v[0] || 0));
  const W = 900, rh = 18, H = rows.length * rh + 30;
  const sv = svgEl('svg', {width: W, height: H, viewBox: `0 0 ${W} ${H}`});
  const ML = 200, MR = 30;
  const mx = Math.max(...rows.map(r => r.v.reduce((a, b) => a + b, 0))) || 1;
  const cstep = labelStep(rows.length, 70);
  rows.forEach((r, k) => { let x = ML; const y = k * rh + 6;
    r.v.forEach((v, q) => { const w = v / mx * (W - ML - MR);
      const rect = svgEl('rect', {x, y, width: Math.max(0, w), height: rh - 3, fill: ccol(clans[q])});
      rect.addEventListener('mousemove', e => showTip(e,
        `<b>${fullName(r.sp)}</b><br>${clans[q]}: ${r.raw[q]} genes` +
        `<br>${(r.raw[q] / r.tot * 100).toFixed(1)}% of ${fam}` +
        `<br>${(r.raw[q] / D.species[r.i].n_proteins * 1e4).toFixed(2)} per 10k ${D.size_unit || 'proteins'}`));
      rect.addEventListener('mouseleave', hideTip); sv.append(rect); x += w; });
    if (k % cstep) return;
    sv.append(svgSpeciesLabel(r.sp, {x: ML - 8, y: y + rh - 6, 'text-anchor': 'end',
      fill: traitColour(D.species[r.i].trait)})); });
  host.append(figure(sv, fam + '_clans'), el('div', {class: 'legend'},
    clans.map(c => el('span', {}, el('i', {class: 'sw', style: 'background:' + ccol(c)}), c)),
    TRAITS.map(t => el('span', {style: 'color:' + traitColour(t)}, '● ' + t))));

  if (TRAITS.length) {
    const agg = {}; TRAITS.forEach(t => agg[t] = clans.map(() => 0));
    rows.forEach(r => { const t = D.species[r.i].trait; if (!t) return;
      r.raw.forEach((v, q) => agg[t][q] += v); });
    const tp = el('div', {class: 'panel'});
    tp.append(el('h2', {}, `${fam} clans by ${D.trait_name}`), el('p', {class: 'hint'},
      'Totals summed across species, then the same table as within-group fractions. If the fractions are near identical, clan balance does not track the trait, which is a real negative result.'));
    const mkT = (title, fmtf) => { const t = el('table', {},
      el('thead', {}, el('tr', {}, [D.trait_name, ...clans, 'total'].map(h => el('th', {}, h)))),
      el('tbody', {}, TRAITS.map(tr => { const tot = agg[tr].reduce((a, b) => a + b, 0);
        return el('tr', {}, el('td', {}, tr), ...agg[tr].map(v => el('td', {}, fmtf(v, tot))),
          el('td', {}, tot)); })));
      return el('div', {}, el('h2', {style: 'font-size:13px;margin:10px 0 4px'}, title), t); };
    tp.append(mkT('counts', v => v), mkT('fraction within group', (v, t) => t ? (v / t).toFixed(3) : '–'));
    host.append(tp);
  }
  const tbl = el('div', {class: 'panel'});
  tbl.append(el('h2', {}, 'Per species'));
  makeTable(tbl, ['species', D.trait_name || 'trait', ...clans, 'total'],
    rows.map(r => [r.sp, D.species[r.i].trait ?? '–', ...r.raw, r.tot]),
    {sortBy: 2, desc: true, name: fam + '_clans'});
  host.append(tbl);
}

/* ---------------------------------------------------------------- TREES */
/* Gene trees are read as evidence about the METHOD comparison: when DToL,
   genome-guided and de novo proteomes each recover a gene, their copies
   should fall together; a clade built from one method alone is either a
   transcript-only gene, a fragment family, or a contaminant, and is worth
   finding. Everything here is computed from the embedded Newick. */
function parseNewick(txt) {
  txt = txt.replace(/\[[^\]]*\]/g, '').trim().replace(/;\s*$/, '');
  let i = 0;
  function node() {
    const n = {name: '', len: 0, kids: [], support: null};
    if (txt[i] === '(') { i++;
      for (;;) { const k = node(); k.parent = n; n.kids.push(k);
        if (txt[i] === ',') i++; else break; }
      if (txt[i] === ')') i++; }
    let s = i;
    if (txt[i] === "'") { i++; while (i < txt.length && txt[i] !== "'") i++; i++; }
    else while (i < txt.length && !',():;'.includes(txt[i])) i++;
    const label = txt.slice(s, i).trim().replace(/^'|'$/g, '');
    if (txt[i] === ':') { i++; s = i;
      while (i < txt.length && /[-\d.eE+]/.test(txt[i])) i++;
      n.len = parseFloat(txt.slice(s, i)) || 0; }
    if (!n.kids.length) n.name = label;
    else if (label !== '' && !isNaN(+label)) n.support = +label;   // FastTree local support
    return n; }
  return node();
}
const walk = r => { const out = [], st = [r];
  while (st.length) { const n = st.pop(); out.push(n); st.push(...n.kids); } return out; };
const leaves = r => walk(r).filter(n => !n.kids.length);
function pruneTree(root, keep) {
  const rec = n => { if (!n.kids.length) return keep.has(n.name) ? n : null;
    const k = n.kids.map(rec).filter(Boolean);
    if (!k.length) return null;
    if (k.length === 1) { k[0].len += n.len; return k[0]; }
    n.kids = k; return n; };
  const r = rec(root); if (r) r.len = 0; return r; }

/* ---- per-tip annotation ------------------------------------------------ */
function tipInfo(name) {
  const bar = name.indexOf('|');
  const sid = bar >= 0 ? name.slice(0, bar) : name;
  const prot = bar >= 0 ? name.slice(bar + 1) : '';
  const s = (D.samples || {})[sid] || {};
  const i = SPI[sid];
  const row = i !== undefined ? D.species[i] : {};
  return {sid, prot,
    species: s.species || row.species_name || sid,
    method: s.method || row.method || '',
    pipeline: s.pipeline || row.pipeline || '',
    clan: (D.tip_clan || {})[name] || '',
    status: prot ? ((D.tip_status || {})[name] || 'COMPLETE') : '',
    trait: row.trait ?? null};
}

/* ---- palettes: every ring gets its own, so no colour means two things --- */
const OKABE = ['#0072B2', '#D55E00', '#009E73', '#CC79A7', '#E69F00', '#56B4E9', '#F0E442', '#000000'];
const RING_FIXED = {
  method: {DToL: '#1d3557', genome_guided: '#457b9d', denovo: '#e63946'},
  clan: {CYP2: '#7b2cbf', CYP3: '#f4a261', CYP4: '#2a9d8f', MITO: '#8d6e63'},
  status: {COMPLETE: '#ced4da', PARTIAL: '#e9a13b', FRAGMENT: '#bc4749', MISSING_DOMAIN: '#343a40'}};
function ringPalette(field, values) {
  const vals = [...new Set(values.filter(v => v !== '' && v != null))].sort();
  const pal = {};
  if (CVD) { vals.forEach((v, i) => pal[v] = OKABE[i % OKABE.length]); return pal; }
  const fixed = RING_FIXED[field] || {};
  let h = field === 'species' ? 12 : 200, extra = 0;
  vals.forEach(v => {
    if (fixed[v]) pal[v] = fixed[v];
    else if (field === 'trait') pal[v] = traitColour(v);
    else { pal[v] = `hsl(${(h + extra * 137.5) % 360},55%,${46 + (extra % 3) * 8}%)`; extra++; }
  });
  return pal;
}
const METHOD_ABBR = {DToL: 'DToL', genome_guided: 'GG', denovo: 'dn'};
/* On a species tree the tip IS the species, so the common name belongs on the
   label; above a couple of hundred tips it doubles the label block and the tree
   itself disappears, so past that the binomial stands alone and the common name
   stays on hover. */
const TREE_COMMON_MAX = 200;
function tipText(t, mode, isSpecies, nTips) {
  if (isSpecies) return (nTips && nTips > TREE_COMMON_MAX) ? nice(t.name)
                                                           : fullName(t.name, 40);
  const i = t.info;
  if (mode === 'prot') return i.prot || t.name;
  if (mode === 'sid') return i.sid;
  if (mode === 'sp') return nice(i.species);
  return `${nice(i.species)} ${METHOD_ABBR[i.method] || i.method || ''}`.trim();
}
const RING_LABEL = {method: 'method', species: 'species', clan: 'CYP clan',
                    status: 'domain completeness', trait: () => D.trait_name || 'trait'};

/* ---- model: parse once, ladderise, index tips in drawing order --------- */
const TREE_CACHE = {};
function treeModel(key) {
  if (TREE_CACHE[key]) return TREE_CACHE[key];
  let root = parseNewick(D.trees[key]);
  const isSpecies = leaves(root).every(n => SPI[n.name] !== undefined);
  if (isSpecies) { const pr = pruneTree(root, new Set(SP)); if (pr) root = pr; }
  (function lad(n) { if (!n.kids.length) { n.nt = 1; return 1; }
    n.nt = n.kids.reduce((a, k) => a + lad(k), 0);
    n.kids.sort((a, b) => a.nt - b.nt); return n.nt; })(root);
  (function par(n, p) { n.parent = p; n.kids.forEach(k => par(k, n)); })(root, null);
  const lv = []; (function dfs(n) { if (!n.kids.length) lv.push(n); else n.kids.forEach(dfs); })(root);
  const all = []; (function po(n) { n.kids.forEach(po); all.push(n); })(root);
  (function dd(n, dep, dist) { n.dep = dep; n.dist = dist;
    n.kids.forEach(k => dd(k, dep + 1, dist + Math.max(0, k.len))); })(root, 0, 0);
  lv.forEach((n, i) => { n.idx = i; n.info = tipInfo(n.name); });
  // clade tips are a contiguous range in DFS order, so a clade is [lo, hi]
  all.forEach(n => { if (!n.kids.length) { n.lo = n.hi = n.idx; }
    else { n.lo = n.kids[0].lo; n.hi = n.kids[n.kids.length - 1].hi; } });
  let maxDep = 0, maxDist = 0;
  lv.forEach(n => { if (n.dep > maxDep) maxDep = n.dep; if (n.dist > maxDist) maxDist = n.dist; });
  const m = {key, root, lv, all, isSpecies, maxDep: maxDep || 1, maxDist: maxDist || 1,
             hasSupport: all.some(n => n.support != null),
             fields: ['method', 'species', 'clan', 'status', 'trait'].filter(f =>
               lv.some(t => t.info[f] !== '' && t.info[f] != null) &&
               !(f === 'trait' && D.trait_name === 'method'))};
  TREE_CACHE[key] = m; return m;
}

function cladeStats(m, n) {
  const tips = m.lv.slice(n.lo, n.hi + 1);
  const by = f => { const o = {}; tips.forEach(t => { const v = t.info[f] || '–'; o[v] = (o[v] || 0) + 1; }); return o; };
  return {n: tips.length, tips, method: by('method'), species: by('species'),
          clan: by('clan'), status: by('status')};
}

/* maximal clades whose every tip shares one value of `field`. For method,
   these are the clades only one annotation route produced. */
function exclusiveClades(m, field, minTips = 3) {
  const val = new Map();
  m.all.forEach(n => {
    if (!n.kids.length) { val.set(n, n.info[field] || '–'); return; }
    const vs = n.kids.map(k => val.get(k));
    val.set(n, vs.every(v => v !== null && v === vs[0]) ? vs[0] : null); });
  const out = [];
  m.all.forEach(n => {
    if (!n.kids.length) return;
    const v = val.get(n); if (v === null) return;
    const size = n.hi - n.lo + 1; if (size < minTips) return;
    if (n.parent && val.get(n.parent) === v) return;           // not maximal
    const sp = new Set(m.lv.slice(n.lo, n.hi + 1).map(t => t.info.species));
    out.push({node: n, value: v, n: size, species: sp.size, support: n.support}); });
  return out.sort((a, b) => b.n - a.n);
}

/* ---- geometry ---------------------------------------------------------- */
const TG = {S: 940, RT: 0.29, RW: 13, GAP: 2, TOP: 30, LEFT: 30, TW: 360, CW: 13};
function place(m, layout, useLen) {
  const nr = n => useLen ? n.dist / m.maxDist : (n.kids.length ? n.dep / m.maxDep : 1);
  const n = m.lv.length;
  if (layout === 'radial') {
    const span = 2 * Math.PI * 0.97, a0 = -Math.PI / 2 + Math.PI * 0.015;
    const step = span / Math.max(n, 1), RT = TG.S * TG.RT;
    m.lv.forEach(t => t.a = a0 + t.idx * step);
    m.all.forEach(x => { if (x.kids.length) x.a = (x.kids[0].a + x.kids[x.kids.length - 1].a) / 2; });
    m.all.forEach(x => { x.r = RT * nr(x); x.x = x.r * Math.cos(x.a); x.y = x.r * Math.sin(x.a); });
    return {layout, step, a0, span, RT, unit: RT / (useLen ? m.maxDist : 1)};
  }
  const H = Math.max(640, Math.min(n * 11, 2400)), dy = (H - 2 * TG.TOP) / Math.max(n - 1, 1);
  m.lv.forEach(t => t.y = TG.TOP + t.idx * dy);
  m.all.forEach(x => { if (x.kids.length) x.y = (x.kids[0].y + x.kids[x.kids.length - 1].y) / 2; });
  m.all.forEach(x => x.x = TG.LEFT + TG.TW * nr(x));
  return {layout, dy, H, unit: TG.TW / (useLen ? m.maxDist : 1)};
}
function branchPath(m, g) {
  let d = '';
  const f = v => v.toFixed(1);
  m.all.forEach(x => {
    if (!x.kids.length) return;
    if (g.layout === 'radial') {
      x.kids.forEach(k => { d += `M${f(x.r * Math.cos(k.a))},${f(x.r * Math.sin(k.a))}L${f(k.x)},${f(k.y)}`; });
      const a0 = x.kids[0].a, a1 = x.kids[x.kids.length - 1].a;
      if (x.r > 0 && a1 > a0) d += `M${f(x.r * Math.cos(a0))},${f(x.r * Math.sin(a0))}A${f(x.r)},${f(x.r)} 0 ${a1 - a0 > Math.PI ? 1 : 0},1 ${f(x.r * Math.cos(a1))},${f(x.r * Math.sin(a1))}`;
    } else {
      d += `M${f(x.x)},${f(x.kids[0].y)}L${f(x.x)},${f(x.kids[x.kids.length - 1].y)}`;
      x.kids.forEach(k => { d += `M${f(x.x)},${f(k.y)}L${f(k.x)},${f(k.y)}`; });
    }
  });
  return d;
}
function ringPaths(m, g, field, k, pal) {
  const byCol = {}, f = v => v.toFixed(1);
  m.lv.forEach(t => {
    const v = field === '__value' ? t.__v : t.info[field];
    const col = field === '__value' ? t.__c : (pal[v] || '#f1f3f5');
    let d;
    if (g.layout === 'radial') {
      const r0 = g.RT + 8 + k * (TG.RW + TG.GAP), r1 = r0 + TG.RW, h = g.step / 2 * 0.96;
      const p = (r, a) => `${f(r * Math.cos(a))},${f(r * Math.sin(a))}`;
      d = `M${p(r0, t.a - h)}A${f(r0)},${f(r0)} 0 0,1 ${p(r0, t.a + h)}L${p(r1, t.a + h)}A${f(r1)},${f(r1)} 0 0,0 ${p(r1, t.a - h)}Z`;
    } else {
      const x0 = TG.LEFT + TG.TW + 8 + k * (TG.CW + TG.GAP), h = Math.max(g.dy / 2 * 0.96, .3);
      d = `M${f(x0)},${f(t.y - h)}h${TG.CW}v${f(2 * h)}h${-TG.CW}Z`;
    }
    (byCol[col] = byCol[col] || []).push(d);
  });
  return byCol;
}

/* ---- the tab ----------------------------------------------------------- */
function tabTrees(host) {
  const keys = Object.keys(D.trees || {});
  const R0 = (typeof ROUTE === 'object' && ROUTE) || {};
  if (!keys.length) { host.append(el('div', {class: 'panel'}, el('h2', {}, 'Trees'),
    el('p', {class: 'hint'}, 'No Newick files embedded. Run `trees`, and pass --species-tree to `dashboard` for a species tree.'))); return; }
  const gene = keys.filter(k => k !== 'species');
  const S = {key: keys.includes(R0.tree) ? R0.tree : gene.includes('CYP') ? 'CYP' : keys[0], layout: 'radial', useLen: false,
             rings: new Set(), labels: 'auto', txt: 'sm', support: false, famRing: '',
             excl: 'method', sel: null, tip: null, k: 1, tx: 0, ty: 0, query: R0.q || ''};

  const p = el('div', {class: 'panel'});
  p.append(el('h2', {}, 'Trees'), el('p', {class: 'hint'},
    'Scroll to zoom, drag to pan, double-click to reset. Click a branch point to inspect that clade; click a tip for its record and a link to the gene page. ' +
    'Where several annotation methods were run, copies of one gene should sit together; a clade built from a single method is a method-specific gene model, a fragment family or a contaminant, and the table below lists them.'));
  const mk = (l, n) => el('div', {class: 'ctl'}, el('label', {}, l), n);
  const selTree = el('select', {style: 'max-width:420px'}, ...keys.map(k => el('option', {value: k}, treeLabel(k))));
  selTree.value = S.key;
  const selLay = el('select', {}, el('option', {value: 'radial'}, 'radial'), el('option', {value: 'rect'}, 'rectangular'));
  const selLen = el('select', {}, el('option', {value: 'clad'}, 'cladogram'), el('option', {value: 'phylo'}, 'phylogram'));
  const selLab = el('select', {}, el('option', {value: 'auto'}, 'labels: when readable'),
    el('option', {value: 'on'}, 'labels: always'), el('option', {value: 'off'}, 'labels: off'));
  const selTxt = el('select', {}, el('option', {value: 'sm'}, 'species + method'),
    el('option', {value: 'prot'}, 'protein ID'), el('option', {value: 'sid'}, 'sample'),
    el('option', {value: 'sp'}, 'species'));
  const search = el('input', {type: 'search', placeholder: 'find protein, sample or species', value: S.query});
  const famSel = el('select', {}, el('option', {value: ''}, '(none)'), ...VARFAM.map(f => el('option', {value: f}, f)));
  const chips = el('div', {class: 'chips'});
  const supChip = el('span', {class: 'chip'}, 'support');
  p.append(el('div', {class: 'row'}, mk('tree', selTree), mk('layout', selLay), mk('branches', selLen),
    mk('labels', selLab), mk('label text', selTxt), mk('search', search),
    mk('value ring (species tree)', famSel)),
    el('div', {class: 'row'}, mk('rings', chips), mk('nodes', supChip)));
  const zbar = el('div', {class: 'row', style: 'margin:0 0 6px'});
  const zb = (t, fn) => el('button', {class: 'act', onclick: fn}, t);
  zbar.append(zb('+', () => zoomBy(1.5)), zb('−', () => zoomBy(1 / 1.5)), zb('reset', () => resetView()),
              zb('zoom to selection', () => S.sel && fitClade(S.sel)), zb('clear selection', () => { S.sel = null; S.tip = null; paint(); }));
  const stage = el('div', {style: 'display:grid;grid-template-columns:minmax(0,1fr) 300px;gap:14px;align-items:start'});
  const figHost = el('div'), side = el('div', {class: 'panel', style: 'margin:0;max-height:900px;overflow:auto'});
  stage.append(figHost, side);
  const legend = el('div', {class: 'legend', style: 'margin-top:8px'});
  const exHost = el('div', {class: 'panel'});
  p.append(zbar, stage, legend); host.append(p, exHost);

  let m, g, svg, vp, gBranch, gSel, gRing, gNode, gHit, gLab, gBar;
  const NS = 'http://www.w3.org/2000/svg';

  function build() {
    if (D.trees[S.key] == null) {      // a large tree lives in its own data file
      figHost.innerHTML = ''; side.innerHTML = ''; exHost.innerHTML = '';
      figHost.append(el('p', {class: 'empty'}, `Loading ${treeLabel(S.key)}…`));
      const want = S.key;
      DEF.load('tree', want).then(o => { D.trees[want] = o.nwk;
        Object.assign(D.tip_clan = D.tip_clan || {}, o.clan || {});
        Object.assign(D.tip_status = D.tip_status || {}, o.status || {});
        if (S.key === want) build(); })
        .catch(() => { figHost.innerHTML = ''; figHost.append(el('p', {class: 'empty'},
          'Tree data not found. Keep the data/ folder beside index.html.')); });
      return;
    }
    m = treeModel(S.key);
    if (!S.rings.size) (m.isSpecies ? ['trait'] : ['method', 'clan']).forEach(f => m.fields.includes(f) && S.rings.add(f));
    chips.innerHTML = '';
    m.fields.forEach(f => chips.append(el('span', {class: 'chip' + (S.rings.has(f) ? ' on' : ''),
      onclick: e => { S.rings.has(f) ? S.rings.delete(f) : S.rings.add(f); e.target.classList.toggle('on'); render(); }},
      typeof RING_LABEL[f] === 'function' ? RING_LABEL[f]() : RING_LABEL[f])));
    supChip.style.display = m.hasSupport ? '' : 'none';
    famSel.disabled = !m.isSpecies;
    S.sel = null; S.tip = null; resetView(false); render(); renderExclusive();
    focusQuery();
  }

  function render() {
    hideTip();
    g = place(m, S.layout, S.useLen);
    const rings = [...S.rings].filter(f => m.fields.includes(f));
    if (S.famRing && m.isSpecies) rings.push('__value');
    if (S.famRing && m.isSpecies) {
      const kind = D.per10k ? 'per10k' : 'raw', z = zcol(kind, S.famRing), v = col(kind, S.famRing);
      m.lv.forEach(t => { const i = SPI[t.name]; t.__v = i === undefined ? null : v[i];
        t.__c = i === undefined ? '#f1f3f5' : rdbu((z[i] + 2.5) / 5); });
    }
    // Reserve room for tip labels when they will be drawn at the starting
    // zoom, so the outermost text is never cut off at the canvas edge.
    const ringsOuter = S.layout === 'radial' ? g.RT + 8 + rings.length * (TG.RW + TG.GAP)
                                            : TG.LEFT + TG.TW + 8 + rings.length * (TG.CW + TG.GAP);
    const sp0 = S.layout === 'radial' ? (ringsOuter + 4) * g.step : g.dy;
    const labelsAt1 = S.labels === 'on' || (S.labels === 'auto' && sp0 >= 6);
    const fs0 = Math.min(11, Math.max(4, sp0 * .82));
    const maxChars = labelsAt1 ? Math.max(...m.lv.map(t => tipText(t, S.txt, m.isSpecies, m.lv.length).length)) : 0;
    const labW = labelsAt1 ? maxChars * fs0 * .56 + 10 : 0;
    let vb;
    if (S.layout === 'radial') { const half = Math.max(TG.S / 2, ringsOuter + 4 + labW + 12); vb = [-half, -half, 2 * half, 2 * half]; }
    else vb = [0, 0, Math.max(TG.S, ringsOuter + 4 + labW + 12), g.H];
    g.vb = vb;
    const W = TG.S, H = Math.round(vb[3] * TG.S / vb[2]);
    svg = svgEl('svg', {width: W, height: H, viewBox: vb.join(' '),
      style: 'background:#fff;cursor:grab;touch-action:none;max-width:100%;height:auto'});
    vp = svgEl('g'); svg.append(vp);
    gBranch = svgEl('path', {d: branchPath(m, g), fill: 'none', stroke: '#5c636a',
      'stroke-width': m.lv.length > 800 ? .45 : .9, 'vector-effect': 'non-scaling-stroke'});
    gSel = svgEl('path', {fill: 'none', stroke: '#e63946', 'stroke-width': 2.2, 'vector-effect': 'non-scaling-stroke'});
    gRing = svgEl('g'); gNode = svgEl('g'); gHit = svgEl('g'); gLab = svgEl('g'); gBar = svgEl('g');
    vp.append(gBranch, gSel, gRing, gNode, gHit, gLab);
    svg.append(gBar);
    const pals = {};
    rings.forEach((f, k) => {
      if (f !== '__value') pals[f] = ringPalette(f, m.lv.map(t => t.info[f]));
      const byCol = ringPaths(m, g, f, k, pals[f] || {});
      Object.entries(byCol).forEach(([c, ds]) => gRing.append(svgEl('path', {d: ds.join(''), fill: c, stroke: 'none'})));
    });
    g.ringsOuter = S.layout === 'radial' ? g.RT + 8 + rings.length * (TG.RW + TG.GAP)
                                        : TG.LEFT + TG.TW + 8 + rings.length * (TG.CW + TG.GAP);
    if (S.support && m.hasSupport) m.all.forEach(x => {
      if (x.support == null || !x.kids.length) return;
      const c = x.support >= .9 ? '#1d3557' : x.support >= .7 ? '#8d99ae' : '#e63946';
      gNode.append(svgEl('circle', {cx: x.x.toFixed(1), cy: x.y.toFixed(1), r: 1.8, fill: c}));
    });
    if (S.useLen) {
      const raw = m.maxDist / 5, pow = Math.pow(10, Math.floor(Math.log10(raw)));
      const L = [1, 2, 5, 10].map(q => q * pow).find(q => q >= raw) || raw, px = L * g.unit;
      const x0 = g.vb[0] + 20, y0 = g.vb[1] + g.vb[3] - 16;
      gBar.append(svgEl('path', {d: `M${x0},${y0}h${px.toFixed(1)}M${x0},${y0 - 4}v8M${(x0 + px).toFixed(1)},${y0 - 4}v8`, stroke: '#333', 'stroke-width': 1.2, fill: 'none'}));
      const tt = svgEl('text', {x: x0 + px / 2, y: y0 - 7, 'text-anchor': 'middle', 'font-size': 10, fill: '#333'});
      tt.textContent = `${+L.toPrecision(2)} substitutions/site`; gBar.append(tt);
    }
    wire(); paint(); labels();
    figHost.innerHTML = '';
    figHost.append(figure(svg, `tree_${S.key}_${S.layout}`));
    legend.innerHTML = '';
    rings.forEach(f => {
      if (f === '__value') { legend.append(colourbar(rdbu, -2.5, 2.5, `${S.famRing}, z-score`, 220, 10)); return; }
      const pal = pals[f], lab = typeof RING_LABEL[f] === 'function' ? RING_LABEL[f]() : RING_LABEL[f];
      const cnt = {}; m.lv.forEach(t => { const v = t.info[f]; if (v) cnt[v] = (cnt[v] || 0) + 1; });
      legend.append(el('span', {style: 'font-weight:600'}, lab + ':'),
        ...Object.keys(pal).slice(0, 14).map(v => el('span', {},
          el('i', {class: 'sw', style: 'background:' + pal[v]}), `${nice(v)} (${cnt[v] || 0})`)));
      if (Object.keys(pal).length > 14) legend.append(el('span', {class: 'hint'}, `+${Object.keys(pal).length - 14} more`));
    });
    if (S.support && m.hasSupport) legend.append(el('span', {style: 'font-weight:600'}, 'support:'),
      el('span', {}, el('i', {class: 'sw', style: 'background:#1d3557'}), '≥0.9'),
      el('span', {}, el('i', {class: 'sw', style: 'background:#8d99ae'}), '0.7–0.9'),
      el('span', {}, el('i', {class: 'sw', style: 'background:#e63946'}), '<0.7'));
    legend.append(el('span', {class: 'hint'}, `${m.lv.length.toLocaleString()} tips`));
  }

  /* labels are rebuilt on zoom so they appear exactly when they become legible */
  function labels() {
    if (!gLab) return;
    gLab.innerHTML = '';
    const spacing = (S.layout === 'radial' ? (g.ringsOuter + 4) * g.step : g.dy) * S.k;
    if (S.labels === 'off' || (S.labels === 'auto' && spacing < 6)) return;
    const fs = Math.min(11, Math.max(4, spacing * .82)) / S.k;
    const q = S.query.toLowerCase();
    m.lv.forEach(t => {
      const txt = tipText(t, S.txt, m.isSpecies, m.lv.length);
      const hit = q && t.name.toLowerCase().includes(q);
      let node;
      if (S.layout === 'radial') {
        const r = g.ringsOuter + 4, deg = t.a * 180 / Math.PI, flip = deg > 90 && deg < 270;
        const x = r * Math.cos(t.a), y = r * Math.sin(t.a);
        node = svgEl('text', {x: x.toFixed(1), y: y.toFixed(1), 'font-size': fs.toFixed(2),
          'text-anchor': flip ? 'end' : 'start', 'dominant-baseline': 'middle',
          transform: `rotate(${(flip ? deg + 180 : deg).toFixed(2)} ${x.toFixed(1)} ${y.toFixed(1)})`,
          fill: hit ? '#e63946' : '#212529', 'font-weight': hit ? 700 : 400});
      } else {
        node = svgEl('text', {x: (g.ringsOuter + 4).toFixed(1), y: t.y.toFixed(1), 'font-size': fs.toFixed(2),
          'dominant-baseline': 'middle', fill: hit ? '#e63946' : '#212529', 'font-weight': hit ? 700 : 400});
      }
      node.textContent = txt; gLab.append(node);
    });
  }

  /* selection highlight, search markers and the side panel */
  function paint() {
    if (!gSel) return;
    let d = '';
    if (S.sel) {
      const f = v => v.toFixed(1), sub = [];
      (function rec(x) { sub.push(x); x.kids.forEach(rec); })(S.sel);
      sub.forEach(x => { if (!x.kids.length) return;
        if (g.layout === 'radial') {
          x.kids.forEach(k => { d += `M${f(x.r * Math.cos(k.a))},${f(x.r * Math.sin(k.a))}L${f(k.x)},${f(k.y)}`; });
          const a0 = x.kids[0].a, a1 = x.kids[x.kids.length - 1].a;
          if (x.r > 0) d += `M${f(x.r * Math.cos(a0))},${f(x.r * Math.sin(a0))}A${f(x.r)},${f(x.r)} 0 ${a1 - a0 > Math.PI ? 1 : 0},1 ${f(x.r * Math.cos(a1))},${f(x.r * Math.sin(a1))}`;
        } else {
          d += `M${f(x.x)},${f(x.kids[0].y)}L${f(x.x)},${f(x.kids[x.kids.length - 1].y)}`;
          x.kids.forEach(k => { d += `M${f(x.x)},${f(k.y)}L${f(k.x)},${f(k.y)}`; }); } });
    }
    gSel.setAttribute('d', d);
    gHit.innerHTML = '';
    // Search the vernacular name too: "zebrafish" is how most readers,
    // including a regulator, will look for Danio rerio.
    const q = S.query.toLowerCase(), hits = q ? m.lv.filter(t => t.name.toLowerCase().includes(q) ||
      String(t.info.species).toLowerCase().includes(q) ||
      commonName(t.name).toLowerCase().includes(q)) : [];
    hits.slice(0, 2000).forEach(t => {
      const r = S.layout === 'radial' ? g.ringsOuter + 2 : 0;
      const x = S.layout === 'radial' ? r * Math.cos(t.a) : g.ringsOuter + 2, y = S.layout === 'radial' ? r * Math.sin(t.a) : t.y;
      gHit.append(svgEl('circle', {cx: x.toFixed(1), cy: y.toFixed(1), r: (3.2 / S.k).toFixed(2), fill: '#ffb703', stroke: '#333', 'stroke-width': (.5 / S.k).toFixed(2)}));
    });
    if (S.tip) { const t = S.tip;
      const x = S.layout === 'radial' ? t.x : t.x, y = t.y;
      gHit.append(svgEl('circle', {cx: x.toFixed(1), cy: y.toFixed(1), r: (5 / S.k).toFixed(2), fill: 'none', stroke: '#e63946', 'stroke-width': (2 / S.k).toFixed(2)})); }
    renderSide(hits);
  }

  function renderSide(hits) {
    side.innerHTML = '';
    if (S.query) side.append(el('div', {class: 'note'}, `${hits.length} match${hits.length === 1 ? '' : 'es'} for “${S.query}”`),
      ...(hits.length ? [el('button', {class: 'act', onclick: () => fitTips(hits)}, 'zoom to matches')] : []));
    if (S.tip) {
      const i = S.tip.info;
      side.append(el('h2', {}, 'Tip'), el('table', {}, el('tbody', {},
        ...[['sample', i.sid], ['species', nice(i.species)], ['method', i.method], ['pipeline', i.pipeline],
            ['protein', i.prot], ['CYP clan', i.clan], ['domains', i.status],
            ['branch length', S.tip.len.toPrecision(3)]]
          .filter(r => r[1]).map(r => el('tr', {}, el('td', {style: 'color:#6b7280'}, r[0]),
            el('td', {class: r[0] === 'protein' ? 'mono' : ''}, String(r[1])))))));
      const sidx = SPKEY[i.sid] !== undefined ? i.sid : (SPKEY[S.tip.name] !== undefined ? S.tip.name : null);
      if (sidx) side.append(el('div', {class: 'row', style: 'margin-top:8px'},
        i.prot ? link(geneHref(sidx, i.prot), 'gene page ›') : null, link(spHref(sidx), 'species page ›')));
    }
    if (!S.sel) {
      if (!S.tip && !S.query) side.append(el('h2', {}, 'Clade'), el('p', {class: 'hint'},
        'Click a branch point to see what that clade is made of: how many tips each method contributed, which species, and which clans.'));
      return;
    }
    const st = cladeStats(m, S.sel);
    side.append(el('h2', {}, `Clade: ${st.n} tips`),
      el('div', {class: 'hint'}, S.sel.support != null ? `support ${S.sel.support}` : ''));
    const bars = (title, obj, field) => {
      const pal = ringPalette(field, Object.keys(obj)), tot = Object.values(obj).reduce((a, b) => a + b, 0);
      const box = el('div', {style: 'margin:8px 0'}, el('div', {style: 'font-size:12px;font-weight:600'}, title));
      Object.entries(obj).sort((a, b) => b[1] - a[1]).slice(0, 12).forEach(([k, v]) => box.append(
        el('div', {style: 'display:flex;align-items:center;gap:6px;font-size:11.5px;margin:2px 0'},
          el('span', {style: 'width:92px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap'}, nice(k)),
          el('span', {style: `height:10px;width:${(v / tot * 120).toFixed(0)}px;background:${pal[k] || '#adb5bd'};border-radius:2px`}),
          el('span', {style: 'color:#6b7280'}, String(v)))));
      return box; };
    side.append(bars('by method', st.method, 'method'));
    if (Object.keys(st.clan).some(k => k !== '–')) side.append(bars('by clan', st.clan, 'clan'));
    side.append(bars('by species', st.species, 'species'));
    if (Object.keys(st.status).some(k => k !== 'COMPLETE' && k !== '–')) side.append(bars('domain completeness', st.status, 'status'));
    const methods = Object.keys(st.method).filter(k => k !== '–');
    if (methods.length === 1 && st.n >= 2) side.append(el('div', {class: 'note', style: 'border-left-color:var(--warm)'},
      `Every tip here comes from ${methods[0]}. If the other methods sampled these species, this is a ${methods[0]}-only gene, a fragment family, or a contaminant.`));
    side.append(el('div', {class: 'row'},
      el('button', {class: 'act', onclick: () => fitClade(S.sel)}, 'zoom to clade'),
      el('button', {class: 'act', onclick: () => dl(`clade_${S.key}_${st.n}tips.tsv`,
        ['sample\tspecies\tmethod\tprotein\tclan\tdomains', ...st.tips.map(t =>
          [t.info.sid, t.info.species, t.info.method, t.info.prot, t.info.clan, t.info.status].join('\t'))].join('\n'))}, 'download tips')));
  }

  function renderExclusive() {
    exHost.innerHTML = '';
    if (m.isSpecies) return;
    const fields = ['method', 'species', 'clan'].filter(f => m.fields.includes(f));
    if (!fields.length) return;
    const sel = el('select', {}, ...fields.map(f => el('option', {value: f}, f)));
    sel.value = fields.includes(S.excl) ? S.excl : fields[0];
    exHost.append(el('h2', {}, 'Single-source clades'), el('p', {class: 'hint'},
      'Maximal clades of three or more tips in which every tip shares one value. By method, these are the genes only one annotation route recovered, which is where fragments and contaminants hide. By species, they are lineage-specific expansions.'),
      el('div', {class: 'row'}, el('div', {class: 'ctl'}, el('label', {}, 'exclusive by'), sel)));
    const body = el('div'); exHost.append(body);
    const draw = () => { S.excl = sel.value; body.innerHTML = '';
      const ex = exclusiveClades(m, S.excl, 3);
      const tot = {}; ex.forEach(e => { tot[e.value] = (tot[e.value] || 0) + e.n; });
      body.append(el('div', {class: 'hint'}, ex.length ? `${ex.length} clades; tips in them by ${S.excl}: ` +
        Object.entries(tot).map(([k, v]) => `${nice(k)} ${v}`).join(', ') : 'none'));
      if (!ex.length) return;
      const wrap = el('div', {class: 'scroll', style: 'max-height:360px'});
      wrap.append(el('table', {}, el('thead', {}, el('tr', {},
        [S.excl, 'tips', 'species', 'support', ''].map(h => el('th', {}, h)))),
        el('tbody', {}, ex.slice(0, 400).map(e => el('tr', {},
          el('td', {}, nice(e.value)), el('td', {}, e.n), el('td', {}, e.species),
          el('td', {}, e.support == null ? '–' : e.support),
          el('td', {}, el('button', {class: 'act', onclick: () => { S.sel = e.node; S.tip = null; paint(); fitClade(e.node); }}, 'show')))))));
      body.append(wrap);
      body.append(el('div', {class: 'row'}, el('button', {class: 'act', onclick: () => dl(`single_source_${S.key}_${S.excl}.tsv`,
        [`${S.excl}\ttips\tspecies\tsupport`, ...ex.map(e => [e.value, e.n, e.species, e.support ?? ''].join('\t'))].join('\n'))}, 'download TSV')));
    };
    sel.onchange = draw; draw();
  }

  /* ---- zoom and pan -------------------------------------------------- */
  function applyView() { if (vp) vp.setAttribute('transform', `translate(${S.tx.toFixed(2)},${S.ty.toFixed(2)}) scale(${S.k.toFixed(4)})`); }
  let labTimer = null;
  function afterZoom() { applyView(); clearTimeout(labTimer); labTimer = setTimeout(() => { labels(); paint(); }, 110); }
  function resetView(redraw = true) { S.k = 1; S.tx = 0; S.ty = 0; if (redraw) { applyView(); labels(); paint(); } }
  function toModel(e) {
    const b = svg.getBoundingClientRect(), vb = svg.getAttribute('viewBox').split(' ').map(Number);
    const sx = vb[2] / (b.width || vb[2]), sy = vb[3] / (b.height || vb[3]);
    const vx = vb[0] + (e.clientX - b.left) * sx, vy = vb[1] + (e.clientY - b.top) * sy;
    return {vx, vy, x: (vx - S.tx) / S.k, y: (vy - S.ty) / S.k};
  }
  function zoomBy(f, at) {
    const c = at || {vx: g.vb[0] + g.vb[2] / 2, vy: g.vb[1] + g.vb[3] / 2};
    const k2 = Math.max(.5, Math.min(60, S.k * f)), r = k2 / S.k;
    S.tx = c.vx - (c.vx - S.tx) * r; S.ty = c.vy - (c.vy - S.ty) * r; S.k = k2; afterZoom();
  }
  function fitBox(x0, y0, x1, y1) {
    const W = g.vb[2], H = g.vb[3], pad = W * .06;
    const k = Math.max(.5, Math.min(60, Math.min((W - pad) / Math.max(x1 - x0, 1), (H - pad) / Math.max(y1 - y0, 1))));
    const cx = (x0 + x1) / 2, cy = (y0 + y1) / 2;
    const vcx = g.vb[0] + W / 2, vcy = g.vb[1] + H / 2;
    S.k = k; S.tx = vcx - cx * k; S.ty = vcy - cy * k; afterZoom();
  }
  function fitTips(tips) {
    const xs = [], ys = [];
    tips.forEach(t => { if (S.layout === 'radial') { const r = g.ringsOuter; xs.push(r * Math.cos(t.a), t.x); ys.push(r * Math.sin(t.a), t.y); }
      else { xs.push(t.x, g.ringsOuter + 120); ys.push(t.y); } });
    fitBox(Math.min(...xs), Math.min(...ys), Math.max(...xs), Math.max(...ys));
  }
  /* Arriving from a gene page with ?q=<protein>: frame the match and, when
     there is exactly one, open its record. */
  function focusQuery() {
    if (!S.query || !m) return;
    const q = S.query.toLowerCase(), hits = m.lv.filter(t => t.name.toLowerCase().includes(q));
    if (hits.length === 1) {       // the tip in context: its smallest clade of 25+ tips
      S.tip = hits[0]; paint();
      let n = hits[0]; const want = Math.min(25, Math.ceil(m.lv.length / 4));
      while (n.parent && n.nt < want) n = n.parent;
      if (n.parent || n.nt < m.lv.length) fitClade(n);
    } else if (hits.length && hits.length < m.lv.length) fitTips(hits);
  }
  tabTrees._focus = () => focusQuery();
  function fitClade(n) { const tips = m.lv.slice(n.lo, n.hi + 1);
    fitTips(tips.length > 3000 ? tips.filter((_, i) => i % Math.ceil(tips.length / 3000) === 0) : tips.concat([n])); }

  /* hit testing is geometric: one listener, no per-tip DOM nodes */
  function tipAt(pt) {
    if (S.layout === 'radial') {
      const r = Math.hypot(pt.x, pt.y);
      if (r < g.RT * .6 || r > g.ringsOuter + 260) return null;
      let a = Math.atan2(pt.y, pt.x); while (a < g.a0) a += 2 * Math.PI;
      const i = Math.round((a - g.a0) / g.step);
      return (i >= 0 && i < m.lv.length) ? m.lv[i] : null;
    }
    if (pt.x < TG.LEFT + TG.TW * .5) return null;
    const i = Math.round((pt.y - TG.TOP) / g.dy);
    return (i >= 0 && i < m.lv.length) ? m.lv[i] : null;
  }
  function nodeAt(pt) {
    let best = null, bd = 9 / S.k;
    m.all.forEach(x => { if (!x.kids.length) return;
      const d = Math.hypot(x.x - pt.x, x.y - pt.y); if (d < bd) { bd = d; best = x; } });
    return best;
  }
  function wire() {
    let drag = null;
    svg.addEventListener('wheel', e => { e.preventDefault(); zoomBy(e.deltaY < 0 ? 1.2 : 1 / 1.2, toModel(e)); }, {passive: false});
    svg.addEventListener('mousedown', e => { drag = {x: e.clientX, y: e.clientY, tx: S.tx, ty: S.ty, moved: false};
      svg.style.cursor = 'grabbing'; });
    window.addEventListener('mouseup', () => { if (svg) svg.style.cursor = 'grab'; setTimeout(() => drag = null, 0); });
    svg.addEventListener('mousemove', e => {
      if (drag && (e.buttons & 1)) {
        const b = svg.getBoundingClientRect(), vb = svg.getAttribute('viewBox').split(' ').map(Number);
        const s = vb[2] / (b.width || vb[2]);
        const dx = (e.clientX - drag.x) * s, dy = (e.clientY - drag.y) * s;
        if (Math.abs(dx) + Math.abs(dy) > 2) drag.moved = true;
        S.tx = drag.tx + dx; S.ty = drag.ty + dy; applyView(); hideTip(); return;
      }
      const pt = toModel(e), n = nodeAt(pt);
      if (n) { const st = cladeStats(m, n);
        showTip(e, `<b>clade, ${st.n} tips</b>` + (n.support != null ? `<br>support ${n.support}` : '') +
          '<br>' + Object.entries(st.method).map(([k, v]) => `${k} ${v}`).join(' · ') + '<br><i>click to inspect</i>'); return; }
      const t = tipAt(pt);
      if (t) { const i = t.info;
        showTip(e, `<b>${nice(i.species)}</b>${i.method ? ' · ' + i.method : ''}` +
          (i.prot ? `<br><span class="mono">${i.prot}</span>` : '') +
          (i.clan ? `<br>clan ${i.clan}` : '') + (i.status && i.status !== 'COMPLETE' ? `<br>${i.status}` : '')); }
      else hideTip();
    });
    svg.addEventListener('mouseleave', hideTip);
    svg.addEventListener('click', e => { if (drag && drag.moved) return;
      const pt = toModel(e), n = nodeAt(pt);
      if (n) { S.sel = n; S.tip = null; } else { const t = tipAt(pt); if (t) { S.tip = t; } else { S.sel = null; S.tip = null; } }
      paint(); });
    svg.addEventListener('dblclick', e => { e.preventDefault(); resetView(); });
    applyView();
  }

  selTree.onchange = () => { S.key = selTree.value; S.rings = new Set(); S.famRing = ''; famSel.value = ''; build(); };
  /* opened from a gene or family page: #/trees/<key>?q=<protein> */
  tabTrees._open = (key, q) => { if (key && keys.includes(key) && key !== S.key) {
      S.key = key; selTree.value = key; S.rings = new Set(); S.famRing = ''; famSel.value = ''; S.query = q || ''; search.value = S.query; build(); }
    else if (q != null && q !== S.query) { S.query = q; search.value = q; if (m) { paint(); labels(); focusQuery(); } } };
  selLay.onchange = () => { S.layout = selLay.value; resetView(false); render(); };
  selLen.onchange = () => { S.useLen = selLen.value === 'phylo'; render(); };
  selLab.onchange = () => { S.labels = selLab.value; render(); };
  selTxt.onchange = () => { S.txt = selTxt.value; render(); };
  famSel.onchange = () => { S.famRing = famSel.value; render(); };
  supChip.onclick = e => { S.support = !S.support; e.target.classList.toggle('on'); render(); };
  let qTimer = null;
  search.oninput = e => { clearTimeout(qTimer); qTimer = setTimeout(() => { S.query = e.target.value.trim(); paint(); labels(); }, 150); };
  build();
  /* exposed for the headless tests */
  tabTrees._state = S; tabTrees._model = () => m; tabTrees._geom = () => g;
  tabTrees._tipAt = tipAt; tabTrees._nodeAt = nodeAt;
  tabTrees._setLayout = l => { S.layout = l; resetView(false); render(); };
  tabTrees._svg = () => svg; tabTrees._select = n => { S.sel = n; paint(); };
  tabTrees._zoomBy = zoomBy; tabTrees._fitClade = fitClade;
}

/* ----------------------------------------------------------- DATA TABLES */
function tabData(host) {
  const p = el('div', {class: 'panel'});
  const tsv = kind => { const m = matrix(kind);
    return [['species', D.trait_name || 'group', 'proteins', ...FAM].join('\t'),
      ...SP.map((s, i) => [s, D.species[i].trait ?? '', D.species[i].n_proteins,
        ...FAM.map(f => m[f] ? m[f][i] : 0)].join('\t'))].join('\n'); };
  p.append(el('h2', {}, 'Copy numbers'), el('div', {class: 'row'},
    el('button', {class: 'act', onclick: () => dl('defensome_counts.tsv', tsv('raw'))}, 'copies per species (TSV)'),
    ...(D.per10k ? [el('button', {class: 'act', onclick: () => dl('defensome_per10k.tsv', tsv('per10k'))},
      `per 10k ${D.size_unit || 'proteins'} (TSV)`)] : []),
    el('span', {class: 'count', style: 'margin:0'}, `${SP.length.toLocaleString()} species in view`)));
  host.append(p);

  const fz = SITE.family_fasta || {};
  if (Object.keys(fz).length) {
    const q = el('div', {class: 'panel'});
    q.append(el('h2', {}, 'Sequences'), el('p', {class: 'hint'},
      'Every called protein, one gzipped FASTA per family; headers are sample|protein. Per-species FASTA is on each species page, per-gene on each gene page.'));
    dataGrid(q, [
      {label: 'family', get: r => r[0], cell: r => link(famHref(r[0]), r[0])},
      {label: 'file', get: r => r[0] + '.faa.gz', cell: r => el('a', {href: `data/fasta/${r[0]}.faa.gz`,
        download: r[0] + '.faa.gz', class: 'lk mono'}, r[0] + '.faa.gz')},
      {label: 'size', num: true, get: r => r[1], cell: r => kb(r[1])}],
      Object.entries(fz).filter(([, b]) => b > 40), {sortBy: 0, pageSize: 100, name: 'fasta_files',
        unit: 'file', units: 'files'});
    host.append(q);
  }
}

/* -------------------------------------------------------------- MAP TAB */
function tabMap(host) {
  const p = el('div', {class: 'panel'});
  p.append(el('h2', {}, 'Defensome map'), el('p', {class: 'hint'},
    'The rules that produced everything else. rule ALL means every listed domain must occur in the same protein; ANY means at least one. min_cov is the fraction of the HMM that must be aligned, summed over non-overlapping segments.'));
  makeTable(p, ['family', 'category', 'tier', 'Pfam', 'rule', 'min_cov', 'min_len',
    'median (raw)', 'notes'],
    D.families.map(f => [f.family, f.category, f.tier, f.pfam, f.rule, f.min_cov, f.min_len,
      +stats(D.counts[f.family]).med.toFixed(1), f.notes]),
    {sortBy: 0, mono: [3], name: 'defensome_map', cap: 200});
  host.append(p);
  if (D.qc_zero && D.qc_zero.length) { const p2 = el('div', {class: 'panel'});
    p2.append(el('h2', {}, 'Why each empty family is empty'), el('p', {class: 'hint'},
      'Every family with no calls at Pfam\u2019s gathering threshold, with a verdict. A zero can mean the HMM was never in the database, that a map rule filtered everything, that real copies score just below the threshold, that the input cannot contain proteins that short, or genuine absence. Those need different responses.'));
    const vc = {ACCESSION_NOT_IN_DATABASE: '#6a040f', FILTERED_BY_MAP: '#e9a13b', FOUND_BELOW_GA: '#2a9d8f',
      BELOW_GA_FAILED_VALIDATION: '#f4a261', INPUT_LACKS_SHORT_PROTEINS: '#7b2cbf',
      COMPOSITION_CANDIDATES_ONLY: '#8d99ae', NOT_DETECTED: '#6c757d', NOT_DETECTED_AT_GA_ONLY: '#6c757d'};
    const list = el('div');
    D.qc_zero.forEach(z => list.append(el('div', {style: 'display:grid;grid-template-columns:220px 1fr;gap:10px;padding:8px 0;border-bottom:1px solid var(--line)'},
      el('div', {}, el('b', {}, z.family), el('div', {class: 'mono', style: 'color:#6b7280'}, z.pfam_ids)),
      el('div', {},
        el('span', {class: 'badge', style: `background:${vc[z.verdict] || '#6c757d'};color:#fff`},
          String(z.verdict).replace(/_/g, ' ').toLowerCase().replace(/ ga\b/g, ' GA')),
        el('div', {style: 'font-size:12.5px;margin-top:4px'}, z.what_it_means || ''),
        el('div', {class: 'hint', style: 'margin:3px 0 0'},
          `GA hits ${z.raw_hits_at_GA ?? '–'} · below-GA hits ${z.raw_hits_below_GA ?? '–'} · rescued ${z.rescued_validated ?? '–'}` +
          (z['proteomes_with_<=5_proteins_under_60aa'] ? ` · proteomes lacking short proteins ${z['proteomes_with_<=5_proteins_under_60aa']}` : ''))))));
    p2.append(list);
    host.append(p2); }
  if (D.rescued && Object.keys(D.rescued).length) { const p3 = el('div', {class: 'panel'});
    p3.append(el('h2', {}, 'Found below the gathering threshold'), el('p', {class: 'hint'},
      'Genes from the rescue pass that passed the map\u2019s validation rules. They are kept out of counts.tsv and every other tab; this is the only place they are added in.'));
    const fams = Object.keys(D.rescued);
    makeTable(p3, ['sample', ...fams.flatMap(f => [f + ' at GA', f + ' rescued'])],
      SP.map((s, i) => [s, ...fams.flatMap(f => [(D.counts[f] || [])[i] ?? 0, D.rescued[f][i] ?? 0])]),
      {name: 'rescued_counts'});
    host.append(p3); }
}

/* ------------------------------------------------- FIGURE EXPORT + LEGENDS */
/* Every figure gets PNG and, where the source is SVG, vector SVG too.
   A publication figure that cannot leave the dashboard is not much use. */
function savePNG(canvas, name) {
  canvas.toBlob(b => { const u = URL.createObjectURL(b);
    const a = el('a', {href: u, download: name + '.png'}); a.click();
    setTimeout(() => URL.revokeObjectURL(u), 1000); });
}
function svgSource(svg) {
  const c = svg.cloneNode(true);
  c.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
  c.setAttribute('style', 'background:#fff');
  return '<?xml version="1.0" encoding="UTF-8"?>\n' + new XMLSerializer().serializeToString(c);
}
function saveSVG(svg, name) {
  const b = new Blob([svgSource(svg)], {type: 'image/svg+xml'});
  const u = URL.createObjectURL(b);
  const a = el('a', {href: u, download: name + '.svg'}); a.click();
  setTimeout(() => URL.revokeObjectURL(u), 1000);
}
function svgToPNG(svg, name, scale = 3) {
  const w = +svg.getAttribute('width') || 800, h = +svg.getAttribute('height') || 600;
  const img = new Image();
  img.onload = () => { const c = document.createElement('canvas');
    c.width = w * scale; c.height = h * scale;
    const g = c.getContext('2d');
    g.fillStyle = '#fff'; g.fillRect(0, 0, c.width, c.height);
    g.setTransform(scale, 0, 0, scale, 0, 0); g.drawImage(img, 0, 0);
    savePNG(c, name); };
  img.src = 'data:image/svg+xml;base64,' + btoa(unescape(encodeURIComponent(svgSource(svg))));
}
/* wrap any figure with a download toolbar */
/* `scrollAt` wraps the figure in its own scroll box once it is taller than
   that many pixels, so a 140-row heatmap does not push the colour bar, the
   legend and the rest of the page thousands of pixels down. */
function figure(node, name, extra, scrollAt) {
  const bar = el('div', {class: 'row', style: 'margin:6px 0 0'});
  if (node.tagName === 'CANVAS') {
    bar.append(el('button', {class: 'act', onclick: () => savePNG(node, name)}, 'PNG'));
  } else {
    bar.append(el('button', {class: 'act', onclick: () => svgToPNG(node, name)}, 'PNG (3×)'),
               el('button', {class: 'act', onclick: () => saveSVG(node, name)}, 'SVG (vector)'));
  }
  if (extra) bar.append(extra);
  const box = scrollAt
    ? el('div', {class: 'figscroll', style: `max-height:${scrollAt}px`}, node)
    : node;
  return el('div', {}, box, bar);
}

/* continuous colour scale legend */
function colourbar(ramp, lo, hi, label, w = 240, h = 12, ticks = 5) {
  const H = h + 34;
  const sv = svgEl('svg', {width: w, height: H, viewBox: `0 0 ${w} ${H}`});
  const N = 120;
  for (let i = 0; i < N; i++) {
    sv.append(svgEl('rect', {x: i * w / N, y: 14, width: w / N + .6, height: h,
      fill: ramp(i / (N - 1)), stroke: 'none'}));
  }
  sv.append(svgEl('rect', {x: 0, y: 14, width: w, height: h, fill: 'none', stroke: '#c9ced4'}));
  for (let k = 0; k < ticks; k++) {
    const x = w * k / (ticks - 1), v = lo + (hi - lo) * k / (ticks - 1);
    sv.append(svgEl('line', {x1: x, x2: x, y1: 14 + h, y2: 14 + h + 4, stroke: '#6b7280'}));
    const t = svgEl('text', {x: Math.min(w - 10, Math.max(10, x)), y: 14 + h + 15,
      'text-anchor': 'middle', 'font-size': 10, fill: '#6b7280'});
    t.textContent = Math.abs(v) >= 100 ? v.toFixed(0) : v.toFixed(Math.abs(v) < 1 ? 2 : 1);
    sv.append(t);
  }
  const lb = svgEl('text', {x: 0, y: 9, 'font-size': 10.5, fill: '#1c1f23'});
  lb.textContent = label; sv.append(lb);
  return sv;
}

/* categorical legend with a colour-vision-safe toggle */
const CVD_SAFE = ['#0072B2','#D55E00','#009E73','#CC79A7','#E69F00','#56B4E9','#F0E442','#000000'];
traitColour = t => { if (t == null) return '#c9ced4';
  const i = TRAITS.indexOf(t);
  return CVD ? CVD_SAFE[i % CVD_SAFE.length] : TCOL[i % TCOL.length]; };

/* ---------------------------------------------------- STATISTICS + METHODS */
function tabStats(host) {
  if (!D.stats || !Object.keys(D.stats).length) {
    host.append(el('div', {class: 'panel'}, el('h2', {}, 'Statistics'),
      el('p', {class: 'hint'}, 'No group tests found. Run report with --metadata and --group-by.')));
  } else {
    for (const [grp, rows] of Object.entries(D.stats)) {
      const p = el('div', {class: 'panel'});
      p.append(el('h2', {}, `Kruskal-Wallis by ${grp}`));
      p.append(el('div', {class: 'note', style: 'border-left-color:var(--warm)'},
        'These tests treat species as independent observations. They are not. ' +
        'Diet breadth is phylogenetically conserved in Lepidoptera, so congeners ' +
        'contribute almost no independent information while inflating the degrees ' +
        'of freedom. Read this table as a screen for what to look at, never as a ' +
        'result. A defensible test needs a phylogenetic comparative model ' +
        '(PGLS, or a phylogenetic GLMM) on a time-calibrated tree, and a count of ' +
        'independent trait transitions rather than a count of tips.'));
      makeTable(p, ['family', 'H', 'p', 'p (BH-adjusted)'],
        rows.map(r => [r.family, +Number(r.H).toFixed(3), +Number(r.p).toFixed(4),
                       +Number(r.p_bh).toFixed(4)]),
        {sortBy: 2, name: 'kruskal_' + grp});
      const sig = rows.filter(r => Number(r.p_bh) < 0.05).length;
      p.append(el('div', {class: 'hint'},
        `${sig} of ${rows.length} families pass a 5% false discovery rate. ` +
        `With ${SP.length} species and strong phylogenetic structure, expect this ` +
        `to shrink under a phylogenetic model.`));
      host.append(p);
    }
  }
  for (const [grp, gm] of Object.entries(D.group_means || {})) {
    const p = el('div', {class: 'panel'});
    p.append(el('h2', {}, `Group means by ${grp}`));
    makeTable(p, ['family', ...gm.groups],
      Object.entries(gm.data).map(([f, v]) => [f, ...v.map(x => x === null ? null : +x.toFixed(2))]),
      {sortBy: 0, name: 'group_means_' + grp});
    host.append(p);
  }
}

function tabMethods(host) {
  const p = el('div', {class: 'panel'});
  p.append(el('h2', {}, 'Provenance'), el('p', {class: 'hint'},
    'Which pipeline outputs this dashboard was built from, and which are absent. Empty rows tell you which tabs will be thin.'));
  if (D.provenance) {
    makeTable(p, ['output', 'path', 'present', 'files', 'last modified'],
      D.provenance.map(r => [r.item, r.path, r.found ? 'yes' : 'no',
        r.n_files ?? '–', r.modified ?? '–']), {sortBy: 0, name: 'provenance'});
  }
  p.append(el('div', {class: 'note'},
    `defensome ${D.tool_version || '?'} · generated ${D.generated} · source ${D.out_dir}` +
    (D.counts_basis ? ` · counts: ${D.counts_basis}` : '')));
  host.append(p);

  const m = el('div', {class: 'panel'});
  m.append(el('h2', {}, 'Methods text'), el('p', {class: 'hint'},
    'A first draft of the methods paragraph, filled in from this run. Check every number before it goes near a thesis.'));
  const nCalls = FAM.reduce((a, f) => a + D.counts[f].reduce((x, y) => x + y, 0), 0);
  const inv = Object.entries(INVARIANT).map(([f, w]) => `${f} (${w})`).join('; ');
  const txt =
`Proteomes for ${SP.length} species were searched against Pfam-A with hmmsearch
using the curated per-family gathering thresholds (--cut_ga), which are used in
preference to a fixed E-value because E-values depend on database size and are
therefore not comparable across proteomes of different sizes. Domain hits were
assigned to ${FAM.length} defensome families using explicit architecture rules:
families defined by a single domain require that domain at or above a
family-specific fraction of the HMM length, and families defined by a domain
combination require every listed domain within the same protein. HMM coverage
was computed by merging non-overlapping segments, so models that match in two
pieces are not scored on their best segment alone. This yielded ${nCalls.toLocaleString()}
gene calls. Counts were normalised to copies per 10,000 predicted proteins to
remove the dependence of raw counts on proteome size and annotation
completeness; both raw and normalised matrices are reported. Low-copy families
were used as an internal calibration check${inv ? ' (' + inv + ')' : ''}.
Analyses were performed with defensome ${D.tool_version || 'x.y.z'}.`;
  m.append(el('div', {class: 'mono',
    style: 'white-space:pre-wrap;background:#f6f8fa;padding:12px;border-radius:6px'}, txt));
  m.append(el('div', {class: 'row'}, el('button', {class: 'act',
    onclick: () => dl('methods_draft.txt', txt)}, 'download methods draft')));
  host.append(m);

  const c = el('div', {class: 'panel'});
  c.append(el('h2', {}, 'Display'), el('p', {class: 'hint'},
    'The default palette is not safe for the commonest forms of colour vision deficiency. Switch to the Okabe-Ito set for figures that go into print.'));
  c.append(el('div', {class: 'row'},
    el('button', {class: 'act' + (CVD ? ' on' : ''), onclick: e => {
      CVD = !CVD; e.target.classList.toggle('on');
      Object.keys(built).forEach(k => { if (k !== 'about') {
        const h = $('#t-' + k); if (h) h.innerHTML = ''; delete built[k]; } });
    }}, 'colour-vision-safe palette')));
  const sw = el('div', {class: 'legend'});
  TRAITS.forEach(t => sw.append(el('span', {},
    el('i', {class: 'sw', style: 'background:' + traitColour(t)}), t)));
  c.append(sw);
  host.append(c);
}

/* ---------------------------------------------------------------- ABOUT */
/* Everything about how the numbers were made, in one place and folded away:
   downloads, methods, rules, quality control. Sections build when opened. */
function tabAbout(host) {
  const secs = [
    ['downloads', 'Downloads', tabData, true],
    ['methods', 'Methods', tabMethods, false],
    ['rules', 'Gene-calling rules', tabMap, false],
    ['qc', 'Quality control', qcPanels, false],
    ...(D.completeness ? [['domains', 'Domain completeness', tabDomains, false]] : []),
    ...(!MULTI_METHOD && ((D.cmp || {}).collapse || []).length ? [['collapse', 'Isoform collapse', tabAMethods, false]] : []),
    ...(D.stats && Object.keys(D.stats).length ? [['stats', 'Group tests', tabStats, false]] : []),
  ];
  const want = (ROUTE && ROUTE.part) || '';
  secs.forEach(([id, title, fn, open]) => {
    const inner = el('div', {class: 'in'});
    const d = el('details', {class: 'sec', id: 'about-' + id}, el('summary', {}, title), inner);
    let done = false;
    const build = () => { if (done) return; done = true;
      try { fn(inner); } catch (e) { inner.append(el('pre', {class: 'mono'}, String(e && e.stack || e))); } };
    d.addEventListener('toggle', () => { if (d.open) build(); });
    if (open || want === id) { d.open = true; build(); }
    host.append(d);
  });
  tabAbout._secs = secs.map(x => x[0]);
}

/* --------------------------------------------------------- SPECIES COMPARE */
/* Mann-Whitney U with a tie-corrected normal approximation. Two-sided.
   Non-parametric because gene-family counts are bounded below, skewed and
   occasionally zero-inflated, none of which a t-test likes. */
function mannWhitney(a, b) {
  const n1 = a.length, n2 = b.length;
  if (n1 < 3 || n2 < 3) return {U: NaN, p: NaN, z: NaN};
  const all = a.map(v => [v, 0]).concat(b.map(v => [v, 1])).sort((x, y) => x[0] - y[0]);
  const rank = new Array(all.length);
  let tie = 0, i = 0;
  while (i < all.length) {
    let j = i; while (j + 1 < all.length && all[j + 1][0] === all[i][0]) j++;
    const r = (i + j) / 2 + 1, t = j - i + 1;
    for (let k = i; k <= j; k++) rank[k] = r;
    tie += t * t * t - t; i = j + 1;
  }
  let R1 = 0; all.forEach((v, k) => { if (v[1] === 0) R1 += rank[k]; });
  const U1 = R1 - n1 * (n1 + 1) / 2, U = Math.min(U1, n1 * n2 - U1);
  const N = n1 + n2, mu = n1 * n2 / 2;
  const sd = Math.sqrt(n1 * n2 / 12 * ((N + 1) - tie / (N * (N - 1))));
  if (!(sd > 0)) return {U, p: 1, z: 0};
  const z = (Math.abs(U1 - mu) - 0.5) / sd;
  // two-sided normal tail, Abramowitz and Stegun 7.1.26
  const t_ = 1 / (1 + 0.2316419 * z);
  const phi = Math.exp(-z * z / 2) / Math.sqrt(2 * Math.PI);
  const tail = phi * t_ * (0.319381530 + t_ * (-0.356563782 + t_ * (1.781477937 +
    t_ * (-1.821255978 + t_ * 1.330274429))));
  return {U, z, p: Math.min(1, 2 * tail)};
}
/* Benjamini-Hochberg, returning adjusted p in the input order. */
function bh(ps) {
  const idx = ps.map((p, i) => [p, i]).filter(x => isFinite(x[0])).sort((x, y) => x[0] - y[0]);
  const m = idx.length, out = ps.map(() => NaN);
  let prev = 1;
  for (let k = m - 1; k >= 0; k--) {
    const q = Math.min(prev, idx[k][0] * m / (k + 1));
    out[idx[k][1]] = q; prev = q;
  }
  return out;
}

function valuesOfTax(rank, value) {
  const out = [];
  D.species.forEach((sp, i) => {
    const v = (rank === D.trait_name && !(sp.tax && rank in sp.tax))
      ? sp.trait : (sp.tax || {})[rank];
    if (String(v ?? '') === value) out.push(i);
  });
  return out;
}
function taxaValues(rank) {
  const s = new Set();
  D.species.forEach(sp => {
    const v = (rank === D.trait_name && !(sp.tax && rank in sp.tax))
      ? sp.trait : (sp.tax || {})[rank];
    if (v != null && v !== '') s.add(String(v));
  });
  return [...s].sort();
}

const REST = '__rest__';
/* One row per family, one column per group: the median, coloured by how it
   sits against the family's other groups. The whole trait at a glance. */
function traitGrid(rank, kind, onPick) {
  const vals = taxaValues(rank).map(v => [v, valuesOfTax(rank, v)]).filter(([, ix]) => ix.length >= 3)
    .sort((a, b) => b[1].length - a[1].length).slice(0, 14);
  if (vals.length < 2) return el('div', {class: 'empty'}, 'Fewer than two groups with three or more species.');
  const fams = [...VARFAM].sort((a, b) => CATS.indexOf(FMETA[a].category) - CATS.indexOf(FMETA[b].category));
  const med = fams.map(f => { const v = col(kind, f); return vals.map(([, ix]) => stats(ix.map(i => v[i])).med); });
  const mg = measurer(11.5);
  const L = Math.ceil(Math.max(...fams.map(f => mg.measureText(f).width)) + 26);
  const heads = vals.map(([v, ix]) => `${nice(v)} (${ix.length})`);
  const T = Math.ceil(Math.min(180, Math.max(...heads.map(h => mg.measureText(h).width)) * .72 + 24));
  const CW = Math.max(44, Math.min(86, Math.floor(820 / vals.length))), CH = 19;
  const W = L + vals.length * CW + Math.ceil(Math.max(...heads.map(h => mg.measureText(h).width)) * .77) + 6;
  const H = T + fams.length * CH + 4;
  const sv = svgEl('svg', {width: W, height: H, viewBox: `0 0 ${W} ${H}`});
  vals.forEach(([v], j) => { const x = L + j * CW + CW / 2;
    const t = svgEl('text', {x, y: T - 8, 'font-size': 11.5, transform: `rotate(-40 ${x} ${T - 8})`,
      style: 'cursor:pointer', fill: '#1b4965'});
    t.textContent = heads[j]; t.addEventListener('click', () => onPick && onPick(v)); sv.append(t); });
  fams.forEach((f, r) => { const y = T + r * CH, row = med[r];
    const m = row.reduce((a, b) => a + b, 0) / row.length;
    const sd = Math.sqrt(row.reduce((a, b) => a + (b - m) ** 2, 0) / Math.max(1, row.length - 1)) || 1;
    sv.append(svgEl('rect', {x: 4, y: y + 4, width: 6, height: CH - 8, rx: 1, fill: catColour(FMETA[f].category)}));
    const lab = svgEl('text', {x: L - 8, y: y + CH / 2 + 4, 'text-anchor': 'end', 'font-size': 11.5,
      fill: '#1b4965', style: 'cursor:pointer'});
    lab.textContent = f; lab.addEventListener('click', () => { location.hash = famHref(f); }); sv.append(lab);
    row.forEach((x, j) => { const z = (x - m) / sd;
      const c = svgEl('rect', {x: L + j * CW, y, width: CW - 2, height: CH - 2, rx: 2,
        fill: rdbu((Math.max(-2, Math.min(2, z)) + 2) / 4)});
      c.addEventListener('mousemove', e => showTip(e, `<b>${f}</b> in ${nice(vals[j][0])}<br>median ${fmt(x, 2)}`));
      c.addEventListener('mouseleave', hideTip); sv.append(c);
      const t = svgEl('text', {x: L + j * CW + CW / 2 - 1, y: y + CH / 2 + 3.5, 'text-anchor': 'middle',
        'font-size': 10.5, fill: Math.abs(z) > 1.3 ? '#fff' : '#1c1f23', 'pointer-events': 'none'});
      t.textContent = x >= 100 ? Math.round(x) : x >= 10 ? x.toFixed(0) : x.toFixed(1); sv.append(t); });
  });
  traitGrid._vals = vals.map(v => v[0]);
  return sv;
}

function taxaComparePanel(groupCols) {
  const p = el('div', {class: 'panel'});
  const mk = (label, node) => el('div', {class: 'ctl'}, el('label', {}, label), node);
  const rankSel = el('select', {}, ...groupCols.map(c => el('option', {value: c}, c)));
  const aSel = el('select', {}), bSel = el('select', {});
  const kindSel = el('select', {},
    ...(D.per10k ? [el('option', {value: 'per10k'}, `per 10k ${D.size_unit || 'proteins'}`)] : []),
    el('option', {value: 'raw'}, 'copies'));
  p.append(el('div', {class: 'row', style: 'align-items:flex-end'}, mk('compare by', rankSel), mk('values', kindSel)));
  const gridBox = el('div', {style: 'overflow:auto'});
  const pair = el('div', {style: 'min-width:0'});
  pair.append(el('div', {class: 'row', style: 'align-items:flex-end'}, mk('group', aSel),
    el('span', {style: 'padding-bottom:8px;color:var(--mut)'}, 'vs'), mk('against', bSel)));
  const body = el('div'); pair.append(body);
  p.append(el('div', {class: 'split'}, gridBox, pair));

  function fillValues() {
    const rank = pick(rankSel, groupCols), vals = taxaValues(rank);
    const lab = v => `${nice(v)} (${valuesOfTax(rank, v).length})`;
    const was = [aSel.value, bSel.value];
    aSel.innerHTML = ''; bSel.innerHTML = '';
    vals.forEach(v => { aSel.append(el('option', {value: v}, lab(v))); bSel.append(el('option', {value: v}, lab(v))); });
    bSel.append(el('option', {value: REST}, 'all others'));
    const big = [...vals].sort((x, y) => valuesOfTax(rank, y).length - valuesOfTax(rank, x).length);
    aSel.value = vals.includes(was[0]) ? was[0] : (big[0] ?? '');
    bSel.value = (vals.includes(was[1]) || was[1] === REST) ? was[1] : REST;
  }
  function drawGrid() {
    gridBox.innerHTML = '';
    gridBox.append(figure(traitGrid(pick(rankSel, groupCols), pick(kindSel, ['per10k', 'raw']),
      v => { aSel.value = v; bSel.value = REST; draw(); }), 'by_' + pick(rankSel, groupCols)));
    gridBox.append(el('div', {class: 'legend'},
      el('span', {}, el('i', {class: 'sw', style: 'background:' + rdbu(.08)}), 'fewer'),
      el('span', {}, el('i', {class: 'sw', style: 'background:' + rdbu(.92)}), 'more'),
      el('span', {style: 'color:var(--mut)'}, 'than the same family in the other groups; numbers are medians')));
  }
  function draw() {
    body.innerHTML = '';
    const rank = pick(rankSel, groupCols), kind = pick(kindSel, ['per10k', 'raw']);
    const A = aSel.value, B = bSel.value;
    const ia = valuesOfTax(rank, A);
    const inA = new Set(ia);
    const ib = B === REST ? taxaValues(rank).filter(v => v !== A).flatMap(v => valuesOfTax(rank, v))
                            .filter(i => !inA.has(i)) : valuesOfTax(rank, B);
    const bName = B === REST ? 'all others' : nice(B);
    if (!A || A === B || ia.length < 2 || ib.length < 2) {
      body.append(el('div', {class: 'empty'}, 'Pick two different groups with at least two species each.')); return; }
    const rows = VARFAM.map(f => {
      const v = col(kind, f);
      const va = ia.map(i => v[i]), vb = ib.map(i => v[i]);
      const sa = stats(va), sb = stats(vb);
      return {f, medA: sa.med, medB: sb.med, diff: sa.med - sb.med,
              lfc: Math.log2((sa.med + 0.5) / (sb.med + 0.5)), p: mannWhitney(va, vb).p};
    });
    const q = bh(rows.map(r => r.p));
    rows.forEach((r, i) => r.q = q[i]);
    const SIG = 0.05;
    const fold = r => { const x = (r.medA + 0.5) / (r.medB + 0.5); return x >= 1 ? `${x.toFixed(1)}×` : `1/${(1 / x).toFixed(1)}`; };
    const TOP = 15;
    const shown = rows.length <= TOP || taxaComparePanel._all ? rows
      : [...rows].sort((x, y) => (x.q < SIG) === (y.q < SIG) ? Math.abs(y.lfc) - Math.abs(x.lfc) : (x.q < SIG ? -1 : 1)).slice(0, TOP);
    const items = [...shown].sort((x, y) => y.lfc - x.lfc).map(r => ({
      label: r.f + (r.q < SIG ? ' *' : ''), value: r.lfc, vlabel: fold(r), href: famHref(r.f),
      colour: r.q < SIG ? (r.lfc >= 0 ? '#1b4965' : '#bc4749') : '#d3d8dd',
      extra: `${nice(A)} ${fmt(r.medA, 2)} vs ${bName} ${fmt(r.medB, 2)} · q ${fmt(r.q, 3)}`}));
    body.append(el('div', {class: 'legend', style: 'margin:4px 0 8px'},
      el('span', {}, el('i', {class: 'sw', style: 'background:#1b4965'}), `more in ${nice(A)}`),
      el('span', {}, el('i', {class: 'sw', style: 'background:#bc4749'}), `more in ${bName}`),
      el('span', {}, el('i', {class: 'sw', style: 'background:#d3d8dd'}), 'no clear difference')));
    body.append(figure(barChart(items, 620, items.length * 18, true), `compare_${rank}_${A}_vs_${B}`));
    if (rows.length > TOP) body.append(el('button', {class: 'act', style: 'margin-top:6px',
      onclick: () => { taxaComparePanel._all = !taxaComparePanel._all; draw(); }},
      taxaComparePanel._all ? `show the ${TOP} largest differences` : `show all ${rows.length} families`));
    const tb = el('div');
    makeTable(tb, ['family', 'category', `median ${nice(A)}`, `median ${bName}`,
      'difference', 'p', 'q (BH)'],
      rows.map(r => [r.f, FMETA[r.f].category, +fmt(r.medA, 2), +fmt(r.medB, 2),
        +r.diff.toFixed(2), r.p, r.q]),
      {sortBy: 6, sci: [5, 6], name: `compare_${rank}_${A}_${B}`});
    body.append(collapsible('table', tb));
  }
  rankSel.onchange = () => { fillValues(); drawGrid(); draw(); };
  kindSel.onchange = () => { drawGrid(); draw(); };
  [aSel, bSel].forEach(x => x.onchange = draw);
  fillValues(); drawGrid(); draw();
  taxaComparePanel._draw = draw; taxaComparePanel._sel = {rankSel, aSel, bSel};
  return p;
}

function tabCompare(host) {
  const cols = [...new Set([...TRAIT_COLS, ...(D.trait_name ? [D.trait_name] : []), ...RANK_COLS, ...TAX_COLS])];
  if (cols.length) host.append(taxaComparePanel(cols));
  /* two genomes, for the specialist question */
  const p = el('div', {class: 'panel'});
  const opts = () => SP.map(s => el('option', {value: s}, fullName(s)));
  const a1 = el('select', {style: 'max-width:340px'}, ...opts()), b1 = el('select', {style: 'max-width:340px'}, ...opts());
  if (SP.length > 1) b1.value = SP[1];
  const body = el('div');
  p.append(el('div', {class: 'row', style: 'align-items:flex-end'},
    el('div', {class: 'ctl'}, el('label', {}, 'species'), a1),
    el('span', {style: 'padding-bottom:8px;color:var(--mut)'}, 'vs'),
    el('div', {class: 'ctl'}, el('label', {}, 'species'), b1)), body);
  host.append(el('details', {class: 'sec'}, el('summary', {}, 'Two species'), el('div', {class: 'in'}, p)));
  const draw = () => { body.innerHTML = '';
    const A = pick(a1, SP), B = pick(b1, SP), kind = D.per10k ? 'per10k' : 'raw';
    const ia = SPI[A], ib = SPI[B];
    const items = VARFAM.map(f => { const v = col(kind, f);
      return {label: f, value: v[ia] - v[ib], href: famHref(f),
        colour: v[ia] >= v[ib] ? '#1b4965' : '#bc4749',
        extra: `${sciName(A)} ${fmt(v[ia])} vs ${sciName(B)} ${fmt(v[ib])}`}; })
      .sort((x, y) => y.value - x.value);
    body.append(figure(barChart(items, 880, Math.max(240, items.length * 16), true), `compare_${A}_vs_${B}`));
  };
  [a1, b1].forEach(s => s.onchange = draw); draw();
}


/* ---------------------------------------------------- ANNOTATION METHODS */
/* The method comparison is the central result of a multi-method run, so it
   gets its own tab: one panel per question, each with the numbers beside it. */
const METHOD_ORDER = ['DToL', 'genome_guided', 'denovo'];
const methodColour = m => (CVD ? {DToL: '#0072B2', genome_guided: '#009E73', denovo: '#D55E00'}
                               : RING_FIXED.method)[m] || '#8d99ae';
const methodsIn = recs => { const s = new Set(recs.map(r => r.method));
  return [...METHOD_ORDER.filter(m => s.has(m)), ...[...s].filter(m => !METHOD_ORDER.includes(m)).sort()]; };

/* Freedman-Diaconis bins, which is the right default when the point is to see
   a long tail rather than to match someone's preferred bin count. */
function histogram(values, W, H, xlab, ylab) {
  const v = values.filter(x => isFinite(x)).sort((a, b) => a - b);
  if (v.length < 2) return el('div', {class: 'note'}, 'Not enough values to plot.');
  const st = stats(v), iqr = st.q3 - st.q1;
  const bw0 = iqr > 0 ? 2 * iqr / Math.cbrt(v.length) : (st.max - st.min) / 20;
  const nb = Math.max(6, Math.min(60, Math.ceil((st.max - st.min) / (bw0 || 1))));
  const lo = st.min, hi = st.max + 1e-9, bw = (hi - lo) / nb;
  const counts = new Array(nb).fill(0);
  v.forEach(x => counts[Math.min(nb - 1, Math.floor((x - lo) / bw))]++);
  const top = Math.max(...counts) * 1.08;
  const ML = 56, MB = 46, MT = 10, MR = 14;
  const sv = svgEl('svg', {width: W, height: H, viewBox: `0 0 ${W} ${H}`});
  const X = x => ML + (W - ML - MR) * (x - lo) / (hi - lo);
  const Y = n => MT + (H - MT - MB) * (1 - n / top);
  for (let k = 0; k <= 4; k++) { const n = top * k / 4;
    sv.append(svgEl('line', {x1: ML, x2: W - MR, y1: Y(n), y2: Y(n), stroke: '#eceff2'}));
    const t = svgEl('text', {x: ML - 8, y: Y(n) + 4, 'text-anchor': 'end',
      'font-size': 10, fill: '#6b7280'}); t.textContent = Math.round(n); sv.append(t); }
  counts.forEach((n, i) => { if (!n) return;
    const x0 = X(lo + i * bw), x1 = X(lo + (i + 1) * bw);
    const r = svgEl('rect', {x: x0, y: Y(n), width: Math.max(1, x1 - x0 - 1),
      height: Math.max(1, Y(0) - Y(n)), fill: '#1b4965', 'fill-opacity': .85});
    r.addEventListener('mousemove', e => showTip(e,
      `${fmt(lo + i * bw, 1)} to ${fmt(lo + (i + 1) * bw, 1)}<br><b>${n}</b> ${ylab || ''}`));
    r.addEventListener('mouseleave', hideTip); sv.append(r); });
  // median line, because the eye reads the mode and the median is the number
  const med = svgEl('line', {x1: X(st.med), x2: X(st.med), y1: MT, y2: Y(0),
    stroke: '#bc4749', 'stroke-width': 1.5, 'stroke-dasharray': '4 3'});
  sv.append(med);
  const ml = svgEl('text', {x: X(st.med) + 5, y: MT + 11, 'font-size': 10.5, fill: '#bc4749'});
  ml.textContent = 'median ' + fmt(st.med, 1); sv.append(ml);
  for (let k = 0; k <= 4; k++) { const x = lo + (hi - lo) * k / 4;
    const t = svgEl('text', {x: X(x), y: H - MB + 16, 'text-anchor': 'middle',
      'font-size': 10, fill: '#6b7280'}); t.textContent = fmt(x, 1); sv.append(t); }
  const xl = svgEl('text', {x: (ML + W - MR) / 2, y: H - 10, 'text-anchor': 'middle',
    'font-size': 11, fill: '#6b7280'}); xl.textContent = xlab; sv.append(xl);
  const yl = svgEl('text', {x: 14, y: (MT + H - MB) / 2, 'font-size': 11, fill: '#6b7280',
    transform: `rotate(-90 14 ${(MT + H - MB) / 2})`, 'text-anchor': 'middle'});
  yl.textContent = ylab; sv.append(yl);
  return sv;
}

function groupedBars(groups, series, valueOf, W, H, ylab) {
  /* groups: x categories; series: bar within group; valueOf(g,s) -> number|null */
  const ML = 64, MB = 86, MT = 14, MR = 12;
  const vals = []; groups.forEach(g => series.forEach(s => { const v = valueOf(g, s); if (v != null) vals.push(v); }));
  const top = (Math.max(...vals, 1)) * 1.08;
  const sv = svgEl('svg', {width: W, height: H, viewBox: `0 0 ${W} ${H}`});
  const Y = v => MT + (H - MT - MB) * (1 - v / top), gw = (W - ML - MR) / groups.length;
  for (let k = 0; k <= 4; k++) { const v = top * k / 4;
    sv.append(svgEl('line', {x1: ML, x2: W - MR, y1: Y(v), y2: Y(v), stroke: '#eceff2'}));
    const t = svgEl('text', {x: ML - 8, y: Y(v) + 4, 'text-anchor': 'end', 'font-size': 10, fill: '#6b7280'});
    t.textContent = v >= 1000 ? (v / 1000).toFixed(v >= 10000 ? 0 : 1) + 'k' : v.toFixed(0); sv.append(t); }
  const yl = svgEl('text', {x: 14, y: (MT + H - MB) / 2, 'font-size': 11, fill: '#6b7280',
    transform: `rotate(-90 14 ${(MT + H - MB) / 2})`, 'text-anchor': 'middle'}); yl.textContent = ylab; sv.append(yl);
  const gstep = labelStep(groups.length, Math.max(8, Math.floor((W - ML - MR) / 15)));
  groups.forEach((g, i) => {
    const bw = gw * .78 / series.length, x0 = ML + gw * i + gw * .11;
    series.forEach((s, j) => { const v = valueOf(g, s); if (v == null) return;
      const r = svgEl('rect', {x: x0 + j * bw, y: Y(v), width: bw - 1, height: Math.max(.5, Y(0) - Y(v)), fill: methodColour(s)});
      r.addEventListener('mousemove', e => showTip(e, `<b>${nice(g)}</b><br>${s}: ${v.toLocaleString()}`));
      r.addEventListener('mouseleave', hideTip); sv.append(r); });
    if (i % gstep) return;
    const t = svgEl('text', {x: x0 + gw * .39, y: H - MB + 14, 'font-size': 10.5, 'font-style': 'italic',
      'text-anchor': 'end', transform: `rotate(-40 ${x0 + gw * .39} ${H - MB + 14})`});
    t.textContent = nice(g); sv.append(t); });
  return sv;
}

function pairedDots(rows, species, methods, W, H, ylab) {
  /* rows: {species, method, value}; one grey line per species joins its methods */
  const ML = 64, MB = 86, MT = 14, MR = 12;
  const v = rows.map(r => r.value), top = Math.max(...v, 1) * 1.1, bot = Math.min(...v, 0);
  const sv = svgEl('svg', {width: W, height: H, viewBox: `0 0 ${W} ${H}`});
  const Y = x => MT + (H - MT - MB) * (1 - (x - bot) / (top - bot)), gw = (W - ML - MR) / species.length;
  for (let k = 0; k <= 4; k++) { const x = bot + (top - bot) * k / 4;
    sv.append(svgEl('line', {x1: ML, x2: W - MR, y1: Y(x), y2: Y(x), stroke: '#eceff2'}));
    const t = svgEl('text', {x: ML - 8, y: Y(x) + 4, 'text-anchor': 'end', 'font-size': 10, fill: '#6b7280'});
    t.textContent = x.toFixed(0); sv.append(t); }
  const yl = svgEl('text', {x: 14, y: (MT + H - MB) / 2, 'font-size': 11, fill: '#6b7280',
    transform: `rotate(-90 14 ${(MT + H - MB) / 2})`, 'text-anchor': 'middle'}); yl.textContent = ylab; sv.append(yl);
  const pstep = labelStep(species.length, Math.max(8, Math.floor((W - ML - MR) / 15)));
  species.forEach((sp, i) => {
    const pts = methods.map((m, j) => { const r = rows.find(q => q.species === sp && q.method === m);
      return r ? {x: ML + gw * i + gw * (.25 + .5 * j / Math.max(methods.length - 1, 1)), y: Y(r.value), m, v: r.value} : null; })
      .filter(Boolean);
    if (pts.length > 1) sv.append(svgEl('path', {d: pts.map((p, k) => `${k ? 'L' : 'M'}${p.x},${p.y}`).join(''),
      stroke: '#c9ced4', 'stroke-width': 1.5, fill: 'none'}));
    pts.forEach(p => { const c = svgEl('circle', {cx: p.x, cy: p.y, r: 5.5, fill: methodColour(p.m), stroke: '#fff', 'stroke-width': 1.2});
      c.addEventListener('mousemove', e => showTip(e, `<b>${fullName(sp)}</b><br>${p.m}: ${p.v.toFixed(0)}`));
      c.addEventListener('mouseleave', hideTip); sv.append(c); });
    if (i % pstep) return;
    const t = svgEl('text', {x: ML + gw * (i + .5), y: H - MB + 14, 'font-size': 10.5, 'font-style': 'italic',
      'text-anchor': 'end', transform: `rotate(-40 ${ML + gw * (i + .5)} ${H - MB + 14})`});
    t.textContent = nice(sp); sv.append(t); });
  return sv;
}

function recoveryHeatmap(recs, methods) {
  /* median recovery per family x method, on a log2 scale centred on parity */
  const fams = [...new Set(recs.map(r => r.family))];
  const med = {}; fams.forEach(f => methods.forEach(m => {
    const v = recs.filter(r => r.family === f && r.method === m && r.recovery != null && isFinite(r.recovery)).map(r => r.recovery);
    if (v.length) { const s = [...v].sort((a, b) => a - b); med[f + '|' + m] = s[Math.floor(s.length / 2)]; } }));
  fams.sort((a, b) => (med[a + '|' + methods[methods.length - 1]] ?? 9) - (med[b + '|' + methods[methods.length - 1]] ?? 9));
  const CW = 110, CH = 20, L = 130, T = 30;
  const W = L + methods.length * CW + 20, H = T + fams.length * CH + 10;
  const sv = svgEl('svg', {width: W, height: H, viewBox: `0 0 ${W} ${H}`});
  const col = r => rdbu(1 - (Math.max(-2, Math.min(2, Math.log2(r))) + 2) / 4);
  methods.forEach((m, j) => { const t = svgEl('text', {x: L + j * CW + CW / 2, y: T - 10, 'text-anchor': 'middle', 'font-size': 11.5, 'font-weight': 600});
    t.textContent = m; sv.append(t); });
  fams.forEach((f, i) => {
    const t = svgEl('text', {x: L - 8, y: T + i * CH + CH * .7, 'text-anchor': 'end', 'font-size': 11}); t.textContent = f; sv.append(t);
    methods.forEach((m, j) => { const r = med[f + '|' + m];
      const x = L + j * CW, y = T + i * CH;
      const rect = svgEl('rect', {x: x + 1, y: y + 1, width: CW - 2, height: CH - 2, fill: r == null ? '#f1f3f5' : col(r)});
      rect.addEventListener('mousemove', e => showTip(e, `<b>${f}</b> · ${m}<br>median recovery ${r == null ? '–' : r.toFixed(2)}`));
      rect.addEventListener('mouseleave', hideTip); sv.append(rect);
      if (r != null) { const tx = svgEl('text', {x: x + CW / 2, y: y + CH * .7, 'text-anchor': 'middle', 'font-size': 10.5,
        fill: Math.abs(Math.log2(r)) > 1.1 ? '#fff' : '#212529'}); tx.textContent = r.toFixed(2); sv.append(tx); } }); });
  return sv;
}

function tabAMethods(host) {
  const C = D.cmp || {}, long = C.long || [], rec = C.recovery || [], col = C.collapse || [];
  if (!long.length && !col.length) {
    host.append(el('div', {class: 'panel'}, el('h2', {}, 'Annotation methods'),
      el('p', {class: 'hint'}, 'No method comparison in this run. It appears when the sample sheet assigns several methods to the same species and `compare` has run.')));
    return; }
  const methods = methodsIn(long.length ? long : col).filter(m => m);
  const multi = methods.length > 1;         // a method comparison, or just one method
  const BAR_SPECIES_MAX = 60;
  const species = [...new Set((long.length ? long : col).map(r => r.species))].sort();
  const ref = (rec.length && !methods.includes(rec[0].method)) ? '' : 'DToL';
  const lab = long.length && long[0] ? (D.cmp_label || 'CORE defensome genes') : '';
  host.append(el('div', {class: 'grid g3'},
    kpi(species.length.toLocaleString(), 'species'),
    kpi(methods.length || null, methods.length === 1 ? 'annotation method' : 'annotation methods'),
    kpi((long.length ? new Set(long.map(r => r.sample_id)).size : col.length).toLocaleString(),
        'proteomes')));
  if (multi) host.append(el('div', {class: 'legend', style: 'margin:8px 0 14px'},
    methods.map(m => el('span', {}, el('i', {class: 'sw', style: 'background:' + methodColour(m)}), m))));

  if (col.length) { const p = el('div', {class: 'panel'});
    p.append(el('h2', {}, (multi ? '1 · ' : '') + 'Genes after isoform collapse'),
      el('p', {class: 'hint'}, multi
        ? 'Proteins become genes once isoforms are merged. Transcriptome "genes" still run above the genome annotation of the same species because assemblies split one real gene across several components; that inflates any per-gene denominator, which is why the comparison below uses counts.'
        : 'Proteins become genes once isoforms are merged. With one annotation method there is nothing to compare across methods, so what matters is the spread: a proteome that collapses far more than its neighbours was listing isoforms, and its gene total is not comparable with theirs untouched.'));
    /* A bar per species is a comparison when there are three species and two
       methods. At catalogue scale it is an 11,200-pixel row of hairlines,
       which is what this used to draw. Chart what the run can actually show. */
    if (multi && species.length <= BAR_SPECIES_MAX) {
      p.append(figure(groupedBars(species, methods,
        (g, s) => { const r = col.find(q => q.species === g && q.method === s);
          return r ? r.n_genes_out : null; },
        Math.max(960, species.length * 80), 380, 'genes'), 'collapse_by_method'));
    } else {
      const pct = col.map(r => +r.collapse_pct).filter(v => isFinite(v));
      p.append(figure(histogram(pct, 880, 300, 'collapse %', 'proteomes'),
        'collapse_distribution'));
      const st = stats(pct), iqr = st.q3 - st.q1;
      p.append(el('div', {class: 'note'}, multi
        ? `${species.length.toLocaleString()} species is too many for a bar each; the distribution is above and every row is in the table below.`
        : `Median ${st.med.toFixed(1)}% collapsed, range ${st.min.toFixed(1)} to ${st.max.toFixed(1)}%. `
          + `${pct.filter(v => v > st.q3 + 1.5 * iqr).length} proteomes collapse far more than the rest; `
          + 'sort the table by collapse % to see which.'));
    }
    makeTable(p, ['species', 'method', 'proteins in', 'genes out', 'collapse %', 'utrorf dropped'],
      col.map(r => [nice(r.species || r.sample_id), r.method, r.n_proteins_in,
        r.n_genes_out, r.collapse_pct, r.dropped_utrorf ?? 0]),
      {sortBy: 4, desc: true, ital: [0], name: 'collapse_stats'});
    host.append(p); }

  if (long.length) { const p = el('div', {class: 'panel'});
    const tot = {}; long.forEach(r => { const k = r.species + '|' + r.method; tot[k] = (tot[k] || 0) + (r.value || 0); });
    const rows = Object.entries(tot).map(([k, v]) => ({species: k.split('|')[0], method: k.split('|')[1], value: v}));
    p.append(el('h2', {}, '2 · The same genome, measured by each method'), el('p', {class: 'hint'},
      'Total CORE defensome per species under each method. A grey line joins one species, so the slope is the method effect with the genome held fixed. Flat lines mean the methods agree.'));
    p.append(figure(pairedDots(rows, species, methods, Math.max(960, species.length * 80), 400, 'CORE defensome genes'), 'paired_totals'));
    host.append(p); }

  if (rec.length) { const others = methods.filter(m => m !== ref && rec.some(r => r.method === m));
    const p = el('div', {class: 'panel'});
    p.append(el('h2', {}, `3 · Which families each method recovers, relative to ${ref || 'the reference'}`),
      el('p', {class: 'hint'}, 'Median over species of (method count ÷ reference count). 1.00 is parity; red is loss, blue is excess. Families are ordered by the last column, so the worst-recovered sit at the top.'));
    const hm = recoveryHeatmap(rec, others);
    p.append(el('div', {style: 'display:flex;gap:18px;align-items:flex-start;flex-wrap:wrap'},
      figure(hm, 'recovery_heatmap'),
      el('div', {}, colourbar(t => rdbu(1 - t), -2, 2, 'log2 recovery (0 = parity)', 220, 12),
        el('div', {class: 'hint', style: 'max-width:260px;margin-top:8px'},
          'Recovery above 1 is rarely extra biology. It is usually residual isoform inflation in the transcriptome, or members the genome annotation missed; check those families in the Trees tab.'))));
    const pipes = [...new Set(rec.map(r => r.ref_pipeline).filter(p => p && p !== 'unknown'))];
    if (pipes.length > 1) {
      const cells = pipes.map(pp => [pp, ...others.map(m => { const v = rec.filter(r => r.ref_pipeline === pp && r.method === m && r.recovery != null).map(r => r.recovery).sort((a, b) => a - b);
        return v.length ? +v[Math.floor(v.length / 2)].toFixed(3) : null; })]);
      p.append(el('h2', {style: 'margin-top:16px'}, 'Split by the pipeline that produced the reference'),
        el('p', {class: 'hint'}, 'If these rows differ, "recovery against DToL" averages over two different denominators and should be reported separately.'));
      makeTable(p, ['reference pipeline', ...others], cells, {name: 'recovery_by_pipeline'}); }
    host.append(p); }

  if ((C.tests || []).length) { const p = el('div', {class: 'panel'});
    p.append(el('h2', {}, '4 · Which differences survive testing'), el('p', {class: 'hint'},
      'Wilcoxon signed-rank tests paired within species, Benjamini-Hochberg adjusted. Pairing within species removes phylogeny from this comparison: the same genome is measured twice.'));
    const nPairs = Math.min(...C.tests.map(r => +r.n_pairs));
    const floor = 2 / Math.pow(2, nPairs);
    p.append(el('div', {class: 'note', style: floor >= 0.05 ? 'border-left-color:var(--warm)' : ''},
      `With ${nPairs} pairs, the smallest p-value a two-sided signed-rank test can return is 2/2^${nPairs} = ${floor.toPrecision(2)}.` +
      (floor >= 0.05 ? ' That is above 0.05, so no difference can reach significance here whatever the data; add species before reading this table.' : '')));
    makeTable(p, ['family', 'method A', 'method B', 'pairs', 'median A', 'median B', 'p', 'p (BH)'],
      C.tests.map(r => [r.family, r.method_a, r.method_b, r.n_pairs, +(+r.median_a).toFixed(2), +(+r.median_b).toFixed(2),
        +(+r.p).toExponential(2), +(+r.p_bh).toExponential(2)]), {sortBy: 7, name: 'paired_tests'});
    host.append(p); }

  /* small proteins: where metallothioneins and MATE went */
  const so = D.short_orfs || [], comp = D.mt_composition || {}, resc = D.rescued || {};
  if (so.length || Object.keys(comp).length || resc.Metallothionein || resc.MATE) {
    const p = el('div', {class: 'panel'});
    p.append(el('h2', {}, '5 · Small proteins, and why they go missing'), el('p', {class: 'hint'},
      'Insect metallothioneins are 40–64 residues and no Pfam model covers them. TransDecoder keeps ORFs of at least 100 residues by default, so a transcriptome proteome cannot contain them. Each column is an independent line of evidence; none of them is merged into counts.tsv.'));
    const idx = s => SPI[s];
    const rows = SP.map((s, i) => { const r = so.find(q => q.sample_id === s);
      return [s, D.species[i].method || '', (D.counts.Metallothionein || [])[i] ?? 0,
        resc.Metallothionein ? (resc.Metallothionein[i] ?? 0) : '–', comp[s] ?? 0,
        r ? r.mt_like_genes : '–', (D.counts.MATE || [])[i] ?? 0, resc.MATE ? (resc.MATE[i] ?? 0) : '–']; });
    makeTable(p, ['sample', 'method', 'MT at Pfam GA', 'MT rescued below GA', 'MT by composition',
      'MT from transcripts (short-orfs)', 'MATE at GA', 'MATE rescued'], rows, {sortBy: 1, name: 'small_proteins'});
    if (so.length) {
      const bm = {}; so.forEach(r => { (bm[r.method] = bm[r.method] || []).push(r.mt_like_genes); });
      p.append(el('div', {class: 'note'}, 'Metallothionein-like genes recovered from transcripts, median by method: ' +
        Object.entries(bm).map(([m, v]) => { const s = [...v].sort((a, b) => a - b); return `${m} ${s[Math.floor(s.length / 2)]}`; }).join(' · ') +
        '. Composition-screen candidates: confirm a few by alignment before citing counts.')); }
    host.append(p); }
}

/* ------------------------------------------------------------------ BOOT */
/* Tabs. The fourth entry says whether the run has anything to put in the tab;
   a tab with nothing in it is left out of the nav rather than shown as a page
   that explains its own emptiness. The Methods tab already lists which
   pipeline outputs are missing and what to run for them, so nothing is lost.
   Related views are one tab with two sections, not two tabs: a reader does not
   think of "ordination" and "correlations" as separate errands. */
/* =====================================================================
   SITE: data loading, routing, search, and the entity pages.
   The landing payload (D) holds the comparative matrices. Everything per gene
   lives in data/: one script per species (genes, domain coordinates,
   sequences), one per family (members), one per large tree. They are loaded
   with <script> tags because that works from file://, where fetch() does not,
   so the folder opens straight off disk and can be hosted unchanged.
   ===================================================================== */
const SITE = D.site || {genes: false, sequences: false, family_fasta: {}, family_bytes: {}};
const safeKey = k => String(k).replace(/[^A-Za-z0-9._-]+/g, '_');
const DEF = (() => {
  const store = {species: {}, family: {}, tree: {}, structure: {}}, wait = {};
  const dir = {species: 'species', family: 'families', tree: 'trees', structure: 'structures'};
  return {
    put(kind, key, obj) { store[kind][key] = obj;
      const w = wait[kind + '/' + key]; delete wait[kind + '/' + key];
      (w || []).forEach(f => f[0](obj)); },
    get: (kind, key) => store[kind][key],
    load(kind, key) {
      if (store[kind][key]) return Promise.resolve(store[kind][key]);
      const id = kind + '/' + key;
      return new Promise((res, rej) => {
        if (wait[id]) { wait[id].push([res, rej]); return; }
        wait[id] = [[res, rej]];
        const head = document.head || document.querySelector('head');
        if (!head || !head.appendChild) { delete wait[id]; rej(new Error('no document')); return; }
        const sc = document.createElement('script');
        sc.src = `data/${dir[kind]}/${safeKey(key)}.js`;
        sc.onerror = () => { const w = wait[id]; delete wait[id]; sc.remove();
          (w || []).forEach(f => f[1](new Error('not found: ' + sc.src))); };
        head.appendChild(sc);
      });
    }};
})();
window.DEF = DEF;

/* ---- names, keys and links ------------------------------------------- */
const SPKEY = {}, KEYSP = {};
D0.species.forEach((s, i) => { const k = s.key || safeKey(s.name);
  SPKEY[s.name] = k; KEYSP[k] = i; });
const spHref = s => '#/species/' + encodeURIComponent(SPKEY[s] || safeKey(s));
const famHref = f => '#/family/' + encodeURIComponent(f);
const geneHref = (s, pid) => `#/gene/${encodeURIComponent(SPKEY[s] || safeKey(s))}/${encodeURIComponent(pid)}`;
const treeHref = (k, q) => '#/trees/' + encodeURIComponent(k) + (q ? '?q=' + encodeURIComponent(q) : '');
const link = (href, text, cls = 'lk') => el('a', {href, class: cls}, text);
const spLink = s => { const a = link(spHref(s), ''); a.append(el('i', {}, sciName(s)));
  return a; };
const extLink = (href, text) => el('a', {href, target: '_blank', rel: 'noopener', class: 'ext'}, text + ' ↗');
const binomial = s => sciName(s).split(' ').slice(0, 2).join(' ');
const interpro = pf => `https://www.ebi.ac.uk/interpro/entry/pfam/${pf}/`;
const pfamLinks = pf => String(pf || '').split(',').filter(x => x).map((x, i) =>
  [i ? ', ' : '', el('a', {href: interpro(x.trim()), target: '_blank', rel: 'noopener',
    class: 'lk mono'}, x.trim())]).flat();
const STATUS_COL = {COMPLETE: '#2a9d8f', PARTIAL: '#e9a13b', FRAGMENT: '#bc4749', MISSING_DOMAIN: '#6c757d'};
const evBadge = ev => ev ? el('span', {class: 'badge', title: 'found below Pfam’s gathering threshold, passed the family’s checks',
  style: 'background:#eef0f3;color:#5a6470;margin-left:4px'}, 'below threshold') : null;
const statusBadge = st => !st ? '' : el('span', {class: 'badge',
  style: `background:${STATUS_COL[st] || '#adb5bd'}22;color:${STATUS_COL[st] || '#555'}`},
  st.toLowerCase().replace('_', ' '));
const RANKS = ['superkingdom', 'kingdom', 'phylum', 'subphylum', 'class', 'subclass', 'order',
  'suborder', 'superfamily', 'family', 'subfamily', 'genus'];
/* Plain names for the default map's families; a custom map falls back to the code. */
const FAMDESC = {CYP: 'cytochrome P450s', CCE: 'carboxyl/cholinesterases', FMO: 'flavin monooxygenases',
  ALDH: 'aldehyde dehydrogenases', AKR: 'aldo-keto reductases', SDR: 'short-chain dehydrogenases/reductases',
  EH: 'epoxide hydrolases and α/β-hydrolases', Cyt_b5: 'cytochrome b5', GST: 'glutathione S-transferases',
  GST_MAPEG: 'microsomal glutathione transferases', UGT: 'UDP-glycosyltransferases', SULT: 'sulfotransferases',
  ABC_full: 'ABC transporters', ABC_NBD: 'ATP-binding cassette domains', MATE: 'multidrug and toxin extrusion transporters',
  MFS: 'major facilitator transporters', Catalase: 'catalases', SOD_CuZn: 'copper-zinc superoxide dismutases',
  SOD_Fe: 'iron/manganese superoxide dismutases', Peroxiredoxin: 'peroxiredoxins', GPx: 'glutathione peroxidases',
  GCL: 'glutamate-cysteine ligase', Thioredoxin: 'thioredoxins', Glutaredoxin: 'glutaredoxins',
  TrxR_GR: 'thioredoxin and glutathione reductases', Hsp70: 'heat shock proteins, 70 kDa', Hsp90: 'heat shock proteins, 90 kDa',
  Hsp20: 'small heat shock proteins', DnaJ: 'DnaJ co-chaperones', Metallothionein: 'metallothioneins', Ferritin: 'ferritins',
  CDF_ZnT: 'zinc exporters', ZIP: 'zinc and iron importers', BTB_BACK_Kelch: 'Kelch-like proteins',
  CncC_bZIP: 'bZIP transcription factors', NR_xeno: 'nuclear receptors', bHLH_PAS: 'bHLH-PAS sensors'};
const RANK_COLS = TAX_COLS.filter(c => RANKS.includes(c));
/* trait columns (habitat, diet...): every grouping column that is not a rank */
const TRAIT_COLS = [...new Set([...TAX_COLS.filter(c => !RANKS.includes(c)),
  ...(D.trait_name && !RANKS.includes(D.trait_name) && D.trait_name !== 'method' ? [D.trait_name] : [])])];
const kb = n => n > 1e6 ? (n / 1e6).toFixed(1) + ' MB' : Math.max(1, Math.round(n / 1e3)) + ' kB';
function crumb(...parts) {
  const c = el('div', {class: 'crumb'});
  parts.filter(x => x).forEach((x, i) => { if (i) c.append(el('span', {class: 'sep'}, '›'));
    c.append(x); });
  return c;
}
/* Clicking a rank anywhere means "show me the species in this taxon". */
function showTaxon(rank, value) {
  FILTER = {[rank]: value}; buildFilterBar(); rerender(false);
  location.hash = '#/species';
}
const taxLink = (rank, value) => el('a', {href: '#/species', title: `all species with ${rank} ${value}`,
  onclick: e => { e.preventDefault(); showTaxon(rank, value); }}, nice(value));

/* ---- a paged, sortable data grid ------------------------------------ */
/* cols: {label, get(row) for sorting and TSV, cell(row) -> node|string,
   num: right-aligned}. One page at a time: a 1,575-row DOM table is slow and
   nobody reads it top to bottom. */
function dataGrid(host, cols, rows, opts = {}) {
  let sortI = opts.sortBy ?? 0, dir = opts.desc ? -1 : 1, q = '', page = 0;
  const size = opts.pageSize || 50;
  const box = el('div');
  const search = el('input', {type: 'search', placeholder: opts.placeholder || 'filter rows…',
    oninput: e => { q = e.target.value.toLowerCase(); page = 0; draw(); }});
  const count = el('span', {class: 'hint', style: 'margin:0'});
  const tools = el('div', {class: 'row', style: 'align-items:center;margin-bottom:8px'},
    search, count, ...(opts.tools || []),
    el('button', {class: 'act', onclick: () => dl((opts.name || 'table') + '.tsv',
      [cols.map(c => c.label).join('\t'),
       ...view().map(r => cols.map(c => { const v = c.get(r); return v == null ? '' : v; }).join('\t'))]
        .join('\n'))}, 'TSV'));
  const wrap = el('div', {style: 'overflow:auto;border:1px solid var(--line);border-radius:6px'});
  const pager = el('div', {class: 'pager'});
  box.append(tools, wrap, pager); host.append(box);
  const text = r => (opts.text ? opts.text(r) : cols.map(c => c.get(r)).join(' ')).toLowerCase();
  const idx = rows.map(text);
  function view() {
    let r = rows.map((x, i) => i);
    if (q) r = r.filter(i => idx[i].includes(q));
    const c = cols[sortI];
    return r.map(i => rows[i]).sort((a, b) => { const x = c.get(a), y = c.get(b);
      if (x == null || x === '') return 1; if (y == null || y === '') return -1;
      if (typeof x === 'number' && typeof y === 'number') return (x - y) * dir;
      return String(x).localeCompare(String(y)) * dir; });
  }
  function draw() {
    const v = view(), pages = Math.max(1, Math.ceil(v.length / size));
    page = Math.min(page, pages - 1);
    count.textContent = `${v.length.toLocaleString()} ${v.length === 1 ? (opts.unit || 'row') : (opts.units || 'rows')}`;
    const t = el('table', {class: 'dg'},
      el('thead', {}, el('tr', {}, cols.map((c, i) => el('th', {class: c.num ? 'num' : '',
        onclick: () => { if (sortI === i) dir = -dir; else { sortI = i; dir = c.num ? -1 : 1; } draw(); }},
        c.label + (sortI === i ? (dir > 0 ? ' ▲' : ' ▼') : ''))))),
      el('tbody', {}, v.slice(page * size, (page + 1) * size).map(r => el('tr', {},
        cols.map(c => { const x = c.cell ? c.cell(r) : c.get(r);
          return el('td', {class: (c.num ? 'num ' : '') + (c.cls || '')},
            x == null ? '' : (x.nodeType ? x : Array.isArray(x) ? x : (typeof x === 'number' ? fmt(x) : String(x)))); })))));
    wrap.innerHTML = ''; wrap.append(t);
    pager.innerHTML = '';
    if (pages > 1) pager.append(
      el('button', {class: 'act', onclick: () => { page = Math.max(0, page - 1); draw(); }}, '‹ prev'),
      el('span', {}, `page ${page + 1} of ${pages}`),
      el('button', {class: 'act', onclick: () => { page = Math.min(pages - 1, page + 1); draw(); }}, 'next ›'));
  }
  draw();
  return {draw, setQuery: x => { q = x.toLowerCase(); search.value = x; page = 0; draw(); },
          view};
}

/* ---- domain architecture, to scale --------------------------------- */
const PFCOL = {};
const pfColour = id => PFCOL[id] || (PFCOL[id] =
  `hsl(${(Object.keys(PFCOL).length * 47 + 200) % 360},52%,${48 + (Object.keys(PFCOL).length % 3) * 7}%)`);
function archSVG(len, doms, W = 220, H = 14, labels = false) {
  const h = labels ? H + 16 : H, sv = svgEl('svg', {width: W, height: h, viewBox: `0 0 ${W} ${h}`});
  const x = v => 1 + (W - 2) * (v / Math.max(len, 1));
  const mid = H / 2;
  sv.append(svgEl('line', {x1: 1, x2: W - 1, y1: mid, y2: mid, stroke: '#9aa1a8', 'stroke-width': 2}));
  (doms || []).forEach(d => { const [name, id, a, b, hf, ht, hl, sc] = d;
    const r = svgEl('rect', {x: x(a), y: 1, width: Math.max(2, x(b) - x(a)), height: H - 2, rx: 3,
      fill: pfColour(id), stroke: '#fff', 'stroke-width': .6});
    r.addEventListener('mousemove', e => showTip(e, `<b>${name}</b> ${id}<br>residues ${a}–${b}<br>` +
      `HMM ${hf}–${ht} of ${hl} (${Math.round((ht - hf + 1) / hl * 100)}%)<br>bit score ${sc}`));
    r.addEventListener('mouseleave', hideTip);
    sv.append(r);
    if (labels && x(b) - x(a) > 24) { const t = svgEl('text', {x: (x(a) + x(b)) / 2, y: H + 12,
      'text-anchor': 'middle', 'font-size': 10.5, fill: '#3d4852'}); t.textContent = name; sv.append(t); }
  });
  return sv;
}

/* ---- species list ---------------------------------------------------- */
function tabSpeciesList(host) {
  const coreF = FAM;
  const tot = i => coreF.reduce((a, f) => a + D.counts[f][i], 0);
  const p10 = i => D.per10k ? coreF.reduce((a, f) => a + D.per10k[f][i], 0) : null;
  const ranks = RANK_COLS.filter(r => ['class', 'order', 'family'].includes(r));
  const shownRanks = ranks.length ? ranks : RANK_COLS.slice(-3);
  const gcol = D.trait_name && !shownRanks.includes(D.trait_name) ? D.trait_name : null;
  const anyCommon = SP.some(s => commonName(s));
  const rows = SP.map((s, i) => ({s, i}));
  const cols = [
    {label: 'species', get: r => sciName(r.s), cell: r => spLink(r.s)},
    ...(anyCommon ? [{label: 'common name', get: r => commonName(r.s)}] : []),
    ...shownRanks.map(k => ({label: k, get: r => (D.species[r.i].tax || {})[k] || ''})),
    ...(gcol ? [{label: gcol, get: r => D.species[r.i].trait ?? ''}] : []),
    {label: 'defensome genes', num: true, get: r => tot(r.i)},
    ...(D.per10k ? [{label: `per 10k ${D.size_unit || 'proteins'}`, num: true, get: r => +p10(r.i).toFixed(1)}] : []),
    {label: `proteome (${D.size_unit || 'proteins'})`, num: true, get: r => D.species[r.i].n_proteins,
     cell: r => D.species[r.i].n_proteins.toLocaleString()},
  ];
  const p = el('div', {class: 'panel'});
  const nf = Object.keys(FILTER).length;
  p.append(el('h2', {}, 'Species'));
  tabSpeciesList._grid = dataGrid(p, cols, rows, {sortBy: 0, name: 'species', unit: 'species', units: 'species',
    placeholder: 'filter by name, taxon…'});
  host.append(p);
}

/* ---- species page --------------------------------------------------- */
function speciesPage(host, key) {
  const i0 = KEYSP[key];
  if (i0 === undefined) { host.append(el('div', {class: 'panel'}, el('h2', {}, 'Species not found'),
    el('p', {}, `No species with key “${key}”.`), link('#/species', 'Browse all species'))); return; }
  const s0 = D0.species[i0], name = s0.name;
  const tax = s0.tax || {};
  const fams = FAM, cnt = f => D0.counts[f][i0];
  const total = fams.reduce((a, f) => a + cnt(f), 0);
  const present = fams.filter(f => cnt(f) > 0).length;
  const per10 = D0.per10k ? fams.reduce((a, f) => a + D0.per10k[f][i0], 0) : null;

  const head = el('div', {class: 'ehead'});
  head.append(crumb(link('#/species', 'Species'), ...RANK_COLS.filter(r => tax[r]).map(r => taxLink(r, tax[r]))));
  head.append(el('h1', {}, el('i', {}, sciName(name))));
  const sub = [commonName(name), qualifier(name)].filter(x => x).join(' · ');
  head.append(el('div', {class: 'sub'}, sub ? sub + ' · ' : '', el('span', {class: 'mono'}, name)));
  const meta = el('div', {class: 'meta'});
  TAX_COLS.filter(c => !RANKS.includes(c) && tax[c]).forEach(c =>
    meta.append(el('span', {class: 'tag', title: c}, `${c}: ${nice(tax[c])}`)));
  const bn = encodeURIComponent(binomial(name));
  meta.append(extLink(`https://en.wikipedia.org/wiki/Special:Search?search=${bn}`, 'Wikipedia'),
    extLink(`https://www.gbif.org/species/search?q=${bn}`, 'GBIF'),
    extLink(`https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?name=${bn}`, 'NCBI Taxonomy'),
    extLink(`https://goat.genomehubs.org/search?query=tax_name(${bn})&result=taxon`, 'GoaT genome'),
    extLink(`https://www.ebi.ac.uk/ena/browser/text-search?query=${bn}`, 'ENA'));
  head.append(meta);
  head.append(el('div', {class: 'jump'}, link('#sp-profile', 'Defensome profile', ''),
    link('#sp-genes', 'Genes', '')));
  head.querySelectorAll('.jump a').forEach(a => a.onclick = e => { e.preventDefault();
    const t = $(a.getAttribute('href')); if (t && t.scrollIntoView) t.scrollIntoView({behavior: 'smooth'}); });
  host.append(head);

  const kp = el('div', {class: 'grid g4'},
    kpi(total.toLocaleString(), 'defensome genes'),
    kpi(per10 == null ? null : per10.toFixed(1), `per 10k ${D.size_unit || 'proteins'}`),
    kpi(`${present} of ${fams.length}`, 'families present'),
    kpi(s0.n_proteins ? s0.n_proteins.toLocaleString() : null, `${D.size_unit || 'proteins'} in the proteome`));
  host.append(kp);

  /* profile: this species against the median of the species in view */
  const p = el('div', {class: 'panel', id: 'sp-profile'});
  const kind = D0.per10k ? 'per10k' : 'raw';
  const val = f => kind === 'per10k' ? D0.per10k[f][i0] : D0.counts[f][i0];
  const med = f => stats(col(kind, f)).med;
  p.append(el('h2', {}, 'Defensome profile'), el('div', {class: 'legend', style: 'margin:0 0 8px'},
    el('span', {}, 'bar: this species'),
    el('span', {}, el('i', {class: 'sw', style: 'background:#1c1f23;width:2px'}), 'typical across species'),
    el('span', {style: 'color:var(--mut)'}, kind === 'per10k' ? `copies per 10k ${D.size_unit || 'proteins'}` : 'copies')));
  p.append(figure(profileChart(fams.map(f => ({f, v: val(f), m: med(f), raw: cnt(f)}))),
    safeKey(name) + '_profile'));
  host.append(p);

  /* genes */
  const g = el('div', {class: 'panel', id: 'sp-genes'});
  g.append(el('h2', {}, 'Genes'));
  host.append(g);
  if (!SITE.genes) { g.append(el('p', {class: 'empty'},
    'Gene-level data was not exported with this dashboard. Regenerate it with `defensome.py dashboard`.')); return; }
  const wait = el('p', {class: 'empty'}, 'Loading genes…'); g.append(wait);
  DEF.load('species', key).then(d => { wait.remove(); speciesGenes(g, name, d); })
    .catch(() => { wait.textContent = 'Gene data not found. Keep the data/ folder beside index.html.'; });
}
function speciesGenes(host, name, d) {
  const famSel = el('select', {}, el('option', {value: ''}, 'all families'),
    ...[...new Set(d.genes.map(x => x[1]))].sort().map(f => el('option', {value: f},
      `${f} (${d.genes.filter(x => x[1] === f).length})`)));
  host.append(el('p', {class: 'count'}, `${d.genes.length.toLocaleString()} genes`));
  const fasta = rows => rows.filter(r => d.seq[r[0]]).map(r => `>${name}|${r[0]} ${r[1]}\n${wrap60(d.seq[r[0]])}`).join('\n') + '\n';
  const cols = [
    {label: 'protein', get: r => r[0], cell: r => link(geneHref(name, r[0]), r[0], 'lk mono')},
    {label: 'family', get: r => r[1], cell: r => link(famHref(r[1]), r[1])},
    {label: 'length', num: true, get: r => r[2]},
    {label: 'domains', get: r => (d.domains[r[0]] || []).map(x => x[0]).join('-'),
     cell: r => archSVG(r[2], d.domains[r[0]], 200, 12)},
    {label: 'status', get: r => r[4] + (r[6] ? ' below threshold' : ''), cell: r => el('span', {}, statusBadge(r[4]), evBadge(r[6]))},
    ...(d.genes.some(r => r[5]) ? [{label: 'CYP clan', get: r => r[5]}] : []),
  ];
  let rows = d.genes;
  const box = el('div');
  const draw = () => { box.innerHTML = '';
    rows = famSel.value ? d.genes.filter(r => r[1] === famSel.value) : d.genes;
    speciesGenes._grid = dataGrid(box, cols, rows, {sortBy: 1, pageSize: 20, name: safeKey(name) + '_genes', unit: 'gene', units: 'genes',
      placeholder: 'filter by protein ID…',
      tools: SITE.sequences ? [el('button', {class: 'act', onclick: () =>
        dlText(`${safeKey(name)}${famSel.value ? '_' + famSel.value : ''}.faa`, fasta(rows))}, 'FASTA')] : []}); };
  famSel.onchange = draw;
  host.append(el('div', {class: 'row'}, el('div', {class: 'ctl'}, el('label', {}, 'family'), famSel)), box);
  draw();
}
const wrap60 = s => s.replace(/(.{60})/g, '$1\n').trim();
const dlText = (name, text) => { const b = new Blob([text], {type: 'text/plain'});
  const u = URL.createObjectURL(b), a = el('a', {href: u, download: name}); a.click();
  setTimeout(() => URL.revokeObjectURL(u), 1000); };

/* One row per family: a bar for this species and a tick for the median of the
   species in view, coloured by defence category. */
function profileChart(items) {
  const W = 880, rh = 18, ML = 150, MR = 70, MT = 6, H = MT + items.length * rh + 8;
  const mx = Math.max(...items.map(d => Math.max(d.v, d.m)), 1e-9);
  const x = v => ML + (W - ML - MR) * v / mx;
  const sv = svgEl('svg', {width: W, height: H, viewBox: `0 0 ${W} ${H}`});
  items.forEach((d, k) => { const y = MT + k * rh;
    const lab = svgEl('text', {x: ML - 8, y: y + rh / 2 + 4, 'text-anchor': 'end', 'font-size': 11.5,
      fill: d.raw ? '#1b4965' : '#9aa1a8', style: 'cursor:pointer'});
    lab.textContent = d.f; lab.addEventListener('click', () => { location.hash = famHref(d.f); });
    sv.append(lab);
    const c = catColour(FMETA[d.f] ? FMETA[d.f].category : CATS[0]);
    const r = svgEl('rect', {x: ML, y: y + 3, width: Math.max(d.v ? 1.5 : 0, x(d.v) - ML), height: rh - 6,
      fill: c, opacity: .85});
    const tipf = e => showTip(e, `<b>${d.f}</b> (${FMETA[d.f] ? FMETA[d.f].category : ''})<br>` +
      `this species ${fmt(d.v, 2)} (${d.raw} genes)<br>median in view ${fmt(d.m, 2)}`);
    r.addEventListener('mousemove', tipf); r.addEventListener('mouseleave', hideTip); sv.append(r);
    sv.append(svgEl('line', {x1: x(d.m), x2: x(d.m), y1: y + 1, y2: y + rh - 1, stroke: '#1c1f23',
      'stroke-width': 1.6}));
    const t = svgEl('text', {x: Math.max(x(d.v), x(d.m)) + 6, y: y + rh / 2 + 4, 'font-size': 10.5, fill: '#6b7280'});
    t.textContent = d.raw; sv.append(t);
  });
  return sv;
}

/* ---- gene page ------------------------------------------------------ */
function genePage(host, key, pid) {
  const i0 = KEYSP[key];
  if (i0 === undefined) { host.append(el('div', {class: 'panel'}, el('h2', {}, 'Gene not found'))); return; }
  const name = D0.species[i0].name;
  const wait = el('div', {class: 'panel'}, el('p', {class: 'empty'}, 'Loading…'));
  host.append(wait);
  DEF.load('species', key).then(d => { wait.remove(); geneBody(host, name, d, pid); })
    .catch(() => { wait.innerHTML = '';
      wait.append(el('h2', {}, pid), el('p', {class: 'empty'},
        'Gene data not found. Keep the data/ folder beside index.html, or regenerate with `defensome.py dashboard`.')); });
}
function geneBody(host, name, d, pid) {
  const rows = d.genes.filter(r => r[0] === pid);
  if (!rows.length) { host.append(el('div', {class: 'panel'}, el('h2', {}, 'Gene not found'),
    el('p', {}, `${pid} is not among the defensome calls of `, el('i', {}, sciName(name)), '.'),
    link(spHref(name), 'Back to the species'))); return; }
  const [, fam, len, cov, st, clan, ev] = rows[0];
  const doms = d.domains[pid] || [], seq = d.seq[pid] || '';
  const fams = rows.map(r => r[1]);
  const head = el('div', {class: 'ehead'});
  head.append(crumb(link('#/species', 'Species'), spLink(name), link(famHref(fam), fam), el('span', {}, pid)));
  head.append(el('h1', {class: 'mono', style: 'font-size:21px'}, pid));
  head.append(el('div', {class: 'sub'}, fams.join(' + ') + ' gene in ', el('i', {}, sciName(name)),
    commonName(name) ? ` (${commonName(name)})` : '', ` · ${len.toLocaleString()} aa`));
  const meta = el('div', {class: 'meta'}, statusBadge(st), evBadge(ev), clan ? el('span', {class: 'tag'}, 'CYP clan ' + clan) : null,
    ...[...new Set(doms.map(x => x[1]))].map(id => extLink(interpro(id), 'InterPro ' + id)));
  if (seq) meta.append(
    extLink('https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&PAGE_TYPE=BlastSearch&DATABASE=swissprot&QUERY=' +
      encodeURIComponent(seq), 'closest characterised proteins'),
    extLink('https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&PAGE_TYPE=BlastSearch&QUERY=' +
      encodeURIComponent(seq), 'BLASTP'));
  head.append(meta); host.append(head);


  const a = el('div', {class: 'panel'});
  a.append(el('h2', {}, 'Domain architecture'), el('p', {class: 'hint'},
    'To scale. Hover a domain for its coordinates and HMM coverage; only defensome Pfam models were searched.'));
  a.append(figure(archSVG(len, doms, 880, 22, true), safeKey(pid) + '_domains'));
  if (doms.length) {
    const t = el('table', {class: 'dg'}, el('thead', {}, el('tr', {},
      ['domain', 'Pfam', 'residues', 'HMM span', 'HMM coverage', 'bit score'].map(h => el('th', {}, h)))),
      el('tbody', {}, doms.map(x => el('tr', {}, el('td', {}, el('i', {class: 'sw',
        style: `background:${pfColour(x[1])};margin-right:6px`}), x[0]),
        el('td', {}, el('a', {href: interpro(x[1]), target: '_blank', rel: 'noopener', class: 'lk mono'}, x[1])),
        el('td', {}, `${x[2]}–${x[3]}`), el('td', {}, `${x[4]}–${x[5]} of ${x[6]}`),
        el('td', {}, Math.round((x[5] - x[4] + 1) / x[6] * 100) + '%'), el('td', {}, x[7])))));
    a.append(el('div', {style: 'margin-top:10px'}, t));
  }
  host.append(a);

  const q = el('div', {class: 'panel'});
  q.append(el('h2', {}, 'Sequence'));
  if (!seq) q.append(el('p', {class: 'empty'}, SITE.sequences
    ? 'This protein was not found in the proteome FASTA (header and ID mismatch?).'
    : 'Sequences were not exported (dashboard --no-sequences).'));
  else {
    const fa = `>${name}|${pid} ${fams.join('+')}\n${wrap60(seq)}\n`;
    q.append(el('div', {class: 'row'},
      el('button', {class: 'act', onclick: e => { (navigator.clipboard ? navigator.clipboard.writeText(fa) : Promise.reject())
        .then(() => { e.target.textContent = 'copied'; }).catch(() => {}); }}, 'copy FASTA'),
      el('button', {class: 'act', onclick: () => dlText(safeKey(pid) + '.faa', fa)}, 'download FASTA'),
      el('span', {class: 'hint', style: 'margin:0'}, `${seq.length} aa`)));
    q.append(el('div', {class: 'seq'}, ...formatSeq(seq)));
  }
  host.append(q);

  host.append(structurePanel(name, pid, fams, seq));

  /* context: paralogues here, and trees that may hold this gene */
  const c = el('div', {class: 'panel'});
  const par = d.genes.filter(r => fams.includes(r[1]) && r[0] !== pid);
  c.append(el('h2', {}, `Other ${fams.join('/')} genes in `, el('i', {}, sciName(name)), ` (${par.length})`));
  if (par.length) dataGrid(c, [
    {label: 'protein', get: r => r[0], cell: r => link(geneHref(name, r[0]), r[0], 'lk mono')},
    {label: 'family', get: r => r[1]}, {label: 'length', num: true, get: r => r[2]},
    {label: 'domains', get: r => (d.domains[r[0]] || []).map(x => x[0]).join('-'),
     cell: r => archSVG(r[2], d.domains[r[0]], 200, 12)},
    {label: 'status', get: r => r[4], cell: r => statusBadge(r[4])}], par,
    {sortBy: 2, desc: true, pageSize: 15, name: safeKey(pid) + '_paralogues', unit: 'gene', units: 'genes'});
  const trees = Object.entries(D.tree_meta || {}).filter(([, m]) => fams.includes(m.family));
  if (trees.length) {
    c.append(el('h2', {style: 'margin-top:18px'}, 'Trees'));
    const tip = `${name}|${pid}`;
    trees.forEach(([k, m]) => { const nwk = D.trees[k];
      const inTree = typeof nwk === 'string' ? nwk.includes(tip) : null;
      c.append(el('div', {style: 'padding:4px 0'}, link(treeHref(k, pid), treeLabel(k)),
        el('span', {class: 'hint'}, inTree === false ? ' · this gene is not in it (subset or subsample)'
          : inTree ? ' · contains this gene' : ' · opens and searches for this gene'))); });
  }
  host.append(c);
}
function formatSeq(seq) {
  const out = [];
  for (let i = 0; i < seq.length; i += 60) {
    const line = seq.slice(i, i + 60).replace(/(.{10})/g, '$1 ').trim();
    out.push(el('span', {class: 'pos'}, String(i + 1).padStart(5, ' ') + '  '), line + '\n');
  }
  return out;
}
function treeLabel(k) {
  const m = (D.tree_meta || {})[k] || {};
  if (k === 'species') return `species tree (${m.n_tips || '?'} tips)`;
  const bits = [];
  Object.entries(m.where || {}).forEach(([c, v]) => bits.push(`${c} ${[].concat(v).join(' or ')}`));
  if (m.clan && m.clan.length) bits.push('clan ' + m.clan.join('+'));
  if (m.complete_only) bits.push('complete models');
  if (m.per_species) bits.push(`≤${m.per_species} per species`);
  if (m.thin) { const [c, k] = String(m.thin).split(':'); bits.push(`${k || 1} species per ${c}`); }
  if (!bits.length && k.includes('__')) bits.push(k.split('__').slice(1).join(' '));
  return `${m.family || k}${bits.length ? ': ' + bits.join(', ') : ''} (${(m.n_tips || 0).toLocaleString()} tips)`;
}

/* ---- structure and binding ------------------------------------------ */
/* Predictions already made for this gene are shown; anything else can be
   set up here and run with Boltz-2 (`structure`) or AlphaFold Server. */
const COFACTOR_NAME = {HEM: 'haem', FAD: 'FAD', NAP: 'NADP+', NAD: 'NAD+', GSH: 'glutathione',
  UPG: 'UDP-glucose', ATP: 'ATP', ADP: 'ADP'};
const AF_SERVER_LIGANDS = new Set(['HEM', 'FAD', 'NAP', 'NAD', 'ATP', 'ADP']);
function structurePanel(name, pid, fams, seq) {
  const p = el('div', {class: 'panel', id: 'gene-structure'});
  p.append(el('h2', {}, 'Structure'));
  const preds = ((D.structures || {})[`${name}|${pid}`] || []);
  const done = preds.filter(x => x.predicted);
  if (done.length) {
    const view = el('div', {class: 'mol', style: 'display:none'});
    const cards = el('div', {class: 'cards', style: 'margin-bottom:10px'});
    done.forEach(x => {
      const bits = [];
      if (x.ptm != null) bits.push(`pTM ${x.ptm.toFixed(2)}`);
      if (x.complex_plddt != null) bits.push(`pLDDT ${Math.round(x.complex_plddt * 100)}`);
      const bind = x.affinity_probability != null
        ? `binder ${Math.round(x.affinity_probability * 100)}% · IC50 ~${fmtIC50(x.affinity_log10_ic50_uM)}` : '';
      cards.append(el('a', {href: '#', class: 'card', style: 'flex-direction:column;align-items:flex-start;gap:3px',
        onclick: e => { e.preventDefault(); showStructure(view, x); }},
        el('b', {}, x.ligand_name ? `with ${x.ligand_name}` : 'apo'),
        el('span', {class: 'n'}, [x.cofactors.map(c => COFACTOR_NAME[c] || c).join(' + '), bits.join(' · ')].filter(t => t).join(' · ')),
        bind ? el('span', {style: 'font-size:12.5px'}, bind) : null));
    });
    p.append(cards, view);
    structurePanel._first = () => showStructure(view, done[0]);
  }
  if (!seq) { p.append(el('p', {class: 'empty'}, 'No sequence exported.')); return p; }

  /* set up a new prediction */
  const presets = [...new Set(fams.flatMap(f => (D.cofactors || {})[f] || []))];
  const all = [...new Set([...presets, 'HEM', 'FAD', 'NAP', 'NAD', 'GSH', 'UPG', 'ATP'])];
  const cof = new Set(presets);
  const chips = el('div', {class: 'chips'}, ...all.map(c => el('span', {class: 'chip' + (cof.has(c) ? ' on' : ''),
    title: c, onclick: e => { cof.has(c) ? cof.delete(c) : cof.add(c); e.target.classList.toggle('on'); }},
    COFACTOR_NAME[c] || c)));
  const lig = el('input', {type: 'text', placeholder: 'ligand SMILES (optional)', style: 'min-width:320px;flex:1'});
  const ligName = el('input', {type: 'text', placeholder: 'ligand name', style: 'min-width:140px;width:160px'});
  const stem = () => `${safeKey(name)}__${safeKey(pid)}${ligName.value.trim() ? '__' + safeKey(ligName.value.trim()) : ''}`;
  const yaml = () => { const ids = 'BCDEFGHIJKLMNOPQRSTUVWXYZ'.split(''); let k = 0;
    const L = ['version: 1', 'sequences:', '  - protein:', '      id: A', `      sequence: ${seq.replace(/U/g, 'C')}`];
    [...cof].forEach(c => L.push('  - ligand:', `      id: ${ids[k++]}`, `      ccd: ${c}`));
    const sm = lig.value.trim();
    if (sm) { const id = ids[k++]; L.push('  - ligand:', `      id: ${id}`, `      smiles: '${sm}'`,
      'properties:', '  - affinity:', `      binder: ${id}`); }
    return L.join('\n') + '\n'; };
  const afjson = () => JSON.stringify([{name: stem(), modelSeeds: [], dialect: 'alphafoldserver', version: 1,
    sequences: [{proteinChain: {sequence: seq.replace(/U/g, 'C'), count: 1}},
      ...[...cof].filter(c => AF_SERVER_LIGANDS.has(c)).map(c => ({ligand: {ligand: 'CCD_' + c, count: 1}}))]}], null, 1);
  const cmd = () => `python3 defensome/defensome.py structure --out ${D.out_dir || 'results'} --gene '${name}|${pid}'` +
    (cof.size ? ` --cofactor ${[...cof].join(',')}` : ' --cofactor none') +
    (lig.value.trim() ? ` --ligand '${lig.value.trim()}' --ligand-name ${safeKey(ligName.value.trim() || 'ligand')}` : '') + ' --run';
  const cmdBox = el('div', {class: 'mono', style: 'display:none;white-space:pre-wrap;background:#f6f8fa;padding:10px;border-radius:6px;margin-top:8px'});
  const pub = el('a', {href: 'https://pubchem.ncbi.nlm.nih.gov/', target: '_blank', rel: 'noopener', class: 'ext'}, 'PubChem ↗');
  ligName.addEventListener('input', () => { const q = ligName.value.trim();
    pub.setAttribute('href', q ? `https://pubchem.ncbi.nlm.nih.gov/#query=${encodeURIComponent(q)}` : 'https://pubchem.ncbi.nlm.nih.gov/'); });
  const box = el('div', {});
  box.append(el('div', {class: 'row', style: 'align-items:center;margin-bottom:8px'}, el('span', {class: 'count', style: 'margin:0'}, 'with'), chips),
    el('div', {class: 'row', style: 'align-items:center'}, ligName, lig, pub),
    el('div', {class: 'row', style: 'margin:4px 0 0'},
      el('button', {class: 'act', onclick: () => dlText(stem() + '.yaml', yaml())}, 'Boltz input'),
      el('button', {class: 'act', onclick: () => dlText(stem() + '.json', afjson())}, 'AlphaFold Server input'),
      el('button', {class: 'act', onclick: () => { cmdBox.textContent = cmd(); cmdBox.style.display = 'block'; }}, 'command'),
      el('a', {class: 'ext', href: 'https://alphafoldserver.com/', target: '_blank', rel: 'noopener'}, 'AlphaFold Server ↗'),
      el('a', {class: 'ext', href: 'https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb',
        target: '_blank', rel: 'noopener'}, 'ColabFold ↗')), cmdBox);
  p.append(done.length ? collapsible('predict another', box) : box);
  structurePanel._yaml = yaml; structurePanel._af = afjson; structurePanel._cmd = cmd;
  structurePanel._set = (smiles, nm) => { lig.value = smiles; ligName.value = nm; };
  return p;
}
const fmtIC50 = v => { if (v == null) return '?'; const x = Math.pow(10, v);
  return x >= 1000 ? `${(x / 1000).toPrecision(2)} mM` : x >= 1 ? `${x.toPrecision(2)} µM` : `${(x * 1000).toPrecision(2)} nM`; };
/* 3Dmol is fetched only when a structure is opened, so the site stays
   offline-capable for everything else. */
const MOL_URL = 'https://cdn.jsdelivr.net/npm/3dmol@2.4.0/build/3Dmol-min.js';
let MOL_P = null;
const loadMol = () => MOL_P || (MOL_P = new Promise((res, rej) => {
  if (window.$3Dmol) return res(window.$3Dmol);
  const sc = document.createElement('script'); sc.src = MOL_URL;
  sc.onload = () => window.$3Dmol ? res(window.$3Dmol) : rej(new Error('3Dmol'));
  sc.onerror = () => { MOL_P = null; rej(new Error('offline')); };
  document.head.appendChild(sc); }));
function showStructure(view, x) {
  view.style.display = 'block'; view.innerHTML = '';
  const st = el('p', {class: 'empty'}, 'Loading structure…'); view.append(st);
  DEF.load('structure', x.id).then(o => {
    const save = el('button', {class: 'act', onclick: () => dlText(x.id + '.cif', o.cif)}, 'download mmCIF');
    loadMol().then(M => {
      view.innerHTML = '';
      const stage = el('div', {style: 'width:100%;height:460px;position:relative;border:1px solid var(--line);border-radius:6px'});
      view.append(stage, el('div', {class: 'row', style: 'margin-top:8px;align-items:center'}, save,
        el('span', {class: 'legend', style: 'margin:0'},
          ...[['#0053d6', 'very high'], ['#65cbf3', 'confident'], ['#ffdb13', 'low'], ['#ff7d45', 'very low']]
            .map(([c, t]) => el('span', {}, el('i', {class: 'sw', style: 'background:' + c}), t)))));
      const v = M.createViewer(stage, {backgroundColor: 'white'});
      v.addModel(o.cif, 'cif');
      const plddt = a => { const b = a.b > 1 ? a.b : a.b * 100;
        return b > 90 ? '#0053d6' : b > 70 ? '#65cbf3' : b > 50 ? '#ffdb13' : '#ff7d45'; };
      v.setStyle({hetflag: false}, {cartoon: {colorfunc: plddt}});
      v.setStyle({hetflag: true}, {stick: {colorscheme: 'greenCarbon', radius: .25}});
      v.zoomTo(); v.render();
      showStructure._viewer = v;
    }).catch(() => { view.innerHTML = '';
      view.append(el('p', {class: 'empty'}, 'The 3D viewer needs an internet connection.'), save); });
  }).catch(() => { st.textContent = 'Structure file not found beside index.html.'; });
}

/* ---- families list --------------------------------------------------- */
function tabFamilies(host) {
  const comp = D.completeness || {};
  const pct = f => { const c = comp[f]; if (!c) return null;
    const t = Object.values(c).reduce((a, b) => a + b, 0); return t ? +(100 * (c.COMPLETE || 0) / t).toFixed(1) : null; };
  const nSp = SP.length;
  const rows = FAM.map(f => { const v = D.counts[f], st = stats(v);
    return {f, m: FMETA[f], pres: v.filter(x => x > 0).length, med: st.med, max: st.max,
            tot: v.reduce((a, b) => a + b, 0), pc: pct(f)}; });
  const p = el('div', {class: 'panel'});
  p.append(el('h2', {}, 'Gene families'), el('p', {class: 'hint'},
    `${FAM.length} families in ${CATS.length} defence categories, counted over the ${nSp.toLocaleString()} species in view. Click a family for its distribution, members, domain architectures and trees.`));
  const bar = (frac, c) => el('span', {style: 'display:inline-flex;align-items:center;gap:6px'},
    el('span', {style: `display:inline-block;width:60px;height:7px;background:#edf0f3;border-radius:4px;overflow:hidden`},
      el('span', {style: `display:block;height:100%;width:${(frac * 100).toFixed(1)}%;background:${c}`})),
    (frac * 100).toFixed(0) + '%');
  dataGrid(p, [
    {label: 'family', get: r => r.f, cell: r => link(famHref(r.f), r.f)},
    {label: '', get: r => FAMDESC[r.f] || '', cls: 'muted'},
    {label: 'category', get: r => r.m.category, cell: r => el('span', {},
      el('i', {class: 'sw', style: `background:${catColour(r.m.category)};margin-right:6px;vertical-align:-1px`}), r.m.category)},
    {label: 'tier', get: r => r.m.tier},
    {label: 'Pfam', get: r => r.m.pfam, cell: r => el('span', {}, ...pfamLinks(r.m.pfam))},
    {label: 'species with ≥1', num: true, get: r => r.pres / Math.max(nSp, 1),
     cell: r => bar(r.pres / Math.max(nSp, 1), catColour(r.m.category))},
    {label: 'median copies', num: true, get: r => r.med},
    {label: 'max', num: true, get: r => r.max},
    {label: 'genes', num: true, get: r => r.tot, cell: r => r.tot.toLocaleString()},
    ...(Object.keys(comp).length ? [{label: '% complete', num: true, get: r => r.pc}] : []),
  ], rows, {sortBy: 1, pageSize: 100, name: 'families', unit: 'family', units: 'families',
    placeholder: 'filter by family, category or Pfam…'});
  host.append(p);
}

/* ---- family page additions: members, trees, downloads ---------------- */
function familyExtras(host, fam) {
  const m = el('div', {class: 'panel', id: 'fam-members'});
  m.append(el('h2', {}, 'Members'));
  host.append(m);
  if (!SITE.genes) m.append(el('p', {class: 'empty'}, 'Gene-level data was not exported with this dashboard.'));
  else {
    const wait = el('p', {class: 'empty'}, 'Loading members…'); m.append(wait);
    DEF.load('family', fam).then(d => { wait.remove();
      const inView = d.members.filter(r => SPI[D0.species[r[0]].name] !== undefined);
      m.append(el('p', {class: 'count'}, `${inView.length.toLocaleString()} genes in ${new Set(inView.map(r => r[0])).size.toLocaleString()} species`));
      const any = inView.some(r => r[4]);
      familyExtras._grid = dataGrid(m, [
        {label: 'species', get: r => sciName(D0.species[r[0]].name), cell: r => spLink(D0.species[r[0]].name)},
        {label: 'protein', get: r => r[1], cell: r => link(geneHref(D0.species[r[0]].name, r[1]), r[1], 'lk mono')},
        {label: 'length', num: true, get: r => r[2]},
        {label: 'status', get: r => r[3] + (r[5] ? ' below threshold' : ''), cell: r => el('span', {}, statusBadge(r[3]), evBadge(r[5]))},
        ...(any ? [{label: 'CYP clan', get: r => r[4]}] : []),
      ], inView, {sortBy: 0, pageSize: 20, name: fam + '_members', unit: 'gene', units: 'genes',
        placeholder: 'filter by species or protein…'});
    }).catch(() => { wait.textContent = 'Member data not found. Keep the data/ folder beside index.html.'; });
  }
  const t = el('div', {class: 'panel', id: 'fam-trees'});
  t.append(el('h2', {}, 'Trees'));
  const trees = Object.entries(D.tree_meta || {}).filter(([k, x]) => x.family === fam);
  if (trees.length) trees.forEach(([k, x]) => t.append(el('div', {style: 'padding:5px 0;border-bottom:1px solid var(--line)'},
    link(treeHref(k), treeLabel(k)),
    el('span', {class: 'hint'}, (x.n_species ? ` · ${x.n_species.toLocaleString()} species` : '') +
      (x.of_family && x.of_family !== x.n_tips ? ` · ${x.n_tips.toLocaleString()} of ${x.of_family.toLocaleString()} sequences` : '')))));
  if (trees.length) host.append(t);
  const dlp = el('div', {class: 'panel', id: 'fam-downloads'});
  dlp.append(el('h2', {}, 'Downloads'));
  const fz = (SITE.family_fasta || {})[fam];
  const row = (a, txt) => el('div', {style: 'padding:4px 0'}, a, el('span', {class: 'hint'}, ' ' + txt));
  if (fz) dlp.append(row(el('a', {href: `data/fasta/${fam}.faa.gz`, download: `${fam}.faa.gz`, class: 'lk'}, `${fam}.faa.gz`),
    `every ${fam} protein in the run, gzipped FASTA, ${kb(fz)}`));
  dlp.append(row(el('a', {href: '#', class: 'lk', onclick: e => { e.preventDefault();
    dl(`${fam}_counts.tsv`, ['species\traw' + (D.per10k ? '\tper10k' : ''),
      ...SP.map((s, i) => [s, D.counts[fam][i], ...(D.per10k ? [D.per10k[fam][i]] : [])].join('\t'))].join('\n')); }},
    `${fam}_counts.tsv`), 'copies per species in view'));
  host.append(dlp);
}

function exampleWhere() {
  const r = ['order', 'class', ...TAX_COLS].find(c => TAX_COLS.includes(c));
  if (!r) return 'order=Lepidoptera';
  const n = {}; D.species.forEach(s => { const v = (s.tax || {})[r]; if (v) n[v] = (n[v] || 0) + 1; });
  const top = Object.entries(n).sort((a, b) => b[1] - a[1])[0];
  return `${r}=${top ? top[0] : 'VALUE'}`;
}

/* ---- global search -------------------------------------------------- */
const SEARCH = (() => {
  const items = [];
  D0.species.forEach(s => { const tax = s.tax || {};
    items.push({type: 'Species', label: sciName(s.name), sub: [commonName(s.name), tax.order || tax.class || ''].filter(x => x).join(' · '),
      text: [sciName(s.name), commonName(s.name), s.name].join(' ').toLowerCase(), href: spHref(s.name), ital: true}); });
  D.families.forEach(f => items.push({type: 'Family', label: f.family, sub: FAMDESC[f.family] || f.category,
    text: [f.family, f.pfam, f.category, FAMDESC[f.family] || ''].join(' ').toLowerCase(), href: famHref(f.family)}));
  RANK_COLS.forEach(r => { const n = {};
    D0.species.forEach(s => { const v = (s.tax || {})[r]; if (v) n[v] = (n[v] || 0) + 1; });
    Object.entries(n).forEach(([v, k]) => items.push({type: 'Taxon', label: nice(v), sub: `${r} · ${k} species`,
      text: String(v).toLowerCase().replace(/_/g, ' '), act: () => showTaxon(r, v)})); });
  return items;
})();
function runSearch(q) {
  const box = $('#qres'); box.innerHTML = '';
  q = q.trim().toLowerCase();
  if (q.length < 2) { box.classList.remove('on'); return; }
  const score = it => it.text.startsWith(q) || it.label.toLowerCase().startsWith(q) ? 0
    : (' ' + it.text).includes(' ' + q) ? 1 : it.text.includes(q) ? 2 : 9;
  const hits = SEARCH.map(it => [score(it), it]).filter(x => x[0] < 9)
    .sort((a, b) => a[0] - b[0] || a[1].label.localeCompare(b[1].label));
  const groups = {};
  hits.forEach(([, it]) => { (groups[it.type] = groups[it.type] || []).length < 7 && groups[it.type].push(it); });
  ['Species', 'Family', 'Taxon'].forEach(t => { if (!groups[t]) return;
    box.append(el('div', {class: 'grp'}, t === 'Family' ? 'Gene families' : t === 'Taxon' ? 'Taxa' : 'Species'));
    groups[t].forEach(it => { const a = el('a', {href: it.href || '#', onclick: e => { closeSearch();
        if (it.act) { e.preventDefault(); it.act(); } }},
      it.ital ? el('i', {}, it.label) : el('span', {}, it.label), el('small', {}, it.sub));
      box.append(a); }); });
  if (SITE.genes && /^[A-Za-z0-9_.:|-]{3,}$/.test(q)) {
    box.append(el('div', {class: 'grp'}, 'Gene IDs'));
    const mb = Object.values(SITE.family_bytes || {}).reduce((a, b) => a + b, 0);
    box.append(el('a', {href: '#', onclick: e => { e.preventDefault(); geneSearch(q); }},
      el('span', {}, `Search protein IDs for “${q}”`), el('small', {}, `loads ${kb(mb)} once`)));
  }
  if (!box.children.length) box.append(el('div', {class: 'empty'}, 'No species, family or taxon matches.'));
  box.classList.add('on');
}
function geneSearch(q) {
  const box = $('#qres'); box.innerHTML = '';
  const st = el('div', {class: 'empty'}, 'Loading family member lists…'); box.append(st);
  const fams = FAM.filter(f => (SITE.family_bytes || {})[f]);
  let done = 0;
  Promise.all(fams.map(f => DEF.load('family', f).then(d => { done++;
    st.textContent = `Searching… ${done} of ${fams.length} families`; return [f, d]; }).catch(() => [f, null])))
    .then(all => { const hits = [];
      all.forEach(([f, d]) => d && d.members.forEach(r => { if (hits.length < 60 && String(r[1]).toLowerCase().includes(q))
        hits.push([f, r]); }));
      box.innerHTML = '';
      box.append(el('div', {class: 'grp'}, `Protein IDs containing “${q}”${hits.length >= 60 ? ' (first 60)' : ''}`));
      hits.forEach(([f, r]) => { const s = D0.species[r[0]].name;
        box.append(el('a', {href: geneHref(s, r[1]), onclick: closeSearch},
          el('span', {class: 'mono'}, r[1]), el('small', {}, `${f} · ${sciName(s)}`))); });
      if (!hits.length) box.append(el('div', {class: 'empty'}, 'No protein ID matches.'));
    });
}
function closeSearch() { const b = $('#qres'); if (b) { b.classList.remove('on'); } }
function wireHelp() {
  const b = $('#help'); if (!b || !b.addEventListener) return;
  let on = false; try { on = localStorage.getItem('defensome.help') === '1'; } catch (e) {}
  const set = v => { on = v; document.body.classList.toggle('explain', on); b.classList.toggle('on', on);
    try { localStorage.setItem('defensome.help', on ? '1' : '0'); } catch (e) {} };
  b.addEventListener('click', () => set(!on)); set(on);
}
function wireSearch() {
  const q = $('#q'); if (!q || !q.addEventListener) return;
  let t = null, hi = -1;
  q.addEventListener('input', () => { clearTimeout(t); hi = -1; t = setTimeout(() => runSearch(q.value), 90); });
  q.addEventListener('keydown', e => {
    const as = [...$('#qres').querySelectorAll('a')];
    if (e.key === 'ArrowDown' || e.key === 'ArrowUp') { e.preventDefault();
      hi = Math.max(0, Math.min(as.length - 1, hi + (e.key === 'ArrowDown' ? 1 : -1)));
      as.forEach((a, i) => a.classList.toggle('hl', i === hi)); }
    else if (e.key === 'Enter' && as.length) { e.preventDefault(); as[Math.max(0, hi)].click(); q.blur(); }
    else if (e.key === 'Escape') { closeSearch(); q.blur(); } });
  document.addEventListener('click', e => { if (!e.target.closest || !e.target.closest('.search')) closeSearch(); });
  document.addEventListener('keydown', e => { if (e.key === '/' && document.activeElement !== q &&
    !/INPUT|SELECT|TEXTAREA/.test((document.activeElement || {}).tagName || '')) { e.preventDefault(); q.focus(); } });
}

/* ---- sections and routes ---------------------------------------------
   Seven top-level destinations; the analysis views sit under Explore with
   their own sub-navigation. A view with nothing to show is left out of the
   navigation rather than shown as a page explaining its own emptiness. */
const MULTI_METHOD = (() => { const C = D.cmp || {};
  return new Set([...(C.long || []), ...(C.collapse || [])].map(r => r.method).filter(m => m)).size > 1; })();
const EXPLORE = [
  ['heatmap',     'Heatmap',      tabHeatmap,     () => FAM.length > 0],
  ['ordination',  'Ordination',   tabOrdination,  () => VARFAM.length >= 3],
  ['correlation', 'Correlations', tabCorrelation, () => VARFAM.length >= 2],
  ['clans',       'CYP clans',    tabClans,       () => Object.keys(D.clans || {}).length > 0],
  ['amethods',    'Annotation methods', tabAMethods, () => MULTI_METHOD],
];
const TABS = [
  ['overview', 'Home',     tabOverview,    () => true],
  ['species',  'Species',  tabSpeciesList, () => D0.species.length > 0],
  ['families', 'Families', tabFamilies,    () => FAM.length > 0],
  ['family',   'Family',   tabFamily,      () => FAM.length > 0],
  ['compare',  'Compare',  tabCompare,     () => D0.species.length > 1],
  ...EXPLORE,
  ['trees',    'Trees',    tabTrees,       () => Object.keys(D.trees || {}).length > 0],
  ['about',    'About',    tabAbout,       () => true],
];
const avail = id => { const t = TABS.find(x => x[0] === id); try { return !!(t && t[3]()); } catch (e) { return false; } };
const SHOWN = TABS.filter(t => avail(t[0]));
const NAV = [['home', 'Home', '#/'], ['species', 'Species', '#/species'], ['families', 'Families', '#/families'],
  ['compare', 'Compare', '#/compare'], ['trees', 'Trees', '#/trees'], ['explore', 'Explore', '#/explore'],
  ['about', 'About', '#/about']].filter(([k]) =>
  k === 'explore' ? EXPLORE.some(e => avail(e[0])) : ['trees', 'compare'].includes(k) ? avail(k) : true);
const LEGACY = {overview: '/', species: '/species', families: '/families', family: '/family',
  compare: '/compare', trees: '/trees', about: '/about', data: '/about', methods: '/about',
  map: '/about', stats: '/about', domains: '/about'};
EXPLORE.forEach(([k]) => { LEGACY[k] = '/explore/' + k; });
if (!MULTI_METHOD) LEGACY.amethods = '/about/collapse';
LEGACY.map = '/about/rules'; LEGACY.stats = '/about/stats'; LEGACY.domains = '/about/domains'; LEGACY.methods = '/about/methods';

function parseRoute(h) {
  h = String(h || '').replace(/^#/, '');
  let [path, qs] = h.split('?');
  if (!path.startsWith('/')) path = LEGACY[path] || '/';
  const q = {}; (qs || '').split('&').filter(x => x).forEach(kv => { const [k, v] = kv.split('=');
    q[decodeURIComponent(k)] = decodeURIComponent(v || ''); });
  const seg = path.split('/').filter(x => x).map(x => { try { return decodeURIComponent(x); } catch (e) { return x; } });
  const firstExplore = (EXPLORE.find(e => avail(e[0])) || ['heatmap'])[0];
  switch (seg[0]) {
    case 'species': return seg[1] ? {sec: 'sp', top: 'species', key: seg[1]} : {sec: 'species', top: 'species'};
    case 'gene': return {sec: 'gene', top: 'species', key: seg[1], pid: seg.slice(2).join('/')};
    case 'families': return {sec: 'families', top: 'families'};
    case 'family': return {sec: 'family', top: 'families', fam: seg[1] || FAM_CUR};
    case 'explore': return {sec: avail(seg[1]) ? seg[1] : firstExplore, top: 'explore'};
    case 'trees': return avail('trees') ? {sec: 'trees', top: 'trees', tree: seg[1], q: q.q} : {sec: 'overview', top: 'home'};
    case 'compare': return avail('compare') ? {sec: 'compare', top: 'compare'} : {sec: 'overview', top: 'home'};
    case 'about': case 'data': case 'methods': return {sec: 'about', top: 'about', part: seg[1]};
    default: return {sec: 'overview', top: 'home'};
  }
}
const NO_FILTER = new Set(['sp', 'gene', 'about']);
let ROUTE = null, SKIP_HASH = null;
const built = {};
function render(r) {
  ROUTE = r;
  if (r.sec === 'family' && r.fam && FAM.includes(r.fam)) FAM_CUR = r.fam;
  const sig = [r.sec, r.key, r.pid, r.sec === 'family' ? FAM_CUR : '', r.part || ''].join('|');
  document.querySelectorAll && document.querySelectorAll('main > section').forEach(n =>
    n.classList.toggle('on', n.id === 't-' + r.sec));
  [...$('#nav').children].forEach(b => b.classList.toggle('on', b.dataset.k === r.top));
  const sub = $('#subnav');
  if (sub) { sub.innerHTML = ''; sub.classList.toggle('on', r.top === 'explore');
    if (r.top === 'explore') EXPLORE.filter(e => avail(e[0])).forEach(([k, label]) =>
      sub.append(el('button', {'data-k': k, class: k === r.sec ? 'on' : '', onclick: () => go(k)}, label))); }
  const fb = $('#filters'); if (fb && fb.classList) fb.classList.toggle('off', NO_FILTER.has(r.sec));
  if (built[r.sec] !== sig) {
    const host = $('#t-' + r.sec); host.innerHTML = '';
    try {
      if (r.sec === 'sp') speciesPage(host, r.key);
      else if (r.sec === 'gene') genePage(host, r.key, r.pid);
      else SHOWN.find(t => t[0] === r.sec)[2](host);
    } catch (e) { host.append(el('div', {class: 'panel'}, el('h2', {}, 'Error'),
      el('pre', {class: 'mono'}, String(e && e.stack || e)))); }
    built[r.sec] = sig;
    if (typeof window.scrollTo === 'function' && (r.sec === 'sp' || r.sec === 'gene' || r.sec === 'family'))
      window.scrollTo(0, 0);
  }
  if (r.sec === 'trees' && tabTrees._open) tabTrees._open(r.tree, r.q);
}
/* Navigate. Accepts a route ('#/species/x') or a legacy section id ('heatmap'). */
function go(id) {
  const h = String(id).startsWith('#') ? id : '#' + (LEGACY[id] || '/');
  const r = parseRoute(h);
  if (!['sp', 'gene'].includes(r.sec) && !SHOWN.some(t => t[0] === r.sec)) return go('#/');
  render(r);
  if (location.hash !== h) { SKIP_HASH = h; location.hash = h; }
}
if (window.addEventListener) window.addEventListener('hashchange', () => {
  if (location.hash === SKIP_HASH) { SKIP_HASH = null; return; }
  SKIP_HASH = null; render(parseRoute(location.hash)); });
NAV.forEach(([k, label, href]) => $('#nav').append(
  el('button', {'data-k': k, onclick: () => go(href)}, label)));

/* ---- taxonomic filter bar -------------------------------------------- */
/* One select per rank, each listing only values still available under the
   other active filters, so the ranks cascade (pick Insecta and the order list
   narrows to insect orders). Changing anything rebuilds the view and re-renders
   the open tab; other tabs rebuild when next opened. */
function valuesFor(rank) {
  const saved = FILTER[rank]; delete FILTER[rank];
  const vals = new Set();
  D0.species.forEach(s => { if (speciesPasses(s)) {
    const v = (s.tax && s.tax[rank]) ?? (rank === D.trait_name ? s.trait : null);
    if (v != null && v !== '') vals.add(String(v)); } });
  if (saved) FILTER[rank] = saved;
  return [...vals].sort();
}
function rerender(redraw = true) {
  const n = rebuildView();
  VARFAM = FAM.filter(f => !INVARIANT[f]);
  for (const k in built) delete built[k];
  if (typeof TREE_CACHE === 'object') for (const k in TREE_CACHE) delete TREE_CACHE[k];
  TABS.forEach(([k]) => { const n_ = $('#t-' + k); if (n_) n_.innerHTML = ''; });
  const fc = $('#fcount');
  if (fc) fc.textContent = n === D0.species.length
    ? `all ${n.toLocaleString()} species` : `${n.toLocaleString()} of ${D0.species.length.toLocaleString()} species`;
  if (redraw) render(parseRoute(location.hash));
  setHeader();
  return n;
}
function buildFilterBar() {
  const bar = $('#filters'); if (!bar) return;
  if (!TAX_COLS.length && !D.trait_name) { bar.style.display = 'none'; return; }
  bar.innerHTML = '';
  const ranks = TAX_COLS.length ? TAX_COLS : [D.trait_name];
  bar.append(el('span', {class: 'flab'}, 'Filter'));
  ranks.forEach(r => {
    const vals = valuesFor(r);
    if (!vals.length) return;
    const sel = el('select', {title: r, onchange: e => {
      if (e.target.value) FILTER[r] = e.target.value; else delete FILTER[r];
      buildFilterBar(); rerender(); }},
      el('option', {value: ''}, `${r} (all)`),
      ...vals.map(v => el('option', {value: v}, `${nice(v)}`)));
    sel.value = FILTER[r] || '';
    if (FILTER[r]) sel.classList.add('active');
    bar.append(sel);
  });
  if (ranks.length > 1 || TAX_COLS.length) {
    bar.append(el('span', {class: 'flab', style: 'margin-left:10px'}, 'Group by'));
    const g = el('select', {onchange: e => { GROUP_BY = e.target.value || null;
      buildFilterBar(); rerender(); }},
      ...(D.trait_name && !TAX_COLS.includes(D.trait_name)
          ? [el('option', {value: D.trait_name}, D.trait_name)] : []),
      ...TAX_COLS.map(r => el('option', {value: r}, r)));
    g.value = GROUP_BY || '';
    bar.append(g);
  }
  bar.append(el('span', {id: 'fcount', class: 'flab'}, ''),
    Object.keys(FILTER).length
      ? el('button', {class: 'act', onclick: () => { FILTER = {}; buildFilterBar(); rerender(); }}, 'clear')
      : el('span', {}));
  const n = D0.species.filter(speciesPasses).length;
  $('#fcount').textContent = n === D0.species.length
    ? `all ${n.toLocaleString()} species` : `${n.toLocaleString()} of ${D0.species.length.toLocaleString()} species`;
}
buildFilterBar();
function setHeader() { const h = $('#hdr'); if (h) h.textContent = ''; }
setHeader();
wireSearch();
wireHelp();
go(location.hash || '#/');
