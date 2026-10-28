/* Real-browser tests for the taxonomic presentation: name display, label
   decluttering, habitat sorting. Run against the 140-species fixture, whose
   known answers are documented in make_taxonomy_fixture.py. */
const { chromium } = require(process.env.PLAYWRIGHT_MODULE || 'playwright');
let bad = 0;
const chk = (n, c, extra) => { console.log((c ? '  OK   ' : '  FAIL ') + n +
  (!c && extra ? '\n         ' + extra : '')); if (!c) bad++; };

(async () => {
  const b = await chromium.launch();
  const p = await b.newPage({ viewport: { width: 1500, height: 1200 } });
  const errs = []; p.on('pageerror', e => errs.push(e.message));
  const url = 'file://' + require('path').resolve(process.argv[2]);
  /* goto() with only a different fragment does not navigate, so the tab would
     never be built; drive the app's own router instead. */
  const goTab = async t => { await p.evaluate(t => go(t), t); await p.waitForTimeout(500); };

  /* ------------------------------------------------- 1. the "pep" suffix */
  await p.goto(url); await p.waitForTimeout(500);
  await goTab('overview');
  await goTab('species');
  const pepText = await p.evaluate(() => {
    const out = [];
    document.querySelectorAll('#t-species td, #t-species th, #t-species option, #t-overview text')
      .forEach(n => { if (/_pep\b/.test(n.textContent)) out.push(n.textContent); });
    return out.slice(0, 3);
  });
  chk('no "_pep" provenance suffix anywhere on screen', pepText.length === 0, pepText.join(' | '));

  /* -------------------------------------- 2. common name beside the name */
  const cn = await p.evaluate(() => {
    const tb = document.querySelector('#t-species table');
    const head = [...tb.querySelectorAll('th')].map(t => t.textContent.trim());
    const j = head.indexOf('common name');
    const rows = [...tb.querySelectorAll('tbody tr')].map(tr =>
      [...tr.querySelectorAll('td')].map(td => td.textContent.trim()));
    return {j, filled: j < 0 ? 0 : rows.filter(r => /number \d+/.test(r[j])).length,
            first: rows[0] ? rows[0].slice(0, 2) : []};
  });
  chk(`species table carries a common-name column (${cn.filled} of the first page filled)`,
    cn.j === 1 && cn.filled >= 40 && /^Genus\w+ species\d+$/.test(cn.first[0]),
    JSON.stringify(cn));

  /* ------------------------------------------------- 3. heatmap renders */
  await goTab('heatmap');
  const cells = await p.evaluate(() => {
    const cv = tabHeatmap._canvas(), g = cv.getContext('2d'), G = tabHeatmap._geom();
    let coloured = 0;
    for (let r = 0; r < 40; r++) for (let j = 0; j < 5; j++) {
      const d = g.getImageData(G.L + j * G.CW + 6, r * G.CH + 4, 1, 1).data;
      if (d[3] > 0 && !(d[0] > 250 && d[1] > 250 && d[2] > 250)) coloured++;
    }
    return coloured;
  });
  chk(`heatmap draws its data cells (${cells}/200 sampled pixels coloured)`, cells === 200);

  /* --------------------------------- 4. row labels thinned, bands named */
  const bands = await p.evaluate(() => tabHeatmap._lastBands);
  const nSp = await p.evaluate(() => D.species.length);
  chk(`140 rows grouped by class get band labels, not 140 names (${bands.length} bands)`,
    nSp === 140 && bands.length === 7, JSON.stringify(bands));
  const bandTot = bands.reduce((a, x) => a + x[1], 0);
  chk(`band sizes account for every row (${bandTot} of ${nSp})`, bandTot === nSp);

  /* ------------------------------------------------ 5. sorting by habitat */
  const hasHab = await p.evaluate(() =>
    [...tabHeatmap._sortSel().options].map(o => o.value).includes('tax:habitat'));
  chk('sort selector offers habitat', hasHab);
  if (hasHab) {
    await p.evaluate(() => tabHeatmap._setSort('tax:habitat'));
    await p.waitForTimeout(300);
    const ord = await p.evaluate(() => tabHeatmap._lastOrder);
    const runs = [];
    ord.forEach(([, h]) => { if (!runs.length || runs[runs.length - 1][0] !== h) runs.push([h, 0]);
      runs[runs.length - 1][1]++; });
    chk(`habitat sort puts each habitat in one contiguous block (${runs.length} runs)`,
      runs.length === 4, JSON.stringify(runs));
    const hb = await p.evaluate(() => tabHeatmap._lastBands.map(x => x[0]));
    chk('bands are relabelled to the habitats',
      JSON.stringify(hb) === JSON.stringify(['estuarine', 'freshwater', 'marine', 'terrestrial']),
      JSON.stringify(hb));
    const cells2 = await p.evaluate(() => {
      const g = tabHeatmap._canvas().getContext('2d'), G = tabHeatmap._geom();
      let n = 0; for (let r = 0; r < 40; r++) {
        const d = g.getImageData(G.L + 6, r * G.CH + 4, 1, 1).data;
        if (d[3] > 0 && !(d[0] > 250 && d[1] > 250 && d[2] > 250)) n++; }
      return n; });
    chk(`cells survive the re-sort (${cells2}/40)`, cells2 === 40);
  }

  /* ------------------------------- 6. ordination labels capped, not soup */
  await goTab('ordination');
  const ord = await p.evaluate(() => {
    const sv = document.querySelectorAll('#t-ordination svg')[0];
    const t = [...sv.querySelectorAll('text')];
    const note = t.find(x => /labelled/.test(x.textContent));
    return {n: t.length, note: note ? note.textContent : ''};
  });
  chk(`ordination labels a readable subset of 140 points (${ord.n} text nodes)`,
    ord.n > 5 && ord.n < 60, JSON.stringify(ord));
  chk('ordination says how many points are labelled', /of 140 labelled/.test(ord.note), ord.note);

  /* ------------------------- 7. family bar chart capped with a toggle */
  await goTab('family');
  const bars0 = await p.evaluate(() => {
    const sv = document.querySelector('#fam-species svg');
    return sv.querySelectorAll('rect').length;
  });
  chk(`family per-species bars show the top 20, not 140 (${bars0})`, bars0 === 20);
  await p.click('#t-family button:has-text("show all 140 species")');
  await p.waitForTimeout(400);
  const bars1 = await p.evaluate(() => {
    const sv = document.querySelector('#fam-species svg');
    return sv.querySelectorAll('rect').length;
  });
  chk(`the toggle opens all 140 (${bars1})`, bars1 === 140);

  /* ------------------------------------------- 8. habitat as a filter */
  await goTab('overview');
  const filtered = await p.evaluate(() => {
    const sel = [...document.querySelectorAll('#filters select')]
      .find(s => s.title === 'habitat');
    if (!sel) return null;
    sel.value = 'marine'; sel.dispatchEvent(new Event('change'));
    return D.species.length;
  });
  chk(`habitat is a filter and marine keeps 35 of 140 (${filtered})`, filtered === 35);

  /* --------------------- 9. empty tabs are not in the nav at all */
  const nav = await p.evaluate(() =>
    [...document.querySelectorAll('#nav button')].map(b => b.textContent));
  chk(`top navigation is the seven destinations less those without data (${nav.join(', ')})`,
    nav.join(',') === 'Home,Species,Families,Compare,Explore,About', nav.join(', '));
  await goTab('heatmap');
  const sub = await p.evaluate(() =>
    [...document.querySelectorAll('#subnav button')].map(b => b.textContent));
  chk(`Explore lists only the analyses with data (${sub.join(', ')})`,
    sub.join(',') === 'Heatmap,Ordination,Correlations', sub.join(', '));
  const empties = await p.evaluate(() => {
    const out = [];
    document.querySelectorAll('main > section').forEach(sec => {
      const t = sec.innerText.replace(/\s+/g, ' ').trim();
      if (/No .* found|No method comparison|No Newick|No group tests/.test(t)) out.push(sec.id);
    });
    return out;
  });
  chk('no tab renders a "nothing here" page', empties.length === 0, empties.join(', '));

  /* --------------------- 10. empty table columns are dropped */
  const cols = await p.evaluate(() => {
    go('species');
    return [...document.querySelectorAll('#t-species table th')].map(t => t.textContent.trim());
  });
  chk('columns with nothing in them are not in the species table',
    !cols.some(c => /QC|worst family|flag|per 10k/.test(c)), cols.join(' | '));
  const dashes = await p.evaluate(() => {
    let n = 0;
    document.querySelectorAll('#t-species table td').forEach(td => {
      const t = td.textContent.trim(); if (t === '\u2013' || t === '\u2014') n++; });
    return n;
  });
  chk(`no column of dashes survives (${dashes} dash cells)`, dashes === 0);

  /* --------------------- 11. the matrix header stays put when rows scroll */
  await p.evaluate(() => { FILTER = {}; buildFilterBar(); rerender(); });  // step 8 left marine on
  await p.waitForTimeout(400);
  await goTab('heatmap');
  const stick = await p.evaluate(async () => {
    const hdr = tabHeatmap._header(), cv = tabHeatmap._canvas();
    const box = cv.parentElement;
    const before = hdr.getBoundingClientRect().top;
    box.scrollTop = 400;
    await new Promise(r => setTimeout(r, 120));
    return {moved: Math.abs(hdr.getBoundingClientRect().top - before),
            scrolled: box.scrollTop, cap: box.clientHeight};
  });
  chk(`column headers hold position while rows scroll (moved ${stick.moved.toFixed(1)}px)`,
    stick.scrolled > 300 && stick.moved < 2, JSON.stringify(stick));
  chk(`the matrix is capped in height, not ${140 * 17}px tall (${stick.cap}px)`,
    stick.cap > 100 && stick.cap <= 640);

  /* ------------- 12. the Overview is a figure, not a table dump */
  await goTab('overview');
  const land = await p.evaluate(() => {
    const h = document.querySelector('#t-overview');
    const big = [...h.querySelectorAll('table')]
      .filter(t => t.querySelectorAll('tbody tr').length > 12).length;
    const segs = h.querySelectorAll('svg rect').length;
    const groups = (tabOverview._rows || []).length;
    return {big, segs, groups};
  });
  chk('the Overview leads with a figure, not a long table',
    land.big === 0 && land.segs > 20, JSON.stringify(land));
  chk(`composition has one row per class (${land.groups})`, land.groups === 7);

  /* ------------- 13. category axis adapts instead of overlapping */
  await goTab('family');
  const axis = await p.evaluate(() => {
    const sv = document.querySelector('#t-family svg');
    const b = sv.getBoundingClientRect();
    const t = [...sv.querySelectorAll('text.axlab')];
    const rot = t.filter(x => /rotate/.test(x.getAttribute('transform') || '')).length;
    let gap = Infinity, overlap = 0;
    if (rot === t.length && t.length > 1) {
      // Rotated labels: their axis-aligned boxes always touch at the corners,
      // so what matters is the distance between baselines measured across them.
      const xs = t.map(x => +x.getAttribute('x')).sort((a, c) => a - c);
      for (let i = 1; i < xs.length; i++)
        gap = Math.min(gap, (xs[i] - xs[i - 1]) * Math.sin(38 * Math.PI / 180));
    } else {
      gap = 99;
      for (let i = 0; i < t.length; i++) for (let j = i + 1; j < t.length; j++) {
        const a = t[i].getBoundingClientRect(), c = t[j].getBoundingClientRect();
        if (a.left < c.right - 1 && a.right - 1 > c.left) overlap++;
      }
    }
    const clipped = t.filter(x => { const r = x.getBoundingClientRect();
      return r.left < b.left - 1 || r.right > b.right + 1; }).length;
    return {n: t.length, rot, gap: Math.round(gap * 10) / 10, overlap, clipped};
  });
  chk(`category axis labels are clear of each other (${axis.n} labels, ` +
      `${axis.rot ? axis.gap + 'px apart' : axis.overlap + ' collisions'})`,
    axis.n >= 5 && axis.overlap === 0 && axis.gap >= 12, JSON.stringify(axis));
  chk('category axis labels are not clipped at the panel edge', axis.clipped === 0,
    JSON.stringify(axis));

  /* ------------- 14. one annotation method: the collapse lives on About, as a distribution */
  const col = await p.evaluate(async () => {
    go('amethods'); await new Promise(r => setTimeout(r, 500));
    const sv = document.querySelector('#about-collapse svg');
    return {hash: location.hash, open: document.querySelector('#about-collapse').open,
            w: sv ? +sv.getAttribute('width') : 0, bars: sv ? sv.querySelectorAll('rect').length : 0};
  });
  chk(`a single-method run charts the collapse distribution on About (${col.w}px wide, ${col.bars} bars)`,
    col.hash === '#/about/collapse' && col.open && col.w > 0 && col.w <= 1000 && col.bars >= 4 && col.bars <= 60, JSON.stringify(col));
  const words = await p.evaluate(async () => { go('#/'); await new Promise(r => setTimeout(r, 400));
    return document.querySelector('#t-overview').innerText; });
  chk('Home says nothing about QC, flags, gene calls or collapse',
    !/QC|flagged|gene calls|collapse|inflated/i.test(words), words.slice(0, 200));

  /* ------------- 15. two-taxa comparison recovers the injected difference */
  const tax = await p.evaluate(async () => {
    go('compare'); await new Promise(r => setTimeout(r, 400));
    const s = taxaComparePanel._sel;
    s.rankSel.value = 'habitat'; s.rankSel.onchange();
    s.aSel.value = 'marine'; s.bSel.value = 'freshwater'; s.aSel.onchange();
    await new Promise(r => setTimeout(r, 300));
    const tb = [...document.querySelectorAll('#t-compare table')][0];
    const head = [...tb.querySelectorAll('th')].map(t => t.textContent.trim());
    const qi = head.findIndex(h => /^q \(BH\)/.test(h));
    const di = head.indexOf('difference');
    const rows = [...tb.querySelectorAll('tbody tr')].map(tr =>
      [...tr.querySelectorAll('td')].map(td => td.textContent.trim()));
    const sig = rows.filter(r => parseFloat(r[qi]) < 0.05);
    return {n: rows.length, sig: sig.map(r => [r[0], r[di], r[qi]])};
  });
  chk(`only the family with an injected habitat difference is significant ` +
      `(${tax.sig.length} of ${tax.n})`,
    tax.sig.length === 1 && tax.sig[0][0] === 'CYP' && tax.sig[0][1] === '20',
    JSON.stringify(tax.sig));

  chk('no page errors', errs.length === 0, errs.join(' | '));
  await b.close();
  console.log(bad ? `\n${bad} FAILED` : '\n  taxonomy UI tests passed');
  process.exitCode = bad ? 1 : 0;
})();
