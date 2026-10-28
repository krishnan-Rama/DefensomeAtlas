/* The linked site: routing, search, lazy loading from file://, and that what
   a gene page shows is what is in the proteome. Run against
   make_taxonomy_fixture.py --genes, exported with a small tree-embed limit so
   the trees load from data/ too. Usage: node site_browser_test.js <site dir> <results dir> */
const { chromium } = require(process.env.PLAYWRIGHT_MODULE || 'playwright');
const fs = require('fs'), path = require('path');
let bad = 0;
const chk = (n, c, x) => { console.log((c ? '  OK   ' : '  FAIL ') + n + (!c && x ? '\n         ' + x : '')); if (!c) bad++; };
const site = path.resolve(process.argv[2]), res = path.resolve(process.argv[3]);

function readFasta(f) { const o = {}; let id = null;
  fs.readFileSync(f, 'utf8').split('\n').forEach(l => { if (l.startsWith('>')) { id = l.slice(1).split(/\s/)[0]; o[id] = ''; }
    else if (id) o[id] += l.trim(); }); return o; }
const tsv = f => { const [h, ...r] = fs.readFileSync(f, 'utf8').trim().split('\n').map(l => l.split('\t'));
  return r.map(x => Object.fromEntries(h.map((k, i) => [k, x[i]]))); };

(async () => {
  const calls = tsv(path.join(res, 'gene_calls_all.tsv'));
  const doms = tsv(path.join(res, 'domains', 'domain_table.tsv'));
  const b = await chromium.launch();
  const ctx = await b.newContext({ acceptDownloads: true, viewport: { width: 1440, height: 1100 } });
  const p = await ctx.newPage();
  const errs = []; p.on('pageerror', e => errs.push(e.message));
  const url = 'file://' + path.join(site, 'index.html');
  const nav = async h => { await p.evaluate(h => { location.hash = h; }, h); await p.waitForTimeout(700); };
  await p.goto(url + '#/'); await p.waitForTimeout(700);

  /* ---------------------------------------------------------- search */
  await p.fill('#q', 'insecta'); await p.waitForTimeout(300);
  const taxa = await p.$$eval('#qres a', as => as.map(a => a.innerText));
  chk(`search offers the taxon (${taxa.length} results)`, taxa.some(t => /^Insecta/.test(t)), taxa.join(' | '));
  await p.click('#qres a:has(span:text-is("Insecta"))'); await p.waitForTimeout(600);
  const ins = await p.evaluate(() => ({hash: location.hash, n: D.species.length,
    rows: document.querySelector('#t-species .count') ? document.querySelector('#t-species .count').innerText : ''}));
  chk(`choosing a taxon opens the species list filtered to it (${ins.n})`, ins.hash === '#/species' && ins.n === 46, JSON.stringify(ins));
  await p.evaluate(() => { FILTER = {}; buildFilterBar(); rerender(); });

  await p.fill('#q', 'PF02798'); await p.waitForTimeout(300);
  await p.keyboard.press('Enter'); await p.waitForTimeout(600);
  chk('a Pfam accession finds its family and Enter opens it',
    await p.evaluate(() => location.hash === '#/family/GST' && /GST/.test(document.querySelector('#t-family h1').innerText)));

  /* ---------------------------------------------------- species page */
  const sid = 'Genusins070_species070_pep';
  await nav('#/species');
  await p.fill('#t-species input[type=search]', 'species070'); await p.waitForTimeout(200);
  await p.click('#t-species a:has-text("Genusins070 species070")'); await p.waitForTimeout(900);
  const sp = await p.evaluate(() => ({hash: location.hash, h1: document.querySelector('#t-sp h1').innerText,
    rows: document.querySelectorAll('#sp-profile svg rect').length,
    count: (document.querySelector('#sp-genes .count') || {}).innerText || '',
    filters: getComputedStyle(document.querySelector('#filters')).display}));
  const nGenes = calls.filter(r => r.species === sid).length;
  chk('clicking a name opens its page', sp.hash === '#/species/' + sid && sp.h1 === 'Genusins070 species070', JSON.stringify(sp));
  chk(`profile has one bar per family (${sp.rows})`, sp.rows === 37);
  chk(`genes load from data/ over file:// and match the calls (${nGenes})`, sp.count === nGenes + ' genes', sp.count);
  chk('the filter bar is hidden on an entity page', sp.filters === 'none');

  const [dlf] = await Promise.all([p.waitForEvent('download'), p.click('#sp-genes button:has-text("FASTA")')]);
  const fa = fs.readFileSync(await dlf.path(), 'utf8');
  const prot = readFasta(path.join(res, 'proteomes', sid + '.faa'));
  const got = readFasta(path.join(await dlf.path()));
  const ids = Object.keys(got).map(k => k.split('|')[1]);
  chk(`species FASTA holds every gene, sequences identical to the proteome (${ids.length})`,
    ids.length === nGenes && ids.every(i => got[sid + '|' + i] === prot[i]));
  chk('decoys are not in it', !/decoy/.test(fa));

  /* breadcrumb: a rank opens the filtered list */
  await p.click('#t-sp .crumb a:has-text("Insecta")'); await p.waitForTimeout(600);
  chk('a rank in the breadcrumb filters the species list', await p.evaluate(() => location.hash === '#/species' && D.species.length === 46));
  await p.goBack(); await p.waitForTimeout(700);
  chk('Back returns to the species page', await p.evaluate(() => location.hash.startsWith('#/species/Genusins070')));
  await p.evaluate(() => { FILTER = {}; buildFilterBar(); rerender(); });

  /* ------------------------------------------------------- gene pages */
  const picks = calls.filter(r => r.species === sid && ['GST', 'CYP', 'ABC_full'].includes(r.family)).slice(0, 40)
    .filter((r, i) => i % 13 === 0);
  let seqOK = 0, domOK = 0;
  for (const r of picks) {
    await nav(`#/gene/${encodeURIComponent(sid)}/${encodeURIComponent(r.protein)}`);
    const g = await p.evaluate(() => ({seq: [...document.querySelectorAll('#t-gene .seq')].map(x => x.innerText)[0] || '',
      boxes: document.querySelector('#t-gene .panel svg') ? document.querySelector('#t-gene .panel svg').querySelectorAll('rect').length : -1,
      h1: (document.querySelector('#t-gene h1') || {}).innerText}));
    const shown = g.seq.split('\n').map(l => l.replace(/^\s*\d+\s+/, '').replace(/\s/g, '')).join('');
    if (shown === prot[r.protein] && g.h1 === r.protein) seqOK++;
    if (g.boxes === doms.filter(d => d.species === sid && d.protein === r.protein).length) domOK++;
  }
  chk(`gene pages show the proteome's sequence exactly (${seqOK}/${picks.length})`, seqOK === picks.length && picks.length >= 3);
  chk(`one domain box per domain-table row (${domOK}/${picks.length})`, domOK === picks.length);

  /* ------------------------------------------- family members, filtered */
  await nav('#/family/CYP');
  const mem = async () => p.evaluate(() => (document.querySelector('#fam-members .count') || {}).innerText || '');
  const all = await mem();
  const cypAll = calls.filter(r => r.family === 'CYP').length;
  chk(`family members load (${all.split(' ')[0]})`, all.startsWith(cypAll.toLocaleString('en-US') + ' genes'), all);
  await p.evaluate(() => { FILTER = {class: 'Insecta'}; buildFilterBar(); rerender(); }); await p.waitForTimeout(700);
  const insN = calls.filter(r => r.family === 'CYP' && /^Genusins/.test(r.species)).length;
  const sub = await mem();
  chk(`members follow the filter bar (${sub.split(' ')[0]} of ${cypAll})`, sub.startsWith(insN.toLocaleString('en-US') + ' genes'), sub);
  await p.evaluate(() => { FILTER = {}; buildFilterBar(); rerender(); });
  const trees = await p.$$eval('#fam-trees a', as => as.map(a => a.getAttribute('href')));
  chk(`the family page lists its trees (${trees.length})`, trees.length === 2 && trees.some(h => /order-insectaformes/.test(h)));

  /* ----------------------------------------------- trees load on demand */
  const pid = calls.find(r => r.species === sid && r.family === 'CYP').protein;
  await nav(`#/trees/CYP__order-insectaformes?q=${encodeURIComponent(pid)}`); await p.waitForTimeout(800);
  const tr = await p.evaluate(() => ({n: tabTrees._model() ? tabTrees._model().lv.length : 0,
    tip: tabTrees._state.tip ? tabTrees._state.tip.info.prot : null,
    embedded: D.tree_meta && D.trees['CYP__order-insectaformes'] !== undefined}));
  const nIns = calls.filter(r => r.family === 'CYP' && /^Genusins/.test(r.species)).length;
  chk(`a large tree loads from data/ when opened (${tr.n} tips)`, tr.n === nIns, JSON.stringify(tr));
  chk('arriving with ?q= selects that gene', tr.tip === pid, tr.tip);
  await p.click('#t-trees a:has-text("gene page")'); await p.waitForTimeout(800);
  chk('the tip record links to the gene page', await p.evaluate(p => location.hash.endsWith(encodeURIComponent(p)), pid));

  /* ------------------------------------------------ protein-ID search */
  await p.fill('#q', 'ins070.gst-g3'); await p.waitForTimeout(300);
  await p.click('#qres a:has-text("Search protein IDs")'); await p.waitForTimeout(1200);
  const hits = await p.$$eval('#qres a', as => as.map(a => a.getAttribute('href')));
  chk(`protein-ID search finds the gene (${hits.length})`, hits.some(h => h.endsWith(encodeURIComponent('ins070.gst-g3.t1'))), hits.join(' '));

  /* --------------------------------- below-threshold genes are counted */
  const resc = tsv(path.join(res, 'rescue_calls.tsv'));
  const pass = resc.filter(r => r.validation === 'PASS'), fail = resc.filter(r => r.validation !== 'PASS');
  const fishSid = pass[0].species;
  await nav('#/species/' + fishSid);
  const fishN = calls.filter(r => r.species === fishSid).length + pass.filter(r => r.species === fishSid).length;
  const fc = await p.evaluate(() => (document.querySelector('#sp-genes .count') || {}).innerText || '');
  chk(`validated below-threshold MATE genes are on the species page (${fc})`, fc === fishN + ' genes', fc);
  await p.selectOption('#sp-genes select', 'MATE'); await p.waitForTimeout(300);
  await p.fill('#sp-genes input[type=search]', 'below threshold'); await p.waitForTimeout(300);
  const badge = await p.evaluate(() => [...document.querySelectorAll('#sp-genes .badge')].filter(b => /below threshold/.test(b.textContent)).length);
  chk(`and marked as below threshold (${badge} of the MATE genes)`, badge === 2);
  const mate = await p.evaluate(s => D.counts.MATE[SPI[s]], fishSid);
  const mateGA = calls.filter(r => r.species === fishSid && r.family === 'MATE').length;
  chk(`dashboard MATE count = GA + validated (${mate} = ${mateGA} + 2)`, mate === mateGA + 2);
  await nav(`#/gene/${encodeURIComponent(fishSid)}/${encodeURIComponent(pass[0].protein)}`);
  const rg = await p.evaluate(() => ({badge: /below threshold/.test(document.querySelector('#t-gene .ehead').innerText),
    boxes: document.querySelector('#t-gene .panel svg').querySelectorAll('rect').length,
    seq: [...document.querySelectorAll('#t-gene .seq')].map(x => x.innerText)[0] || ''}));
  const fprot = readFasta(path.join(res, 'proteomes', fishSid + '.faa'));
  chk('a below-threshold gene page has its badge, both MatE repeats and its sequence',
    rg.badge && rg.boxes === 2 && rg.seq.split('\n').map(l => l.replace(/^\s*\d+\s+/, '').replace(/\s/g, '')).join('') === fprot[pass[0].protein],
    JSON.stringify({badge: rg.badge, boxes: rg.boxes}));
  chk('a gene that failed validation has no page entry',
    await p.evaluate(([s, pid]) => DEF.get('species', SPKEY[s]).genes.every(g => g[0] !== pid), [fail[0].species, fail[0].protein]));

  /* ------------------------------------------------- structure panel */
  const gpid = 'ins070.cyp-g3.t1';
  await nav(`#/gene/${encodeURIComponent(sid)}/${encodeURIComponent(gpid)}`);
  const sp_ = await p.evaluate(() => ({card: (document.querySelector('#gene-structure .card') || {}).innerText || '',
    chips: [...document.querySelectorAll('#gene-structure .chip.on')].map(c => c.textContent)}));
  chk(`a finished prediction shows its confidence and binding (${sp_.card.replace(/\n/g, ' | ')})`,
    /permethrin/.test(sp_.card) && /pTM 0\.77/.test(sp_.card) && /binder 73%/.test(sp_.card) && /3\.2 µM/.test(sp_.card));
  chk(`CYP comes preset with haem (${sp_.chips.join(', ')})`, sp_.chips.join() === 'haem');
  const pyYaml = fs.readFileSync(fs.readdirSync(path.join(res, 'structures')).map(dn =>
    path.join(res, 'structures', dn)).filter(dn => dn.includes('permethrin')).map(dn =>
    path.join(dn, fs.readdirSync(dn).find(f => f.endsWith('.yaml'))))[0], 'utf8');
  const js = await p.evaluate(() => { structurePanel._set('CC1(C)C(C=C(Cl)Cl)C1C(=O)OCc1cccc(Oc2ccccc2)c1', 'permethrin');
    return {yaml: structurePanel._yaml(), af: JSON.parse(structurePanel._af()), cmd: structurePanel._cmd()}; });
  chk('the page writes the same Boltz input as the `structure` command', js.yaml === pyYaml, js.yaml + '\n---\n' + pyYaml);
  chk('AlphaFold Server input: one protein chain and haem as CCD_HEM',
    js.af[0].dialect === 'alphafoldserver' && js.af[0].sequences.length === 2 && js.af[0].sequences[1].ligand.ligand === 'CCD_HEM');
  chk('the command names the gene, cofactor and ligand',
    js.cmd.includes(`--gene '${sid}|${gpid}'`) && js.cmd.includes('--cofactor HEM') && js.cmd.includes('--ligand-name permethrin'));
  await p.route(/3Dmol-min\.js/, r => r.abort());
  await p.click('#gene-structure .card'); await p.waitForTimeout(800);
  chk('offline, the viewer says so and still offers the mmCIF',
    await p.evaluate(() => /internet/.test(document.querySelector('#gene-structure .mol').innerText) &&
      !!document.querySelector('#gene-structure .mol button')));

  /* --------------------------------------------------- plain pages */
  await nav('#/');
  const home = await p.evaluate(() => document.querySelector('#t-overview').innerText);
  chk('Home carries no QC, flag or gene-call wording', !/QC|flagged|gene calls|inflated|collapse/i.test(home));
  const hidden = await p.evaluate(() => [...document.querySelectorAll('main p.hint')].filter(x => x.checkVisibility()).length);
  chk(`no explanatory text is visible by default (${hidden})`, hidden === 0);
  await p.click('#help'); await p.waitForTimeout(200);
  await nav('#/families');
  chk('the ? button brings the explanations back', await p.evaluate(() =>
    [...document.querySelectorAll('#t-families p.hint')].some(x => x.checkVisibility())));
  await p.click('#help');

  chk(`no page errors (${errs.length})`, errs.length === 0, errs.join(' | '));

  /* ------------------------------- index.html alone, without data/ */
  const lone = fs.mkdtempSync(path.join(require('os').tmpdir(), 'lone'));
  fs.copyFileSync(path.join(site, 'index.html'), path.join(lone, 'index.html'));
  const q = await ctx.newPage(); const e2 = []; q.on('pageerror', e => e2.push(e.message));
  await q.goto('file://' + path.join(lone, 'index.html') + '#/species/' + sid); await q.waitForTimeout(900);
  const msg = await q.evaluate(() => document.querySelector('#sp-genes').innerText);
  chk('index.html on its own still opens; gene sections say where the data went', /data\/ folder/.test(msg) && e2.length === 0, msg + e2.join('|'));
  await q.goto('file://' + path.join(lone, 'index.html') + '#/explore/heatmap'); await q.waitForTimeout(700);
  chk('the comparative views need nothing but index.html', await q.evaluate(() => !!tabHeatmap._canvas()));
  fs.rmSync(lone, {recursive: true});

  await b.close();
  console.log(bad ? `\n${bad} FAILED` : '\n  site tests passed');
  process.exitCode = bad ? 1 : 0;
})().catch(e => { console.error('FAILED', e); process.exit(1); });
