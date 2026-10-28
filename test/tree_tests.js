
/* ---- tree tab: known-answer tests (run with the method-comparison fixture) */
function approxT(a,b,t=1e-6){ return Math.abs(a-b)<t; }
let _tb=0; const tchk=(n,c)=>{ console.log((c?'  OK   ':'  FAIL ')+n); if(!c) _tb++; };
if (D.trees && D.trees.CYP && D.samples && D.samples.X_alpha_pep_DToL) {   // known-answer fixture only
  const m = treeModel('CYP');
  tchk('tips parsed: 12 genes x 9 + 12 injected = 120', m.lv.length === 120);
  tchk('every tip annotated with a method from the sample sheet',
       m.lv.every(t => ['DToL','genome_guided','denovo'].includes(t.info.method)));
  tchk('every tip carries its CYP clan', m.lv.every(t => ['CYP2','CYP3','CYP4','MITO'].includes(t.info.clan)));
  tchk('fragment status carried only on injected fragments',
       m.lv.filter(t => t.info.status === 'FRAGMENT').length === 12 &&
       m.lv.filter(t => t.info.status === 'FRAGMENT').every(t => t.info.prot.includes('frag')));
  tchk('support values parsed from internal labels', m.hasSupport &&
       m.all.filter(n => n.support === 0.95).length === 12 && m.all.filter(n => n.support === 0.61).length === 3);
  tchk('clade tips are a contiguous DFS range', m.all.every(n => n.hi - n.lo + 1 === n.nt));
  const ex = exclusiveClades(m, 'method', 3);
  tchk(`single-source clades by method: exactly the 3 injected (${ex.length} found)`, ex.length === 3);
  tchk('  all three are denovo', ex.every(e => e.value === 'denovo'));
  tchk('  sizes are 5, 4, 3', ex.map(e => e.n).join(',') === '5,4,3');
  tchk('  none of the 12 real gene clades is reported', !ex.some(e => e.n === 9));
  const byClan = exclusiveClades(m, 'clan', 3);
  tchk('by clan, each 9-tip gene clade is single-clan', byClan.filter(e => e.n >= 9).length >= 12);
  const gene = m.all.find(n => n.support === 0.95);
  const st = cladeStats(m, gene);
  tchk('clade stats: a gene clade has 3 tips per method',
       st.n === 9 && st.method.DToL === 3 && st.method.genome_guided === 3 && st.method.denovo === 3);
  tchk('clade stats: 3 species x 3', Object.values(st.species).every(v => v === 3) && Object.keys(st.species).length === 3);
  // geometry and hit testing
  const h = document.createElement('div'); tabTrees(h);
  const g = tabTrees._geom(), mm = tabTrees._model();
  const t7 = mm.lv[37];
  const hit = tabTrees._tipAt({x: (g.RT + 4) * Math.cos(t7.a), y: (g.RT + 4) * Math.sin(t7.a)});
  tchk('radial hit test returns the tip under the pointer', hit === t7);
  const inner = mm.all.find(n => n.kids.length && n.support === 0.61);
  tchk('node hit test finds a clicked branch point', tabTrees._nodeAt({x: inner.x + .5, y: inner.y}) === inner);
  const S = tabTrees._state, k0 = S.k;
  tabTrees._zoomBy(2);
  tchk('zoom doubles the scale', approxT(S.k, 2 * k0));
  tabTrees._fitClade(inner);
  tchk('zoom-to-clade magnifies a small clade', S.k > 2);
  tchk('tip lookup outside the tree returns nothing', tabTrees._tipAt({x: 0, y: 0}) === null);
}
console.log(_tb ? '\n' + _tb + ' TREE FAILURES' :
  (D.samples && D.samples.X_alpha_pep_DToL) ? '\ntree tests passed' : '\n(tree known-answer tests skipped: not the fixture)');
if (_tb) process.exitCode = 1;
