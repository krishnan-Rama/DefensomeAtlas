const { chromium } = require(process.env.PLAYWRIGHT_MODULE || 'playwright');
let bad = 0; const chk = (n, c) => { console.log((c ? '  OK   ' : '  FAIL ') + n); if (!c) bad++; };
(async () => {
  const b = await chromium.launch();
  const p = await b.newPage({ viewport: { width: 1400, height: 1200 } });
  const errs = []; p.on('pageerror', e => errs.push(e.message));
  await p.goto('file://' + require('path').resolve(process.argv[2]) + '#trees'); await p.waitForTimeout(500);

  // screen position of tip i, from the browser's own transform matrix
  const tipScreen = i => p.evaluate(i => {
    const m = tabTrees._model(), g = tabTrees._geom(), svg = tabTrees._svg(), vp = svg.querySelector('g');
    const t = m.lv[i], pt = svg.createSVGPoint();
    if (g.layout === 'radial') { pt.x = (g.RT + 6) * Math.cos(t.a); pt.y = (g.RT + 6) * Math.sin(t.a); }
    else { pt.x = g.ringsOuter - 6; pt.y = t.y; }
    const s = pt.matrixTransform(vp.getScreenCTM()), r = svg.getBoundingClientRect();
    return { x: s.x, y: s.y, prot: t.info.prot, inside: s.x > Math.max(r.left, 0) && s.x < Math.min(r.right, innerWidth) && s.y > Math.max(r.top, 0) && s.y < Math.min(r.bottom, innerHeight) };
  }, i);
  const tipUnderMouse = async (i) => {
    const s = await tipScreen(i); if (!s.inside) return null;
    await p.mouse.move(s.x, s.y); await p.waitForTimeout(40);
    const txt = await p.$eval('#tip', e => e.innerText);
    return txt.includes(s.prot);
  };
  const probe = async (label) => {
    const n = await p.evaluate(() => tabTrees._model().lv.length);
    let tested = 0, ok = 0;
    for (const i of [2, 17, 40, 63, 88, 101, 119]) { const r = await tipUnderMouse(i);
      if (r === null) continue; tested++; if (r) ok++; }
    chk(`${label}: tooltip names the tip under the real mouse (${ok}/${tested})`, tested > 0 && ok === tested);
  };
  await p.$eval('#t-trees svg', e => e.scrollIntoView({block: 'center'}));
  await probe('radial, unzoomed');

  // labels must sit inside the canvas at the starting zoom
  const clipped = await p.evaluate(() => { const svg = tabTrees._svg(), r = svg.getBoundingClientRect();
    return [...svg.querySelectorAll('text')].filter(t => { const b = t.getBoundingClientRect();
      return b.width && (b.left < r.left - 1 || b.right > r.right + 1 || b.top < r.top - 1 || b.bottom > r.bottom + 1); }).length; });
  const nlab = await p.evaluate(() => tabTrees._svg().querySelectorAll('text').length);
  chk(`radial: no tip label clipped at the canvas edge (${clipped} of ${nlab})`, clipped === 0 && nlab >= 120);

  // real wheel zoom around a point, then re-probe
  const bb = await (await p.$('#t-trees svg')).boundingBox();
  await p.mouse.move(bb.x + bb.width * .62, bb.y + bb.height * .40);
  for (let k = 0; k < 4; k++) { await p.mouse.wheel(0, -200); await p.waitForTimeout(60); }
  await p.waitForTimeout(200);
  const k1 = await p.evaluate(() => tabTrees._state.k);
  chk(`wheel zoom raises the scale (k = ${k1.toFixed(2)})`, k1 > 1.5);
  await probe('radial, after wheel zoom');

  // real drag pan, then re-probe
  const cx = bb.x + bb.width / 2, cy = bb.y + bb.height / 2;
  const t0 = await p.evaluate(() => ({ ...tabTrees._state }));
  await p.mouse.move(cx, cy); await p.mouse.down(); await p.mouse.move(cx - 140, cy + 90, { steps: 8 }); await p.mouse.up();
  await p.waitForTimeout(150);
  const t1 = await p.evaluate(() => ({ ...tabTrees._state }));
  chk('drag pans the view', Math.abs(t1.tx - t0.tx) > 20 && Math.abs(t1.ty - t0.ty) > 20);
  await probe('radial, after drag pan');

  // clicking a branch point selects that clade
  await p.dblclick('#t-trees svg'); await p.waitForTimeout(200);
  const node = await p.evaluate(() => { const m = tabTrees._model(), svg = tabTrees._svg(), vp = svg.querySelector('g');
    const n = m.all.find(x => x.support === 0.61), pt = svg.createSVGPoint(); pt.x = n.x; pt.y = n.y;
    const s = pt.matrixTransform(vp.getScreenCTM()); return { x: s.x, y: s.y, n: n.hi - n.lo + 1 }; });
  await p.mouse.click(node.x, node.y); await p.waitForTimeout(150);
  const side = await p.$eval('#t-trees', e => e.innerText);
  chk(`clicking a branch point opens that clade (${node.n} tips)`, side.includes(`Clade: ${node.n} tips`));
  chk('a single-method clade is flagged in the side panel', side.includes('Every tip here comes from denovo'));

  // zoom-to-clade must frame every tip of the clade inside the canvas
  await p.click('text=zoom to clade'); await p.waitForTimeout(300);
  const framed = await p.evaluate(() => { const m = tabTrees._model(), svg = tabTrees._svg(), vp = svg.querySelector('g');
    const S = tabTrees._state, n = S.sel, r = svg.getBoundingClientRect(), g = tabTrees._geom();
    return m.lv.slice(n.lo, n.hi + 1).every(t => { const pt = svg.createSVGPoint();
      pt.x = (g.RT + 6) * Math.cos(t.a); pt.y = (g.RT + 6) * Math.sin(t.a);
      const s = pt.matrixTransform(vp.getScreenCTM());
      return s.x > r.left && s.x < r.right && s.y > r.top && s.y < r.bottom; }); });
  chk('zoom to clade keeps every clade tip on screen', framed);

  // rectangular layout
  await p.selectOption('#t-trees select >> nth=1', 'rect'); await p.waitForTimeout(300);
  await p.$eval('#t-trees svg', e => e.scrollIntoView({block: 'start'}));
  const stale = await p.$eval('#tip', e => getComputedStyle(e).opacity);
  chk('switching layout clears any stale tooltip', stale === '0');
  await probe('rectangular, unzoomed');
  for (let k = 0; k < 3; k++) { await p.mouse.wheel(0, -200); await p.waitForTimeout(60); }
  await p.waitForTimeout(200);
  await probe('rectangular, after wheel zoom');

  // search
  await p.fill('#t-trees input[type=search]', 'frag1'); await p.waitForTimeout(400);
  const found = await p.$eval('#t-trees', e => e.innerText);
  chk('search reports its matches', /3 matches for/.test(found));

  chk(`no page errors (${errs.length})`, errs.length === 0); if (errs.length) console.log(errs);
  await b.close();
  console.log(bad ? `\n${bad} BROWSER FAILURES` : '\nbrowser tests passed');
  process.exitCode = bad ? 1 : 0;
})().catch(e => { console.error('FAILED', e); process.exit(1); });
