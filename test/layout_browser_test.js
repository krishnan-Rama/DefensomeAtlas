/* Catalogue-scale layout: the views that only clutter at 30+ groups and
   1,000+ species. Run against make_taxonomy_fixture.py --catalogue. */
const { chromium } = require(process.env.PLAYWRIGHT_MODULE || 'playwright');
let bad = 0;
const chk = (n, c, x) => { console.log((c ? '  OK   ' : '  FAIL ') + n + (!c && x ? '\n         ' + x : '')); if (!c) bad++; };
(async () => {
  const b = await chromium.launch();
  const p = await b.newPage({ viewport: { width: 1440, height: 1200 } });
  const errs = []; p.on('pageerror', e => errs.push(e.message));
  await p.goto('file://' + require('path').resolve(process.argv[2])); await p.waitForTimeout(800);
  const go = async t => { await p.evaluate(t => go(t), t); await p.waitForTimeout(900); };

  await go('family');
  const f = await p.evaluate(() => {
    const sv = document.querySelector('#t-family svg');
    return {rows: sv.querySelectorAll('text.axlab').length,
            rotated: [...sv.querySelectorAll('text.axlab')].filter(t => t.getAttribute('transform')).length,
            dots: sv.querySelectorAll('circle').length,
            fills: new Set([...sv.querySelectorAll('rect')].map(r => r.getAttribute('fill'))).size,
            note: (sv.parentElement.parentElement.innerText.match(/\d+ groups with fewer than 3 species not shown/) || [''])[0],
            height: document.querySelector('#t-family').scrollHeight,
            openTables: [...document.querySelectorAll('#t-family table')]
              .filter(t => t.checkVisibility() && t.querySelectorAll('tbody tr').length > 12).length};
  });
  chk(`32 classes become horizontal rows, none rotated (${f.rows} rows)`, f.rows >= 25 && f.rotated === 0, JSON.stringify(f));
  chk(`no cloud of 1,236 dots over the boxes (${f.dots})`, f.dots === 0);
  chk(`one colour once groups outnumber the palette (${f.fills})`, f.fills === 1);
  chk(`tiny groups are named in a line, not drawn (${f.note})`, /^\d+ groups/.test(f.note));
  chk(`no long table open by default on Families (${f.openTables})`, f.openTables === 0);
  chk(`Families tab fits in ${f.height}px (budget 2,200)`, f.height < 2200);

  chk('no page errors', errs.length === 0, errs.join(' | '));
  await b.close();
  console.log(bad ? `\n${bad} FAILED` : '\n  catalogue layout tests passed');
  process.exitCode = bad ? 1 : 0;
})();
