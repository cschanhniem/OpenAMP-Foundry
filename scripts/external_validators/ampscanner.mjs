// AMP Scanner v2 driver — file upload, inline results table (ID, Class, Probability).
import { chromium } from '@playwright/test';
import fs from 'fs';

const OUT = process.env.OUT_CSV;
const SHOT = process.env.SHOT_PATH || 'screenshot.png';

const b = await chromium.launch({ headless: true });
const p = await b.newPage();
await p.goto('https://www.dveltri.com/ascan/v2/ascan.html', { timeout: 45000 });
await p.setInputFiles('input[name="seqInputFile"]', process.env.FASTA_PATH);
await Promise.all([
  p.waitForLoadState('networkidle', { timeout: 180000 }).catch(() => {}),
  p.click('input[type=submit], button[type=submit]'),
]);
await p.waitForTimeout(8000);
await p.screenshot({ path: SHOT, fullPage: true }).catch(() => {});

const rows = await p.evaluate(() => {
  const out = [];
  for (const tr of document.querySelectorAll('table tr')) {
    const c = [...tr.querySelectorAll('td,th')].map(td => td.innerText.trim());
    if (c.length >= 3) out.push(c);
  }
  return out;
});
// Each data row: ["XPRT_0001\nSEQ", "Non-AMP", "0.1306"]
const parsed = [];
for (const r of rows) {
  const id = (r[0].match(/XPRT_\d+/) || [])[0];
  if (!id) continue;
  parsed.push({ candidate_id: id, call: r[r.length - 2], prob: r[r.length - 1] });
}
await b.close();

const hdr = 'candidate_id,ampscanner_call,ampscanner_score,is_amp_positive\n';
const body = parsed.map(r =>
  `${r.candidate_id},${r.call},${r.prob},${/^AMP$/i.test(r.call) || parseFloat(r.prob) > 0.5}`
).join('\n');
fs.writeFileSync(OUT, hdr + body + '\n');
const pos = parsed.filter(r => parseFloat(r.prob) > 0.5).length;
console.log(`WROTE ${parsed.length} rows -> ${OUT}`);
console.log(`AMPScanner AMP+ (prob>0.5): ${pos}/${parsed.length}`);
