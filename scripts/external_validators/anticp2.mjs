// AntiCP2 driver — submits candidates in chunks, parses the disp.php result table.
import { chromium } from '@playwright/test';
import fs from 'fs';

const FASTA = process.env.FASTA_PATH;
const OUT = process.env.OUT_CSV;
const SHOT = process.env.SHOT_PATH || 'screenshot.png';
const CHUNK = 150;

function readFasta(path) {
  const txt = fs.readFileSync(path, 'utf8').trim().split('\n');
  const out = [];
  for (let i = 0; i < txt.length; i += 2) out.push([txt[i].slice(1), txt[i + 1]]);
  return out;
}

const seqs = readFasta(FASTA);
const chunks = [];
for (let i = 0; i < seqs.length; i += CHUNK) chunks.push(seqs.slice(i, i + CHUNK));

const b = await chromium.launch({ headless: true });
const rows = [];
for (let ci = 0; ci < chunks.length; ci++) {
  const fasta = chunks[ci].map(([id, s]) => `>${id}\n${s}`).join('\n');
  const p = await b.newPage();
  await p.goto('https://webs.iiitd.edu.in/raghava/anticp2/predict.php', { timeout: 45000 });
  await p.fill('textarea[name="seq"]', fasta);
  await Promise.all([
    p.waitForLoadState('networkidle', { timeout: 90000 }).catch(() => {}),
    p.click('input[type="submit"][value="Submit"]'),
  ]);
  await p.waitForTimeout(3000);
  if (ci === 0) await p.screenshot({ path: SHOT, fullPage: true }).catch(() => {});
  const data = await p.evaluate(() => {
    const trs = [...document.querySelectorAll('table tr')];
    return trs.map(tr => [...tr.querySelectorAll('td,th')].map(td => td.innerText.trim()));
  });
  for (const r of data) {
    // expected: ID, Seq, Score, Prediction, ...
    if (r.length >= 4 && /^XPRT_/.test(r[0])) {
      rows.push({ candidate_id: r[0], anticp_score: r[2], anticp_call: r[3] });
    }
  }
  console.log(`chunk ${ci + 1}/${chunks.length}: parsed ${rows.length} rows so far`);
  await p.close();
}
await b.close();

const header = 'candidate_id,anticp_score,anticp_call,is_non_anticp\n';
const body = rows.map(r =>
  `${r.candidate_id},${r.anticp_score},${r.anticp_call},${/non-?anticp/i.test(r.anticp_call)}`
).join('\n');
fs.writeFileSync(OUT, header + body + '\n');
console.log(`WROTE ${rows.length} rows -> ${OUT}`);
const nonanti = rows.filter(r => /non-?anticp/i.test(r.anticp_call)).length;
console.log(`Non-AntiCP (off-target clear): ${nonanti}/${rows.length}`);
