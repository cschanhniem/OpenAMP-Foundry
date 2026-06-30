// CAMP-R4 driver — 4 classifiers (SVM/RF/ANN/DA), results mapped by input order.
import { chromium } from '@playwright/test';
import fs from 'fs';

const FASTA = '/tmp/pw-validate/candidates.fasta';
const OUT = '/Volumes/SSD/openamp-foundry/outputs/external_validation/camp4_results.csv';
const SHOT = '/Volumes/SSD/openamp-foundry/outputs/external_validation/screenshots/camp4.png';
const CHUNK = 100;

function readFasta(p) {
  const t = fs.readFileSync(p, 'utf8').trim().split('\n');
  const o = []; for (let i = 0; i < t.length; i += 2) o.push([t[i].slice(1), t[i + 1]]); return o;
}
const seqs = readFasta(FASTA);
const chunks = []; for (let i = 0; i < seqs.length; i += CHUNK) chunks.push(seqs.slice(i, i + CHUNK));

const b = await chromium.launch({ headless: true });
const all = [];
for (let ci = 0; ci < chunks.length; ci++) {
  const fasta = chunks[ci].map(([id, s]) => `>${id}\n${s}`).join('\n');
  const p = await b.newPage();
  await p.goto('https://camp3.bicnirrh.res.in/predict/', { timeout: 45000 });
  await p.fill('textarea[name="S1"]', fasta);
  for (const v of ['svm', 'rf', 'ann', 'da']) await p.check(`input[name="algo[]"][value="${v}"]`).catch(() => {});
  await p.click('input[name="B1"]');
  await p.waitForTimeout(Math.min(8000 + chunks[ci].length * 250, 60000));
  if (ci === 0) await p.screenshot({ path: SHOT, fullPage: true }).catch(() => {});
  // Collect compact table rows: [idx, class, prob?]
  const rows = await p.evaluate(() => {
    const out = [];
    for (const tr of document.querySelectorAll('table tr')) {
      const c = [...tr.querySelectorAll('td,th')].map(td => td.innerText.trim());
      if (c.length >= 2 && c.length <= 3) out.push(c);
    }
    return out;
  });
  // Split into 4 sections by header rows "Seq. ID."
  const sections = []; let cur = null;
  for (const r of rows) {
    if (/seq\.?\s*id/i.test(r[0])) { cur = []; sections.push(cur); }
    else if (cur && /^\d+$/.test(r[0])) cur.push(r);
  }
  const [svm = [], rf = [], ann = [], da = []] = sections;
  const get = (sec, idx) => (sec.find(r => +r[0] === idx) || [, '?', '']);
  chunks[ci].forEach(([id], k) => {
    const idx = k + 1;
    const s = get(svm, idx), r = get(rf, idx), a = get(ann, idx), d = get(da, idx);
    all.push({ candidate_id: id, svm: s[1], svm_p: s[2] || '', rf: r[1], rf_p: r[2] || '',
      ann: a[1], da: d[1], da_p: d[2] || '' });
  });
  console.log(`chunk ${ci + 1}/${chunks.length}: ${all.length} rows`);
  await p.close();
}
await b.close();

const ampCall = r => {
  const votes = [r.svm, r.rf, r.ann, r.da].filter(x => /^AMP$/i.test(x)).length;
  return { votes, call: votes >= 2 ? 'AMP' : 'non-AMP' };
};
const hdr = 'candidate_id,camp_svm,camp_svm_prob,camp_rf,camp_rf_prob,camp_ann,camp_da,camp_da_prob,camp_amp_votes,camp_call,is_amp_positive\n';
const body = all.map(r => {
  const { votes, call } = ampCall(r);
  return `${r.candidate_id},${r.svm},${r.svm_p},${r.rf},${r.rf_p},${r.ann},${r.da},${r.da_p},${votes},${call},${call === 'AMP'}`;
}).join('\n');
fs.writeFileSync(OUT, hdr + body + '\n');
const pos = all.filter(r => ampCall(r).call === 'AMP').length;
console.log(`WROTE ${all.length} rows -> ${OUT}`);
console.log(`CAMP AMP+ (>=2/4 classifiers): ${pos}/${all.length}`);
