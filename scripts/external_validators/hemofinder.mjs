// HemoFinder (dbAMP/CUHK) driver — async result page, poll until table populates.
import { chromium } from '@playwright/test';
import fs from 'fs';

const FASTA = process.env.FASTA_PATH;
const OUT = process.env.OUT_CSV;
const SHOT = process.env.SHOT_PATH || 'screenshot.png';
const CHUNK = 100;

function readFasta(p){const t=fs.readFileSync(p,'utf8').trim().split('\n');const o=[];for(let i=0;i<t.length;i+=2)o.push([t[i].slice(1),t[i+1]]);return o;}
const seqs = readFasta(FASTA);
const chunks=[];for(let i=0;i<seqs.length;i+=CHUNK)chunks.push(seqs.slice(i,i+CHUNK));

const b = await chromium.launch({ headless: true });
const rows = [];
for (let ci=0; ci<chunks.length; ci++){
  const fasta = chunks[ci].map(([id,s])=>`>${id}\n${s}`).join('\n');
  const tmp = `/tmp/pw-validate/hemo_chunk_${ci}.fasta`; fs.writeFileSync(tmp, fasta+'\n');
  const p = await b.newPage();
  await p.goto('https://ycclab.cuhk.edu.cn/dbAMP/HemoFinder.php', { timeout:120000, waitUntil:'domcontentloaded' });
  await p.setInputFiles('input[name="SEQFILE"]', tmp).catch(async()=>{ await p.fill('textarea[name="SEQTEXT"]', fasta); });
  await p.click('#startfinder');
  let got=null;
  for (let i=0;i<30;i++){
    await p.waitForTimeout(5000);
    got = await p.evaluate(()=>{
      if(!/Low-hemolysis|High-hemolysis/i.test(document.body.innerText)) return null;
      return [...document.querySelectorAll('table tr')].map(tr=>[...tr.querySelectorAll('td,th')].map(td=>td.innerText.trim())).filter(r=>r.length>=3);
    });
    if(got) break;
  }
  if (ci===0 && got) await p.screenshot({ path:SHOT, fullPage:true }).catch(()=>{});
  let n=0;
  for (const r of (got||[])){
    if (/^XPRT_/.test(r[0])){ rows.push({candidate_id:r[0], hemo_call:r[2], half_life:r[3]||''}); n++; }
  }
  console.log(`chunk ${ci+1}/${chunks.length}: +${n} (total ${rows.length})`);
  await p.close();
}
await b.close();

const hdr='candidate_id,hemo_call,half_life,is_nonhemolytic\n';
const body=rows.map(r=>`${r.candidate_id},${r.hemo_call},${r.half_life},${/low-?hemolysis/i.test(r.hemo_call)}`).join('\n');
fs.writeFileSync(OUT, hdr+body+'\n');
const low=rows.filter(r=>/low-?hemolysis/i.test(r.hemo_call)).length;
console.log(`WROTE ${rows.length} rows -> ${OUT}`);
console.log(`HemoFinder Low-hemolysis: ${low}/${rows.length}`);
