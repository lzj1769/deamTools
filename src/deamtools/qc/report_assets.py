"""Offline dashboard styles and progressive enhancement (no external resources)."""

CSS = """
:root{--teal:#16767a;--ink:#203438;--muted:#607277;--line:#dae5e6;--bg:#f3f7f7;--good:#27683d;--amber:#865600;--red:#a72d32}
*{box-sizing:border-box}html{scroll-behavior:smooth;scroll-padding-top:86px}body{margin:0;background:var(--bg);color:var(--ink);font:15px/1.55 system-ui,-apple-system,Segoe UI,sans-serif;font-variant-numeric:tabular-nums}main{max-width:1200px;margin:auto;padding:28px 24px 60px}h1{font-size:30px;letter-spacing:-.7px;margin:0;color:var(--teal)}h2{font-size:21px;margin:0 0 16px;color:var(--teal)}h3{font-size:16px;margin:8px 0}p{margin:8px 0 16px}.eyebrow{font-size:12px;font-weight:700;letter-spacing:.1em;text-transform:uppercase;color:var(--muted)}.muted,.intro{color:var(--muted)}.small{font-size:12px}.header-line,.row{display:flex;align-items:center;justify-content:space-between;gap:16px;flex-wrap:wrap}.sample{font-size:20px;margin:8px 0;overflow-wrap:anywhere}.meta-table{margin:14px 0}.meta-table th,.meta-table td{padding:5px 10px}.subnav{display:flex;gap:12px;flex-wrap:wrap;font-size:12px}.subnav a{color:var(--teal)}.meta-table th{width:160px}nav{position:sticky;top:0;z-index:5;background:#fffefc;border-block:1px solid var(--line);box-shadow:0 2px 6px #20343808}nav .nav-inner{max-width:1200px;margin:auto;display:flex;gap:5px;padding:10px 24px;overflow:auto;white-space:nowrap}nav a{padding:8px 12px;border-radius:6px;color:var(--teal);font-size:13px;font-weight:650;text-decoration:none}nav a:hover,nav a:focus{background:#e9f3f3}section,.panel{background:white;border:1px solid var(--line);border-radius:12px;padding:24px;margin:22px 0;scroll-margin-top:90px}.overall{border-left:5px solid var(--teal)}.overall.pass{border-left-color:var(--good)}.overall.fail{border-left-color:var(--red)}.badge{display:inline-block;border-radius:20px;padding:4px 10px;font-size:12px;font-weight:650;background:#eef1f2;color:var(--muted)}.badge.observation{background:#e5f2f3;color:var(--teal)}.badge.good{background:#e7f4e9;color:var(--good)}.badge.advisory{background:#fff2d8;color:var(--amber)}.badge.fail{background:#fce7e8;color:var(--red)}.cards{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:12px;margin:20px 0}.card{background:white;border:1px solid var(--line);border-radius:12px;padding:20px;min-width:0}.card-value{font-size:30px;line-height:1.3;font-weight:650;color:var(--teal);margin:10px 0}.card-label{font-size:13px;color:var(--muted);display:flex;justify-content:space-between;gap:8px}.grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:20px}.bias-grid{display:grid;grid-template-columns:1fr 1.2fr;gap:18px;align-items:start}.grid>*{min-width:0}.flow{display:flex;align-items:center;gap:18px;flex-wrap:wrap;margin:18px 0}.flow .number{font-size:28px;font-weight:650;color:var(--teal)}.flow .arrow{color:#8fa6a7;font-size:24px}.mini-grid{display:grid;grid-template-columns:1fr 1fr;gap:14px}.mini-grid strong{display:block;font-size:21px}.retention{height:12px;background:#e5eeee;border-radius:8px;overflow:hidden;margin:10px 0 18px}.retention span{display:block;height:100%;background:var(--teal)}table{width:100%;border-collapse:collapse;font-size:13px}th{background:#edf4f4;color:#36565a;font-weight:650}th,td{text-align:left;padding:10px 12px;border-bottom:1px solid var(--line);vertical-align:top;overflow-wrap:anywhere}td.value{font-variant-numeric:tabular-nums;white-space:nowrap}.table-scroll{overflow:auto}details{border-top:1px solid var(--line);padding-top:12px;margin-top:16px}summary{cursor:pointer;color:var(--teal);font-weight:650;padding:4px 0 10px}details details{margin-left:8px}.info{position:relative;display:inline-block;cursor:help;border:0;background:none;color:var(--teal);padding:0 4px;font:inherit}.info .tip{display:none;position:absolute;z-index:8;right:0;top:100%;width:min(290px,70vw);background:#173f43;color:white;padding:12px;border-radius:7px;font:12px/1.5 system-ui;text-align:left;box-shadow:0 4px 12px #0002}.info:hover .tip,.info:focus .tip{display:block}.chart{border:1px solid var(--line);border-radius:9px;padding:12px;margin:12px 0;background:white;min-width:0}.chart h3{margin:0}.chart svg{width:100%;height:auto;display:block}.chart-tools{display:flex;align-items:center;gap:8px;flex-wrap:wrap;font-size:12px;margin-top:8px}button,.download{border:1px solid #bcd1d2;border-radius:5px;padding:5px 9px;background:white;color:var(--teal);font:inherit;text-decoration:none;cursor:pointer}button:hover,.download:hover{background:#eef7f7}button:focus-visible,a:focus-visible,summary:focus-visible,.info:focus-visible{outline:3px solid #daae46;outline-offset:3px}.chart-tools input{max-width:100px}.plot-point:focus{outline:none;stroke:#b76d00;stroke-width:3}.method{background:#f6f9f9;border-left:3px solid #a9c8c9;padding:12px 16px;font-size:13px}.downloads{display:flex;gap:7px;flex-wrap:wrap;margin:12px 0}.downloads a{font-size:12px}.reference{border-left:4px solid var(--teal);padding:8px 14px;background:#eef6f6}.recommendations li{margin:10px 0}img{max-width:100%;height:auto}pre{white-space:pre-wrap;overflow-wrap:anywhere;background:#edf4f4;padding:14px;border-radius:6px;font-size:12px}code{font-size:12px;overflow-wrap:anywhere}.empty{color:var(--muted);padding:20px;text-align:center;background:#f7f9f9;border-radius:8px}.footer{color:var(--muted);font-size:12px;margin-top:26px}.no-js-note{font-size:12px}.print-only{display:none}
@media(max-width:850px){.bias-grid{grid-template-columns:1fr 1fr}.cards{grid-template-columns:repeat(2,minmax(0,1fr))}}
@media(max-width:600px){main{padding:18px 12px}section{padding:18px 14px}.grid,.bias-grid{grid-template-columns:1fr}.bias-grid>:last-child{grid-column:auto}.card{padding:14px}.card-value{font-size:25px}nav .nav-inner{padding:8px 12px}.flow{gap:12px}.flow .number{font-size:23px}th,td{padding:8px}.meta-table th{width:110px}}
@media print{body{background:white;font-size:10pt}main{padding:0;max-width:none}nav,button,.chart-tools,.downloads,.no-print{display:none!important}section{break-inside:avoid;box-shadow:none;margin:12px 0}.cards{grid-template-columns:repeat(3,1fr)}details> :not(summary){display:block!important}details{content-visibility:visible}.info .tip{display:block;position:static;background:none;color:var(--muted);width:auto;padding:0;box-shadow:none}.info{font-size:0}.tip{font-size:9pt!important}.print-only{display:block}.meta-table td{word-break:break-all}a{color:inherit;text-decoration:none}}
"""

JS = r"""
(() => {
'use strict';
const reveal = target => {
  let parent = target.parentElement;
  while (parent) { if (parent.tagName === 'DETAILS') parent.open = true; parent = parent.parentElement; }
};
document.querySelectorAll('a[href^="#"]').forEach(link => link.addEventListener('click', () => {
  const target = document.getElementById(link.hash.slice(1)); if (target) reveal(target);
}));
if (location.hash) {const t=document.getElementById(location.hash.slice(1)); if(t){reveal(t);setTimeout(()=>t.scrollIntoView(),0);}}
document.querySelectorAll('.info').forEach(button=>{
 const tip=button.querySelector('.tip');
 const place=()=>{const box=button.getBoundingClientRect();tip.style.position='fixed';tip.style.right='auto';tip.style.top=(box.bottom+8)+'px';tip.style.left=Math.max(8,Math.min(window.innerWidth-300,box.right-280))+'px';};
 button.addEventListener('mouseenter',place);button.addEventListener('focus',place);
});
document.querySelectorAll('[data-chart]').forEach(card => {
 const svg=card.querySelector('svg'), viewport=card.querySelector('[data-viewport]');
 let zoom=1;
 const pan=card.querySelector('input'), labels=card.querySelectorAll('[data-tick]');
 const min=Number(card.dataset.min), max=Number(card.dataset.max), names=JSON.parse(card.dataset.labels || '[]');
 const show=()=>{
   const width=540/zoom, start=(540-width)*Number(pan.value)/100;
   viewport.setAttribute('viewBox',`${start} 0 ${width} 220`);
   labels.forEach((label,i)=>{
     const value=min+(max-min)*(start+width*i/4)/540;
     label.textContent=names.length ? (names[Math.max(0,Math.min(names.length-1,Math.round(value)))] || '') : (card.dataset.scale==='symlog-rate' ? (100*(value<=.01?value:.01*Math.exp(value/.01-1))).toFixed(value<.01?2:1)+'%' : Number(value.toPrecision(4)).toLocaleString());
   });
   card.querySelector('[data-zoom-label]').textContent=zoom+'×';
 };
 card.querySelectorAll('[data-zoom]').forEach(button=>button.addEventListener('click',()=>{
   zoom=button.dataset.zoom==='reset'?1:Math.max(1,Math.min(16,zoom*(button.dataset.zoom==='in'?2:.5)));
   if(zoom===1)pan.value=0; show();
 }));
 pan.addEventListener('input',show);
 card.querySelector('[data-export]').addEventListener('click',()=>{
   const copy=svg.cloneNode(true);copy.setAttribute('xmlns','http://www.w3.org/2000/svg');
   const blob=new Blob([new XMLSerializer().serializeToString(copy)],{type:'image/svg+xml'});
   const url=URL.createObjectURL(blob), a=document.createElement('a'); a.href=url;a.download=card.dataset.chart+'.svg';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);
 });
});
let saved=[];
window.addEventListener('beforeprint',()=>{saved=[...document.querySelectorAll('details')].map(d=>[d,d.open]);saved.forEach(([d])=>d.open=true);});
window.addEventListener('afterprint',()=>saved.forEach(([d,open])=>d.open=open));
})();
"""
