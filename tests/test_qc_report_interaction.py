"""Exercise the offline enhancement script with a small DOM harness in Node."""

import shutil
import subprocess

import pytest

from deamtools.qc.report_assets import JS


def test_navigation_zoom_reset_and_print_state(tmp_path):
    node = shutil.which("node")
    if node is None:
        pytest.skip("Node is optional; Python report tests still run")
    harness = r"""
const assert=require('node:assert/strict');
function element(extra={}) {return Object.assign({listeners:{},addEventListener(event,fn){this.listeners[event]=fn;}},extra);}
const fold=element({tagName:'DETAILS',open:false,parentElement:null});
const target=element({parentElement:fold});
const link=element({hash:'#target'});
const viewport={setAttribute(k,v){this[k]=v;}};
const pan=element({value:'0'});
const ticks=Array.from({length:5},()=>({textContent:''}));
const zoomLabel={};
const buttons=['in','out','reset'].map(zoom=>element({dataset:{zoom}}));
const exportButton=element();
const chart={dataset:{min:'0',max:'1',scale:'number'},querySelector(selector){
 if(selector==='svg')return {};
 if(selector==='[data-viewport]')return viewport;
 if(selector==='input')return pan;
 if(selector==='[data-zoom-label]')return zoomLabel;
 if(selector==='[data-export]')return exportButton;
 throw new Error(selector);
},querySelectorAll(selector){return selector==='[data-tick]'?ticks:buttons;}};
global.document={querySelectorAll(selector){
 if(selector==='a[href^="#"]')return [link];
 if(selector==='[data-chart]')return [chart];
 if(selector==='details')return [fold];
 if(selector==='.info')return [];
 throw new Error(selector);
},getElementById(id){return id==='target'?target:null;}};
global.location={hash:''};
global.window=element();
"""
    assertions = r"""
link.listeners.click();assert.equal(fold.open,true);
buttons[0].listeners.click();assert.equal(zoomLabel.textContent,'2×');
assert.equal(viewport.viewBox,'0 0 270 220');
pan.value='100';pan.listeners.input();assert.equal(viewport.viewBox,'270 0 270 220');
assert.equal(ticks[0].textContent,'0.5');assert.equal(ticks[4].textContent,'1');
buttons[2].listeners.click();assert.equal(viewport.viewBox,'0 0 540 220');assert.equal(pan.value,0);
fold.open=false;window.listeners.beforeprint();assert.equal(fold.open,true);
window.listeners.afterprint();assert.equal(fold.open,false);
"""
    path = tmp_path / "interaction.cjs"
    path.write_text(harness + JS + assertions)
    subprocess.run([node, str(path)], check=True, capture_output=True, text=True)
