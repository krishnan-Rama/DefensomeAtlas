// Minimal DOM shim: enough to execute every tab builder and catch real errors.
class N {
  constructor(t){ this.tag=t; this.nodeType=1; this.tagName=String(t).toUpperCase(); this.children=[]; this.attrs={}; this.style={}; this._txt='';
    this.dataset={}; this.classList={ _s:new Set(),
      add:(...c)=>c.forEach(x=>this.classList._s.add(x)),
      remove:(...c)=>c.forEach(x=>this.classList._s.delete(x)),
      toggle:(c,f)=>{ f===undefined? (this.classList._s.has(c)?this.classList._s.delete(c):this.classList._s.add(c)) : (f?this.classList._s.add(c):this.classList._s.delete(c)); },
      contains:c=>this.classList._s.has(c) }; }
  append(...k){ k.flat().forEach(x=>{ this.children.push(x);
      if(this.tag==='select' && x && x.tag==='option' && this.value===undefined)
        this.value=x.attrs.value; }); }
  setAttribute(k,v){ this.attrs[k]=v; }
  getAttribute(k){ return this.attrs[k]; }
  addEventListener(){}
  querySelectorAll(){ return []; }
  remove(){}
  closest(){ return null; }
  getBoundingClientRect(){ return {left:0,top:0,width:900,height:600}; }
  getContext(){ return { clearRect(){},fillRect(){},fillText(){},save(){},restore(){},
    translate(){},rotate(){},measureText(){return{width:10}}, set font(v){}, set fillStyle(v){},
    set textAlign(v){} }; }
  get textContent(){ return this._txt; } set textContent(v){ this._txt=v; }
  get innerHTML(){ return ''; } set innerHTML(v){ this.children=[]; }
  get offsetWidth(){ return 100; }
  get firstChild(){ return this.children[0]; }
  dispatchEvent(){}
  cloneNode(){ return this; }
  toBlob(cb){ cb(new Blob()); }
  click(){}
}
const nodes={};
global.document={ createElement:t=>new N(t), createElementNS:(ns,t)=>new N(t),
  createTextNode:v=>({nodeType:3,textContent:v}),
  querySelector:s=>{ if(!nodes[s]) nodes[s]=new N('div'); return nodes[s]; },
  querySelectorAll:()=>[], getElementById:()=>null, addEventListener(){}, body:new N('body'), head:null };
global.window={innerWidth:1400,scrollX:0,addEventListener(){}};
global.localStorage={getItem(){return null;},setItem(){}};
global.location={hash:'',};
global.Blob=function(){}; global.URL={createObjectURL:()=>'x',revokeObjectURL:()=>{}};
global.Event=function(){};
