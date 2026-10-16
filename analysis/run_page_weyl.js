// Load the page's Weyl detector (unchanged source) into Node and run it on parameter sets.
const fs=require('fs');
const html=fs.readFileSync(process.argv[2],'utf8');
const a=html.indexOf('/* ------------------------------------------------------------------ Weyl spectral block */');
const b=html.indexOf('function bandsOnly(');
const src=html.slice(a,b);
const pre=`const TAU=2*Math.PI,SQ3=Math.sqrt(3),D2R=Math.PI/180;const clamp=(x,a,b)=>x<a?a:x>b?b:x;let P={},ON={weyl:true};`;
const run=new Function('cases',pre+src+`
  const out=[];
  for(const c of cases){P={t:1,breath:0,flux:90,mA:0,tz:0,dz:.8,gapSel:'12',Nk:18,kzCut:0,...c};const t0=Date.now();weylCompute();
    out.push({case:c,nodes:WEYL.nodes.map(n=>({chi:n.chi,kz:+n.kz.toFixed(3),gap:n.gap})),neutral:WEYL.neutral.length,flicker:WEYL.flicker,extended:WEYL.extended,
      oldStyleCount:(()=>{let s=0;for(let q=0;q<48;q++)if(WEYL.Cs[q]!==WEYL.Cs[q+1]&&WEYL.gaps[q]>=GAPTOL&&WEYL.gaps[q+1]>=GAPTOL)s+=Math.abs(WEYL.Cs[q+1]-WEYL.Cs[q]);return s;})(),ms:Date.now()-t0});}
  return out;`);
const cases=JSON.parse(fs.readFileSync(process.argv[3],'utf8'));
fs.writeFileSync(process.argv[4],JSON.stringify(run(cases)));
