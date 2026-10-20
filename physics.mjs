// Units: nm, eV, effective mass / electron rest mass.
export const HBAR2_2ME = 0.0380998212;
export const HC = 1239.841984;
export const COULOMB = 1.43996455;
export function dotEmission(radius, material = 'CdSe') {
  const p = material === 'InP' ? {gap:1.344,me:0.077,mh:0.60,eps:12.4} : {gap:1.74,me:0.13,mh:0.45,eps:9.5};
  const confinement = HBAR2_2ME*Math.PI**2/radius**2*(1/p.me+1/p.mh);
  const attraction = 1.8*COULOMB/(p.eps*radius);
  const energy = p.gap+confinement-attraction;
  return {energy,wavelength:HC/energy,confinement,attraction,...p};
}
export function wavelengthRGB(w) {
  let r=0,g=0,b=0;
  if(w>=380&&w<440){r=(440-w)/60;b=1;}else if(w<490&&w>=440){g=(w-440)/50;b=1;}else if(w<510&&w>=490){g=1;b=(510-w)/20;}else if(w<580&&w>=510){r=(w-510)/70;g=1;}else if(w<645&&w>=580){r=1;g=(645-w)/65;}else if(w<=780&&w>=645){r=1;}
  const a=w<420?0.3+0.7*(w-380)/40:w>700?0.3+0.7*(780-w)/80:1;
  return [r,g,b].map(v=>Math.round(255*Math.max(0,v*Math.max(0,a))**0.8));
}
export function parsePotential(input) {
  if(input.length>320) throw Error('expression_length');
  const tokens=input.match(/(?:\d*\.?\d+(?:[eE][+-]?\d+)?)|[A-Za-z_]+|[+\-*/^(),]/g)||[];
  if(tokens.join('')!==input.replace(/\s/g,'')) throw Error('expression_syntax');
  let at=0;
  const funcs={sin:Math.sin,cos:Math.cos,exp:Math.exp,abs:Math.abs,sqrt:Math.sqrt,tanh:Math.tanh,min:Math.min,max:Math.max};
  function primary(){const t=tokens[at++];if(t==='('){const a=sum();if(tokens[at++]!==')')throw Error('expression_syntax');return a;}if(t==='x')return x=>x;if(t==='pi')return ()=>Math.PI;if(t==='e')return ()=>Math.E;if(t in funcs){if(tokens[at++]!=='(')throw Error('expression_syntax');const args=[sum()];while(tokens[at]===','){at++;args.push(sum());}if(tokens[at++]!==')'||args.length>2||((t==='min'||t==='max')?args.length!==2:args.length!==1))throw Error('expression_syntax');return x=>funcs[t](...args.map(a=>a(x)));}if(t!==undefined&&/^\d*\.?\d/.test(t)){const n=Number(t);return ()=>n;}throw Error('expression_syntax');}
  function power(){let a=primary();if(tokens[at]==='^'){at++;const b=unary();const old=a;a=x=>old(x)**b(x);}return a;}
  function unary(){if(tokens[at]==='+'){at++;return unary();}if(tokens[at]==='-'){at++;const a=unary();return x=>-a(x);}return power();}
  function product(){let a=unary();while(['*','/'].includes(tokens[at])){const op=tokens[at++],b=unary(),old=a;a=x=>op==='*'?old(x)*b(x):old(x)/b(x);}return a;}
  function sum(){let a=product();while(['+','-'].includes(tokens[at])){const op=tokens[at++],b=product(),old=a;a=x=>op==='+'?old(x)+b(x):old(x)-b(x);}return a;}
  const fn=sum();if(at!==tokens.length)throw Error('expression_syntax');return fn;
}
// Symmetric tridiagonal Hamiltonian. Sturm bisection + inverse iteration.
export function solveWell(potential,{mass=0.13,span=20,states=5}={}) {
  const n=potential.length,dx=span/(n+1),t=HBAR2_2ME/(mass*dx*dx);
  if(n<8||!Number.isFinite(t)||mass<=0||potential.some(v=>!Number.isFinite(v)))throw Error('invalid_potential');
  const diag=potential.map(v=>v+2*t),off=-t,small=1e-12;
  const count=value=>{let q=diag[0]-value,c=q<0?1:0;for(let i=1;i<n;i++){if(Math.abs(q)<small)q=-small;q=diag[i]-value-off*off/q;if(q<0)c++;}return c;};
  const min=Math.min(...potential),max=Math.max(...potential)+4*t;
  const energies=[];for(let k=0;k<states;k++){let lo=min,hi=max;for(let j=0;j<65;j++){const mid=(lo+hi)/2;if(count(mid)<=k)lo=mid;else hi=mid;}energies.push((lo+hi)/2);}
  const waves=[];
  for(let k=0;k<states;k++){
    let v=Array.from({length:n},(_,i)=>Math.sin((i+1)*(k+1)*Math.PI/(n+1))+0.13*Math.cos(i*1.37));
    const shift=energies[k]+1e-8;
    for(let step=0;step<7;step++){
      const c=new Float64Array(n),d=new Float64Array(n);let pivot=diag[0]-shift;if(Math.abs(pivot)<small)pivot=small;c[0]=off/pivot;d[0]=v[0]/pivot;
      for(let i=1;i<n;i++){pivot=diag[i]-shift-off*c[i-1];if(Math.abs(pivot)<small)pivot=small;c[i]=off/pivot;d[i]=(v[i]-off*d[i-1])/pivot;}
      v[n-1]=d[n-1];for(let i=n-2;i>=0;i--)v[i]=d[i]-c[i]*v[i+1];
      for(const w of waves){let dot=0;for(let i=0;i<n;i++)dot+=v[i]*w[i]*dx;for(let i=0;i<n;i++)v[i]-=dot*w[i];}
      const norm=Math.sqrt(v.reduce((s,a)=>s+a*a*dx,0));v=v.map(a=>a/norm);
    }
    const peak=v.find(a=>Math.abs(a)>0.05);if(peak<0)v=v.map(a=>-a);waves.push(v);
  }
  return {energies,waves,dx,x:potential.map((_,i)=>-span/2+(i+1)*dx),t};
}
