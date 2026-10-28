import {HBAR2_2ME} from './physics.mjs';
// Block imaginary-time iteration of the 5-point finite-difference Hamiltonian.
// Returns discrete-normalized states; physical psi = vector / dx.
export function solve2D(V,{n=24,span=16,mass=.13,states=4,iterations=3500,tolerance=1e-5}={}){
 if(V.length!==n*n||V.some(v=>!Number.isFinite(v))||mass<=0||span<=0)throw Error('Invalid grid');
 const dx=span/(n+1),t=HBAR2_2ME/(mass*dx*dx),vmin=Math.min(...V),vmax=Math.max(...V),dt=.9/(8*t+vmax-vmin);
 const H=(v,out)=>{for(let y=0;y<n;y++)for(let x=0;x<n;x++){const i=y*n+x;out[i]=(4*t+V[i])*v[i]-t*((x?v[i-1]:0)+(x<n-1?v[i+1]:0)+(y?v[i-n]:0)+(y<n-1?v[i+n]:0));}};
 let seed=4291;const random=()=>{seed=(Math.imul(seed,1664525)+1013904223)>>>0;return seed/4294967296-.5;};
 const waves=Array.from({length:states},()=>Float64Array.from({length:n*n},random)),work=new Float64Array(n*n);
 function ortho(){for(let k=0;k<states;k++){const a=waves[k];for(let pass=0;pass<2;pass++)for(let j=0;j<k;j++){const b=waves[j];let dot=0;for(let i=0;i<a.length;i++)dot+=a[i]*b[i];for(let i=0;i<a.length;i++)a[i]-=dot*b[i];}let norm=0;for(const v of a)norm+=v*v;norm=Math.sqrt(norm);for(let i=0;i<a.length;i++)a[i]/=norm;}}
 ortho();let energies=[],residuals=[],steps=0;
 for(let step=0;step<iterations;step++){
  for(const a of waves){H(a,work);for(let i=0;i<a.length;i++)a[i]-=dt*(work[i]-vmin*a[i]);}ortho();steps=step+1;
  if(step%50===49||step===iterations-1){energies=[];residuals=[];for(const a of waves){H(a,work);let e=0;for(let i=0;i<a.length;i++)e+=a[i]*work[i];let r=0;for(let i=0;i<a.length;i++)r+=(work[i]-e*a[i])**2;energies.push(e);residuals.push(Math.sqrt(r));}if(Math.max(...residuals)<tolerance)break;}
 }
 const order=energies.map((e,k)=>k).sort((a,b)=>energies[a]-energies[b]);
 return {energies:order.map(k=>energies[k]),waves:order.map(k=>Array.from(waves[k],v=>v/dx)),residuals:order.map(k=>residuals[k]),dx,n,span,steps,converged:Math.max(...residuals)<tolerance};
}
export function ledModel({drive=1,electron=1,hole=1,traps=.2,outcoupling=.2}={}){
 // Dimensionless teaching model, not a calibrated device simulation.
 const Je=drive*electron,Jh=drive*hole,total=Je+Jh,balance=total?2*Math.min(Je,Jh)/total:0;
 const trapLoss=traps,highDriveLoss=.12*drive*drive;
 const radiative=1/(1+trapLoss+highDriveLoss),eqe=balance*radiative*outcoupling;
 const flux=total/2*eqe;
 return {Je,Jh,balance,radiative,eqe,flux,trapLoss,highDriveLoss};
}
