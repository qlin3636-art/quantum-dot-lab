import {solve2D} from './physics2d.mjs?v=20261028-ritz8';
self.onmessage=({data})=>{try{self.postMessage({id:data.id,result:solve2D(data.V,data.options)});}catch(e){self.postMessage({id:data.id,error:e.message});}};
