"""Artifact-first scientific agent ecosystem. Requires numpy; no API key."""
import json,hashlib,argparse
from pathlib import Path
import numpy as np
from dataclasses import dataclass,asdict
@dataclass
class Artifact:
    agent:str
    kind:str
    payload:dict
    parents:list
    sha256:str=''
    def finalize(self):
        body=json.dumps(dict(agent=self.agent,kind=self.kind,payload=self.payload,parents=self.parents),sort_keys=True)
        self.sha256=hashlib.sha256(body.encode()).hexdigest()
        return self
class Registry:
    def __init__(self):self.items={}
    def put(self,a):
        a.finalize();self.items[a.sha256]=a;return a.sha256
    def save(self,path):Path(path).write_text(json.dumps([asdict(a) for a in self.items.values()],indent=2))
FAMILIES={'linear':[0,3],'polynomial':[0,1,3],'harmonic':[0,2,3],'hybrid':[0,1,2,3]}
def features(x):return np.stack([x,x*x,np.sin(3*x),np.ones_like(x)],axis=1)
def predict(x,ix,coef):return features(x)[:,ix]@np.asarray(coef)
def dataset(seed=9,n=240):
    rng=np.random.default_rng(seed);x=rng.uniform(-2,2,n)
    y=.75*x+.65*x*x+.8*np.sin(3*x)+rng.normal(0,.07,n)
    return x,y
class Agent:
    def __init__(self,name):self.name=name
    def propose(self,x,y,r,parent):
        ix=FAMILIES[self.name];coef=np.linalg.lstsq(features(x)[:,ix],y,rcond=None)[0]
        return r.put(Artifact(self.name,'model',dict(family=self.name,features=ix,coefficients=coef.tolist()),[parent]))
class Verifier:
    def __init__(self,x,y):self.x=x;self.y=y
    def check(self,a,r):
        p=a.payload;err=float(np.mean((self.y-predict(self.x,p['features'],p['coefficients']))**2))
        vid=r.put(Artifact('independent_verifier','evaluation',dict(candidate=a.sha256,heldout_mse=err,verified=bool(np.isfinite(err))),[a.sha256]))
        return vid,err
class Coordinator:
    def __init__(self,r,verifier,mode):self.r=r;self.v=verifier;self.mode=mode;self.trace=[]
    def evaluate(self,name,x,y,parent):
        aid=Agent(name).propose(x,y,self.r,parent);vid,err=self.v.check(self.r.items[aid],self.r)
        self.trace.append(dict(need='explain observations' if name!='hybrid' else 'unexplained residual: expand representation',agent=name,artifact=aid,verification=vid,mse=err))
        return err,aid,vid,name
    def run(self,x,y,parent,threshold=.03):
        candidates=[self.evaluate(name,x,y,parent) for name in ('linear','polynomial','harmonic')]
        if min(candidates)[0]>threshold and self.mode=='need':
            candidates.append(self.evaluate('hybrid',x,y,parent))
        return min(candidates)
def experiment(seed=9,mode='need',out='artifacts.json'):
    x,y=dataset(seed);train=np.arange(len(x))%3!=0;test=~train
    r=Registry();root=r.put(Artifact('lab','dataset',dict(seed=seed,n_train=int(sum(train)),n_test=int(sum(test)),ground_truth='hidden from agents'),[]))
    c=Coordinator(r,Verifier(x[test],y[test]),mode);best=c.run(x[train],y[train],root)
    r.put(Artifact('coordinator','decision',dict(winner=best[3],heldout_mse=best[0],mode=mode,trace=c.trace),[best[1],best[2]]))
    if out:r.save(out)
    return dict(mode=mode,winner=best[3],heldout_mse=best[0],candidate_count=len(c.trace),artifact_count=len(r.items))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--seed',type=int,default=9);p.add_argument('--mode',choices=['need','fixed'],default='need');p.add_argument('--out',default='artifacts.json');args=p.parse_args()
    print(json.dumps(experiment(args.seed,args.mode,args.out),indent=2))
