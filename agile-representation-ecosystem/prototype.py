import numpy as np
from dataclasses import dataclass

@dataclass
class Specialist:
    name: str
    projection: np.ndarray
    reliability: float
    cost: float

class Ecosystem:
    def __init__(self, seed=42, mode='adaptive'):
        self.rng=np.random.default_rng(seed); self.mode=mode
        self.specialists=[Specialist('causal',np.array([1,0,0,0]),.84,1),Specialist('spectral',np.array([0,1,0,0]),.84,1),Specialist('experimental',np.array([0,0,1,0]),.84,1),Specialist('structural',np.array([0,0,0,1]),.84,1)]
        self.counts=np.zeros(4); self.beliefs=np.zeros(4); self.history=[]; self.artifacts=[]
    def select(self, need, t):
        if self.mode=='fixed': return 0
        if self.mode=='round_robin': return t%4
        if self.mode=='random': return int(self.rng.integers(4))
        scores=need + .2/np.sqrt(1+self.counts)-.12*self.counts/max(1,t) + .08*np.cos(.8*t+np.arange(4)*np.pi/2)
        temp=.14 if self.mode=='adaptive' else .001
        probs=np.exp((scores-scores.max())/temp);probs/=probs.sum()
        return int(self.rng.choice(4,p=probs))
    def run(self, tasks):
        for t,(need,truth) in enumerate(tasks):
            i=self.select(need,t);self.counts[i]+=1
            signal=truth[i]+self.rng.normal(0,.17)
            self.beliefs[i]=.75*self.beliefs[i]+.25*signal
            # Artifact: independently verifiable observation, retained with provenance
            artifact={'step':t,'specialist':self.specialists[i].name,'dimension':i,'observation':round(float(signal),3),'verified':abs(signal-truth[i])<.35}
            self.artifacts.append(artifact)
            coverage=sum((self.counts>0).astype(float)*need)/sum(need)
            self.history.append({'step':t,'active':self.specialists[i].name,'need':need.tolist(),'coverage':coverage,'verified':artifact['verified']})
        return {'coverage':float(np.mean([x['coverage'] for x in self.history])),'verified':sum(x['verified'] for x in self.history),'specialists_used':int(sum(self.counts>0))}

def tasks(seed=123,n=100):
    rng=np.random.default_rng(seed); out=[]
    for t in range(n):
        need=np.full(4,.08); active=(t//15)%4;need[active]=1
        if t%19 in (0,1,2):need[(active+2)%4]=.65
        truth=rng.normal(size=4)
        out.append((need,truth))
    return out

if __name__=='__main__':
    ts=tasks()
    for mode in ['fixed','round_robin','random','context','adaptive']:
        e=Ecosystem(mode=mode);print(mode,e.run(ts))
