from pathlib import Path
import json,collections
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parent/'saved'
out=root/'tradeoff_preview';out.mkdir(exist_ok=True)
rows=json.loads((root/'experiment_tables/all_runs_full_precision.json').read_text(encoding='utf-8'))
pairs=collections.defaultdict(dict)
for r in rows:
    if r['group']=='השוואה ראשית':pairs[(r['domain'],r['pair'])][r['method']]=r
good={k:v for k,v in pairs.items() if len(v)==4 and all(r['validated']=='True' for r in v.values())}
fig,axes=plt.subplots(1,2,figsize=(11,4.6),sharex=True,sharey=True,layout='constrained')
spec=[('Independent continuous','IC','#2878b5','x'),('CPS-3F','CPS-3F','#d48826','s'),('CPS-2F','CPS-2F','#238448','o')]
points=[]
for ax,domain,title in zip(axes,['protein','hurricane'],['Proteins (12 pairs)','Hurricanes (31 pairs)']):
    ax.axhspan(0,1,color='#238448',alpha=.045)
    ax.axhline(1,color='#353535',linestyle='--',linewidth=1.2,zorder=1)
    for method,label,color,marker in spec:
        data=[]
        for (dom,pair),v in good.items():
            if dom!=domain:continue
            r=v[method];x=100*(1-(float(r['kA'])+float(r['kB']))/(float(r['nA'])+float(r['nB'])));y=float(r['coupling_distance'])/float(r['delta3']);data.append((x,y));points.append(dict(domain=dom,pair=pair,method=method,compression_percent=x,coupling_ratio=y))
        x,y=zip(*data)
        if marker=='o':ax.scatter(x,y,s=62,facecolors='none',edgecolors=color,marker=marker,linewidths=1.6,label=label,zorder=4)
        else:ax.scatter(x,y,s=40,c=color,marker=marker,alpha=.85,label=label,zorder=3)
    ax.set_title(title,fontsize=13);ax.set_xlabel('Vertex reduction (%)',fontsize=12);ax.set_xlim(-3,100);ax.set_ylim(0,3);ax.grid(alpha=.17);ax.tick_params(labelsize=11)
axes[0].set_ylabel('Output coupling distance / δ₃',fontsize=12)
axes[0].legend(loc='upper left',frameon=False,fontsize=11)
axes[1].text(98,1.05,'Coupling threshold',ha='right',fontsize=10,color='#444')
fig.suptitle('Compression and coupling on the common comparison set · α = 0.5',fontsize=14)
for ext in ['png','pdf']:fig.savefig(out/f'compression_coupling.{ext}',dpi=180)
plt.close(fig)
(out/'plotted_values.json').write_text(json.dumps(points,indent=2),encoding='utf-8')
print('Saved plot for',len(good),'pairs; 129 method observations. No jitter or clipping; overlapping points retained.')
