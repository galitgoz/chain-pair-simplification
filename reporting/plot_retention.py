from pathlib import Path
import json,collections
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
root=Path(__file__).resolve().parents[1];out=root/'reproduced/figures';out.mkdir(parents=True,exist_ok=True)
rows=json.loads((root/'results/plot_data.json').read_text(encoding='utf-8'))
pairs=collections.defaultdict(dict)
for r in rows:
    if r['group']=='השוואה ראשית':pairs[(r['domain'],r['pair'])][r['method']]=r
good={k:v for k,v in pairs.items() if len(v)==4 and all(r['validated']=='True' for r in v.values())}
data=[]
for (domain,pair),v in good.items():
    def retention(method):
        r=v[method];return 100*(float(r['kA'])+float(r['kB']))/(float(r['nA'])+float(r['nB']))
    a,b=retention('CPS-2F'),retention('CPS-3F');data.append(dict(domain=domain,pair=pair,alias=v['CPS-2F'].get('pair_alias') or pair,continuous=a,discrete=b,gap=b-a))
fig=plt.figure(figsize=(10.5,10),layout='constrained');gs=fig.add_gridspec(1,2,width_ratios=[1,1]);left=gs[0].subgridspec(31,1)
axp=fig.add_subplot(left[:12,0]);axh=fig.add_subplot(gs[1])
for ax,domain,title in [(axp,'protein','Proteins · 12 pairs'),(axh,'hurricane','Hurricanes · 31 pairs')]:
    a=sorted([r for r in data if r['domain']==domain],key=lambda r:(-r['gap'],r['pair']))
    for i,r in enumerate(a):
        ax.plot([r['continuous'],r['discrete']],[i,i],color='#b5bec7',lw=2.4,zorder=1)
        ax.scatter(r['discrete'],i,color='#d48826',marker='s',s=40,zorder=2)
        ax.scatter(r['continuous'],i,facecolors='white',edgecolors='#238448',marker='o',s=55,linewidths=1.8,zorder=3)
    ax.set_yticks(range(len(a)),[r['alias'] for r in a]);ax.set_ylim(len(a)-.35,-.65);ax.set_xlim(0,103);ax.set_xticks([0,20,40,60,80,100]);ax.grid(axis='x',alpha=.18);ax.set_axisbelow(True);ax.set_title(title,fontsize=14,pad=12);ax.set_xlabel('Input vertices retained (%)',fontsize=11);ax.tick_params(labelsize=10);ax.spines[['top','right','left']].set_visible(False);ax.tick_params(axis='y',length=0)
fig.suptitle('Compression under continuous and discrete fidelity',fontsize=17)
legend=fig.add_subplot(left[15:23,0]);legend.axis('off')
legend.legend(handles=[Line2D([],[],marker='o',markerfacecolor='white',markeredgecolor='#238448',markeredgewidth=1.8,color='none',label='CPS-2F · continuous fidelity',markersize=8),Line2D([],[],marker='s',color='none',markerfacecolor='#d48826',markeredgecolor='#d48826',label='CPS-3F · discrete fidelity',markersize=7)],loc='upper left',frameon=False,fontsize=11)
legend.text(.02,.55,'Further left = fewer vertices retained\nLonger line = greater difference between methods\n\nOne row per pair; sorted by percentage-point gap.\nBoth outputs satisfy their respective constraints.\n\nRetained = 100 × (kA + kB) / (nA + nB)\nCommon comparison set; α = 0.5',transform=legend.transAxes,va='top',fontsize=10,linespacing=1.65,color='#405265')
for ext in ['png','pdf']:fig.savefig(out/f'vertices_retained.{ext}',dpi=175)
plt.close(fig)
(out/'pair_mapping_and_values.json').write_text(json.dumps(data,indent=2),encoding='utf-8')
print('Created figure for 43 pairs, with original pair aliases and exact values saved.')
