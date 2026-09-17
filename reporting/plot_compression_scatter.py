from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parent/'saved'/'compression_preview'
data=json.loads((root/'pair_mapping_and_values.json').read_text(encoding='utf-8'))
assert len(data)==43
fig,ax=plt.subplots(figsize=(6.6,6.4),layout='constrained')
ax.plot([0,100],[0,100],color='#697683',linestyle='--',linewidth=1.3,zorder=1,label='Equal compression')
for domain,label,color,marker in [('hurricane','Hurricanes (31 pairs)','#167e96','o'),('protein','Proteins (12 pairs)','#8753a6','^')]:
    a=[r for r in data if r['domain']==domain]
    ax.scatter([100-r['discrete'] for r in a],[100-r['continuous'] for r in a],s=65,c=color,marker=marker,alpha=.8,edgecolors='white',linewidths=.6,label=label,zorder=3)
ax.set_xlim(0,100);ax.set_ylim(0,100);ax.set_aspect('equal',adjustable='box')
ax.set_xticks(range(0,101,20));ax.set_yticks(range(0,101,20));ax.grid(alpha=.18);ax.set_axisbelow(True)
ax.set_xlabel('Vertices removed by CPS-3F (%)',fontsize=12,labelpad=9)
ax.set_ylabel('Vertices removed by CPS-2F (%)',fontsize=12,labelpad=9)
ax.set_title('Compression comparison · α = 0.5',fontsize=15,pad=14)
ax.legend(loc='lower right',frameon=True,facecolor='white',edgecolor='#ddd',fontsize=10)
ax.text(52,57,'Equal compression',rotation=45,color='#697683',fontsize=10,ha='center',rotation_mode='anchor')
ax.spines[['top','right']].set_visible(False)
for ext in ['png','pdf']:fig.savefig(root/f'compression_scatter.{ext}',dpi=190)
plt.close(fig)
print('Saved scatter plot for 43 pairs; no jitter, identical scales, exact saved values.')
