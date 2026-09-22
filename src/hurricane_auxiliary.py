"""Explain the CPS-2F auxiliary construction using actual hurricane coordinates."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.ticker import MaxNLocator
from cps_paper_algorithms import extended_curve, EPS
from hurricane_cps import save_figure


def construction_events(P, delta):
    """Record each accepted circle/edge root, before parameter deduplication."""
    rows=[]
    for edge in range(len(P)-1):
        start,end=P[edge:edge+2];v=end-start;aa=float(v@v)
        if aa<1e-24:continue
        for center,c in enumerate(P):
            q=start-c;t0=-float(q@v)/aa
            slack=delta**2-float((q+t0*v)@(q+t0*v))
            if slack < -1e-10:continue
            h=np.sqrt(max(0.,slack)/aa)
            for branch,t in [('minus',t0-h),('plus',t0+h)]:
                if EPS<t<1-EPS:
                    point=start+t*v
                    rows.append(dict(edge=edge,center_vertex=center,branch=branch,t=t,
                        curve_parameter=edge+t,x_km=point[0],y_km=point[1],
                        radius_km=delta,distance_to_center_km=float(np.linalg.norm(point-c)),
                        quadratic_a=aa,quadratic_b=2*float(q@v),quadratic_c=float(q@q)-delta**2))
    events=pd.DataFrame(rows,columns=['edge','center_vertex','branch','t','curve_parameter',
        'x_km','y_km','radius_km','distance_to_center_km','quadratic_a','quadratic_b','quadratic_c'])
    X,parameters=extended_curve(P,delta)
    # Identical curve-parameter deduplication policy to the actual solver.
    raw=sorted(list(map(float,range(len(P))))+events.curve_parameter.tolist())
    unique=[]
    for t in raw:
        if not unique or t-unique[-1]>1e-10:unique.append(t)
    assert np.allclose(unique,parameters,rtol=0,atol=1e-10)
    mask=np.abs(parameters-np.round(parameters))>1e-10
    if len(events):assert np.allclose(events.distance_to_center_km,delta,rtol=1e-9,atol=1e-7)
    assert int(mask.sum())==len(X)-len(P)
    return events,X,parameters,mask


def point_catalog(data,parameters,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    points=[];summaries=[];event_tables=[]
    for r in parameters.itertuples():
        item=data['instances'][r.pair+'@1']
        for curve,P,delta in [('A',item['A'],r.delta1),('B',item['B_xyz'],r.delta2)]:
            events,X,t,aux=construction_events(P,delta)
            for i in range(len(X)):
                points.append(dict(pair=r.pair,alpha=r.alpha,curve=curve,extended_index=i,
                    curve_parameter=t[i],point_type='auxiliary' if aux[i] else 'original',x_km=X[i,0],y_km=X[i,1]))
            event_tables.append(events.assign(pair=r.pair,alpha=r.alpha,curve=curve))
            summaries.append(dict(pair=r.pair,alpha=r.alpha,curve=curve,n_original=len(P),
                delta1=r.delta1,delta2=r.delta2,delta3=r.delta3,
                accepted_root_events=len(events),duplicate_events_removed=len(events)-int(aux.sum()),
                added_auxiliary=int(aux.sum()),extended_size=len(X),auxiliary_upper=2*len(P)*(len(P)-1),
                utilization_pct=100*int(aux.sum())/(2*len(P)*(len(P)-1))))
    summary=pd.DataFrame(summaries)
    summary.to_csv(out/'auxiliary_construction_summary.csv',index=False)
    pd.DataFrame(points).to_csv(out/'extended_curve_points.csv',index=False)
    pd.concat(event_tables,ignore_index=True).to_csv(out/'auxiliary_intersection_events.csv',index=False)
    return summary


def worked_intersection(item,delta,out):
    P=item['A'];events,X,t,aux=construction_events(P,delta)
    event=events.iloc[(events.t-.5).abs().argmin()]
    edge=int(event.edge);center=int(event.center_vertex);q=P[edge]-P[center];v=P[edge+1]-P[edge]
    report=pd.DataFrame([dict(curve=item['pair'].split('__')[0],edge=f'{edge} → {edge+1}',
        center_vertex=center,radius_km=delta,a=event.quadratic_a,b=event.quadratic_b,c=event.quadratic_c,
        discriminant=event.quadratic_b**2-4*event.quadratic_a*event.quadratic_c,
        accepted_t=event.t,curve_parameter=event.curve_parameter,x_km=event.x_km,y_km=event.y_km,
        distance_to_center_km=event.distance_to_center_km)])
    report.to_csv(Path(out)/'worked_intersection.csv',index=False)
    fig,axes=plt.subplots(1,2,figsize=(12,5))
    ax=axes[0];ax.plot(P[:,0],P[:,1],'-',color='.7',lw=1)
    ax.plot(P[edge:edge+2,0],P[edge:edge+2,1],'-o',color='#4477AA',lw=3,label='Chosen input edge')
    ax.add_patch(Circle(P[center],delta,fill=False,color='#AA3377',ls='--',lw=1.5))
    ax.scatter(*P[center],marker='*',s=130,color='#AA3377',label=f'Circle centre: vertex {center}')
    ax.scatter(event.x_km,event.y_km,s=90,marker='x',color='#CC6600',lw=2,label='Accepted auxiliary point')
    ax.plot([P[center,0],event.x_km],[P[center,1],event.y_km],':',color='#AA3377')
    selected=np.vstack([P[edge:edge+2],P[center]+delta,P[center]-delta]);mid=.5*(selected.min(0)+selected.max(0));half=.6*np.ptp(selected,axis=0).max()
    ax.set_xlim(mid[0]-half,mid[0]+half);ax.set_ylim(mid[1]-half,mid[1]+half)
    ax.set_aspect('equal');ax.set_title('1. Intersect an edge with one radius-δ₁ circle')
    ax.set_xlabel('Projected easting (km)');ax.set_ylabel('Projected northing (km)');ax.legend(fontsize=8)
    ax=axes[1]
    region=(t>=edge)&(t<=edge+1)
    ax.hlines(0,edge,edge+1,color='.65',lw=2)
    ax.scatter(t[region&aux],np.zeros((region&aux).sum()),marker='x',color='#CC6600',s=65,label='All added points on this edge')
    ax.scatter([edge,edge+1],[0,0],color='#4477AA',s=65,label='Original endpoints')
    ax.annotate(f'Chosen root: t={event.t:.4f}',(event.curve_parameter,0),xytext=(edge+.1,.2),arrowprops=dict(arrowstyle='->'))
    ax.set_ylim(-.1,.35);ax.set_yticks([]);ax.set_xlabel('Curve parameter s = edge index + t')
    ax.set_title('2. Sort and deduplicate along the curve');ax.legend(loc='upper left',fontsize=8)
    fig.suptitle(f"Worked example: {item['pair'].replace('__',' / ')}; A, δ₁={delta:.3f} km")
    fig.tight_layout();save_figure(fig,Path(out)/'figures','auxiliary_worked_example')
    return report


def auxiliary_curve_plots(data,parameters,out,alpha=1):
    """Show all points plus a circle-centred local view for A and B of each pair."""
    for r in parameters[parameters.alpha.eq(alpha)].itertuples():
        item=data['instances'][r.pair+'@1']
        fig,axes=plt.subplots(2,2,figsize=(12,10))
        for col,curve,P,delta in [(0,'A',item['A'],r.delta1),(1,'B',item['B_xyz'],r.delta2)]:
            events,X,t,aux=construction_events(P,delta)
            center=int(events.center_vertex.value_counts().index[0]) if len(events) else 0
            roots=events[events.center_vertex.eq(center)]
            for ax in axes[:,col]:
                ax.plot(P[:,0],P[:,1],color='.55',lw=1,zorder=1)
                ax.scatter(P[:,0],P[:,1],s=26,facecolors='white',edgecolors='#4477AA',lw=1.2,zorder=3,label='Original vertices')
                ax.scatter(X[aux,0],X[aux,1],s=25,marker='x',color='#CC6600',lw=1.1,zorder=4,label='Added auxiliary points')
                ax.set_aspect('equal',adjustable='box');ax.grid(alpha=.15)
                ax.xaxis.set_major_locator(MaxNLocator(4));ax.yaxis.set_major_locator(MaxNLocator(5))
                ax.set_xlabel('Projected easting (km)');ax.set_ylabel('Projected northing (km)')
            axes[0,col].set_title(f'Curve {curve}: all {len(P)} original + {int(aux.sum())} auxiliary points\nExtended size {len(X)}; bound {2*len(P)*(len(P)-1)} added points',fontsize=10)
            mid=.5*(P.min(0)+P.max(0));half=.55*np.ptp(P,axis=0).max()
            axes[0,col].set_xlim(mid[0]-half,mid[0]+half);axes[0,col].set_ylim(mid[1]-half,mid[1]+half)
            ax=axes[1,col]
            ax.add_patch(Circle(P[center],delta,fill=False,color='#AA3377',ls='--',lw=1.5))
            ax.scatter(*P[center],marker='*',s=120,color='#AA3377',zorder=5)
            ax.scatter(roots.x_km,roots.y_km,s=95,facecolors='none',edgecolors='#AA3377',lw=1.4,zorder=6)
            ax.set_xlim(P[center,0]-1.5*delta,P[center,0]+1.5*delta)
            ax.set_ylim(P[center,1]-1.5*delta,P[center,1]+1.5*delta)
            ax.set_title(f'Local construction: centre vertex {center}, radius {delta:.2f} km\nPurple rings: intersections generated by this circle',fontsize=10)
        handles,labels=axes[0,0].get_legend_handles_labels()
        fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,.93),ncol=2,frameon=False)
        fig.suptitle(f"{r.pair.replace('__',' / ')} — α={alpha}\nδ₁={r.delta1:.2f}, δ₂={r.delta2:.2f}, δ₃={r.delta3:.2f} km",fontsize=13)
        fig.text(.5,.012,'Auxiliary points lie on the existing track; they subdivide edges without changing the curve.\nOnly CPS-2F uses them as owner positions. They are not new storm observations or output candidates.',ha='center',fontsize=10)
        fig.subplots_adjust(left=.08,right=.98,bottom=.09,top=.85,hspace=.40,wspace=.26)
        save_figure(fig,Path(out)/'figures',f'auxiliary_on_curves_{r.pair}_alpha{alpha}')
