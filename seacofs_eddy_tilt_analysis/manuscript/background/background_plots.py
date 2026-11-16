"""Original tagged Figure 11 presentation; local manuscript map helper only."""
import numpy as np
import matplotlib.pyplot as plt


def plot(planetary,summary,grid,threshold,settings):
    import census_plots
    fig,axes=plt.subplots(1,3,figsize=settings['size_inches'],constrained_layout=True)
    colors={'AE':'firebrick','CE':'royalblue'}
    levels={'surface':'Surface','200':'0–200 m','200_500':'200–500 m'}
    ax=axes[0];p=summary[summary.component.eq('east')]
    for offset,cyc in zip((-.08,.08),('AE','CE')):
        q=p[p.Cyc.eq(cyc)].set_index('level').reindex(levels);m=q['median']
        ax.errorbar(np.arange(3)+offset,m,[m-q.ci_low,q.ci_high-m],fmt='o-',capsize=3,color=colors[cyc],label=cyc)
    ax.axhline(0,color='k',lw=.8)
    ax.set(xticks=range(3),ylabel='Zonal background flow (m s$^{-1}$)')
    ax.set_xticklabels(levels.values(),rotation=15);ax.legend()
    # Preserve the original panel-label position without colliding with a top tick.
    ax.margins(y=settings['summary_y_margin'])
    col='surface_minus_deep_east_ms'
    limit=np.nanpercentile(abs(planetary[col]),settings['colour_percentile'])
    if not np.isfinite(limit) or limit<=0:raise ValueError('Cannot establish a positive map colour scale')
    for k,(ax,cyc) in enumerate(zip(axes[1:],('AE','CE'))):
        q=planetary[planetary.Cyc.eq(cyc)].sort_values(col,key=abs)
        sc=ax.scatter(q.xc,q.yc,c=q[col],s=9,alpha=.7,cmap='PiYG',vmin=-limit,vmax=limit)
        ax.contour(grid.X_grid,grid.Y_grid,grid.h,[threshold],colors='k',linewidths=.7)
        census_plots.lat_lon_contours(ax,grid,land_mask=True,land_clr=colors[cyc])
        ax.text(30,830,cyc,fontsize=20,weight='bold',color='w')
        ax.set(xlabel='x (km)',ylabel='y (km)' if k==0 else '')
        if k:ax.set_yticklabels([])
    fig.colorbar(sc,ax=axes[1:],label='Surface − 200–500 m zonal background flow (m s$^{-1}$)',shrink=.7)
    for ax,label in zip(axes,['a)','b)','c)']):
        ax.text(-.1,1.01,label,transform=ax.transAxes,ha='left',va='top',fontsize=12,fontweight='bold')
    return fig,float(limit)
