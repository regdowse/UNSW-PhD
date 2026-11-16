"""Original six-panel Figure 17 layout, colours and confidence-band style."""
REGIMES=['planetary','mixed','topographic']
COLORS={'AE':'firebrick','CE':'royalblue'}

def plot(stats, max_depth_m=1000., min_eddies=2, sparse_eddies=20, figsize=(13,6)):
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 3, figsize=figsize, sharey=True, constrained_layout=True)
    for i, (metric, label) in enumerate([('mean_along', 'Along PV gradient'),
                                       ('mean_perp', 'Perpendicular (+CCW)')]):
        for j, regime in enumerate(REGIMES):
            ax = axes[i, j]; found = False
            for cyc, color in COLORS.items():
                if stats.empty: continue
                s = stats.loc[stats.regime.eq(regime) & stats.Cyc.eq(cyc)].sort_values('Depth')
                if s.empty: continue
                ok = s.n_eddies.ge(min_eddies); found |= bool(ok.any())
                ax.plot(s[metric].where(ok), s.Depth, color=color, label=cyc)
                ax.fill_betweenx(s.Depth, s[metric+'_ci_low'].where(ok),
                                 s[metric+'_ci_high'].where(ok), color=color, alpha=.15)
                sparse = ok & s.n_eddies.lt(sparse_eddies)
                ax.scatter(s.loc[sparse, metric], s.loc[sparse, 'Depth'],
                           facecolors='none', edgecolors=color, s=18)
            ax.axvline(0, color='.4', lw=.7); ax.grid(alpha=.15)
            ax.set(xlabel=label+' displacement (km)', ylabel='Depth (m)' if j==0 else '',
                   title={'planetary':'Deep planetary','mixed':'Mixed','topographic':'Shallow topographic'}[regime] if i==0 else '', ylim=(max_depth_m, 0))
            if found & (i==1) & (j==2): ax.legend()
            # else: ax.text(.5, .5, 'Insufficient contributors', ha='center', transform=ax.transAxes)
        limit = max(1e-6, max(abs(v) for ax in axes[i] for v in ax.get_xlim()))
        for ax in axes[i]: ax.set_xlim(-limit, limit)
    for ax,l in zip(axes.flatten(),['a)','b)','c)','d)','e)','f)']):
        ax.text(-.1,1.09,l,transform=ax.transAxes,ha='left',va='top',
                fontsize=12,fontweight='bold')
    # fig.suptitle('Mean constituent displacement relative to environmental PV gradient\n'
    #              + f'95% track-bootstrap CI; open circles <{sparse_eddies} eddies')
    return fig, axes
