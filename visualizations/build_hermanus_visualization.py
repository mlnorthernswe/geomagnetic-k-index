"""Build an offline interactive LRNS explainer and export PNG/SVG figures.

Only this optional visualization builder needs matplotlib. Calculators remain
standard-library Python. Run from any directory with a Python containing it.
"""
import json
import sys
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
import calculate_k_index as ih
import calculate_hermanus_k_index as hermanus

OUTPUT = Path(__file__).resolve().parent
COLORS = dict(raw='#a4b1c2', quiet='#13a88a', forward='#537be0', reverse='#b668cc', average='#df9c30', residual='#cc5562')


def first_pass(signal, mu, qstart=300, qend=1009, seed=1):
    """Trace the first L1 cascade only, using the reference's M-window settings."""
    f=hermanus._f32
    values=[f(v) for v in signal]
    total=0.
    for value in values: total=f(total+value)
    mean=f(total/1440)
    raw=[f(v-mean) for v in values]
    total=0.
    for v in raw: total=f(total+v)
    w=f(total/1440)
    guide=[];state=seed
    for _ in raw:
        state=16807*state%2147483647
        guide.append(f(1+f(.005)*(2*state/2147483647-1)))
    forward=[];weights=[];next_weights=[];steps=[]
    for i,(value,g) in enumerate(zip(raw,guide),1):
        step=f(f(2*f(mu))/3) if i<=qstart or i>=qend else f(2*f(mu))
        steps.append(step);weights.append(w)
        estimate=f(w*g);forward.append(estimate)
        error=f(value-estimate)
        sign=0 if abs(error)<f(.001) else 1 if error>0 else -1
        w=f(w+f(step*sign*g));next_weights.append(w)
    reverse=[0.]*1440
    # M mode carries the forward end weight into the reverse operation.
    for i in range(1439,-1,-1):
        reverse[i]=f(w*guide[i])
        error=f(raw[i]-reverse[i])
        sign=0 if abs(error)<f(.001) else 1 if error>0 else -1
        w=f(w+f(steps[i]*sign*guide[i]))
    average=[f(f(a+b)/2) for a,b in zip(forward,reverse)]
    return dict(raw=raw,forward=forward,reverse=reverse,average=average,
                weights=weights,nextWeights=next_weights,guide=guide,steps=steps)


def build_data():
    record=ih.parse_month_file(ROOT/'KIR/KIR23DEC.BIN')[0]
    day=datetime.strptime(str(record[1]),'%Y%j').date()
    k9=record[10]
    source={c:[v/10 for v in record[start:start+1440]] for c,start in ih.COMPONENT_START_WORDS.items()}
    data=dict(date=day.isoformat(),station='Kiruna',k9=k9,qfactor=k9*.01,mu={})
    for mu in (.08,.17,.20):
        results={};component_k=[]
        for c,signal in source.items():
            trace=first_pass(signal,mu)
            residual,smooth=hermanus.Hermanus().calculate(signal,mu=mu,qfactor=k9*.01)
            assert trace['raw']==[a+b for a,b in zip(residual,smooth)]
            trace.update(quiet=smooth,residual=residual)
            ranges=[max(residual[i:i+180])-min(residual[i:i+180]) for i in range(0,1440,180)]
            trace['ranges']=ranges
            trace['k']=[ih.classify_range_pt(v*1000,k9) for v in ranges]
            component_k.append(trace['k']);results[c]=trace
        results['k']=[max(a,b) for a,b in zip(*component_k)]
        # Cross-check bar values against the public calculator API.
        rows=hermanus.calculate_days({day:(list(record[16:1456]),list(record[1456:2896]),record,None)},mu=mu)
        assert results['k']==[row['k'] for row in rows]
        data['mu'][f'{mu:.2f}']=results
    return data


def export_figure(data):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                         'axes.edgecolor':'#d6dde6','axes.labelcolor':'#4a596b','xtick.color':'#657588','ytick.color':'#657588',
                         'text.color':'#243348','axes.titleweight':'bold','svg.fonttype':'none'})
    fig,axes=plt.subplots(4,1,figsize=(13,12),gridspec_kw={'height_ratios':[1.2,1.2,1,.7]},layout='constrained')
    fig.set_facecolor('#f8fafc')
    fig.suptitle('How Hermanus LRNS separates daily variation from magnetic activity',fontsize=17,fontweight='bold',x=.06,ha='left')
    t=[i/60 for i in range(1440)];result=data['mu']['0.17'];trace=result['x']
    titles=['1  Measured X and the final quiet-day estimate',
            '2  Inside L1: the first forward and backward pass',
            '3  Subtract the final quiet curve to obtain the disturbance',
            '4  Measure each 3-hour residual range and classify K']
    for ax,title in zip(axes,titles):
        ax.set_facecolor('white');ax.set_title(title,loc='left',pad=13)
        ax.set_xlim(0,24);ax.set_xticks(range(0,25,3));ax.grid(axis='y',color='#edf1f5');ax.set_axisbelow(True)
        for hour in range(3,24,3): ax.axvline(hour,color='#e6ebf1',linewidth=.8,zorder=0)
    axes[0].axvspan(5,16.8,color=COLORS['quiet'],alpha=.06)
    axes[0].plot(t,trace['raw'],color=COLORS['raw'],lw=.8,label='Measured X, daily mean removed')
    axes[0].plot(t,trace['quiet'],color=COLORS['quiet'],lw=2,label='Final quiet curve after L1 + L2')
    axes[0].set_ylabel('Departure (nT)');axes[0].legend(loc='upper right',frameon=False,fontsize=9)
    axes[1].plot(t,trace['raw'],color=COLORS['raw'],lw=.5,alpha=.32)
    for key,label in [('forward','Forward →'),('reverse','Backward ←'),('average','Average of this first pass')]:
        axes[1].plot(t,trace[key],color=COLORS[key],lw=1.5 if key=='average' else 1,label=label)
    axes[1].scatter([0],[trace['forward'][0]],color=COLORS['forward'],s=35,zorder=5)
    axes[1].annotate('First estimate starts near zero',xy=(0,trace['forward'][0]),xytext=(2.3,20),
                     textcoords='data',fontsize=9,arrowprops={'arrowstyle':'->','color':'#657588'})
    axes[1].set_ylabel('Departure (nT)');axes[1].legend(loc='upper right',frameon=False,fontsize=9)
    axes[1].text(.01,.03,'First rough pass only; a second L1 cascade and two L2 cascades follow.',transform=axes[1].transAxes,fontsize=9,color='#657588')
    axes[2].plot(t,trace['residual'],color=COLORS['residual'],lw=.9)
    axes[2].axhline(0,color='#657588',linewidth=.7);axes[2].set_ylabel('Residual (nT)')
    # Show one amplitude explicitly to connect the curve to the K bars.
    i=2;section=trace['residual'][i*180:(i+1)*180];lo,hi=min(section),max(section)
    axes[2].vlines(7.6,lo,hi,color='#243348',linewidth=1.3)
    axes[2].annotate(f'Range = {hi-lo:.1f} nT',xy=(7.6,(hi+lo)/2),xytext=(8.1,(hi+lo)/2),fontsize=9,va='center')
    centers=[1.5+3*i for i in range(8)]
    axes[3].bar(centers,result['k'],width=2.35,color='#405673',alpha=.95)
    for x,k,kx,ky in zip(centers,result['k'],result['x']['k'],result['y']['k']):
        axes[3].text(x,k+.15,str(k),ha='center',fontweight='bold')
        axes[3].text(x,.25,f'X {kx} / Y {ky}',ha='center',va='bottom',fontsize=8,color='white' if k else '#657588')
    axes[3].set_ylim(0,10);axes[3].set_yticks([0,3,6,9]);axes[3].set_ylabel('Combined K');axes[3].set_xlabel('Time (UT)')
    fig.supxlabel('Kiruna · 1 December 2023 · minute means · μ = 0.17 · M window 300–1009 · K9 = 1500 nT\nExample Hermanus parameters require calibration for Kiruna; plotted K values are estimates.',fontsize=10,color='#657588')
    fig.savefig(OUTPUT/'hermanus.png',dpi=160,facecolor=fig.get_facecolor())
    fig.savefig(OUTPUT/'hermanus.svg',facecolor=fig.get_facecolor())
    plt.close(fig)


if __name__=='__main__':
    data=build_data()
    export_figure(data)
    template=(OUTPUT/'hermanus-template.html').read_text()
    (OUTPUT/'hermanus.html').write_text(template.replace('__DATA__',json.dumps(data,separators=(',',':'))))
    print('Created hermanus.html, hermanus.png, and hermanus.svg')
