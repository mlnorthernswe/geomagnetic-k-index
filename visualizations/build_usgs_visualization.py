"""Build an offline USGS explainer and PNG/SVG figures from real Kiruna data.

Only optional figure export needs matplotlib; the calculator remains standard-library Python.
"""
import json
import math
import sys
from datetime import timedelta
from pathlib import Path

ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT))
import calculate_k_index as ih
import calculate_usgs_k_index as usgs

OUTPUT=Path(__file__).resolve().parent
COLORS=dict(raw='#a4b1c2',accepted='#13a88a',rejected='#d98536',quiet='#13a88a',segment='#8c75c5',residual='#cc5562')


def build_data():
    december=usgs.iaf_days(ih.parse_month_file(ROOT/'KIR/KIR23DEC.BIN'))
    days=usgs.iaf_days(ih.parse_month_file(ROOT/'KIR/KIR23NOV.BIN'))
    days.update(december)
    dates=sorted(december)
    chosen=[dates[i] for i in (0,1,10)]
    data=dict(station='Kiruna',month='December 2023',k9=1500,configurations={})
    for nppm in (30,60):
        params=usgs.bootstrap_parameters(days,dates,nppm=nppm)
        encoded=usgs._encode(usgs._polar(days),params['base_h_nT'],params['base_d_arcmin'])
        means={day:[usgs._mean_day(values,params,c,day.isoformat() in params['quiet_days'])
                    for c,values in enumerate(arrays)] for day,arrays in encoded.items()}
        config=dict(parameters=params,histograms=[],days={})
        for c in (0,1):
            scale=usgs.F(params['xd']) if c else usgs.F(1)
            all_means=[v for d in dates for v in usgs._blocks(encoded[d][c],nppm)[0] if v is not None]
            low,high=min(all_means),max(all_means);width=max(1,(high-low+9)//10)
            counts=[0]*10
            for value in all_means: counts[min((value-low)//width,9)]+=1
            limits=params['components'][c]
            config['histograms'].append(dict(edges=[float(usgs.F(low+i*width)*.1*scale) for i in range(11)],counts=counts,
                low=float(usgs.F(limits['low'])*.1*scale),high=float(usgs.F(limits['high'])*.1*scale),
                rangeLimit=float(usgs.F(limits['range'])*.1*scale)))
        for day in chosen:
            results=[]
            rows=usgs.calculate_days(days,[day],parameters=params)
            for c in (0,1):
                values=encoded[day][c];scale=usgs.F(params['xd']) if c else usgs.F(1)
                limits=params['components'][c];quiet_day=day.isoformat() in params['quiet_days']
                block_means,block_ranges=usgs._blocks(values,nppm)
                accepted=means[day][c][0]
                blocks=[]
                for block,(calib_mean,span,selected) in enumerate(zip(block_means,block_ranges,accepted)):
                    good=[v for v in values[block*nppm:(block+1)*nppm] if v is not None]
                    mean=int(usgs._sum(good)/len(good)) if calib_mean is not None else None
                    reason=('Too few valid minute means' if mean is None else 'Within-block range exceeds limit'
                            if span is not None and span>limits['range'] else 'Mean exceeds upper limit'
                            if mean>limits['high'] else 'Mean is below lower limit'
                            if not quiet_day and mean<limits['low'] else 'Accepted quiet mean')
                    blocks.append(dict(time=(block+.5)*nppm/60,mean=None if mean is None else float(usgs.F(mean)*.1*scale),
                        range=None if span is None else float(usgs.F(span)*.1*scale),accepted=selected is not None,
                        reason=reason,coverage=len(good)/nppm))
                neighbors=[means.get(day+timedelta(days=d),[(None,None)]*2)[c] for d in (-1,0,1)]
                nm=1440//nppm
                triplet=[v if v is not None else [None]*nm for v,_ in neighbors]
                dailies=[v if v is not None else limits['end'] for _,v in neighbors]
                x,y=usgs._coordinates(triplet,dailies,nppm)
                knots=usgs._segments(x,y,nppm,params['segment_means'])
                curve=usgs._spline(*knots,nm)
                raw=[None if v is None else float(usgs.F(v)*.1*scale) for v in values]
                quiet=[float(v*.1*scale) for v in curve]
                residual=[None if v is None else float(.1*(usgs.F(v)-q)*scale) for v,q in zip(values,curve)]
                intervals=[]
                for i in range(8):
                    section=[v for v in residual[i*180:(i+1)*180] if v is not None]
                    rounded=[usgs._nint(usgs.F(v)*.1*scale) for v in values[i*180:(i+1)*180] if v is not None]
                    residual_range=float(usgs.F(max(section))-usgs.F(min(section)));raw_range=max(rounded)-min(rounded)
                    amplitude=rows[i]['x_range_nT' if c==0 else 'y_range_nT']
                    if amplitude is None: reason='Unavailable: insufficient coverage'
                    elif abs(amplitude-raw_range)<1e-5:
                        reason=('Residual amplitude capped at raw range' if residual_range>raw_range+1e-5
                                else 'Raw and residual ranges agree; no amplitude reduction' if abs(residual_range-raw_range)<1e-5
                                else 'Raw range selected by confidence tests')
                    elif amplitude<residual_range-1e-5: reason='Significant trend adjustment reduced the residual amplitude'
                    else: reason='Residual amplitude retained by confidence tests'
                    intervals.append(dict(rawRange=raw_range,residualRange=residual_range,amplitude=amplitude,reason=reason,
                        k=rows[i]['k_x' if c==0 else 'k_y'],decimal=rows[i]['k_x_decimal' if c==0 else 'k_y_decimal']))
                # Validate the chart's final amplitudes with the numerical kernel and public API.
                actual=usgs._amplitudes(values,curve,scale,params['source_confidence'][c],params['adjustment_confidence'][c])
                assert actual==[r['x_range_nT' if c==0 else 'y_range_nT'] for r in rows]
                results.append(dict(raw=raw,quiet=quiet,residual=residual,blocks=blocks,quietDay=quiet_day,
                    supportX=[float(v) for v in x],supportY=[float(v*.1*scale) for v in y],
                    knotX=[float(v) for v in knots[0]],knotY=[float(v*.1*scale) for v in knots[1]],intervals=intervals))
            config['days'][day.isoformat()]=dict(components=results,k=[r['k'] for r in rows],decimal=[r['k_decimal'] for r in rows])
        data['configurations'][str(nppm)]=config
    return data


def export_figure(data):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,
        'axes.edgecolor':'#d6dde6','axes.labelcolor':'#4a596b','xtick.color':'#657588','ytick.color':'#657588',
        'text.color':'#243348','axes.titleweight':'bold','svg.fonttype':'none'})
    fig,axes=plt.subplots(5,1,figsize=(13,14.5),gridspec_kw={'height_ratios':[.7,1.1,1,.95,.65]},layout='constrained')
    fig.set_facecolor('#f8fafc')
    fig.suptitle('How USGS turns selected quiet means into geomagnetic K indices',fontsize=17,fontweight='bold',x=.06,ha='left')
    config=data['configurations']['60'];result=config['days']['2023-12-01'];trace=result['components'][0];hist=config['histograms'][0]
    t=[(i+1)/60 for i in range(1440)]
    for ax in axes:
        ax.set_facecolor('white');ax.grid(axis='y',color='#edf1f5');ax.set_axisbelow(True)
    axes[0].set_title('1  Calibrate quiet-mean limits from December’s data',loc='left',pad=12)
    axes[0].bar(hist['edges'][:-1],hist['counts'],width=hist['edges'][1]-hist['edges'][0],align='edge',color='#839bb8',edgecolor='white')
    axes[0].axvline(hist['low'],color=COLORS['accepted'],lw=2,label='Mean-selection limits')
    axes[0].axvline(hist['high'],color=COLORS['accepted'],lw=2)
    axes[0].set_xlabel('H block mean relative to calibration baseline (nT)');axes[0].set_ylabel('Count');axes[0].legend(frameon=False,fontsize=9)
    titles=['2  Accept quiet hourly means and reject disturbed blocks','3  Fit line segments, then a cubic spline through their intersections',
            '4  Subtract the spline; confidence tests select or adjust the amplitude','5  Convert final H/D amplitudes to digital K']
    for ax,title in zip(axes[1:],titles):
        ax.set_title(title,loc='left',pad=12);ax.set_xlim(0,24);ax.set_xticks(range(0,25,3))
        for hour in range(3,24,3):ax.axvline(hour,color='#e6ebf1',lw=.8,zorder=0)
    axes[1].plot(t,trace['raw'],color=COLORS['raw'],lw=.8,label='Centered, quantized H minute means')
    for keep,color,marker,label in [(True,COLORS['accepted'],'o','Accepted block mean'),(False,COLORS['rejected'],'x','Rejected block mean')]:
        blocks=[b for b in trace['blocks'] if b['accepted']==keep and b['mean'] is not None]
        axes[1].scatter([b['time'] for b in blocks],[b['mean'] for b in blocks],c=color,marker=marker,s=40,label=label,zorder=4)
    axes[1].set_ylabel('H departure (nT)');axes[1].legend(loc='lower left',frameon=False,fontsize=9)
    axes[2].plot(t,trace['raw'],color=COLORS['raw'],lw=.5,alpha=.35)
    axes[2].plot(trace['knotX'],trace['knotY'],color=COLORS['segment'],linestyle='--',marker='D',ms=4,lw=1.2,label='Piecewise-linear fit / intersections')
    axes[2].plot(t,trace['quiet'],color=COLORS['quiet'],lw=2,label='Cubic-spline regular curve')
    axes[2].set_ylabel('H departure (nT)');axes[2].legend(frameon=False,fontsize=9)
    axes[3].plot(t,trace['residual'],color=COLORS['residual'],lw=.9)
    axes[3].axhline(0,color='#657588',lw=.7);axes[3].set_ylabel('H residual (nT)')
    sample=trace['intervals'][2]
    axes[3].text(.02,.06,f"06–09 UT: raw range {sample['rawRange']:.1f} nT · residual range {sample['residualRange']:.1f} nT · final amplitude {sample['amplitude']:.1f} nT\n{sample['reason']}",
                 transform=axes[3].transAxes,fontsize=9,bbox={'facecolor':'white','edgecolor':'none','alpha':.9})
    centers=[1.5+3*i for i in range(8)]
    axes[4].bar(centers,result['k'],width=2.3,color='#405673')
    for i,(x,k,decimal) in enumerate(zip(centers,result['k'],result['decimal'])):
        axes[4].text(x,k+.15,f'{decimal:.1f}',ha='center',fontweight='bold')
        axes[4].text(x,.25,f"H {trace['intervals'][i]['k']} / D {result['components'][1]['intervals'][i]['k']}",ha='center',va='bottom',fontsize=8,color='white' if k else '#657588')
    axes[4].set_ylim(0,10);axes[4].set_yticks([0,3,6,9]);axes[4].set_ylabel('Integer K');axes[4].set_xlabel('Time (UT)')
    axes[4].text(.01,.89,'Labels show decimal digital K; bars show truncated integer K.',transform=axes[4].transAxes,fontsize=9,color='#657588')
    fig.supxlabel('Kiruna · 1 December 2023 · H shown above; combined K uses H and D · 60-minute means / 3 means per segment\nMonthly bootstrap calibration; H reference estimated from December. These K values require observatory validation.',fontsize=10,color='#657588')
    fig.savefig(OUTPUT/'usgs.png',dpi=160,facecolor=fig.get_facecolor())
    fig.savefig(OUTPUT/'usgs.svg',facecolor=fig.get_facecolor());plt.close(fig)


if __name__=='__main__':
    data=build_data();export_figure(data)
    (OUTPUT/'usgs.html').write_text((OUTPUT/'usgs-template.html').read_text().replace('__DATA__',json.dumps(data,separators=(',',':'))))
    print('Created usgs.html, usgs.png, and usgs.svg')
