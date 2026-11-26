"""Plot USGS monthly H calibration with the actual ten-bin statistics.

Optional matplotlib is required only to export this explanatory figure.
"""
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT))
import calculate_k_index as ih
import calculate_usgs_k_index as usgs


def histogram(values, width=None, low=None, count=10):
    low=min(values) if low is None else low
    width=max(1,(max(values)-low+9)//10) if width is None else width
    counts=[0]*count
    for value in values:
        counts[min(max((value-low)//width,0),count-1)]+=1
    return [low+i*width for i in range(count+1)],counts


def main():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    days=usgs.iaf_days(ih.parse_month_file(ROOT/'KIR/KIR23DEC.BIN'))
    params=usgs.bootstrap_parameters(days)
    encoded=usgs._encode(usgs._polar(days),params['base_h_nT'],params['base_d_arcmin'])
    means,ranges=[],[]
    for values in encoded.values():
        m,r=usgs._blocks(values[0],60)
        means.extend(v for v in m if v is not None)
        ranges.extend(v for v in r if v is not None)
    cutoff=usgs._nint(usgs.F(1500*.14)*10)
    included=[r for r in ranges if r<=cutoff];excluded=[r for r in ranges if r>cutoff]
    mean,spread,_=usgs._histogram(means)
    rmean,rspread,_=usgs._histogram(included)
    limits=params['components'][0]
    low,high,rlimit=(limits[key]/10 for key in ('low','high','range'))
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False,
        'axes.edgecolor':'#d6dde6','text.color':'#243348','axes.labelcolor':'#4a596b',
        'xtick.color':'#657588','ytick.color':'#657588','svg.fonttype':'none'})
    fig,axes=plt.subplots(2,1,figsize=(13,8.8),layout='constrained')
    fig.set_facecolor('#f8fafc')
    fig.suptitle('How USGS chooses monthly quiet-mean limits',fontsize=20,fontweight='bold',x=.06,ha='left')
    for ax in axes:
        ax.set_facecolor('white');ax.set_axisbelow(True);ax.grid(axis='y',color='#edf1f5');ax.set_ylabel('Hourly blocks')
    edges,counts=histogram(means)
    axes[0].set_title('1  Estimate the monthly mean level and spread from 744 block means',loc='left',fontweight='bold',pad=14)
    axes[0].axvspan(low,high,color='#13a88a',alpha=.07)
    axes[0].bar([e/10 for e in edges[:-1]],counts,width=(edges[1]-edges[0])/10,align='edge',color='#839bb8',edgecolor='white')
    axes[0].axvline(low,color='#13a88a',lw=2,label=f'Lower limit: {low:.1f} nT')
    axes[0].axvline(high,color='#13a88a',lw=2,ls='--',label=f'Upper limit: +{high:.1f} nT')
    axes[0].axvline(float(mean)/10,color='#405673',ls=':',lw=1.5,label=f'Histogram mean: {float(mean)/10:.2f} nT')
    axes[0].legend(loc='upper left',frameon=False,fontsize=10)
    axes[0].text(.04,.64,
        f'Mean = {float(mean)/10:.4f} nT\nSpread = {float(spread)/10:.4f} nT\n\nLower = mean − 2.576 × spread = {low:.1f} nT\nUpper = mean + 2.576 × spread = +{high:.1f} nT',
        transform=axes[0].transAxes,va='top',fontsize=10.5,bbox=dict(facecolor='white',alpha=.94,edgecolor='none',pad=10))
    axes[0].set_xlabel('H block mean relative to the calibration baseline (nT)')
    range_edges,range_counts=histogram(included)
    width=range_edges[1]-range_edges[0];range_low=range_edges[0]
    n=max(10,(max(ranges)-range_low)//width+1)
    all_edges,used_counts=histogram(included,width,range_low,n)
    _,excluded_counts=histogram(excluded,width,range_low,n)
    axes[1].set_title('2  Estimate a separate range ceiling after excluding strongly disturbed blocks',loc='left',fontweight='bold',pad=14)
    axes[1].bar([v/10 for v in all_edges[:-1]],used_counts,width=width/10,align='edge',color='#839bb8',edgecolor='white',label=f'{len(included)} ranges used for calibration')
    axes[1].bar([v/10 for v in all_edges[:-1]],excluded_counts,width=width/10,align='edge',color='#d98536',edgecolor='white',label=f'{len(excluded)} ranges excluded (> {cutoff/10:g} nT)')
    axes[1].axvline(rlimit,color='#13a88a',lw=2,label=f'Final range ceiling: {rlimit:.1f} nT')
    axes[1].axvline(cutoff/10,color='#d98536',lw=1.7,ls='--',label='Initial K5-based prefilter')
    axes[1].legend(loc='upper right',frameon=False,fontsize=10)
    axes[1].text(.31,.48,
        f'Mean range = {float(rmean)/10:.4f} nT\nRange spread = {float(rspread)/10:.4f} nT\n\nCeiling = mean range + 1.96 × spread\n             = {rlimit:.1f} nT',
        transform=axes[1].transAxes,va='top',fontsize=11,bbox=dict(facecolor='white',alpha=.94,edgecolor='none',pad=10))
    axes[1].set_xlabel('Within-hour H range: maximum − minimum (nT)')
    fig.supxlabel('Kiruna · December 2023 · 60-minute blocks · H calibration baseline 10503.0 nT · K9 = 1500 nT\nThe reference uses ten-bin histogram estimates and rounds limits to 0.1 nT. Spread is not divided by √N.',fontsize=10,color='#657588')
    output=Path(__file__).resolve().parent
    fig.savefig(output/'usgs-monthly-limits.png',dpi=160,facecolor=fig.get_facecolor())
    fig.savefig(output/'usgs-monthly-limits.svg',facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f'Mean limits: {low}, {high} nT; range ceiling: {rlimit} nT; {len(included)} of {len(ranges)} ranges retained.')


if __name__=='__main__':main()
