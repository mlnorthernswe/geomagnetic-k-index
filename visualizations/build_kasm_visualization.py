"""Build an offline KASM explainer; matplotlib is needed only for PNG/SVG export."""
import json
import math
import sys
from datetime import date, timedelta
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import calculate_k_index as ih
import calculate_kasm_k_index as kasm

OUTPUT = Path(__file__).resolve().parent
EXAMPLES = [('KIR', date(2023, 12, 21), 1500, 'Quiet day · Kiruna · 21 Dec 2023'),
            ('KIR', date(2024, 5, 10), 1500, 'Storm onset · Kiruna · 10 May 2024'),
            ('BEL', date(2024, 5, 13), 450, 'Final K cap · Belsk · 13 May 2024')]


def build_data():
    examples = []
    for station, day, k9, label in EXAMPLES:
        records = []
        for path in sorted((ROOT / station).glob('*.[bB][iI][nN]')):
            records.extend(ih.parse_month_file(path))
        days = kasm.iaf_days(records)
        variants = {}
        for smoothing in (.5, 1., 2.):
            result = {'components': {}}
            public = kasm.calculate_days(days, [day], k9_nT=k9, smoothing=smoothing)
            for ci, component in enumerate(('x', 'y')):
                observed = (days[day-timedelta(days=1)][ci][-180:] + days[day][ci]
                            + days[day+timedelta(days=1)][ci][:180])
                signal = [kasm._UNDEF if v is None else v/10 for v in observed]
                flags = kasm._fill_gaps(signal)
                low = min(signal)
                scale = 500/k9
                normalized = [(v-low)*scale for v in signal]
                first = kasm._solve(normalized, [math.sqrt(5e-7*smoothing)]*len(signal))
                weights = []
                for i in range(0, len(first), 60):
                    q = min(kasm._range(first[i:i+60])/6.66-1, 50)
                    weights.extend([math.sqrt(7e-7*smoothing)*math.exp(-q/2)]*60)
                final = kasm._solve(normalized, weights)
                section = slice(180, 1620)
                center = sum(signal[section])/1440
                raw = [v-center for v in signal[section]]
                initial_residual = [v/scale for v in first[section]]
                residual = [v/scale for v in final[section]]
                ranges = [kasm._range(residual[i:i+180]) for i in range(0,1440,180)]
                raw_ranges = [kasm._range(signal[i:i+180]) for i in range(180,1620,180)]
                component_k = [kasm._classify(v*scale) for v in ranges]
                for i, row in enumerate(public):
                    assert flags[i+1] == 0 and row['k_'+component] == component_k[i]
                    assert math.isclose(ranges[i],row[component+'_range_nT'],abs_tol=1e-8)
                result['components'][component] = dict(
                    raw=raw, initial=[v-r for v,r in zip(raw,initial_residual)],
                    quiet=[v-r for v,r in zip(raw,residual)], residual=residual,
                    initialResidual=initial_residual,
                    # The squared solver weight is the coefficient in the objective.
                    influence=[w*w/(7e-7*smoothing) for w in weights[section]],
                    hourlyRange=[kasm._range(first[i:i+60]) for i in range(180,1620,60)],
                    ranges=ranges, rawRanges=raw_ranges, k=component_k)
            components = result['components']
            filtered = [max(components['x']['k'][i],components['y']['k'][i]) for i in range(8)]
            raw_k = [kasm._classify(max(components[c]['rawRanges'][i] for c in ('x','y'))*500/k9) for i in range(8)]
            combined = [min(a,b) for a,b in zip(filtered,raw_k)]
            assert combined == [row['k'] for row in public]
            result.update(k=combined,filteredK=filtered,rawK=raw_k)
            variants[str(smoothing)] = result
        examples.append(dict(station=station,date=str(day),label=label,k9=k9,variants=variants))
    return examples


def export_figure(data):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'svg.fonttype':'none',
                         'axes.spines.top':False,'axes.spines.right':False})
    example=data[1]; r=example['variants']['1.0']; c=r['components']['x']
    fig,axes=plt.subplots(4,1,figsize=(13,13),layout='constrained',
                          gridspec_kw={'height_ratios':[1.3,1,.9,.8]})
    fig.set_facecolor('#f5f7fa')
    fig.suptitle('Inside KASM: let disturbed hours influence the quiet curve less',
                 fontsize=17,fontweight='bold',x=.06,ha='left')
    t=[i/60 for i in range(1440)]
    titles=['1  Fit once with constant weights, then fit again with adaptive weights',
            '2  First-pass hourly disturbance controls the second fit',
            '3  Subtract the final quiet curve; measure three-hour ranges',
            '4  Combine X/Y classifications and apply the unfiltered K cap']
    for ax,title in zip(axes,titles):
        ax.set_title(title,loc='left',pad=12);ax.set_xlim(0,24);ax.set_xticks(range(0,25,3))
        ax.set_facecolor('white');ax.grid(color='#e7ecf1',linewidth=.6)
    for key,color,label in [('raw','#9aa8ba','Observed X, daily mean removed'),
                            ('initial','#5574d8','First quiet estimate'),('quiet','#079c82','Final quiet estimate')]:
        axes[0].plot(t,c[key],color=color,lw=1.5 if key!='raw' else .6,label=label)
    axes[0].legend(frameon=False,loc='upper left');axes[0].set_ylabel('nT')
    axes[1].bar([i+.5 for i in range(24)],c['hourlyRange'],width=.85,color='#5574d8',label='Hourly range, normalized')
    ax2=axes[1].twinx();ax2.plot([i+.5 for i in range(24)],c['influence'][::60],color='#d08b24',marker='.',label='Relative data influence')
    ax2.set_ylabel('Weight² / (7 × 10⁻⁷ × smoothing)',color='#b37317');axes[1].set_ylabel('Normalized range')
    axes[2].plot(t,c['residual'],color='#c34f69',lw=.8);axes[2].axhline(0,color='#9aa8ba',lw=.5);axes[2].set_ylabel('Residual nT')
    for i in range(8):
        a=c['residual'][i*180:(i+1)*180];axes[2].vlines(i*3+1.5,min(a),max(a),color='#243348',lw=1)
    centers=[i*3+1.5 for i in range(8)]
    axes[3].bar(centers,r['k'],width=2.4,color='#344d6c')
    for i,x in enumerate(centers):axes[3].text(x,r['k'][i]+.15,str(r['k'][i]),ha='center',weight='bold')
    axes[3].set_ylim(0,10);axes[3].set_ylabel('Final K');axes[3].set_xlabel('UTC hour')
    fig.supxlabel('Kiruna · 10 May 2024 · X component shown · K9 = 1500 nT · smoothing = 1\n'
                  '1,440 minute means + 180 minutes from each adjacent day. Final K uses both X and Y.',fontsize=10)
    for suffix in ('png','svg'):fig.savefig(OUTPUT/f'kasm.{suffix}',dpi=160,facecolor=fig.get_facecolor())
    plt.close(fig)


if __name__ == '__main__':
    data=build_data()
    template=(OUTPUT/'kasm-template.html').read_text()
    (OUTPUT/'kasm.html').write_text(template.replace('__DATA__',json.dumps(data,separators=(',',':'),allow_nan=False)))
    export_figure(data)
    print('Created kasm.html, kasm.png and kasm.svg; all displayed ranges and K indices checked against calculate_days.')
