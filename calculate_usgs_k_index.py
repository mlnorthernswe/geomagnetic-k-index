"""Standard-library port of Wilson's USGS mean-selection/spline K method."""
from __future__ import annotations

import json
import math
import struct
import statistics
from bisect import bisect_left
from datetime import date, timedelta
from pathlib import Path

import calculate_k_index as ih
from calculate_kasm_k_index import iaf_days as _iaf_days, iaga_days as _iaga_days

_PACK = struct.Struct('=f')


class F(float):
    """Fortran REAL*4 arithmetic, including rounding of real literals."""
    def __new__(cls, value=0):
        return float.__new__(cls, _PACK.unpack(_PACK.pack(value))[0])
    def __add__(self, other): return F(float(self)+float(F(other)))
    __radd__ = __add__
    def __sub__(self, other): return F(float(self)-float(F(other)))
    def __rsub__(self, other): return F(float(F(other))-float(self))
    def __mul__(self, other): return F(float(self)*float(F(other)))
    __rmul__ = __mul__
    def __truediv__(self, other): return F(float(self)/float(F(other)))
    def __rtruediv__(self, other): return F(float(F(other))/float(self))
    def __neg__(self): return F(-float(self))
    def __abs__(self): return F(abs(float(self)))
    def __pow__(self, power):
        result = F(1)
        for _ in range(power): result = result*self
        return result


def _nint(value):
    return int(math.copysign(math.floor(abs(value)+.5), value))


def _rnd(value):
    """USGS RND (ties to even, including its historical negative-zero quirk)."""
    value = F(value)
    if abs(value) > 32767: return 32767
    integer = int(value)
    magnitude = abs(integer)
    fraction = abs(value-F(integer))*10
    if fraction > 5 or (fraction == 5 and magnitude % 2): magnitude += 1
    return -magnitude if integer < 0 else magnitude


def _sum(values):
    total = F(0)
    for value in values: total = total+value
    return total


def _histogram(values):
    """SELMVS/MSTAT's ten-bin estimates, with defined constant-data handling."""
    values = sorted(v for v in values if v is not None)
    if not values: raise ValueError('USGS calibration has no usable samples')
    low, high = values[0], values[-1]
    if high == low or len(values) == 1: return F(low), F(0), F(low)
    width = (high-low+9)//10
    if width*10 > 32767:
        raise ValueError("USGS histogram exceeds the reference 16-bit range")
    counts = [0]*10
    # The archive overruns the histogram when max is exactly the last upper
    # edge. Include that endpoint in bin 10 rather than reading past it.
    for value in values: counts[min((value-low)//width,9)] += 1
    sg = _sum(counts)
    sfd = _sum(F(i-4)*count for i,count in enumerate(counts))
    sfd2 = _sum(F(i-4)*(F(i-4)*count) for i,count in enumerate(counts))
    hint = F(width-1)*.5
    average = F(low)+hint+4*F(width)+(sfd/sg)*width
    variance = (sg*sfd2-sfd*sfd)/(sg*(sg-1))
    se = F(math.sqrt(max(0,variance)))*width
    modebin = counts.index(max(counts))
    maximum = counts[modebin]
    a = F(maximum-(counts[modebin-1] if modebin else 0))
    b = F(maximum-(counts[modebin+1] if modebin < 9 else 0))
    mode = F(low+modebin*width)+(a/(a+b))*width
    return average,se,mode


def _blocks(values, nppm):
    means, ranges = [], []
    for i in range(0,1440,nppm):
        good = [v for v in values[i:i+nppm] if v is not None]
        if len(good) < nppm//2:
            means.append(None); ranges.append(None)
        else:
            means.append(_nint(_sum(good)/len(good)))
            span = max(good)-min(good)
            ranges.append(span if span < 32767 else None)
    return means,ranges


def iaf_days(records):
    try: return _iaf_days(records)
    except ValueError as error: raise ValueError(str(error).replace('KASM','USGS')) from error


def iaga_days(path, sample_step=1):
    try:
        return _iaga_days(path, sample_step=sample_step)
    except ValueError as error:
        raise ValueError(str(error).replace('KASM', 'USGS')) from error


def _polar(days):
    result = {}
    for day,data in days.items():
        x,y,record,lon=data[:4]
        if len(x) != 1440 or len(y) != 1440:
            raise ValueError('USGS requires 1440 minute samples per component day')
        h,d = [],[]
        for vx,vy in zip(x,y):
            if vx is None or vy is None:
                h.append(None);d.append(None)
            else:
                h.append(math.hypot(vx,vy)/10)
                d.append(math.atan2(vy,vx)*180/math.pi*60)  # arcminutes
        result[day] = (h,d)
    return result


def _encode(polar, base_h, base_d):
    result = {}
    for day,(h,d) in polar.items():
        # Unwrap D around the calibration baseline before quantization.
        arrays = [[None if v is None else _nint((v-base_h)*10) for v in h],
                  [None if v is None else _nint(((v-base_d+10800)%21600-10800)*10) for v in d]]
        if any(v is not None and not -32766 <= v <= 32766 for a in arrays for v in a):
            raise ValueError('USGS centered H/D data exceed the reference 16-bit range; choose a suitable calibration baseline')
        result[day] = arrays
    return result


def bootstrap_parameters(days, target_days=None, *, nppm=60, h_reference=None, k9_nT=None):
    targets = sorted(days if target_days is None else target_days)
    if not targets: raise ValueError('USGS calibration requires input days')
    if nppm < 12 or nppm > 720 or 1440%nppm: raise ValueError('USGS minutes per mean must divide 1440 and be between 12 and 720')
    polar = _polar({d:days[d] for d in targets})
    good_h = [v for day in targets for v in polar[day][0] if v is not None]
    good_d = [v for day in targets for v in polar[day][1] if v is not None]
    if not good_h: raise ValueError('USGS calibration requires valid X/Y pairs')
    base_h = _nint(statistics.median(good_h)*10)/10
    anchor = good_d[0]
    unwrapped = [anchor+(v-anchor+10800)%21600-10800 for v in good_d]
    base_d = _nint(statistics.median(unwrapped)*10)/10
    h_reference = statistics.mean(good_h) if h_reference is None else h_reference
    if not math.isfinite(h_reference) or h_reference <= 0: raise ValueError('USGS H reference must be finite and positive')
    xd = F(F(h_reference)/3438)
    encoded = _encode(polar,base_h,base_d)
    components = []
    day_ranges = []
    for c in (0,1):
        all_means,all_ranges = [],[]
        for day in targets:
            means,ranges = _blocks(encoded[day][c],nppm)
            all_means.extend(means)
            record = days[day][2]
            limit = record[10] if record else None
            # Actual K9 is supplied below by calculate_days for IAGA/overrides.
            k9 = k9_nT if k9_nT is not None else limit
            if k9 is None or not math.isfinite(k9) or k9 < 100:
                raise ValueError('USGS calibration requires a finite K9 limit >= 100 nT')
            cutoff_value = F(k9*.14)*10*F(nppm)/60
            if c: cutoff_value = cutoff_value*34380/F(h_reference)  # SELMVS's published D gate.
            cutoff = _nint(cutoff_value)
            all_ranges.extend(v for v in ranges if v is not None and v <= cutoff)
            if c == 0:
                valid = [32767 if span is None else span for mean,span in zip(means,ranges) if mean is not None]
                if len(valid) >= (1440//nppm)//2:
                    day_ranges.append((float(_sum(valid)/len(valid)),day))
        av,se,mode = _histogram(all_means)
        rav,rse,_ = _histogram(all_ranges)
        components.append(dict(high=_nint(av+F(2.576)*se), low=_nint((av+F(2.576)*se)-2*(F(2.576)*se)),
                               end=_nint(mode), range=_nint(rav+F(1.96)*rse)))
    quiet = [day.isoformat() for _,day in sorted(day_ranges)[:10]]
    return dict(version=1, nppm=nppm, segment_means=3, base_h_nT=base_h,
                base_d_arcmin=base_d, h_reference_nT=h_reference, xd=float(xd),
                components=components,quiet_days=quiet,source_confidence=[[2.58]*8,[2.58]*8],
                adjustment_confidence=[[2.58]*8,[2.58]*8])


def _mean_day(values, params, component, quiet):
    nppm = params['nppm']; limits = params['components'][component]
    means,ranges = _blocks(values,nppm)
    accepted = []
    for block,(mean,span) in enumerate(zip(means,ranges)):
        # KMEAN truncates its mean, unlike SELMVS's rounded calibration means.
        good = [v for v in values[block*nppm:(block+1)*nppm] if v is not None]
        value = int(_sum(good)/len(good)) if mean is not None else None
        if mean is not None and span is None: span = 0  # KMEAN's native wide-range behavior.
        accepted.append(value if value is not None and span is not None and span <= limits['range']
                        and value <= limits['high'] and (quiet or value >= limits['low']) else None)
    good = [v for v in accepted if v is not None]
    daily = int(_sum(good)/len(good)) if len(good) >= len(means)//2 else limits['end']
    return accepted,daily


def _coordinates(triplet, dailies, nppm):
    nm = 1440//nppm; mhv = sum(triplet,[])
    a,b = F(nppm)/60,F(nppm)/120
    x,y = [],[]
    def add(xx,yy): x.append(F(xx));y.append(F(yy))
    if mhv[nm-2] is None or mhv[nm-1] is None:
        old = next((i for i in range(nm-3,-1,-1) if mhv[i] is not None),None)
        current = next((v for v in mhv[nm:2*nm] if v is not None),None)
        add(a*(old+1)-(24+b) if old is not None else -11.5,
            mhv[old] if old is not None else current if current is not None else dailies[0])
    if mhv[nm-2] is not None: add(-3*b,mhv[nm-2])
    if mhv[nm-1] is not None: add(-b,mhv[nm-1])
    for i in range(nm,2*nm):
        if mhv[i] is not None: add(a*(i+1)-(24+b),mhv[i])
    if mhv[2*nm] is not None: add(24+b,mhv[2*nm])
    if mhv[2*nm+1] is not None: add(24+3*b,mhv[2*nm+1])
    if mhv[2*nm] is None or mhv[2*nm+1] is None:
        following = next((i for i in range(2*nm+2,3*nm) if mhv[i] is not None),None)
        current = next((v for v in reversed(mhv[nm:2*nm]) if v is not None),None)
        add(a*(following+1)-(24+b) if following is not None else 24+5*b if current is not None else 35.5,
            mhv[following] if following is not None else current if current is not None else dailies[2])
    if len(x) < 3: return list(map(F,[-11.5,12,35.5])),[F(dailies[2])]*3
    return x,y


def _segments(x0,y0,nppm,ncp):
    """GS's sliding regressions and continuity corrections (1-based indices)."""
    x,y = [F(0)]+x0,[F(0)]+y0
    size = 130
    b0,b1,z,yhat = [[F(0)]*size for _ in range(4)]
    npts = len(x0); nc = min(ncp-1,npts-1)
    k,l,m,n = 0,1,1,1+nc
    guard = 0
    while True:
        guard += 1
        if guard > 10000: raise ValueError('USGS segment fit did not converge')
        sx,sy,sxx,sxy = (F(0) for _ in range(4))
        for j in range(m,n+1):
            sx=sx+x[j];sy=sy+y[j];sxy=sxy+x[j]*y[j];sxx=sxx+x[j]*x[j]
        k += 1; km1=k-1;ce=F(n-m+1)
        b1[k]=(ce*sxy-sx*sy)/(ce*sxx-sx*sx)
        b0[k]=(sy-b1[k]*sx)/ce
        proceed=True
        if k != 1:
            if b1[k] == b1[km1]: k=km1;n+=1
            else:
                l += 1
                z[l]=(b0[km1]-b0[k])/(b1[k]-b1[km1])
                yhat[l]=b0[k]+b1[k]*z[l]
                lost=False
                if z[2] < x[1] or z[l] > x[n]:
                    if l > 2: lost=True
                    else:
                        m -= 1
                        if m > 0: k=km1;l-=1
                        else: k=0;m=2;n=m+nc
                        proceed=False
                elif l > 2 and (z[l] <= z[l-1] or z[l] >= x[n]): lost=True
                if lost:
                    k=km1
                    a=(x[n-1]+x[n])/2-z[l-1]
                    if a <= 0: raise ValueError('USGS segment fit lost continuity')
                    b1[k]=((y[n-1]+y[n])/2-yhat[l-1])/a
                    b0[k]=yhat[l-1]-b1[k]*z[l-1]
                    l-=1
        if not proceed: continue
        m=n-1;n=m+nc
        if n > npts: break
    aa=F(.5)*F(nppm)/60
    if z[2] < 0:
        l-=1
        for i in range(1,l+1): z[i],yhat[i],b1[i],b0[i]=z[i+1],yhat[i+1],b1[i+1],b0[i+1]
    else:
        z[1]=-3*aa;yhat[1]=b0[1]-3*aa*b1[1]
    if z[l] < 24:
        z[l+1]=24+3*aa
        yhat[l+1]=yhat[l]+(24+3*aa-z[l])*b1[l]
        l+=1
    return z[1:l+1],yhat[1:l+1],b1[1:l+1]


def _spline(x,u,slopes,nm):
    """NTRP0 and SPLN2: fill long knot gaps, then clamped cubic spline."""
    x,u,slopes=x.copy(),u.copy(),slopes.copy()
    if len(x) != nm+4:
        i=0
        while i < len(x)-1 and len(x) < nm+4:
            if x[i+1]-x[i] >= 6:
                x.insert(i+1,x[i]+(x[i+1]-x[i])*.5)
                u.insert(i+1,u[i]+(u[i+1]-u[i])*.5)
                slopes.insert(i+1,slopes[i+1])
            else: i+=1
    n=len(x)
    if n < 5: raise ValueError('USGS spline needs at least five knots')
    q1,qn=slopes[0],slopes[n-2]
    delta=[F(0)]+[x[i]-x[i-1] for i in range(1,n)]
    if any(v <= 0 for v in delta[1:]): raise ValueError('USGS spline knots must be increasing')
    if delta[1] == q1: q1=F(0)  # Reference special case.
    v=[F(0)]*n;a=[F(0)]*n;s=[F(0)]*n
    v[0]=6*((u[1]-u[0])/delta[1]-q1)
    for i in range(1,n-1):
        c1=u[i-1]/delta[i]
        c2=u[i]*(1/delta[i]+1/delta[i+1])
        c3=u[i+1]/delta[i+1]
        v[i]=(c1-c2+c3)*6
    v[-1]=(qn+(u[-2]-u[-1])/delta[-1])*6
    a[0]=2*delta[1];a[1]=1.5*delta[1]+2*delta[2];v[1]=v[1]-.5*v[0]
    for i in range(2,n-1):
        c=delta[i]/a[i-1]
        a[i]=2*(delta[i]+delta[i+1])-c*delta[i];v[i]=v[i]-c*v[i-1]
    c=delta[-1]/a[-2];a[-1]=2*delta[-1]-c*delta[-1];v[-1]=v[-1]-c*v[-2]
    s[-1]=v[-1]/a[-1]
    for i in range(n-2,-1,-1): s[i]=(v[i]-delta[i+1]*s[i+1])/a[i]
    curve=[]
    for minute in range(1,1441):
        t=F(minute)/60;k=max(1,bisect_left(x,t));left,right=x[k]-t,t-x[k-1]
        ff1=s[k-1]*(left**3);ff2=s[k]*(right**3);ff3=1/(6*delta[k])
        f1=(ff1+ff2)*ff3
        f2=right*(u[k]/delta[k]-(s[k]*delta[k])/6)
        f3=left*(u[k-1]/delta[k]-(s[k-1]*delta[k])/6)
        curve.append(f1+f2+f3)
    return curve


def _amplitudes(values,curve,xd,sep,sepa):
    """GKR: range/amplitude selection and statistically significant slope adjustment."""
    output=[]
    for interval in range(8):
        start=interval*180
        valid=[(j+1,start+j,values[start+j]) for j in range(180) if values[start+j] is not None]
        ni=len(valid)
        if ni <= 120: output.append(None);continue
        sx,sy,sxx,sxy,syy=(F(0) for _ in range(5))
        residual=[];rounded=[]
        for l,i,value in valid:
            d=.1*(F(value)-curve[i])*xd
            sx=sx+F(l);sxx=sxx+F(l)*F(l);sy=sy+d;syy=syy+d*d;sxy=sxy+d*F(l)
            residual.append((d,i+1))
            rounded.append(_nint(F(value)*.1*xd))
        slope=(F(ni)*sxy-sx*sy)/(F(ni)*sxx-sx*sx)
        v1=syy-sy*sy/F(ni);v2=slope*(sxy-sx*sy/F(ni))
        va=(v1-v2)/F(ni-2)
        if va <= 0: va=F(999)
        vr=syy/F(ni-2)
        # Preserve the first timestamp on ties, as in GKR.
        maximum=max(residual,key=lambda p:p[0]);minimum=min(residual,key=lambda p:p[0])
        amplitude=maximum[0]-minimum[0]
        rawrange=max(rounded)-min(rounded)
        if amplitude > rawrange: amplitude=F(rawrange)
        else:
            jj=maximum[1]//180
            nb=maximum[1]-jj*180;nl=minimum[1]-jj*180
            vxb=(F(nb)-90)*(F(nb)-90);vxl=(F(nl)-90)*(F(nl)-90)
            vx=sxx-sx*sx/F(ni)
            vb=vr/F(ni)+vxb*vr/vx;vl=vr/F(ni)+vr*vxl/vx
            v2=vb+vl;cl1=F(math.sqrt(v2))*F(sep[interval]) if v2 > 0 else F(0)
            if F(rawrange)-amplitude > cl1 and amplitude > cl1:
                if va != 999:
                    adjustment=slope*F(nb)-slope*F(nl)
                    vb=va/F(ni)+va*vxb/vx;vl=va/F(ni)+va*vxl/vx
                    v2=vb+vl;cl2=F(math.sqrt(v2))*F(sepa[interval]) if v2 > 0 else F(0)
                    if abs(adjustment) >= cl2 and abs(adjustment) <= .95*amplitude:
                        amplitude=amplitude-abs(adjustment)
            else: amplitude=F(rawrange)
        output.append(float(amplitude) if 0 <= amplitude < 32767 else None)
    return output


def digital_k(amplitude,k9):
    """KTAB/CURVK/RKV decimal digital K; integer K is its truncated value."""
    magnitude=_rnd(amplitude)
    if magnitude > 6000: return None
    scale=F(k9)/2500
    threshold=_rnd(25*scale)
    if threshold <= 0: raise ValueError('USGS K9 limit is too small for integer-nT digital K')
    value=F(magnitude)
    if magnitude < threshold: vk=value/F(threshold)
    else:
        value=F(math.log10(value/scale))-F(1.39794)
        vk=F(1.003)+F(3.232717)*value+F(.190954)*value*value*value
    vk=F(_rnd(vk*10))*.1
    return round(min(float(vk),float(F(9.9))),1)


def validate_parameters(params):
    if not isinstance(params,dict): raise ValueError('USGS calibration must be a JSON object')
    if params.get('version') != 1: raise ValueError('unsupported USGS calibration version')
    n=params.get('nppm')
    if not isinstance(n,int) or not 12 <= n <= 720 or 1440%n: raise ValueError('USGS nppm must divide 1440 and be between 12 and 720')
    if not isinstance(params.get('segment_means'),int) or not 3 <= params['segment_means'] <= 5: raise ValueError('USGS segment_means must be 3 to 5')
    for key in ('base_h_nT','base_d_arcmin','xd','h_reference_nT'):
        if not isinstance(params.get(key),(int,float)) or not math.isfinite(params[key]): raise ValueError('invalid USGS '+key)
    if params['xd'] <= 0 or params['h_reference_nT'] <= 0: raise ValueError('USGS H reference and XD must be positive')
    if not isinstance(params.get('components'),list) or len(params['components']) != 2: raise ValueError('USGS calibration needs H and D limits')
    for limits in params['components']:
        if not isinstance(limits,dict): raise ValueError('invalid USGS component limits')
        if any(not isinstance(limits.get(k),int) for k in ('high','low','end','range')) or limits['low'] > limits['high'] or limits['range'] < 0:
            raise ValueError('invalid USGS component limits')
    quiet=params.get('quiet_days')
    if not isinstance(quiet,list) or any(not isinstance(v,str) for v in quiet):
        raise ValueError('USGS quiet_days must be a list of ISO dates')
    try:
        for value in quiet: date.fromisoformat(value)
    except ValueError as error:
        raise ValueError('USGS quiet_days must contain ISO dates') from error
    for key in ('source_confidence','adjustment_confidence'):
        rows=params.get(key,[])
        if not isinstance(rows,list) or len(rows) != 2 or any(not isinstance(r,list) or len(r) != 8 or any(not isinstance(v,(int,float)) or not math.isfinite(v) or v <= 0 for v in r) for r in rows):
            raise ValueError('USGS confidence parameters need 2 by 8 positive numbers')


def read_parameters(path):
    params=json.loads(Path(path).read_text())
    validate_parameters(params)
    return params


def calculate_days(days,target_days=None,*,k9_nT=None,minimum_coverage=.9,parameters=None,nppm=60,h_reference=None):
    if not 0 < minimum_coverage <= 1: raise ValueError('minimum coverage must be greater than 0 and at most 1')
    targets=sorted(days if target_days is None else target_days)
    # Each calendar month has its own calibration, as recommended by SELMVS.
    if parameters is None:
        groups={}
        for day in targets: groups.setdefault((day.year,day.month),[]).append(day)
        rows=[]
        for group in groups.values():
            params=bootstrap_parameters(days,group,nppm=nppm,h_reference=h_reference,k9_nT=k9_nT)
            rows.extend(calculate_days(days,group,k9_nT=k9_nT,minimum_coverage=minimum_coverage,parameters=params))
        return rows
    params=read_parameters(parameters) if isinstance(parameters,(str,Path)) else parameters
    validate_parameters(params)
    required={day+timedelta(days=offset) for day in targets for offset in (-1,0,1)}
    selected={day:data for day,data in days.items() if day in required}
    encoded=_encode(_polar(selected),params['base_h_nT'],params['base_d_arcmin'])
    nm=1440//params['nppm'];modes=[c['end'] for c in params['components']]
    means={}
    for day,arrays in encoded.items():
        means[day]=[_mean_day(a,params,c,day.isoformat() in params['quiet_days']) for c,a in enumerate(arrays)]
    rows=[]
    for day in targets:
        record=days[day][2];limit=k9_nT if k9_nT is not None else record[10] if record else None
        if limit is None or not math.isfinite(limit) or limit < 100: raise ValueError('USGS requires K9 >= 100 nT for its rounded integer-nT scale')
        results=[]
        try:
            for c in (0,1):
                neighbors=[means.get(day+timedelta(days=d),[(None,None)]*2)[c] for d in (-1,0,1)]
                triplet=[v if v is not None else [None]*nm for v,_ in neighbors]
                daily=[v if v is not None else modes[c] for _,v in neighbors]
                x,y=_coordinates(triplet,daily,params['nppm'])
                knots=_segments(x,y,params['nppm'],params['segment_means'])
                curve=_spline(*knots,nm)
                results.append(_amplitudes(encoded[day][c],curve,F(params['xd']) if c else F(1),
                    params['source_confidence'][c],params['adjustment_confidence'][c]))
        except (ArithmeticError,OverflowError) as error: raise ValueError('USGS numerical calculation failed') from error
        for i,(start,end) in enumerate(ih.INTERVALS):
            coverage=[sum(v is not None for v in a[i*180:(i+1)*180])/180 for a in encoded[day]]
            amplitudes=[results[c][i] if coverage[c] >= minimum_coverage else None for c in (0,1)]
            decimals=[digital_k(v,limit) if v is not None else None for v in amplitudes]
            indices=[int(v) if v is not None else None for v in decimals]
            k=max(indices) if all(v is not None for v in indices) else None
            decimal=max(decimals) if k is not None else None
            stored_raw=record[ih.K_START_WORD+i] if record else None
            stored=None if stored_raw in (None,ih.MISSING_K_VALUE) else stored_raw/10
            rows.append(dict(date=day.isoformat(),interval_utc=f'{start:02d}:00-{end:02d}:00',
                x_range_nT=amplitudes[0],y_range_nT=amplitudes[1],x_coverage=coverage[0],y_coverage=coverage[1],
                k_x=indices[0],k_y=indices[1],k=k,stored_k_raw=stored_raw,stored_k=stored,
                matches_stored='' if k is None or stored is None else int(k == stored),
                k_decimal=decimal,k_x_decimal=decimals[0],k_y_decimal=decimals[1],components='HD'))
    return rows
