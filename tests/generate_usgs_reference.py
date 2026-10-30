"""Regenerate standalone USGS Fortran fixtures; GNU Fortran needed only here."""
from pathlib import Path
import subprocess, os, shlex, tempfile, json, hashlib, sys
ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT))
import calculate_usgs_k_index as usgs
import calculate_k_index as ih


def build_reference(p):
    files=['KMEAN','MOVDAT','XYSET','GS','SPLN2','NTRP0','CSPLN2','GKR','MPR','RND','RKV','CURVK']
    paths=[]
    for n in files:
        s=(ROOT / 'vendor/usgs' / (n+'.FOR')).read_text().replace('\x1a','').expandtabs(8)
        s=s.replace('IABS(','ABS(').replace('IMOD(','MOD(')
        s=s.replace('FLOATI(','FLOAT(').replace('FLOATJ(','FLOAT(').replace('ININT(','NINT(')
        s=s.replace('REAL FUNCTION CSPLN2*4','REAL*4 FUNCTION CSPLN2').replace('INTEGER FUNCTION MPR*2','INTEGER*2 FUNCTION MPR').replace('INTEGER FUNCTION RND*2','INTEGER*2 FUNCTION RND')
        if n=='XYSET': s=s.replace('INTEGER*2 H,E,DV,I,J,NM,NPPM','INTEGER*2 H,E,DV,I,J,NM,NPPM,NPTS')
        if n=='NTRP0': s=s.replace('IF(I.GT.NR)GOTO 15','IF(I.GE.NR)GOTO 15')
        if n=='GKR': s=s.replace("',DEL='F5.1)","',DEL=',F5.1)")
        f=p/(n.lower()+'.f');f.write_text(s);paths.append(str(f))
    command=shlex.split(os.environ.get('FC','gfortran'))
    result=subprocess.run(command+['-O0','-std=legacy','-ffixed-line-length-132',
        '-ffp-contract=off','-fno-align-commons','-fallow-argument-mismatch',
        *paths,str(ROOT/'vendor/usgs/reference_driver.f90'),'-o',str(p/'reference')],capture_output=True,text=True)
    if result.returncode: raise RuntimeError(result.stderr)
    return p/'reference'


def build_histogram(p):
    paths=[]
    for name in ('MSTAT','SORT0','RND'):
        text=(ROOT/'vendor/usgs'/(name+'.FOR')).read_text().replace('\x1a','').expandtabs(8)
        text=text.replace('IABS(','ABS(').replace('IMOD(','MOD(').replace('INTEGER FUNCTION RND*2','INTEGER*2 FUNCTION RND')
        if name=='MSTAT': text=text.replace('IF(IX(I).GE.IT)GOTO 8','IF(IX(I).GE.IT.AND.J.LT.NTV)GOTO 8')
        path=p/(name.lower()+'.f');path.write_text(text);paths.append(str(path))
    command=shlex.split(os.environ.get('FC','gfortran'))
    result=subprocess.run(command+['-O0','-std=legacy','-ffixed-line-length-132','-ffp-contract=off',
        '-fno-align-commons',*paths,str(ROOT/'vendor/usgs/histogram_driver.f90'),'-o',str(p/'histogram')],capture_output=True,text=True)
    if result.returncode: raise RuntimeError(result.stderr)
    return p/'histogram'


def reference_data(days, params, target, k9):
    from datetime import timedelta
    encoded=usgs._encode(usgs._polar(days),params['base_h_nT'],params['base_d_arcmin'])
    neighbors=[target+timedelta(days=d) for d in (-1,0,1)]
    data=f"{params['nppm']} {params['segment_means']-1} {params['xd']} {k9}\n"
    data+=' '.join('1' if d.isoformat() in params['quiet_days'] else '0' for d in neighbors)+'\n'
    for c in (0,1):
        limits=params['components'][c]
        data+=' '.join(str(limits[k]) for k in ('high','low','end','range'))+'\n'
        data+=' '.join(map(str,params['source_confidence'][c]))+'\n'
        data+=' '.join(map(str,params['adjustment_confidence'][c]))+'\n'
    for d in neighbors:
        for values in encoded.get(d,[[None]*1440]*2):
            data+='\n'.join('32767' if v is None else str(v) for v in values)+'\n'
    return data


def main():
    cases=[]
    histograms=[]
    with tempfile.TemporaryDirectory(prefix='usgs-reference-') as directory:
        executable=build_reference(Path(directory))
        histogram=build_histogram(Path(directory))
        for values in ([i*i-300 for i in range(40)], list(range(-10,31)), [3]*20+[4]*30+[8]*10+[13]*4):
            data=str(len(values))+'\n'+'\n'.join(map(str,values))+'\n'
            result=subprocess.run([str(histogram)],input=data,text=True,capture_output=True,check=True)
            histograms.append(dict(values=values,expected=list(map(float,result.stdout.split()))))
        for station in ('KIR','BEL','WIC','HYB'):
            path=next(p for p in (ROOT/station).iterdir() if p.name.lower()==station.lower()+'23dec.bin')
            records=ih.parse_month_file(path);days=usgs.iaf_days(records);dates=sorted(days)
            configs=[(1,60,3),(10,60,3),(21,60,3),(0,60,3),(30,60,3)]
            if station=='WIC': configs += [(10,30,4),(10,12,3)]
            for index,nppm,ncp in configs:
                params=usgs.bootstrap_parameters(days,nppm=nppm)
                params['segment_means']=ncp
                target=dates[index];k9=days[target][2][10]
                result=subprocess.run([str(executable)],input=reference_data(days,params,target,k9),
                    capture_output=True,text=True,check=True)
                expected=[list(map(float,line.split())) for line in result.stdout.splitlines() if len(line.split())==2]
                if len(expected)!=16: raise RuntimeError(result.stdout+result.stderr)
                cases.append(dict(path=str(path.relative_to(ROOT)),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                    day_index=index,parameters=params,expected=expected))
    hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT/'vendor/usgs').glob('*.FOR')}
    (ROOT/'tests/data/usgs-reference.json').write_text(json.dumps(dict(source_hashes=hashes,cases=cases,histograms=histograms),indent=2)+'\n')


if __name__=='__main__': main()
