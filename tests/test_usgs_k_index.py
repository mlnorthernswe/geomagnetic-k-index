"""USGS standalone Fortran, calibration, input, and integration checks."""
import csv
import hashlib
import io
import json
import math
import subprocess
import sys
import tempfile
import unittest
from datetime import date, timedelta
from pathlib import Path
from unittest.mock import patch

import calculate_k_index as ih
import calculate_usgs_k_index as usgs
import run_pipeline as pipeline

ROOT=Path(__file__).resolve().parent.parent


class TestUSGSReference(unittest.TestCase):
    def test_native_histograms_and_amplitudes(self):
        fixture=json.loads((ROOT/'tests/data/usgs-reference.json').read_text())
        for name,digest in fixture['source_hashes'].items():
            self.assertEqual(hashlib.sha256((ROOT/'vendor/usgs'/name).read_bytes()).hexdigest(),digest)
        for case in fixture['histograms']:
            actual=usgs._histogram(case['values'])
            for a,b in zip(actual,case['expected']): self.assertAlmostEqual(a,b,places=5)
        cache={}
        for case in fixture['cases']:
            with self.subTest(path=case['path'],day=case['day_index'],nppm=case['parameters']['nppm']):
                path=ROOT/case['path']
                self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(),case['sha256'])
                if path not in cache: cache[path]=usgs.iaf_days(ih.parse_month_file(path))
                days=cache[path];target=sorted(days)[case['day_index']]
                rows=usgs.calculate_days(days,[target],parameters=case['parameters'],minimum_coverage=.01)
                for c,field in enumerate(('x_range_nT','y_range_nT')):
                    for i,row in enumerate(rows):
                        want,decimal=case['expected'][c*8+i]
                        if want >= 32767:
                            self.assertIsNone(row[field])
                        else:
                            self.assertAlmostEqual(row[field],want,places=5)
                            key='k_x_decimal' if c==0 else 'k_y_decimal'
                            self.assertEqual(row[key],decimal)

    def test_constant_histogram_and_digital_scale(self):
        self.assertEqual(usgs._histogram([7]*20),(7,0,7))
        self.assertEqual(usgs.digital_k(0,500),0)
        self.assertEqual(usgs.digital_k(5,500),1)
        self.assertEqual(usgs.digital_k(500,500),9)
        self.assertIsNone(usgs.digital_k(6001,500))
        self.assertEqual(usgs._nint(-.5),-1)
        self.assertEqual(usgs._rnd(2.5),2)
        self.assertEqual(usgs._rnd(3.5),4)


class TestUSGSIntegration(unittest.TestCase):
    @staticmethod
    def constant_days():
        day=date(2023,1,2)
        return day,{day:([100000]*1440,[20000]*1440,None,None)}

    def test_polar_conversion_unwrap_and_missing_pairs(self):
        day,days=self.constant_days()
        values=usgs._polar(days)[day]
        self.assertAlmostEqual(values[0][0],math.hypot(10000,2000))
        self.assertAlmostEqual(values[1][0],math.degrees(math.atan2(2000,10000))*60)
        days[day][0][0]=None
        h,d=usgs._polar(days)[day]
        self.assertIsNone(h[0]);self.assertIsNone(d[0])
        encoded=usgs._encode({day:([10000]*1440,[-10799]*1440)},10000,10799)
        self.assertEqual(encoded[day][1][0],20)  # two arcminutes across the branch cut

    def test_no_compiler_boundary_days_and_coverage(self):
        day,days=self.constant_days()
        with patch('subprocess.run',side_effect=AssertionError('compiler invocation')):
            rows=usgs.calculate_days(days,k9_nT=500)
        self.assertEqual([r['k'] for r in rows],[0]*8)
        days[day][0][:60]=[None]*60
        rows=usgs.calculate_days(days,k9_nT=500,minimum_coverage=.01)
        self.assertIsNone(rows[0]['k'])  # native GKR needs at least 121 valid minutes
        self.assertEqual(rows[1]['k'],0)
        self.assertAlmostEqual(rows[0]['x_coverage'],2/3)

    def test_calibration_cli_csv_and_comprehensive(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)/'result.csv';params=Path(tmp)/'parameters.json';comp=Path(tmp)/'comprehensive.txt'
            command=[sys.executable,str(ROOT/'calculate_k_index.py'),str(ROOT/'KIR/KIR23DEC.BIN'),
                '--method','usgs','--usgs-parameters-output',str(params),'-o',str(out),
                '--comprehensive-output',str(comp)]
            result=subprocess.run(command,text=True,capture_output=True)
            self.assertEqual(result.returncode,0,result.stderr)
            self.assertIn('bootstrapped monthly',result.stderr)
            calibration=usgs.read_parameters(params)
            self.assertEqual(calibration['nppm'],60)
            with out.open() as source: rows=list(csv.DictReader(source))
            self.assertEqual(len(rows),248)
            self.assertIn('h_amplitude_nT',rows[0]);self.assertIn('d_amplitude_nT',rows[0])
            self.assertIn('k_decimal',rows[0]);self.assertNotIn('x_range_nT',rows[0])
            self.assertEqual(len([line for line in comp.read_text().splitlines() if line.startswith('2023')]),31)
            original=out.read_bytes()
            result=subprocess.run([sys.executable,str(ROOT/'calculate_k_index.py'),str(ROOT/'KIR/KIR23DEC.BIN'),
                '--method','usgs','--usgs-parameters',str(params),'-o',str(out)],text=True,capture_output=True)
            self.assertEqual(result.returncode,0,result.stderr)
            self.assertEqual(out.read_bytes(),original)

    def test_iaga_minute_cli(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'samples.min';out=Path(tmp)/'result.csv'
            path.write_text('DATE TIME DOY TESTX TESTY\n'+''.join(
                f'2023-01-02 {i//60:02d}:{i%60:02d}:00 002 10000 2000\n' for i in range(1440)))
            result=subprocess.run([sys.executable,str(ROOT/'calculate_k_index.py'),str(path),
                '--method','usgs','--input-format','iaga','--iaga-sample-step','60','--k9-limit','500',
                '-o',str(out)],text=True,capture_output=True)
            self.assertEqual(result.returncode,0,result.stderr)
            with out.open() as source: rows=list(csv.DictReader(source))
            self.assertEqual([r['k'] for r in rows],['0']*8)
            self.assertTrue(all(r['stored_k']=='' for r in rows))

    def test_pipeline_configuration(self):
        with tempfile.TemporaryDirectory() as tmp:
            summary=pipeline.run_station_pipeline('KIR',{'directory':'KIR','method':'usgs'},
                {},ROOT,output_dir=Path(tmp),dry_run=True)
            self.assertEqual(len(summary),4)
            self.assertTrue(all(s.csv_output.name.endswith('-USGS.csv') for s in summary))
            self.assertTrue(all(s.sq_output is None for s in summary))


class TestUSGSValidation(unittest.TestCase):
    def test_invalid_calibration_and_cli(self):
        with self.assertRaisesRegex(ValueError,'version'): usgs.validate_parameters({})
        for args,message in [(['--sq-output','unused.txt'],'own daily quiet curve'),
                             (['--usgs-mean-minutes','13'],'divide 1440'),
                             (['--usgs-h-reference','nan'],'H reference'),
                             (['--longitude','20'],'no longitude parameter')]:
            result=subprocess.run([sys.executable,str(ROOT/'calculate_k_index.py'),
                str(ROOT/'KIR/KIR23DEC.BIN'),'--method','usgs',*args],text=True,capture_output=True)
            self.assertNotEqual(result.returncode,0)
            self.assertIn(message,result.stderr)

    def test_calibration_output_cannot_overwrite_input(self):
        result=subprocess.run([sys.executable,str(ROOT/'calculate_k_index.py'),str(ROOT/'KIR/KIR23DEC.BIN'),
            '--method','usgs','--usgs-parameters-output',str(ROOT/'KIR/KIR23DEC.BIN')],capture_output=True,text=True)
        self.assertNotEqual(result.returncode,0)
        self.assertIn('separate path',result.stderr)
