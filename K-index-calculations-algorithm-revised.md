# K-index Calculation Algorithm

## Purpose and definition

This document describes the algorithm developed by Dr. Ingemar Haggstrom to compute local K-indices for the Kiruna magnetometer (IAGA code: KIR).

The K-index quantifies disturbances in the horizontal magnetic field components, X and Y, on a scale from 0 to 9. For each three-hour interval, calculate the range of the deflection from the presumed solar-quiet (SQ) value. The range is the maximum deflection minus the minimum deflection.

## Kiruna K-index thresholds

| K-index | Deflection range (nT) |
|---:|---:|
| 0 | 0 to less than 15 |
| 1 | 15 to less than 30 |
| 2 | 30 to less than 60 |
| 3 | 60 to less than 120 |
| 4 | 120 to less than 210 |
| 5 | 210 to less than 360 |
| 6 | 360 to less than 600 |
| 7 | 600 to less than 990 |
| 8 | 990 to 1500, inclusive |
| 9 | greater than 1500 |

The source threshold table does not state how to handle shared boundaries. The intervals above make each shared boundary belong to the higher index, while honoring the source's explicit statement that K9 is greater than 1500 nT. Do not round the calculated range before classification.

## Solar-quiet curve

The SQ curve represents the regular daily variation of the geomagnetic field during low solar activity. It consists of 24 values per component, one for each hour of the day. The curve is updated monthly.

The magnetometer records data every second. Interpolate the 24 hourly SQ values for X, Y, and Z to obtain an SQ value for every second.

### Monthly SQ update

1. Read the previous month's 24 hourly SQ mean values for X, Y, and Z.
2. Interpolate the hourly values to one SQ value per second.
3. Select quiet intervals by comparing observations with the interpolated SQ values at the corresponding times. A time qualifies for the first pass when the difference is within 100 nT (the K3 threshold). Form a solar-quiet interval file from the qualifying times.
4. For each hour and component, take the median of the qualifying data to produce a new 24-value SQ file.
5. Interpolate this new hourly curve to one value per second.
6. Repeat the quiet-interval selection using a tighter difference of 10 nT (the K0 threshold).
7. For each hour and component, take the median of the qualifying data to produce the updated 24-value SQ file.
8. Use this file as the SQ curve for the next month.

The two passes use progressively tighter thresholds because the SQ curve is expected to vary by less than 10 nT from month to month. The first pass allows deviations up to 100 nT; the second refines the result using a 10 nT limit.

## First-run bootstrap when no previous SQ file exists

The monthly update above assumes that a previous month's SQ file is available. On the first run, there is no such file, so first create a provisional seed curve from the observations. The original algorithm does not specify this initialization; the following is a recommended bootstrap procedure.

1. Collect an initial calibration window. A full month is a practical default because the SQ curve is updated monthly. Use the GFZ International Quiet and Disturbed Days list to identify candidate dates: Q1 is the quietest day of the month, through Q10, and D1-D5 are the most disturbed days. Prefer Q1-Q5 when they provide enough valid measurements; otherwise include additional Q days up to Q10. Exclude D1-D5 and data already identified as invalid.
2. Group the measurements by hour of day, using the same time basis and timestamp alignment that will be used for the operational SQ curve.
3. For each of the 24 hours and each component (X, Y, and Z), take the median of the valid measurements in that hour across the calibration window. These 72 medians form the provisional 24-hour SQ seed.
4. Interpolate the seed to one value per second and run the standard monthly SQ update above on the calibration data: first with the 100 nT criterion, then with the 10 nT criterion. Use the resulting curve as the initial operational SQ curve.
5. If an hour lacks enough valid observations to estimate a reliable median, do not fill it with zero. Extend the calibration window or obtain a suitable SQ curve from another trusted source. Until all hourly values are available, mark the SQ curve and any K-indices that depend on it as unavailable or provisional.

The GFZ quiet/disturbed-day classification is based on planetary Kp criteria, not Kiruna-specific measurements. Treat it as a way to choose candidate dates, not as a substitute for the local 100 nT and 10 nT filters in the SQ update. The dataset is available from the [GFZ Kp data page](https://kp.gfz.de/daten) and is described as CC BY 4.0; cite Matzka et al. (2021), [Geomagnetic Kp index, version 1.0](https://doi.org/10.5880/Kp.0001), when using it.

The calibration-window length and minimum number of observations per hour are not specified by the source algorithm and should be documented by the implementation.

## K-index calculation

For each day, calculate K-indices in eight consecutive three-hour intervals:

- 00:00-03:00
- 03:00-06:00
- 06:00-09:00
- 09:00-12:00
- 12:00-15:00
- 15:00-18:00
- 18:00-21:00
- 21:00-24:00

For each interval, subtract the interpolated SQ value from each observed horizontal component (X and Y). Find the maximum and minimum deflection in the interval and subtract the minimum from the maximum. Classify the resulting deflection range using the threshold table above.

The source describes the X and Y calculations but does not explicitly state how their two ranges are combined into a single reported K-index. That combination rule must be specified by the implementation.

\newpage

## Example solar-quiet file

The source document's example shows 24 hourly mean values. The columns are X, Y, and Z in pT.

| Hour | X (pT) | Y (pT) | Z (pT) |
|---:|---:|---:|---:|
| 00 | 10638894 | -91957 | 51817732 |
| 01 | 10641294 | -92162 | 51817152 |
| 02 | 10643882 | -92377 | 51817457 |
| 03 | 10643891 | -92356 | 51817470 |
| 04 | 10641757 | -92354 | 51818139 |
| 05 | 10642735 | -93584 | 51817495 |
| 06 | 10641265 | -95282 | 51817208 |
| 07 | 10641714 | -95960 | 51814396 |
| 08 | 10638799 | -94869 | 51817234 |
| 09 | 10637983 | -95847 | 51817535 |
| 10 | 10636815 | -99162 | 51818190 |
| 11 | 10636732 | -99002 | 51818487 |
| 12 | 10641875 | -93074 | 51812544 |
| 13 | 10641863 | -93029 | 51812523 |
| 14 | 10641310 | -98623 | 51820069 |
| 15 | 10644902 | -97646 | 51816755 |
| 16 | 10640593 | -91926 | 51807756 |
| 17 | 10646601 | -96004 | 51818344 |
| 18 | 10646277 | -95996 | 51818349 |
| 19 | 10646277 | -95632 | 51817617 |
| 20 | 10644188 | -94023 | 51816800 |
| 21 | 10644110 | -90536 | 51810423 |
| 22 | 10637333 | -91737 | 51806640 |
| 23 | 10637400 | -91740 | 51806645 |

## Example solar-quiet interval file

| Start time (YYYYMMDDHH) | End time (YYYYMMDDHH) |
|---:|---:|
| 2016030223 | 2016030300 |
| 2016030317 | 2016030321 |
| 2016030418 | 2016030501 |
| 2016030502 | 2016030504 |
| 2016030600 | 2016030601 |
| 2016030603 | 2016030604 |
| 2016031302 | 2016031305 |
| 2016032217 | 2016032219 |
| 2016032316 | 2016032317 |
| 2016032517 | 2016032518 |
| 2016032521 | 2016032600 |
| 2016032616 | 2016032617 |
| 2016032816 | 2016032817 |
| 2016033122 | 2016033123 |
