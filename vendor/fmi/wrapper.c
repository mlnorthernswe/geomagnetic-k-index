/* Adapter for FMI 1.3. Keep the reference algorithm and threshold semantics. */
#include "K_index.h"

void fmi_details(long *x, long *y, long step, long limit, long longitude,
                 long missing, long *ks, long *kx, long *ky, long *rx, long *ry)
{
    long i, j, lo, hi;
    InitGlobals(step, limit, longitude, missing);
    FillGaps(x);
    FillGaps(y);
    FillKTable(ks, x, y, 0);
    for (i = 0; i < 2; i++) {
        FindLengths(ks);
        ComputeHourAverages(x, y);
        ComputeHarmonicFit();
        FillKTable(ks, x, y, 1);
    }
    for (i = 0; i < 8; i++) {
        kx[i] = K_MaxMin(X_temp, i);
        ky[i] = K_MaxMin(Y_temp, i);
        lo = hi = X_temp[i * PointsInDay / 8];
        /* Match the reference K_MaxMin sample selection exactly. */
        for (j = i * PointsInDay / 8; j < (i + 1) * PointsInDay / 8 - 1; j++) {
            if (X_temp[j] < lo) lo = X_temp[j];
            if (X_temp[j] > hi) hi = X_temp[j];
        }
        rx[i] = hi - lo;
        lo = hi = Y_temp[i * PointsInDay / 8];
        for (j = i * PointsInDay / 8; j < (i + 1) * PointsInDay / 8 - 1; j++) {
            if (Y_temp[j] < lo) lo = Y_temp[j];
            if (Y_temp[j] > hi) hi = Y_temp[j];
        }
        ry[i] = hi - lo;
    }
    FreeArrays();
}
