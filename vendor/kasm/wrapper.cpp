// Adapter to the unmodified author's KASM 1.11 smoothing/gap routines.
#include "casm.h"
#include <algorithm>
#include <cmath>
#include <vector>
#include <exception>
extern const Real UNDEF = 1e30;
int fits(int n, Real* x, int* flags, int gap);
int kid(const Real& range);
Real ext(int n, const Real* x);

extern "C" int kasm_details(int n, int offset, double* x, double* y,
                            double k9, double smoothing, int* ks,
                            int* kx, int* ky, double* rx, double* ry,
                            int* invalid_x, int* invalid_y) {
    if ((n != 1440 && n != 1620 && n != 1800) ||
        (offset != 0 && offset != 180) || offset + 1440 > n ||
        !std::isfinite(k9) || k9 <= 0 ||
        !std::isfinite(smoothing) || smoothing <= 0) return 1;
    try {
        int flags_x[10] = {}, flags_y[10] = {};
        fits(n, x, flags_x, 6);
        fits(n, y, flags_y, 6);
        double scale = 500. / k9;
        double xmin = *std::min_element(x, x+n);
        double ymin = *std::min_element(y, y+n);
        for (int i=0; i<n; ++i) {
            x[i] = scale * (x[i] - xmin);
            y[i] = scale * (y[i] - ymin);
        }
        std::vector<Real> dx(n), dy(n);
        Real* raw[2] = {x, y};
        Real* residual[2] = {dx.data(), dy.data()};
        rs::rs1(n, raw, residual, smoothing);
        for (int i=0; i<8; ++i) {
            int start = offset + i*180, flag = offset/180 + i;
            rx[i] = ext(180, dx.data()+start);
            ry[i] = ext(180, dy.data()+start);
            if (!std::isfinite(rx[i]) || !std::isfinite(ry[i])) return 2;
            kx[i] = kid(rx[i]); ky[i] = kid(ry[i]);
            // KASM wwsk caps the final K at the unfiltered combined K.
            ks[i] = std::min(std::max(kx[i], ky[i]),
                         kid(std::max(ext(180, x+start), ext(180, y+start))));
            invalid_x[i] = flags_x[flag]; invalid_y[i] = flags_y[flag];
            rx[i] /= scale; ry[i] /= scale;
        }
    } catch (const std::exception&) { return 3; }
    return 0;
}
