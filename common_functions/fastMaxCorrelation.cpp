#include "mex.h"
#include <cmath>
#include <vector>
#include <algorithm>

mwSize idx(mwSize row, mwSize col, mwSize num_rows) {
    return col * num_rows + row;
}

double computeCorrelation(const std::vector<double>& a, const std::vector<double>& b) {
    mwSize n = a.size();
    if (n == 0 || a.size() != b.size()) return 0.0;

    double meanA = 0.0, meanB = 0.0;
    for (mwSize i = 0; i < n; ++i) {
        meanA += a[i];
        meanB += b[i];
    }
    meanA /= n;
    meanB /= n;

    double covariance = 0.0, varianceA = 0.0, varianceB = 0.0;
    for (mwSize i = 0; i < n; ++i) {
        double da = a[i] - meanA;
        double db = b[i] - meanB;
        covariance += da * db;
        varianceA += da * da;
        varianceB += db * db;
    }

    if (varianceA == 0 || varianceB == 0) return 0.0;  // To prevent division by zero

    return covariance / (std::sqrt(varianceA) * std::sqrt(varianceB));
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Check for proper number of arguments.
    if (nrhs != 2) {
        mexErrMsgIdAndTxt("mexFunction:input", "Two input arguments required.");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("mexFunction:output", "One output argument required.");
    }

    double* N = mxGetPr(prhs[0]);
    double* P = mxGetPr(prhs[1]);

    mwSize num_rows = mxGetM(prhs[0]);
    mwSize num_cols_P = mxGetN(prhs[1]);

    std::vector<double> nUpperTri;
    for (mwSize i = 0; i < num_rows; ++i) {
        for (mwSize j = i + 1; j < num_rows; ++j) {
            nUpperTri.push_back(N[idx(i, j, num_rows)]);
        }
    }

    std::vector<double> correlations(num_cols_P);
    for (mwSize col = 0; col < num_cols_P; ++col) {
        std::vector<double> qUpperTri;
        for (mwSize i = 0; i < num_rows; ++i) {
            for (mwSize j = i + 1; j < num_rows; ++j) {
                qUpperTri.push_back(std::abs(P[idx(i, col, num_rows)] - P[idx(j, col, num_rows)]));
            }
        }
        correlations[col] = computeCorrelation(nUpperTri, qUpperTri);
    }

    double maxCorrelation = *std::max_element(correlations.begin(), correlations.end());

    plhs[0] = mxCreateDoubleScalar(maxCorrelation);
}
