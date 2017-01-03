#include <mex.h>
#include <cassert>

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    const mwSize* imDims = mxGetDimensions(prhs[0]);
    float* imData = (float*)mxGetData(prhs[0]);
    double* searchSize = mxGetPr(prhs[1]);
    double* neighSize = mxGetPr(prhs[2]);

    int dims[2] = { 25 * (imDims[2] - 2 * neighSize[2]) * (imDims[1] - 2 * neighSize[1]) * (imDims[0] - 2 * neighSize[0]), 1 },
        minc[3], maxc[3], nSize = neighSize[0] / 2, coord[3], ii;
    
    long int offset;

    plhs[0] = mxCreateNumericArray(2, dims, mxSINGLE_CLASS, mxREAL);
    float* output = (float*)mxGetData(plhs[0]);

    for (int k = 0; k < (imDims[2] - 2 * neighSize[2]); ++k) {
        coord[2] = k + neighSize[2];
        for (int j = 0; j < (imDims[1] - 2 * neighSize[1]); ++j) {
            coord[1] = j + neighSize[1];
            for (int i = 0; i < (imDims[0] - 2 * neighSize[0]); ++i) {
                coord[0] = i + neighSize[0];

                if (imData[coord[2] * imDims[0] * imDims[1] + coord[1] * imDims[0] + coord[0]] == 0)
                    continue;
                
                offset = 25 * (coord[2] * imDims[0] * imDims[1] + coord[1] * imDims[0] + coord[0]);

                for (ii = 0; ii < 3; ++ii) {
                    minc[ii] = coord[ii] + 1 - searchSize[ii] / 2;
                    maxc[ii] = minc[ii] + searchSize[ii];
                    
                    if (minc[ii] < nSize)
                        minc[ii] = nSize;
                    
                    if (maxc[ii] > imDims[ii] - nSize)
                        maxc[ii] = imDims[ii] - nSize;
                }

                float diff, val, dist[25];
                
                for (ii = 0; ii < 25; ii++)
                    dist[ii] = 9999999;

                for (int z = minc[2]; z < maxc[2]; ++z) {
                    for (int y = minc[1]; y < maxc[1]; ++y) {
                        for (int x = minc[0]; x < maxc[0]; ++x) {
                            
                            if (x == coord[0] && y == coord[1] && z == coord[2])
                                continue;

                            long int bi = z * imDims[0] * imDims[1] + y * imDims[0] + x;
                            long int bj = coord[2] * imDims[0] * imDims[1] + coord[1] * imDims[0] + coord[0];

                            diff = 0;
                            long int rel1 = -nSize * imDims[0] * imDims[1] - nSize * imDims[0] - nSize;
                            
                            for (int dz = -nSize; dz <= nSize; ++dz) {
                                long int rel2 = rel1;
                                for (int dy = -nSize; dy <= nSize; ++dy) {
                                    long int rel3 = rel2;
                                    for (int dx = -nSize; dx <= nSize; ++dx) {
                            
                                        val = imData[bi + rel3] - imData[bj + rel3];
                                        
                                        if (dz || dy || dx) {
                                            diff += (val * val);
                                        }
                                        
                                        rel3++;
                                    }
                                    rel2 += imDims[0];
                                }
                                rel1 += imDims[0] * imDims[1];
                            }
                            
                            // Minor fix for the crash here.
                            for (ii = 24; ii > 0 && diff < dist[ii]; ii--) {
                                dist[ii] = dist[ii - 1];
                                output[offset + ii] = output[offset + ii - 1];
                            }

                            if (ii == 0) {
                                dist[0] = diff;
                                output[offset] = imData[bi];
                            }
                            else if (ii != 24) {
                                dist[ii + 1] = diff;
                                output[offset + ii + 1] = imData[bi];
                            }
                        }
                    }
                }
            }
        }
    }
}
