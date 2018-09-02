#include <mex.h>
#include <cuda.h>
#include <cuda_runtime.h>

#define cudaCheckError()                                                                     \
    {                                                                                        \
        cudaError_t e = cudaGetLastError();                                                  \
        if (e != cudaSuccess) {                                                              \
            printf("Cuda failure %s:%d: '%s'\n", __FILE__, __LINE__, cudaGetErrorString(e)); \
        }                                                                                    \
    }

__global__ void cudaLaunch(float* imData, float* output, int nSize, int searchSize, mwSize dx, mwSize dy, mwSize dz)
{
    int x = threadIdx.x + blockIdx.x * blockDim.x,
        y = threadIdx.y + blockIdx.y * blockDim.y,
        z = threadIdx.z + blockIdx.z * blockDim.z,
        coord[3] = { x + nSize, y + nSize, z + nSize },
        imDims[3] = { dx, dy, dz };

    if (x >= (imDims[0] - 2 * nSize) || y >= (imDims[1] - 2 * nSize) || z >= (imDims[2] - 2 * nSize) || 
        (int)imData[coord[2] * imDims[0] * imDims[1] + coord[1] * imDims[0] + coord[0]] == 0) {
        
        return;
    }

    int minc[3], maxc[3];
    long int offset = 25 * (x + y * (dx - 2 * nSize) + z * (dx - 2 * nSize) * (dy - 2 * nSize));
    nSize >>= 1;

    int i;
    for (i = 0; i < 3; ++i) {
        minc[i] = coord[i] + 1 - searchSize / 2;
        maxc[i] = minc[i] + searchSize;
        if (minc[i] < nSize)
            minc[i] = nSize;
        
        if (maxc[i] > imDims[i] - nSize)
            maxc[i] = imDims[i] - nSize;
    }

    float diff, val, dist[25];
    // Find 25 closest neighbors. Initialize initial distances to INF.
    for (i = 0; i < 25; i++)
        dist[i] = 1E+37; 

    for (z = minc[2]; z < maxc[2]; ++z) {
        for (y = minc[1]; y < maxc[1]; ++y) {
            for (x = minc[0]; x < maxc[0]; ++x) {
                if (x == coord[0] && y == coord[1] && z == coord[2])
                    continue;

                long int bi = z * imDims[0] * imDims[1] + y * imDims[0] + x;
                long int bj = coord[2] * imDims[0] * imDims[1] + coord[1] * imDims[0] + coord[0];

                diff = 0;
                long int rel1 = -nSize * imDims[0] * imDims[1] - nSize * imDims[0] - nSize;
                for (dz = -nSize; dz <= nSize; ++dz) {
                    long int rel2 = rel1;
                    for (dy = -nSize; dy <= nSize; ++dy) {
                        long int rel3 = rel2;
                        for (dx = -nSize; dx <= nSize; ++dx) {
                            val = imData[bi + rel3] - imData[bj + rel3];
                            if (dz || dy || dx)
                                diff += (val * val);
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

void mexFunction(int nlhs, mxArray* plhs[],
    int nrhs, const mxArray* prhs[])
{
    const mwSize* imDims = mxGetDimensions(prhs[0]);
    float* imData = (float*)mxGetData(prhs[0]);
    double* searchSize = mxGetPr(prhs[1]);
    double* neighSize = mxGetPr(prhs[2]);

    int dims[2] = { 25 * (imDims[2] - 2 * neighSize[2]) * (imDims[1] - 2 * neighSize[1]) * (imDims[0] - 2 * neighSize[0]), 1 };
    plhs[0] = mxCreateNumericArray(2, dims, mxSINGLE_CLASS, mxREAL);
    float* output = (float*)mxGetData(plhs[0]);

    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    float *d_imData, *d_output, executionTime;

    cudaEventRecord(start, 0);
    cudaMalloc((void**)&d_output, dims[0] * sizeof(float));
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&executionTime, start, stop);
    mexPrintf("\n1. Execution time: gpu mem allocation for output gpuArray %f", executionTime);

    //	EXECUTION TIME-> MEMORY ALLOCATION ON GPU + TRANSFER OF DATA FROM RAM TO GPU MEMORY.
    cudaEventRecord(start, 0);
    cudaMalloc((void**)&d_imData, imDims[0] * imDims[1] * imDims[2] * sizeof(float));
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&executionTime, start, stop);
    mexPrintf("\n2. Execution time: gpu mem allocation for imageData on GPU memory %f", executionTime);

    cudaEventRecord(start, 0);
    cudaMemcpy(d_imData, imData, imDims[0] * imDims[1] * imDims[2] * sizeof(float), cudaMemcpyHostToDevice);
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&executionTime, start, stop);
    mexPrintf("\n3. Execution time: ImageData transfer from RAM to GPU memory %f", executionTime);
    // EXECUTION TIME-> MEMORY ALLOCATION ON GPU + TRANSFER OF DATA FROM RAM TO GPU MEMORY.

    dim3 grid((imDims[0] - 2 * neighSize[0] + 7) / 8, (imDims[1] - 2 * neighSize[1] + 7) / 8, (imDims[2] - 2 * neighSize[2] + 3) / 4),
        block(8, 8, 4);
    cudaLaunch<<<grid, block> > >(d_imData, d_output, (int)neighSize[0], (int)searchSize[0], imDims[0], imDims[1], imDims[2]);
    cudaDeviceSynchronize();

    //	EXECUTION TIME-> COMPUTED RESULTS FROM GPU MEMORY TO RAM
    cudaEventRecord(start, 0);
    cudaMemcpy(output, d_output, dims[0] * sizeof(float), cudaMemcpyDeviceToHost);
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&executionTime, start, stop);
    mexPrintf("\n4. Execution time: gpuArray transfer from GPU memory to RAM %f\n", executionTime);
    //	EXECUTION TIME-> COMPUTED RESULTS FROM GPU MEMORY TO RAM

    cudaCheckError();
    cudaFree(d_output);
    cudaFree(d_imData);
}
