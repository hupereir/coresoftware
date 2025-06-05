#ifndef PHOOL_ONNXLIB_H
#define PHOOL_ONNXLIB_H

#include <onnxruntime_c_api.h>
#include <onnxruntime_cxx_api.h>
// This is a stub for some ONNX code refactoring

Ort::Session *onnxSession(std::string &modelfile);

std::vector<float> onnxInference(Ort::Session *session, std::vector<float> &input, int N, int Nsamp, int Nreturn);

std::vector<float> onnxInference(Ort::Session *session, std::vector<float> &input, int N, int Nx, int Ny, int Nz, int Nreturn);

#endif
