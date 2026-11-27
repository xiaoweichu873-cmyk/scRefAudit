
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include "sinkhorn_solvers.h"
namespace py = pybind11;
using namespace pybind11::literals;
using Sinkhorn::SinkhornResult;
using Sinkhorn::SinkhornSolverOpts;
using Vec = Eigen::VectorXd;
using Mat = Eigen::MatrixXd;
using RMat = Eigen::Ref<const Mat>;
using RVec = Eigen::Ref<const Vec>;
SinkhornResult run(RMat M,RVec a,RVec b,double eps,RVec x0,double tol,int max_iter,int method,bool splr){
  SinkhornResult res;SinkhornSolverOpts opts;opts.x0=x0;opts.method=method;
  if(splr)Sinkhorn::sinkhorn_splr_internal(res,M,a,b,eps,opts,tol,max_iter,0,std::cout);
  else Sinkhorn::sinkhorn_ssns_internal(res,M,a,b,eps,opts,tol,max_iter,0,std::cout);
  return res;
}
PYBIND11_MODULE(regot_followup,m){
  py::class_<SinkhornResult>(m,"Result")
    .def_readonly("niter",&SinkhornResult::niter)
    .def_readonly("dual",&SinkhornResult::dual)
    .def_readonly("plan",&SinkhornResult::plan)
    .def_readonly("mar_errs",&SinkhornResult::mar_errs);
  m.def("solve",&run,"C"_a,"a"_a,"b"_a,"eps"_a,"x0"_a,"tol"_a=1e-7,"max_iter"_a=10000,"method"_a=1,"splr"_a=false);
  m.def("eigen_threads",[](){return Eigen::nbThreads();});
}
