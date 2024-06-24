#include <python_interface.h>
#include <rtepack.h>

#include "enums.h"
#include "hpy_numpy.h"
#include "matpack_data.h"
#include "nanobind/nanobind.h"
#include "py_macros.h"

namespace Python {
template <typename T, Index M, size_t... N>
void rtepack_array(py::class_<matpack::matpack_data<T, M>> &c) {
  c.def(
      "__array__",
      [](matpack::matpack_data<T, M> &v) {
        constexpr auto n = M + sizeof...(N);
        std::array<size_t, n> shape{};
        std::ranges::copy(v.shape(), shape.begin());
        std::ranges::copy(std::array{N...}, shape.begin() + M);
        return py::ndarray<py::numpy, Numeric, py::ndim<n>, py::c_contig>(
            v.data_handle(), n, shape.data(), py::handle());
      },
      py::rv_policy::reference_internal);

  c.def_prop_rw("value",
                    [](matpack::matpack_data<T, M> &x) {
                      py::object np = py::module_::import_("numpy");
                      return np.attr("array")(x, py::arg("copy") = false);
                    },
                [](matpack::matpack_data<T, M> &x,
                   matpack::matpack_data<T, M> &y) { x = y; });

  common_ndarray(c);
}

void py_rtepack(py::module_ &m) try {
  py::class_<Stokvec> sv(m, "Stokvec");
  sv.def(py::init_implicit<Numeric>())
      .def("__init__",
           [](Stokvec *s, const PolarizationChoice p) {
             new (s) Stokvec{rtepack::to_stokvec(p)};
           })
      .def("__init__",
           [](Stokvec *s, const String &p) {
             new (s) Stokvec{rtepack::to_stokvec(to<PolarizationChoice>(p))};
           })
      .def_static(
          "linpol",
          [](const Numeric angle) {
            return Stokvec{1.0,
                           Conversion::cosd(2.0 * angle),
                           Conversion::sind(2.0 * angle),
                           0.0};
          },
          "Returns [1.0, cos(2*angle), sin(2*angle), 0.0], the linear polarization vector for a given angle",
          py::arg("angle"))
      .def_static(
          "cirpol",
          [](const Numeric angle) {
            return Stokvec{1.0, 0.0, 0.0, Conversion::sind(angle)};
          },
          "Returns [1.0, 0.0, 0.0, sin(angle)], the circular polarization vector for a given phase delay angle",
          py::arg("angle"))
      .def(
          "__array__",
          [](Stokvec &v) {
            std::array<size_t, 1> shape = {4};
            return py::ndarray<py::numpy, Numeric, py::shape<4>, py::c_contig>(
                v.data.data(), 1, shape.data(), py::handle());
          },
          py::rv_policy::reference_internal)
      .def_prop_rw("value",
                       [](Stokvec &x) {
                         py::object np = py::module_::import_("numpy");
                         return np.attr("array")(x, py::arg("copy") = false);
                       },
                   [](Stokvec &x, Stokvec &y) { x = y; });
  common_ndarray(sv);
  py::implicitly_convertible<PolarizationChoice, Stokvec>();
  py::implicitly_convertible<String, Stokvec>();

  py::class_<StokvecVector> vsv(m, "StokvecVector");
  vsv.def(py::init_implicit<std::vector<Numeric>>())
      .def(py::init_implicit<std::vector<Stokvec>>());
  rtepack_array<Stokvec, 1, 4>(vsv);

  py::class_<StokvecMatrix> msv(m, "StokvecMatrix");
  rtepack_array<Stokvec, 2, 4>(msv);

  py::class_<StokvecTensor3> t3sv(m, "StokvecTensor3");
  rtepack_array<Stokvec, 3, 4>(t3sv);

  py::class_<StokvecTensor4> t4sv(m, "StokvecTensor4");
  rtepack_array<Stokvec, 4, 4>(t4sv);

  py::class_<StokvecTensor5> t5sv(m, "StokvecTensor5");
  rtepack_array<Stokvec, 5, 4>(t5sv);

  py::class_<StokvecTensor6> t6sv(m, "StokvecTensor6");
  rtepack_array<Stokvec, 6, 4>(t6sv);

  py::class_<Propmat> pm(m, "Propmat");
  pm.def(py::init_implicit<Numeric>())
      .def(py::init<Numeric,
                    Numeric,
                    Numeric,
                    Numeric,
                    Numeric,
                    Numeric,
                    Numeric>())
      .def("__array__",
           [](Propmat &x) {
             std::array<size_t, 1> shape = {7};
             return py::ndarray<py::numpy, Numeric, py::shape<7>, py::c_contig>(
                 x.data.data(), 1, shape.data(), py::handle());
           })
      .def_prop_rw("value",
                       [](Propmat &x) {
                         py::object np = py::module_::import_("numpy");
                         return np.attr("array")(x, py::arg("copy") = false);
                       },
                   [](Propmat &x, Propmat &y) { x = y; });
  common_ndarray(pm);

  py::class_<PropmatVector> vpm(m, "PropmatVector");
  vpm.def(py::init_implicit<std::vector<Numeric>>())
      .def(py::init_implicit<std::vector<Propmat>>());
  rtepack_array<Propmat, 1, 7>(vpm);

  py::class_<PropmatMatrix> mpm(m, "PropmatMatrix");
  rtepack_array<Propmat, 2, 7>(mpm);

  py::class_<Muelmat> mm(m, "Muelmat");
  mm.def(py::init_implicit<Numeric>())
      .def(py::init<Numeric,
                    Numeric,
                    Numeric,
                    Numeric,
                    Numeric,
                    Numeric,
                    Numeric,
                    Numeric,
                    Numeric,
                    Numeric,
                    Numeric,
                    Numeric,
                    Numeric,
                    Numeric,
                    Numeric,
                    Numeric>())
      .def("__array__",
           [](Muelmat &x) {
             std::array<size_t, 2> shape = {4, 4};
             return py::
                 ndarray<py::numpy, Numeric, py::shape<4, 4>, py::c_contig>(
                     x.data.data(), 2, shape.data(), py::handle());
           })
      .def_prop_rw("value",
                       [](Muelmat &x) {
                         py::object np = py::module_::import_("numpy");
                         return np.attr("array")(x, py::arg("copy") = false);
                       },
                   [](Muelmat &x, Muelmat &y) { x = y; });
  common_ndarray(mm);

  py::class_<MuelmatVector> vmm(m, "MuelmatVector");
  vmm.def(py::init_implicit<std::vector<Numeric>>())
      .def(py::init_implicit<std::vector<Muelmat>>());
  rtepack_array<Muelmat, 1, 4, 4>(vmm);

  py::class_<MuelmatMatrix> mmm(m, "MuelmatMatrix");
  rtepack_array<Muelmat, 2, 4, 4>(mmm);
} catch (std::exception &e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize rtepack\n", e.what()));
}
}  // namespace Python
