#include <parameters.h>
#include <python_interface.h>

#include <algorithm>
#include <memory>
#include <mutex>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "auto_wsv.h"
#include "workspace_agenda_class.h"
#include "workspace_class.h"

extern Parameters parameters;
namespace Python {
py::tuple pickle_method(const Method& m) {
  return py::make_tuple(
      m.get_name(), m.get_ins(), m.get_outs(), m.get_setval(), m.overwrite());
}

Method unpickle_method(const std::tuple<std::string, std::vector<std::string>, std::vector<std::string>, std::optional<Wsv>, bool>& state) {
  auto name   = std::get<0>(state);
  auto ins    = std::get<1>(state);
  auto outs   = std::get<2>(state);
  auto setval = std::get<3>(state);
  auto ow     = std::get<4>(state);

  return Method{name, ins, outs, setval, ow};
}

std::tuple<std::string, std::vector<Method>, std::vector<std::string>, std::vector<std::string>, bool> pickle_agenda(const Agenda& ag) {
  return std::make_tuple(ag.get_name(),
                        ag.get_methods(),
                        ag.get_share(),
                        ag.get_copy(),
                        ag.is_checked());
}

Agenda unpickle_agenda(const std::tuple<std::string, std::vector<Method>, std::vector<std::string>, std::vector<std::string>, bool>& state) {
  auto n  = std::get<0>(state);
  auto m  = std::get<1>(state);
  auto s  = std::get<2>(state);
  auto c  = std::get<3>(state);
  auto ch = std::get<4>(state);

  return Agenda{n, m, s, c, ch};
}

std::filesystem::path correct_include_path(
    const std::filesystem::path& path_copy);

void py_auto_wsv(py::class_<Workspace>& ws);
void py_auto_wsm(py::class_<Workspace>& ws);

void py_workspace(py::class_<Workspace>& ws) try {
  ws.def("__init__", [](Workspace* w, bool with_defaults) {
           if (with_defaults)
             new (w) Workspace{WorkspaceInitialization::FromGlobalDefaults};
              else new (w) Workspace{WorkspaceInitialization::Empty};
         },
         py::arg("with_defaults") = true)
      .def(py::init<Workspace>())
      .def("__copy__", [](Workspace& w) { return w; })
      .def("__deepcopy__", [](Workspace& w, py::dict&) { return w.deepcopy(); })
      .def(
          "get",
          [](Workspace& w, const std::string& n) { return w.share(n) -> value; },
          py::rv_policy::reference_internal,
          py::keep_alive<0, 1>())
      .def("init", &Workspace::init)
      .def("set",
           [](Workspace& w, const std::string& n, const WsvValue& _x) {
             std::visit([&](auto a){w.set(n, std::make_shared<Wsv>(std::move(a)));}, _x);

             auto& ptr = w.share(n);
             if (ptr->holds<Agenda>()) {
               auto& ag = ptr->get_unsafe<Agenda>();
               if (not ag.is_checked()) {
                 ag.set_name(n);
                 ag.finalize();
               }
             } else if (ptr->holds<ArrayOfAgenda>()) {
               auto& ags = ptr->get_unsafe<ArrayOfAgenda>();
               for (auto& ag : ags) {
                 if (not ag.is_checked()) {
                   ag.set_name(n);
                   ag.finalize();
                 }
               }
             }
           })
      .def("has",
           [](Workspace& w, const std::string& n) { return w.contains(n); })
      .def(
          "swap",
          [](Workspace& w1, Workspace& w2) {
            using std::swap;
            swap(w1, w2);
          },
          py::arg("other"));

  ws.def("__str__", [](const Workspace& w) { return var_string(w); });

  py_auto_wsv(ws);
  py_auto_wsm(ws);
} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize workspace\n", e.what()));
}
}  // namespace Python
