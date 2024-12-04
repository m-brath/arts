#include "workspace_method_class.h"

#include <auto_wsm.h>

#include <algorithm>
#include <exception>
#include <iomanip>
#include <ranges>
#include <stdexcept>
#include <string_view>

#include "auto_wsv.h"
#include "callback.h"
#include "format_tags.h"
#include "workspace_agenda_class.h"
#include "workspace_class.h"

std::ostream& operator<<(std::ostream& os, const Method& m) {
  if (m.setval) {
    const Wsv& wsv = m.setval.value();
    if (wsv.holds<CallbackOperator>()) {
      const CallbackOperator& cb = wsv.get_unsafe<CallbackOperator>();
      for (auto& var : cb.outputs) {
        os << var << ", ";
      }
      os << " := " << m.name << '(';
      for (auto& var : cb.inputs) {
        os << var << ", ";
      }
      os << ')';
    } else {
      os << m.name;
      std::string var = std::visit(
          [](auto& v) { return std::format("{}", *v); }, wsv.value());
      constexpr std::size_t maxsize = 50;
      if (var.size() > maxsize) {
        var = std::string(var.begin(), var.begin() + maxsize) + "...";
      }
      std::replace(var.begin(), var.end(), '\n', ' ');

      os << " = " << var;
    }
  } else {
    os << m.name;
    os << '(';

    os << "o : ";
    bool first = true;
    for (auto& o : m.outargs) {
      if (not first) os << ", ";
      first = false;
      os << o;
    }

    os << '\n' << std::string(m.name.size() + 1, ' ') << "i : ";
    first = true;
    for (auto& o : m.inargs) {
      if (not first) os << ", ";
      first = false;
      os << o;
    }
    os << ')';
  }

  return os;
}

Method::Method() : name("this-is-not-a-method") {}

Method::Method(const std::string& n,
               const std::vector<std::string>& a,
               const std::unordered_map<std::string, std::string>& kw) try
    : name(n), outargs(workspace_methods().at(name).out), inargs(workspace_methods().at(name).in) {
  const std::size_t nargout = outargs.size();
  const std::size_t nargin  = inargs.size();

  // FIXME: IN C++23, USE ZIP HERE INSTEAD AS WE CAN REMOVE LATER CODE DOING THAT
  std::vector<std::pair<std::string, bool>> outargs_set(nargout);
  std::vector<std::pair<std::string, bool>> inargs_set(nargin);
  for (std::size_t i = 0; i < nargout; ++i)
    outargs_set[i] = {outargs[i], false};
  for (std::size_t i = 0; i < nargin; ++i) inargs_set[i] = {inargs[i], false};

  // Common filter
  const auto unset =
      std::views::filter([](const auto& p) { return not p.second; });

  // Common G-name
  const auto is_gname = [](const auto& str1, auto& str2) {
    return str1.front() == internal_prefix and
           std::string_view(str1.begin() + 1, str1.end()) == str2;
  };

  // Positional arguments
  {
    const auto fuzzy_equals = [is_gname](auto& arg) {
      return std::views::filter([&arg, is_gname](const auto& p) {
        return p.first == arg or is_gname(p.first, arg);
      });
    };

    std::transform(a.begin(),
                   a.begin() + std::min(nargout, a.size()),
                   outargs.begin(),
                   outargs_set.begin(),
                   [&](auto& arg, auto& orig) {
                     if (auto x = inargs_set | unset | fuzzy_equals(orig);
                         bool(x)) {
                       auto ptr    = x.begin();
                       ptr->first  = arg;
                       ptr->second = true;
                     }
                     return std::pair<std::string, bool>{arg, true};
                   });

    auto unset_filt = inargs_set | unset;
    std::transform(
        a.begin() + std::min(nargout, a.size()),
        a.end(),
        unset_filt.begin(),
        [](auto& arg) { return std::pair<std::string, bool>{arg, true}; });
  }

  // Named arguments
  {
    for (auto& [key, val] : kw) {
      bool any = false;

      for (auto& arg : outargs_set | unset) {
        if (arg.first == key) {
          arg.first  = val;
          arg.second = true;
          any        = true;
        }

        if (is_gname(arg.first, key)) {
          arg.first  = val;
          arg.second = true;
          any        = true;
        }
      }

      for (auto& arg : inargs_set | unset) {
        if (arg.first == key) {
          arg.first  = val;
          arg.second = true;
          any        = true;
        }

        if (is_gname(arg.first, key)) {
          arg.first  = val;
          arg.second = true;
          any        = true;
        }
      }

      if (not any) {
        throw std::runtime_error(std::format("No named argument \"{}\"", key));
      }
    }
  }

  // FIXME: REMOVE THESE TWO IN C++23 WITH ZIP
  for (std::size_t i = 0; i < nargout; ++i) {
    if (outargs_set[i].second) {
      outargs[i] = outargs_set[i].first;
    }
  }
  for (std::size_t i = 0; i < nargin; ++i) {
    if (inargs_set[i].second) {
      inargs[i] = inargs_set[i].first;
    }
  }

  // Check that all non-defaulted GINS are set
  for (std::size_t i = 0; i < nargin; i++) {
    if (inargs[i].front() == internal_prefix and
        not workspace_methods().at(n).defs.contains(inargs[i])) {
      throw std::runtime_error(std::format(
          "Missing required generic input argument \"{}\"",
          std::string_view(inargs[i].begin() + 1, inargs[i].end())));
    }
  }
} catch (std::out_of_range&) {
  throw std::runtime_error(std::format("No method named \"{}\"", n));
} catch (std::exception& e) {
  throw std::runtime_error(
      std::format("Error in method construction for \"{}\"\n{}",
                  n,
                  std::string_view(e.what())));
}

Method::Method(std::string n, const Wsv& wsv, bool overwrite)
    : name(std::move(n)), setval(wsv), overwrite_setval(overwrite) {
  if (not setval) {
    throw std::runtime_error(std::format(
        "Cannot set workspace variable \"{}\" to empty value", name));
  }

  if (wsv.holds<CallbackOperator>()) {
    const auto& cb = wsv.get_unsafe<CallbackOperator>();
    outargs        = cb.outputs;
    inargs         = cb.inputs;
  } else {
    outargs = {name};
  }
}

void Method::operator()(Workspace& ws) const try {
  if (setval) {
    if (const Wsv& wsv = setval.value(); wsv.holds<CallbackOperator>()) {
      wsv.get_unsafe<CallbackOperator>()(ws);
    } else {
      if (overwrite_setval) {
        ws.overwrite(name, wsv.copy());
      } else {
        ws.set(name, wsv.copy());
      }
    }
  } else {
    workspace_methods().at(name).func(ws, outargs, inargs);
  }
} catch (std::out_of_range&) {
  throw std::runtime_error(std::format("No method named \"{}\"", name));
} catch (std::exception& e) {
  throw std::runtime_error(
      std::format("Error in method {}\n{}", *this, std::string_view(e.what())));
}

void Method::add_defaults_to_agenda(Agenda& agenda) const {
  if (not setval) {
    const auto& map = workspace_methods().at(name).defs;
    for (auto& arg : inargs) {
      if (arg.front() == internal_prefix and map.contains(arg)) {
        agenda.add(Method{arg, map.at(arg), true});
      }
    }
  }
}

Method::Method(const std::string& n,
               const std::vector<std::string>& ins,
               const std::vector<std::string>& outs,
               const std::optional<Wsv>& wsv,
               bool overwrite)
    : name(n),
      outargs(outs),
      inargs(ins),
      setval(wsv),
      overwrite_setval(overwrite) {}

std::string std::formatter<Wsv>::to_string(const Wsv& wsv) const {
  return std::visit(
      []<typename T>(const std::shared_ptr<T>& val) {
        return std::format("{}", *val);
      },
      wsv.value());
}

std::string wsv_format(const std::string& x) {
  std::string_view n = x;

  while (n.size() > 1 and
         (n.front() == internal_prefix or n.front() == named_input_prefix))
    n.remove_prefix(1);

  std::string_view sep =
      workspace_variables().contains(std::string{n}) ? "*" : "";

  return std::format("{}{}{}", sep, n, sep);
}

struct SetvalHelper {
  std::string method_name;
  std::string wsv_name;
};

template <>
struct std::formatter<std::vector<SetvalHelper>> {
  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return ctx.begin();
  }

  template <class FmtContext>
  FmtContext::iterator format(const std::vector<SetvalHelper>& vec,
                              FmtContext& ctx) const {
    if (vec.empty()) return ctx.out();

    std::format_to(ctx.out(), ", using");
    std::string_view s = ": ";
    for (auto& v : vec) {
      std::format_to(ctx.out(),
                     "{}{} = *{}*",
                     std::exchange(s, ", "sv),
                     wsv_format(v.method_name),
                     v.wsv_name);
    }

    return ctx.out();
  }
};

std::string Method::sphinx_list_item() const {
  if (setval.has_value()) {
    return std::format(
        "{} = {}",
        wsv_format(std::string{name}),
        std::visit(
            []<typename T>(const std::shared_ptr<T>& x) -> std::string {
              if constexpr (std::same_as<T, std::string>) {
                return std::format("\"{}\"", *x);
              } else if constexpr (std::same_as<Numeric, T> or
                                   std::same_as<Size, T> or
                                   std::same_as<Index, T>) {
                return std::format("{}", *x);
              } else {
                return std::format("{:sqNB,}", *x);
              }
            },
            setval->value()));
  }

  std::vector<SetvalHelper> setvals;
  const WorkspaceMethodRecord& wsm = workspace_methods().at(name);
  for (Size i = 0; i < outargs.size(); i++) {
    if (outargs[i] != wsm.out[i] and outargs[i].front() != named_input_prefix) {
      setvals.push_back({wsm.out[i], outargs[i]});
    }
  }

  for (Size i = 0; i < inargs.size(); i++) {
    if (inargs[i] != wsm.in[i] and inargs[i].front() != named_input_prefix and
        wsm.out.end() == std::ranges::find(wsm.out, wsm.in[i])) {
      setvals.push_back({wsm.in[i], inargs[i]});
    }
  }

  return std::format("*{}*{}", name, setvals);
}
