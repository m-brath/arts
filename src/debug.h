/* Copyright (C) 2013
   Oliver Lemke <olemke@core-dump.info>

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA. */

/**
  \file  debug.h

  Helper macros for debugging.

  \author Oliver Lemke
  \date 2013-04-25 */

#ifndef debug_h
#define debug_h

#ifndef NDEBUG
#include <exception>
#include <iostream>
#include <sstream>

// Use this macro around function parameter names and variable definitions
// which are only used in assertions
#define DEBUG_ONLY(...) __VA_ARGS__

// Use this macro to output a counter value everytime a
// certain place is reached
#define DEBUG_COUNTER(n)                                    \
  {                                                         \
    static Index n = 0;                                     \
    std::cerr << "DBG: " << #n << ": " << ++n << std::endl; \
  }

// Print value of expression for debugging
#define DEBUG_PRINT(e) \
  { std::cerr << "DBG: " << (e) << std::endl; }

// Print expression and value for debugging
#define DEBUG_VAR(e) \
  { std::cerr << "DBG: " << #e << ": " << (e) << std::endl; }

// Print expression and value with the given fp precision for debugging
#define DEBUG_VAR_FLT(p, e)                                           \
  {                                                                   \
    std::streamsize old_p = std::cerr.precision();                    \
    std::cerr << "DBG: " << #e << ": " << std::setprecision(p) << (e) \
              << std::endl                                            \
              << std::setprecision(old_p);                            \
  }

/*! Turn off noexcept */
#define ARTS_NOEXCEPT

/*! Take all arguments and turn to string by their operator<<() */
template <typename ... Args>
std::string variadic_to_string(Args ... args [[maybe_unused]]) {
  if constexpr (sizeof...(Args) not_eq 0) {
    std::ostringstream os;
    os << std::boolalpha;
    (os << ... << args);
    return os.str();
  } else {
    return "";
  }
}

/*! Condition should be true to pass */
#define ARTS_ASSERT_FALSE(condition, ...) { \
  if (condition) {                          \
    throw std::runtime_error(               \
      std::string("Failed Assert False: "   \
                  #condition "\n") +        \
      variadic_to_string(__VA_ARGS__));     \
  } }

/*! Condition should be false to pass */
#define ARTS_ASSERT_TRUE(condition, ...) {  \
  if (not (condition)) {                    \
    throw std::runtime_error(               \
      std::string("Failed Assert True: "    \
                  #condition "\n") +        \
      variadic_to_string(__VA_ARGS__));     \
  } }

#else

#define DEBUG_ONLY(...)

#define DEBUG_COUNTER(n)

#define DEBUG_PRINT(e)

#define DEBUG_VAR(e)

#define DEBUG_VAR_FLT(p, e)

/*! Turn on noexcept explicitly */
#define ARTS_NOEXCEPT noexcept

/*! Condition should be true to pass */
#define ARTS_ASSERT_TRUE(condition, ...)

/*! Condition should be false to pass */
#define ARTS_ASSERT_FALSE(condition, ...)

#endif /* NDEBUG */

#endif /* debug_h */
