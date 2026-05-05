// @file comin_keyval.cpp
// @brief Key value storage backend based on std::unordered_map and std::variant
//
// @authors 08/2024 :: ICON Community Interface  <comin@icon-model.org>
//
// SPDX-License-Identifier: BSD-3-Clause
//
// See LICENSES for license information.
// Where software is supplied by third parties, it is indicated in the
// headers of the routines.

#include <cstdint>
#include <iostream>
#include <string>
#include <string_view>
#include <tuple>
#include <unordered_map>
#include <variant>

namespace {
using VariantType = std::variant<int, double, std::string, bool>;
using MapType     = std::unordered_map<std::string, VariantType>;

using VarMapType = std::unordered_map<std::tuple<std::string, int>, void*>;
} // namespace

// Specialization of std::hash for tuples.
template <typename... Ts> struct std::hash<std::tuple<Ts...>> {

  static uint64_t fmix64(uint64_t k) {
    k ^= k >> 33;
    k *= 0xff51afd7ed558ccdULL;
    k ^= k >> 33;
    k *= 0xc4ceb9fe1a85ec53ULL;
    k ^= k >> 33;

    return k;
  }

  static uint64_t rotl(uint64_t x, unsigned n) {
    return (x << n) | (x >> (64 - n));
  }

  static void murmur_round(uint64_t& h, size_t k) {
    constexpr uint64_t c1 = 0x87c37b91114253d5ULL;
    constexpr uint64_t c2 = 0x4cf5ad432745937fULL;

    uint64_t k1 = k;

    k1 *= c1;
    k1 = rotl(k1, 31);
    k1 *= c2;

    h ^= k1;
    h = rotl(h, 27);
    h = h * 5 + 0x52dce729;
  }

  size_t operator()(const std::tuple<Ts...>& tpl) const noexcept {
    return std::apply(
        [](Ts const&... args) {
          size_t h = 0;

          // Do a simplified murmur3 round to combine hashes.
          (..., murmur_round(h, std::hash<Ts>{}(args)));
          return h;
        },
        tpl);
  }
};

extern "C" {
  void comin_keyval_set_int_c(const char* ckey, int val, MapType* map) {
    (*map)[ckey] = val;
  };

  void comin_keyval_get_int_c(const char* ckey, int* val, MapType* map) {
    *val = std::get<int>(map->at(ckey));
  };

  void comin_keyval_set_double_c(const char* ckey, double val, MapType* map) {
    (*map)[ckey] = val;
  };

  void comin_keyval_get_double_c(const char* ckey, double* val, MapType* map) {
    *val = std::get<double>(map->at(ckey));
  };

  void comin_keyval_set_char_c(const char* ckey, char* val, MapType* map) {
    (*map)[ckey] = std::string(val);
  };

  void comin_keyval_get_char_c(const char* ckey, const char** val,
                               MapType* map) {
    *val = std::get<std::string>(map->at(ckey)).data();
  };

  void comin_keyval_set_bool_c(const char* ckey, bool val, MapType* map) {
    (*map)[ckey] = val;
  };

  void comin_keyval_get_bool_c(const char* ckey, bool* val, MapType* map) {
    *val = std::get<bool>(map->at(ckey));
  };

  void comin_keyval_create_c(MapType** map) { *map = new MapType(); };

  void comin_keyval_delete_c(MapType* map) { delete map; };

  void comin_keyval_query_c(const char* ckey, int* idx, MapType* map) {
    auto it = map->find(ckey);
    if (it == map->end())
      *idx = -1;
    else
      *idx = it->second.index();
  };

  void comin_keyval_iterator_begin_c(MapType* map, MapType::iterator** it) {
    *it = new MapType::iterator(map->begin());
  }

  void comin_keyval_iterator_end_c(MapType* map, MapType::iterator** it) {
    *it = new MapType::iterator(map->end());
  }

  const char* comin_keyval_iterator_get_key_c(MapType::iterator* it) {
    return ((*it)->first).c_str();
  }

  bool comin_keyval_iterator_compare_c(MapType::iterator* it1,
                                       MapType::iterator* it2) {
    return (*it1 == *it2);
  }

  void comin_keyval_iterator_next_c(MapType::iterator* it) { (*it)++; }

  void comin_keyval_iterator_delete_c(MapType::iterator* it) { delete it; }

  void comin_varmap_new_c(VarMapType** map) { *map = new VarMapType(); }
  void comin_varmap_delete_c(VarMapType* map) { delete map; }
  void* comin_varmap_get_c(const VarMapType* map, const char* name, size_t len,
                           int id) {
    auto it = map->find({{name, len}, id});
    if (it == map->end())
      return nullptr;

    return it->second;
  }

  void comin_varmap_put_c(VarMapType* map, const char* name, size_t len, int id,
                          void* ptr) {
    map->insert({{{name, len}, id}, ptr});
  }

  void comin_varmap_iterator_begin_c(VarMapType* map,
                                     VarMapType::iterator** it) {
    *it = new VarMapType::iterator(map->begin());
  }
  void comin_varmap_iterator_delete_c(VarMapType::iterator* it) { delete it; }
  void comin_varmap_iterator_next_c(VarMapType::iterator* it) { ++*it; }
  void comin_varmap_iterator_value_c(const VarMapType::iterator* it,
                                     void** ptr) {
    *ptr = (*it)->second;
  }
  bool comin_varmap_iterator_is_end_c(const VarMapType* map,
                                      const VarMapType::iterator* it) {
    return *it == map->end();
  }
}
