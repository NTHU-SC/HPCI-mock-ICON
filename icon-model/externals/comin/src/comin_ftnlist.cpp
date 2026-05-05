// @file comin_ftnlist.cpp
// @brief List implementation based on std::list
//
// @authors 10/2024 :: ICON Community Interface  <comin@icon-model.org>
//
// SPDX-License-Identifier: BSD-3-Clause
//
// See LICENSES for license information.
// Where software is supplied by third parties, it is indicated in the
// headers of the routines.

#include <iostream>
#include <list>

namespace comin_ftnlist {
using ListType     = std::list<void *>;
using IteratorType = std::list<void *>::iterator;
} // namespace comin_ftnlist

extern "C" {
  void comin_ftnlist_new(comin_ftnlist::ListType **listptr) {
    *listptr = new comin_ftnlist::ListType{};
  }

  void comin_ftnlist_delete(comin_ftnlist::ListType **listptr) {
    delete *listptr;
  }

  void comin_ftnlist_push_back(comin_ftnlist::ListType *listptr, void *ptr) {
    listptr->push_back(ptr);
  }

  void comin_ftnlist_iterator_begin(comin_ftnlist::ListType *listptr,
                                    comin_ftnlist::IteratorType **it_begin) {
    *it_begin = new comin_ftnlist::IteratorType(listptr->begin());
  }

  void comin_ftnlist_iterator_delete(comin_ftnlist::IteratorType **it) {
    delete *it;
  }

  void comin_ftnlist_iterator_next(comin_ftnlist::IteratorType *it) { (*it)++; }

  void comin_ftnlist_iterator_value(comin_ftnlist::IteratorType *it,
                                    void **val) {
    *val = (void *)(**it);
  }

  bool comin_ftnlist_is_end(comin_ftnlist::ListType *listptr,
                            comin_ftnlist::IteratorType *it) {
    return (*it) == listptr->end();
  }
}
