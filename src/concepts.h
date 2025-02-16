/*
 * Created by Jack Tyler on 07/06/2021.
 * Copyright (c) 2021 University of Southampton. All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>. */

#ifndef STRAINLINES_CONCEPTS_H
#define STRAINLINES_CONCEPTS_H

#include <type_traits>

/* @brief Containts concepts (C++20 feature) that enforces certain properties on templated types.
 */

/* @brief Requires that a type supports arithmetic, or is a DACE::DA.*/
template <class T>
concept arithmetic = std::is_arithmetic_v<T> or std::is_same_v<DACE::DA, T>;

/* @brief Requires that a type supports being indexed.*/
template <class T>
concept is_indexable = requires (T A) {
    {A[0]};
};

/* @brief Requires that a type has a size method. */
template <class T>
concept has_size = requires (T A) {
    {A.size()};
};

#endif //STRAINLINES_CONCEPTS_H