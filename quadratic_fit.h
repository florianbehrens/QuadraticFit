/*
 * Copyright (c) 2014 Florian Behrens
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef QUADRATIC_FIT_H
#define QUADRATIC_FIT_H

#include <array>
#include <vector>
#include <cmath>
#include <cstddef>

/** @brief Computes a least squares quadratic fit curve.
 *
 * Calculate the three coefficients in y = a * x^2 + b * x + c.
 * See @a http://mathforum.org/library/drmath/view/72047.html for more
 * information about the algorithm. */
template<typename T = double>
class quadratic_fit
{
public:
    struct point {
        point(T _x = 0.0, T _y = 0.0) : x(_x), y(_y) {}
        T x, y;
    };

    /// Default constructor.
    quadratic_fit()
    {}

    /// @brief Constructor with initial memory reservation.
    /// @param n Number of points to initially reserve memory for.
    explicit quadratic_fit(std::size_t n)
    {
        samples_.reserve(n);
    }

    /// @brief Returns a reference to the point at specified location.
    /// If idx is not within the valid range, an exception of type
    /// std::out_of_range is thrown.
    point& operator[](std::size_t idx) {
        return samples_.at(idx);
    }

    /// @brief Returns a const reference to the point at specified location.
    /// If idx is not within the valid range, an exception of type
    /// std::out_of_range is thrown.
    point const& operator[](std::size_t idx) const {
        return samples_.at(idx);
    }

    /// Add a new point to the algorithm.
    void add(T x, T y) {
        samples_.push_back(point(x, y));
    }

    /// Clear all points of the algorithm.
    void clear() {
        samples_.clear();
    }

    /// Compute the three coefficients.
    std::array<T, 3> compute() const
    {
        std::array<T, 3> retval;

        // Compute coefficient a
        retval[0] = (sj1(0) * sj0(1) * sj0(3)
                   - sj1(1) * sj0(0) * sj0(3)
                   - sj1(0) * std::pow(sj0(2), 2)
                   + sj1(1) * sj0(1) * sj0(2)
                   + sj1(2) * sj0(0) * sj0(2)
                   - sj1(2) * std::pow(sj0(1), 2))
                   / (sj0(0) * sj0(2) * sj0(4)
                   - std::pow(sj0(1), 2) * sj0(4)
                   - sj0(0) * std::pow(sj0(3), 2)
                   + 2 * sj0(1) * sj0(2) * sj0(3)
                   - std::pow(sj0(2), 3));

        // Compute coefficient b
        retval[1] = (sj1(1) * sj0(0) * sj0(4)
                   - sj1(0) * sj0(1) * sj0(4)
                   + sj1(0) * sj0(2) * sj0(3)
                   - sj1(2) * sj0(0) * sj0(3)
                   - sj1(1) * std::pow(sj0(2), 2)
                   + sj1(2) * sj0(1) * sj0(2) )
                   / (sj0(0) * sj0(2) * sj0(4)
                   - std::pow(sj0(1), 2) * sj0(4)
                   - sj0(0) * std::pow(sj0(3), 2)
                   + 2 * sj0(1) * sj0(2) * sj0(3)
                   - std::pow(sj0(2), 3));

        // Compute coefficient c
        retval[2] = (sj1(0) * sj0(2) * sj0(4)
                   - sj1(1) * sj0(1) * sj0(4)
                   - sj1(0) * std::pow(sj0(3), 2)
                   + sj1(1) * sj0(2) * sj0(3)
                   + sj1(2) * sj0(1) * sj0(3)
                   - sj1(2) * std::pow(sj0(2), 2))
                   / (sj0(0) * sj0(2) * sj0(4)
                   - std::pow(sj0(1), 2) * sj0(4)
                   - sj0(0) * std::pow(sj0(3), 2)
                   + 2 * sj0(1) * sj0(2) * sj0(3)
                   - std::pow(sj0(2), 3));

        return retval;
    }

private:
    /// @brief Compute sj0 sum.
    /// Computes the sum sj0 = sum(i = 0, n, samples_[i].x^j).
    inline T sj0(std::size_t j) const {
        T retval = 0.0;
        for (auto &val : samples_)
            retval += std::pow(val.x, j);
        return retval;
    }

    /// @brief Compute sj1 sum.
    /// Computes the sum sj1 = sum(i = 0, n, samples_[i].x^j * samples_[i].y).
    inline T sj1(std::size_t j) const {
        T retval = 0.0;
        for (auto &val : samples_)
            retval += (std::pow(val.x, j) * val.y);
        return retval;
    }

    std::vector<point> samples_;
};

#endif // QUADRATIC_FIT_H
