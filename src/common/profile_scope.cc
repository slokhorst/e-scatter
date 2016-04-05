/*!
 * @file src/common/profile_scope.cc
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 */

#include <common/profile_scope.hh>
#include <chrono>

size_t profile_scope(std::function<void()> callback) {
    const auto chrono_start = std::chrono::high_resolution_clock::now();
    callback();
    const auto chrono_stop = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::milliseconds>(chrono_stop-chrono_start).count();
}
