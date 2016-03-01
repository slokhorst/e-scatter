/*
 * Copyright 2014-2016 Thomas Verduin <T.Verduin@tudelft.nl>
 *
 * This program is part of the electron-matter interaction program (SCATTER).
 *
 * SCATTER is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 */

/**
 * @file src/cuda_safe_call.cuh
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#include "cuda_safe_call.cuh"
#include <sstream>
#include <stdexcept>

__host__ void cuda_safe_call(const char* file, int line, std::function<void(void)> callback) {
	try {
		callback();
	} catch(...) {
		cuda_safe_call(file, line, [&](){});
		throw;
	}
	if(cudaPeekAtLastError() != cudaSuccess) {
		std::ostringstream oss;
		oss << cudaGetErrorString(cudaGetLastError());
		oss << " in '" << file << "' at line " << line;
		throw std::runtime_error(oss.str());
	}
}
