#pragma once

#include <vector>
#include <Eigen/Sparse>

namespace Eigen {
	template <typename T>
	std::vector<Eigen::Triplet<T>> sparseMatrixToTriplets(const SparseMatrix<T> & M) {
		std::vector<Eigen::Triplet<T>> triplets;
		for (int i = 0; i < M.outerSize(); i++) {
			for (typename Eigen::SparseMatrix<T>::InnerIterator it(M, i); it; ++it) {
				triplets.push_back(Eigen::Triplet<T>(it.row(), it.col(), it.value()));
			}
		}
		return triplets;
	}

	template <typename T>
	void sparseMatrixToFile(const SparseMatrix<T> & M, const std::string & filename) {
		std::cout << "Write matrix to file: " << filename << std::endl;
		std::ofstream file(filename);
		for (int i = 0; i < M.outerSize(); i++) {
			for (typename Eigen::SparseMatrix<T>::InnerIterator it(M, i); it; ++it) {
				file << it.row() << "," << it.col() << "," << it.value() << std::endl;
			}
		}
	}
}

