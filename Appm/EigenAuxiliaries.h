#pragma once

#include <vector>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <iostream>
#include <iomanip>
#include <fstream>

namespace Eigen {
	/**
	* @return Eigen Library Version as given by the macro file at compile time.
	*/
	inline
	std::string showVersion() {
		std::stringstream ss;
		ss << "Eigen Library Version ";
		ss << EIGEN_WORLD_VERSION << ".";
		ss << EIGEN_MAJOR_VERSION << ".";
		ss << EIGEN_MINOR_VERSION;
		return ss.str();
	}

	/**
	* @return triplets of non-zero entries of a sparse matrix.
	*/
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

	/**
	* Write sparse matrix to file with specified filename.
	*/
	template <typename T>
	void sparseMatrixToFile(const SparseMatrix<T> & M, const std::string & filename) {
		std::cout << "Write matrix to file: " << filename << std::endl;
		std::ofstream file(filename);
		for (int i = 0; i < M.outerSize(); i++) {
			for (typename Eigen::SparseMatrix<T>::InnerIterator it(M, i); it; ++it) {
				file << it.row() << "," << it.col() << "," << std::scientific << std::setprecision(20) << it.value() << std::endl;
			}
		}
		if (M.rows() > 0 && M.cols() > 0) {
			if (M.coeff(M.rows() - 1, M.cols() - 1) != 0) {
				file << M.rows() - 1 << "," << M.cols() - 1 << "," << 0 << std::endl;
			}
		}
	}

	/** 
	* @return if the sparse matrix is symmetric.
	*/
	template <typename T>
	bool isSymmetric(const SparseMatrix<T> & M, const bool showInfo = false) {
		if (!M.isApprox(M.transpose())) {
			if (showInfo) {
				std::cout << "Matrix is not symmetric" << std::endl;
			}
			return false;
		}
		return true;
	}

	/**
	* @return if the sparse matrix is positive definite.
	*/
	template <typename T>
	bool isPositiveDefinite(const SparseMatrix<T> & M, const bool showInfo = false) {
		Eigen::SimplicialLLT<Eigen::SparseMatrix<T>> llt(M);
		if (llt.info() != Eigen::Success) {
			if (showInfo) {
				std::cout << "Matrix is not positive definite (according to Cholesky factorization)" << std::endl;
			}
			return false;
		}
		return true;
	}

	/**
	* @return if the sparse matrix is symmetric positive definite. 
	* 
	* See also: https://stackoverflow.com/a/35230714
	*/
	template <typename T>
	bool isSymmetricPositiveDefinite(const SparseMatrix<T> & M, const bool showInfo = false) {
		return isSymmetric(M, showInfo) && isPositiveDefinite(M, showInfo);
	}



	/** 
	* Get a sparse identitiy matrix of non-square shape. (See also speye() in Matlab)
	*/
	//template <typename T>
	//Eigen::SparseMatrix<double> speye(const int rows, const int cols) {
	//	assert(rows > 0);
	//	assert(cols > 0);
	//	const int n = std::min(rows, cols);

	//	typedef Eigen::Triplet<double> TripletType;
	//	std::vector<TripletType> triplets;
	//	triplets.reserve(n);
	//	for (int i = 0; i < n; i++) {
	//		triplets.emplace_back(TripletType(i, i, 1));
	//	}
	//	Eigen::SparseMatrix<double> M(rows, cols);
	//	M.setFromTriplets(triplets.begin(), triplets.end());
	//	return M;
	//}


	/** 
	* @return Concatenate two sparse matrices vertically.
	*/
	template <typename T>
	Eigen::SparseMatrix<T> vertcat(const Eigen::SparseMatrix<T> & U, const Eigen::SparseMatrix<T> & L) {
		assert(U.cols() == L.cols());
		
		typedef Eigen::Triplet<T> TripletType;
		std::vector<TripletType> triplets;
		triplets.reserve(U.nonZeros() + L.nonZeros());

		// get triplets from upper matrix
		

		for (int k = 0; k < U.outerSize(); ++k) {
			for (typename Eigen::SparseMatrix<T>::InnerIterator it(U, k); it; ++it) {
				triplets.emplace_back(it.row(), it.col(), it.value());
			}
		}
		// get triplets from lower matrix
		for (int k = 0; k < L.outerSize(); ++k) {
			for (typename Eigen::SparseMatrix<T>::InnerIterator it(L, k); it; ++it) {
				triplets.emplace_back(U.rows() + it.row(), it.col(), it.value());
			}
		}
		// set matrix from triplets
		Eigen::SparseMatrix<T> M(U.rows() + L.rows(), U.cols());
		M.setFromTriplets(triplets.begin(), triplets.end());
		return M;
	}
}

