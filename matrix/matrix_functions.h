#pragma once


#ifndef __MATRIX_FUNCTIONS_H__
#define __MATRIX_FUNCTIONS_H__  

#include "matrix_base.h"

namespace la {

	template<
		typename T, unsigned int _RowsT, unsigned int _ColumnsT,
		typename U, unsigned int _RowsU, unsigned int _ColumnsU>
	auto operator * (
		Matrix<T, _RowsT, _ColumnsT> _Left,
		Matrix<U, _RowsU, _ColumnsU> _Right
	) {
		static_assert((_ColumnsT == _RowsU), 
			"product of the matrices is not defined");

		using res_type = decltype(_Left.m_Data[0][0] * _Right.m_Data[0][0]);
		Matrix<res_type, _RowsT, _ColumnsU> res;

		for (int l_line = 0; l_line < _RowsT; ++l_line) {
			for (int r_col = 0; r_col < _ColumnsU; ++r_col) {
				res_type sum = (res_type)0;

				for (int r_line = 0; r_line < _RowsU; ++r_line) {
					sum += (res_type)_Left.at(l_line, r_line) * (res_type)_Right.at(r_line, r_col);
				}

				res.at(l_line, r_col) = sum;
			}
		}

		return res;
	}

	template<typename T, typename U, unsigned int _Rows, unsigned int _Columns>
	auto operator + (
		Matrix<T, _Rows, _Columns> _Left,
		Matrix<U, _Rows, _Columns> _Right
		) {
		using res_type = decltype(_Left.m_Data[0][0] + _Right.m_Data[0][0]);
		Matrix<res_type, _Rows, _Columns> res;

		for (unsigned int i = 0; i < _Rows; ++i) {
			for (unsigned int j = 0; j < _Columns; ++j) {
				res.at(i, j) = ((res_type)_Left.at(i, j) + (res_type)_Right.at(i, j));
			}
		}

		return res;
	}


	template<typename T, typename U, unsigned int _Rows, unsigned int _Columns>
	auto operator - (
		Matrix<T, _Rows, _Columns> _Left,
		Matrix<U, _Rows, _Columns> _Right
		) {
		using res_type = decltype(_Left.m_Data[0][0] - _Right.m_Data[0][0]);
		Matrix<res_type, _Rows, _Columns> res;

		for (unsigned int i = 0; i < _Rows; ++i) {
			for (unsigned int j = 0; j < _Columns; ++j) {
				res.at(i, j) = ((res_type)_Left.at(i, j) - (res_type)_Right.at(i, j));
			}
		}

		return res;
	}

	template<typename T, typename U, unsigned int _Rows, unsigned int _Columns>
	auto operator * (
		Matrix<T, _Rows, _Columns> _Left,
		U _RightScalar
		) {
		using res_type = decltype(_Left.m_Data[0][0] * _RightScalar);
		Matrix<res_type, _Rows, _Columns> res;

		for (unsigned int i = 0; i < _Rows; ++i) {
			for (unsigned int j = 0; j < _Columns; ++j) {
				res.at(i, j) = ((res_type)_Left.at(i, j) * (res_type)_RightScalar);
			}
		}

		return res;
	}

	template<typename T, typename U, unsigned int _Rows, unsigned int _Columns>
	auto operator * (
		U _LeftScalar,
		Matrix<T, _Rows, _Columns> _Right
		) {
		using res_type = decltype(_Right.m_Data[0][0] * _LeftScalar);
		Matrix<res_type, _Rows, _Columns> res;

		for (unsigned int i = 0; i < _Rows; ++i) {
			for (unsigned int j = 0; j < _Columns; ++j) {
				res.at(i, j) = ((res_type)_Right.at(i, j) * (res_type)_LeftScalar);
			}
		}

		return res;
	}

	template<typename T, typename U, unsigned int _Rows, unsigned int _Columns>
	auto operator / (
		Matrix<T, _Rows, _Columns> _Left,
		U _RightScalar
		) {
		using res_type = decltype(_Left.m_Data[0][0] / _RightScalar);
		Matrix<res_type, _Rows, _Columns> res;

		for (unsigned int i = 0; i < _Rows; ++i) {
			for (unsigned int j = 0; j < _Columns; ++j) {
				res.at(i, j) = ((res_type)_Left.at(i, j) / (res_type)_RightScalar);
			}
		}

		return res;
	}

	// The matrix determinant of matrix 1x1
	template<typename T>
	inline T dt(Matrix<T, 1, 1> _Mat) {
		return _Mat.at(0, 0);
	}

	// The matrix determinant
	// 
	// The method of minor and algebraic addition
	// @param _Mat sqaure matrix n*n
	template<typename T, int _N>
	T dt(Matrix<T, _N, _N> _Mat) {
		static_assert(_N > 1);

		T res{};

		for (int row = 0; row < _N; ++row) {
			Matrix<T, _N - 1, _N - 1> tmp_mat;

			for (int i = 1, tline = 0, trow = 0; i < _N; ++i) {
				for (int j = 0; j < _N; ++j) {
					if (j != row)
						tmp_mat.at(tline, trow++) = _Mat.at(i, j);
				}

				tline += 1;
				trow = 0;
			}

			T sign = (T)(1 - 2 * ((2 + row) % 2));

			res += _Mat.at(0, row) * sign * dt<T>(tmp_mat);
		}

		return res;
	}
}

#endif // !__MATRIX_FUNCTIONS_H__