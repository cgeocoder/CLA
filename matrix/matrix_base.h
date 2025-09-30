#pragma once
#include <exception>

#ifndef __MATRIX_BASE_H__
#define __MATRIX_BASE_H__

/*
	Cpp Linear Algebra Library
*/
namespace cppml {

	/*
		Matrix template base class
	*/
	template <typename T, unsigned int _Rows, unsigned int _Columns>
	class Matrix {
	private:
		void alloc_data() {
			m_Data = new T*[_Rows];
			for (unsigned int i = 0; i < _Rows; ++i) {
					m_Data[i] = new T[_Columns];
			}
		}

	public:
		Matrix() {
			static_assert((_Rows >= 1) && (_Columns >= 1));

			alloc_data();

			for (unsigned int i = 0; i < _Rows; ++i) {
				for (unsigned int j = 0; j < _Columns; ++j) {
					m_Data[i][j] = (T)0;
				}
			}
		}

		Matrix(const Matrix<T, _Rows, _Columns>& _Copy) {
			static_assert((_Rows >= 1) && (_Columns >= 1));

			alloc_data();

			for (unsigned int i = 0; i < _Rows; ++i) {
				std::memcpy(m_Data[i], _Copy.m_Data[i], sizeof(T) * _Columns);
			}
		}

		template<size_t rl, size_t cl>
		Matrix(const T(&_Mat)[rl][cl]) {
			static_assert((_Rows >= 1) && (_Columns >= 1) && (rl == _Rows) && (cl == _Columns));

			alloc_data();

			for (unsigned int i = 0; i < _Rows; ++i) {
				std::memcpy(m_Data[i], _Mat[i], sizeof(T) * _Columns);
			}
		}

		template<size_t rl, size_t cl>
		Matrix(T (&_Mat)[rl][cl]) {
			static_assert((_Rows >= 1) && (_Columns >= 1) && (rl == _Rows) && (cl == _Columns));

			alloc_data();

			for (unsigned int i = 0; i < _Rows; ++i) {
				std::memcpy(m_Data[i], _Mat[i], sizeof(T) * _Columns);
			}
		}

		~Matrix() {
			for (unsigned int i = 0; i < _Rows; ++i) {
				delete[] m_Data[i];
			}

			delete[] m_Data;
		}

		inline bool is_square() const { return _Rows == _Columns; }

		inline T** data() const { return m_Data; }

		inline T* operator[](unsigned int _Index) const {
			if (_Index >= _Rows)
				throw std::exception(__FUNCTION__ "(): _Index out of range (_Rows)");

			return m_Data[_Index];
		}

		inline T& at(unsigned int _Row, unsigned int _Column) {
			if (_Row >= _Rows)
				throw std::exception(__FUNCTION__ "(): _Row out of range (_Rows)");

			if (_Column >= _Columns)
				throw std::exception(__FUNCTION__ "(): _Column out of range (_Columns)");

			return m_Data[_Row][_Column];
		}

		Matrix<T, _Columns, _Rows> transpose(void) const {
			Matrix<T, _Columns, _Rows> res;

			for (unsigned int i = 0; i < _Rows; ++i) {
				for (unsigned int j = 0; j < _Columns; ++j) {
					res.m_Data[j][i] = m_Data[i][j];
				}
			}

			return res;
		}

		template<typename U>
		Matrix<U, _Rows, _Columns> cast_to(void) const {
			Matrix<U, _Rows, _Columns> res;

			for (unsigned int i = 0; i < _Rows; ++i) {
				for (unsigned int j = 0; j < _Columns; ++j) {
					res.m_Data[i][j] = (U)m_Data[i][j];
				}
			}

			return res;
		}

		Matrix<T, _Rows, _Columns>& operator = (
			const Matrix<T, _Rows, _Columns>& _Right
		) {
			for (unsigned int i = 0; i < _Rows; ++i) {
				std::memcpy(m_Data[i], _Right.m_Data[i], sizeof(T) * _Columns);
			}

			return *this;
		}

		bool operator == (
			Matrix<T, _Rows, _Columns> _Right
		) const {
			for (unsigned int i = 0; i < _Rows; ++i) {
				for (unsigned int j = 0; j < _Columns; ++j) {
					if (m_Data[i][j] != (T)_Right.m_Data[i][j])
						return false;
				}
			}

			return true;
		}

		bool operator != (
			Matrix<T, _Rows, _Columns> _Right
		) const {
			for (unsigned int i = 0; i < _Rows; ++i) {
				for (unsigned int j = 0; j < _Columns; ++j) {
					if (m_Data[i][j] == (T)_Right.m_Data[i][j])
						return false;
				}
			}

			return true;
		}

		Matrix<T, _Rows, _Columns>& operator = (
			Matrix<T, _Rows, _Columns>& _Right
		) {
			for (unsigned int i = 0; i < _Rows; ++i) {
				std::memcpy(m_Data[i], _Right.m_Data[i], sizeof(T) * _Columns);
			}

			return *this;
		}

		template<typename U>
		Matrix<T, _Rows, _Columns>& operator += (
			Matrix<U, _Rows, _Columns>& _Right
			) {
			for (unsigned int i = 0; i < _Rows; ++i) {
				for (unsigned int j = 0; j < _Columns; ++j) {
					m_Data[i][j] += (T)_Right.m_Data[i][j];
				}
			}

			return *this;
		}

		template<typename U>
		Matrix<T, _Rows, _Columns>& operator -= (
			Matrix<U, _Rows, _Columns>& _Right
			) {
			for (unsigned int i = 0; i < _Rows; ++i) {
				for (unsigned int j = 0; j < _Columns; ++j) {
					m_Data[i][j] -= (T)_Right.m_Data[i][j];
				}
			}

			return *this;
		}

		template<typename U>
		Matrix<T, _Rows, _Columns>& operator *= (
			const U& _Right
			) {
			for (unsigned int i = 0; i < _Rows; ++i) {
				for (unsigned int j = 0; j < _Columns; ++j) {
					m_Data[i][j] *= (T)_Right;
				}
			}

			return *this;
		}

		template<typename U>
		Matrix<T, _Rows, _Columns>& operator /= (
			const U& _Right
			) {
			for (unsigned int i = 0; i < _Rows; ++i) {
				for (unsigned int j = 0; j < _Columns; ++j) {
					m_Data[i][j] /= (T)_Right;
				}
			}

			return *this;
		}

	public:
		T **m_Data;
	};

	template<
		typename T, unsigned int _RowsT, unsigned int _ColumnsT,
		typename U, unsigned int _RowsU, unsigned int _ColumnsU>
	auto operator * (
		Matrix<T, _RowsT, _ColumnsT>
		_Left, Matrix<U, _RowsU, _ColumnsU> _Right
	) {
		static_assert((_ColumnsT == _RowsU) && "product of the matrices is not defined");

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

	// 2x2
	template<typename T>
	inline T det(Matrix<T, 2, 2> _Mat) {
		return _Mat.at(0, 0) * _Mat.at(1, 1) - _Mat.at(0, 1) * _Mat.at(1, 0);
	}

	template<typename T, int n_src, int n_res>
	T matrix_high_det(Matrix<T, n_src, n_src> _Mat) {
		T res{};

		for (int row = 0; row < n_src; ++row) {
			Matrix<T, n_res, n_res> tmp_mat;

			for (int i = 1, tline = 0, trow = 0; i < n_src; ++i) {
				for (int j = 0; j < n_src; ++j) {
					if (j != row)
						tmp_mat.at(tline, trow++) = _Mat.at(i, j);
				}

				tline += 1;
				trow = 0;
			}

			res += _Mat.at(0, row) * (T)std::pow(-1, 2 + row) * det<T>(tmp_mat);
		}

		return res;
	}

	// 3x3
	template<typename T>
	inline T det(Matrix<T, 3, 3> _Mat) { return matrix_high_det<T, 3, 2>(_Mat); }

	// 4x4
	template<typename T>
	inline T det(Matrix<T, 4, 4> _Mat) { return matrix_high_det<T, 4, 3>(_Mat); }

	// 5x5
	template<typename T>
	inline T det(Matrix<T, 5, 5> _Mat) { return matrix_high_det<T, 5, 4>(_Mat); }

	// 6x6
	template<typename T>
	inline T det(Matrix<T, 6, 6> _Mat) { return matrix_high_det<T, 6, 5>(_Mat); }

	// 7x7
	template<typename T>
	inline T det(Matrix<T, 7, 7> _Mat) { return matrix_high_det<T, 7, 6>(_Mat); }

	// 8x8
	template<typename T>
	inline T det(Matrix<T, 8, 8> _Mat) { return matrix_high_det<T, 8, 7>(_Mat); }

	// 9x9
	template<typename T>
	inline T det(Matrix<T, 9, 9> _Mat) { return matrix_high_det<T, 9, 8>(_Mat); }

	// 10x10
	template<typename T>
	inline T det(Matrix<T, 10, 10> _Mat) { return matrix_high_det<T, 10, 9>(_Mat); }
}

#endif // !__MATRIX_BASE_H__