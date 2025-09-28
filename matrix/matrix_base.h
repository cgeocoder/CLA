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
	public:
		Matrix() {
			static_assert((_Rows >= 1) && (_Columns >= 1));

			for (unsigned int i = 0; i < _Rows; ++i) {
				for (unsigned int j = 0; j < _Columns; ++j) {
					m_Data[i][j] = (T)0;
				}
			}
		}

		Matrix(const Matrix<T, _Rows, _Columns>& _Copy) {
			for (unsigned int i = 0; i < _Rows; ++i) {
				std::memcpy(m_Data[i], _Copy.m_Data[i], sizeof(T) * _Columns);
			}
		}

		template<size_t rl, size_t cl>
		Matrix(const T(&_Mat)[rl][cl]) {
			static_assert((_Rows >= 1) && (_Columns >= 1) && (rl == _Rows) && (cl == _Columns));

			for (unsigned int i = 0; i < _Rows; ++i) {
				std::memcpy(m_Data[i], _Mat[i], sizeof(T) * _Columns);
			}
		}

		~Matrix() = default;

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
			Matrix<T, _Rows, _Columns>& _Right
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
			Matrix<T, _Rows, _Columns>& _Right
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
		T m_Data[_Rows][_Columns];
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

}

#endif // !__MATRIX_BASE_H__