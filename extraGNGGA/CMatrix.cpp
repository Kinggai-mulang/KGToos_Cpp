#include"CMatrix.h"
constexpr int MAXMATRIXSIZE = 50; //设置最大矩阵维数


//矩阵乘法
Matrix operator*(const Matrix& a, const Matrix& b) {
	if (a.m_cols != b.m_rows) {
		throw std::invalid_argument("Matrix dimensions do not match for multiplication.");
	}
	Matrix result(a.m_rows, b.m_cols);

	for (size_t i = 0; i < a.m_rows; ++i) {
		for (size_t j = 0; j < b.m_cols; ++j) {
			for (size_t k = 0; k < a.m_cols; ++k) {
				result.m_data[i * b.m_cols + j] += a.m_data[i * a.m_cols + k] * b.m_data[k * b.m_cols + j];
			}
		}
	}
	return result;
}


//矩阵加法
Matrix operator+(const Matrix& a, const Matrix& b) {
	if (a.m_cols != b.m_cols || a.m_rows != b.m_rows) {
		throw std::invalid_argument("Matrix dimensions do not match for add.");
	}
	Matrix result(a.m_rows, a.m_cols);
	for (size_t i = 0; i < a.m_rows * a.m_cols; ++i) {
		result.m_data[i] = a.m_data[i] + b.m_data[i];
	}
	return result;
}

//减法
Matrix operator-(const Matrix& a, const Matrix& b) {
	if (a.m_cols != b.m_cols && a.m_rows != b.m_rows) {
		throw std::invalid_argument("Matrix dimensions do not match for sub.");
	}
	Matrix result(a.m_rows, a.m_cols);
	for (size_t i = 0; i < a.m_rows * a.m_cols; ++i) {
		result.m_data[i] = a.m_data[i] - b.m_data[i];
	}
	return result;
}

/**
 * @file         "CMatrix.cpp"
 * @brief        输出,比较耗时，除了Debug，慎用
 * @param[in]    out
 * @param[in]    A
 * @return       std::ostream &
 * @note         
 * @par History
 *               2024/12/4, Ken Ga  new \n
 */
std::ostream& operator<<(std::ostream& out, const Matrix& A) {
	// 第一步，找到最大宽度
	char precflo = MatrixFixPres;   // 设置浮点小数点后精度

	int maxWidth = 0;   //采用动态宽度
	for (size_t i = 0; i < A.m_rows; ++i) {
		for (size_t j = 0; j < A.m_cols; ++j) {
			// 计算当前元素的字符长度
			int width = (A(i, j) > 10) ? log10(A(i, j)) + precflo+1 : int(log10(A(i, j))) + 2 + precflo;
			if (A(i, j) < 0) { // 负数需要多一个字符用于负号
				width++;
			}
			maxWidth = std::max(maxWidth, width);
		}
	}

	// 第二步，输出矩阵
	for (size_t i = 0; i < A.m_rows; ++i) {
		for (size_t j = 0; j < A.m_cols; ++j) {
			out << std::setw(maxWidth + 2) << std::fixed << std::setprecision(precflo) << A(i, j);
		}
		out << std::endl;
	}
	printf("\n");
	return out;
}

/**
 * @file         "CMatrix.cpp"
 * @brief        
 * @param[in]    a
 * @param[in]    b
 * @return       bool
 * @note         判断是否两个矩阵相等
 * @par History
 *               2024/11/15, Ken Ga  new \n
 */
bool operator==(const Matrix& a, const Matrix& b)
{
	if (a.m_cols==b.m_cols && a.m_rows==b.m_rows)
	{
		for (size_t i = 0; i < a.m_cols * a.m_rows; i++) {
			if (abs(a.m_data[i]-b.m_data[i])>Epsilon_math)
			{
				return false;
			}
		}
		return true;
	}
	else
	{
		return false;
	}
}


const bool Matrix::isNULL()const
{
	if (m_cols == 0 && m_rows == 0 && m_data.empty())return true;
	return false;
}


//矩阵转置
const Matrix Matrix::Tra() const {
	Matrix result(m_cols, m_rows);  // 转置后的矩阵尺寸是原矩阵的 cols x rows

	for (size_t i = 0; i < m_rows; ++i) {
		for (size_t j = 0; j < m_cols; ++j) {
			result.m_data[j * m_rows + i] = m_data[i * m_cols + j];
		}
	}

	return result;
}

/**
 * @file         "CMatrix.cpp"
 * @brief        LU分解
 * @param[in]    L
 * @param[in]    U
 * @return       void
 * @note         
 * @par History
 *               2024/12/4, Ken Ga  new \n
 */
void Matrix::luDecompose(Matrix& L, Matrix& U) const {
	size_t n = m_rows;

	for (size_t i = 0; i < n; ++i) {
		// 构造上三角矩阵 U
		for (size_t k = i; k < n; ++k) {
			double sum = 0;
			for (size_t j = 0; j < i; ++j) {
				sum += (L(i, j) * U(j, k));
			}
			U(i, k) = m_data[i * m_cols + k] - sum;
		}

		// 构造下三角矩阵 L
		for (size_t k = i; k < n; ++k) {
			if (i == k) {
				L(i, i) = 1; // 对角线元素为1
			}
			else {
				double sum = 0;
				for (size_t j = 0; j < i; ++j) {
					sum += (L(k, j) * U(j, i));
				}
				L(k, i) = (m_data[k * m_cols + i] - sum) / U(i, i);
				if (fabs(m_data[k * m_cols + i] - sum) < Epsilon_math && fabs(U(i, i) < Epsilon_math))L(k, i) = 1;
			}
		}
	}
}

/**
 * @file         "CMatrix.cpp"
 * @brief        
 * @return       bool
 * @note         是否可逆判断
 * @par History
 *               2024/12/4, Ken Ga  new \n
 */
bool Matrix::isInvertible() const {
	if (m_rows != m_cols) {
		throw std::invalid_argument("Matrix must be square to check invertibility.");
	}

	Matrix L(m_rows, m_cols), U(m_rows, m_cols);
	try {
		luDecompose(L, U);
	}
	catch (const std::runtime_error&) {
		return false;  // LU 分解失败，矩阵不可逆
	}

	// 如果 U 的对角线元素存在 0，则不可逆
	for (size_t i = 0; i < m_rows; ++i) {
		if (U(i, i) == 0) {
			return false;
		}
	}
	return true;
}


/**
 * @file         "CMatrix.cpp"
 * @brief        
 * @return       Matrix
 * @note         正定对称矩阵求逆
 * @par History
 *               2024/12/4, Ken Ga  new \n
 */
Matrix Matrix::inv_sp() const
{
	if (m_rows != m_cols) {
		throw std::invalid_argument("Matrix must be square to compute inverse.");
	}
	Matrix M(m_rows, m_cols, m_data);

	int j, k;

	for (k = 0; k < m_rows; k++)
	{
		if (m_data[k * m_rows + k] <= 0.0) return NULLMATRIX;
		if (M.m_data[k * m_rows + k] < DBL_EPSILON) M.m_data[k * m_rows + k] += 1e-15; // M矩阵接近奇异时,主对角线加上一个小量

		for (int i = 0; i < m_rows; i++)
		{
			if (i != k)   M.m_data[i * m_rows + k] = -M.m_data[i * m_rows + k] / M.m_data[k * m_rows + k];
		}
		M.m_data[k * m_rows + k] = 1.0 / M.m_data[k * m_rows + k];

		for (int i = 0; i < m_rows; i++)
		{
			if (i != k)
			{
				for (j = 0; j < m_rows; j++)
				{
					if (j != k)  M.m_data[i * m_rows + j] += M.m_data[k * m_rows + j] * M.m_data[i * m_rows + k];
				}
			}
		}
		for (j = 0; j < m_rows; j++)
		{
			if (j != k)  M.m_data[k * m_rows + j] *= M.m_data[k * m_rows + k];
		}
	}
	return M;
}

/**
 * @file         "CMatrix.cpp"
 * @brief        
 * @param[in]    excludeRow
 * @param[in]    excludeCol
 * @return       Matrix
 * @note         获得余子式
 * @par History
 *               2024/12/4, Ken Ga  new \n
 */
Matrix Matrix::getCofactor(int excludeRow, int excludeCol) const
{
	if (excludeRow>m_rows || excludeCol>m_cols) {
		throw std::invalid_argument("越界了.");
	}

	Matrix subMatrix(m_rows - 1, m_cols - 1);
	int sub_i = 0;
	for (int i = 0; i < m_rows; ++i) {
		if (i == excludeRow) continue;
		int sub_j = 0;
		for (int j = 0; j < m_cols; ++j) {
			if (j == excludeCol) continue;
			subMatrix(sub_i, sub_j) = (*this)(i, j);
			++sub_j;
		}
		++sub_i;
	}
	return subMatrix;
}

/**
 * @file         "CMatrix.cpp"
 * @brief        计算矩阵的行列式，超体积
 * @return       double
 * @note         计算矩阵的行列式，超体积
 * @par History
 *               2024/12/4, Ken Ga  new \n
 */
double Matrix::determinant() const
{
	if (m_rows != m_cols) throw std::runtime_error("行列式仅适用于方阵！");

	// 1x1 矩阵
	if (m_rows == 1) return (*this)(0, 0);

	// 2x2 矩阵
	if (m_rows == 2) return (*this)(0, 0) * (*this)(1, 1) - (*this)(0, 1) * (*this)(1, 0);

	// NxN 矩阵
	double det = 0;
	for (int col = 0; col < m_cols; ++col) {
		Matrix subMatrix = getCofactor(0, col);  // 获取子矩阵,余子式
		double cofactor = ((col % 2 == 0) ? 1 : -1) * (*this)(0, col);  // 余子式符号
		det += cofactor * subMatrix.determinant();  // 递归计算行列式
	}
	return det;
}


/**
 * @file         "CMatrix.cpp"
 * @brief        获得马氏距离
 * @param[in]    P  权阵  V是矩阵本身
 * @return       double
 * @note         使用动态数组，更底层的实现计算，简化运算
 * @par History
 *               2024/12/13, Ken Ga  new \n
 */
double Matrix::getMaDist(const Matrix P)
{
	return 0.0;
}


bool Matrix::choleskyDecompose(Matrix& L) const {
	if (m_rows != m_cols) {
		throw std::invalid_argument("Matrix must be square for Cholesky decomposition.");
	}

	L.resize(m_rows, m_cols);
	for (size_t i = 0; i < m_rows; ++i) {
		for (size_t j = 0; j <= i; ++j) {
			double sum = 0;

			for (size_t k = 0; k < j; ++k) {
				sum += L(i, k) * L(j, k);
			}

			if (i == j) {
				double value = m_data[i * m_cols + i] - sum;
				if (value <= 0) {
					return false;  // 不是正定矩阵
				}
				L(i, j) = std::sqrt(value);
			}
			else {
				L(i, j) = (m_data[i * m_cols + j] - sum) / L(j, j);
			}
		}
	}
	return true;
}



/**
 * @file         "CMatrix.cpp"
 * @brief        获得其逆矩阵
 * @return       Matrix
 * @note         
 * @par History
 *               2025/03/30, Ken Ga  new \n
 */
Matrix Matrix::inv() {

	if (m_rows != m_cols) {
		throw std::invalid_argument("Matrix must be square for inversion.");
	}

	Matrix result(m_rows, m_cols);

	// 创建一个临时矩阵 result.data 用于存储逆矩阵
	std::vector<int>is(m_rows), js(m_rows);   // 主元行列交换的索引
	int i, j, k, l, u, v;
	double d, p;

BEGIN:
	for (i = 0; i < m_rows; i++)
	{
		for (j = 0; j < m_rows; j++)
		{
			result.m_data[i * m_rows + j] = m_data[i * m_rows + j];
		}
	}

	for (k = 0; k < m_rows; k++)
	{
		d = 0.0;
		for (i = k; i < m_rows; i++)
		{
			for (j = k; j < m_rows; j++)
			{
				l = m_rows * i + j;
				p = fabs(result.m_data[l]);
				if (p > d)
				{
					d = p;
					is[k] = i;
					js[k] = j;
				}
			}
		}

		if (d < DBL_EPSILON)
		{
			if (d < 1e-30) {
				//throw std::invalid_argument("Matrix not be inv.");
				return NULLMATRIX;
			}

			for (i = 0; i < m_rows; i++)   m_data[i * (m_rows + 1)] += 1E-15; // data矩阵接近奇异时,主对角线加上一个小量
			goto BEGIN;
		}

		if (is[k] != k)
		{
			for (j = 0; j < m_rows; j++)
			{
				u = k * m_rows + j;
				v = is[k] * m_rows + j;
				p = result.m_data[u];
				result.m_data[u] = result.m_data[v];
				result.m_data[v] = p;
			}
		}

		if (js[k] != k)
		{
			for (i = 0; i < m_rows; i++)
			{
				u = i * m_rows + k;
				v = i * m_rows + js[k];
				p = result.m_data[u];
				result.m_data[u] = result.m_data[v];
				result.m_data[v] = p;
			}
		}

		l = k * m_rows + k;
		result.m_data[l] = 1.0 / result.m_data[l];
		for (j = 0; j < m_rows; j++)
		{
			if (j != k)
			{
				u = k * m_rows + j;
				result.m_data[u] = result.m_data[u] * result.m_data[l];
			}
		}
		for (i = 0; i < m_rows; i++)
		{
			if (i != k)
			{
				for (j = 0; j < m_rows; j++)
				{
					if (j != k)
					{
						u = i * m_rows + j;
						result.m_data[u] = result.m_data[u] - result.m_data[i * m_rows + k] * result.m_data[k * m_rows + j];
					}
				}
			}
		}
		for (i = 0; i < m_rows; i++)
		{
			if (i != k)
			{
				u = i * m_rows + k;
				result.m_data[u] = -result.m_data[u] * result.m_data[l];
			}
		}
	}

	for (k = m_rows - 1; k >= 0; k--)
	{
		if (js[k] != k)
		{
			for (j = 0; j < m_rows; j++)
			{
				u = k * m_rows + j;
				v = js[k] * m_rows + j;
				p = result.m_data[u];
				result.m_data[u] = result.m_data[v];
				result.m_data[v] = p;
			}
		}
		if (is[k] != k)
		{
			for (i = 0; i < m_rows; i++)
			{
				u = i * m_rows + k;
				v = is[k] + i * m_rows;
				p = result.m_data[u];
				result.m_data[u] = result.m_data[v];
				result.m_data[v] = p;
			}
		}
	}

	return result;
}

bool Matrix::LDLTcom(Matrix& L, Matrix& D) const
{
	if (m_rows != m_cols) { printf("%s", "非方阵，无法进行LDLT分解"); return false; }
	int row = m_rows;
	
	double Qahat_[MAXMATRIXSIZE * MAXMATRIXSIZE] = { 0.0 };
	toArr(Qahat_, row, row);
	double a = 0.0;

	for (int i = row - 1; i >= 0; i--)
	{
		if ((D.m_data[i] = Qahat_[i * row + i]) <= 0)
		{
			return false;
		}
		a = sqrt(Qahat_[i * row + i]);
		for (int k = 0; k <= i; k++) L.m_data[i * row + k] = Qahat_[i * row + k] / a;
		for (int j = 0; j <= i - 1; j++) for (int k = 0; k <= j; k++) Qahat_[j * row + k] = Qahat_[j * row + k] - L.m_data[i * row + k] * L.m_data[i * row + j];
		for (int k = 0; k <= i; k++) L.m_data[i * row + k] = L.m_data[i * row + k] / L.m_data[i * row + i];
	}

	for (int i = 0; i < row; i++)
	{
		if (D.m_data[i] < 1E-10)
		{
			return false;
		}
	}

	return true;
}


double Matrix::Dot(const Matrix& A, const Matrix& B)
{
	if (A.m_cols != 1 || B.m_cols != 1 || A.m_rows != B.m_rows)
	{
		throw std::out_of_range("Dot dimensions exceed matrix dimensions.");
		return 0;
	}
	return (A.Tra() * B).m_data[0];
}


// 提取子矩阵的方法

Matrix Matrix::getSubMatrix(size_t startRow, size_t startCol, size_t numRows, size_t numCols) const {
	if (startRow + numRows > m_rows || startCol + numCols > m_cols) {
		throw std::out_of_range("Submatrix dimensions exceed matrix dimensions.");
	}

	Matrix subMatrix(numRows, numCols);
	for (size_t i = 0; i < numRows; ++i) {
		for (size_t j = 0; j < numCols; ++j) {
			subMatrix(i, j) = m_data[(startRow + i) * m_cols + startCol + j];
		}
	}
	return subMatrix;
}

// 插入子矩阵的方法

void Matrix::setSubMatrix(size_t startRow, size_t startCol, const Matrix& subMatrix) {
	if (startRow + subMatrix.m_rows > m_rows || startCol + subMatrix.m_cols > m_cols) {
		throw std::out_of_range("Submatrix dimensions exceed matrix dimensions.");
	}
	for (size_t i = 0; i < subMatrix.m_rows; ++i) {
		for (size_t j = 0; j < subMatrix.m_cols; ++j) {
			m_data[(startRow + i) * m_cols + startCol + j] = subMatrix(i, j);
		}
	}
}


Matrix Matrix::exp(int n) const
{
	if (this->m_cols != this->m_rows) {
		throw std::out_of_range("matrix is not a square(非方阵).");
	}

	Matrix I = Identity(this->m_rows);  // 单位矩阵
	Matrix result = I;                // 结果矩阵，初始化为单位矩阵
	Matrix term = I;                  // 用于存储每一项的矩阵
	double factorial = 1;                     // 阶乘初值为1

	for (int i = 1; i <= n; ++i) {
		term = term * (*this);           // 计算矩阵的 i 次方
		factorial *= i;                  // 计算 i!
		result = result + term * (1 / factorial); // 累加到结果矩阵中
	}

	return result;  // 返回放缩的结果
}

//矩阵对角化
Matrix Matrix::symm() const {
	Matrix M(this->m_rows, this->m_cols);
	int i, j;
	int r = this->m_rows;
	int c = this->m_cols;
	for (i = 0; i < r; i++) for (j = 0; j < i; j++) M.m_data[i * c + j] = M.m_data[j * r + i] = 0.5 * (M.m_data[i * c + j] + M.m_data[j * r + i]);
	return M;
}

// 判断是否对角线元素都大于0
bool Matrix::tr_positive() const
{
	if (m_cols!=m_rows)
	{
		return false;   //说明非方阵
	}
	for (size_t i = 0; i < m_cols; i++)
	{
		if (m_data[i * m_cols + i] < 0) {
			return false;
		}
	}

	return true;
}

/**
 * @file         "CMatrix.cpp"
 * @brief        转成数组
 * @param[in]    Arr	接受的数组
 * @param[in]    rows	要截取的行列数
 * @param[in]    cols	
 * @return       bool	
 * @note         
 * @par History
 *               2025/01/17, Ken Ga  new \n
 */
bool Matrix::toArr(double Arr[],const size_t rows, const size_t cols)const
{
	if (rows <= 0 || cols <= 0)return false;
	for (size_t i = 0; i < rows; i++)
	{
		for (size_t j = 0; j < cols; j++) {
			Arr[i * cols + j] = m_data[i * cols + j];
		}
	}
	return true;
}
/**
 * @file         "CMatrix.cpp"
 * @brief        取某行的元素转为数列
 * @param[in]    Arr
 * @param[in]    row
 * @return       bool
 * @note         
 * @par History
 *               2025/01/17, Ken Ga  new \n
 */
bool Matrix::rowtoArr(double Arr[], const size_t row)
{
	for (size_t i = 0; i < m_cols; i++)
	{
		Arr[i] = at(row, i);
	}
	return true;
}

/**
 * @file         "CMatrix.cpp"
 * @brief        取某列元素转为数组
 * @param[in]    Arr
 * @param[in]    row
 * @return       bool
 * @note         
 * @par History
 *               2025/01/17, Ken Ga  new \n
 */
bool Matrix::coltoArr(double Arr[], const size_t row)
{
	for (size_t i = 0; i < m_rows; i++)
	{
		Arr[i] = at(i, m_cols);
	}
	return true;
}
/**
 * @file         "CMatrix.cpp"
 * @brief        直接对矩阵的内部元素进行更改，将矩阵降为一维数组，对数组内元素的固定位置开始填充固定大小的数组元素
 * @param[in]    Arr	填充数组
 * @param[in]    begin	填充的开始位置
 * @param[in]    size	填充的大小
 * @return       bool	返回真值表示填充成功，其他为失败
 * @note         
 * @par History
 *               2025/01/17, Ken Ga  new \n
 */
bool Matrix::Arr2data(const double Arr[], const int begin, const int size)
{
	if (begin + size > m_rows * m_cols)return false;
	if (begin == -1 && size == -1) {
		for (size_t i = 0; i < m_rows*m_cols; i++)
		{
			m_data[i] = Arr[i];
		}
		return true;
	}
	else if (begin < 0 || size < 0)return false;
	for (size_t i = begin; i < begin+size; i++)
	{
		m_data[i] = Arr[i];
	}
	return false;
}

/**
 * @file         "CMatrix.cpp"
 * @brief        直接对矩阵的内部元素进行更改，将矩阵降为一维数组，对数组内元素的固定位置开始填充固定大小的数组元素
 * @param[in]    Arr
 * @param[in]    begin
 * @param[in]    size
 * @return       bool
 * @note         
 * @par History
 *               2025/01/17, Ken Ga  new \n
 */
bool Matrix::Arr2data(const int Arr[], const int begin, const int size)
{
	if (begin + size > m_rows * m_cols)return false;
	if (begin == -1 && size == -1) {
		for (size_t i = 0; i < m_rows * m_cols; i++)
		{
			m_data[i] = Arr[i];
		}
		return true;
	}
	else if (begin < 0 || size < 0)return false;
	for (size_t i = begin; i < begin + size; i++)
	{
		m_data[i] = Arr[i];
	}
	return false;
}



/**
 * @file         "CMatrix.cpp"
 * @brief        重设矩阵大小
 * @param[in]    newRows 新矩阵的列
 * @param[in]    newCols 新矩阵的行
 * @param[in]    reserve 用来进行判定原来矩阵中的值是否保留
 * @return       void
 * @note         
 * @par History
 *               2024/12/4, Ken Ga  new \n
 */
void Matrix::resize(const size_t newRows, const size_t newCols, bool reserve) {
	if (newRows == m_rows && newCols == m_cols && reserve) {
		return;  // 如果大小没变，且保留原来值，直接返回
	}

	// 创建新的矩阵数据，并全部初始化为 0
	std::vector<double> newData(newRows * newCols, double{});

	if (reserve) {
		// 计算要拷贝的最小行数和列数，防止越界
		size_t minRows = std::min(newRows, m_rows);
		size_t minCols = std::min(newCols, m_cols);

		// 将旧数据拷贝到新矩阵中
		for (size_t i = 0; i < minRows; ++i) {
			for (size_t j = 0; j < minCols; ++j) {
				newData[i * newCols + j] = m_data[i * m_cols + j];
			}
		}
	}

	// 更新矩阵的行列和数据
	m_data = std::move(newData);  // 使用移动语义避免额外的拷贝
	m_rows = newRows;
	m_cols = newCols;
}

/**
 * @file         "CMatrix.cpp"
 * @brief        加入了用vector可变数组进行创建矩阵
 * @param[in]    r
 * @param[in]    c
 * @param[in]    data
 * @return       void
 * @note         
 * @par History
 *               2025/01/15, Ken Ga  new \n
 */
Matrix::Matrix(size_t r, size_t c, std::vector<double> vec)
{
	m_data.resize(r * c);
	m_cols = c;
	m_rows = r;
	for (size_t i = 0; i < r*c; i++)
	{
		m_data[i] = vec[i];
	}
}

/**
 * @file         "CMatrix.cpp"
 * @brief        利用列表初始化进行对矩阵的赋值 类似A={{1,23},{3,2}};
 * @param[in]    init
 * @return       void
 * @note         
 * @par History
 *               2024/12/4, Ken Ga  new \n
 */
Matrix::Matrix(std::initializer_list<std::initializer_list<double>> init)
{
	m_rows = init.size();
	m_cols = init.begin()->size();
	m_data.reserve(m_rows * m_cols);

	for (auto& row : init) {
		m_data.insert(m_data.end(), row.begin(), row.end());
	}
}

/**
 * @file         "CMatrix.cpp"
 * @brief        利用双参数和列表进行初始化，类似A(3,1,{1,2,3});
 * @param[in]    r
 * @param[in]    c
 * @param[in]    init
 * @return       void
 * @note         
 * @par History
 *               2025/01/16, Ken Ga  new \n
 */
Matrix::Matrix(size_t r, size_t c, std::initializer_list<double> init)
{
	m_rows = r;
	m_cols = c;

	if (init.size()>0 && init.size() != r * c) {
		throw std::invalid_argument("Initializer list size does not match specified dimensions.");
	}

	m_data.reserve(r * c);
	if (init.size()>0) {
		m_data.insert(m_data.end(), init.begin(), init.end()); // 如果提供了列表数据，插入数据
	}
	else {
		m_data.resize(r * c, 0.0); // 否则用默认值填充
	}

}


/**
 * @file         "CMatrix.cpp"
 * @brief        拷贝构造获得矩阵的子矩阵,任意维度
 * @param[in]    a
 * @param[in]    srow	start行
 * @param[in]    scol	start列
 * @param[in]    nrow	子矩阵的行
 * @param[in]    ncol	子矩阵的列
 * @return       void	构造函数
 * @note         
 * @par History
 *               2024/12/4, Ken Ga  new \n
 */
Matrix::Matrix(const Matrix& a, const size_t srow, const size_t scol, const size_t nrow, const size_t ncol)
{
	if (srow + nrow > a.getRows() || scol + ncol > a.getCols())
	{
		printf("submatrix exceeds the original matrix\n");
		return;
	}
	m_rows = nrow;
	m_cols = ncol;
	m_data.resize(nrow * ncol, 0.0);
	size_t startIdx = srow * a.getCols() + scol;
	for (size_t i = 0; i < nrow; i++)
	{
		for (size_t j = 0; j < ncol; j++) {
			m_data[i * ncol + j] = a(srow + i, scol + j);
		}
	}
}

/**
 * @file         "CMatrix.cpp"
 * @brief        拷贝构造获得矩阵的子方阵,固定为方阵
 * @param[in]    a
 * @param[in]    srow
 * @param[in]    scol
 * @param[in]    ndimensin
 * @return       void
 * @note         
 * @par History
 *               2024/12/4, Ken Ga  new \n
 */
Matrix::Matrix(const Matrix& a, const size_t srow, const size_t scol, const size_t ndimensin)
{
	if (srow + ndimensin > a.getRows() || scol + ndimensin > a.getCols())
	{
		printf("submatrix exceeds the original matrix\n");
		return;
	}
	m_rows = ndimensin;
	m_cols = ndimensin;
	m_data.resize(m_rows * m_cols, 0.0);
	size_t startIdx = srow * a.getCols() + scol;
	for (size_t i = 0; i < ndimensin; i++)
	{
		for (size_t j = 0; j < ndimensin; j++) {
			m_data[i * ndimensin + j] = a(srow + i, scol + j);
		}
	}
}

Matrix& Matrix::operator=(const Matrix& p)
{
	if (this != &p) {  // 防止自赋值
		m_rows = p.m_rows;
		m_cols = p.m_cols;
		m_data = p.m_data;   //这里因为vector自带深拷贝，所以两个矩阵指向不同的内存，不会出现A发生修改，影响B的情况！
	}
	return *this;
}

Matrix& Matrix::operator=(std::initializer_list<double> init)
{
	if (init.size() != m_rows * m_cols) {
		throw std::invalid_argument("Initializer list size does not match matrix dimensions.");
	}
	m_data.clear();
	m_data.insert(m_data.end(), init.begin(), init.end());
	return *this;
}

double& Matrix::at(const size_t i, const size_t j)
{
	if (i >= m_rows || j >= m_cols) {
		throw std::invalid_argument("Matrix dimensions do not match for at.");
	}
	return m_data[i * m_cols + j];
}

double& Matrix::operator()(const size_t i, const size_t j)
{
	// 参数检查，防止越界访问
	if (i >= m_rows || j >= m_cols) {
		throw std::out_of_range("Matrix index out of bounds.");
	}
	return m_data[i * m_cols + j];  // 直接访问 data 数组
}

const double& Matrix::operator()(const size_t i, const size_t j) const
{
	// 参数检查，防止越界访问
	if (i >= m_rows || j >= m_cols) {
		throw std::out_of_range("Matrix index out of bounds.");
	}
	return m_data[i * m_cols + j];  // 直接访问 data 数组
}

/**
 * @file         "CMatrix.cpp"
 * @brief        返回矩阵大小
 * @return       int
 * @note         
 * @par History
 *               2025/01/16, Ken Ga  new \n
 */
int Matrix::size() const
{
	return m_rows * m_cols;
}

/**
 * @file         "CMatrix.cpp"
 * @brief        打印矩阵
 * @param[in]    precflo	小数点后精度
 * @return       void		
 * @note         
 * @par History
 *               2025/01/16, Ken Ga  new \n
 */
void Matrix::print(int precflo) const
{
	// 第一步，找到最大宽度

	int maxWidth = 0;   //采用动态宽度
	for (size_t i = 0; i < m_rows; ++i) {
		for (size_t j = 0; j < m_cols; ++j) {
			// 计算当前元素的字符长度
			int width = (m_data[i*m_cols+j] > 10) ? log10(m_data[i*m_cols+j]) + precflo + 1 : int(log10(m_data[i*m_cols+j])) + 2 + precflo;
			if (m_data[i*m_cols+j] < 0) { // 负数需要多一个字符用于负号
				width++;
			}
			maxWidth = std::max(maxWidth, width);
		}
	}

	// 第二步，输出矩阵
	for (size_t i = 0; i < m_rows; ++i) {
		for (size_t j = 0; j < m_cols; ++j) {
			std::cout << std::setw(maxWidth + 2) << std::fixed << std::setprecision(precflo) << m_data[i*m_cols+j];
		}
		std::cout << std::endl;
	}
	printf("\n");  //防止多个输出矩阵连起来
}

/**
 * @file         "CMatrix.cpp"
 * @brief        创建一个零矩阵
 * @return       void
 * @note         
 * @par History
 *               2025/01/16, Ken Ga  new \n
 */
void Matrix::zero()
{
	std::fill(m_data.begin(), m_data.end(), 0.0);
}

/**
 * @file         "CMatrix.cpp"
 * @brief        创建一个填充指定值的矩阵
 * @param[in]    size
 * @param[in]    value
 * @return       Matrix
 * @note         
 * @par History
 *               2025/01/16, Ken Ga  new \n
 */
Matrix Matrix::fill(const int size, const double value)
{
	Matrix Fill(size, size);
	std::fill(Fill.m_data.begin(), Fill.m_data.end(), value);
	return Fill;
}

/**
 * @file         "CMatrix.cpp"
 * @brief        矩阵所有值赋值1
 * @param[in]    value
 * @return       void
 * @note         
 * @par History
 *               2025/01/16, Ken Ga  new \n
 */
void Matrix::fill(const double value)
{
	std::fill(m_data.begin(), m_data.end(), value);
}
/**
 * @file         "CMatrix.cpp"
 * @brief        判断是否为0矩阵，如果为0矩阵(元素全为0)，返回真值
 * @return       bool
 * @note         
 * @par History
 *               2025/01/17, Ken Ga  new \n
 */
bool Matrix::isZero()
{
	for (size_t i = 0; i < m_rows*m_cols; i++)
	{
		if (!fabs(m_data[i]) < Epsilon_math)return false;
	}
	return true;
}

/**
 * @file         "CMatrix.cpp"
 * @brief        对角矩阵构建1
 * @param[in]    size	矩阵维数
 * @param[in]    values	等角值
 * @return       Matrix
 * @note
 * @par History
 *               2025/01/16, Ken Ga  new \n
 */
Matrix Matrix::diag(const int size, const double value)
{
	Matrix diag(size, size);
	for (size_t i = 0; i < size; i++) {
		diag.m_data[i * size + i] = value;
	}
	return diag;
}
/**
 * @file         "CMatrix.cpp"
 * @brief        对角矩阵构建2
 * @param[in]    size
 * @param[in]    values
 * @return       Matrix
 * @note         
 * @par History
 *               2025/01/16, Ken Ga  new \n
 */
Matrix Matrix::diag(const int size, const double values[])
{
	Matrix diag(size, size);
	for (size_t i = 0; i < size; i++) {
		diag.m_data[i * size + i] = values[i];
	}
	return diag;
}
/**
 * @file         "CMatrix.cpp"
 * @brief        创建单位阵，
 * @param[in]    size	输入的是维数 size*size
 * @return       Matrix
 * @note         
 * @par History
 *               2025/01/16, Ken Ga  new \n
 */
Matrix Matrix::Identity(const int size)
{
	Matrix Identity(size, size);
	for (size_t i = 0; i < size; i++) {
		Identity.m_data[i * size + i] = 1;
	}
	return Identity;
}

/**
 * @file         "CMatrix.cpp"
 * @brief        矩阵求迹
 * @return       double
 * @note         
 * @par History
 *               2025/01/16, Ken Ga  new \n
 */
double Matrix:: tr() const {
	if (m_rows != m_cols) {
		throw std::invalid_argument("Matrix dimensions do not match for tr.");
	}
	double result{};
	for (size_t i = 0; i < m_rows; i++) {
		result += m_data[i * m_rows + i];
	}
	return result;
}

/**
 * @file         "CMatrix.cpp"
 * @brief        
 * @return       int
 * @note         获取列
 * @par History
 *               2025/01/16, Ken Ga  new \n
 */
int Matrix::getCols() const
{
	return m_cols;
}

/**
 * @file         "CMatrix.cpp"
 * @brief        获取行
 * @return       int
 * @note         
 * @par History
 *               2025/01/16, Ken Ga  new \n
 */
int Matrix::getRows() const
{
	return m_rows;
}

/**
 * @file         "CMatrix.cpp"
 * @brief        计算最小二乘解，（BTPB）^-1*BTPL
 * @param[in]    B
 * @param[in]    L
 * @param[in]    P
 * @return       Matrix
 * @note         
 * @par History
 *               2025/01/16, Ken Ga  new \n
 */
Matrix LeastSqure(const Matrix& B, const Matrix& L, const Matrix* P)
{
	if (!P)return (B.Tra() * B).inv_sp() * B.Tra() * L;
	return (B.Tra()*(*P)* B).inv_sp() * B.Tra()*(*P) * L;;
}

/**
 * @file         "CMatrix.cpp"
 * @brief        实现C语言(数组形式的)矩阵打印
 * @param[in]    data
 * @param[in]    rows
 * @param[in]    cols
 * @return       void
 * @note         
 * @par History
 *               2025/01/16, Ken Ga  new \n
 */
void PrintMatrix(const double* data, const int rows, const int cols)
{
	char precision = 3;   //小数点后精度
	int maxWidth = 0;
	// 第一步，找到最大宽度
	for (size_t i = 0; i < rows; ++i) {
		for (size_t j = 0; j < cols; ++j) {
			char buffer[50]; // 临时缓冲区存储格式化字符串
			snprintf(buffer, sizeof(buffer), "%.*f", precision, data[i * cols + j]);
			int width = (int)strlen(buffer);
			if (data[i * cols + j] < 0) {
				width++; // 负号需要额外宽度
			}
			maxWidth = width > maxWidth ? width : maxWidth;
		}
	}
	// 第二步，打印矩阵
	for (size_t i = 0; i < rows; ++i) {
		for (size_t j = 0; j < cols; ++j) {
			printf("%*.*f ", maxWidth + 2, precision, data[i * cols + j]);
		}
		printf("\n");  //防止多个输出矩阵连起来
	}
	printf("\n");  //防止多个输出矩阵连起来
}


bool PrintMatrixToFile(const double* data, const int rows, const int cols,char mode,const char precision)
{
	int maxWidth = 0;
	// 第一步，找到最大宽度
	for (size_t i = 0; i < rows; ++i) {
		for (size_t j = 0; j < cols; ++j) {
			char buffer[50]; // 临时缓冲区存储格式化字符串
			snprintf(buffer, sizeof(buffer), "%.*f", precision, data[i * cols + j]);
			int width = (int)strlen(buffer);
			if (data[i * cols + j] < 0) {
				width++; // 负号需要额外宽度
			}
			maxWidth = width > maxWidth ? width : maxWidth;
		}
	}

	static FILE* fpmat = NULL;
	if (!fpmat) {
		if (mode == 'w')fopen_s(&fpmat, "PirntMatrix.txt", "w");
		if (mode == 'a')fopen_s(&fpmat, "PirntMatrix.txt", "a");
	}
	if (!fpmat)return false;

	// 第二步，打印矩阵
	for (size_t i = 0; i < rows; ++i) {
		for (size_t j = 0; j < cols; ++j) {
			fprintf(fpmat,"%*.*f ", maxWidth + 2, precision, data[i * cols + j]);
		}
		fprintf(fpmat,"\n");  //防止多个输出矩阵连起来
	}
	fprintf(fpmat,"\n");  //防止多个输出矩阵连起来
	return true;
}

//void ArrAlgo::MN1EQU(const double M[], double M2[],const int len)
//{
//	for (size_t i = 0; i < len; i++)
//	{
//		M2[i] = M[i];
//	}
//}
//
//void ArrAlgo::MN1Scale(const double a, double M[], const int length)
//{
//	for (size_t i = 0; i < length; i++)
//	{
//		M[i] = a * M[i];
//	}
//}
//
//void ArrAlgo::MN1Scale(const double a, const double M[], double M2[], const int length)
//{
//	for (size_t i = 0; i < length; i++)
//	{
//		M2[i] = a * M[i];
//	}
//}
//
///**
// * @file         "CMatrix.cpp"
// * @brief        将数组所有元素取倒数
// * @param[in]    M
// * @param[in]    M2
// * @param[in]    length
// * @return       bool	如果数组有元素为 0，返回false
// * @note         
// * @par History
// *               2025/01/17, Ken Ga  new \n
// */
//bool ArrAlgo::MN1Reci(const double M[], double M2[], const int length)
//{
//	bool status = false;
//	for (size_t i = 0; i < length; i++)
//	{
//		if(isZero(M[i]))status=false;
//		M2[i] = 1.0 / M[i];
//	}
//	return true;
//}
//
///**
// * @file         "BaseGT.cpp"
// * @brief        计算二维距离
// * @param[in]    A
// * @param[in]    B
// * @return       double
// * @note
// * @par History
// *               2024/12/12, Ken Ga  new \n
// */
//double ArrAlgo::L2d(const double A[], const double B[])
//{
//	double sum = 0;
//	for (size_t i = 0; i < 2; i++)
//	{
//		sum += pow(A[i] - B[i], 2);
//	}
//	return sqrt(sum);
//}
//
///**
// * @file         "BaseGT.cpp"
// * @brief        计算三维距离
// * @param[in]    A
// * @param[in]    B
// * @return       double
// * @note
// * @par History
// *               2024/12/12, Ken Ga  new \n
// */
//double ArrAlgo::L3d(const double A[], const double B[])
//{
//	double sum = 0;
//	for (size_t i = 0; i < 3; i++)
//	{
//		sum +=pow(A[i] - B[i], 2);
//	}
//	return sqrt(sum);
//}
//
//
///**
// * @file         "BaseGT.cpp"
// * @brief        二次型，计算马氏距离，vTPv
// * @param[in]    v	残差
// * @param[in]    P	协方差
// * @param[in]    n	维数
// * @return       double
// * @note         未检验!!!
// * @par History
// *               2024/12/12, Ken Ga  new \n
// */
//double ArrAlgo::Md(const double v[], const double P[], const int n)
//{
//	// 计算 v^T P
//	double* vTp = new double[n];
//	for (int i = 0; i < n; ++i) {
//		vTp[i] = 0;
//		for (int j = 0; j < n; ++j) {
//			vTp[i] += v[j] * P[i * n + j];  // P[i][j] 是按行主序存储的
//		}
//	}
//
//	// 计算 (v^T P) v
//	double result = 0.0;
//	for (int i = 0; i < n; ++i) {
//		result += vTp[i] * v[i];
//	}
//
//	// 释放内存
//	delete[] vTp;
//
//	return result;
//}
//
//double ArrAlgo::mean(const double A[], const int length)
//{
//	double sum = 0.0;
//	for (size_t i = 0; i < length; i++)
//	{
//		sum += A[i];
//	}
//	return sum/length;
//}
//
//double ArrAlgo::std(const double A[], const int length)
//{
//	double m = 0.0, sum = 0.0;
//	m = mean(A, length);
//	for (size_t i = 0; i < length; i++)
//	{
//		sum += pow((A[i] - m), 2);
//	}
//	return sqrt(sum/length);
//}
//
//
//double ArrAlgo::Dot(const double v1[], const double v2[],const int n)
//{
//	double c = 0.0;
//	for (int i = 0; i < n;i++) c += v1[i] * v2[i];
//	return c;
//}
//
///**
// * @file         "CMatrix.cpp"
// * @brief        统计浮点的两个数组中，有多少个元素相等
// * @param[in]    arr
// * @param[in]    arr2
// * @param[in]    n
// * @return       int	返回有多个元素相等
// * @note         统计浮点的两个数组中，有多少个元素相等
// * @par History
// *               2025/01/17, Ken Ga  new \n
// */
//int ArrAlgo::countIntArrEQU(const int arr[], const int arr2[], size_t n)
//{
//	int sum = 0;
//	for (size_t i = 0; i < n; i++)
//	{
//		if (arr[i] == arr2[i])sum++;
//	}
//	return sum;
//}

/**
 * @file         "CMatrix.cpp"
 * @brief        判断一个浮点数是否为0，精度由第二个参数决定
 * @param[in]    x
 * @param[in]    epsilon
 * @return       bool
 * @note         
 * @par History
 *               2025/01/17, Ken Ga  new \n
 */
bool isZero(const double x, const double epsilon)
{
	if (fabs(x) < epsilon)return true;
	return false;
}

/**
 * @file         "CMatrix.cpp"
 * @brief        获取输入数字a与最近整数的距离绝对值
 * @param[in]    x
 * @return       bool
 * @note         值只在0-0.5之间
 * @par History
 *               2025/01/17, Ken Ga  new \n
 */
bool getFrac(const double x)
{
	return fabs(x - round(x));
}

/********************************
 IPS矩阵运算库
/********************************

/** @brief  Initialize assignment (double) */
void SetValsD(double* dst, double val, int n)
{
	if (!dst)return;
	for (int i = 0; i < n; i++) dst[i] = val;
}


/** @brief  Initialize assignment (int) */
void SetValsI(int* dst, int val, int n)
{
	if (!dst)return;
	for (int i = 0; i < n; i++) dst[i] = val;
}


/** @brief  Initialize assignment (bool) */
void SetValsB(bool* dst, bool val, int n)
{
	if (!dst)return;
	for (int i = 0; i < n; i++) dst[i] = val;
}


/** @brief  Initialize assignment (char) */
void SetValsSC(signed char* dst, signed char val, int n)
{
	if (!dst)return;
	for (int i = 0; i < n; i++) dst[i] = val;
}


/** @brief  Initialize assignment (char) */
void SetValsUC(unsigned char* dst, unsigned char val, int n)
{
	if (!dst)return;
	for (int i = 0; i < n; i++) dst[i] = val;
}


/** @brief  Assigns matrix : M2 = M */
void M31EQU(const double M[3], double M2[3])
{
	for (int i = 0; i < 3; i++) M2[i] = M[i];
}


/** @brief  Assigns matrix : M2 = -M */
void M31EQUM(const double M[3], double M2[3])
{
	for (int i = 0; i < 3; i++) M2[i] = -M[i];
}


/** @brief  Assigns matrix : M2 = M */
void M41EQU(const double M[4], double M2[4])
{
	for (int i = 0; i < 4; i++) M2[i] = M[i];
}


/** @brief Assigns matrix : M2 = M */
void M33EQU(const double M[9], double M2[9])
{
	for (int i = 0; i < 9; i++) M2[i] = M[i];
}


/** @brief Assigns Vector (double) : V2 = V */
void VecEQUD(const int n, const double V[], double V2[])
{
	for (int i = 0; i < n; i++) V2[i] = V[i];
}


/** @brief Assigns Vector (int) : V2 = V */
void VecEQUI(const int n, const int V[], int V2[])
{
	for (int i = 0; i < n; i++) V2[i] = V[i];
}


/** @brief Assigns Vector (bool) : V2 = V */
void VecEQUB(const int n, const bool V[], bool V2[])
{
	for (int i = 0; i < n; i++) V2[i] = V[i];
}


/** @brief Set M to 0 */
void M31Zero(double M[3])
{
	for (int i = 0; i < 3; i++) M[i] = 0.0;
}


/** @brief Set M to 0 */
void M33Zero(double M[9])
{
	for (int i = 0; i < 9; i++) M[i] = 0.0;
}


/** @brief  Scale matrix : M2 = a*M */
void M31Scale(double a, double M[3])
{
	for (int i = 0; i < 3; i++) M[i] *= a;
}


/** @brief  Scale matrix : M2 = a*M */
void M33Scale(double a, double M[9])
{
	for (int i = 0; i < 9; i++) M[i] *= a;
}


/** @brief  Adds two matrices : M3 = s1*M1 + s2*M2 */
void M31M311(const double M1[3], const double M2[3], double M3[3], double s1, double s2)
{
	for (int i = 0; i < 3; i++) M3[i] = s1 * M1[i] + s2 * M2[i];
}


/** @brief  Adds two matrices : M2 = s1*M1 + s2*M2 */
void M31M312(const double M1[3], double M2[3], double s1, double s2)
{
	for (int i = 0; i < 3; i++) M2[i] = s1 * M1[i] + s2 * M2[i];
}


/** @brief  Adds two matrices : M3 = s1*M1 + s2*M2 */
void M33M331(const double M1[9], const double M2[9], double M3[9], double s1, double s2)
{
	for (int i = 0; i < 9; i++) M3[i] = s1 * M1[i] + s2 * M2[i];
}


/** @brief  Adds two matrices : M2 = s1*M1 + s2*M2 */
void M33M332(const double M1[9], double M2[9], double s1, double s2)
{
	for (int i = 0; i < 9; i++) M2[i] = s1 * M1[i] + s2 * M2[i];
}


/** @brief  Subtracts two matrices : M3 = M1 - M2 */
void M31_M31(const double M1[3], const double M2[3], double M3[3])
{
	for (int i = 0; i < 3; i++) M3[i] = M1[i] - M2[i];
}


/** @brief  Multiplys two matrices : M31_ = M33 * M31 */
void M33XM31(const double M33[9], const double M31[3], double M31_[3])
{
	M31_[0] = M33[0] * M31[0] + M33[1] * M31[1] + M33[2] * M31[2];
	M31_[1] = M33[3] * M31[0] + M33[4] * M31[1] + M33[5] * M31[2];
	M31_[2] = M33[6] * M31[0] + M33[7] * M31[1] + M33[8] * M31[2];
}


/** @brief  Multiplys two matrices : M33_3 = M33_1 * M33_2 */
void M33XM33(const double M33_1[9], const double M33_2[9], double M33_3[9])
{
	int i, j, k;
	double Sum;
	for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) {
		Sum = 0.0; for (k = 0; k < 3; k++) Sum = Sum + M33_1[i * 3 + k] * M33_2[k * 3 + j]; M33_3[i * 3 + j] = Sum;
	}
}


/** @brief  Cross two vectors : v3 = v1 X v2 */
void CrossM3(const double v1[3], const double v2[3], double v3[3])
{
	v3[0] = v1[1] * v2[2] - v1[2] * v2[1];
	v3[1] = v1[2] * v2[0] - v1[0] * v2[2];
	v3[2] = v1[0] * v2[1] - v1[1] * v2[0];
}


/** @brief  dot product of two vectors */
double MatrixDot1(int n, const double v1[], const double v2[])
{
	double c = 0.0;
	while (--n >= 0) c += v1[n] * v2[n];
	return c;
}


/*
* @brief        Get Skew-Symmetric Matrix from vector
* @param[in]    v    double[3]    Vector
* @param[out]   M    double[9]    Skew-Symmetric Matrix, 3 x 3
* @return       void
* @note         无
* @internals    无
*/
void MatrixSkewSymmetric(const double v[3], double M[9])
{
	M[0] = 0.0;  M[1] = -v[2]; M[2] = v[1];
	M[3] = v[2]; M[4] = 0.0;  M[5] = -v[0];
	M[6] = -v[1]; M[7] = v[0]; M[8] = 0.0;
}


/*
* @brief        Get vector from Skew-Symmetric Matrix
* @param[in]    M    double[9]    Skew-Symmetric Matrix, 3 x 3
* @param[out]   v    double[3]    Vector
* @return       void
* @note         无
* @internals    无
*/
void MatrixSkewSymmetric_T(const double M[9], double v[3])
{
	v[0] = M[7]; v[1] = M[2]; v[2] = M[3];
}


/** @brief  Matrix Symmetrization */
void MatrixSymmetrization(int r, int c, double M[])
{
	int i, j;
	for (i = 0; i < r; i++) for (j = 0; j < i; j++) M[i * c + j] = M[j * r + i] = 0.5 * (M[i * c + j] + M[j * r + i]);
}


/**
* @brief        Copy Matrix
* @param[in]    r      int        The number of rows in the matrix M
* @param[in]    c      int        The number of cols in the matrix M
* @param[in]    src    double*    Src matrix
* @param[out]   dst    double*    Dest matrix
* @return       void
* @note         无
* @internals    无
*/
void MatrixCopy(int r, int c, const double* src, double* dst)
{
	int n = r * c, i;

	for (i = 0; i < n; i++) dst[i] = src[i];
}


/**
* @brief        Copy a sub-matrix of the matrix
* @param[in]    s1        int        (s1[0],s1[1]) is the starting point of M1
* @param[in]    s2        int        (s2[0],s2[1]) is the starting point of M2
* @param[in]    nRow      int        The number of rows that is copied
* @param[in]    nCol      int        The number of cols that is copied
* @param[in]    r1        int        The number of rows in the matrix M1
* @param[in]    c1        int        The number of cols in the matrix M1
* @param[in]    M1        double*    Src matrix
* @param[in]    r2        int        The number of rows in the matrix M2
* @param[in]    c2        int        The number of cols in the matrix M2
* @param[out]   M2        double*    Dst matrix
* @param[in]    bMinus    bool       true = the sub-matrix will take a negative sign
* @return       void
* @note         The new matrix is given by M2[(s2[0] + i)*c2 + s2[1] + j] = M1[(s1[0] + i)*c1 + s1[1] + j];
* @internals
*/
void MatrixCopySub(const int s1[2], const int s2[2], int nRow, int nCol, int r1, int c1, const double* M1, int r2, int c2, double* M2, bool bMinus)
{
	if (s1[0] + nRow > r1) return;
	if (s1[1] + nCol > c1) return;
	if (s2[0] + nRow > r2) return;
	if (s2[1] + nCol > c2) return;

	int i, j;

	if (bMinus)
	{
		for (i = 0; i < nRow; i++)
		{
			for (j = 0; j < nCol; j++)
			{
				M2[(s2[0] + i) * c2 + s2[1] + j] = -M1[(s1[0] + i) * c1 + s1[1] + j];
			}
		}
	}
	else
	{
		for (i = 0; i < nRow; i++)
		{
			for (j = 0; j < nCol; j++)
			{
				M2[(s2[0] + i) * c2 + s2[1] + j] = M1[(s1[0] + i) * c1 + s1[1] + j];
			}
		}
	}

}


/**
* @brief        Swaps two matrix
* @param[in]    r      int        The number of rows in the matrix M
* @param[in]    c      int        The number of cols in the matrix M
* @param[in]    src    double*    Src matrix
* @param[out]   dst    double*    Dst matrix
* @return       void
* @note         无
* @internals    无
*/
void MatrixSwap(int r, int c, double* src, double* dst)
{
	int n = r * c, i;
	double tmp;

	for (i = 0; i < n; i++)
	{
		tmp = dst[i];
		dst[i] = src[i];
		src[i] = tmp;
	}
}


/**
* @brief          Sets all elements of a matrix equal to zero
* @param[in]      r    int         The number of rows in the matrix M
* @param[in]      c    int         The number of cols in the matrix M
* @param[in,out]  M    double[]    Matrix
* @return         void
* @note           无
* @internals      无
*/
void MatrixZero(int r, int c, double M[])
{
	if (M && r > 0 && c > 0)
	{
		memset(M, 0, sizeof(double) * r * c);
	}
}


/**
* @brief          Sets the diagonal elements of a matrix equal to zero
* @param[in]      r    int         The number of rows in the matrix M
* @param[in]      c    int         The number of cols in the matrix M
* @param[in,out]  M    double[]    Matrix
* @return         void
* @note           its main diagonal elements will be equal to 1.0 and the rest zero
* @internals      无
*/
void MatrixZero_diag(int r, int c, double M[])
{
	for (int i = 0; i < r; i++)
	{
		M[(r + 1) * i] = 0.0;
	}
}


/**
* @brief          Sets the matrix as a unit matrix
* @param[in]      r    int         The number of rows in the matrix M
* @param[in]      c    int         The number of cols in the matrix M
* @param[in,out]  M    double[]    Matrix
* @return         void
* @note           its main diagonal elements will be equal to 1.0 and the rest zero
* @internals      无
*/
void MatrixUnit(int r, int c, double M[])
{
	if (M && r > 0 && c > 0 && r == c)
	{
		memset(M, 0, sizeof(double) * r * c);
		for (int i = 0; i < r; i++) M[i * r + i] = 1.0;
	}
}


/**
* @brief        Matrix Trace
* @param[in]    r       int         The number of rows in the matrix M
* @param[in]    c       int         The number of cols in the matrix M
* @param[in]    M       double[]    Matrix
* @return       double  Returns the sum of diagonal elements of a matrix
* @note         无
* @internals    无
*/
double MatrixTrace(int r, int c, const double M[])
{
	int i;
	double b = 0;
	for (i = 0; i < r; i++) b += M[i * r + i];
	return b;
}


/**
* @brief        Matrix Norm2
* @param[in]    r       int         The number of rows in the matrix M
* @param[in]    c       int         The number of cols in the matrix M
* @param[in]    M       double[]    Matrix
* @return       double  Returns the two norm of a matrix
* @note         Either r or c must be 1, otherwise it is invalid
* @internals    无
*/
double MatrixNorm2(int r, int c, const double M[])
{
	int i;
	double val = 0;

	int n = (r == 1 ? c : (c == 1 ? r : 0));

	for (i = 0; i < n; i++)
	{
		val += M[i] * M[i];
	}

	return sqrt(val);
}


/**
* @brief        Matrix Dot
* @param[in]    r       int         The number of rows in the matrix M
* @param[in]    c       int         The number of cols in the matrix M
* @param[in]    M1      double[]    Matrix
* @param[in]    M2      double[]    Matrix
* @return       double  Returns the scalar product of two matrice
* @note         Either r or c must be 1, otherwise it is invalid
* @internals    无
*/
double MatrixDot2(int r, int c, const double M1[], const double M2[])
{
	int i;
	double val = 0;

	int n = (r == 1 ? c : (c == 1 ? r : 0));

	for (i = 0; i < n; i++)
	{
		val += M1[i] * M2[i];
	}

	return val;
}


/**
* @brief        Matrix is Symmetric?
* @param[in]    r       int         The number of rows in the matrix M
* @param[in]    c       int         The number of cols in the matrix M
* @param[in]    M       double[]    M is self multiplied n times
* @return       bool    true: symmetric, false: non-symmetric
* @note         Here M must be a square matrix, i.e. r == c
* @internals    无
*/
bool MatrixIsSymmetric(int r, int c, const double M[])
{
	if (r != c) return false;
	int i, j;
	for (i = 0; i < r; i++)
		for (j = 0; j < c; j++)
		{
			if (fabs(M[i * c + j] - M[j * r + i]) > 1e-8)
			{
				//z = fabs(M[i*c + j] - M[j*r + i]);
				return false;
			}
		}
	return true;
}


/**
* @brief        Matrix is SkewSymmetric?
* @param[in]    r       int         The number of rows in the matrix M
* @param[in]    c       int         The number of cols in the matrix M
* @param[in]    M       double[]    M is self multiplied n times
* @return       bool    true: SkewSymmetric, false: non-SkewSymmetric
* @note         Here M must be a square matrix, i.e. r == c
* @internals    无
*/
bool MatrixIsSkewSymmetric(int r, int c, const double M[])
{
	if (r != c) return false;
	int i, j;
	for (i = 0; i < r; i++)
		for (j = 0; j < c; j++)
		{
			if (fabs(M[i * c + j] + M[j * r + i]) > MATRIX_EPSILON) return false;
		}
	return true;
}


/**
* @brief        Adds two matrix : M3 = s1*M1 + s2*M2
* @param[in]    r     int         The number of rows in the matrix M1, M2 or M3
* @param[in]    c     int         The number of cols in the matrix M1, M2 or M3
* @param[in]    M1    double[]    Src matrix
* @param[in]    M2    double[]    Src matrix
* @param[out]   M3    double[]    Dst matrix
* @param[in]    s1    double      The scale factor of M1
* @param[in]    s2    double      The scale factor of M2
* @return       void
* @note         These matrices must have the same dimensions
* @internals    无
*/
void MatrixAddition1(int r, int c, const double M1[], const double M2[], double M3[], double s1, double s2)
{
	int i, j;

	for (i = 0; i < r; i++)
	{
		for (j = 0; j < c; j++)
		{
			M3[i * c + j] = s1 * M1[i * c + j] + s2 * M2[i * c + j];
		}
	}
}


/**
* @brief          Adds two matrix : M2 = s1*M1 + s2*M2
* @param[in]      r     int         The number of rows in the matrix M1 or M2
* @param[in]      c     int         The number of cols in the matrix M1 or M2
* @param[in]      M1    double[]    Src matrix
* @param[in,out]  M2    double[]    Dst matrix
* @param[in]      s1    double      The scale factor of M1
* @param[in]      s2    double      The scale factor of M2
* @return         void
* @note           These matrices must have the same dimensions
* @internals      无
*/
void MatrixAddition2(int r, int c, const double M1[], double M2[], double s1, double s2)
{
	int i, j;

	for (i = 0; i < r; i++)
	{
		for (j = 0; j < c; j++)
		{
			M2[i * c + j] = s1 * M1[i * c + j] + s2 * M2[i * c + j];
		}
	}
}


/**
* @brief        Subtracts two matrices : M3 = s1*M1 - s2*M2
* @param[in]    r     int         The number of rows in the matrix M1, M2 or M3
* @param[in]    c     int         The number of cols in the matrix M1, M2 or M3
* @param[in]    M1    double[]    Src matrix
* @param[in]    M2    double[]    Src matrix
* @param[out]   M3    double[]    Dst matrix
* @param[in]    s1    double      The scale factor of M1
* @param[in]    s2    double      The scale factor of M2
* @return       void
* @note         These matrices must have the same dimensions
* @internals    无
*/
void MatrixSubtraction1(int r, int c, const double M1[], const double M2[], double M3[], double s1, double s2)
{
	int i, j;

	for (i = 0; i < r; i++)
	{
		for (j = 0; j < c; j++)
		{
			M3[i * c + j] = s1 * M1[i * c + j] - s2 * M2[i * c + j];
		}
	}
}


/**
* @brief          Subtracts two matrix : left:M2 = M1 - M2, right: M2 = s2*M2 - s1*M1
* @param[in]      r        int         The number of rows in the matrix M1 or M2
* @param[in]      c        int         The number of cols in the matrix M1 or M2
* @param[in]      bleft    bool        left minus or right minus
* @param[in]      M1       double[]    Src matrix
* @param[in,out]  M2       double[]    Dst matrix
* @param[in]      s1       double      The scale factor of M1
* @param[in]      s2       double      The scale factor of M2
* @return         void
* @note           These matrices must have the same dimensions
*                 left:M2 = M1 - M2, right: M2 = M2 - M1
* @internals      无
*/
void MatrixSubtraction2(int r, int c, bool bleft, const double M1[], double M2[], double s1, double s2)
{
	int i, j;

	if (bleft)
	{
		for (i = 0; i < r; i++)
		{
			for (j = 0; j < c; j++)
			{
				M2[i * c + j] = s1 * M1[i * c + j] - s2 * M2[i * c + j];
			}
		}
	}
	else
	{
		for (i = 0; i < r; i++)
		{
			for (j = 0; j < c; j++)
			{
				M2[i * c + j] = s2 * M2[i * c + j] - s1 * M1[i * c + j];
			}
		}
	}
}


/**
* @brief        Multiplys two matrices : M3 = scale* (M1 * M2)
* @param[in]    r1       int         The number of rows in the matrix M1
* @param[in]    c1       int         The number of cols in the matrix M1
* @param[in]    M1       double[]    Src matrix
* @param[in]    r2       int         The number of rows in the matrix M2
* @param[in]    c2       int         The number of cols in the matrix M2
* @param[in]    M2       double[]    Src matrix
* @param[out]   M3       double[]    Dst matrix
* @param[in]    scale    double      M3 = scale*M3
* @return       void
* @note         The number of M1 cols and M2 rows must be same
* @internals    无
*/
void MatrixMultiply(int r1, int c1, const double M1[], int r2, int c2, const double M2[], double M3[], double scale)
{
	int i, j, k;
	double Sum;

	for (i = 0; i < r1; i++)
	{
		for (j = 0; j < c2; j++)
		{
			Sum = 0.0;

			for (k = 0; k < c1; k++)
			{
				Sum = Sum + *(M1 + i * c1 + k) * *(M2 + k * c2 + j);
			}

			*(M3 + i * c2 + j) = Sum * scale;
		}
	}
}


/**
* @brief        convolution of two Matrices
* @param[in]    r1,r2,r3    int         The number of rows in the matrix M1, M2 or M3
* @param[in]    M1          double[]    Src matrix
* @param[in]    M2          double[]    Src matrix
* @param[out]   M3          double[]    Dst matrix
* @return       void
* @note         r3 = r1 + r2
* @internals    无
*/
void MatrixConv(int r1, const double M1[], int r2, const double M2[], int r3, double M3[])
{
	const double* signal;
	const double* filter;

	int sigLength = 0, filLength = 0;

	if (r1 > r2)
	{
		sigLength = r1;
		filLength = r2;
		signal = M1;
		filter = M2;
	}
	else
	{
		sigLength = r2;
		filLength = r1;
		signal = M2;
		filter = M1;
	}

	int length = sigLength + filLength - 1;

	for (int i = 1; i <= length; ++i)
	{
		M3[i] = 0;
		if (i < filLength)
			for (int j = 1; j <= i; ++j)
				M3[i] += filter[j] * signal[i - j + 1];
		else if (i <= sigLength)
			for (int j = 1; j <= filLength; ++j)
				M3[i] += filter[j] * signal[i - j + 1];
		else
			for (int j = i - sigLength + 1; j <= filLength; ++j)
				M3[i] += filter[j] * signal[i - j + 1];
	}
}


/**
* @brief        Multiplys two matrices : M3 = a*M1*M2 + b*M3
* @param[in]    tr       char*       "NN" or "NT" or "TN" or "TT"
* @param[in]    r,c,m    int         The number of cols or rows in the matrix M1, M2 or M3
* @param[in]    M1       double[]    Src matrix
* @param[in]    M2       double[]    Src matrix
* @param[out]   M3       double[]    Dst matrix
* @param[in]    a,b      double      The scale factors
* @return       void
* @note         The number of M1 cols and M2 rows must be same
*               NN: M1(r,m) * M2(m,c) = M3(r,c)
*               NT: M1(r,m) * M2(c,m) = M3(r,c)
*               TN: M1(m,r) * M2(m,c) = M3(r,c)
*               TT: M1(m,r) * M2(c,m) = M3(r,c)
* @internals    无
*/
void MatrixMultiplyEx(const char* tr, int r, int c, int m, const double M1[], const double M2[], double M3[], double a, double b)
{
	double d;
	int i, j, k;
	int f = tr[0] == 'N' ? (tr[1] == 'N' ? 1 : 2) : (tr[1] == 'N' ? 3 : 4);

	for (i = 0; i < r; i++) for (j = 0; j < c; j++) {
		d = 0.0;
		switch (f) {
		case 1: for (k = 0; k < m; k++) d += M1[i * m + k] * M2[k * c + j]; break;
		case 2: for (k = 0; k < m; k++) d += M1[i * m + k] * M2[j * m + k]; break;
		case 3: for (k = 0; k < m; k++) d += M1[k * r + i] * M2[k * c + j]; break;
		case 4: for (k = 0; k < m; k++) d += M1[k * r + i] * M2[j * m + k]; break;
		}
		if (b == 0.0) M3[i * c + j] = a * d; else M3[i * c + j] = a * d + b * M3[i * c + j];
	}
}


/**
* @brief        Multiplys three matrices : BTPB = BT * P * B; if P == NULL, BTPB = BT * B
* @param[in]    r             int         The number of rows
* @param[in]    c             int         The number of cols
* @param[in]    B             double[]    Src matrix, r x c
* @param[in]    P             double*     Src matrix, r x r, if there is no P, P == NULL!!
* @param[in]    bPDiagonal    bool        The P is diagonal?
* @param[out]   BTPB          double[]    Dst matrix, c x c
* @return       void
* @note         无
* @internals    无
*/
void MatrixMultiply_BTPB(int r, int c, const double B[], const double* P, bool bPDiagonal, double BTPB[])
{
	int i, j, k;
	double Sum;

	if (!P)
	{
		for (i = 0; i < c; i++)
		{
			for (j = 0; j < c; j++)
			{
				Sum = 0.0;

				for (k = 0; k < r; k++)
				{
					Sum += B[k * c + i] * B[k * c + j];
				}

				BTPB[i * c + j] = Sum;
			}
		}
	}
	else if (bPDiagonal)
	{
		for (i = 0; i < c; i++)
		{
			for (j = 0; j < c; j++)
			{
				Sum = 0.0;

				for (k = 0; k < r; k++)
				{
					Sum += B[k * c + i] * B[k * c + j] * P[k * r + k];
				}

				BTPB[i * c + j] = Sum;
			}
		}
	}
	else
	{
		//double* vec = (double*)malloc(r, sizeof(double));
		double vec[MAXMATRIXSIZE] = { 0.0 };

		for (i = 0; i < c; i++)
		{
			for (j = 0; j < r; j++)
			{
				Sum = 0.0;

				for (k = 0; k < r; k++)
				{
					Sum = Sum + B[i + k * c] * P[k * r + j];
				}

				vec[j] = Sum;
			}

			for (j = 0; j < c; j++)
			{
				Sum = 0.0;

				for (k = 0; k < r; k++)
				{
					Sum = Sum + vec[k] * B[j + k * c];
				}

				BTPB[i * c + j] = Sum;
			}
		}
	}
}


/**
* @brief        Multiplys three matrices : BTPL = BT * P * L; if P == NULL, BTPL = BT * L
* @param[in]    r       int         The number of rows
* @param[in]    c       int         The number of cols
* @param[in]    B       double[]    Src matrix, r x c
* @param[in]    P       double*     Src matrix, r x r, if there is no P, P == NULL!!
* @param[in]    L       double[]    Src matrix, r x 1
* @param[out]   BTPL    double[]    Dst matrix, c x 1
* @return       void
* @note         The number of B rows and L rows must be same
* @internals    无
*/
void MatrixMultiply_BTPL(int r, int c, const double B[], const double* P, const double L[], bool bPDiagonal, double BTPL[])
{
	int i, k, j;
	double Sum;

	if (!P)
	{
		for (i = 0; i < c; i++)
		{
			Sum = 0;

			for (k = 0; k < r; k++)
			{
				Sum += B[k * c + i] * L[k];
			}

			BTPL[i] = Sum;
		}
	}
	else if (bPDiagonal)
	{
		for (i = 0; i < c; i++)
		{
			Sum = 0;

			for (k = 0; k < r; k++)
			{
				Sum += B[k * c + i] * L[k] * P[k * r + k];
			}

			BTPL[i] = Sum;
		}
	}
	else
	{
		//double* vec = (double*)malloc(r, sizeof(double));
		double vec[MAXMATRIXSIZE] = { 0.0 };

		for (i = 0; i < c; i++)
		{
			for (j = 0; j < r; j++)
			{
				Sum = 0.0;

				for (k = 0; k < r; k++)
				{
					Sum = Sum + B[i + k * c] * P[k * r + j];
				}

				vec[j] = Sum;
			}

			Sum = 0.0;

			for (k = 0; k < r; k++)
			{
				Sum = Sum + vec[k] * L[k];
			}

			BTPL[i] = Sum;
		}
	}
}


/**
* @brief        Multiplys three matrices : HPHT = H * P * HT; if P == NULL, HPHT = H * HT
* @param[in]    r             int         The number of rows
* @param[in]    c             int         The number of cols
* @param[in]    H             double[]    Src matrix, r x c
* @param[in]    P             double*     Src matrix, c x c, if there is no P, P == NULL!!
* @param[in]    bPDiagonal    bool        The P is diagonal?
* @param[out]   HPHT          double[]    Dst matrix, r x r
* @return       void
* @note         无
* @internals    无
*/
void MatrixMultiply_HPHT(int r, int c, const double H[], const double* P, bool bPDiagonal, double HPHT[])
{
	int i, j, k;
	double Sum;

	if (!P)
	{
		for (i = 0; i < r; i++)
		{
			for (j = 0; j < r; j++)
			{
				Sum = 0.0;

				for (k = 0; k < c; k++)
				{
					Sum += H[i * c + k] * H[j * c + k];
				}

				HPHT[i * r + j] = Sum;
			}
		}
	}
	else if (bPDiagonal || (!P))
	{
		for (i = 0; i < r; i++)
		{
			for (j = 0; j < r; j++)
			{
				Sum = 0.0;

				for (k = 0; k < c; k++)
				{
					Sum += H[i * c + k] * H[j * c + k] * P[k * c + k];
				}

				HPHT[i * r + j] = Sum;
			}
		}
	}
	else
	{
		//double* vec = (double*)malloc(c, sizeof(double));
		double vec[MAXMATRIXSIZE] = { 0.0 };

		for (i = 0; i < r; i++)
		{
			for (j = 0; j < c; j++)
			{
				Sum = 0.0;

				for (k = 0; k < c; k++)
				{
					Sum = Sum + H[i * c + k] * P[k * c + j];
				}

				vec[j] = Sum;
			}

			for (j = 0; j < r; j++)
			{
				Sum = 0.0;

				for (k = 0; k < c; k++)
				{
					Sum = Sum + vec[k] * H[j * c + k];
				}

				HPHT[i * r + j] = Sum;
			}
		}
	}
}


/**
* @brief          Scale multiplys a matrix : M = scale * M
* @param[in]      scale    double      The scale factor
* @param[in]      r1       int         The number of rows in the matrix M
* @param[in]      c1       int         The number of cols in the matrix M
* @param[in,out]  M        double[]    Src/Dst matrix
* @return         void
* @note           无
* @internals      无
*/
void MatrixScaleMultiply(double scale, int r1, int c1, double M[])
{
	int i, j;

	for (i = 0; i < r1; i++)
	{
		for (j = 0; j < c1; j++)
		{
			M[i * c1 + j] *= scale;
		}
	}
}


/**
* @brief        Transposes matrix, MT = transpose(M)
* @param[in]    r     int         The number of rows in the matrix M or MT
* @param[in]    c     int         The number of cols in the matrix M or MT
* @param[in]    M     double[]    Src matrix
* @param[out]   MT    double[]    Transpose of matrix M
* @return       void
* @note         无
* @internals    无
*/
void MatrixTranspose1(int r, int c, const double M[], double MT[])
{
	int i, j;

	for (i = 0; i < r; i++)
	{
		for (j = 0; j < c; j++)
		{
			MT[j * r + i] = M[i * c + j];
		}
	}
}


/**
* @brief          Transpose of matrix, M = transpose(M)
* @param[in]      r     int         The number of rows in the matrix M
* @param[in]      c     int         The number of cols in the matrix M
* @param[in,out]  M     double[]    Transpose of matrix M
* @return         void
* @note           The original matrix M is replaced by the transpose of matrix
* @internals      无
*/
void MatrixTranspose2(int r, int c, double M[])
{
	//double* MT = (double*)malloc(r * c, sizeof(double));
	double MT[MAXMATRIXSIZE * MAXMATRIXSIZE];
	MatrixTranspose1(r, c, M, MT);
	MatrixCopy(r, c, MT, M);
}


/**
* @brief        Inverse of matrix, M_Inv = Inverse(M)
* @param[in]    r        int         The number of rows in the matrix M or M_Inv
* @param[in]    c        int         The number of cols in the matrix M or M_Inv
* @param[in]    M        double[]    Src matrix
* @param[out]   M_Inv    double[]    Inverse of matrix M
* @return       bool     success or failure
* @note         无
* @internals    无
*/
bool MatrixInv1(int r, int c, double M[], double M_Inv[])
{
	//int* is = (int*)malloc(r, sizeof(int));
	//int* js = (int*)malloc(r, sizeof(int));
	int is[MAXMATRIXSIZE] = { 0 };
	int js[MAXMATRIXSIZE] = { 0 };
	int i, j, k, l, u, v;
	double d, p;

BEGIN:
	for (i = 0; i < r; i++)
	{
		for (j = 0; j < r; j++)
		{
			M_Inv[i * r + j] = M[i * r + j];
		}
	}

	for (k = 0; k < r; k++)
	{
		d = 0.0;
		for (i = k; i < r; i++)
		{
			for (j = k; j < r; j++)
			{
				l = r * i + j;
				p = fabs(M_Inv[l]);
				if (p > d)
				{
					d = p;
					is[k] = i;
					js[k] = j;
				}
			}
		}

		if (d < MATRIX_EPSILON)
		{
			if (d < 1e-30) {
				return false;
			}

			for (i = 0; i < r; i++)   M[i * (r + 1)] += 1E-15; // M矩阵接近奇异时,主对角线加上一个小量
			goto BEGIN;
		}

		if (is[k] != k)
		{
			for (j = 0; j < r; j++)
			{
				u = k * r + j;
				v = is[k] * r + j;
				p = M_Inv[u];
				M_Inv[u] = M_Inv[v];
				M_Inv[v] = p;
			}
		}

		if (js[k] != k)
		{
			for (i = 0; i < r; i++)
			{
				u = i * r + k;
				v = i * r + js[k];
				p = M_Inv[u];
				M_Inv[u] = M_Inv[v];
				M_Inv[v] = p;
			}
		}

		l = k * r + k;
		M_Inv[l] = 1.0 / M_Inv[l];
		for (j = 0; j < r; j++)
		{
			if (j != k)
			{
				u = k * r + j;
				M_Inv[u] = M_Inv[u] * M_Inv[l];
			}
		}
		for (i = 0; i < r; i++)
		{
			if (i != k)
			{
				for (j = 0; j < r; j++)
				{
					if (j != k)
					{
						u = i * r + j;
						M_Inv[u] = M_Inv[u] - M_Inv[i * r + k] * M_Inv[k * r + j];
					}
				}
			}
		}
		for (i = 0; i < r; i++)
		{
			if (i != k)
			{
				u = i * r + k;
				M_Inv[u] = -M_Inv[u] * M_Inv[l];
			}
		}
	}

	for (k = r - 1; k >= 0; k--)
	{
		if (js[k] != k)
		{
			for (j = 0; j < r; j++)
			{
				u = k * r + j;
				v = js[k] * r + j;
				p = M_Inv[u];
				M_Inv[u] = M_Inv[v];
				M_Inv[v] = p;
			}
		}
		if (is[k] != k)
		{
			for (i = 0; i < r; i++)
			{
				u = i * r + k;
				v = is[k] + i * r;
				p = M_Inv[u];
				M_Inv[u] = M_Inv[v];
				M_Inv[v] = p;
			}
		}
	}

	return true;
}


/**
* @brief          Inverse of matrix, M = Inverse(M)
* @param[in]      r     int         The number of rows in the matrix M
* @param[in]      c     int         The number of cols in the matrix M
* @param[in,out]  M     double[]    Inverse of matrix M
* @return         bool  success or failure
* @note           The original matrix M is replaced by the inverse of matrix
* @internals      无
*/
bool MatrixInv2(int r, int c, double M[])
{
	int i;
	double b[MAXMATRIXSIZE * MAXMATRIXSIZE];
	bool flag = MatrixInv1(r, c, M, b);

	r = r * r;
	for (i = 0; i < r; i++) M[i] = b[i];

	return flag;
}


/**
* @brief          Inverse of the real symmetric positive definite matrix
* @param[in]      r     int         The number of rows in the matrix M
* @param[in]      c     int         The number of cols in the matrix M
* @param[in,out]  M     double[]    Inverse of matrix M
* @return         bool  success or failure
* @note           The real symmetric positive definite matrix has the form: M = ATPA. \n
*                 This function for inverse of SP matrix is efficient.
* @internals      无
*/
bool MatrixInvSP(int r, int c, double M[])
{
	int j, k;

	for (k = 0; k < r; k++)
	{
		if (M[k * r + k] <= 0.0) return false;
		if (M[k * r + k] < MATRIX_EPSILON) M[k * r + k] += 1e-15; // M矩阵接近奇异时,主对角线加上一个小量

		for (int i = 0; i < r; i++)
		{
			if (i != k)   M[i * r + k] = -M[i * r + k] / M[k * r + k];
		}
		M[k * r + k] = 1.0 / M[k * r + k];

		for (int i = 0; i < r; i++)
		{
			if (i != k)
			{
				for (j = 0; j < r; j++)
				{
					if (j != k)  M[i * r + j] += M[k * r + j] * M[i * r + k];
				}
			}
		}
		for (j = 0; j < r; j++)
		{
			if (j != k)  M[k * r + j] *= M[k * r + k];
		}
	}

	return true;
}


/**
* @brief        Inverse of the real symmetric positive definite matrix
* @param[in]    index     int*       List of indices to be copied in the src matrix
* @param[in]    nindex    int        The size of index
* @param[in]    src       double*    Src matrix
* @param[out]   dst       double*    Dst matrix
* @param[in]    col       int        The number of cols in the matrix src or dst
* @return       void
* @note         The function only supports a vector (col = 1) or square matrix
* @internals    无
*/
void EQU(int* index, int nindex, const double* src, double* dst, int col)
{
	if (!index || !src || !dst || nindex < 1) return;

	int i, j;
	if (col == 1)
	{
		for (i = 0; i < nindex; i++)
			dst[index[i]] = src[index[i]];
	}
	else
	{
		for (i = 0; i < nindex; i++) for (j = 0; j < nindex; j++)
			dst[index[i] * col + index[j]] = src[index[i] * col + index[j]];
	}
}


/**
* @brief        状态向量与状态方差复制
* @param[in]    X_Src    double*    Src matrix of StateX, row*1
* @param[in]    P_Src    double*    Src matrix of StateP, row*row
* @param[out]   X_Dst    double*    Dst matrix of StateX, row*1
* @param[out]   P_Dst    double*    Dst matrix of StateP, row*row
* @param[in]    row      int        The number of rows in the square matrix src or dst
* @note         The function only serves for GNSS_GLSINFO.StateX/StateXa/MKFStateX and GNSS_GLSINFO.StateP/StatePa/MKFStateP \n
				The purpose is to accelerate the copying process of large matrices
* @par History:
*               2024/10/06, Ying Liu, modify \n
* @internals    无
*/
void EQU_StateXP(const double* X_Src, const double* P_Src, double* X_Dst, double* P_Dst, int row)
{
	if (!X_Src || !P_Src || !X_Dst || !P_Dst)return;

	int StateIndex[MAXMATRIXSIZE], nStateIndex = 0;
	for (int i = 0; i < row; i++)
	{
		if (P_Src[i * row + i] > 0.0) StateIndex[nStateIndex++] = i;
	}
	EQU(StateIndex, nStateIndex, X_Src, X_Dst, row);
	EQU(StateIndex, nStateIndex, P_Src, P_Dst, row);
}


/**
* @brief        状态方差复制
* @param[in]    P_Src    double*    Src matrix of StateP, row*row
* @param[out]   P_Dst    double*    Dst matrix of StateP, row*row
* @param[in]    row      int        The number of rows in the square matrix src or dst
* @note         The function only serves for GNSS_GLSINFO.StateP/StatePa/MKFStateP \n
				The purpose is to accelerate the copying process of large matrices
* @par History:
*               2024/10/06, Ying Liu, modify \n
* @internals    无
*/
void EQU_StateP(const double* P_Src, double* P_Dst, int row)
{
	if (!P_Src || !P_Dst)return;

	int StateIndex[MAXMATRIXSIZE], nStateIndex = 0;
	for (int i = 0; i < row; i++)
	{
		if (P_Src[i * row + i] > 0.0) StateIndex[nStateIndex++] = i;
	}
	EQU(StateIndex, nStateIndex, P_Src, P_Dst, row);
}


/**
* @brief          删除向量或者矩阵的第dr行和第dc列
* @param[in]      r     int         The number of rows in the matrix M
* @param[in]      c     int         The number of cols in the matrix M
* @param[in]      dr    int         Delete row number, starting from 1
* @param[in]      dc    int         Delete col number, starting from 1
* @param[in,out]  M     double[]    Src/dst matrix
* @note           无
* @par History:
*                 2024/10/06, Ying Liu, new \n
* @internals      无
*/
void MatrixDelRC(int r, int c, int dr, int dc, double M[])
{
	if (dr > r || dc > c || dr <= 0 || dc <= 0) return;

	int i = 0, j = 0;
	double MCopy[MAXMATRIXSIZE * MAXMATRIXSIZE] = { 0.0 };
	MatrixCopy(r, c, M, MCopy);
	MatrixZero(r, c, M);

	if (r == 1)
	{
		for (i = 0; i < dc - 1; i++)        M[i] = MCopy[i];
		for (i = dc - 1; i < c - 1; i++)    M[i] = MCopy[i + 1];
	}
	else if (c == 1)
	{
		for (i = 0; i < dr - 1; i++)        M[i] = MCopy[i];
		for (i = dr - 1; i < r - 1; i++)    M[i] = MCopy[i + 1];
	}
	else
	{
		for (i = 0; i < dr - 1; i++)
		{
			for (j = 0; j < dc - 1; j++)    M[i * c + j] = MCopy[i * c + j];
		}
		for (i = 0; i < dr - 1; i++)
		{
			for (j = dc - 1; j < c - 1; j++)M[i * c + j] = MCopy[i * c + j + 1];
		}
		for (i = dr - 1; i < r - 1; i++)
		{
			for (j = 0; j < dc - 1; j++)    M[i * c + j] = MCopy[(i + 1) * c + j];
		}
		for (i = dr - 1; i < r - 1; i++)
		{
			for (j = dc - 1; j < c - 1; j++)M[i * c + j] = MCopy[(i + 1) * c + j + 1];
		}
	}
}


/**
* @brief        LtDL分解
* @param[in]    Qahat    const double*    协方差矩阵，维数为row*row
* @param[out]   L        double*          Qahat = LT * D * L，L为单位下三角矩阵
* @param[out]   D        double*          Qahat = LT * D * L，D为对角阵（输出按row*1的数组输出）
* @param[in]    nfixed   int              模糊度维数
* @return       bool     true = 分解成功
* @note         得到协方差矩阵的条件方差和条件方差的系数阵
* @internals    无
*/
bool LTDLcom(const double Qahat[], double L[], double D[],const int row)
{
	double Qahat_[MAXMATRIXSIZE * MAXMATRIXSIZE] = { 0.0 };
	MatrixCopy(row, row, Qahat, Qahat_);
	double a = 0.0;

	for (int i = row - 1; i >= 0; i--)
	{
		if ((D[i] = Qahat_[i * row + i]) <= 0)
		{
			return false;
		}
		a = sqrt(Qahat_[i * row + i]);
		for (int k = 0; k <= i; k++) L[i * row + k] = Qahat_[i * row + k] / a;
		for (int j = 0; j <= i - 1; j++) for (int k = 0; k <= j; k++) Qahat_[j * row + k] = Qahat_[j * row + k] - L[i * row + k] * L[i * row + j];
		for (int k = 0; k <= i; k++) L[i * row + k] = L[i * row + k] / L[i * row + i];
	}

	for (int i = 0; i < row; i++)
	{
		if (D[i] < 1E-10)
		{
			return false;
		}
	}

	return true;
}


/**
* @brief        模糊度降相关
* @param[in]    Qahat    const double*    浮点模糊度方差，维数为nfixed*nfixed
* @param[in]    ahat     const double*    浮点模糊度，维数为nfixed
* @param[out]   Qzhat    double*          降相关浮点模糊度方差，Qzhat = ZT * Qahat * Z
* @param[out]   Z        double*          降相关矩阵（整数变换矩阵）
* @param[out]   L        double*          Qzhat = LT * D * L，L为单位下三角矩阵
* @param[out]   D        double*          Qzhat = LT * D * L，D为对角阵
* @param[out]   zhat     double*          降相关浮点模糊度，zhat = ZT * ahat
* @param[out]   iZt      double*          逆降相关矩阵，ahat = iZt * zhat
* @param[in]    nfixed   int              模糊度维数
* @return       bool     true = 正常
* @note         无
* @internals    无
*/
bool decorrel(const double Qahat[], const double ahat[], double Qzhat[], double Z[], double L[], double D[], double zhat[], double iZt[], int nfixed)
{
	//--- Initialisations ---
	int row = nfixed;
	int r = 0;
	int il = row - 1;
	int sw = 1;
	int i = 0, j = 0;
	int mu = 0;
	double delta = 0.0;
	double lambda = 0.0;
	double eta = 0.0;
	double as[2 * 2] = { 0.0 };
	double tmp[2 * MAXMATRIXSIZE] = { 0.0 };
	double as_tmp[2 * MAXMATRIXSIZE] = { 0.0 };
	MatrixUnit(row, row, iZt);

	// --- LtDL decomposition ---
	if (LTDLcom(Qahat, L, D, row) == false)
	{
		return false;
	}

	while (sw)
	{
		i = row;
		sw = 0;

		while ((!sw) && (i > 1))
		{
			i = i - 1;

			if (i <= il)
			{
				for (int j = i; j < row; j++)
				{
					mu = (int)(round(L[j * row + i - 1]));
					if (mu)
					{
						for (int k = j; k < row; k++) L[k * row + i - 1] = L[k * row + i - 1] - mu * L[k * row + j];
						for (int k = 0; k < row; k++) iZt[k * row + j] = iZt[k * row + j] + mu * iZt[k * row + i - 1];
					}
				}
			}

			delta = D[i - 1] + (L[i * row + i - 1]) * (L[i * row + i - 1]) * D[i];

			if (delta < D[i])
			{
				lambda = D[i] * L[i * row + i - 1] / delta;
				eta = D[i - 1] / delta;
				D[i - 1] = eta * D[i];
				D[i] = delta;

				as[0] = -L[i * row + i - 1]; as[1] = 1;
				as[2] = eta;                 as[3] = lambda;

				if (i > 1)
				{
					int s1[2] = { i - 1, 0 };
					int s2[2] = { 0, 0 };
					MatrixCopySub(s1, s2, 2, i - 1, row, row, L, 2, i - 1, tmp, false);
					MatrixMultiply(2, 2, as, 2, i - 1, tmp, as_tmp, 1.0);
					int s10[2] = { 0, 0 };
					int s20[2] = { i - 1, 0 };
					MatrixCopySub(s10, s20, 2, i - 1, 2, i - 1, as_tmp, row, row, L, false);
				}

				L[i * row + i - 1] = lambda;

				// i 和 i-1列调换
				for (r = i + 1; r < row; r++)
				{
					lambda = L[r * row + i];
					L[r * row + i] = L[r * row + i - 1];
					L[r * row + i - 1] = lambda;
				}

				for (r = 0; r < row; r++)
				{
					lambda = iZt[r * row + i];
					iZt[r * row + i] = iZt[r * row + i - 1];
					iZt[r * row + i - 1] = lambda;
				}

				il = i;
				sw = 1;
			}
		}
	}

	double ZT[MAXMATRIXSIZE * MAXMATRIXSIZE] = { 0.0 };
	MatrixTranspose1(row, row, iZt, ZT);
	MatrixInv2(row, row, ZT);
	for (i = 0; i < row; i++)
		for (j = 0; j < row; j++)
		{
			Z[row * i + j] = round(ZT[row * i + j]);
		}

	MatrixMultiply_BTPB(row, row, Z, Qahat, false, Qzhat);
	MatrixMultiplyEx("TN", row, 1, row, Z, ahat, zhat, 1.0, 0.0);

	return true;
}




// PCA分解函数：输入协方差矩阵P，输出Q（特征向量）和Lambda（特征值）
void pca_decomposition(int n, double* P, double* Q, double* Lambda) {
	double A[MAXMATRIXSIZE * MAXMATRIXSIZE] = { 0.0 };
	// 初始化工作矩阵和特征向量矩阵
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			A[i * n + j] = P[i * n + j];
			Q[i * n + j] = (i == j) ? 1.0 : 0.0;  // Q初始化为单位矩阵
		}
	}

	// 雅可比迭代
	for (int iter = 0; iter < 1000; iter++) {  // 1000为最大迭代次数
		double max_offdiag = 0.0;
		int p = 0, q = 0;

		// 找到最大的非对角元素
		for (int i = 0; i < n; i++) {
			for (int j = i + 1; j < n; j++) {
				if (fabs(A[i * n + j]) > max_offdiag) {
					max_offdiag = fabs(A[i * n + j]);
					p = i;
					q = j;
				}
			}
		}

		// 收敛检查
		if (max_offdiag < 1e-10) break;

		// 计算旋转角度
		double theta = 0.5 * atan2(2 * A[p * n + q], A[q * n + q] - A[p * n + p]);
		double c = cos(theta);
		double s = sin(theta);

		// 更新矩阵A
		double App = A[p * n + p];
		double Aqq = A[q * n + q];
		double Apq = A[p * n + q];

		A[p * n + p] = c * c * App + s * s * Aqq - 2 * s * c * Apq;
		A[q * n + q] = s * s * App + c * c * Aqq + 2 * s * c * Apq;
		A[p * n + q] = A[q * n + p] = 0.0;

		// 更新其他行和列
		for (int j = 0; j < n; j++) {
			if (j != p && j != q) {
				// 更新第p行和第q行
				double Ajp = A[j * n + p];
				double Ajq = A[j * n + q];
				A[j * n + p] = c * Ajp - s * Ajq;
				A[p * n + j] = A[j * n + p];  // 保持对称性
				A[j * n + q] = s * Ajp + c * Ajq;
				A[q * n + j] = A[j * n + q];
			}
		}

		// 更新特征向量矩阵Q
		for (int j = 0; j < n; j++) {
			double Q_jp = Q[j * n + p];
			double Q_jq = Q[j * n + q];
			Q[j * n + p] = c * Q_jp - s * Q_jq;
			Q[j * n + q] = s * Q_jp + c * Q_jq;
		}
	}

	// 提取对角线元素作为特征值
	for (int i = 0; i < n; i++) {
		Lambda[i] = A[i * n + i];
	}

	// 对特征值和特征向量按降序排序
	for (int i = 0; i < n - 1; i++) {
		int max_idx = i;
		for (int j = i + 1; j < n; j++) {
			if (Lambda[j] > Lambda[max_idx]) {
				max_idx = j;
			}
		}
		if (max_idx != i) {
			// 交换特征值
			swap(&Lambda[i], &Lambda[max_idx]);
			// 交换特征向量列
			for (int k = 0; k < n; k++) {
				swap(&Q[k * n + i], &Q[k * n + max_idx]);
			}
		}
	}
}


// 取A-B=C,即A中有，但是B中没有的元素集合
static int findComplement(const int A[], const int B[], const int size, int C[]) {
	bool status = false;
	int k = 0;
	for (int i = 0; i < size; i++)
	{
		status = true;
		for (int j = 0; j < size; j++)
		{
			if (A[i] == B[j]) { status = false; break; }
		}
		if (status) { C[k] = A[i]; k++; }

	}
	return  k;
}