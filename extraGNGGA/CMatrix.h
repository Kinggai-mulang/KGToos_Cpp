#pragma once
#include<vector>
#include <iostream>  //用于print输出，或者<<输出！！！
#include <iomanip>   //主要用于输出小数的控制


constexpr double Epsilon_math = 2.20e-16;
const char MatrixFixPres = 4;  //输出精度设置
#define MATRIX_EPSILON (2.2204460492503131e-016)	///< == DBL_EPSILON

class Matrix {
private:
    size_t m_rows, m_cols;
    std::vector<double> m_data;
public:

    /**      构造函数            **/
    Matrix() : m_rows(0), m_cols(0), m_data() {}                                                  // 无参构造
    Matrix(size_t r, size_t c) : m_rows(r), m_cols(c), m_data(r* c, double{}) {}                  // 双参构造
    Matrix(size_t r, size_t c, std::vector<double>vec);                                     // 使用vec实现构造
    Matrix(size_t r, size_t c, double const a[]) : m_rows(r), m_cols(c), m_data(a, a + r * c) {}  // 三参构造（便于赋值）
    Matrix(const Matrix& a) : m_rows(a.m_rows), m_cols(a.m_cols), m_data(a.m_data) {}                   // 拷贝构造函数，`std::vector` 自带深拷贝
    Matrix(std::initializer_list<std::initializer_list<double>> init);                      // 列表进行初始化1
    Matrix(size_t r, size_t c, std::initializer_list<double> init);                         // 列表初始化2,如果加入默认构造如{}，可以和双参构造合并
    Matrix(const Matrix& a, const size_t srow, const size_t scol, const size_t nrow, const size_t ncol);  //拷贝构造获得矩阵的子矩阵,任意维度
    Matrix(const Matrix& a, const size_t srow, const size_t scol, const size_t ndimensin);  // 拷贝构造获得矩阵的子方阵,固定为方阵
    
    /**      运算符重载          **/
    Matrix& operator=(const Matrix& p);                                                     // 赋值运算符重载，`std::vector` 自动处理深拷贝和内存管理
    Matrix& operator=(std::initializer_list<double> init);                                  // A(3,3)={1,2,3,4,5,6,7,8,9};
    double& at(const size_t i, const size_t j);                                             // 获取元素 (i, j)      // 更方便的取元素操作符重载
    double& operator()(const size_t i, const size_t j);                                     // 非 const 版本，允许修改元素
    const double& operator()(const size_t i, const size_t j) const;                         // const 版本，允许读取元素但不允许修改
    
    //利用友元函数进行封装
    friend Matrix operator*(const Matrix& a, const Matrix& b);              //矩阵乘法
    friend Matrix operator+(const Matrix& a, const Matrix& b);              //矩阵加法
    friend Matrix operator-(const Matrix& a, const Matrix& b);              //减法
    friend std::ostream& operator<<(std::ostream& out, const Matrix& A);    //输出
    friend bool operator==(const Matrix& a, const Matrix& b);               //矩阵相等
   

    //使用模板进行数乘
    template<class G>
    friend Matrix operator*(const Matrix& a, const G& b);  //数乘的重载(数右)
    template<class G>
    friend Matrix operator*(const Matrix& a, const G& b);  //数乘的重载(数左)

    /**      静态类函数          **/  
    static Matrix fill(const int size, const double value);     // 创建一个填充指定值的矩阵
    static Matrix diag(const int size, const double value);     // 创建一个对角矩阵
    static Matrix diag(const int size, const double values[]);  // 使用数组创建对角矩阵
    static Matrix Identity(const int size);                     // 创建单位矩阵,这里的size应该表示维数
    static double Dot(const Matrix& A, const Matrix& B);        // 点乘积;;


    double tr() const;                                          // 矩阵的迹
    int getCols() const;                                        // 矩阵列
    int getRows()const;                                         // 矩阵行
    void resize(const size_t newRows, const size_t newCols, bool reserve = false);   //重新改正矩阵大小
    int size() const;                                           // 返回矩阵的总大小
    void print(int size= MatrixFixPres) const;                                         // 打印矩阵
    void zero();                                                // 将矩阵所有元素设置为0
    void fill(const double value);                              // 将矩阵所有元素设置为指定值
    bool isZero();                                              // 检测是否为0矩阵

    Matrix inv();
    bool LDLTcom(Matrix &L,Matrix &D)const;
    // 求逆
    const Matrix Tra() const;                                   // 转置
    bool isInvertible() const;                                  // 是否可逆
    Matrix inv_sp() const;                                      // 正定矩阵求逆
    Matrix getCofactor(int excludeRow, int excludeCol) const;   // 获得余子式
    double determinant() const;                                 // 计算行列式
    static double getMaDist(const Matrix P);           /******     尚未满意！！！       *****//// 获得马氏距离 vTPv


    //功能性
    Matrix getSubMatrix(size_t startRow, size_t startCol, size_t numRows, size_t numCols) const;  // 取出对应行列子矩阵
    void setSubMatrix(size_t startRow, size_t startCol, const Matrix& subMatrix);                 // 将子矩阵赋值给当前矩阵
    Matrix exp(int n) const;                                                                      // exp是展开成指数
    Matrix symm() const;                                                                          // 矩阵对角化
    bool tr_positive()const;                                                                      // 判断是否对角线元素都大于0
    const bool isNULL()const;                                                                                // 判断是否为空矩阵

    bool toArr(double Arr[], const size_t rows, const size_t cols)const;                               // 转为数组类型
    bool rowtoArr(double Arr[], const size_t row);
    bool coltoArr(double Arr[], const size_t row);
    bool Arr2data(const double Arr[], const int begin = -1, const int size = -1);
    bool Arr2data(const int Arr[], const int begin = -1, const int size = -1);
private:
    void luDecompose(Matrix& L, Matrix& U)const;    // LU分解,用于分解矩阵，实现一般函数的矩阵求逆
    bool choleskyDecompose(Matrix& L) const;        // cholesky分解，一般用于正定矩阵求逆的中间过程
};

static const Matrix NULLMATRIX(0, 0);  //设置空矩阵
Matrix LeastSqure(const Matrix& B, const Matrix& L, const Matrix* P = nullptr);  // 最小二乘


//模板类函数需要在头文件中声明并实现
template<class G>  //数乘的重载(数左),注意数乘时候，两个类型需要相等，比如10*C，C如果是double则改为10.0*C
Matrix operator*(const G& a, const Matrix& b) {

    Matrix result(b.m_rows, b.m_cols);

    for (size_t i = 0; i < b.m_rows * b.m_cols; ++i) {
        result.m_data[i] = b.m_data[i] * a;
    }
    return result;
}
template<class G>  //数乘的重载(数右)
Matrix operator*(const Matrix& a, const G& b) {

    Matrix result(a.m_rows, a.m_cols);

    for (size_t i = 0; i < a.m_rows * a.m_cols; ++i) {
        result.m_data[i] = a.m_data[i] * b;
    }
    return result;
}

template<class T>
inline void swap(T* a, T* b)
{
    T temp = *a;
    *a = *b;
    *b = temp;
}



///* 这里数组类实现数组算法相关运算*/
//class ArrAlgo
//{
//public:
//    // 赋值函数
//    static void M31EQU(const double M[3], double M2[3]) {
//        for (int i = 0; i < 3; i++) M2[i] = M[i];
//    }
//    // 数组三维向量加减法
//    static void M31_M31(const double M1[3], const double M2[3], double M3[3],bool isAdd=false) {
//        if(!isAdd)for (int i = 0; i < 3; i++) M3[i] = M1[i] - M2[i];
//        else for (int i = 0; i < 3; i++) M3[i] = M1[i] + M2[i];
//    }
//    // 数组三维向量扩大
//    static void M31Scale(double a, double M[3]) {
//        for (int i = 0; i < 3; i++) M[i] *= a;
//    }
//    static void MN1EQU(const double M[], double M2[],const int len);
//    static void MN1Scale(const double a, double M[],const int length);
//    static bool MN1Reci(const double M[], double M2[], const int length);
//    static void MN1Scale(const double a, const double M[], double M2[], const int length);
//
//    static double L2d(const double A[], const double B[]);
//    static double L3d(const double A[], const double B[]);
//    static double Md(const double v[], const double P[], const int n);
//
//    static double mean(const double A[], const int length);
//    static double std(const double A[], const int length);
//
//    template<class T>
//    static typename std::enable_if<std::is_arithmetic<T>::value, T>::type
//    SumSquare(const T A[], const int arrlen);
//
//    static double Dot(const double v1[], const double v2[],const int n=3);
//    
//    static int countIntArrEQU(const int arr[],const int arr2[], size_t n);
//};
//
//
//template<class T>
//typename std::enable_if<std::is_arithmetic<T>::value, T>::type
//ArrAlgo::SumSquare(const T A[], const int arrlen)
//{
//    T sum = 0;
//    for (size_t i = 0; i < arrlen; i++)
//    {
//        sum += A[i]*A[i];
//    }
//    return sum;
//}











// 一维矩阵运算
extern void SetValsD(double* dst, double val, int n);
extern void SetValsI(int* dst, int val, int n);
extern void SetValsB(bool* dst, bool val, int n);
extern void SetValsSC(signed char* dst, signed char val, int n);
extern void SetValsUC(unsigned char* dst, unsigned char val, int n);
extern void M31EQU(const double M[3], double M2[3]);
extern void M31EQUM(const double M[3], double M2[3]);
extern void M41EQU(const double M[4], double M2[4]);
extern void M33EQU(const double M[9], double M2[9]);
extern void VecEQUD(const int n, const double V[], double V2[]);
extern void VecEQUI(const int n, const int V[], int V2[]);
extern void VecEQUB(const int n, const bool V[], bool V2[]);
extern void M31Zero(double M[3]);
extern void M33Zero(double M[9]);
extern void M31Scale(double a, double M[3]);
extern void M33Scale(double a, double M[9]);
extern void M31M311(const double M1[3], const double M2[3], double M3[3], double s1, double s2);
extern void M31M312(const double M1[3], double M2[3], double s1, double s2);
extern void M33M331(const double M1[9], const double M2[9], double M3[9], double s1, double s2);
extern void M33M332(const double M1[9], double M2[9], double s1, double s2);
extern void M31_M31(const double M1[3], const double M2[3], double M3[3]);
extern void M33XM31(const double M33[9], const double M31[3], double M31_[3]);
extern void M33XM33(const double M33_1[9], const double M33_2[9], double M33_3[9]);
extern void CrossM3(const double v1[3], const double v2[3], double v3[3]);
extern double MatrixDot1(int n, const double v1[], const double v2[]);
extern void MatrixSkewSymmetric(const double v[3], double M[9]);
extern void MatrixSkewSymmetric_T(const double M[9], double v[3]);
extern void MatrixSymmetrization(int r, int c, double M[]);
extern void MatrixCopy(int r, int c, const double* src, double* dst);
extern void MatrixCopySub(const int s1[2], const int s2[2], int nRow, int nCol, int r1, int c1, const double* M1, int r2, int c2, double* M2, bool bMinus);
extern void MatrixSwap(int r, int c, double* src, double* dst);
extern void MatrixZero(int r, int c, double M[]);
extern void MatrixZero_diag(int r, int c, double M[]);
extern void MatrixUnit(int r, int c, double M[]);
extern double MatrixTrace(int r, int c, const double M[]);
extern double MatrixNorm2(int r, int c, const double M[]);
double MatrixDot2(int r, int c, const double M1[], const double M2[]);
extern bool MatrixIsSymmetric(int r, int c, const double M[]);
extern bool MatrixIsSkewSymmetric(int r, int c, const double M[]);
extern void MatrixAddition1(int r, int c, const double M1[], const double M2[], double M3[], double s1, double s2);
extern void MatrixAddition2(int r, int c, const double M1[], double M2[], double s1, double s2);
extern void MatrixSubtraction1(int r, int c, const double M1[], const double M2[], double M3[], double s1, double s2);
extern void MatrixSubtraction2(int r, int c, bool bleft, const double M1[], double M2[], double s1, double s2);
extern void MatrixMultiply(int r1, int c1, const double M1[], int r2, int c2, const double M2[], double M3[], double scale=1.);
extern void MatrixConv(int r1, const double M1[], int r2, const double M2[], int r3, double M3[]);
extern void MatrixMultiplyEx(const char* tr, int r, int c, int m, const double M1[], const double M2[], double M3[], double a, double b);
extern void MatrixMultiply_BTPB(int r, int c, const double B[], const double* P, bool bPDiagonal, double BTPB[]);
extern void MatrixMultiply_BTPL(int r, int c, const double B[], const double* P, const double L[], bool bPDiagonal, double BTPL[]);
extern void MatrixMultiply_HPHT(int r, int c, const double H[], const double* P, bool bPDiagonal, double HPHT[]);
extern void MatrixScaleMultiply(double scale, int r1, int c1, double M3[]);
extern void MatrixTranspose1(int r, int c, const double M1[], double MT[]);
extern void MatrixTranspose2(int r, int c, double M[]);
extern bool MatrixInv1(int r, int c, double M[], double M_Inv[]);
extern bool MatrixInv2(int r, int c, double M[]);
extern bool MatrixInvSP(int r, int c, double M[]);
extern void EQU(int* index, int nindex, const double* src, double* dst, int col);
extern void EQU_StateXP(const double* X_Src, const double* P_Src, double* X_Dst, double* P_Dst, int row);
extern void EQU_StateP(const double* P_Src, double* P_Dst, int row);
extern void MatrixDelRC(int r, int c, int dr, int dc, double M[]);
bool LTDLcom(const double Qahat[], double L[], double D[], const int row);
extern bool decorrel(const double Qahat[], const double ahat[], double Qzhat[], double Z[], double L[], double D[], double zhat[], double iZt[], int nfixed);

void pca_decomposition(int n, double* P, double* Q, double* Lambda);

/**************************************** 其他自建函数************************************/

extern void PrintMatrix(const double* data, const int rows, const int cols);
extern bool PrintMatrixToFile(const double* data, const int rows, const int cols, char mode='w', const char precision=3);
extern bool isZero(const double x, const double epsilon = Epsilon_math);
extern bool getFrac(const double x);


