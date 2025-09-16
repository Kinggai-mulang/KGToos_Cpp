#pragma once
#include<cmath>
#include<iostream>
//#include<string>

constexpr double PI = 3.1415926535897932384626433832795;		///< pi
constexpr double PI_2 = 6.2831853071795864769252867665590;      ///< 2*pi
constexpr double PI_12 = 1.5707963267948966192313216916398;     ///< 0.5*pi
constexpr double D2R = 0.0174532925199432957692369076849;       ///< deg to rad
constexpr double R2D = 57.295779513082322864647721871734;       ///< rad to deg
constexpr double CLIGHT = 299792458.0;                      ///< speed of light (m/s)
constexpr double Epsilon = (2.2204460492503131e-016);			///< minunm

constexpr double SECOND_WEEK = 604800.0;
constexpr double HALF_SECOND_WEEK = SECOND_WEEK / 2.0;




//**************** 坐标转化 *********************/
static constexpr double gs_WGS84_a = 6378137.0;                     ///< earth semimajor axis (WGS84) (m)
static constexpr double gs_WGS84_FE = 1.0 / 298.257223563;             ///< earth flattening (WGS84)
static const double gs_WGS84_OMGE = 7.2921151467E-5;                ///< earth angular velocity (IS-GPS) (rad/s)

static constexpr double gs_CGCS2000_a = 6378137.0;                      ///< earth semimajor axis (WGS84) (m)
static constexpr double gs_CGCS2000_FE = 1.0 / 298.257222101;           ///< earth flattening (WGS84)
static const double gs_CGCS2000_OMGE = 7.2921150E-5;                ///< earth angular velocity (IS-GPS) (rad/s




// 椭球体参数结构体
struct Ellipsoid {
	double a;    // 长半轴
	double f;    // 扁率
	double b;    // 椭球短半轴
	double c;    // 子午线曲率半径
	double e2;   // 第一偏心率平方
	double e2_prime; // 第二偏心率平方

	Ellipsoid(double a, double f)
		: a(a), f(f), b(a* (1 - f)), c(a* a / b), e2(2 * f - f * f), e2_prime(e2 / (1 - e2)) {}
};

// 默认 WGS84 参数  // 示例 CGCS2000 参数
const Ellipsoid WGS84(gs_WGS84_a, gs_WGS84_FE);
const Ellipsoid CGCS2000(gs_CGCS2000_a, gs_CGCS2000_FE);


class CoordinateTransform {
private:
	static Ellipsoid defaultEllipsoid;
public:

	// 静态成员变量，默认使用 WGS84
	static void XYZ2BLH(const double xyz[3], double blh[3]);
	static void BLH2XYZ(const double blh[3], double xyz[3],bool isRad=true);
	static void BLH2ENU(const double blh[3], const double refblh[3], double enu[3]);
	static void XYZ2ENU(const double xyz[3], const double refXyz[3], double enu[3]);
	static void Rxyz2enu(const double LLH[3], double R_E2L[9]);
	static void Vxyz2enu(const double LLH[3], const double v_xyz[3], double v_enu[3],bool isrefBLH=true);

	static double UAngel(const double enu[3]);
	static double EAngel(const double enu[3]);

	static double  Getyita2(const double& B);
	static double  Gett(double B);
	static double  GetM(double B);
	static double  GetN(double B);
	static double  GetMLength(double B);
	static double  GetBf(double X);

	static double getDeltaBLH(const double B,const double dxyz[], double dBLH[]);

};


//*********************** 时间转化 *********************/
/*============================ start const =======================================*/
// 每个月的天数（非闰年）
const int days_in_month[12] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
const static long    ls_JAN61980 = 44244;    ///< gps time reference
const static long    ls_JAN11901 = 15385;
const static double  ls_SECPERDAY = 86400.0;

const static long    ls_MonthDay[2][13] = {
	{0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365},
	{0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366}
};

const static double ls_Leaps[][7] = { ///< leap seconds {y,m,d,h,m,s,utc-gpst,...}
	{2017,1,1,0,0,0,-18},
	{2015,7,1,0,0,0,-17},
	{2012,7,1,0,0,0,-16},
	{2009,1,1,0,0,0,-15},
	{2006,1,1,0,0,0,-14},
	{1999,1,1,0,0,0,-13},
	{1997,7,1,0,0,0,-12},
	{1996,1,1,0,0,0,-11},
	{1994,7,1,0,0,0,-10},
	{1993,7,1,0,0,0, -9},
	{1992,7,1,0,0,0, -8},
	{1991,1,1,0,0,0, -7},
	{1990,1,1,0,0,0, -6},
	{1988,1,1,0,0,0, -5},
	{1985,7,1,0,0,0, -4},
	{1983,7,1,0,0,0, -3},
	{1982,7,1,0,0,0, -2},
	{1981,7,1,0,0,0, -1}
};

//GPS时间
class GPSTIME
{
public:
	int week;
	int secWeek;
	double fracSec;
	GPSTIME()
	{
		week = -1;
		secWeek = 0;
		fracSec = 0.0;
	}
	GPSTIME(const int Week, const int Sec_week, const double Frac)
	{
		week = Week;
		secWeek = Sec_week;
		fracSec = Frac;
	}
	GPSTIME(const int Week, const double sec) {

		week = Week;
		secWeek = (int)(sec);
		fracSec = sec - secWeek;
		if (fracSec < Epsilon)
		{
			fracSec = 0.0;
		}
	}

	void Init()
	{
		week = -1;
		secWeek = 0;
		fracSec = 0.0;
	}

	double getFracSec() const;
	double getSec()const;
	double getSow()const;

	bool isBeyond(const double Otime, const double precision = 0.0) const;
	bool isLess(const double Otime, const double precision = 0.0) const;
	bool isBetween(const double start, const double end, const double precision = 0.0) const;
	void Mutator(const int Iweek,const double Isec);
	bool isSpecify(const double Isec, const double pres = 0.1);

	static GPSTIME adjweek(GPSTIME t, GPSTIME t0);
};
//运算发重载
GPSTIME operator+(GPSTIME& a, GPSTIME& b);
GPSTIME operator-(const GPSTIME& a,const  GPSTIME& b);
GPSTIME operator-(const GPSTIME& a, const double sec);
GPSTIME operator+(const GPSTIME& a, const double sec);
std::ostream& operator<<(std::ostream& out, GPSTIME t);
bool operator==(GPSTIME& a, GPSTIME& b);
bool operator!=(GPSTIME& a, GPSTIME& b);
bool operator<=(GPSTIME& a, GPSTIME& b);
bool operator>=(GPSTIME& a, GPSTIME& b);
bool operator<(GPSTIME& a, GPSTIME& b);
bool operator>(GPSTIME& a, GPSTIME& b);


//公历时间结构体
class UTCTIME
{
public:
	int  Year;
	short Month;
	short Day;
	short Hour;
	short Minute;
	double Second;
	short WeekN; //即周几,0算周日，6是周六，其他对应

	UTCTIME() {
		Year = 0;
		Month = 0;
		Hour = 0;
		Minute = 0;
		Second = 0.0;
		WeekN = 0;
	}
	static short calcuWeekD(int year, int month, int day);
	short getWeekD() const;
	static bool IsLeapYear(int year);
	int getDoy(const UTCTIME& common);
	void Mutator(const int &year, const int& month, const int& day,const int &hour,const int &min,const double &second);

private:
	int DaysSince1980(int year);
};

UTCTIME GPSTtoUTC(const GPSTIME& gps);
GPSTIME UTCtoGPST(const UTCTIME& utc);
GPSTIME BDST2GPST(GPSTIME bdst, bool bSecCorr = false);
GPSTIME TOD2GPST(const int day,const double tod ,const int leapSec=18);


double Distance3D(const double A[], const double B[]);

void TransformNtoB(const double n[3], const double angles[3], double b[3]);
void ins_to_gnss(const double insBLH[3], const double hpr[3], const double lever[3], double newBLH[3]);
void gnss_to_ins(const double gnssBLH[3], const double hpr[3], const double lever[3], double insBLH[3]);