#include "BaseGT.h"     //基础时间、地理库


// 在类外定义静态成员变量，初始化为 WGS84
Ellipsoid CoordinateTransform::defaultEllipsoid = WGS84;   //必须在CPP文件中定义

/**
 * @file         "BaseGT.cpp"
 * @brief        实现对操作符号+号的重载
 * @param[in]    a
 * @param[in]    b
 * @return       GPSTIME
 * @note         
 * @par History
 *               2024/12/4, Ken Ga  new \n
 */
GPSTIME operator+(GPSTIME& a, GPSTIME& b)
{
	GPSTIME temp;
	temp.week = a.week + b.week;
	temp.secWeek = a.secWeek + b.secWeek;
	temp.fracSec = a.fracSec + b.fracSec;

	// 处理秒的小数部分的进位
	if (temp.fracSec >= 1.0)
	{
		temp.fracSec -= 1.0;
		temp.secWeek += 1.0;
	}

	// 处理周内秒的进位 (一周有604800秒)
	if (temp.secWeek >= SECOND_WEEK)
	{
		temp.secWeek -= SECOND_WEEK;
		temp.week += 1;
	}

	return temp;
}


/**
 * @file         "BaseGT.cpp"
 * @brief        运算符号重载，实现GPS时相减
 * @param[in]    a
 * @param[in]    b
 * @return       GPSTIME
 * @note         
 * @par History
 *               2024/12/4, Ken Ga  new \n
 */
GPSTIME operator-(const GPSTIME& a,const GPSTIME& b)
{
	GPSTIME temp;
	temp.week = a.week - b.week;
	temp.secWeek = a.secWeek - b.secWeek;
	temp.fracSec = a.fracSec - b.fracSec;

	// 处理秒的小数部分的借位
	if (temp.fracSec < 0.0)
	{
		temp.fracSec += 1.0;
		temp.secWeek -= 1.0;
	}

	// 处理周内秒的借位
	if (temp.secWeek < 0.0)
	{
		temp.secWeek += SECOND_WEEK;
		temp.week -= 1;
	}

	return temp;
}

/**
 * @file         "BaseGT.cpp"
 * @brief        GPS时间-double秒数
 * @param[in]    gt1
 * @param[in]    sec
 * @return       GPSTIME
 * @note         
 * @par History
 *               2024/12/11, Ken Ga  new \n
 */
GPSTIME operator-(const GPSTIME& gt1, const double sec)
{
	GPSTIME gt2;
	double t = gt1.secWeek + gt1.fracSec - sec;

	long   lt = (int)t;

	double dt = (gt1.secWeek - lt) + gt1.fracSec - sec;

	if (t >= 0 && t < 604800)
	{
		gt2.week = gt1.week;
		gt2.secWeek = lt;
		gt2.fracSec = dt;
	}
	else if (t < 0)
	{
		int i = (int)(-lt / 604800) + 1;
		gt2.week = gt1.week - i;
		if (dt == 0)
		{
			gt2.secWeek = i * 604800 + lt;
			gt2.fracSec = dt;
		}
		else
		{
			gt2.secWeek = i * 604800 + lt - 1;
			gt2.fracSec = dt + 1;
		}
	}
	else if (t >= 604800)
	{
		int i = (int)(lt / 604800);
		gt2.week = gt1.week + i;
		gt2.secWeek = lt - i * 604800;
		gt2.fracSec = dt;
	}

	if (gt2.fracSec < Epsilon)
	{
		gt2.fracSec = 0.0;
	}

	return gt2;
}

GPSTIME operator+(const GPSTIME& gt1, const double sec)
{
	return (gt1 - -sec);
}



std::ostream& operator<<(std::ostream& out, GPSTIME t)  //重载<<，使得能直接输出GPSTIME，为了保证是链式编程，使用返回类型修改
{
	out << t.week << "  " << t.secWeek + t.fracSec;
	return out;
}

bool operator==(GPSTIME& a, GPSTIME& b)  //小数部分小于1e-14待定
{
	if (a.week != b.week || a.secWeek != b.secWeek)return false;

	return (((a.week - b.week) * 604800 + a.secWeek + a.fracSec - b.secWeek + b.fracSec) < Epsilon);
}

bool operator!=(GPSTIME& a, GPSTIME& b)  //小数部分小于1e-14待定
{
	return !(a == b);
}

bool operator>(GPSTIME& a, GPSTIME& b)
{
	return ((a.week - b.week) * 604800 + a.secWeek + a.fracSec - b.secWeek + b.fracSec) > Epsilon;
}

bool operator<(GPSTIME& a, GPSTIME& b)
{
	return ((a.week - b.week) * 604800 + a.secWeek + a.fracSec - b.secWeek + b.fracSec) < Epsilon;
}

bool operator>=(GPSTIME& a, GPSTIME& b)
{
	return !(a < b);
}

bool operator<=(GPSTIME& a, GPSTIME& b)
{
	return !(a > b);
}


UTCTIME GPSTtoUTC(const GPSTIME& gps)
{
	UTCTIME utc_t;
	// 1. GPS时间起始日期是 1980年1月6日
	int base_year = 1980;
	int base_month = 1;
	int base_day = 6;

	// 2. 计算从GPS起点到当前的总天数和秒数
	int total_days = gps.week * 7;   // 周数转化为天数
	double total_seconds = gps.secWeek + gps.fracSec; // 周内的总秒数

	// 3. 将秒数分离为时、分、秒
	int hours = static_cast<int>(total_seconds) / 3600;
	total_seconds -= hours * 3600;
	int minutes = static_cast<int>(total_seconds) / 60;
	double seconds = total_seconds - minutes * 60;

	// 4. 处理天数的累加
	total_days += hours / 24;        // 将小时数中的整天部分加到天数上
	hours %= 24;                     // 剩下的小时数

	// 5. 从1980年1月6日起，累加这些天数，找到对应的年、月、日
	int year = base_year;
	int month = base_month;
	int day = base_day + total_days;

	// 6. 调整日期到合理范围
	while (day > (UTCTIME::IsLeapYear(year) ? (month == 2 ? 29 : days_in_month[month - 1]) : days_in_month[month - 1]))
	{
		day -= (UTCTIME::IsLeapYear(year) ? (month == 2 ? 29 : days_in_month[month - 1]) : days_in_month[month - 1]);
		month++;
		if (month > 12)
		{
			month = 1;
			year++;
		}
	}

	// 7. 计算周几 (0 是周日, 1-6 对应周一到周六)
	int total_weeks = gps.week;
	utc_t.WeekN = (total_weeks * 7 + static_cast<int>(gps.secWeek) / 86400) % 7;

	// 8. 填充 COMMONTIME 结构体
	utc_t.Year = year;
	utc_t.Month = month;
	utc_t.Day = day;
	utc_t.Hour = hours;
	utc_t.Minute = minutes;
	utc_t.Second = seconds;

	return utc_t;
}

// 补充跳秒 18s
GPSTIME UTCtoGPST(const UTCTIME& utc)
{
	int year = utc.Year;
	int month = utc.Month;
	int mday = utc.Day;
	int hour = utc.Hour;
	int minute = utc.Minute;
	double second = utc.Second;

	long   leap = (year % 4 == 0);
	long   yday = ls_MonthDay[leap][month - 1] + mday;
	long   mjd = ((year - 1901) / 4) * 1461 + ((year - 1901) % 4) * 365 + yday - 1 + ls_JAN11901;
	int    gps_week = (mjd - ls_JAN61980) / 7;
	double fmjd = second + minute * 60 + hour * 3600;
	double sec_of_week = ((mjd - ls_JAN61980) - gps_week * 7) * ls_SECPERDAY + fmjd + 18;  //18

	GPSTIME ot(gps_week, sec_of_week);
	return ot;
}

/*
* @brief        BDS时转GPS时,默认只进行14s改正
* @param[in]    bdst        GPSTIME    BDS时间
* @param[in]    bSecCorr    bool       true = 只进行14s改正, false = BDS周改正
* @return       GPSTIME     GPS时间
* @note         无
* @internals    无
*/
GPSTIME BDST2GPST(GPSTIME bdst, bool bSecCorr)
{
	if (bSecCorr == false) bdst.week += 1356;
	return (bdst + 14.0);
}

/**
 * @file         "BaseGT.cpp"
 * @brief        天内秒转周内秒
 * @param[in]    day
 * @param[in]    leapSec
 * @return       GPSTIME
 * @note         
 * @par History
 *               2025/07/2, Ken Ga  new \n
 */
GPSTIME TOD2GPST(const int day, const double tod ,const int leapSec)
{
	return GPSTIME(0, day * 86400.0 + tod + leapSec);
}

double Distance3D(const double A[], const double B[])
{
	double sum = 0.0;
	for (size_t i = 0; i < 3; i++)
	{
		sum += (A[i] - B[i]) * (A[i] - B[i]);
	}
	return sqrt(sum);
}



///**
// * @file         "BaseGT.cpp"
// * @brief        输入字符串，转为double类型的数组
// * @param[in]    input
// * @param[in]    data
// * @param[in]    symbol
// * @return       void
// * @note         
// * @par History
// *               2025/02/5, Ken Ga  new \n
// */
//int String2Dptr(std::string& input, double* data, std::string symbol)
//{
//	// 使用字符串流分割字符串
//	std::istringstream iss(input);
//	std::string token;
//	int index = 0;
//
//	while (std::getline(iss, token, ',')) {
//		// 将每个子字符串转换为 double
//		data[index++] = std::stod(token);
//	}
//	return index;
//}
//
//std::string D2str(double X, const int pres)
//{
//	// 创建一个字符串流
//	std::ostringstream oss;
//	// 设置固定的小数位数
//	oss << std::fixed << std::setprecision(pres) << X;
//	// 返回流中的字符串
//	return oss.str();
//}


/**
 * @file         "BaseGT.cpp"
 * @brief        获得小数点后秒数
 * @return       double
 * @note         
 * @par History
 *               2024/12/11, Ken Ga  new \n
 */
double GPSTIME::getFracSec() const
{
	return this->fracSec;
}

/**
 * @file         "BaseGT.cpp"
 * @brief        获得自1986年1月6日到现在的所有秒数
 * @return       double
 * @note         获得所有秒数 Week*604800+sec
 * @par History
 *               2024/12/11, Ken Ga  new \n
 */
double GPSTIME::getSec() const
{
	return week * SECOND_WEEK + secWeek + fracSec;
}

/**
 * @file         "BaseGT.cpp"
 * @brief        获得周内秒
 * @return       double
 * @note         获得卫星周内秒
 * @par History
 *               2024/12/11, Ken Ga  new \n
 */
double GPSTIME::getSow() const
{
	return secWeek + fracSec;
}

/**
 * @file         "BaseGT.cpp"
 * @brief        超过目标时间多少秒
 * @param[in]    Otime
 * @param[in]    precision
 * @return       bool
 * @note         
 * @par History
 *               2024/12/4, Ken Ga  new \n
 */
bool GPSTIME::isBeyond(const double Otime, const double precision) const
{
	return (this->getFracSec() - Otime) > precision;
}

/**
 * @file         "BaseGT.cpp"
 * @brief        少于目标时间多少秒
 * @param[in]    Otime
 * @param[in]    precision
 * @return       bool
 * @note         
 * @par History
 *               2024/12/4, Ken Ga  new \n
 */
bool GPSTIME::isLess(const double Otime, const double precision) const
{
	return (this->getFracSec() - Otime) < precision;
}

/**
 * @file         "BaseGT.cpp"
 * @brief        在目标时间多少秒以内
 * @param[in]    start
 * @param[in]    end
 * @param[in]    precision
 * @return       bool
 * @note         
 * @par History
 *               2024/12/4, Ken Ga  new \n
 */
bool GPSTIME::isBetween(const double start, const double end, const double precision) const
{
	return ((this->getSow() - start) > precision) && ((this->getSow() - end) < precision);
}

/**
 * @file         "BaseGT.cpp"
 * @brief        赋值函数
 * @param[in]    week
 * @param[in]    sec
 * @return       GPSTIME
 * @note         
 * @par History
 *               2024/12/11, Ken Ga  new \n
 */
void GPSTIME::Mutator(const int Iweek, const double Isec)
{
	week = Iweek;
	secWeek = (int)(Isec);
	fracSec = Isec - secWeek;
	if (fracSec < Epsilon)
	{
		fracSec = 0.0;
	}
}

/**
 * @file         "BaseGT.cpp"
 * @brief        选择是否为特点的时间
 * @param[in]    Isec	定位的时间(周内秒)
 * @param[in]    pres	定位的精度(在该精度下的时间会被定位)
 * @return       bool
 * @note         选择是否为特点的时间
 * @par History
 *               2025/01/16, Ken Ga  new \n
 */
bool GPSTIME::isSpecify(const double Isec, const double pres)
{
	return (fabs(secWeek + fracSec - Isec)) < pres;
}

/**
* @brief        调整时间,和参考时间差不超过半周
* @param[in]    t        GPSTIME    待调整时间
* @param[in]    t0       GPSTIME    参考时间
* @return       GPSTIME  调整后的时间
* @note         无
* @internals    无
*/
GPSTIME GPSTIME::adjweek(GPSTIME t, GPSTIME t0)
{
	double tt = (t - t0).getSec();
	if (tt < -HALF_SECOND_WEEK) return (t + SECOND_WEEK);
	if (tt > HALF_SECOND_WEEK) return (t - SECOND_WEEK);
	return t;
}


/*****************************************************
坐标转化模块
*****************************************************/
/**
* @brief       WGS84坐标转换到CGCS2000坐标,坐标对应的是框架历元
* @param[in]   XYZ_WGS84         double  WGS84坐标(m)
* @param[out]  XYZ_CGCS2000      double  CGCS2000坐标(m)
* @return      void
* @note        WGS84对准到ITRF2000框架,CGCS2000对准到ITRF1997框架,无站速坐标框架转换
* @par History:
*              2018/01/17,Feng Zhu, new \n
* @internals
*/
void WGS84_CGCS2000(const double XYZ_WGS84[3], double XYZ_CGCS2000[3])
{
	double T[3] = { 0.0 };
	T[0] = 0.0067;
	T[1] = 0.0061;
	T[2] = -0.0185;

	const double D = 1.55 * 1E-9;

	// 采用简化写法
	for (int i = 0; i < 3; i++)
	{
		XYZ_CGCS2000[i] = (1.0 + D) * XYZ_WGS84[i] + T[i];
	}
}


/**
	* @brief       CGCS2000坐标转换到WGS84坐标,坐标对应的是框架历元
	* @param[in]   XYZ_CGCS2000      double  CGCS2000坐标(m)
	* @param[out]  XYZ_WGS84         double  WGS84坐标(m)
	* @return      void
	* @note        WGS84对准到ITRF2000框架,CGCS2000对准到ITRF1997框架,无站速坐标框架转换
	* @par History:
	*              2018/01/17,Feng Zhu, new \n
	* @internals
	*/
static void CGCS2000_WGS84(const double XYZ_CGCS2000[3], double XYZ_WGS84[3])
{
	double T[3] = { 0.0 };
	T[0] = 0.0067;
	T[1] = 0.0061;
	T[2] = -0.0185;

	double D = 1.55 * 1E-9;
	// 采用简化写法
	for (int i = 0; i < 3; i++)
	{
		XYZ_WGS84[i] = (XYZ_CGCS2000[i] - T[i]) / (1.0 + D);
	}
}



/**
* @brief        WGS84坐标系: XYZ坐标转LLH坐标
* @param[in]    XYZ        double    XYZ(m) 坐标
* @param[in]    CoorSys    int       坐标系统: 0:WGS84, 1:CGCS2000
* @param[out]   LLH        double    弧度
* @return       void
* @note         无
* @internals    无
*/
void CoordinateTransform::XYZ2BLH(const double xyz[3], double LLH[3]) {
	double a = defaultEllipsoid.a;
	double e2 = defaultEllipsoid.e2;


	double X = xyz[0];
	double Y = xyz[1];
	double Z = xyz[2];

	double r2 = X * X + Y * Y;
	double z = 0.0;
	double zk = 0.0;
	double v = a;
	double sinp = 0.0;

	for (z = Z, zk = 0.0; fabs(z - zk) >= 1E-4;)
	{
		zk = z;
		sinp = z / sqrt(r2 + z * z);
		v = a / sqrt(1.0 - e2 * sinp * sinp);
		z = Z + v * e2 * sinp;
	}
	LLH[0] = r2 > 1E-12 ? atan(z / sqrt(r2)) : (Z > 0.0 ? PI / 2.0 : -PI / 2.0);
	LLH[1] = r2 > 1E-12 ? atan2(Y, X) : 0.0;
	LLH[2] = sqrt(r2 + z * z) - v;
}

/**
* @brief        WGS84坐标系: LLH坐标转XYZ坐标
* @param[in]    LLH        double    Lat[-90,90],Lon[-180,180],Hgt (deg,m),输入应为弧度
* @param[in]    CoorSys    int       坐标系统: 0:WGS84, 1:CGCS2000
* @param[out]   XYZ        double    XYZ(m) 坐标
* @return       void
* @note         输入应为弧度
* @internals    无
*/
void  CoordinateTransform::BLH2XYZ(const double blh[3], double xyz[3], bool isRad) {
	double lat, lon, h;
	if (isRad) {
		 lat = blh[0];
		 lon = blh[1];
		 h = blh[2];
	}
	else
	{
		lat = blh[0]*D2R;
		lon = blh[1]*D2R;
		h = blh[2];
	}

	double N = defaultEllipsoid.a / sqrt(1 - defaultEllipsoid.e2 * sin(lat) * sin(lat));
	xyz[0] = (N + h) * cos(lat) * cos(lon);
	xyz[1] = (N + h) * cos(lat) * sin(lon);
	xyz[2] = (N * (1 - defaultEllipsoid.e2) + h) * sin(lat);

}

// BLH 到 ENU 转换
void CoordinateTransform::BLH2ENU(const double blh[3], const double refblh[3], double enu[3]) {
	double refxyz[3] = { 0.0 }, xyz[3] = { 0.0 };
	BLH2XYZ(refblh, refxyz);
	BLH2XYZ(blh, xyz);

	double dx = xyz[0] - refxyz[0];
	double dy = xyz[1] - refxyz[1];
	double dz = xyz[2] - refxyz[2];

	double sinLat = sin(refblh[0]);
	double cosLat = cos(refblh[0]);
	double sinLon = sin(refblh[1]);
	double cosLon = cos(refblh[1]);

	enu[0] = -sinLon * dx + cosLon * dy;
	enu[1] = -sinLat * cosLon * dx - sinLat * sinLon * dy + cosLat * dz;
	enu[2] = cosLat * cosLon * dx + cosLat * sinLon * dy + sinLat * dz;

}

// XYZ 到 ENU 转换
void  CoordinateTransform::XYZ2ENU(const double xyz[3], const double refxyz[3], double enu[3]) {
	double refblh[3] = { 0.0 };
	XYZ2BLH(refxyz, refblh);
	double dx = xyz[0] - refxyz[0];
	double dy = xyz[1] - refxyz[1];
	double dz = xyz[3] - refxyz[3];

	double sinLat = sin(refblh[0]);
	double cosLat = cos(refblh[0]);
	double sinLon = sin(refblh[1]);
	double cosLon = cos(refblh[1]);

	enu[0] = -sinLon * dx + cosLon * dy;
	enu[1] = -sinLat * cosLon * dx - sinLat * sinLon * dy + cosLat * dz;
	enu[2] = cosLat * cosLon * dx + cosLat * sinLon * dy + sinLat * dz;
}

/**
* @brief        ENU转向ECEF系的旋转矩阵
* @param[in]    LLH      double    纬度,经度,大地高
* @param[out]   R_L2E    double    ENU->ECEF的旋转矩阵
* @return       void
* @note			这里是将参考坐标的BLH传入，获得旋转矩阵
* @internals    无
*/
void CoordinateTransform::Rxyz2enu(const double LLH[3], double R_E2L[9])
{
	double sinp = sin(LLH[0]), cosp = cos(LLH[0]), sinl = sin(LLH[1]), cosl = cos(LLH[1]);
	R_E2L[0] = -sinl;        R_E2L[1] = cosl;         R_E2L[2] = 0.0;
	R_E2L[3] = -sinp * cosl; R_E2L[4] = -sinp * sinl; R_E2L[5] = cosp;
	R_E2L[6] = cosp * cosl;  R_E2L[7] = cosp * sinl;  R_E2L[8] = sinp;
}

/**
* @brief        ECEF的向量转到ENU系下
* @param[in]    LLH      double    纬度,经度,大地高,弧度制
* @param[in]    v_xyz    double    ECEF系下向量
* @param[out]   v_enu    double    ENU系下向量
* @return		将ECEF坐标转为ENU，LLH为参考站的坐标,bool为true时为BLH，为false时为XYZ
* @note         无
* @internals    无
*/
void CoordinateTransform::Vxyz2enu(const double LLH[3], const double v_xyz[3], double v_enu[3], bool isrefBLH)
{
	double R[9] = { 0.0 };
	if (!isrefBLH) {
		double blh[3];
		XYZ2BLH(LLH, blh);
		Rxyz2enu(blh, R);
	}
	else Rxyz2enu(LLH, R);
	v_enu[0] = R[0] * v_xyz[0] + R[1] * v_xyz[1] + R[2] * v_xyz[2];
	v_enu[1] = R[3] * v_xyz[0] + R[4] * v_xyz[1] + R[5] * v_xyz[2];
	v_enu[2] = R[6] * v_xyz[0] + R[7] * v_xyz[1] + R[8] * v_xyz[2];
}


// 计算 ENU 坐标的上仰角
double CoordinateTransform::UAngel(const double enu[3]) {
	return atan2(enu[2], sqrt(enu[0] * enu[0] + enu[1] * enu[1]));
}


//计算方位角
double CoordinateTransform::EAngel(const double enu[3])
{
	return atan2(enu[0], enu[1]);
}

double CoordinateTransform::Getyita2(const double& B)
{
	return defaultEllipsoid.e2_prime * cos(B) * cos(B);
}

double CoordinateTransform::Gett(double B)
{
	return tan(B);
}

double CoordinateTransform::GetM(double B)
{
	double M = 0.0;
	M = defaultEllipsoid.a * (1 - defaultEllipsoid.e2) * pow((1 - defaultEllipsoid.e2 * sin(B) * sin(B)), -1.5);
	return M;
}

double CoordinateTransform::GetN(double B)
{
	double N = 0.0;
	N = defaultEllipsoid.a / sqrt(1 - defaultEllipsoid.e2 * sin(B) * sin(B));
	return N;
}

double CoordinateTransform::GetMLength(double B)
{
	double m0 = defaultEllipsoid.a * (1 - defaultEllipsoid.e2);
	double m2 = 1.5 * defaultEllipsoid.e2 * m0;
	double m4 = 5.0 / 4.0 * defaultEllipsoid.e2 * m2;
	double m6 = 7.0 / 6.0 * defaultEllipsoid.e2 * m4;
	double m8 = 9.0 / 8.0 * defaultEllipsoid.e2 * m6;

	double a0 = m0 + m2 / 2.0 + 3.0 / 8.0 * m4 + 5.0 / 16.0 * m6 + 35.0 / 128.0 * m8;
	double a2 = m2 / 2.0 + m4 / 2.0 + 15.0 / 32.0 * m6 + 7.0 / 16.0 * m8;
	double a4 = m4 / 8.0 + 3.0 / 16.0 * m6 + 7.0 / 32.0 * m8;
	double a6 = m6 / 32.0 + m8 / 16.0;
	double a8 = m8 / 128.0;

	double X;
	X = a0 * B - a2 / 2.0 * sin(2 * B) + a4 / 4.0 * sin(4 * B) - a6 / 6.0 * sin(6 * B) + a8 / 8.0 * sin(8 * B);
	return X;
}

double CoordinateTransform::GetBf(double X)
{
	double m0 = defaultEllipsoid.a * (1 - defaultEllipsoid.e2);
	double m2 = 1.5 * defaultEllipsoid.e2 * m0;
	double m4 = 5.0 / 4.0 * defaultEllipsoid.e2 * m2;
	double m6 = 7.0 / 6.0 * defaultEllipsoid.e2 * m4;
	double m8 = 9.0 / 8.0 * defaultEllipsoid.e2 * m6;

	double a0 = m0 + m2 / 2.0 + 3.0 / 8.0 * m4 + 5.0 / 16.0 * m6 + 35.0 / 128.0 * m8;
	double a2 = m2 / 2.0 + m4 / 2.0 + 15.0 / 32.0 * m6 + 7.0 / 16.0 * m8;
	double a4 = m4 / 8.0 + 3.0 / 16.0 * m6 + 7.0 / 32.0 * m8;
	double a6 = m6 / 32.0 + m8 / 16.0;
	double a8 = m8 / 128.0;

	double B = 0;
	double Btemp = 0;
	Btemp = B;
	B = X / a0;
	do
	{
		B = Btemp;
		double Fun;
		Fun = -a2 / 2.0 * sin(2 * B) + a4 / 4.0 * sin(4 * B) - a6 / 6.0 * sin(6 * B) + a8 / 8.0 * sin(8 * B);
		double X_Fun;
		X_Fun = X - Fun;
		Btemp = X_Fun / a0;
	} while (fabs(Btemp - B) > 0.0000001);

	return Btemp;
}

/**
 * @file         "BaseGT.cpp"
 * @brief        实现给BLH坐标增加dxyz 东北天 ENU 平移后的 Delta BLH量
 * @param[in]    B
 * @param[in]    dxyz	东北天 平移量  m
 * @param[in]    dBLH	平移了 BLH
 * @return       double	
 * @note         
 * @par History
 *               2025/01/21, Ken Ga  new \n
 */
double CoordinateTransform::getDeltaBLH(const double B,const double dxyz[], double dBLH[])
{
	dBLH[0] = dxyz[1] / GetM(B);
	dBLH[1] = dxyz[0] / GetN(B) * cos(B);
	dBLH[2] = dxyz[2];
	return 0.0;
}

short UTCTIME::calcuWeekD(int year, int month, int day)
{
	if (month < 3) {
		month += 12;
		year -= 1;
	}
	int k = year % 100;
	int j = year / 100;
	int h = (day + 13 * (month + 1) / 5 + k + k / 4 + j / 4 - 2 * j) % 7;
	return (h + 5) % 7;  // 转换为 0（周日）到 6（周六）
}

short UTCTIME::getWeekD() const
{
	return UTCTIME::calcuWeekD(Year, Month, Day);
}

/**
 * @file         "BaseGT.cpp"
 * @brief        判断闰年
 * @param[in]    year
 * @return       bool
 * @note         
 * @par History
 *               2024/12/11, Ken Ga  new \n
 */
bool UTCTIME::IsLeapYear(int year)
{
	return (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0);
}

int UTCTIME::getDoy(const UTCTIME& common)
{
	int day_of_year = 0;

	// 根据年份，判断是否为闰年
	bool is_leap_year = IsLeapYear(common.Year);

	// 累加月份的天数
	for (int i = 0; i < common.Month - 1; i++)
	{
		if (i == 1 && is_leap_year)  // 如果是2月并且是闰年
		{
			day_of_year += 29;
		}
		else
		{
			day_of_year += days_in_month[i];
		}
	}

	// 加上当前月的天数
	day_of_year += common.WeekN + 1;  // 1月1日算第1天

	return day_of_year;
}

/**
 * @file         "BaseGT.cpp"
 * @brief        赋值函数
 * @param[in]    year
 * @param[in]    month
 * @param[in]    day
 * @param[in]    hour
 * @param[in]    min
 * @param[in]    second
 * @return       void
 * @note         
 * @par History
 *               2025/04/7, Ken Ga  new \n
 */
void UTCTIME::Mutator(const int& year, const int& month,const int &day ,const int& hour, const int& min, const double& second)
{
	Year = year;
	if (month > 12)return;
	Month = month;
	if (day > 31)return;
	Day = day;
	if (hour > 23)return;
	Hour = hour;
	if (min > 59)return;
	Minute = min;
	if (second >= 60.0)return;
	Second = second;
}

// 计算某年到1980年1月6日的天数差
int UTCTIME::DaysSince1980(int year)
{
	int days = 0;
	for (int y = 1980; y < year; y++)
	{
		days += UTCTIME::IsLeapYear(y) ? 366 : 365;
	}
	return days;
}


// 杆臂（INS到GNSS天线相位中心），单位：米
// 坐标顺序：前、右、下(x,y,z)
// 如果是：右前上，应该变为：y,x,-z
//static  double lever[3] = { 0.3, 0.3 ,-1.3 };//{ 0.3, 0.3 ,-1.3 };  1.317,0.385,-1.338 
// 将INS中心位置（BLH，单位：度、米）转换为GNSS天线相位中心位置（BLH）
/**
 * @file         "BaseGT.cpp"
 * @brief        
 * @param[in]    insBLH		INS惯导中心的大地坐标(度)
 * @param[in]    hpr		航向、俯仰、横滚(度)
 * @param[in]    newBLH		输出的天线相位中心的大地坐标(度)
 * @return       void
 * @note         
 * @par History
 *               2025/03/11, Ken Ga  new \n
 */
void ins_to_gnss(const double insBLH[3], const double hpr[3],const double lever[3], double newBLH[3])
{
	//lever[0] = 1.317;
	//lever[1] = 0.385;
	//lever[2] = -1.338;

	// WGS84椭球参数
	double a = 6378137.0;               // 长半轴，单位：米
	double f = 1.0 / 298.257223563;       // 扁率
	double e2 = 2 * f - f * f;            // 第一偏心率平方

	// 输入的INS经纬度（单位：度）和高度（单位：米）
	double lat_deg = insBLH[0];
	double lon_deg = insBLH[1];
	double h = insBLH[2];

	// 转换为弧度
	double lat = lat_deg * PI / 180.0;
	double lon = lon_deg * PI / 180.0;

	// 杆臂（INS到GNSS天线相位中心），单位：米
	// 坐标顺序：前、右、下
	//const double lever[3] = { 0.3, -0.3 ,1.3 };//{ 0.3, 0.3 ,-1.3 };  1.317,0.385,-1.338 

	// 姿态角（单位：度）：航向角H、俯仰角P、横滚角R
	double H = hpr[0] * PI / 180.0;   // 航向角（yaw）
	double P = hpr[1] * PI / 180.0;   // 俯仰角（pitch）
	double R = hpr[2] * PI / 180.0;   // 横滚角（roll）

	// 计算从INS机体系到NED坐标系的旋转矩阵（R_body_to_ned）
	// 假设INS机体系中：x向前，y向右，z向下
	double Rbn[3][3];
	Rbn[0][0] = cos(P) * cos(H);
	Rbn[0][1] = sin(R) * sin(P) * cos(H) - cos(R) * sin(H);
	Rbn[0][2] = cos(R) * sin(P) * cos(H) + sin(R) * sin(H);

	Rbn[1][0] = cos(P) * sin(H);
	Rbn[1][1] = sin(R) * sin(P) * sin(H) + cos(R) * cos(H);
	Rbn[1][2] = cos(R) * sin(P) * sin(H) - sin(R) * cos(H);

	Rbn[2][0] = -sin(P);
	Rbn[2][1] = sin(R) * cos(P);
	Rbn[2][2] = cos(R) * cos(P);

	// 将杆臂从INS机体系转换到NED坐标系，得到偏移量（单位：米）
	double offset_n = Rbn[0][0] * lever[0] + Rbn[0][1] * lever[1] + Rbn[0][2] * lever[2];  // 北向分量
	double offset_e = Rbn[1][0] * lever[0] + Rbn[1][1] * lever[1] + Rbn[1][2] * lever[2];  // 东向分量
	double offset_d = Rbn[2][0] * lever[0] + Rbn[2][1] * lever[1] + Rbn[2][2] * lever[2];  // 地向分量（向下为正）

	// 计算当地的卯酉圈曲率半径和子午圈曲率半径
	double sinLat = sin(lat);
	double cosLat = cos(lat);
	double N_val = a / sqrt(1 - e2 * sinLat * sinLat);               // 卯酉圈曲率半径
	double M_val = a * (1 - e2) / pow(1 - e2 * sinLat * sinLat, 1.5);    // 子午圈曲率半径

	// 利用小量近似计算经纬度的修正量（单位：弧度）
	double dLat = offset_n / (M_val + h);
	double dLon = offset_e / ((N_val + h) * cosLat);
	// 高度修正：NED中“向下”为正，所以新高度 = 原高度 - offset_d
	double new_lat = lat + dLat;
	double new_lon = lon + dLon;
	double new_h = h - offset_d;

	// 转换回度制输出
	newBLH[0] = new_lat * 180.0 / PI;
	newBLH[1] = new_lon * 180.0 / PI;
	newBLH[2] = new_h;
}


// 
/**
 * @file         "BaseGT.cpp"
 * @brief        
 * @param[in]    gnssBLH  BLH，单位：度、米
 * @param[in]    hpr
 * @param[in]    lever	  INS-天线相位中心的矢量
 * @param[in]    insBLH   BLH，单位：度、米
 * @return       void
 * @note         将GNSS天线相位中心（BLH，单位：度、米）转换到INS中心位置（BLH）（INS-天线相位中心的矢量）
 * @par History
 *               2025/03/12, Ken Ga  new \n
 */
void gnss_to_ins(const double gnssBLH[3], const double hpr[3], const double lever[3], double insBLH[3])
{
	//lever[0] = -0.077;
	//lever[1] = 0.028;
	//lever[2] = 0.052;
	double Rlever[3]{};

	Rlever[0] = -lever[0];
	Rlever[1] = -lever[1];
	Rlever[2] = -lever[2];


	ins_to_gnss(gnssBLH, hpr,Rlever ,insBLH);
}
static void MatrixMultiply(int r1, int c1, const double M1[], int r2, int c2, const double M2[], double M3[], double scale)
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


static void TransMatCbn(const double angles[3], double Cnb[]) {
	double psi = angles[0] * D2R;    // 航向角 (Yaw)
	double theta = angles[1] * D2R;  // 俯仰角 (Pitch)
	double phi = angles[2] * D2R;    // 横滚角 (Roll)

	double cpsi = cos(psi), spsi = sin(psi);
	double ctheta = cos(theta), stheta = sin(theta);
	double cphi = cos(phi), sphi = sin(phi);

	// 构造方向余弦矩阵 C_b^n
	
	Cnb[0] = cpsi * ctheta; Cnb[1] = spsi * ctheta; Cnb[2] = -stheta;
	Cnb[3] = cpsi * stheta * sphi - spsi * cphi; Cnb[4] = spsi * stheta * sphi + cpsi * cphi; Cnb[5] = ctheta * sphi;
	Cnb[6] = cpsi * stheta * cphi + spsi * sphi; Cnb[7] = spsi * stheta * cphi - cpsi * sphi; Cnb[8] = ctheta * cphi;
}



// b系转n系的矩阵，输入是航向，俯仰，横滚，输出转移矩阵R
static void rotationMatrix(const double att[3], double R[9]) {
	double pitch = att[0];
	double roll = att[1];
	double yaw = att[2];

	// 绕Z轴（yaw）的旋转矩阵
	double R_yaw[9] = {
		cos(yaw), -sin(yaw), 0,
		sin(yaw), cos(yaw), 0,
		0, 0, 1
	};

	// 绕Y轴（pitch）的旋转矩阵
	double R_pitch[9] = {
		1, 0, 0,
		0, cos(pitch), -sin(pitch),
		0, sin(pitch), cos(pitch)
	};

	// 绕X轴（roll）的旋转矩阵
	double R_roll[9] = {
		cos(roll), 0, sin(roll),
		0, 1, 0,
		-sin(roll), 0, cos(roll)
	};
	double temp[9];
	MatrixMultiply(3, 3, R_yaw, 3, 3, R_pitch, temp, 1.0);
	MatrixMultiply(3, 3, temp, 3, 3, R_roll, R, 1.0);
}


static bool MatrixInv(int r, int c, double M[], double M_Inv[])
{
	//int* is = (int*)malloc(r, sizeof(int));
	//int* js = (int*)malloc(r, sizeof(int));
	int is[50] = { 0 };
	int js[50] = { 0 };
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

		if (d < 1e-10)
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
 * @file         "BaseGT.cpp"
 * @brief        输出ENU，输出右前上，输入NED，输出前右下
 * @param[in]    n
 * @param[in]    angles
 * @param[in]    b
 * @return       void
 * @note         
 * @par History
 *               2025/03/11, Ken Ga  new \n
 */
void TransformNtoB(const double n[3], const double angles[3], double b[3]) {
	double C_b_n[9], Cnb[9];
	//TransMatCbn(angles, C_b_n);
	rotationMatrix(angles, C_b_n);
	MatrixInv(3, 3, C_b_n, Cnb);
	// 计算 b = C_b^n * n
	MatrixMultiply(3, 3, Cnb, 3, 1, n, b,1.0);
}

