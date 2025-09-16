#include "BaseGT.h"     //����ʱ�䡢�����


// �����ⶨ�徲̬��Ա��������ʼ��Ϊ WGS84
Ellipsoid CoordinateTransform::defaultEllipsoid = WGS84;   //������CPP�ļ��ж���

/**
 * @file         "BaseGT.cpp"
 * @brief        ʵ�ֶԲ�������+�ŵ�����
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

	// �������С�����ֵĽ�λ
	if (temp.fracSec >= 1.0)
	{
		temp.fracSec -= 1.0;
		temp.secWeek += 1.0;
	}

	// ����������Ľ�λ (һ����604800��)
	if (temp.secWeek >= SECOND_WEEK)
	{
		temp.secWeek -= SECOND_WEEK;
		temp.week += 1;
	}

	return temp;
}


/**
 * @file         "BaseGT.cpp"
 * @brief        ����������أ�ʵ��GPSʱ���
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

	// �������С�����ֵĽ�λ
	if (temp.fracSec < 0.0)
	{
		temp.fracSec += 1.0;
		temp.secWeek -= 1.0;
	}

	// ����������Ľ�λ
	if (temp.secWeek < 0.0)
	{
		temp.secWeek += SECOND_WEEK;
		temp.week -= 1;
	}

	return temp;
}

/**
 * @file         "BaseGT.cpp"
 * @brief        GPSʱ��-double����
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



std::ostream& operator<<(std::ostream& out, GPSTIME t)  //����<<��ʹ����ֱ�����GPSTIME��Ϊ�˱�֤����ʽ��̣�ʹ�÷��������޸�
{
	out << t.week << "  " << t.secWeek + t.fracSec;
	return out;
}

bool operator==(GPSTIME& a, GPSTIME& b)  //С������С��1e-14����
{
	if (a.week != b.week || a.secWeek != b.secWeek)return false;

	return (((a.week - b.week) * 604800 + a.secWeek + a.fracSec - b.secWeek + b.fracSec) < Epsilon);
}

bool operator!=(GPSTIME& a, GPSTIME& b)  //С������С��1e-14����
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
	// 1. GPSʱ����ʼ������ 1980��1��6��
	int base_year = 1980;
	int base_month = 1;
	int base_day = 6;

	// 2. �����GPS��㵽��ǰ��������������
	int total_days = gps.week * 7;   // ����ת��Ϊ����
	double total_seconds = gps.secWeek + gps.fracSec; // ���ڵ�������

	// 3. ����������Ϊʱ���֡���
	int hours = static_cast<int>(total_seconds) / 3600;
	total_seconds -= hours * 3600;
	int minutes = static_cast<int>(total_seconds) / 60;
	double seconds = total_seconds - minutes * 60;

	// 4. �����������ۼ�
	total_days += hours / 24;        // ��Сʱ���е����첿�ּӵ�������
	hours %= 24;                     // ʣ�µ�Сʱ��

	// 5. ��1980��1��6�����ۼ���Щ�������ҵ���Ӧ���ꡢ�¡���
	int year = base_year;
	int month = base_month;
	int day = base_day + total_days;

	// 6. �������ڵ�����Χ
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

	// 7. �����ܼ� (0 ������, 1-6 ��Ӧ��һ������)
	int total_weeks = gps.week;
	utc_t.WeekN = (total_weeks * 7 + static_cast<int>(gps.secWeek) / 86400) % 7;

	// 8. ��� COMMONTIME �ṹ��
	utc_t.Year = year;
	utc_t.Month = month;
	utc_t.Day = day;
	utc_t.Hour = hours;
	utc_t.Minute = minutes;
	utc_t.Second = seconds;

	return utc_t;
}

// �������� 18s
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
* @brief        BDSʱתGPSʱ,Ĭ��ֻ����14s����
* @param[in]    bdst        GPSTIME    BDSʱ��
* @param[in]    bSecCorr    bool       true = ֻ����14s����, false = BDS�ܸ���
* @return       GPSTIME     GPSʱ��
* @note         ��
* @internals    ��
*/
GPSTIME BDST2GPST(GPSTIME bdst, bool bSecCorr)
{
	if (bSecCorr == false) bdst.week += 1356;
	return (bdst + 14.0);
}

/**
 * @file         "BaseGT.cpp"
 * @brief        ������ת������
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
// * @brief        �����ַ�����תΪdouble���͵�����
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
//	// ʹ���ַ������ָ��ַ���
//	std::istringstream iss(input);
//	std::string token;
//	int index = 0;
//
//	while (std::getline(iss, token, ',')) {
//		// ��ÿ�����ַ���ת��Ϊ double
//		data[index++] = std::stod(token);
//	}
//	return index;
//}
//
//std::string D2str(double X, const int pres)
//{
//	// ����һ���ַ�����
//	std::ostringstream oss;
//	// ���ù̶���С��λ��
//	oss << std::fixed << std::setprecision(pres) << X;
//	// �������е��ַ���
//	return oss.str();
//}


/**
 * @file         "BaseGT.cpp"
 * @brief        ���С���������
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
 * @brief        �����1986��1��6�յ����ڵ���������
 * @return       double
 * @note         ����������� Week*604800+sec
 * @par History
 *               2024/12/11, Ken Ga  new \n
 */
double GPSTIME::getSec() const
{
	return week * SECOND_WEEK + secWeek + fracSec;
}

/**
 * @file         "BaseGT.cpp"
 * @brief        ���������
 * @return       double
 * @note         �������������
 * @par History
 *               2024/12/11, Ken Ga  new \n
 */
double GPSTIME::getSow() const
{
	return secWeek + fracSec;
}

/**
 * @file         "BaseGT.cpp"
 * @brief        ����Ŀ��ʱ�������
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
 * @brief        ����Ŀ��ʱ�������
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
 * @brief        ��Ŀ��ʱ�����������
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
 * @brief        ��ֵ����
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
 * @brief        ѡ���Ƿ�Ϊ�ص��ʱ��
 * @param[in]    Isec	��λ��ʱ��(������)
 * @param[in]    pres	��λ�ľ���(�ڸþ����µ�ʱ��ᱻ��λ)
 * @return       bool
 * @note         ѡ���Ƿ�Ϊ�ص��ʱ��
 * @par History
 *               2025/01/16, Ken Ga  new \n
 */
bool GPSTIME::isSpecify(const double Isec, const double pres)
{
	return (fabs(secWeek + fracSec - Isec)) < pres;
}

/**
* @brief        ����ʱ��,�Ͳο�ʱ����������
* @param[in]    t        GPSTIME    ������ʱ��
* @param[in]    t0       GPSTIME    �ο�ʱ��
* @return       GPSTIME  �������ʱ��
* @note         ��
* @internals    ��
*/
GPSTIME GPSTIME::adjweek(GPSTIME t, GPSTIME t0)
{
	double tt = (t - t0).getSec();
	if (tt < -HALF_SECOND_WEEK) return (t + SECOND_WEEK);
	if (tt > HALF_SECOND_WEEK) return (t - SECOND_WEEK);
	return t;
}


/*****************************************************
����ת��ģ��
*****************************************************/
/**
* @brief       WGS84����ת����CGCS2000����,�����Ӧ���ǿ����Ԫ
* @param[in]   XYZ_WGS84         double  WGS84����(m)
* @param[out]  XYZ_CGCS2000      double  CGCS2000����(m)
* @return      void
* @note        WGS84��׼��ITRF2000���,CGCS2000��׼��ITRF1997���,��վ��������ת��
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

	// ���ü�д��
	for (int i = 0; i < 3; i++)
	{
		XYZ_CGCS2000[i] = (1.0 + D) * XYZ_WGS84[i] + T[i];
	}
}


/**
	* @brief       CGCS2000����ת����WGS84����,�����Ӧ���ǿ����Ԫ
	* @param[in]   XYZ_CGCS2000      double  CGCS2000����(m)
	* @param[out]  XYZ_WGS84         double  WGS84����(m)
	* @return      void
	* @note        WGS84��׼��ITRF2000���,CGCS2000��׼��ITRF1997���,��վ��������ת��
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
	// ���ü�д��
	for (int i = 0; i < 3; i++)
	{
		XYZ_WGS84[i] = (XYZ_CGCS2000[i] - T[i]) / (1.0 + D);
	}
}



/**
* @brief        WGS84����ϵ: XYZ����תLLH����
* @param[in]    XYZ        double    XYZ(m) ����
* @param[in]    CoorSys    int       ����ϵͳ: 0:WGS84, 1:CGCS2000
* @param[out]   LLH        double    ����
* @return       void
* @note         ��
* @internals    ��
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
* @brief        WGS84����ϵ: LLH����תXYZ����
* @param[in]    LLH        double    Lat[-90,90],Lon[-180,180],Hgt (deg,m),����ӦΪ����
* @param[in]    CoorSys    int       ����ϵͳ: 0:WGS84, 1:CGCS2000
* @param[out]   XYZ        double    XYZ(m) ����
* @return       void
* @note         ����ӦΪ����
* @internals    ��
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

// BLH �� ENU ת��
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

// XYZ �� ENU ת��
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
* @brief        ENUת��ECEFϵ����ת����
* @param[in]    LLH      double    γ��,����,��ظ�
* @param[out]   R_L2E    double    ENU->ECEF����ת����
* @return       void
* @note			�����ǽ��ο������BLH���룬�����ת����
* @internals    ��
*/
void CoordinateTransform::Rxyz2enu(const double LLH[3], double R_E2L[9])
{
	double sinp = sin(LLH[0]), cosp = cos(LLH[0]), sinl = sin(LLH[1]), cosl = cos(LLH[1]);
	R_E2L[0] = -sinl;        R_E2L[1] = cosl;         R_E2L[2] = 0.0;
	R_E2L[3] = -sinp * cosl; R_E2L[4] = -sinp * sinl; R_E2L[5] = cosp;
	R_E2L[6] = cosp * cosl;  R_E2L[7] = cosp * sinl;  R_E2L[8] = sinp;
}

/**
* @brief        ECEF������ת��ENUϵ��
* @param[in]    LLH      double    γ��,����,��ظ�,������
* @param[in]    v_xyz    double    ECEFϵ������
* @param[out]   v_enu    double    ENUϵ������
* @return		��ECEF����תΪENU��LLHΪ�ο�վ������,boolΪtrueʱΪBLH��ΪfalseʱΪXYZ
* @note         ��
* @internals    ��
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


// ���� ENU �����������
double CoordinateTransform::UAngel(const double enu[3]) {
	return atan2(enu[2], sqrt(enu[0] * enu[0] + enu[1] * enu[1]));
}


//���㷽λ��
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
 * @brief        ʵ�ָ�BLH��������dxyz ������ ENU ƽ�ƺ�� Delta BLH��
 * @param[in]    B
 * @param[in]    dxyz	������ ƽ����  m
 * @param[in]    dBLH	ƽ���� BLH
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
	return (h + 5) % 7;  // ת��Ϊ 0�����գ��� 6��������
}

short UTCTIME::getWeekD() const
{
	return UTCTIME::calcuWeekD(Year, Month, Day);
}

/**
 * @file         "BaseGT.cpp"
 * @brief        �ж�����
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

	// ������ݣ��ж��Ƿ�Ϊ����
	bool is_leap_year = IsLeapYear(common.Year);

	// �ۼ��·ݵ�����
	for (int i = 0; i < common.Month - 1; i++)
	{
		if (i == 1 && is_leap_year)  // �����2�²���������
		{
			day_of_year += 29;
		}
		else
		{
			day_of_year += days_in_month[i];
		}
	}

	// ���ϵ�ǰ�µ�����
	day_of_year += common.WeekN + 1;  // 1��1�����1��

	return day_of_year;
}

/**
 * @file         "BaseGT.cpp"
 * @brief        ��ֵ����
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

// ����ĳ�굽1980��1��6�յ�������
int UTCTIME::DaysSince1980(int year)
{
	int days = 0;
	for (int y = 1980; y < year; y++)
	{
		days += UTCTIME::IsLeapYear(y) ? 366 : 365;
	}
	return days;
}


// �˱ۣ�INS��GNSS������λ���ģ�����λ����
// ����˳��ǰ���ҡ���(x,y,z)
// ����ǣ���ǰ�ϣ�Ӧ�ñ�Ϊ��y,x,-z
//static  double lever[3] = { 0.3, 0.3 ,-1.3 };//{ 0.3, 0.3 ,-1.3 };  1.317,0.385,-1.338 
// ��INS����λ�ã�BLH����λ���ȡ��ף�ת��ΪGNSS������λ����λ�ã�BLH��
/**
 * @file         "BaseGT.cpp"
 * @brief        
 * @param[in]    insBLH		INS�ߵ����ĵĴ������(��)
 * @param[in]    hpr		���򡢸��������(��)
 * @param[in]    newBLH		�����������λ���ĵĴ������(��)
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

	// WGS84�������
	double a = 6378137.0;               // �����ᣬ��λ����
	double f = 1.0 / 298.257223563;       // ����
	double e2 = 2 * f - f * f;            // ��һƫ����ƽ��

	// �����INS��γ�ȣ���λ���ȣ��͸߶ȣ���λ���ף�
	double lat_deg = insBLH[0];
	double lon_deg = insBLH[1];
	double h = insBLH[2];

	// ת��Ϊ����
	double lat = lat_deg * PI / 180.0;
	double lon = lon_deg * PI / 180.0;

	// �˱ۣ�INS��GNSS������λ���ģ�����λ����
	// ����˳��ǰ���ҡ���
	//const double lever[3] = { 0.3, -0.3 ,1.3 };//{ 0.3, 0.3 ,-1.3 };  1.317,0.385,-1.338 

	// ��̬�ǣ���λ���ȣ��������H��������P�������R
	double H = hpr[0] * PI / 180.0;   // ����ǣ�yaw��
	double P = hpr[1] * PI / 180.0;   // �����ǣ�pitch��
	double R = hpr[2] * PI / 180.0;   // ����ǣ�roll��

	// �����INS����ϵ��NED����ϵ����ת����R_body_to_ned��
	// ����INS����ϵ�У�x��ǰ��y���ң�z����
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

	// ���˱۴�INS����ϵת����NED����ϵ���õ�ƫ��������λ���ף�
	double offset_n = Rbn[0][0] * lever[0] + Rbn[0][1] * lever[1] + Rbn[0][2] * lever[2];  // �������
	double offset_e = Rbn[1][0] * lever[0] + Rbn[1][1] * lever[1] + Rbn[1][2] * lever[2];  // �������
	double offset_d = Rbn[2][0] * lever[0] + Rbn[2][1] * lever[1] + Rbn[2][2] * lever[2];  // �������������Ϊ����

	// ���㵱�ص�î��Ȧ���ʰ뾶������Ȧ���ʰ뾶
	double sinLat = sin(lat);
	double cosLat = cos(lat);
	double N_val = a / sqrt(1 - e2 * sinLat * sinLat);               // î��Ȧ���ʰ뾶
	double M_val = a * (1 - e2) / pow(1 - e2 * sinLat * sinLat, 1.5);    // ����Ȧ���ʰ뾶

	// ����С�����Ƽ��㾭γ�ȵ�����������λ�����ȣ�
	double dLat = offset_n / (M_val + h);
	double dLon = offset_e / ((N_val + h) * cosLat);
	// �߶�������NED�С����¡�Ϊ���������¸߶� = ԭ�߶� - offset_d
	double new_lat = lat + dLat;
	double new_lon = lon + dLon;
	double new_h = h - offset_d;

	// ת���ض������
	newBLH[0] = new_lat * 180.0 / PI;
	newBLH[1] = new_lon * 180.0 / PI;
	newBLH[2] = new_h;
}


// 
/**
 * @file         "BaseGT.cpp"
 * @brief        
 * @param[in]    gnssBLH  BLH����λ���ȡ���
 * @param[in]    hpr
 * @param[in]    lever	  INS-������λ���ĵ�ʸ��
 * @param[in]    insBLH   BLH����λ���ȡ���
 * @return       void
 * @note         ��GNSS������λ���ģ�BLH����λ���ȡ��ף�ת����INS����λ�ã�BLH����INS-������λ���ĵ�ʸ����
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
	double psi = angles[0] * D2R;    // ����� (Yaw)
	double theta = angles[1] * D2R;  // ������ (Pitch)
	double phi = angles[2] * D2R;    // ����� (Roll)

	double cpsi = cos(psi), spsi = sin(psi);
	double ctheta = cos(theta), stheta = sin(theta);
	double cphi = cos(phi), sphi = sin(phi);

	// ���췽�����Ҿ��� C_b^n
	
	Cnb[0] = cpsi * ctheta; Cnb[1] = spsi * ctheta; Cnb[2] = -stheta;
	Cnb[3] = cpsi * stheta * sphi - spsi * cphi; Cnb[4] = spsi * stheta * sphi + cpsi * cphi; Cnb[5] = ctheta * sphi;
	Cnb[6] = cpsi * stheta * cphi + spsi * sphi; Cnb[7] = spsi * stheta * cphi - cpsi * sphi; Cnb[8] = ctheta * cphi;
}



// bϵתnϵ�ľ��������Ǻ��򣬸�������������ת�ƾ���R
static void rotationMatrix(const double att[3], double R[9]) {
	double pitch = att[0];
	double roll = att[1];
	double yaw = att[2];

	// ��Z�ᣨyaw������ת����
	double R_yaw[9] = {
		cos(yaw), -sin(yaw), 0,
		sin(yaw), cos(yaw), 0,
		0, 0, 1
	};

	// ��Y�ᣨpitch������ת����
	double R_pitch[9] = {
		1, 0, 0,
		0, cos(pitch), -sin(pitch),
		0, sin(pitch), cos(pitch)
	};

	// ��X�ᣨroll������ת����
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

			for (i = 0; i < r; i++)   M[i * (r + 1)] += 1E-15; // M����ӽ�����ʱ,���Խ��߼���һ��С��
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
 * @brief        ���ENU�������ǰ�ϣ�����NED�����ǰ����
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
	// ���� b = C_b^n * n
	MatrixMultiply(3, 3, Cnb, 3, 1, n, b,1.0);
}

