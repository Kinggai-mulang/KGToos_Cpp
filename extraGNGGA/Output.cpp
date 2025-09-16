#include"BaseGT.h"
// 仅用于校验码，作为静态函数服务于 void PrintNEMA()
static int NMEACRC(char* buff, int len)
{
	int checknum = 0;
	char* pdata = NULL;
	short i;

	if (!buff)return 0;

	pdata = buff;
	for (i = 0; i < len; i++)
	{
		checknum ^= *(pdata + i);
	}
	return checknum;
}

/**
 * @brief        将结果用NEMA格式进行输出，一次输出一个历元的结果
 * @param[in]    gt			当前历元GPS时间结构体
 * @param[in]    isFixed	传入的坐标是否为固定解坐标
 * @param[in]    XYZ		当前解算的位置坐标
 * @param[in]    nsat		共视卫星数  BDS+GPS
 * @param[in]    OutputFilePath	  输出的文件路径
 * @return       void
 * @note         注意这里GPS周秒转成通用时，ECEF坐标转LLH使用自己的所写的代码实现
 * @par History
 *               2025/03/12, 卫星导航算法课程  new \n
 */
static void PrintNEMA(GPSTIME gt, bool isFixed, const double XYZ[], const int nsat, const char OutputFilePath[]) {
	if (!OutputFilePath)return;
	char   sol[1024] = "";
	memset(sol, 0, 1024);
	double BLH[3]{}, min;
	int    checkNum = 0, checkStart = 0, checkEnd = 0, deg = 0, solLen = 0, QNMEA = 5;
	static bool bHead1 = true;
	UTCTIME ymd;			// 年月日时分秒的转换
	ymd = GPSTtoUTC(gt);
	CoordinateTransform::XYZ2BLH(XYZ, BLH);      //XYZ是解算出流动站ecef坐标,默认输出的BLH是弧度制
	isFixed ? QNMEA = 4 : QNMEA = 5;   //这里是针对RTK，如果是SPP结果，设QNMEA为1即可

	if (bHead1)
	{
		solLen += sprintf(sol + solLen, "$ACCUR,1,12,0,%04d%02d%02d,,123206.00,,,,,,,,,,*", ymd.Year, ymd.Month, ymd.Day);
		checkNum = NMEACRC(sol + 1, solLen - 2); // 校验和
		solLen += sprintf(sol + solLen, "%X\n", checkNum);
		bHead1 = false;
	}

	checkStart = solLen + 1; //起始下标
	checkEnd = solLen + 1; //结束下标

	solLen += sprintf(sol + solLen, "$GNGGA,");

	// 时间: 时分秒
	solLen += sprintf(sol + solLen, "%02d%02d%05.2f,", ymd.Hour, ymd.Minute, ymd.Second);

	// 纬度 ddmm.mmmmmm 度分 
	deg = (int)(fabs(BLH[0]) * R2D);
	min = (fabs(BLH[0] * R2D) - deg) * 60;
	if (BLH[0] > 0) solLen += sprintf(sol + solLen, "%d%09.6f,N,", deg, min);
	else solLen += sprintf(sol + solLen, "%d%09.6f,S,", deg, min);

	// 经度 dddmm.mmmmmm 度分
	deg = (int)(fabs(BLH[1]) * R2D);
	min = (fabs(BLH[1] * R2D) - deg) * 60;
	if (BLH[1] > 0) solLen += sprintf(sol + solLen, "%d%09.6f,E,", deg, min);
	else solLen += sprintf(sol + solLen, "%d%09.6f,W,", deg, min);

	// 解状态
	solLen += sprintf(sol + solLen, "%d,", QNMEA);

	// 卫星数
	solLen += sprintf(sol + solLen, "%2d,", nsat);

	// HDOP
	solLen += sprintf(sol + solLen, "%3.1f,", 0.0);

	// 高程
	solLen += sprintf(sol + solLen, "%5.3f,M,", BLH[2]);

	// 高程异常
	solLen += sprintf(sol + solLen, "%5.3f,M,", 0.0);

	solLen += sprintf(sol + solLen, ",*");

	// 校验和：$与*之间所有字符ASCII码的校验和(各字节做异或运算，得到校验和后，再转换16进制格式的ASCII字符。)
	checkEnd = solLen - 2; //-2:*与下标
	checkNum = NMEACRC(sol + checkStart, checkEnd - checkStart + 1);
	solLen += sprintf(sol + solLen, "%X\n", checkNum);


	static FILE* fpNEMA = NULL;
	if (fpNEMA == NULL)
	{
		fopen_s(&fpNEMA, OutputFilePath, "w");
	}
	if (fpNEMA != NULL) {
		fprintf(fpNEMA, sol);
	}
}