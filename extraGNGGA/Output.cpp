#include"BaseGT.h"
// ������У���룬��Ϊ��̬���������� void PrintNEMA()
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
 * @brief        �������NEMA��ʽ���������һ�����һ����Ԫ�Ľ��
 * @param[in]    gt			��ǰ��ԪGPSʱ��ṹ��
 * @param[in]    isFixed	����������Ƿ�Ϊ�̶�������
 * @param[in]    XYZ		��ǰ�����λ������
 * @param[in]    nsat		����������  BDS+GPS
 * @param[in]    OutputFilePath	  ������ļ�·��
 * @return       void
 * @note         ע������GPS����ת��ͨ��ʱ��ECEF����תLLHʹ���Լ�����д�Ĵ���ʵ��
 * @par History
 *               2025/03/12, ���ǵ����㷨�γ�  new \n
 */
static void PrintNEMA(GPSTIME gt, bool isFixed, const double XYZ[], const int nsat, const char OutputFilePath[]) {
	if (!OutputFilePath)return;
	char   sol[1024] = "";
	memset(sol, 0, 1024);
	double BLH[3]{}, min;
	int    checkNum = 0, checkStart = 0, checkEnd = 0, deg = 0, solLen = 0, QNMEA = 5;
	static bool bHead1 = true;
	UTCTIME ymd;			// ������ʱ�����ת��
	ymd = GPSTtoUTC(gt);
	CoordinateTransform::XYZ2BLH(XYZ, BLH);      //XYZ�ǽ��������վecef����,Ĭ�������BLH�ǻ�����
	isFixed ? QNMEA = 4 : QNMEA = 5;   //���������RTK�������SPP�������QNMEAΪ1����

	if (bHead1)
	{
		solLen += sprintf(sol + solLen, "$ACCUR,1,12,0,%04d%02d%02d,,123206.00,,,,,,,,,,*", ymd.Year, ymd.Month, ymd.Day);
		checkNum = NMEACRC(sol + 1, solLen - 2); // У���
		solLen += sprintf(sol + solLen, "%X\n", checkNum);
		bHead1 = false;
	}

	checkStart = solLen + 1; //��ʼ�±�
	checkEnd = solLen + 1; //�����±�

	solLen += sprintf(sol + solLen, "$GNGGA,");

	// ʱ��: ʱ����
	solLen += sprintf(sol + solLen, "%02d%02d%05.2f,", ymd.Hour, ymd.Minute, ymd.Second);

	// γ�� ddmm.mmmmmm �ȷ� 
	deg = (int)(fabs(BLH[0]) * R2D);
	min = (fabs(BLH[0] * R2D) - deg) * 60;
	if (BLH[0] > 0) solLen += sprintf(sol + solLen, "%d%09.6f,N,", deg, min);
	else solLen += sprintf(sol + solLen, "%d%09.6f,S,", deg, min);

	// ���� dddmm.mmmmmm �ȷ�
	deg = (int)(fabs(BLH[1]) * R2D);
	min = (fabs(BLH[1] * R2D) - deg) * 60;
	if (BLH[1] > 0) solLen += sprintf(sol + solLen, "%d%09.6f,E,", deg, min);
	else solLen += sprintf(sol + solLen, "%d%09.6f,W,", deg, min);

	// ��״̬
	solLen += sprintf(sol + solLen, "%d,", QNMEA);

	// ������
	solLen += sprintf(sol + solLen, "%2d,", nsat);

	// HDOP
	solLen += sprintf(sol + solLen, "%3.1f,", 0.0);

	// �߳�
	solLen += sprintf(sol + solLen, "%5.3f,M,", BLH[2]);

	// �߳��쳣
	solLen += sprintf(sol + solLen, "%5.3f,M,", 0.0);

	solLen += sprintf(sol + solLen, ",*");

	// У��ͣ�$��*֮�������ַ�ASCII���У���(���ֽ���������㣬�õ�У��ͺ���ת��16���Ƹ�ʽ��ASCII�ַ���)
	checkEnd = solLen - 2; //-2:*���±�
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